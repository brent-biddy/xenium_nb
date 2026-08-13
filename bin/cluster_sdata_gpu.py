#!/usr/bin/env python3
"""
cluster_sdata_gpu.py - GPU-accelerated QC, normalize, and cluster a SpatialData zarr.

Mirrors cluster_sdata.py but uses rapids-singlecell for the compute-heavy steps
(QC, normalization, PCA, neighbors, UMAP, and a sweep of Leiden clusterings).
Data is moved back to CPU before zarr I/O.

Requires an Apptainer/Docker image with rapids-singlecell and a CUDA-capable GPU.

Writes clustered.zarr into the current working directory.

Usage:
    cluster_sdata_gpu.py --sample ROI1_A --path /data/ROI1_A.zarr
    cluster_sdata_gpu.py --sample ROI1_A --path /data/ROI1_A.zarr --resolutions 0.5 1.0
"""

import argparse

import rapids_singlecell as rsc
import spatialdata

from sdata_io import write_table_only
from timer import timer, timing_summary

DEFAULT_RESOLUTIONS = [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0]

# Must stay in step with cluster_sdata.py — the CPU and GPU steps are meant to be
# interchangeable, so a different default cut or recipe constant here would make them
# silently disagree. See that file for why these values were chosen.
DEFAULT_MIN_COUNTS = 20
DEFAULT_MIN_CELLS = 100
# See cluster_sdata.py: off by default, because a top-percentile cut on ovary
# removes oocytes, which are the most transcript-rich cells in the tissue.
DEFAULT_MAX_COUNTS_QUANTILE = 0

HVG_N_TOP_GENES = 2000
SCALE_MAX_VALUE = 10
N_PCS = 30
NEIGHBORS_METRIC = "cosine"
# n_iterations is deliberately not set; see cluster_sdata.py's Leiden call.


def parse_args():
    parser = argparse.ArgumentParser(
        description="GPU-accelerated clustering of a SpatialData zarr"
    )
    parser.add_argument("--sample", required=True, help="Sample identifier")
    parser.add_argument("--path", required=True, help="Path to input SpatialData zarr")
    parser.add_argument(
        "--resolutions",
        type=float,
        nargs="+",
        default=DEFAULT_RESOLUTIONS,
        metavar="RES",
        help="Leiden resolutions to sweep; one obs column is written per value "
        f"(default: {' '.join(str(r) for r in DEFAULT_RESOLUTIONS)}).",
    )
    parser.add_argument(
        "--min_counts",
        type=int,
        default=DEFAULT_MIN_COUNTS,
        help="Drop cells with fewer than this many transcripts. Read the cut off the "
        f"retention curves in the qc_report deck (default: {DEFAULT_MIN_COUNTS}).",
    )
    parser.add_argument(
        "--min_cells",
        type=int,
        default=DEFAULT_MIN_CELLS,
        help="Drop genes detected in fewer than this many cells "
        f"(default: {DEFAULT_MIN_CELLS}).",
    )
    parser.add_argument(
        "--max_counts_quantile",
        type=float,
        default=DEFAULT_MAX_COUNTS_QUANTILE,
        help="Drop cells above this quantile of transcript_counts, to remove doublets "
        f"and segmentation merges. 0 disables (default: {DEFAULT_MAX_COUNTS_QUANTILE}).",
    )
    return parser.parse_args()


def main():
    args = parse_args()

    resolutions = sorted(args.resolutions)

    output_path = "clustered.zarr"

    print(f"Sample:  {args.sample}")
    print(f"Input:   {args.path}")
    print(f"Output:  {output_path}")
    print(f"Res:     {', '.join(f'{r:g}' for r in resolutions)}")
    print(f"Filter:  min_counts={args.min_counts}, min_cells={args.min_cells}, "
          f"max_counts_quantile={args.max_counts_quantile:g}")

    with timer("Read zarr"):
        sdata = spatialdata.read_zarr(args.path)

    table_key = "table"
    with timer("Extract table"):
        adata = sdata.tables[table_key].copy()

    print(f"Table:   {adata.n_obs:,} cells × {adata.n_vars:,} genes  (key: '{table_key}')")

    with timer("Move to GPU"):
        rsc.get.anndata_to_GPU(adata)

    with timer("Filter"):
        # No calculate_qc_metrics call here, matching cluster_sdata — create_sdata
        # already annotated every per-cell and per-gene metric, and recomputing them
        # overwrote the Xenium-native obs["total_counts"] with a plain row sum. See
        # cluster_sdata.py's Filter block for the full argument.
        #
        # Dropping it is safe: rsc.pp.filter_cells and filter_genes both derive their
        # thresholds from X via _basic_qc and never read obs, so the cut is unchanged.
        n_before = adata.n_obs
        # Quantile on the UNFILTERED table, before min_counts removes the low tail —
        # see cluster_sdata.py's Filter block.
        max_counts = (
            float(adata.obs["transcript_counts"].quantile(args.max_counts_quantile))
            if args.max_counts_quantile else None
        )
        rsc.pp.filter_cells(adata, min_counts=args.min_counts)
        if max_counts is not None:
            rsc.pp.filter_cells(adata, max_counts=max_counts)
            print(f"Upper cut at q{args.max_counts_quantile:g} = {max_counts:,.0f} transcripts.")
        # Masks rather than removes — see cluster_sdata.py for why a detection cut on a
        # targeted panel deletes exactly the rare-cell-type markers annotation needs.
        passing, _ = rsc.pp.filter_genes(adata, min_cells=args.min_cells, inplace=False)
        adata.var["passes_min_cells"] = passing

    print(f"Filtered {n_before - adata.n_obs:,} low-quality cells.")
    print(f"Retained {adata.n_obs:,} cells × {adata.n_vars:,} genes.")

    with timer("HVG"):
        # seurat_v3 on raw counts, before normalization — see cluster_sdata.py.
        rsc.pp.highly_variable_genes(
            adata, flavor="seurat_v3", n_top_genes=HVG_N_TOP_GENES)
        # Re-select among the genes passing min_cells — see cluster_sdata.py.
        ranked = adata.var.loc[passing, "variances_norm"].nlargest(HVG_N_TOP_GENES).index
        adata.var["highly_variable"] = adata.var_names.isin(ranked)

    with timer("Normalize"):
        adata.layers["counts"] = adata.X.copy()
        # Normalization deliberately takes NO target_sum, so scanpy's default applies: the
        # median pre-normalization cell total. That default is principled here, not arbitrary.
        # A fixed target far from the data's own scale — CP10K, say, against a ~200 transcript
        # median — multiplies every cell by a factor inversely proportional to its depth, and
        # because log1p follows, that factor leaks straight back into the values. Measured on
        # both test ROIs, CP10K drove PC1's correlation with depth to 0.97 (from 0.31-0.53),
        # lost one sample's immune cluster entirely, and cut every reference-correlation margin
        # by half or more (endothelial 0.457 -> 0.168). The median target keeps the multiplier
        # near 1 and log1p behaving consistently across cells.
        #
        # The cost is that "normalized" is a per-sample scale (226 and 198 on the two test
        # ROIs), so layers["lognorm"] is not comparable across samples or against an external
        # atlas. Anything needing that comparability normalizes to CP10K itself from
        # layers["counts"] — sample_summary's centroids already do.
        rsc.pp.normalize_total(adata, inplace=True)
        rsc.pp.log1p(adata)
        # Preserved before scaling overwrites X; downstream annotation reads this.
        adata.layers["lognorm"] = adata.X.copy()

    with timer("Scale"):
        rsc.pp.scale(adata, zero_center=False, max_value=SCALE_MAX_VALUE)

    with timer("PCA"):
        rsc.pp.pca(adata, n_comps=N_PCS, random_state=0)

    with timer("Neighbors"):
        rsc.pp.neighbors(adata, metric=NEIGHBORS_METRIC, random_state=0)

    with timer("UMAP"):
        rsc.tl.umap(adata, random_state=0)

    # Sweep resolutions rather than committing to one: the neighbour graph is
    # already built, so each extra resolution only re-runs the (comparatively
    # cheap) community detection on it. Downstream picks a column by name.
    for res in resolutions:
        with timer(f"Leiden res={res:g}"):
            rsc.tl.leiden(
                adata,
                resolution=res,
                key_added=f"leiden_res_{res:.2f}",
                random_state=0,
            )

    # Move back to CPU for zarr I/O — rapids-singlecell keeps arrays on GPU.
    with timer("Move to CPU"):
        rsc.get.anndata_to_CPU(adata)

    with timer("Write zarr"):
        # Only the table changed, so only the table is written; the images are
        # hardlinked from the input store. See sdata_io.write_table_only.
        sdata.tables[table_key] = adata
        write_table_only(sdata, args.path, output_path, table_key)

    print(f"Written to {output_path}")

    timing_summary(path="cluster_sdata_gpu_timing.tsv")


if __name__ == "__main__":
    main()
