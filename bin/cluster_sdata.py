#!/usr/bin/env python3
"""
cluster_sdata.py - QC, normalize, and cluster a SpatialData zarr.

Reads an existing SpatialData zarr, filters, selects highly variable genes,
normalizes, scales, runs PCA, neighbours, UMAP, and a sweep of Leiden
clusterings, then writes a new self-contained SpatialData zarr whose table
carries the embeddings, neighbour graph, and one cluster-label column per
resolution alongside the original spatial elements.

The recipe follows 10x's Xenium Prime workshop vignette, which targets the same
5K panel this pipeline is used with. Only the three filter thresholds are
adjustable; the embedding steps are constants, because they define what a
cluster means here and all three cluster_sdata* steps must agree on them.

The written table carries three expression forms: X (scaled, what the embedding
was built from), layers["counts"] (raw), and layers["lognorm"] (log-normalized).
Downstream marker and annotation work wants lognorm — after scaling, X no longer
holds expression values.

Writes clustered.zarr into the current working directory.

Usage:
    cluster_sdata.py --sample ROI1_A --path /data/ROI1_A.zarr
    cluster_sdata.py --sample ROI1_A --path /data/ROI1_A.zarr --resolutions 0.5 1.0
"""

import argparse

import scanpy as sc
import spatialdata

from timer import timer, timing_summary

DEFAULT_RESOLUTIONS = [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0]

# The filtering cut. These live here rather than in nextflow.config for the same reason
# DEFAULT_RESOLUTIONS does: the modules omit the flag when the param is unset, so the
# default is defined in exactly one place and the three cluster_sdata* steps cannot
# drift apart on what "unfiltered" means.
#
# The values are 10x's, from the Xenium Prime workshop vignette, which is built on the
# same 5K panel this pipeline targets. They were adopted on evidence rather than
# authority: on both downsampled test ROIs, clusters from these filters annotate
# measurably better against a reference than the previous 10/5 cut did — mean
# correlation margin 0.140 -> 0.263 on one ROI, with all four reference cell types
# resolved instead of three. See HVG_N_TOP_GENES for the embedding half of the same
# comparison.
DEFAULT_MIN_COUNTS = 20
DEFAULT_MIN_CELLS = 100

# Upper bound as a quantile of transcript_counts, not an absolute number: cells vary
# ten-fold in depth between samples, so a fixed ceiling would be a different cut on
# every one. Targets doublets and segmentation merges — two cells called as one — which
# sit between populations and blur the boundary the clustering is trying to find.
# Set to 0 to disable.
DEFAULT_MAX_COUNTS_QUANTILE = 0.98

# The embedding recipe, from the same 10x Xenium Prime vignette. Constants rather than
# CLI flags on purpose: these define what "a cluster" means in this pipeline, and the
# three cluster_sdata* steps have to agree on them or their labels stop being
# comparable. Changing one is a change to the pipeline, made here, not a per-run knob.
#
# Measured against the previous recipe (no HVG, no scaling, 50 PCs, euclidean) on both
# downsampled test ROIs, holding filtering constant: on one ROI the old recipe produced
# 3 clusters, one confident reference cell type, and marker-score separations of 0.010
# — it did not isolate oocytes or granulosa at all — while this one produced 7 clusters,
# 3 types, and separations 5-6x higher. Oocyte separation on the other ROI more than
# doubled, 0.079 -> 0.176.
HVG_N_TOP_GENES = 2000
SCALE_MAX_VALUE = 10
N_PCS = 30
NEIGHBORS_METRIC = "cosine"
# The vignette uses -1, scanpy's "run to convergence" sentinel, but rapids-singlecell
# cannot: it forwards the value to cuGraph as a size_t and raises OverflowError on a
# negative. 100 is rapids' own default, is far above the handful of passes Leiden
# actually needs to converge here, and is accepted by both backends — which matters
# more than matching the sentinel, since the CPU and GPU steps have to agree.
LEIDEN_N_ITERATIONS = 100


def parse_args():
    parser = argparse.ArgumentParser(
        description="Cluster a SpatialData zarr and write results back"
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

    with timer("Filter"):
        # No calculate_qc_metrics call here. create_sdata already annotated every
        # per-cell and per-gene metric, at the one point X was in memory; recomputing
        # them produced nothing new and actively did harm.
        #
        # With the default expr_type="counts" the call overwrote the Xenium-native
        # obs["total_counts"] — transcripts plus the five control/codeword counters,
        # which the reader takes from cells.parquet — with a plain X.sum(axis=1). The
        # column therefore meant one thing on create_sdata's zarrs and another on
        # these, silently. Running before filter_genes, it did not even leave a correct
        # row sum: the value described the full panel while the object shipped only the
        # genes that survived (4,753 of 5,101 on one ROI, wrong for 802 cells).
        #
        # It also carried percent_top=(10, 20, 50, 150), the setting create_sdata moved
        # off for good reason — at a median 174 genes per cell, top_150 reads ~100% for
        # a large share of cells and inverts the metric.
        #
        # Dropping it is safe: filter_cells and filter_genes derive their thresholds
        # from X directly and never read obs, so the cut is unchanged.
        n_before = adata.n_obs
        # The quantile is taken on the UNFILTERED table, before min_counts removes the
        # low tail — matching the vignette, and the more reproducible definition: a cut
        # derived from the surviving cells would move when the lower threshold moved.
        # Read from obs["transcript_counts"] rather than a row sum: the two are equal
        # here, but transcript_counts is the column create_sdata guarantees excludes the
        # control/codeword counters.
        max_counts = (
            float(adata.obs["transcript_counts"].quantile(args.max_counts_quantile))
            if args.max_counts_quantile else None
        )
        sc.pp.filter_cells(adata, min_counts=args.min_counts)
        if max_counts is not None:
            sc.pp.filter_cells(adata, max_counts=max_counts)
            print(f"Upper cut at q{args.max_counts_quantile:g} = {max_counts:,.0f} transcripts.")
        sc.pp.filter_genes(adata, min_cells=args.min_cells)

    print(f"Filtered {n_before - adata.n_obs:,} low-quality cells.")
    print(f"Retained {adata.n_obs:,} cells × {adata.n_vars:,} genes.")

    with timer("HVG"):
        # seurat_v3 runs on RAW counts, before normalization — that is what the flavor
        # expects, and it is why this sits above the normalize block rather than after
        # it. Genes are flagged, not subset: scanpy's pca reads var["highly_variable"]
        # and restricts itself, while every other gene stays available to downstream
        # marker and annotation work.
        sc.pp.highly_variable_genes(
            adata, flavor="seurat_v3", n_top_genes=HVG_N_TOP_GENES)

    with timer("Normalize"):
        adata.layers["counts"] = adata.X.copy()
        sc.pp.normalize_total(adata, inplace=True)
        sc.pp.log1p(adata)
        # Keep the log-normalized values before scaling overwrites X. Everything
        # downstream that reads expression rather than a reduced representation —
        # rank_genes_groups, marker scoring, the reference centroid correlation — needs
        # this layer, because after scale() X holds z-like values, not expression.
        adata.layers["lognorm"] = adata.X.copy()

    with timer("Scale"):
        # zero_center=False keeps X sparse, which matters at whole-slide cell counts;
        # PCA centers internally, so the centering is not lost. max_value clips the
        # tail so a handful of extreme cells cannot dominate a component.
        sc.pp.scale(adata, zero_center=False, max_value=SCALE_MAX_VALUE)

    with timer("PCA"):
        sc.pp.pca(adata, n_comps=N_PCS, random_state=0)

    with timer("Neighbors"):
        # Cosine rather than the euclidean default: it compares the shape of an
        # expression profile and ignores its magnitude, so cells of one type separate
        # from another type rather than from shallower cells of the same type.
        sc.pp.neighbors(adata, metric=NEIGHBORS_METRIC, random_state=0)

    with timer("UMAP"):
        sc.tl.umap(adata, random_state=0)

    # Sweep resolutions rather than committing to one: the neighbour graph is
    # already built, so each extra resolution only re-runs the (comparatively
    # cheap) community detection on it. Downstream picks a column by name.
    for res in resolutions:
        with timer(f"Leiden res={res:g}"):
            # Use the igraph backend (orders of magnitude faster than the legacy
            # leidenalg default). flavor="igraph" requires an undirected graph.
            sc.tl.leiden(
                adata,
                resolution=res,
                key_added=f"leiden_res_{res:.2f}",
                flavor="igraph",
                n_iterations=LEIDEN_N_ITERATIONS,
                directed=False,
                random_state=0,
            )

    with timer("Write zarr"):
        sdata.tables[table_key] = adata
        sdata.write(output_path)

    print(f"Written to {output_path}")

    timing_summary(path="cluster_sdata_timing.tsv")


if __name__ == "__main__":
    main()
