#!/usr/bin/env python3
"""
cluster_sdata_gpu_ooc.py - Out-of-core GPU QC, normalize, and cluster a SpatialData zarr.

Mirrors cluster_sdata_gpu.py but streams the table's X matrix through Dask
instead of loading it into GPU memory whole, so tables too large for VRAM
(e.g. a merged cohort from concat_sdata) can still be processed on a single
GPU. Only QC/normalization/HVG/PCA run against the full lazy matrix —
rapids-singlecell has no Dask-native neighbors/UMAP/Leiden, so those run on
the already-computed, much smaller PCA embedding once PCA has reduced the
working set.

Requires an Apptainer/Docker image with rapids-singlecell, dask, and zarr,
plus a CUDA-capable GPU.

Writes clustered.zarr into the current working directory.

Usage:
    cluster_sdata_gpu_ooc.py --sample cohort --path /data/cohort.zarr --chunk-size 20000
    cluster_sdata_gpu_ooc.py --sample cohort --path /data/cohort.zarr --resolutions 0.5 1.0
"""

import argparse

import anndata as ad
import cupy as cp
import numpy as np
import rapids_singlecell as rsc
import rmm
import spatialdata
import zarr
from rmm.allocators.cupy import rmm_cupy_allocator

from timer import timer, timing_summary

try:
    from anndata.experimental import read_elem_lazy as read_dask
except ImportError:  # older anndata: same functionality under the old name
    from anndata.experimental import read_elem_as_dask as read_dask

DEFAULT_RESOLUTIONS = [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0]

# Must stay in step with cluster_sdata.py and cluster_sdata_gpu.py. Note this step
# applies the cut as a boolean mask over obs/var rather than through filter_cells /
# filter_genes (see the QC block), so the flags mean the same thing but reach the data
# by a different route.
DEFAULT_MIN_COUNTS = 20
DEFAULT_MIN_CELLS = 100
DEFAULT_MAX_COUNTS_QUANTILE = 0.98

# The embedding recipe — see cluster_sdata.py for what these are and why.
HVG_N_TOP_GENES = 2000
SCALE_MAX_VALUE = 10
N_PCS = 30
NEIGHBORS_METRIC = "cosine"
LEIDEN_N_ITERATIONS = 100


def parse_args():
    parser = argparse.ArgumentParser(
        description="Out-of-core GPU clustering of a SpatialData zarr via Dask"
    )
    parser.add_argument("--sample", required=True, help="Sample identifier")
    parser.add_argument("--path", required=True, help="Path to input SpatialData zarr")
    parser.add_argument(
        "--table-key", default="table", help="Table key within the zarr (default: table)"
    )
    parser.add_argument(
        "--chunk-size",
        type=int,
        default=20_000,
        help="Row chunk size for the lazy X array (default: 20000)",
    )
    parser.add_argument(
        "--n-top-genes",
        type=int,
        default=None,
        help=f"Hard-subset to this many highly variable genes, replacing the "
        f"{HVG_N_TOP_GENES} all three steps otherwise select and flag. Off by default, "
        "so results match cluster_sdata/cluster_sdata_gpu. Set it only when the "
        "materialized X (n_obs x n_vars) will not fit — it shrinks that matrix at the "
        "cost of clustering on a different feature space and dropping the other genes "
        "from the output entirely.",
    )
    parser.add_argument(
        "--resolutions",
        type=float,
        nargs="+",
        default=DEFAULT_RESOLUTIONS,
        metavar="RES",
        help="Leiden resolutions to sweep; one obs column is written per value "
        f"(default: {' '.join(str(r) for r in DEFAULT_RESOLUTIONS)}).",
    )
    # Underscored, unlike this file's --table-key / --chunk-size / --n-top-genes, to
    # match cluster_sdata.py and cluster_sdata_gpu.py: the three steps take the same
    # cut and one module pattern passes it, so the flag spelling has to agree.
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
        help="Drop cells above this quantile of the count column, to remove doublets "
        f"and segmentation merges. 0 disables (default: {DEFAULT_MAX_COUNTS_QUANTILE}).",
    )
    return parser.parse_args()


def read_table_lazy(path, table_key, chunk_size):
    """Build an AnnData with X as a lazy Dask array from a SpatialData table group.

    obs/var/uns are read eagerly — they're small relative to X and rapids-singlecell's
    QC/filtering steps need them as plain pandas objects, not Dask-backed.
    """
    store = zarr.open(f"{path}/tables/{table_key}", mode="r")
    obs = ad.io.read_elem(store["obs"])
    var = ad.io.read_elem(store["var"])
    uns = ad.io.read_elem(store["uns"]) if "uns" in store else {}
    # ad.io.read_elem returns spatialdata_attrs.region as an ndarray, unlike
    # spatialdata's own reader — TableModel.validate() requires a list.
    region = uns.get("spatialdata_attrs", {}).get("region")
    if isinstance(region, np.ndarray):
        uns["spatialdata_attrs"]["region"] = region.tolist()
    X = read_dask(store["X"], (chunk_size, var.shape[0]))
    return ad.AnnData(X=X, obs=obs, var=var, uns=uns)


def main():
    args = parse_args()

    resolutions = sorted(args.resolutions)

    output_path = "clustered.zarr"

    print(f"Sample:     {args.sample}")
    print(f"Input:      {args.path}")
    print(f"Output:     {output_path}")
    print(f"Chunk size: {args.chunk_size:,} rows")
    print(f"Res:        {', '.join(f'{r:g}' for r in resolutions)}")
    print(f"Filter:     min_counts={args.min_counts}, min_cells={args.min_cells}, "
          f"max_counts_quantile={args.max_counts_quantile:g}")

    # Managed memory lets chunks spill to host RAM instead of OOM-ing when the
    # dataset (or an intermediate) doesn't fit in VRAM — the whole point of an
    # out-of-core run. Trades some throughput for headroom vs. the plain pool
    # allocator cluster_sdata_gpu.py relies on implicitly.
    rmm.reinitialize(managed_memory=True, pool_allocator=False)
    cp.cuda.set_allocator(rmm_cupy_allocator)

    with timer("Read spatial elements"):
        # selection excludes tables: the table's X is streamed lazily below
        # instead of materialized whole by spatialdata's default AnnData reader.
        sdata = spatialdata.read_zarr(
            args.path, selection=("images", "labels", "points", "shapes")
        )

    with timer("Read table (lazy)"):
        adata = read_table_lazy(args.path, args.table_key, args.chunk_size)

    print(f"Table:      {adata.n_obs:,} cells x {adata.n_vars:,} genes  (key: '{args.table_key}')")

    with timer("Move to GPU (lazy)"):
        # anndata_to_GPU understands Dask arrays: it map_blocks-converts each
        # chunk's meta from scipy/numpy to cupy without forcing a compute.
        # rapids-singlecell's ops require GPU-backed input (even when Dask-lazy),
        # so this has to happen before the first pp call, not after PCA.
        rsc.get.anndata_to_GPU(adata)

    with timer("QC"):
        # Lazy — does not force computation of the underlying Dask array.
        rsc.pp.calculate_qc_metrics(adata)
        n_before = adata.n_obs
        # Boolean-index + .copy() rather than filter_cells/filter_genes: filtering
        # through views is incompatible with — and much slower on — Dask-backed X.
        #
        # Both bounds read obs["total_counts"], the column rsc.pp.calculate_qc_metrics
        # just overwrote with a plain row sum. That is deliberate and is why the
        # recompute above is kept: on a merged cohort from concat_sdata the row sum is
        # the quantity these thresholds are defined against. The other two steps reach
        # the same numbers through transcript_counts; see CLAUDE.md.
        counts = adata.obs["total_counts"].to_numpy()
        cell_mask = counts >= args.min_counts
        if args.max_counts_quantile:
            # Quantile over the unfiltered column, matching the other two steps.
            max_counts = float(np.quantile(counts, args.max_counts_quantile))
            cell_mask &= counts <= max_counts
            print(f"Upper cut at q{args.max_counts_quantile:g} = {max_counts:,.0f} counts.")
        adata = adata[cell_mask].copy()
        # Recompute before the gene cut. n_cells_by_counts from the pass above counted
        # detections in cells the mask has since dropped, so genes would be judged on a
        # population this run does not keep — cluster_sdata and cluster_sdata_gpu run
        # filter_genes after filter_cells and so count only survivors. Left stale, this
        # step retained 3,328 genes where the other two retained 3,187 on one test ROI.
        rsc.pp.calculate_qc_metrics(adata)
        gene_mask = adata.var["n_cells_by_counts"].to_numpy() >= args.min_cells
        adata = adata[:, gene_mask].copy()

    print(f"Filtered {n_before - adata.n_obs:,} low-quality cells.")
    print(f"Retained {adata.n_obs:,} cells x {adata.n_vars:,} genes.")

    with timer("HVG"):
        # Now mandatory and no longer an escape hatch: all three cluster_sdata* steps
        # select the same 2000 genes, so their labels stay comparable. --n_top_genes
        # remains only to narrow further when a merged cohort will not fit on the GPU.
        #
        # Genes are FLAGGED, not subset, matching the other two steps — pca reads
        # var["highly_variable"] and restricts itself, while lognorm below keeps every
        # surviving gene available to downstream annotation. filter_highly_variable is
        # only called when --n_top_genes forces a hard narrowing for memory.
        n_top = args.n_top_genes or HVG_N_TOP_GENES
        rsc.pp.highly_variable_genes(adata, flavor="seurat_v3", n_top_genes=n_top)
        if args.n_top_genes:
            rsc.pp.filter_highly_variable(adata)
            print(f"HVG subset: {adata.n_vars:,} genes retained.")

    with timer("Normalize"):
        rsc.pp.normalize_total(adata, inplace=True)
        rsc.pp.log1p(adata)
        # Preserved before scaling overwrites X; downstream annotation reads this.
        adata.layers["lognorm"] = adata.X.copy()

    with timer("Scale"):
        rsc.pp.scale(adata, zero_center=False, max_value=SCALE_MAX_VALUE)

    with timer("PCA"):
        rsc.pp.pca(adata, n_comps=N_PCS, random_state=0)
        # pca() only auto-syncs the covariance/mean computation; the resulting
        # embedding is still a lazy Dask array and needs an explicit compute().
        adata.obsm["X_pca"] = adata.obsm["X_pca"].compute()

    # From here on the working set is the PCA embedding (n_obs x n_comps), not
    # the full gene matrix — small enough to finish in-memory on GPU. X still has to
    # be materialized for the zarr write; it's already GPU-resident (cupy chunks),
    # just still wrapped in a Dask array. This is the step --n-top-genes shrinks, and
    # the reason managed_memory is on above: an oversized X spills to host rather
    # than OOM-ing.
    with timer("Materialize"):
        adata.X = adata.X.compute()

    with timer("Neighbors"):
        rsc.pp.neighbors(adata, random_state=0)

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

    with timer("Move to CPU"):
        rsc.get.anndata_to_CPU(adata)

    with timer("Write zarr"):
        sdata.tables[args.table_key] = adata
        sdata.write(output_path)

    print(f"Written to {output_path}")

    timing_summary(path="cluster_sdata_gpu_ooc_timing.tsv")


if __name__ == "__main__":
    main()
