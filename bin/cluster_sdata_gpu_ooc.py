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
"""

import argparse

import anndata as ad
import cupy as cp
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
        default=2000,
        help="Number of highly variable genes to keep (default: 2000)",
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
    X = read_dask(store["X"], (chunk_size, var.shape[0]))
    return ad.AnnData(X=X, obs=obs, var=var, uns=uns)


def main():
    args = parse_args()

    output_path = "clustered.zarr"

    print(f"Sample:     {args.sample}")
    print(f"Input:      {args.path}")
    print(f"Output:     {output_path}")
    print(f"Chunk size: {args.chunk_size:,} rows")

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
        cell_mask = adata.obs["total_counts"].to_numpy() >= 10
        adata = adata[cell_mask].copy()
        gene_mask = adata.var["n_cells_by_counts"].to_numpy() >= 5
        adata = adata[:, gene_mask].copy()

    print(f"Filtered {n_before - adata.n_obs:,} low-quality cells.")
    print(f"Retained {adata.n_obs:,} cells x {adata.n_vars:,} genes.")

    with timer("Normalize"):
        rsc.pp.normalize_total(adata, inplace=True)
        rsc.pp.log1p(adata)

    with timer("HVG"):
        rsc.pp.highly_variable_genes(adata, n_top_genes=args.n_top_genes)
        rsc.pp.filter_highly_variable(adata)

    with timer("PCA"):
        rsc.pp.pca(adata, random_state=0)
        # pca() only auto-syncs the covariance/mean computation; the resulting
        # embedding is still a lazy Dask array and needs an explicit compute().
        adata.obsm["X_pca"] = adata.obsm["X_pca"].compute()

    # From here on the working set is the PCA embedding (n_obs x n_comps), not
    # the full gene matrix — small enough to finish in-memory on GPU. Materialize
    # the (now HVG-subset) X too, since it still needs to be written to the zarr;
    # it's already GPU-resident (cupy chunks), just still wrapped in a Dask array.
    with timer("Materialize"):
        adata.X = adata.X.compute()

    with timer("Neighbors"):
        rsc.pp.neighbors(adata, random_state=0)

    with timer("UMAP"):
        rsc.tl.umap(adata, random_state=0)

    with timer("Leiden"):
        rsc.tl.leiden(adata, random_state=0)

    print(f"Leiden clustering: {adata.obs['leiden'].nunique()} clusters")

    with timer("Move to CPU"):
        rsc.get.anndata_to_CPU(adata)

    with timer("Write zarr"):
        sdata.tables[args.table_key] = adata
        sdata.write(output_path)

    print(f"Written to {output_path}")

    timing_summary(path="cluster_sdata_gpu_ooc_timing.tsv")


if __name__ == "__main__":
    main()
