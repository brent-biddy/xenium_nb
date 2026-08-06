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

from timer import timer, timing_summary

DEFAULT_RESOLUTIONS = [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0]


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
    return parser.parse_args()


def main():
    args = parse_args()

    resolutions = sorted(args.resolutions)

    output_path = "clustered.zarr"

    print(f"Sample:  {args.sample}")
    print(f"Input:   {args.path}")
    print(f"Output:  {output_path}")
    print(f"Res:     {', '.join(f'{r:g}' for r in resolutions)}")

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
        rsc.pp.filter_cells(adata, min_counts=10)
        rsc.pp.filter_genes(adata, min_cells=5)

    print(f"Filtered {n_before - adata.n_obs:,} low-quality cells.")
    print(f"Retained {adata.n_obs:,} cells × {adata.n_vars:,} genes.")

    with timer("Normalize"):
        adata.layers["counts"] = adata.X.copy()
        rsc.pp.normalize_total(adata, inplace=True)
        rsc.pp.log1p(adata)

    with timer("PCA"):
        rsc.pp.pca(adata, random_state=0)

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

    # Move back to CPU for zarr I/O — rapids-singlecell keeps arrays on GPU.
    with timer("Move to CPU"):
        rsc.get.anndata_to_CPU(adata)

    with timer("Write zarr"):
        sdata.tables[table_key] = adata
        sdata.write(output_path)

    print(f"Written to {output_path}")

    timing_summary(path="cluster_sdata_gpu_timing.tsv")


if __name__ == "__main__":
    main()
