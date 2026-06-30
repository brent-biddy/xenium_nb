#!/usr/bin/env python3
"""
cluster_sdata_gpu.py - GPU-accelerated QC, normalize, and cluster a SpatialData zarr.

Mirrors cluster_sdata.py but uses rapids-singlecell for the compute-heavy steps
(QC, normalization, PCA, neighbors, UMAP, Leiden). After clustering, data is moved
back to CPU for squidpy spatial steps and zarr I/O.

Requires an Apptainer/Docker image with rapids-singlecell and a CUDA-capable GPU.

Writes clustered.zarr into the current working directory.

Usage:
    cluster_sdata_gpu.py --sample ROI1_A --path /data/ROI1_A.zarr
"""

import argparse

import rapids_singlecell as rsc
import scanpy as sc
import squidpy as sq
import spatialdata

from timer import timer, timing_summary


def parse_args():
    parser = argparse.ArgumentParser(
        description="GPU-accelerated clustering of a SpatialData zarr"
    )
    parser.add_argument("--sample", required=True, help="Sample identifier")
    parser.add_argument("--path", required=True, help="Path to input SpatialData zarr")
    return parser.parse_args()


def main():
    args = parse_args()

    output_path = "clustered.zarr"

    print(f"Sample:  {args.sample}")
    print(f"Input:   {args.path}")
    print(f"Output:  {output_path}")

    with timer("Read zarr"):
        sdata = spatialdata.read_zarr(args.path)

    table_key = "table"
    with timer("Extract table"):
        adata = sdata.tables[table_key].copy()

    print(f"Table:   {adata.n_obs:,} cells × {adata.n_vars:,} genes  (key: '{table_key}')")

    with timer("Move to GPU"):
        rsc.get.anndata_to_GPU(adata)

    with timer("QC"):
        rsc.pp.calculate_qc_metrics(adata, percent_top=(10, 20, 50, 150), inplace=True)
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
        rsc.pp.pca(adata)

    with timer("Neighbors"):
        rsc.pp.neighbors(adata)

    with timer("UMAP"):
        rsc.tl.umap(adata)

    with timer("Leiden"):
        rsc.tl.leiden(adata)

    print(f"Leiden clustering: {adata.obs['leiden'].nunique()} clusters")

    # Move back to CPU before squidpy spatial steps and zarr I/O —
    # sq.gr has no GPU-accelerated equivalents.
    with timer("Move to CPU"):
        rsc.get.anndata_to_CPU(adata)

    with timer("Spatial neighbors"):
        sq.gr.spatial_neighbors(adata, coord_type="generic", delaunay=True)

    with timer("Nhood enrichment"):
        sq.gr.nhood_enrichment(adata, cluster_key="leiden")

    with timer("Write zarr"):
        sdata.tables[table_key] = adata
        sdata.write(output_path)

    print(f"Written to {output_path}")

    timing_summary(path="cluster_sdata_gpu_timing.tsv")


if __name__ == "__main__":
    main()
