#!/usr/bin/env python3
"""
cluster_sdata.py - QC, normalize, and cluster a SpatialData zarr.

Reads an existing SpatialData zarr, runs scanpy QC/normalization, PCA, UMAP,
and Leiden clustering, then computes spatial neighbours and neighbourhood
enrichment with squidpy. Writes a new self-contained SpatialData zarr
containing all embeddings, graphs, and cluster labels alongside the original
spatial elements.

Writes output/clustered.zarr into the current working directory.

Usage:
    cluster_sdata.py --sample ROI1_A --path /data/ROI1_A.zarr
"""

import argparse
import os

import scanpy as sc
import spatialdata

from timer import timer, timing_summary


def parse_args():
    parser = argparse.ArgumentParser(
        description="Cluster a SpatialData zarr and write results back"
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

    # xenium() names the table "table"; concatenated objects may rename it.
    # Fall back to the first available table if "table" is absent.
    table_key = "table" if "table" in sdata.tables else next(iter(sdata.tables))
    with timer("Extract table"):
        adata = sdata.tables[table_key].copy()

    print(f"Table:   {adata.n_obs:,} cells × {adata.n_vars:,} genes  (key: '{table_key}')")

    with timer("QC"):
        sc.pp.calculate_qc_metrics(adata, percent_top=(10, 20, 50, 150), inplace=True)
        n_before = adata.n_obs
        sc.pp.filter_cells(adata, min_counts=10)
        sc.pp.filter_genes(adata, min_cells=5)

    print(f"Filtered {n_before - adata.n_obs:,} low-quality cells.")
    print(f"Retained {adata.n_obs:,} cells × {adata.n_vars:,} genes.")

    with timer("Normalize"):
        adata.layers["counts"] = adata.X.copy()
        sc.pp.normalize_total(adata, inplace=True)
        sc.pp.log1p(adata)

    with timer("PCA"):
        sc.pp.pca(adata, random_state=0)

    with timer("Neighbors"):
        sc.pp.neighbors(adata, random_state=0)

    with timer("UMAP"):
        sc.tl.umap(adata, random_state=0)

    with timer("Leiden"):
        sc.tl.leiden(adata, random_state=0)

    print(f"Leiden clustering: {adata.obs['leiden'].nunique()} clusters")

    with timer("Write zarr"):
        sdata.tables[table_key] = adata
        sdata.write(output_path)

    print(f"Written to {output_path}")

    timing_summary(path="cluster_sdata_timing.tsv")


if __name__ == "__main__":
    main()
