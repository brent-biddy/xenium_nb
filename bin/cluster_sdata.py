#!/usr/bin/env python3
"""
cluster_sdata.py - QC, normalize, and cluster a SpatialData zarr.

Reads an existing SpatialData zarr, runs scanpy QC/normalization, PCA,
neighbours, UMAP, and a sweep of Leiden clusterings, then writes a new
self-contained SpatialData zarr whose table carries the embeddings, neighbour
graph, and one cluster-label column per resolution alongside the original
spatial elements.

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

    # Sweep resolutions rather than committing to one: the neighbour graph is
    # already built, so each extra resolution only re-runs the (comparatively
    # cheap) community detection on it. Downstream picks a column by name.
    for res in resolutions:
        with timer(f"Leiden res={res:g}"):
            # Use the igraph backend (orders of magnitude faster than the legacy
            # leidenalg default). flavor="igraph" requires an undirected graph, and
            # n_iterations=2 reproduces the old default's convergence behavior.
            sc.tl.leiden(
                adata,
                resolution=res,
                key_added=f"leiden_res_{res:.2f}",
                flavor="igraph",
                n_iterations=2,
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
