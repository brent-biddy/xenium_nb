#!/usr/bin/env python3
"""
downsample_sdata.py - Randomly subsample cells from a SpatialData zarr.

Reads an existing SpatialData zarr, subsamples the cell table to a target
fraction or fixed count, and writes a new self-contained SpatialData zarr.
Useful for reducing dataset size for local clustering runs.

Exactly one of --fraction or --n_cells must be provided.

Writes downsampled.zarr into the current working directory.

Usage:
    downsample_sdata.py --sample ROI1_A --path /data/ROI1_A.zarr --fraction 0.1
    downsample_sdata.py --sample ROI1_A --path /data/ROI1_A.zarr --n_cells 50000
"""

import argparse

import scanpy as sc
import spatialdata

from timer import timer, timing_summary


def parse_args():
    parser = argparse.ArgumentParser(
        description="Subsample cells from a SpatialData zarr"
    )
    parser.add_argument("--sample", required=True, help="Sample identifier")
    parser.add_argument("--path", required=True, help="Path to input SpatialData zarr")

    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument(
        "--fraction", type=float,
        help="Fraction of cells to retain (0 < fraction ≤ 1)"
    )
    group.add_argument(
        "--n_cells", type=int,
        help="Number of cells to retain"
    )

    return parser.parse_args()


def main():
    args = parse_args()

    output_path = "downsampled.zarr"

    print(f"Sample:  {args.sample}")
    print(f"Input:   {args.path}")
    print(f"Output:  {output_path}")

    with timer("Read zarr"):
        sdata = spatialdata.read_zarr(args.path)

    table_key = "table"
    with timer("Extract table"):
        adata = sdata.tables[table_key].copy()

    n_before = adata.n_obs
    print(f"Table:   {n_before:,} cells × {adata.n_vars:,} genes  (key: '{table_key}')")

    if args.fraction is not None:
        if not 0 < args.fraction <= 1:
            raise ValueError(f"--fraction must be between 0 and 1, got {args.fraction}")
        n_target = int(n_before * args.fraction)
    else:
        n_target = args.n_cells

    if n_target >= n_before:
        print(f"Requested {n_target:,} cells ≥ available {n_before:,} — skipping subsample.")
    else:
        with timer("Subsample"):
            sc.pp.subsample(adata, n_obs=n_target, random_state=0, copy=False)
        print(f"Subsampled {n_before:,} → {adata.n_obs:,} cells.")

    with timer("Write zarr"):
        sdata.tables[table_key] = adata
        sdata.write(output_path)

    print(f"Written to {output_path}")

    timing_summary(path="downsample_sdata_timing.tsv")


if __name__ == "__main__":
    main()
