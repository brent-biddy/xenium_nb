#!/usr/bin/env python3
"""
downsample_sdata.py - Randomly subsample cells from a SpatialData zarr.

Reads an existing SpatialData zarr, subsamples the cell table to a target
fraction or fixed count, and writes a new self-contained SpatialData zarr.
Useful for reducing dataset size for local clustering runs.

Subsampling is cluster-aware when the table carries cluster labels. Each
cluster keeps its proportional share of the target, but never fewer than
--min_cells_per_cluster cells, so rare populations survive a deep downsample
instead of being rounded away. Inputs with no label column (e.g. a raw
create_sdata zarr) fall back to uniform subsampling.

Exactly one of --fraction or --n_cells must be provided.

Writes downsampled.zarr into the current working directory.

Usage:
    downsample_sdata.py --sample ROI1_A --path /data/ROI1_A.zarr --fraction 0.1
    downsample_sdata.py --sample ROI1_A --path /data/ROI1_A.zarr --n_cells 50000
    downsample_sdata.py --sample ROI1_A --path /data/ROI1_A.zarr --fraction 0.1 \
        --stratify_by cell_type --min_cells_per_cluster 100
"""

import argparse

import numpy as np
import spatialdata

from timer import timer, timing_summary

# obs column looked for when --stratify_by is not given. Matches the key
# cluster_sdata / cluster_sdata_gpu / cluster_sdata_gpu_ooc write (sc.tl.leiden
# with no key_added), so a clustered input stratifies with no extra flags and a
# raw one falls back to uniform.
DEFAULT_STRATIFY_COL = "leiden"


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

    parser.add_argument(
        "--stratify_by", default=None,
        help=f"obs column to stratify on. Default: use '{DEFAULT_STRATIFY_COL}' if "
             f"present, else subsample uniformly. An explicitly named column must exist."
    )
    parser.add_argument(
        "--min_cells_per_cluster", type=int, default=50,
        help="Floor on cells retained per cluster when stratifying (default: 50). "
             "Clusters smaller than the floor are kept whole."
    )
    parser.add_argument(
        "--seed", type=int, default=0,
        help="Random seed for subsampling (default: 0)"
    )

    return parser.parse_args()


def resolve_stratify_col(adata, requested):
    """Return the obs column to stratify on, or None to subsample uniformly.

    An explicitly requested column is required to exist — a typo should fail the
    run, not silently degrade to uniform. The default column is only a hint, so
    its absence just means the input was never clustered.
    """
    if requested is None:
        if DEFAULT_STRATIFY_COL in adata.obs.columns:
            return DEFAULT_STRATIFY_COL
        print(
            f"No '{DEFAULT_STRATIFY_COL}' column in table.obs — "
            f"subsampling uniformly."
        )
        return None

    if requested not in adata.obs.columns:
        raise KeyError(
            f"--stratify_by '{requested}' not found in table.obs. "
            f"Available columns: {list(adata.obs.columns)}"
        )
    return requested


def allocate_stratified(sizes, n_target, floor):
    """Split n_target cells across clusters proportionally, with a per-cluster floor.

    Each cluster's share of n_target is its share of the total, raised to `floor`
    if proportionality would drop it below that, and capped at the cells it
    actually has. Those three constraints rarely sum to exactly n_target, so the
    remainder is settled by trimming the largest allocations (or topping them up),
    which keeps the rounding error off the rare clusters the floor protects.

    Returns a dict of cluster -> number of cells to keep.
    """
    keys = list(sizes)
    n_total = sum(sizes.values())

    # A cluster can never give up more cells than it has, so the floor is capped
    # by the cluster's own size.
    floors = {k: min(floor, sizes[k]) for k in keys}
    alloc = {
        k: min(sizes[k], max(floors[k], round(n_target * sizes[k] / n_total)))
        for k in keys
    }

    diff = sum(alloc.values()) - n_target
    if diff > 0:
        # One pass over descending allocations clears the excess unless the floors
        # alone already exceed n_target, which is reported below.
        for k in sorted(keys, key=lambda k: alloc[k], reverse=True):
            if diff == 0:
                break
            take = min(diff, alloc[k] - floors[k])
            alloc[k] -= take
            diff -= take
        if diff > 0:
            print(
                f"WARNING: --min_cells_per_cluster {floor} across {len(keys)} clusters "
                f"needs {sum(floors.values()):,} cells, more than the {n_target:,} "
                f"requested. Honoring the floor; output will exceed the target."
            )
    elif diff < 0:
        # Headroom always covers the shortfall here, since n_target <= n_total.
        for k in sorted(keys, key=lambda k: sizes[k], reverse=True):
            if diff == 0:
                break
            give = min(-diff, sizes[k] - alloc[k])
            alloc[k] += give
            diff += give

    return alloc


def select_indices(adata, col, n_target, floor, rng):
    """Return sorted positional indices of the cells to keep."""
    if col is None:
        return np.sort(rng.choice(adata.n_obs, size=n_target, replace=False))

    counts = adata.obs[col].value_counts()
    # Categorical obs columns report unused categories at zero; drop them so they
    # do not show up as empty clusters in the allocation and the summary.
    sizes = {k: int(v) for k, v in counts.items() if v > 0}
    print(f"Stratifying on '{col}': {len(sizes)} clusters")

    alloc = allocate_stratified(sizes, n_target, floor)

    picks = []
    for k in sorted(sizes, key=lambda k: (-sizes[k], str(k))):
        positions = np.flatnonzero((adata.obs[col] == k).to_numpy())
        picks.append(rng.choice(positions, size=alloc[k], replace=False))
        print(f"  {str(k):>12}: {sizes[k]:>8,} → {alloc[k]:>7,}")

    # Sorting keeps the output in the input's cell order, so the table's rows stay
    # aligned with the instance_key order the other sdata elements join on.
    return np.sort(np.concatenate(picks))


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
        col = resolve_stratify_col(adata, args.stratify_by)
        rng = np.random.default_rng(args.seed)
        with timer("Subsample"):
            idx = select_indices(adata, col, n_target, args.min_cells_per_cluster, rng)
            adata = adata[idx].copy()
        print(f"Subsampled {n_before:,} → {adata.n_obs:,} cells.")

    with timer("Write zarr"):
        sdata.tables[table_key] = adata
        sdata.write(output_path)

    print(f"Written to {output_path}")

    timing_summary(path="downsample_sdata_timing.tsv")


if __name__ == "__main__":
    main()
