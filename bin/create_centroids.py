#!/usr/bin/env python3
"""
create_centroids.py - Reduce a clustered SpatialData zarr to per-group centroids.

Sums both expression layers over an obs column, writing one row per group of cells
instead of one per cell:

    X                 each cell's CP10K profile, summed over the group
    layers["counts"]  raw counts, summed over the group

obs describes each row with `grouping`, `group`, `n_cells` and `sample`, plus
`resolution` when the column grouped on is one of the sweep's leiden_res_<r> columns.
Requires layers["counts"], so not cluster_sdata_gpu_ooc output.

Groups by every leiden_res_<r> column in the sweep by default, or by one named obs column
with --group_by. Writes <sample>_centroids.h5ad into the current working directory (or
<sample>_centroids_<column>.h5ad for a --group_by run), alongside timing and session info
files.

See the create_centroids entry in CLAUDE.md for the design rationale.

Usage:
    create_centroids.py --sample ROI1_A --path clustered.zarr
    create_centroids.py --sample ROI1_A --path clustered.zarr --group_by cell_type
"""

import argparse

import anndata as ad
import pandas as pd
import scanpy as sc
import session_info
import spatialdata

from timer import timer, timing_summary

# The prefix cluster_sdata* writes its sweep under.
LEIDEN_PREFIX = "leiden_res_"


def parse_args():
    parser = argparse.ArgumentParser(
        description="Reduce a clustered SpatialData zarr to its per-group centroids"
    )
    parser.add_argument("--sample", required=True, help="Sample identifier")
    parser.add_argument(
        "--path", required=True,
        help="Clustered SpatialData zarr from cluster_sdata or cluster_sdata_gpu. "
             "NOT cluster_sdata_gpu_ooc, which omits layers['counts']."
    )
    parser.add_argument(
        "--group_by",
        default=None,
        help="obs column to sum over. Omit to sum over every leiden_res_<r> in the "
        "sweep, which is what the report decks read.",
    )
    return parser.parse_args()


def reduce_by(adata, column):
    """Sum both layers over one obs column, returning a groups x genes AnnData.

    Two aggregate calls, not one: sc.get.aggregate takes a single `layer`, and its
    `func` list gives several statistics over that one layer rather than one statistic
    over several. Both group on the same column, so the rows line up.
    """
    cp10k = sc.get.aggregate(adata, by=column, func="sum", layer="cp10k")
    counts = sc.get.aggregate(adata, by=column, func="sum", layer="counts")

    out = ad.AnnData(
        X=cp10k.layers["sum"],
        obs=pd.DataFrame(index=cp10k.obs_names),
        var=pd.DataFrame(index=cp10k.var_names),
        layers={"counts": counts.layers["sum"]},
    )

    out.obs["grouping"] = column
    out.obs["group"] = out.obs_names
    out.obs["n_cells"] = cp10k.obs["n_obs_aggregated"].to_numpy()

    if column.startswith(LEIDEN_PREFIX):
        out.obs["resolution"] = column.removeprefix(LEIDEN_PREFIX)

    return out


def main():
    args = parse_args()

    stem = (
        f"{args.sample}_centroids_{args.group_by}"
        if args.group_by
        else f"{args.sample}_centroids"
    )
    output_path = f"{stem}.h5ad"

    print(f"Sample:  {args.sample}")
    print(f"Input:   {args.path}")
    print(f"Output:  {output_path}")
    print(f"Group:   {args.group_by or 'leiden sweep'}")

    with timer("Read zarr"):
        sdata = spatialdata.read_zarr(args.path)

    with timer("Extract table"):
        adata = sdata.tables["table"]

    print(f"Table:   {adata.n_obs:,} cells × {adata.n_vars:,} genes")

    # Checked before the CP10K pass, so a zarr that cannot be reduced fails on its own
    # terms rather than on a bare KeyError deep inside it. cluster_sdata_gpu_ooc omits
    # this layer on purpose. There is deliberately no fallback to layers["lognorm"]:
    # it is on each sample's own median scale, which is the one thing these centroids
    # exist not to be.
    if "counts" not in adata.layers:
        raise KeyError(
            f"{args.sample}: the table has no layers['counts'], so centroids cannot be "
            f"built. Centroids are summed CP10K and need the raw matrix; the layers "
            f"present are {sorted(adata.layers)}. cluster_sdata and cluster_sdata_gpu "
            f"write this layer, but cluster_sdata_gpu_ooc omits it on purpose — if this "
            f"zarr came from that step, re-cluster the sample with cluster_sdata or "
            f"cluster_sdata_gpu if it fits."
        )

    # What to sum over. The sweep is whatever the cluster run that produced this zarr
    # chose, so it is read off the column names; sorted numerically because a sweep
    # including "1.00" does not sort lexically.
    if args.group_by is None:
        columns = sorted(
            (col for col in adata.obs if col.startswith(LEIDEN_PREFIX)),
            key=lambda col: float(col.removeprefix(LEIDEN_PREFIX)),
        )
        if not columns:
            raise ValueError(
                f"no {LEIDEN_PREFIX}<r> columns in obs — was this zarr written by "
                f"cluster_sdata or cluster_sdata_gpu? Found: {sorted(adata.obs.columns)}"
            )
    else:
        # Validated before the CP10K pass, which is the expensive part.
        if args.group_by not in adata.obs:
            raise ValueError(
                f"--group_by {args.group_by!r} is not an obs column of this zarr. This "
                f"step groups by a column some upstream step already wrote; it does not "
                f"compute labels. Found: {sorted(adata.obs.columns)}"
            )
        columns = [args.group_by]

    print(f"Reducing over {len(columns)} column(s): {', '.join(columns)}")

    with timer("CP10K"):
        # Normalized once and reused for every grouping: the per-cell scaling does not
        # depend on which column the cells are summed over. From layers["counts"], not
        # expm1(X) — scaling has overwritten X, so it no longer holds expression.
        adata.layers["cp10k"] = adata.layers["counts"].copy()
        sc.pp.normalize_total(adata, layer="cp10k", target_sum=1e4)

    blocks = []
    for column in columns:
        with timer(f"Reduce {column}"):
            blocks.append(reduce_by(adata, column))

    with timer("Assemble"):
        # index_unique=None keeps the group labels as they are; rows are made unique by
        # the (grouping, group) pair in obs, not by the index, and a suffixed index would
        # break the decks' join back onto per-cell obs. anndata still requires a unique
        # index, so the sweep's repeated labels are renumbered explicitly.
        centroids = ad.concat(blocks, axis=0, join="outer", index_unique=None)
        centroids.obs_names = [str(i) for i in range(centroids.n_obs)]
        centroids.obs["sample"] = args.sample

        # Pin the identity columns to categorical rather than letting the dtype fall out
        # of the data: anndata's writer and ad.concat convert object string columns only
        # when the cardinality makes it pay off, so the same column comes back categorical
        # from a sweep store and object from a single-resolution one.
        for col in ("grouping", "group", "resolution", "sample"):
            if col in centroids.obs:
                centroids.obs[col] = pd.Categorical(centroids.obs[col].astype(str))

    print(f"Wrote {centroids.n_obs:,} group rows × {centroids.n_vars:,} genes.")

    with timer("Write h5ad"):
        centroids.write_h5ad(output_path)

    timing_summary(f"{stem}_timing.tsv")

    session_info_path = f"{stem}_session_info.txt"
    session_info.show(write_req_file=True, req_file_name=session_info_path)
    print(f"Session info written to {session_info_path}")


if __name__ == "__main__":
    main()
