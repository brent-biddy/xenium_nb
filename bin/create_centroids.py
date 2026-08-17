#!/usr/bin/env python3
"""
create_centroids.py - Reduce a clustered SpatialData zarr to per-group centroids.

Sums both expression layers over an obs column, writing one row per group of cells
instead of one per cell:

    X                 each cell's CP10K profile, summed over the group
    layers["counts"]  raw counts, summed over the group

obs describes each row with `grouping`, `group`, `n_cells` and `sample`, plus
`resolution` and `cluster_v1` when the column grouped on is one of the sweep's
leiden_res_<r> columns. Requires layers["counts"], so not cluster_sdata_gpu_ooc output.

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

# The prefix cluster_sdata* writes its sweep under, `leiden_res_<r>` formatted to two
# decimals. Matched on the prefix rather than taken from a --resolutions flag: the sweep
# is whatever the cluster run that produced this zarr chose, and the zarr is the only
# thing that knows.
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


def leiden_resolution(column):
    """`leiden_res_0.60` -> `"0.60"`; None for any column that is not a sweep column.

    A row's resolution is a property OF THE COLUMN, not of how that column was chosen.
    Naming one explicitly (`--group_by leiden_res_0.60`) has to produce the same schema
    the sweep produces for it, or a single-resolution store would carry the same numbers
    as the sweep's rows and still be unreadable by a deck that selects on
    obs["resolution"].

    The float() check keeps an unrelated column that happens to start with the prefix
    from being read as a resolution. The string is returned rather than the float: the
    decks join it against assets/chosen_resolutions.csv, and "0.60" round-trips through
    a CSV where 0.6 does not.
    """
    if not column.startswith(LEIDEN_PREFIX):
        return None
    res = column.removeprefix(LEIDEN_PREFIX)
    try:
        float(res)
    except ValueError:
        return None
    return res


def sweep_columns(adata):
    """The sweep's cluster columns, as [(resolution_str, column_name), ...].

    Sorted numerically rather than lexically: at ten resolutions "0.10" sorts before
    "0.90" either way, but a sweep including "1.00" does not.
    """
    pairs = [
        (res, col)
        for col in adata.obs.columns
        if (res := leiden_resolution(col)) is not None
    ]
    if not pairs:
        raise ValueError(
            f"no {LEIDEN_PREFIX}<r> columns in obs — was this zarr written by "
            f"cluster_sdata, cluster_sdata_gpu or cluster_sdata_gpu_ooc? "
            f"Found: {sorted(adata.obs.columns)}"
        )
    return sorted(pairs, key=lambda pair: float(pair[0]))


def groupings(adata, group_by):
    """What to sum over, as [(column_name, resolution_str_or_None), ...].

    The one place the two modes differ; everything downstream of it just reduces by a
    column. `resolution` rides beside the column because it is only recoverable by
    parsing the leiden name — and it is parsed from the NAME in both modes, so
    `--group_by leiden_res_0.60` yields exactly the sweep's row for that resolution. A
    column with no resolution in its name gets None, which is what leaves
    `resolution`/`cluster_v1` absent from the store rather than placeheld with something
    that would sort and compare as if it meant anything.

    A --group_by column is validated here and not later: a typo'd column name should
    fail before the CP10K pass, which is the expensive part.
    """
    if group_by is None:
        return [(col, res) for res, col in sweep_columns(adata)]
    if group_by not in adata.obs:
        raise ValueError(
            f"--group_by {group_by!r} is not an obs column of this zarr. This step "
            f"groups by a column some upstream step already wrote; it does not compute "
            f"labels. Found: {sorted(adata.obs.columns)}"
        )
    return [(group_by, leiden_resolution(group_by))]


def size_rank(counts):
    """Group label -> size rank as a string, 1 being the largest group.

    Leiden's own ids are assignment order, so cluster 3 in one sample has nothing to do
    with cluster 3 in another, or with cluster 3 at the next resolution up. Ranking by
    cell count at least makes the id mean the same KIND of thing everywhere.

    value_counts already orders by descending count; the sort only settles ties, on the
    original numeric label. Pandas sorts value_counts with a non-stable algorithm, so two
    equal-sized clusters could otherwise swap ids between runs over the same data. int()
    because the labels are strings, where "10" < "2".
    """
    ordered = sorted(counts.index, key=lambda level: (-counts[level], int(level)))
    return {level: str(rank) for rank, level in enumerate(ordered, start=1)}


def reduce_by(adata, column, resolution):
    """Sum both layers over one obs column, returning a groups x genes AnnData.

    Two aggregate calls, not one: sc.get.aggregate takes a single `layer`, so each sum
    is its own pass over the matrix. What guarantees the two are over identical row sets
    is that both group on the SAME materialized `by_key` column — grouping twice on the
    stored categorical could disagree if its unused categories differed between calls.
    """
    # sc.get.aggregate returns its results in .layers, keyed by the func name, and
    # ignores .X entirely. `by` must NAME an obs column — passing the Series itself
    # raises a KeyError on its values — so the string cast goes through a scratch
    # column. Cast to str rather than grouped as the stored categorical so the row
    # labels are plain strings on both sides of the join the decks make onto obs.
    by_key = "_centroid_group"
    adata.obs[by_key] = adata.obs[column].astype(str)
    agg = sc.get.aggregate(adata, by=by_key, func="sum", layer="cp10k")
    counts = sc.get.aggregate(adata, by=by_key, func="sum", layer="counts")
    by = adata.obs.pop(by_key)

    out = ad.AnnData(
        X=agg.layers["sum"],
        obs=pd.DataFrame(index=agg.obs_names.astype(str)),
        var=pd.DataFrame(index=agg.var_names),
    )
    out.layers["counts"] = counts.layers["sum"]

    # The schema every row carries, whatever the grouping.
    group_counts = by.value_counts()
    out.obs["grouping"] = column
    out.obs["group"] = out.obs_names
    out.obs["n_cells"] = group_counts.reindex(out.obs_names).to_numpy()

    # The two fields a leiden column earns and nothing else does. Absent rather than
    # placeheld for a non-leiden grouping: a cell type has no resolution, and a NaN says
    # so where a 0 or an empty string would sort and compare as if it meant something.
    if resolution is not None:
        out.obs["resolution"] = resolution
        out.obs["cluster_v1"] = out.obs_names.map(size_rank(group_counts))

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

    table_key = "table"
    with timer("Extract table"):
        adata = sdata.tables[table_key].copy()

    print(f"Table:   {adata.n_obs:,} cells × {adata.n_vars:,} genes  (key: '{table_key}')")

    # Checked here, before the CP10K pass and before any grouping work, because the
    # alternative is a bare KeyError from deep inside that pass on a zarr that has
    # already been read.
    #
    # cluster_sdata and cluster_sdata_gpu both stash the raw matrix; cluster_sdata_gpu_ooc
    # deliberately does NOT, since holding a second copy of a merged cohort's matrix is
    # the thing that step exists to avoid. So its output cannot be reduced here, and the
    # message says so rather than leaving the reader to find that out.
    #
    # There is no fallback to layers["lognorm"], and that is deliberate rather than
    # unfinished. cluster_sdata* normalizes with no target_sum, so lognorm is on each
    # sample's own median scale — the one thing these centroids exist NOT to be, since
    # they are built to be compared across samples and against an external reference.
    # CP10K is in fact recoverable from it (every cell sums to the same median M after
    # normalize_total, so CP10K is expm1(lognorm) * 1e4 / M), but that round-trips a log
    # transform over the whole matrix and adds a second normalization path that has to
    # stay in step with the first. A clear error beats a subtle numerical difference
    # between samples that nothing would report.
    if "counts" not in adata.layers:
        raise KeyError(
            f"{args.sample}: the table has no layers['counts'], so centroids cannot be "
            f"built. Centroids are summed CP10K and need the raw matrix; the layers "
            f"present are {sorted(adata.layers)}. cluster_sdata and cluster_sdata_gpu "
            f"write this layer, but cluster_sdata_gpu_ooc omits it on purpose — if this "
            f"zarr came from that step, re-cluster the sample with cluster_sdata or "
            f"cluster_sdata_gpu if it fits."
        )

    to_reduce = groupings(adata, args.group_by)
    print(f"Reducing over {len(to_reduce)} column(s): "
          f"{', '.join(col for col, _ in to_reduce)}")

    with timer("CP10K"):
        # Normalized once, then reused for every grouping: the per-cell scaling does not
        # depend on which column the cells are summed over, and on a whole-slide sample
        # this is the pass worth not doing ten times.
        #
        # From layers["counts"], the raw matrix cluster_sdata stashes before normalizing.
        # Not expm1(X): scaling has overwritten X, so it no longer holds expression at
        # all, and not layers["lognorm"], which is on this sample's own median scale.
        adata.layers["cp10k"] = adata.layers["counts"].copy()
        sc.pp.normalize_total(adata, layer="cp10k", target_sum=1e4)

    blocks = []
    for column, resolution in to_reduce:
        with timer(f"Reduce {column}"):
            blocks.append(reduce_by(adata, column, resolution))

    with timer("Assemble"):
        # index_unique=None keeps the group labels as they are; the rows are made unique
        # by the (grouping, group) pair in obs, not by the index, and a suffixed index
        # would break the decks' join back onto per-cell obs. anndata still requires the
        # index to be unique, so the sweep's repeated labels are renumbered explicitly.
        centroids = ad.concat(blocks, axis=0, join="outer", index_unique=None)
        centroids.obs_names = [str(i) for i in range(centroids.n_obs)]
        # The sample id, so a cohort of these stores can be concatenated and the rows
        # still say where they came from. Every deck reading a fan-in of centroid stores
        # takes the sample from here rather than from the staged file name.
        centroids.obs["sample"] = args.sample

        # Pin the identity columns to categorical rather than letting the dtype fall out
        # of the data. Both anndata's writer and ad.concat convert object string columns
        # to categorical only when the cardinality is low enough to pay off, so `group`
        # and `cluster_v1` come back categorical from a sweep store (77 rows, 7 unique at
        # one resolution) and object from a single-resolution one (7 rows, 7 unique) —
        # same values, different dtype, decided by how many resolutions were run.
        #
        # That breaks the promise above that `--group_by leiden_res_0.60` is readable by
        # the code that reads the sweep: a categorical and an object column compare,
        # merge and groupby differently, and a categorical is what round-trips a missing
        # value as NA instead of the string "None" (the same trap create_sdata documents
        # for var["codeword_category"]). Cast explicitly so the schema is this script's
        # decision and not the cohort's shape.
        for col in ("grouping", "group", "resolution", "cluster_v1", "sample"):
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
