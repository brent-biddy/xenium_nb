#!/usr/bin/env python3
"""
annotate_sdata.py - Label Leiden clusters with cell types from a marker YAML.

Reads a clustered SpatialData zarr and a species-specific marker YAML (see
assets/ovarian_markers_human.yaml), scores every cell against each cell type's
marker set with scanpy's score_genes, then assigns each Leiden cluster the cell
type with the highest mean score across its cells.

Labels are assigned per cluster, not per cell: single-cell scores on a targeted
Xenium panel are too sparse to argmax reliably, and clusters are the unit the
downstream reports already work in.

Because the cluster_sdata* steps sweep resolutions and record no plain `leiden`
column, this annotates every leiden_res_* column it finds unless --resolutions
narrows it. Each one yields a matching cell_type_res_* column.

Writes annotated.zarr and annotate_sdata_scores.tsv into the working directory.

Usage:
    annotate_sdata.py --sample ROI1_A --path clustered.zarr \
        --markers assets/ovarian_markers_human.yaml
    annotate_sdata.py --sample ROI1_A --path clustered.zarr \
        --markers assets/ovarian_markers_human.yaml --resolutions 0.4 0.6 \
        --exclude_nonspecific
"""

import argparse

import pandas as pd
import scanpy as sc
import spatialdata
import yaml

from timer import timer, timing_summary

UNASSIGNED = "unassigned"


def parse_args():
    parser = argparse.ArgumentParser(
        description="Annotate Leiden clusters in a SpatialData zarr from a marker YAML"
    )
    parser.add_argument("--sample", required=True, help="Sample identifier")
    parser.add_argument(
        "--path", required=True, help="Path to input clustered SpatialData zarr"
    )
    parser.add_argument(
        "--markers", required=True, help="Path to the marker YAML for this species"
    )
    parser.add_argument(
        "--resolutions",
        type=float,
        nargs="+",
        default=None,
        metavar="RES",
        help="Leiden resolutions to annotate; defaults to every leiden_res_* "
        "column present in the table.",
    )
    parser.add_argument(
        "--exclude_nonspecific",
        action="store_true",
        help="Drop markers the YAML flags as also labelling other cell types. "
        "Off by default: the flag is unevenly annotated in the source sheet, so "
        "excluding shrinks some marker sets far more than others.",
    )
    parser.add_argument(
        "--min_genes",
        type=int,
        default=3,
        help="Skip cell types with fewer than this many markers on the panel "
        "(default: 3). Scores over one or two genes are noise.",
    )
    parser.add_argument(
        "--min_score",
        type=float,
        default=0.0,
        help=f"Clusters whose winning cell type has a raw mean score at or below "
        f"this are labelled '{UNASSIGNED}' (default: 0.0). score_genes centres a "
        "random gene set near zero, so a non-positive score means no enrichment.",
    )
    parser.add_argument(
        "--scale",
        choices=("zscore", "none"),
        default="zscore",
        help="How to compare cell types within a cluster. 'zscore' (default) "
        "standardizes each cell type's mean scores across clusters before "
        "ranking, so the label answers 'which type is this cluster most enriched "
        "for relative to other clusters'. 'none' ranks raw means, which lets a "
        "cell type with broadly-expressed markers win nearly every cluster.",
    )
    return parser.parse_args()


def main():
    args = parse_args()

    output_path = "annotated.zarr"
    scores_path = "annotate_sdata_scores.tsv"

    print(f"Sample:  {args.sample}")
    print(f"Input:   {args.path}")
    print(f"Markers: {args.markers}")
    print(f"Output:  {output_path}")

    with timer("Read markers"):
        spec = yaml.safe_load(open(args.markers))
        cell_types = spec["cell_types"]

    print(f"Species: {spec.get('species', 'unspecified')}")
    print(f"Loaded {len(cell_types)} cell types from {args.markers}")

    with timer("Read zarr"):
        sdata = spatialdata.read_zarr(args.path)

    table_key = "table"
    with timer("Extract table"):
        adata = sdata.tables[table_key].copy()

    print(f"Table:   {adata.n_obs:,} cells x {adata.n_vars:,} genes  (key: '{table_key}')")

    # score_genes compares each set against a background of similarly-expressed
    # genes, which only means anything on normalized data. cluster_sdata* leaves
    # X normalized and log1p'd; warn rather than fail if that changes, since a
    # caller may legitimately pass a differently-prepared table.
    if "log1p" not in adata.uns:
        print(
            "WARNING: table has no 'log1p' entry in uns - X may be raw counts. "
            "score_genes expects normalized, log-transformed data."
        )

    # Resolve which cluster columns to annotate.
    available = sorted(c for c in adata.obs.columns if c.startswith("leiden_res_"))
    if not available:
        raise SystemExit(
            f"No leiden_res_* columns in {args.path}. Run a cluster_sdata* step first."
        )
    if args.resolutions is None:
        cluster_keys = available
    else:
        cluster_keys = [f"leiden_res_{r:.2f}" for r in sorted(args.resolutions)]
        missing = [k for k in cluster_keys if k not in available]
        if missing:
            raise SystemExit(
                f"Requested resolutions not present: {', '.join(missing)}. "
                f"Available: {', '.join(available)}"
            )

    print(f"Annotating {len(cluster_keys)} resolution(s): {', '.join(cluster_keys)}")

    # Restrict each marker set to genes actually on the panel. score_genes drops
    # missing genes silently, so an unchecked set can collapse to a handful of
    # genes - or none - without any visible sign in the output.
    panel = set(adata.var_names)
    scored = {}
    print(f"\n{'cell type':38} {'on panel':>10}")
    for name, entry in sorted(cell_types.items()):
        genes = entry["markers"]
        if args.exclude_nonspecific:
            genes = [g for g in genes if g not in (entry.get("nonspecific") or {})]
        present = [g for g in genes if g in panel]
        flag = "" if len(present) >= args.min_genes else "  SKIPPED"
        print(f"  {name:36} {len(present):4}/{len(genes):<4}{flag}")
        if len(present) >= args.min_genes:
            scored[name] = present

    if not scored:
        raise SystemExit(
            "No cell type has enough markers on this panel. Check that the YAML "
            "species matches the data, or lower --min_genes."
        )
    print(f"\nScoring {len(scored)} of {len(cell_types)} cell types.")

    with timer("Score cell types"):
        for name, genes in scored.items():
            sc.tl.score_genes(adata, genes, score_name=f"score_{name}", random_state=0)

    score_cols = [f"score_{name}" for name in scored]

    # Per cluster: mean of each cell type's score across the cluster's cells,
    # then pick the winner. Raw score_genes values are not comparable between
    # gene sets - a set whose markers are broadly expressed scores high in every
    # cluster and wins them all - so by default rank on each cell type's mean
    # standardized across clusters, and keep the raw mean only as the
    # is-this-enriched-at-all gate.
    summaries = []
    for cluster_key in cluster_keys:
        with timer(f"Assign {cluster_key}"):
            means = adata.obs.groupby(cluster_key, observed=True)[score_cols].mean()

            if args.scale == "zscore" and means.shape[0] > 1:
                # ddof=0: these clusters are the whole population, not a sample.
                spread = means.std(ddof=0).replace(0.0, pd.NA)
                ranked = (means - means.mean()) / spread
            else:
                if args.scale == "zscore":
                    print(
                        f"  {cluster_key}: only one cluster, cannot standardize "
                        "across clusters; ranking raw means instead."
                    )
                ranked = means

            best = ranked.idxmax(axis=1).str.removeprefix("score_")
            # Gate on the raw mean, not the ranked value: z-scores are centred by
            # construction, so every cell type is the relative winner somewhere.
            winning_raw = means.to_numpy()[
                range(len(means)), means.columns.get_indexer(ranked.idxmax(axis=1))
            ]
            best[winning_raw <= args.min_score] = UNASSIGNED

            label_key = cluster_key.replace("leiden_res_", "cell_type_res_")
            adata.obs[label_key] = (
                adata.obs[cluster_key].map(best).astype("category")
            )

        counts = adata.obs[label_key].value_counts()
        n_unassigned = int(counts.get(UNASSIGNED, 0))
        print(
            f"{cluster_key}: {means.shape[0]} clusters -> "
            f"{best.nunique()} labels, {n_unassigned:,} cells {UNASSIGNED}"
        )

        # Long-form so every resolution stacks into one inspectable table. Both
        # the raw mean and the value actually ranked on are kept, so a surprising
        # label can be traced back to which of the two drove it.
        tidy = means.reset_index().melt(
            id_vars=cluster_key, var_name="cell_type", value_name="mean_score"
        )
        tidy["ranked_score"] = ranked.reset_index().melt(
            id_vars=cluster_key, var_name="cell_type", value_name="v"
        )["v"]
        tidy["cell_type"] = tidy["cell_type"].str.removeprefix("score_")
        tidy = tidy.rename(columns={cluster_key: "cluster"})
        tidy.insert(0, "resolution", cluster_key.removeprefix("leiden_res_"))
        tidy.insert(0, "sample", args.sample)
        summaries.append(tidy)

    with timer("Write scores"):
        pd.concat(summaries, ignore_index=True).to_csv(
            scores_path, sep="\t", index=False
        )

    with timer("Write zarr"):
        sdata.tables[table_key] = adata
        sdata.write(output_path)

    print(f"\nWritten to {output_path} and {scores_path}")

    timing_summary(path="annotate_sdata_timing.tsv")


if __name__ == "__main__":
    main()
