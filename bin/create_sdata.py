#!/usr/bin/env python3
"""
create_sdata.py - Convert raw Xenium output to a SpatialData zarr store.

Reads a Xenium output directory and writes a SpatialData zarr store containing
cells, transcripts, segmentation masks, morphology images, and optionally an
aligned H&E image and DAPI z-stack. The zarr is the primary artifact consumed
by downstream analysis notebooks.

Writes <sample>.zarr into the current working directory, alongside
timing and session info files.

Usage:
    create_sdata.py --sample ROI1_A --path /data/ROI1_A --n_jobs 4
    create_sdata.py --sample ROI1_A --path /data/ROI1_A \
        --he_image /data/he.ome.tif --he_alignment /data/he_imagealignment.csv
"""

import argparse
from pathlib import Path

import numpy as np
import scanpy as sc
import spatialdata_io
from spatialdata_io import xenium_aligned_image
from dask_image.imread import imread as dask_imread
from spatialdata.models import Image3DModel
from spatialdata.transformations import Identity
import session_info

from timer import timer, timing_summary

# The five Xenium control/codeword counters, as spatialdata_io's reader names them in
# obs. They are per-cell counts taken from cells.parquet, not features in X — which is
# why the negative-control fraction below is summed by hand rather than passed to
# calculate_qc_metrics as qc_vars.
CONTROL_COLS = [
    "control_probe_counts",
    "genomic_control_counts",
    "control_codeword_counts",
    "unassigned_codeword_counts",
    "deprecated_codeword_counts",
]


def parse_args():
    parser = argparse.ArgumentParser(
        description="Convert Xenium output to a SpatialData zarr store"
    )
    parser.add_argument("--sample", required=True, help="Sample identifier")
    parser.add_argument("--path", required=True, dest="path", help="Xenium output directory")
    parser.add_argument("--n_jobs", type=int, default=1, help="Parallel reader threads")
    parser.add_argument("--he_image", default="", help="Path to H&E OME-TIFF (optional)")
    parser.add_argument("--he_alignment", default="", help="Path to H&E alignment CSV (optional)")
    args = parser.parse_args()
    # Both H&E args must be provided together
    if bool(args.he_image) != bool(args.he_alignment):
        parser.error("--he_image and --he_alignment must be provided together")
    return args


def main():
    args = parse_args()

    output_path = f"{args.sample}.zarr"
    # morphology.ome.tif is the full DAPI z-stack (all focal planes). It is separate
    # from morphology_focus/, which contains the single best-focus plane per channel.
    morphology_3d_path = Path(args.path) / "morphology.ome.tif"

    print(f"Sample:  {args.sample}")
    print(f"Input:   {args.path}")
    print(f"Output:  {output_path}")

    # spatialdata_io.xenium() reads the standard Xenium output files: cells,
    # transcripts, segmentation masks (cell/nucleus labels), and morphology_focus
    # images. n_jobs parallelises transcript and cell reading.
    with timer("Read Xenium"):
        sdata = spatialdata_io.xenium(
            path=args.path,
            n_jobs=args.n_jobs,
        )

    sdata.tables["table"].obs["sample"] = args.sample

    # Annotate per-cell QC metrics but filter nothing: create_sdata produces the raw
    # artifact and leaves thresholds to downstream analysis, mirroring create_adata.
    # Computing them here rather than in a QC notebook means X is read once, at the
    # only point it is already in memory.
    #
    # expr_type="transcripts" is what keeps the outputs (total_transcripts,
    # n_genes_by_transcripts, pct_transcripts_in_top_N_genes) from clobbering the
    # Xenium-native total_counts that the reader takes from cells.parquet. The two
    # are NOT the same quantity: native total_counts = transcript_counts + the five
    # control/codeword counters, whereas calculate_qc_metrics would overwrite it with
    # a plain row sum of X. Only transcript_counts (== X.sum(axis=1)) is the panel
    # signal, so thresholds belong on it, not on total_counts.
    #
    # No qc_vars: the Xenium panel carries no MT- genes, so pct_transcripts_mt would
    # read 0.0 for every cell — a metric that looks real and silently passes every
    # threshold. create_adata warns about this same trap on the scRNA-seq side, where
    # the genes do exist. Ribosomal/hemoglobin sets are likewise omitted: a targeted
    # panel holds only a partial, arbitrary subset of each, so the percentages are not
    # comparable to their scRNA-seq counterparts.
    #
    # percent_top reports the share of a cell's transcripts coming from its N most
    # expressed genes. The values are lower than create_adata's (10, 20, 50, 150)
    # because a cell cannot have more genes than the panel targets: at a median 174
    # genes per cell, top_150 pins 40% of cells at exactly 100% (any cell with <= 150
    # genes detected is 100% by definition), which inverts the metric — the emptiest
    # cells score highest. These four stay below the per-cell gene count.
    #
    # log1p=False — the raw totals are what downstream thresholds are set on.
    with timer("QC metrics"):
        sc.pp.calculate_qc_metrics(
            sdata.tables["table"],
            expr_type="transcripts",
            percent_top=(5, 10, 20, 50),
            log1p=False,
            inplace=True,
        )

    # obs["total_transcripts"] is X.sum(axis=1), which for Xenium is exactly the
    # reader's transcript_counts. Drop it rather than leave obs holding three
    # near-identically named totals — two identical, one (total_counts) meaning
    # something else — for a later reader to pick the wrong one from. The per-gene
    # var["total_transcripts"] is a different quantity (a gene's total across cells)
    # and is kept.
    qc_obs = sdata.tables["table"].obs
    qc_obs.drop(columns="total_transcripts", inplace=True)

    # Negative-control burden, the one per-cell QC metric calculate_qc_metrics cannot
    # produce here. qc_vars needs its variables to be features in X, and the five
    # control/codeword counters are not — the reader lands them in obs, from
    # cells.parquet. So the sum and the fraction are written explicitly.
    #
    # Computed here for the same reason as everything above: this is the one point X
    # and the reader's own columns are together in memory, so every consumer can read
    # the metric instead of re-deriving it. Both cluster_report and qc_report used to
    # derive it themselves, which is two definitions of one number.
    #
    # The denominator is transcript_counts + controls rather than obs["total_counts"].
    # On this object the two are identical by construction, but total_counts is the
    # column cluster_sdata_gpu_ooc's calculate_qc_metrics overwrites with a plain row
    # sum, whereas transcript_counts always equals X.sum(axis=1). Building the
    # denominator from the parts is what keeps the metric meaning one thing everywhere.
    control_counts = qc_obs[CONTROL_COLS].sum(axis=1)
    total = qc_obs["transcript_counts"] + control_counts
    qc_obs["control_counts"] = control_counts
    # A cell with no transcripts and no controls has no fraction to report — 0/0 is a
    # NaN here rather than a 0 that would read as a clean cell.
    qc_obs["pct_control"] = 100 * control_counts / total.replace(0, np.nan)

    # Mean transcripts per detected gene — how concentrated a cell's counts are across
    # the genes it detects. Transcripts and genes rise together, so this is the residual
    # after that trend: a cell high on it has its counts piled into few genes, which is
    # either a very specialised cell or a segmentation artifact that swept one bright
    # neighbour's transcripts into an otherwise empty mask.
    #
    # Written the "transcripts per gene" way round rather than its reciprocal because it
    # is unbounded above, so those concentrated cells separate instead of being squashed
    # against a ceiling of 1. Not the log10(genes)/log10(transcripts) "novelty score"
    # from scRNA-seq either: that convention assumes a whole transcriptome, and on a 5K
    # panel it compresses into ~0.92-0.98 where nothing is legible.
    #
    # NaN, not 0, where a cell detects no genes: 0 would place the emptiest cells at the
    # bottom of the range next to genuinely diffuse ones.
    qc_obs["transcripts_per_gene"] = (
        qc_obs["transcript_counts"] / qc_obs["n_genes_by_transcripts"].replace(0, np.nan)
    )

    # Share of a cell's area taken up by its nucleus. A segmentation check more than an
    # expression one: the ratio sits around 0.5 on these ROIs, and cells far below it are
    # mostly cytoplasm — either genuinely large cells or a boundary that swept in
    # neighbouring space — while cells near 1 are nucleus with almost no cytoplasm around
    # it, which is what an over-tight boundary or a nucleus-expansion fallback looks like.
    #
    # The reader leaves nucleus_area as NaN, NOT 0, for a cell segmented without a
    # nucleus (70 of 21,724 on one ROI, matching nucleus_count == 0 exactly). That NaN
    # propagates here on purpose: those cells have no ratio to report, and a 0 would put
    # them at the bottom of the range beside genuinely cytoplasm-heavy cells.
    qc_obs["nucleus_ratio"] = (
        qc_obs["nucleus_area"] / qc_obs["cell_area"].replace(0, np.nan)
    )

    # Which class of codebook entry each gene's probe belongs to — predesigned catalogue
    # panel vs. separately-designed custom probes, on the gene axis. That distinction is
    # the one question the per-gene abundance and detection metrics cannot answer: whether
    # custom probes land in the same abundance-detection cloud as the catalogue ones, or
    # off in a corner, which would be a probe-design problem rather than biology.
    #
    # This is a LOOKUP, not a summary, which is what makes it safe to compute here.
    # feature_name -> codeword_category is strictly one-to-one (verified: 0 of 7,415
    # features on a test ROI carry two categories), so collapsing the molecule table to
    # one row per feature loses nothing and cannot disagree with a downstream derivation.
    #
    # drop_duplicates rather than groupby().first(): both give one row per feature, but
    # drop_duplicates is a single pass and stays cheap under Dask, which matters when a
    # whole-slide sample has hundreds of millions of molecules. Only the two columns are
    # read — the coordinates and ids are most of the table's width.
    #
    # Only the gene categories survive onto var, since controls are not features in X:
    # the reader lands them in obs as the five counter columns summed above. So this
    # column takes two values plus the unmapped state below.
    with timer("Codeword category"):
        tx = sdata.points["transcripts"][["feature_name", "codeword_category"]]
        lookup = tx.drop_duplicates().compute()
        lookup = lookup.astype({"feature_name": str, "codeword_category": str})
        category = (lookup.set_index("feature_name")["codeword_category"]
                          .reindex(sdata.tables["table"].var.index))

    # A gene with no transcript ANYWHERE in the sample has no row in the molecule table
    # to read a category from, so reindex leaves a missing value. That is a real state —
    # the probe is on the panel and detected nothing — not a lookup failure, and it is
    # kept rather than filled.
    #
    # Normalised to object dtype with None HERE, once, rather than left as pandas
    # "string" dtype: pd.NA == "predesigned_gene" evaluates to pd.NA rather than False,
    # so every comparison-built mask downstream raises "boolean value of NA is ambiguous"
    # instead of selecting rows. Do not reintroduce the string dtype on this column.
    category = category.astype(object).where(category.notna(), None)
    sdata.tables["table"].var["codeword_category"] = category

    n_unmapped = int(category.isna().sum())

    print(f"Cells:              {len(qc_obs):,}")
    print(f"Median transcripts: {qc_obs['transcript_counts'].median():,.0f}")
    print(f"Median genes:       {qc_obs['n_genes_by_transcripts'].median():,.0f}")
    print(f"Gene categories:    "
          f"{dict(category.dropna().value_counts())}, unmapped {n_unmapped:,}")
    print(f"Cells w/ control:   {(control_counts > 0).sum():,} "
          f"({100 * (control_counts > 0).mean():.2f}%)")
    print(f"Mean control %:     {qc_obs['pct_control'].mean():.4f}")

    # spatialdata_io auto-detects an H&E image if one is named with the expected
    # Xenium suffix alongside the data. If not auto-detected, load it explicitly
    # using the provided image path and alignment matrix.
    if "he_image" not in sdata.images and args.he_image and args.he_alignment:
        with timer("Load H&E"):
            # imread reads only the base level of the OME-TIFF pyramid; scale_factors
            # rebuilds it in the zarr. 4 halvings reaches a screen-sized resolution.
            he = xenium_aligned_image(
                image_path=args.he_image,
                alignment_file=args.he_alignment,
                image_models_kwargs={
                    "scale_factors": [2, 2, 2, 2],
                },
            )
        sdata.images["he_image"] = he
        print(f"Loaded H&E from {args.he_image}")
    elif "he_image" in sdata.images:
        print("H&E auto-detected by spatialdata_io.")
    else:
        print("No H&E image found.")

    # morphology.ome.tif is not loaded by the xenium() reader — it adds the full
    # z-stack separately so downstream notebooks can inspect individual focal planes.
    if morphology_3d_path.exists():
        with timer("Add DAPI z-stack"):
            dapi_3d = dask_imread(str(morphology_3d_path))
            sdata.images["dapi_3d"] = Image3DModel.parse(
                dapi_3d[None],  # imread returns (z, y, x); [None] adds the required c axis → (c, z, y, x)
                dims=("c", "z", "y", "x"),
                c_coords=["DAPI"],
                transformations={"global": Identity()},
                # imread reads only the base level; scale_factors rebuilds the pyramid.
                # y/x only — z is not downsampled. 4 halvings reaches a screen-sized resolution.
                scale_factors=[{"y": 2, "x": 2}, {"y": 2, "x": 2}, {"y": 2, "x": 2}, {"y": 2, "x": 2}],
                # Without explicit chunks, building the pyramid above falls back to
                # multiscale_spatial_image's default_chunks=64, so every level is written
                # as 64x64 tiles — one ~90 KB zarr chunk file each, ~10k files for a
                # single image, enough to threaten an HPC inode quota across a cohort.
                # One chunk per z-plane matches how notebooks read this element (a focal
                # plane at a time) and mirrors the chunking spatialdata_io already applies
                # to morphology_focus. Each value is capped at that axis's extent, so
                # smaller images degrade to one chunk per plane rather than padding.
                chunks=(1, 1, 4096, 4096),
            )
        print(f"Loaded DAPI z-stack from {morphology_3d_path}")
    else:
        print("Skipping DAPI z-stack (morphology.ome.tif not found).")

    with timer("Write zarr"):
        sdata.write(output_path, overwrite=True)
    print(f"Written to {output_path}")

    # Print every element in the sdata object
    print("\nElements:")
    for group_name in ("images", "labels", "points", "shapes", "tables"):
        group = getattr(sdata, group_name, {})
        for name, element in group.items():
            print(f"  {name}: {type(element).__name__} [{group_name}]")

    timing_summary(path=f"{args.sample}_timing.tsv")

    session_info_path = f"{args.sample}_session_info.txt"
    session_info.show(write_req_file=True, req_file_name=session_info_path)
    print(f"Session info written to {session_info_path}")


if __name__ == "__main__":
    main()
