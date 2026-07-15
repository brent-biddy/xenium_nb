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

import scanpy as sc
import spatialdata_io
from spatialdata_io import xenium_aligned_image
from dask_image.imread import imread as dask_imread
from spatialdata.models import Image3DModel
from spatialdata.transformations import Identity
import session_info

from timer import timer, timing_summary


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

    print(f"Cells:              {len(qc_obs):,}")
    print(f"Median transcripts: {qc_obs['transcript_counts'].median():,.0f}")
    print(f"Median genes:       {qc_obs['n_genes_by_transcripts'].median():,.0f}")

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
