#!/usr/bin/env python3
"""
create_follicle_sdata.py - Subset a sample-level SpatialData zarr into one zarr per cell ID.

For each cell listed in cell_ids_file that belongs to the given sample, queries a
bounding box of ± radius around the cell centroid and writes a follicle-level zarr
to output/<cell_id>.zarr.

Usage:
    create_follicle_sdata.py --sample ROI1_A --path ROI1_A.zarr \
        --cell_ids_file assets/stage_quality_area_all_rois.csv --radius 100
"""

import argparse
import os
from pathlib import Path

import numpy as np
import pandas as pd
import spatialdata
from spatialdata import transform
from spatialdata.transformations import get_transformation
import session_info

from timer import timer, timing_summary


def normalize_sample_id(value: str) -> str:
    return "".join(ch for ch in str(value).upper() if ch.isalnum())


def parse_args():
    parser = argparse.ArgumentParser(
        description="Subset a sample-level SpatialData zarr into per-cell follicle zarrs"
    )
    parser.add_argument("--sample", required=True, help="Sample identifier")
    parser.add_argument("--path", required=True, dest="path", help="Path to sample-level zarr store")
    parser.add_argument("--cell_ids_file", required=True, help="CSV mapping cell IDs to samples")
    parser.add_argument("--radius", type=float, default=100.0, help="Bounding box radius in µm")
    args = parser.parse_args()
    return args


def load_cells(cell_ids_file: str, sample: str, default_radius: float) -> pd.DataFrame:
    """Read cell_ids_file and return rows matching sample, with radius filled."""
    df = pd.read_csv(cell_ids_file)
    normalized_sample = normalize_sample_id(sample)
    df["_normalized_roi"] = df["Donor.ROI"].map(normalize_sample_id)
    matching_rois = df.loc[df["_normalized_roi"] == normalized_sample, "Donor.ROI"].drop_duplicates().tolist()
    if not matching_rois:
        raise ValueError(
            f"No follicle cells found for sample '{sample}'. "
            "Check naming in the samplesheet and cell_ids_file Donor.ROI column."
        )
    if len(matching_rois) > 1:
        raise ValueError(
            f"Sample '{sample}' matched multiple Donor.ROI values after normalization: {matching_rois}"
        )
    matched_roi = matching_rois[0]
    df = df.loc[df["Donor.ROI"] == matched_roi].copy()
    if "radius" not in df.columns:
        df["radius"] = default_radius
    else:
        df["radius"] = df["radius"].fillna(default_radius)
    cells = df.drop(columns=["Donor.ROI", "_normalized_roi"]).reset_index(drop=True)
    print(f"{len(cells)} cell(s) found for sample '{sample}' (matched Donor.ROI='{matched_roi}'):")
    print(cells[["cell_id", "radius"]].to_string(index=False))
    return cells


def main():
    args = parse_args()

    with timer("Setup"):
        zarr_path = args.path
        default_radius = float(args.radius)

    with timer("Read zarr"):
        sdata = spatialdata.read_zarr(zarr_path)

    with timer("Load cell IDs"):
        cells = load_cells(args.cell_ids_file, args.sample, default_radius)

    with timer("Load cell circles"):
        circles = transform(sdata["cell_circles"], to_coordinate_system="global")
        affine = (
            get_transformation(sdata["cell_circles"], "global")
            .to_affine_matrix(input_axes=("x", "y"), output_axes=("x", "y"))
        )
        radius_scale = float(np.mean(np.abs([affine[0, 0], affine[1, 1]])))

    os.makedirs("output", exist_ok=True)

    for _, row in cells.iterrows():
        cell_id = row["cell_id"]
        radius = float(row["radius"]) * radius_scale

        with timer(f"Subset {cell_id}"):
            if cell_id not in circles.index:
                print(f"  WARNING: {cell_id} not found in cell_circles — skipping")
                continue
            centroid = circles.loc[cell_id, "geometry"]
            cx, cy = centroid.x, centroid.y
            min_coordinate = [cx - radius, cy - radius]
            max_coordinate = [cx + radius, cy + radius]

            sdata_fov = sdata.query.bounding_box(
                axes=("x", "y"),
                min_coordinate=min_coordinate,
                max_coordinate=max_coordinate,
                target_coordinate_system="global",
            )

        with timer(f"Write {cell_id}"):
            out = os.path.join("output", f"{cell_id}.zarr")
            if "table" in sdata_fov.tables:
                obs = sdata_fov["table"].obs
                obs["follicle_id"] = cell_id
                if cell_id not in obs.index:
                    print(f"  WARNING: {cell_id} not in table.obs — per-cell metadata not embedded")
                else:
                    meta = row.drop(labels=["cell_id", "radius"]).to_dict()
                    for col, val in meta.items():
                        obs.loc[cell_id, col] = val
            sdata_fov.write(out, overwrite=True)

        print(f"  {cell_id}: centroid=({cx:.1f}, {cy:.1f})  radius={radius:.1f}  →  {out}")

    print(f"\nSample:        {args.sample}")
    print(f"Input zarr:    {zarr_path}")
    print(f"Cells written: {len(cells)}")
    for _, row in cells.iterrows():
        print(f"  output/{row['cell_id']}.zarr  (radius={row['radius']}µm)")

    timing_summary()
    session_info.show()


if __name__ == "__main__":
    main()
