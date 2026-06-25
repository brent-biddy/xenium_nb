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

import pandas as pd
import spatialdata
from spatialdata import transform
import session_info

from timer import timer, timing_summary

# Xenium pixel size in µm/pixel — a fixed instrument constant used to convert
# the bounding box radius from µm (user-facing) to pixels (coordinate system units).
XENIUM_PIXEL_SIZE_UM = 0.2125


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
    cells = df.loc[df["Donor.ROI"] == sample].copy().reset_index(drop=True)
    if cells.empty:
        raise ValueError(
            f"No follicle cells found for sample '{sample}'. "
            "Check naming in the samplesheet and cell_ids_file Donor.ROI column."
        )
    # radius column is optional in the CSV; fall back to the CLI default when absent or blank.
    if "radius" not in cells.columns:
        cells["radius"] = default_radius
    else:
        cells["radius"] = cells["radius"].fillna(default_radius)
    print(f"{len(cells)} cell(s) found for sample '{sample}':")
    print(cells[["cell_id", "radius"]].to_string(index=False))
    return cells


def main():
    args = parse_args()

    zarr_path = args.path
    default_radius = float(args.radius)

    print(f"Sample:   {args.sample}")
    print(f"Input:    {zarr_path}")
    print(f"Output:   output/")

    # read_zarr opens the store lazily — array data is not loaded into memory
    # until accessed. This keeps startup fast even for large whole-sample zarrs.
    with timer("Read zarr"):
        sdata = spatialdata.read_zarr(zarr_path)

    with timer("Load cell IDs"):
        cells = load_cells(args.cell_ids_file, args.sample, default_radius)

    with timer("Load cell circles"):
        # cell_circles is a GeoDataFrame of Shapely Point geometries, one per cell,
        # stored in the native Xenium coordinate system (µm). transform() applies
        # the registered transformation to bring all coordinates into the shared
        # "global" space (pixels), so centroid lookups are consistent across all elements.
        circles = transform(sdata["cell_circles"], to_coordinate_system="global")
        # Global coordinates are in pixels; radius is in µm. Convert using the fixed
        # Xenium pixel size so the bounding box window has the correct physical size.
        radius_px_per_um = 1.0 / XENIUM_PIXEL_SIZE_UM

    os.makedirs("output", exist_ok=True)

    for _, row in cells.iterrows():
        cell_id = row["cell_id"]
        radius = float(row["radius"]) * radius_px_per_um

        with timer(f"Subset {cell_id}"):
            if cell_id not in circles.index:
                print(f"  WARNING: {cell_id} not found in cell_circles — skipping")
                continue
            centroid = circles.loc[cell_id, "geometry"]
            cx, cy = centroid.x, centroid.y
            min_coordinate = [cx - radius, cy - radius]
            max_coordinate = [cx + radius, cy + radius]

            # bounding_box() returns a new SpatialData object containing only the
            # elements (images, labels, points, shapes, table rows) that overlap
            # the query window. Images and labels are spatially cropped; points and
            # shapes are filtered to those within the box.
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
                # Tag every cell in this follicle's table with the follicle ID so
                # downstream notebooks can identify which follicle a cell belongs to.
                obs["follicle_id"] = cell_id
                if cell_id not in obs.index:
                    print(f"  WARNING: {cell_id} not in table.obs — per-cell metadata not embedded")
                else:
                    # Embed any extra columns from the cell_ids_file (e.g. stage,
                    # quality score) directly into the index cell's obs row so the
                    # metadata travels with the zarr artifact.
                    meta = row.drop(labels=["cell_id", "radius"]).to_dict()
                    for col, val in meta.items():
                        obs.loc[cell_id, col] = val
            sdata_fov.write(out, overwrite=True)

        print(f"  {cell_id}: centroid=({cx:.1f}, {cy:.1f})  radius={radius:.1f}  →  {out}")

    print(f"\nCells written: {len(cells)}")
    for _, row in cells.iterrows():
        print(f"  output/{row['cell_id']}.zarr  (radius={row['radius']}µm)")

    timing_summary(path=f"output/{args.sample}_timing.tsv")

    session_info_path = f"output/{args.sample}_session_info.txt"
    session_info.show(write_req_file=True, req_file_name=session_info_path)
    print(f"Session info written to {session_info_path}")


if __name__ == "__main__":
    main()
