#!/usr/bin/env python3
"""
Generate a follicle-level samplesheet for notebook 03_plot_follicle.qmd.

Each row corresponds to one annotated follicle cell. The sample_id is
constructed as <Donor.ROI>_<cell_id> (e.g. ROI1_aaaaimck-1) so that
notebook 03 can split on the first underscore to recover the ROI ID and
cell ID, then derive the follicle zarr path as:
    <data_path>/<roi_id>/02_subset_follicle/output/<cell_id>.zarr

data_path is set to the pipeline outdir so it is the same for all rows
belonging to the same pipeline run. All other columns from the cell IDs
file (Stage.Labels, Quality, etc.) are carried over as notebook params.

Usage:
    python bin/make_follicle_samplesheet.py \\
        --cell-ids assets/stage_quality_area_all_rois.csv \\
        --outdir   results \\
        --output   assets/follicle_samplesheet.csv
"""

import argparse
import pandas as pd


def parse_args():
    p = argparse.ArgumentParser()
    p.add_argument("--cell-ids", required=True, help="Cell IDs CSV (Donor.ROI, cell_id, ...)")
    p.add_argument("--outdir",   required=True, help="Pipeline outdir used as data_path (e.g. results)")
    p.add_argument("--output",   required=True, help="Path to write follicle samplesheet")
    return p.parse_args()


def main():
    args = parse_args()

    cell_ids = pd.read_csv(args.cell_ids)

    # Combine ROI and cell ID into a single sample_id that encodes both
    cell_ids["sample_id"] = cell_ids["Donor.ROI"] + "_" + cell_ids["cell_id"]

    # Explicit roi_id and cell_id columns let notebook 03 use them directly
    # (no fragile string splitting) and let main.nf group outputs by ROI.
    cell_ids["roi_id"] = cell_ids["Donor.ROI"]

    # data_path is the base results directory — notebook 03 derives the
    # zarr path from roi_id + cell_id + data_path at runtime
    cell_ids["data_path"] = args.outdir

    cell_ids = cell_ids.drop(columns=["Donor.ROI"])

    # Ensure key columns come first
    cols = ["sample_id", "roi_id", "cell_id", "data_path"] + [
        c for c in cell_ids.columns
        if c not in ("sample_id", "roi_id", "cell_id", "data_path")
    ]
    cell_ids[cols].to_csv(args.output, index=False)
    print(f"Wrote {len(cell_ids)} rows to {args.output}")


if __name__ == "__main__":
    main()
