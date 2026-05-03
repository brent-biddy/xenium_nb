#!/usr/bin/env python3
"""
Generate a follicle-level samplesheet for notebook 03_plot_follicle.qmd.

Each row corresponds to one annotated follicle cell. The 'sample' column is
constructed as <Donor.ROI>_<cell_id> (e.g. ROI1_aaaaimck-1) for unique
identification, and 'roi_id' / 'cell_id' columns are emitted explicitly so
notebook 03 can derive the follicle zarr path as:
    <path>/<roi_id>/02_subset_follicle/output/<cell_id>.zarr

'path' is set to the pipeline outdir so it is the same for all rows
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
    p.add_argument("--outdir",   required=True, help="Pipeline outdir used as the 'path' column (e.g. results)")
    p.add_argument("--output",   required=True, help="Path to write follicle samplesheet")
    return p.parse_args()


def main():
    args = parse_args()

    cell_ids = pd.read_csv(args.cell_ids)

    # Combine ROI and cell ID into a unique 'sample' identifier
    cell_ids["sample"] = cell_ids["Donor.ROI"] + "_" + cell_ids["cell_id"]

    # Explicit roi_id and cell_id columns let notebook 03 use them directly
    # (no fragile string splitting) and let main.nf group outputs by ROI.
    cell_ids["roi_id"] = cell_ids["Donor.ROI"]

    # 'path' is the base results directory — notebook 03 derives the
    # follicle zarr path from roi_id + cell_id + path at runtime
    cell_ids["path"] = args.outdir

    cell_ids = cell_ids.drop(columns=["Donor.ROI"])

    # Ensure key columns come first
    cols = ["sample", "roi_id", "cell_id", "path"] + [
        c for c in cell_ids.columns
        if c not in ("sample", "roi_id", "cell_id", "path")
    ]
    cell_ids[cols].to_csv(args.output, index=False)
    print(f"Wrote {len(cell_ids)} rows to {args.output}")


if __name__ == "__main__":
    main()
