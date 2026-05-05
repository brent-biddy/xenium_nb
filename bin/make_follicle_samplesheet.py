#!/usr/bin/env python3
"""
Generate a follicle-level samplesheet for manual or legacy workflows.

The main Nextflow workflow no longer requires a follicle-level samplesheet for
downstream follicle notebooks; it derives per-cell work items from the sample
sheet plus the cell IDs file. This script remains useful when exporting a
standalone follicle-level sheet for ad hoc runs or external tooling.

Each row corresponds to one annotated follicle cell. The output is a two-column
samplesheet with:

    sample = <Donor.ROI>_<cell_id> (e.g. ROI1_aaaaimck-1)
    path   = path to the upstream artifact for that follicle

The 'path' column points at the exact upstream zarr each row needs, so
Nextflow stages it directly into the notebook work dir:

    --upstream create_sdata
        path = <outdir>/<roi_id>/create_sdata/output/<roi_id>.zarr
        (use this to feed notebook subset_follicle)

    --upstream create_follicle_sdata
        path = <outdir>/<roi_id>/create_follicle_sdata/output/<cell_id>.zarr
        (use this to feed notebook plot_follicle)

All other columns from the cell IDs file (Stage.Labels, Quality, etc.) are
carried over as notebook params.

Usage:
    python bin/make_follicle_samplesheet.py \\
        --cell-ids assets/stage_quality_area_all_rois.csv \\
        --outdir   results \\
        --upstream create_sdata \\
        --output   assets/follicle_samplesheet.csv
"""

import argparse
from pathlib import Path
import pandas as pd


UPSTREAM_CHOICES = (
    "create_sdata",
    "create_follicle_sdata",
    # Backward-compatible aliases for older workflow naming.
    "create_spatialdata",
    "subset_follicle",
)


def parse_args():
    p = argparse.ArgumentParser()
    p.add_argument("--cell-ids", required=True, help="Cell IDs CSV (Donor.ROI, cell_id, ...)")
    p.add_argument("--outdir",   required=True, help="Pipeline outdir (e.g. results)")
    p.add_argument("--upstream", required=True, choices=UPSTREAM_CHOICES,
                   help="Which upstream notebook's zarr each row should point at")
    p.add_argument("--output",   required=True, help="Path to write follicle samplesheet")
    return p.parse_args()


def zarr_path(outdir: str, upstream: str, roi_id: str, cell_id: str) -> str:
    if upstream in ("create_sdata", "create_spatialdata"):
        base = Path(outdir) / roi_id / "create_sdata" / "output"
        return str(base / f"{roi_id}.zarr")
    if upstream in ("create_follicle_sdata", "subset_follicle"):
        base = Path(outdir) / roi_id / "create_follicle_sdata" / "output"
        return str(base / f"{cell_id}.zarr")
    raise ValueError(f"Unknown upstream: {upstream}")


def main():
    args = parse_args()

    cell_ids = pd.read_csv(args.cell_ids)

    out = pd.DataFrame({
        "sample": cell_ids["Donor.ROI"] + "_" + cell_ids["cell_id"],
        "path": [
            zarr_path(args.outdir, args.upstream, roi, cid)
            for roi, cid in zip(cell_ids["Donor.ROI"], cell_ids["cell_id"])
        ],
    })

    out.to_csv(args.output, index=False)
    print(f"Wrote {len(out)} rows to {args.output}")


if __name__ == "__main__":
    main()
