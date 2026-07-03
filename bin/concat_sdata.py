#!/usr/bin/env python3
"""
concat_sdata.py - Concatenate multiple SpatialData zarr stores into one.

Reads each input zarr with spatialdata.read_zarr(), merges them with
spatialdata.concatenate(), and writes the result to merged.zarr.

Usage:
    concat_sdata.py --paths ROI1_A.zarr ROI1_B.zarr ROI2_A.zarr
"""

import argparse
import os

import spatialdata
import session_info

from timer import timer, timing_summary


def parse_args():
    parser = argparse.ArgumentParser(
        description="Concatenate multiple SpatialData zarr stores into one"
    )
    parser.add_argument("--paths", required=True, nargs="+", help="Input zarr store paths")
    return parser.parse_args()


def main():
    args = parse_args()

    output_path = "merged.zarr"

    print(f"Output:      {output_path}")
    print(f"Inputs ({len(args.paths)}):")
    for p in args.paths:
        print(f"  {p}")

    # Keys become element-name prefixes, so identical element names across samples
    # (e.g. morphology_focus, cell_labels) don't collide.
    sdata_dict = {}
    for path in args.paths:
        with timer(f"Read {path}"):
            sdata = spatialdata.read_zarr(path)
        key = sdata.tables["table"].obs["sample"].iloc[0]
        sdata_dict[key] = sdata

    with timer("Concatenate"):
        # concatenate_tables merges the per-sample "table" elements into a
        # single "table" (obs distinguished by the existing "sample" column)
        # instead of suffixing them into table-<key>, so cluster_sdata's
        # hardcoded table_key = "table" can read the merged object.
        merged = spatialdata.concatenate(sdata_dict, concatenate_tables=True)

    with timer("Write zarr"):
        merged.write(output_path, overwrite=True)
    print(f"Written to {output_path}")

    timing_summary(path="concat_sdata_timing.tsv")

    session_info_path = "concat_sdata_session_info.txt"
    session_info.show(write_req_file=True, req_file_name=session_info_path)
    print(f"Session info written to {session_info_path}")


if __name__ == "__main__":
    main()
