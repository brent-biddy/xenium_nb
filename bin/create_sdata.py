#!/usr/bin/env python3
"""
create_sdata.py - Convert raw Xenium output to a SpatialData zarr store.

Writes output/<sample>.zarr into the current working directory.
Timing and session info are printed to stdout.

Usage:
    create_sdata.py --sample ROI1_A --path /data/ROI1_A --n_jobs 4
    create_sdata.py --sample ROI1_A --path /data/ROI1_A \
        --he_image /data/he.ome.tif --he_alignment /data/he_imagealignment.csv
"""

import argparse
import os
from pathlib import Path

import spatialdata_io
from dask_image.imread import imread as dask_imread
from spatialdata.models import Image3DModel
from spatialdata.transformations import Identity
import session_info

from timer import timer, timing_summary


def prepare_xenium_input(path_str: str) -> str:
    """Return a spatialdata_io-compatible path, building a compat symlink tree if needed.

    Some downsampled Xenium outputs use channel-prefixed filenames (ch0_, ch1_, …)
    rather than the sequentially numbered names spatialdata_io expects. This
    function detects that case and builds a sibling directory with the expected
    naming via symlinks, leaving the original untouched.
    """
    path = Path(path_str).resolve()
    focus_dir = path / "morphology_focus"
    if not focus_dir.is_dir():
        return str(path)

    files = sorted(p for p in focus_dir.glob("*.ome.tif"))
    if not files:
        return str(path)

    expected = {f"morphology_focus_{i:04d}.ome.tif" for i in range(len(files))}
    actual = {p.name for p in files}
    if actual == expected:
        return str(path)

    channel_files = sorted(p for p in files if p.name.startswith("ch"))
    if not channel_files:
        return str(path)

    compat_root = Path("input_compat")
    compat_root.mkdir(exist_ok=True)

    for child in path.iterdir():
        target = compat_root / child.name
        if target.exists():
            continue
        if child.name != "morphology_focus":
            target.symlink_to(child.resolve(), target_is_directory=child.is_dir())
            continue

        target.mkdir()
        for idx, image in enumerate(channel_files):
            alias = target / f"morphology_focus_{idx:04d}.ome.tif"
            if not alias.exists():
                alias.symlink_to(image.resolve())

    return str(compat_root.resolve())


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
    if bool(args.he_image) != bool(args.he_alignment):
        parser.error("--he_image and --he_alignment must be provided together")
    return args


def main():
    args = parse_args()

    with timer("Setup"):
        xenium_path = prepare_xenium_input(args.path)
        output_path = os.path.join("output", f"{args.sample}.zarr")

    with timer("Read Xenium"):
        sdata = spatialdata_io.xenium(
            path=xenium_path,
            n_jobs=args.n_jobs,
            cells_as_circles=True,
        )

    if "he_image" not in sdata.images and args.he_image and args.he_alignment:
        from spatialdata_io import xenium_aligned_image
        with timer("Load H&E"):
            he = xenium_aligned_image(
                image_path=args.he_image,
                alignment_file=args.he_alignment,
                image_models_kwargs={
                    "chunks": {"y": 2048, "x": 2048, "c": -1},
                    "scale_factors": [2, 2, 2, 2],
                },
            )
        sdata.images["he_image"] = he
        print(f"Loaded H&E from {args.he_image}")
    elif "he_image" in sdata.images:
        print("H&E auto-detected by spatialdata_io.")
    else:
        print("No H&E image found.")

    morphology_3d_path = Path(xenium_path) / "morphology.ome.tif"
    if morphology_3d_path.exists():
        with timer("Add DAPI z-stack"):
            dapi_3d = dask_imread(str(morphology_3d_path))
            sdata.images["dapi_3d"] = Image3DModel.parse(
                dapi_3d[None],
                dims=("c", "z", "y", "x"),
                c_coords=["DAPI"],
                transformations={"global": Identity()},
                scale_factors=[{"y": 2, "x": 2}, {"y": 2, "x": 2}, {"y": 2, "x": 2}, {"y": 2, "x": 2}],
            )
        print(sdata.images["dapi_3d"])
    else:
        print("Skipping DAPI z-stack (morphology.ome.tif not found).")

    with timer("Write zarr"):
        os.makedirs("output", exist_ok=True)
        sdata.write(output_path, overwrite=True)
    print(f"Written to {output_path}")

    print(f"\nSample:      {args.sample}")
    print(f"Input path:  {args.path}")
    print(f"Reader path: {xenium_path}")
    print(f"Output zarr: {output_path}")
    print("\nElements:")
    for group_name in ("images", "labels", "points", "shapes", "tables"):
        group = getattr(sdata, group_name, {})
        for name, element in group.items():
            print(f"  {name}: {type(element).__name__} [{group_name}]")

    timing_summary()
    session_info.show()


if __name__ == "__main__":
    main()
