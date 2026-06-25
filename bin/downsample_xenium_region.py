#!/usr/bin/env python3
"""
downsample_xenium_region.py - Crop a Xenium output directory to one or more regions.

The original downsample_xenium.py keeps the full morphology field of view. This
utility instead selects cells inside a bounding box, keeps transcripts inside the
same box, rebases spatial coordinates to the crop origin, and crops morphology
images and cells.zarr masks so the output is physically smaller.

Usage:
    python bin/downsample_xenium_region.py /path/to/xenium_output \
        --bbox 1000 2000 2500 3500 --region_name region_a

    python bin/downsample_xenium_region.py /path/to/xenium_output \
        --regions_csv regions.csv

regions.csv must contain: region,xmin,ymin,xmax,ymax
"""

import argparse
import gc
import gzip
import json
import shutil
import sys
import tempfile
import time
import zipfile
from dataclasses import dataclass
from pathlib import Path

import h5py
import numpy as np
import pandas as pd
import pyarrow as pa
import pyarrow.compute as pc
import pyarrow.parquet as pq
import zarr

from downsample_xenium import (
    archive_directory_as_zip,
    close_zarr_group,
    create_zarr_dataset,
    open_directory_group,
    open_zip_read_store,
    process_analysis,
    process_analysis_zarr,
    process_cell_feature_matrix_dir,
    process_cell_feature_matrix_h5,
    process_cell_feature_matrix_tar,
    process_cfm_zarr,
    validate_output,
)


@dataclass(frozen=True)
class Region:
    name: str
    xmin: float
    ymin: float
    xmax: float
    ymax: float

    @property
    def width(self):
        return self.xmax - self.xmin

    @property
    def height(self):
        return self.ymax - self.ymin


def parse_args():
    parser = argparse.ArgumentParser(
        description="Crop/downsample Xenium output files to one or more bounding boxes."
    )
    parser.add_argument("input_dir", type=Path, help="Path to Xenium output directory")
    parser.add_argument(
        "--output_dir",
        type=Path,
        default=None,
        help="Output directory root; defaults to <input_dir>_region_downsampled",
    )
    parser.add_argument(
        "--bbox",
        nargs=4,
        type=float,
        metavar=("XMIN", "YMIN", "XMAX", "YMAX"),
        help="Single crop bounding box in Xenium coordinate units.",
    )
    parser.add_argument(
        "--region_name",
        default="region",
        help="Name for the single --bbox output directory.",
    )
    parser.add_argument(
        "--regions_csv",
        type=Path,
        help="CSV with columns region,xmin,ymin,xmax,ymax for multiple crops.",
    )
    parser.add_argument(
        "--proportion",
        type=float,
        default=1.0,
        help="Fraction of in-region cells to keep after cropping (default: 1.0).",
    )
    parser.add_argument(
        "--grid_size",
        type=float,
        default=100.0,
        help="Grid size for optional spatial cell downsampling (default: 100.0).",
    )
    parser.add_argument(
        "--pixel_size",
        type=float,
        default=None,
        help=(
            "Coordinate units per image pixel. Defaults to the cells.zarr mask "
            "transform when available, otherwise 1.0."
        ),
    )
    parser.add_argument(
        "--skip_validation",
        action="store_true",
        help="Skip consistency validation after writing each region.",
    )
    parser.add_argument(
        "--he_image",
        type=Path,
        default=None,
        help="Path to H&E OME-TIFF to crop alongside the Xenium region.",
    )
    parser.add_argument(
        "--he_alignment",
        type=Path,
        default=None,
        help="Path to alignment matrix CSV (3x3, H&E pixels -> Xenium pixels).",
    )
    args = parser.parse_args()

    if bool(args.bbox) == bool(args.regions_csv):
        parser.error("Provide exactly one of --bbox or --regions_csv.")
    if bool(args.he_image) != bool(args.he_alignment):
        parser.error("--he_image and --he_alignment must be provided together.")
    if args.proportion <= 0 or args.proportion > 1:
        parser.error("--proportion must be > 0 and <= 1.")
    return args


def load_regions(args):
    if args.bbox:
        xmin, ymin, xmax, ymax = args.bbox
        regions = [Region(args.region_name, xmin, ymin, xmax, ymax)]
    else:
        df = pd.read_csv(args.regions_csv)
        required = {"region", "xmin", "ymin", "xmax", "ymax"}
        missing = required - set(df.columns)
        if missing:
            raise ValueError(f"{args.regions_csv} missing columns: {sorted(missing)}")
        regions = [
            Region(str(r.region), float(r.xmin), float(r.ymin), float(r.xmax), float(r.ymax))
            for r in df.itertuples(index=False)
        ]

    for region in regions:
        if region.xmax <= region.xmin or region.ymax <= region.ymin:
            raise ValueError(f"Invalid bbox for {region.name}: {region}")
    return regions


def safe_name(name):
    return "".join(c if c.isalnum() or c in "._-" else "_" for c in name)


def bbox_mask(x, y, region):
    return (
        (x >= region.xmin)
        & (x < region.xmax)
        & (y >= region.ymin)
        & (y < region.ymax)
    )


def select_cells_in_region(input_dir, region, proportion, grid_size):
    rng = np.random.default_rng(42)

    cells_table = pq.read_table(
        input_dir / "cells.parquet", columns=["cell_id", "x_centroid", "y_centroid"]
    )
    n_total = len(cells_table)
    cell_ids = cells_table.column("cell_id").to_pylist()
    x = cells_table.column("x_centroid").to_numpy()
    y = cells_table.column("y_centroid").to_numpy()
    in_region = bbox_mask(x, y, region)

    if proportion < 1.0 and in_region.any():
        grid_x = ((x - region.xmin) / grid_size).astype(int)
        grid_y = ((y - region.ymin) / grid_size).astype(int)
        selected_mask = np.zeros(len(x), dtype=bool)
        grid_dict = {}
        for i in np.where(in_region)[0]:
            grid_dict.setdefault((int(grid_x[i]), int(grid_y[i])), []).append(i)
        for indices in grid_dict.values():
            n = max(1, round(len(indices) * proportion))
            selected_mask[rng.choice(indices, size=n, replace=False)] = True
    else:
        selected_mask = in_region

    selected_indices = np.where(selected_mask)[0]
    selected_ids = {cell_ids[i] for i in selected_indices}
    del cells_table
    gc.collect()
    return selected_indices, selected_ids, n_total


def shift_spatial_columns(df, region):
    x_columns = [
        "x_centroid",
        "x_location",
        "vertex_x",
        "x",
        "x0",
        "x1",
    ]
    y_columns = [
        "y_centroid",
        "y_location",
        "vertex_y",
        "y",
        "y0",
        "y1",
    ]
    for column in x_columns:
        if column in df.columns:
            df[column] = df[column] - region.xmin
    for column in y_columns:
        if column in df.columns:
            df[column] = df[column] - region.ymin
    return df


def subset_spatial_parquet_and_csv(input_dir, output_dir, filename, selected_ids, region):
    table = pq.read_table(input_dir / f"{filename}.parquet")
    mask = pc.is_in(table.column("cell_id"), value_set=pa.array(list(selected_ids)))
    sub = table.filter(mask).to_pandas()
    sub = shift_spatial_columns(sub, region)
    sub.to_parquet(output_dir / f"{filename}.parquet", index=False)
    sub.to_csv(output_dir / f"{filename}.csv.gz", index=False, compression="gzip")
    n = len(sub)
    del table, sub, mask
    gc.collect()
    return n


def process_transcripts_in_region(input_dir, output_dir, selected_ids, region, proportion):
    rng = np.random.default_rng(42)
    selected_ids_arrow = pa.array(list(selected_ids))
    parquet_file = pq.ParquetFile(input_dir / "transcripts.parquet")
    writer = None
    csv_path = output_dir / "transcripts.csv"
    first_chunk = True
    total_kept = 0

    for batch in parquet_file.iter_batches(batch_size=1_000_000):
        table = pa.Table.from_batches([batch])
        x = table.column("x_location")
        y = table.column("y_location")
        region_mask = pc.and_(
            pc.and_(pc.greater_equal(x, region.xmin), pc.less(x, region.xmax)),
            pc.and_(pc.greater_equal(y, region.ymin), pc.less(y, region.ymax)),
        )

        cell_id_col = table.column("cell_id")
        assigned_keep = pc.is_in(cell_id_col, value_set=selected_ids_arrow)
        unassigned = pc.equal(cell_id_col, "UNASSIGNED")
        if proportion < 1.0:
            sampled = pa.array(rng.random(len(table)) < proportion)
            unassigned = pc.and_(unassigned, sampled)

        keep = pc.and_(region_mask, pc.or_(assigned_keep, unassigned))
        chunk_sub = table.filter(keep)

        if len(chunk_sub) > 0:
            df = shift_spatial_columns(chunk_sub.to_pandas(), region)
            out_table = pa.Table.from_pandas(df, preserve_index=False)
            if writer is None:
                writer = pq.ParquetWriter(output_dir / "transcripts.parquet", out_table.schema)
            writer.write_table(out_table)
            df.to_csv(csv_path, mode="a", index=False, header=first_chunk)
            first_chunk = False
            total_kept += len(df)

        del table, chunk_sub
        gc.collect()

    if writer:
        writer.close()
    else:
        # Keep the expected files present even if the region has no transcripts.
        empty = parquet_file.schema_arrow.empty_table().to_pandas()
        empty.to_parquet(output_dir / "transcripts.parquet", index=False)
        empty.to_csv(csv_path, index=False)

    with open(csv_path, "rb") as f_in:
        with gzip.open(output_dir / "transcripts.csv.gz", "wb") as f_out:
            shutil.copyfileobj(f_in, f_out)
    csv_path.unlink()
    return total_kept


def infer_pixel_size(input_dir):
    def pixel_size_from_transform(transform):
        if transform.shape[0] >= 2 and transform.shape[1] >= 2:
            pixels_per_unit = float(np.mean(np.abs([transform[0, 0], transform[1, 1]])))
            if np.isfinite(pixels_per_unit) and pixels_per_unit > 0:
                return 1.0 / pixels_per_unit
        return None

    try:
        import numcodecs

        with zipfile.ZipFile(input_dir / "cells.zarr.zip") as zf:
            meta = json.loads(zf.read("masks/homogeneous_transform/.zarray"))
            raw = zf.read("masks/homogeneous_transform/0.0")
            compressor = meta.get("compressor")
            if compressor and compressor.get("id") == "blosc":
                raw = numcodecs.Blosc().decode(raw)
            transform = np.frombuffer(raw, dtype=np.dtype(meta["dtype"])).reshape(meta["shape"])
            pixel_size = pixel_size_from_transform(transform)
            if pixel_size:
                return pixel_size
    except Exception:
        pass

    try:
        store = open_zip_read_store(input_dir / "cells.zarr.zip", mode="r")
        z = zarr.open(store, mode="r")
        transform = np.asarray(z["masks/homogeneous_transform"][:], dtype=float)
        store.close()
        pixel_size = pixel_size_from_transform(transform)
        if pixel_size:
            return pixel_size
    except Exception:
        pass

    exp = input_dir / "experiment.xenium"
    if exp.exists():
        try:
            data = json.loads(exp.read_text())
            candidates = []

            def walk(obj):
                if isinstance(obj, dict):
                    for key, value in obj.items():
                        lowered = str(key).lower()
                        if "pixel" in lowered and isinstance(value, (int, float)):
                            candidates.append(float(value))
                        walk(value)
                elif isinstance(obj, list):
                    for value in obj:
                        walk(value)

            walk(data)
            candidates = [v for v in candidates if np.isfinite(v) and v > 0]
            if candidates:
                return candidates[0]
        except Exception:
            pass

    return 1.0


def crop_bounds_pixels(region, pixel_size):
    x0 = max(0, int(np.floor(region.xmin / pixel_size)))
    y0 = max(0, int(np.floor(region.ymin / pixel_size)))
    x1 = max(x0 + 1, int(np.ceil(region.xmax / pixel_size)))
    y1 = max(y0 + 1, int(np.ceil(region.ymax / pixel_size)))
    return x0, y0, x1, y1


def crop_array_last_yx(arr, x0, y0, x1, y1):
    key = [slice(None)] * arr.ndim
    key[-2] = slice(y0, min(y1, arr.shape[-2]))
    key[-1] = slice(x0, min(x1, arr.shape[-1]))
    return arr[tuple(key)]


def scaled_crop_bounds(base_shape, level_shape, x0, y0, x1, y1):
    base_y, base_x = base_shape[-2], base_shape[-1]
    level_y, level_x = level_shape[-2], level_shape[-1]
    scale_y = base_y / level_y
    scale_x = base_x / level_x
    lx0 = max(0, int(np.floor(x0 / scale_x)))
    ly0 = max(0, int(np.floor(y0 / scale_y)))
    lx1 = min(level_x, max(lx0 + 1, int(np.ceil(x1 / scale_x))))
    ly1 = min(level_y, max(ly0 + 1, int(np.ceil(y1 / scale_y))))
    return lx0, ly0, lx1, ly1


def crop_tiff_level(level, base_shape, x0, y0, x1, y1):
    lx0, ly0, lx1, ly1 = scaled_crop_bounds(base_shape, level.shape, x0, y0, x1, y1)
    crop_height = ly1 - ly0
    crop_width = lx1 - lx0
    pages = list(level.pages)

    if len(pages) == 1 and not pages[0].is_tiled:
        return crop_array_last_yx(level.asarray(), lx0, ly0, lx1, ly1)

    if not pages or not pages[0].is_tiled:
        return crop_array_last_yx(level.asarray(), lx0, ly0, lx1, ly1)

    sample_shape = (len(pages), crop_height, crop_width)
    cropped = np.zeros(sample_shape, dtype=pages[0].dtype)
    keyframe = pages[0].keyframe
    tile_width = keyframe.tilewidth
    tile_height = keyframe.tilelength

    for page_index, page in enumerate(pages):
        page_width = getattr(page, "imagewidth", page.shape[-1])
        page_height = getattr(page, "imagelength", page.shape[-2])
        tiles_x = int(np.ceil(page_width / tile_width))
        tile_x0 = lx0 // tile_width
        tile_x1 = (lx1 - 1) // tile_width
        tile_y0 = ly0 // tile_height
        tile_y1 = (ly1 - 1) // tile_height

        for tile_y in range(tile_y0, tile_y1 + 1):
            for tile_x in range(tile_x0, tile_x1 + 1):
                tile_index = tile_y * tiles_x + tile_x
                fh = page.parent.filehandle
                fh.open()
                fh.seek(page.dataoffsets[tile_index])
                encoded = fh.read(page.databytecounts[tile_index])
                tile, _, _ = page.decode(encoded, tile_index)
                tile = np.squeeze(tile)

                src_x0 = tile_x * tile_width
                src_y0 = tile_y * tile_height
                src_x1 = min(src_x0 + tile.shape[-1], page_width)
                src_y1 = min(src_y0 + tile.shape[-2], page_height)

                overlap_x0 = max(lx0, src_x0)
                overlap_y0 = max(ly0, src_y0)
                overlap_x1 = min(lx1, src_x1)
                overlap_y1 = min(ly1, src_y1)
                if overlap_x0 >= overlap_x1 or overlap_y0 >= overlap_y1:
                    continue

                out_x0 = overlap_x0 - lx0
                out_y0 = overlap_y0 - ly0
                out_x1 = overlap_x1 - lx0
                out_y1 = overlap_y1 - ly0
                in_x0 = overlap_x0 - src_x0
                in_y0 = overlap_y0 - src_y0
                in_x1 = overlap_x1 - src_x0
                in_y1 = overlap_y1 - src_y0

                cropped[page_index, out_y0:out_y1, out_x0:out_x1] = tile[
                    in_y0:in_y1, in_x0:in_x1
                ]

    if len(level.shape) == 2:
        return cropped[0]
    return cropped


def _read_tiff_level_crop(level, y0, x0, y1, x1):
    """Read a (y0:y1, x0:x1) crop from a tifffile series level.

    Handles tiled and non-tiled pages and both greyscale (YX) and RGB (YXS)
    images. Only tiles that overlap the crop are read.
    """
    pages = list(level.pages)
    if not pages:
        return np.zeros((y1 - y0, x1 - x0), dtype=np.uint8)

    page = pages[0]
    lh, lw = level.shape[0], level.shape[1]
    has_channels = level.ndim == 3
    n_ch = level.shape[2] if has_channels else 1

    if not page.is_tiled:
        arr = page.asarray()
        return arr[y0:y1, x0:x1]

    result_shape = (y1 - y0, x1 - x0, n_ch) if has_channels else (y1 - y0, x1 - x0)
    result = np.zeros(result_shape, dtype=page.dtype)

    keyframe = page.keyframe
    tile_h = keyframe.tilelength
    tile_w = keyframe.tilewidth
    tiles_x = int(np.ceil(lw / tile_w))

    for ty in range(y0 // tile_h, (y1 - 1) // tile_h + 1):
        for tx in range(x0 // tile_w, (x1 - 1) // tile_w + 1):
            tile_idx = ty * tiles_x + tx
            fh = page.parent.filehandle
            fh.open()
            fh.seek(page.dataoffsets[tile_idx])
            encoded = fh.read(page.databytecounts[tile_idx])
            tile, _, _ = page.decode(encoded, tile_idx)
            tile = np.squeeze(tile)

            src_y0 = ty * tile_h
            src_x0 = tx * tile_w
            src_y1 = min(src_y0 + tile.shape[0], lh)
            src_x1 = min(src_x0 + tile.shape[1], lw)

            ov_y0, ov_y1 = max(y0, src_y0), min(y1, src_y1)
            ov_x0, ov_x1 = max(x0, src_x0), min(x1, src_x1)
            if ov_y0 >= ov_y1 or ov_x0 >= ov_x1:
                continue

            result[ov_y0 - y0:ov_y1 - y0, ov_x0 - x0:ov_x1 - x0] = (
                tile[ov_y0 - src_y0:ov_y1 - src_y0, ov_x0 - src_x0:ov_x1 - src_x0]
            )

    return result


def crop_he_ome_tiff(src, dst, x0, y0, x1, y1):
    """Crop an RGB OME-TIFF (YXS axes) to the given pixel bounds.

    Uses (x0, y0, x1, y1) to match the convention of crop_ome_tiff.
    Each pyramid level is cropped proportionally; only overlapping tiles
    are read so the full base image is never loaded into memory.
    """
    import tifffile

    with tifffile.TiffFile(src) as tif:
        series = tif.series[0]
        levels = list(series.levels) if getattr(series, "levels", None) else [series]
        base_h, base_w = levels[0].shape[0], levels[0].shape[1]

        crops = []
        for level in levels:
            lh, lw = level.shape[0], level.shape[1]
            scale_y = base_h / lh
            scale_x = base_w / lw
            ly0 = max(0, int(np.floor(y0 / scale_y)))
            lx0 = max(0, int(np.floor(x0 / scale_x)))
            ly1 = min(lh, max(ly0 + 1, int(np.ceil(y1 / scale_y))))
            lx1 = min(lw, max(lx0 + 1, int(np.ceil(x1 / scale_x))))
            crops.append(_read_tiff_level_crop(level, ly0, lx0, ly1, lx1))

    with tifffile.TiffWriter(dst, bigtiff=True, ome=True) as tif_out:
        tif_out.write(
            crops[0],
            subifds=max(0, len(crops) - 1),
            photometric="rgb",
            compression="deflate",
            metadata={"axes": "YXS"},
        )
        for crop in crops[1:]:
            tif_out.write(
                crop,
                subfiletype=1,
                photometric="rgb",
                compression="deflate",
                metadata=None,
            )


def crop_ome_tiff(src, dst, x0, y0, x1, y1):
    import re
    import tifffile

    with tifffile.TiffFile(src) as tif:
        series = tif.series[0]
        axes = series.axes
        levels = list(series.levels) if getattr(series, "levels", None) else [series]
        base_shape = levels[0].shape
        crops = [
            crop_tiff_level(level, base_shape, x0, y0, x1, y1)
            for level in levels
        ]
        # Read channel names from OME metadata before closing the file handle.
        # tifffile generates a fresh OME header on write, so without this the
        # names are lost and spatialdata_io's v4 reader raises "channel without a name".
        channel_names = None
        if tif.is_ome:
            try:
                import ome_types
                ome = ome_types.from_xml(tifffile.tiffcomment(str(src)), validate=False)
                names = [ch.name for ch in ome.images[0].pixels.channels]
                if any(n is not None for n in names):
                    channel_names = names
            except Exception:
                pass

    metadata = {"axes": axes} if axes and len(axes) == crops[0].ndim else None
    if metadata is not None and channel_names:
        metadata["Channel"] = {"Name": channel_names}
    with tifffile.TiffWriter(dst, bigtiff=True, ome=True) as tif:
        tif.write(
            crops[0],
            subifds=max(0, len(crops) - 1),
            photometric="minisblack",
            compression="deflate",
            metadata=metadata,
        )
        for crop in crops[1:]:
            tif.write(
                crop,
                subfiletype=1,
                photometric="minisblack",
                compression="deflate",
                metadata=None,
            )

    # tifffile auto-generates channel IDs as "Channel:<image_index>:<channel_index>"
    # (standard OME spec), but spatialdata_io's v4 xenium reader expects Xenium's
    # native format "Channel:<channel_index>" (no image prefix). Patch the OME-XML
    # comment in-place after writing so the output is readable without workarounds.
    ome_xml = tifffile.tiffcomment(str(dst))
    ome_xml_fixed = re.sub(r'Channel:(\d+):(\d+)', r'Channel:\2', ome_xml)
    if ome_xml_fixed != ome_xml:
        tifffile.tiffcomment(str(dst), ome_xml_fixed)


def load_he_alignment(alignment_path):
    """Read a 3x3 affine matrix (H&E pixels -> Xenium pixels) from a CSV."""
    rows = []
    with open(alignment_path) as f:
        for line in f:
            line = line.strip()
            if line:
                rows.append([float(v) for v in line.split(",")])
    return np.array(rows)


def he_crop_bounds(region, alignment_matrix, pixel_size):
    """Return H&E pixel (x0, y0, x1, y1) for a Xenium coordinate bbox.

    alignment_matrix maps H&E pixels -> Xenium pixels; we invert it to get
    Xenium pixels -> H&E pixels. Region coords are in µm so we first divide
    by pixel_size to convert to Xenium pixels before applying the inverse.
    """
    xenium_to_he = np.linalg.inv(alignment_matrix)
    corners = np.array([
        [region.xmin / pixel_size, region.ymin / pixel_size, 1],
        [region.xmax / pixel_size, region.ymin / pixel_size, 1],
        [region.xmin / pixel_size, region.ymax / pixel_size, 1],
        [region.xmax / pixel_size, region.ymax / pixel_size, 1],
    ])
    transformed = (xenium_to_he @ corners.T).T
    px = transformed[:, 0]
    py = transformed[:, 1]
    x0 = max(0, int(np.floor(px.min())))
    y0 = max(0, int(np.floor(py.min())))
    x1 = max(x0 + 1, int(np.ceil(px.max())))
    y1 = max(y0 + 1, int(np.ceil(py.max())))
    return x0, y0, x1, y1


def copy_and_crop_images(input_dir, output_dir, region, pixel_size, he_image=None, he_alignment=None):
    copy_files = [
        "experiment.xenium",
        "gene_panel.json",
        "metrics_summary.csv",
        "analysis_summary.html",
    ]
    for fname in copy_files:
        src = input_dir / fname
        if src.exists():
            shutil.copy2(src, output_dir / fname)

    x0, y0, x1, y1 = crop_bounds_pixels(region, pixel_size)
    src = input_dir / "morphology.ome.tif"
    if src.exists():
        print(f"    morphology.ome.tif crop pixels x={x0}:{x1}, y={y0}:{y1}...")
        crop_ome_tiff(src, output_dir / "morphology.ome.tif", x0, y0, x1, y1)

    focus_in = input_dir / "morphology_focus"
    if focus_in.exists():
        focus_out = output_dir / "morphology_focus"
        focus_out.mkdir(exist_ok=True)
        for image in sorted(focus_in.glob("*.ome.tif")):
            print(f"    morphology_focus/{image.name} crop...")
            crop_ome_tiff(image, focus_out / image.name, x0, y0, x1, y1)

    if he_image is not None and he_alignment is not None:
        hx0, hy0, hx1, hy1 = he_crop_bounds(region, he_alignment, pixel_size)
        print(f"    he_image.ome.tif crop pixels x={hx0}:{hx1}, y={hy0}:{hy1}...")
        crop_he_ome_tiff(he_image, output_dir / "he_image.ome.tif", hx0, hy0, hx1, hy1)

        # Build the alignment matrix for the cropped H&E in spatialdata global space
        # (Xenium pixel units). The input matrix maps H&E pixels → Xenium pixels; we
        # adjust only the translation to account for the H&E crop origin (hx0, hy0)
        # and the rebased Xenium origin (region.xmin / pixel_size, region.ymin / pixel_size).
        # No row scaling is needed — the matrix already outputs Xenium pixels.
        cropped_matrix = he_alignment.copy()
        cropped_matrix[0, 2] = (
            he_alignment[0, 0] * hx0 + he_alignment[0, 1] * hy0 + he_alignment[0, 2]
            - region.xmin / pixel_size
        )
        cropped_matrix[1, 2] = (
            he_alignment[1, 0] * hx0 + he_alignment[1, 1] * hy0 + he_alignment[1, 2]
            - region.ymin / pixel_size
        )
        cropped_matrix[2] = [0.0, 0.0, 1.0]
        # spatialdata_io.xenium() auto-detects alignment files named
        # <image_stem>alignment.csv, so he_image.ome.tif → he_imagealignment.csv.
        alignment_out = output_dir / "he_imagealignment.csv"
        with open(alignment_out, "w") as f:
            for row in cropped_matrix:
                f.write(",".join(f"{v}" for v in row) + "\n")
        print(f"    he_imagealignment.csv written")


def shift_zarr_dataset_if_spatial(key, data, region):
    lowered = key.lower()
    if data.ndim == 1 and ("x" in lowered or "y" in lowered):
        if lowered in {"x", "vertex_x", "vertices_x", "x_location", "x_centroid"}:
            return data - region.xmin
        if lowered in {"y", "vertex_y", "vertices_y", "y_location", "y_centroid"}:
            return data - region.ymin
    if data.ndim == 2 and data.shape[1] >= 2 and lowered in {"vertices", "points"}:
        data = data.copy()
        data[:, 0] -= region.xmin
        data[:, 1] -= region.ymin
    return data


def process_cells_zarr_region(input_dir, output_dir, selected_indices, n_total, region, pixel_size):
    store_in = open_zip_read_store(input_dir / "cells.zarr.zip", mode="r")
    z_in = zarr.open(store_in, mode="r")
    x0, y0, x1, y1 = crop_bounds_pixels(region, pixel_size)

    with tempfile.TemporaryDirectory() as tmpdir:
        tmpdir_path = Path(tmpdir)
        z_out = open_directory_group(tmpdir, mode="w")

        create_zarr_dataset(z_out, "cell_id", data=z_in["cell_id"][:][selected_indices])
        cell_summary = z_in["cell_summary"][:][selected_indices]
        summary_attrs = z_in["cell_summary"].attrs.asdict()
        summary_columns = (
            summary_attrs.get("columns")
            or summary_attrs.get("column_names")
            or summary_attrs.get("_ARRAY_DIMENSIONS")
            or []
        )
        if cell_summary.ndim == 2 and len(summary_columns) == cell_summary.shape[1]:
            cell_summary = cell_summary.copy()
            for idx, column in enumerate(summary_columns):
                if column in {"x_centroid", "x"}:
                    cell_summary[:, idx] -= region.xmin
                if column in {"y_centroid", "y"}:
                    cell_summary[:, idx] -= region.ymin
        create_zarr_dataset(z_out, "cell_summary", data=cell_summary)
        z_out["cell_summary"].attrs.update(summary_attrs)

        selected_label_ids = selected_indices + 1
        lookup = np.zeros(n_total + 1, dtype=bool)
        lookup[selected_label_ids] = True

        masks_out = z_out.create_group("masks")
        transform = np.asarray(z_in["masks/homogeneous_transform"][:], dtype=float)
        if transform.shape[0] >= 2 and transform.shape[1] >= 3:
            transform = transform.copy()
            transform[0, 2] = transform[0, 2] + x0 * transform[0, 0] - region.xmin
            transform[1, 2] = transform[1, 2] + y0 * transform[1, 1] - region.ymin
        create_zarr_dataset(masks_out, "homogeneous_transform", data=transform)

        for mask_name in ["0", "1"]:
            print(f"    mask {mask_name} crop...")
            mask_in = z_in[f"masks/{mask_name}"]
            crop = np.asarray(mask_in[y0:min(y1, mask_in.shape[0]), x0:min(x1, mask_in.shape[1])])
            nonzero = crop > 0
            if nonzero.any():
                in_range = crop <= n_total
                keep = np.zeros_like(crop, dtype=bool)
                valid = nonzero & in_range
                keep[valid] = lookup[crop[valid]]
                crop[nonzero & ~keep] = 0
            create_zarr_dataset(masks_out, mask_name, data=crop, chunks=mask_in.chunks)

        old_to_new_idx = {old: new for new, old in enumerate(selected_indices)}
        for ps_name in ["0", "1"]:
            ps_in = z_in[f"polygon_sets/{ps_name}"]
            cell_index = ps_in["cell_index"][:]
            keep_mask = np.isin(cell_index, selected_indices)
            ps_out = z_out.create_group(f"polygon_sets/{ps_name}")
            for key in ps_in.keys():
                data = ps_in[key][:][keep_mask]
                if key == "cell_index":
                    data = np.array([old_to_new_idx[i] for i in data], dtype=data.dtype)
                else:
                    data = shift_zarr_dataset_if_spatial(key, data, region)
                create_zarr_dataset(ps_out, key, data=data)
            if ps_in.attrs:
                ps_out.attrs.update(ps_in.attrs.asdict())

        close_zarr_group(z_out)
        archive_directory_as_zip(tmpdir_path, output_dir / "cells.zarr.zip")

    store_in.close()


def run_region(input_dir, output_dir, region, proportion, grid_size, pixel_size, skip_validation, he_image=None, he_alignment=None):
    output_dir.mkdir(parents=True, exist_ok=True)
    print(f"\nRegion:     {region.name}")
    print(f"Output:     {output_dir}")
    print(f"BBox:       {region.xmin}, {region.ymin}, {region.xmax}, {region.ymax}")
    print(f"Pixel size: {pixel_size}")

    print("Step 1: Selecting cells in region...")
    selected_indices, selected_ids, n_total = select_cells_in_region(
        input_dir, region, proportion, grid_size
    )
    print(
        f"  Selected {len(selected_ids):,} / {n_total:,} cells "
        f"({len(selected_ids) / n_total * 100:.1f}%)"
    )
    if not selected_ids:
        raise ValueError(f"No cells selected for region {region.name}")

    print("\nStep 2: Subsetting and rebasing tabular cell files...")
    for name in ["cells", "cell_boundaries", "nucleus_boundaries"]:
        print(f"  {name}...")
        subset_spatial_parquet_and_csv(input_dir, output_dir, name, selected_ids, region)

    print("\nStep 3: Subsetting and rebasing transcripts...")
    n_transcripts = process_transcripts_in_region(
        input_dir, output_dir, selected_ids, region, proportion
    )
    print(f"  Kept {n_transcripts:,} transcripts")

    print("\nStep 4: Subsetting cell feature matrix...")
    print("  cell_feature_matrix/ directory...")
    process_cell_feature_matrix_dir(input_dir, output_dir, selected_ids)
    print("  cell_feature_matrix.h5...")
    process_cell_feature_matrix_h5(input_dir, output_dir, selected_ids)
    print("  cell_feature_matrix.tar.gz...")
    process_cell_feature_matrix_tar(output_dir)

    print("\nStep 5: Subsetting analysis results...")
    process_analysis(input_dir, output_dir, selected_ids)

    print("\nStep 6: Cropping cells.zarr.zip...")
    process_cells_zarr_region(input_dir, output_dir, selected_indices, n_total, region, pixel_size)

    print("\nStep 7: Subsetting analysis.zarr.zip...")
    process_analysis_zarr(input_dir, output_dir, selected_ids)

    print("\nStep 8: Subsetting cell_feature_matrix.zarr.zip...")
    process_cfm_zarr(input_dir, output_dir, selected_ids)

    print("\nStep 9: Copying metadata and cropping images...")
    copy_and_crop_images(input_dir, output_dir, region, pixel_size, he_image=he_image, he_alignment=he_alignment)

    if not skip_validation:
        print("\nStep 10: Validating output...")
        if not validate_output(output_dir):
            raise RuntimeError(f"Validation failed for {region.name}")


def main():
    args = parse_args()
    input_dir = args.input_dir.resolve()
    output_root = (
        args.output_dir.resolve()
        if args.output_dir
        else input_dir.parent / f"{input_dir.name}_region_downsampled"
    )
    regions = load_regions(args)
    pixel_size = args.pixel_size if args.pixel_size else infer_pixel_size(input_dir)
    he_image = args.he_image.resolve() if args.he_image else None
    he_alignment = load_he_alignment(args.he_alignment) if args.he_alignment else None

    print(f"Input:      {input_dir}")
    print(f"Output root:{output_root}")
    print(f"Regions:    {len(regions)}")
    print(f"Proportion: {args.proportion}")
    print(f"Grid size:  {args.grid_size}")
    if he_image:
        print(f"H&E image:  {he_image}")

    t0 = time.time()
    try:
        for region in regions:
            run_region(
                input_dir=input_dir,
                output_dir=output_root / safe_name(region.name),
                region=region,
                proportion=args.proportion,
                grid_size=args.grid_size,
                pixel_size=pixel_size,
                skip_validation=args.skip_validation,
                he_image=he_image,
                he_alignment=he_alignment,
            )
    except Exception as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        sys.exit(1)

    elapsed = time.time() - t0
    print(f"\nDone in {elapsed / 60:.1f} minutes.")
    print(f"Output root: {output_root}")


if __name__ == "__main__":
    main()
