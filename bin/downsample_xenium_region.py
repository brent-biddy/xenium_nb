#!/usr/bin/env python3
"""
downsample_xenium_region.py - Crop a Xenium output directory to one or more regions.

Selects cells inside a bounding box, keeps transcripts inside the same box,
rebases spatial coordinates to the crop origin, and crops morphology images
and cells.zarr masks so the output is physically smaller.

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
import shutil
import sys
import tarfile
import tempfile
import time
from dataclasses import dataclass
from pathlib import Path

import h5py
import numpy as np
import pandas as pd
import pyarrow as pa
import pyarrow.compute as pc
import pyarrow.parquet as pq
import zarr
from scipy.io import mmread, mmwrite
from scipy.sparse import csc_matrix

# Xenium pixel size in µm/pixel — a fixed instrument constant used to convert
# bounding box coordinates from µm (user-facing) to pixels (image space).
XENIUM_PIXEL_SIZE_UM = 0.2125


# ---------------------------------------------------------------------------
# Zarr utilities
# ---------------------------------------------------------------------------

def open_zip_read_store(path, mode):
    return zarr.storage.ZipStore(str(path), mode=mode)


def create_zarr_dataset(group, name, data=None, **kwargs):
    if data is not None:
        kwargs.setdefault("shape", np.shape(data))
        kwargs.setdefault("dtype", getattr(data, "dtype", np.asarray(data).dtype))
        return group.create_dataset(name, data=data, **kwargs)
    return group.create_dataset(name, **kwargs)


def open_directory_group(path, mode):
    open_group = getattr(zarr, "open_group", None)
    if open_group is not None:
        return open_group(store=str(path), mode=mode)
    return zarr.open(str(path), mode=mode)


def archive_directory_as_zip(source_dir, output_zip_path):
    archive_base = output_zip_path.with_suffix("")
    shutil.make_archive(str(archive_base), "zip", root_dir=str(source_dir))


def close_zarr_group(group):
    store = getattr(group, "store", None)
    close = getattr(store, "close", None)
    if callable(close):
        close()


# Region is immutable so it can be freely passed between functions without risk
# of accidental mutation.
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
        help="Coordinate units per image pixel. Defaults to XENIUM_PIXEL_SIZE_UM (0.2125).",
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
    """Parse --bbox or --regions_csv into a list of Region objects."""
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
    """Sanitize a region name for use as a filesystem directory name."""
    return "".join(c if c.isalnum() or c in "._-" else "_" for c in name)


# ---------------------------------------------------------------------------
# Cell selection
# ---------------------------------------------------------------------------

def bbox_mask(x, y, region):
    """Return a boolean array that is True for points inside the region bbox."""
    return (
        (x >= region.xmin)
        & (x < region.xmax)
        & (y >= region.ymin)
        & (y < region.ymax)
    )


def select_cells_in_region(input_dir, region, proportion, grid_size):
    """Return (selected_indices, selected_ids, n_total) for cells inside region.

    When proportion < 1, cells are spatially downsampled: the region is divided
    into a grid and `proportion` of cells are randomly kept per grid cell. This
    preserves spatial coverage better than random global subsampling.
    Seed is fixed (42) for reproducibility.
    """
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
        # Divide the region into grid cells and sample `proportion` of cells
        # from each grid cell independently to preserve spatial coverage.
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


# ---------------------------------------------------------------------------
# Tabular files
# ---------------------------------------------------------------------------

def shift_spatial_columns(df, region):
    """Rebase spatial coordinates in a DataFrame to the crop origin.

    Subtracts region.xmin/ymin from all recognised spatial column names so
    that coordinates in the output are relative to the crop's top-left corner
    rather than the full slide origin.
    """
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


def subset_spatial_parquet_and_csv(input_dir, output_dir, filename, selected_ids, region,
                                    label_id_remap=None):
    """Filter a parquet file to selected_ids, rebase coordinates, write parquet + csv.gz.

    label_id_remap: when provided and the file has a 'label_id' column, renumber
    label_id values according to the dict (old_label -> new_label). Used to keep
    boundary parquets consistent with the cropped cells.zarr.zip polygon_sets, which
    spatialdata-io 0.7.0 requires to be sequential 1..M.
    """
    table = pq.read_table(input_dir / f"{filename}.parquet")
    cell_mask = pc.is_in(table.column("cell_id"), value_set=pa.array(list(selected_ids)))
    sub = table.filter(cell_mask).to_pandas()
    if label_id_remap is not None and "label_id" in sub.columns:
        sub["label_id"] = sub["label_id"].map(label_id_remap).astype(sub["label_id"].dtype)
    sub = shift_spatial_columns(sub, region)
    sub.to_parquet(output_dir / f"{filename}.parquet", index=False)
    sub.to_csv(output_dir / f"{filename}.csv.gz", index=False, compression="gzip")
    n = len(sub)
    del table, sub
    gc.collect()
    return n


def process_transcripts_in_region(input_dir, output_dir, selected_ids, region, proportion):
    """Filter transcripts.parquet to the region, writing parquet + csv.gz.

    Reads in 1M-row batches to avoid loading the full transcript table into
    memory (whole-slide transcript files can exceed several GB).

    Assigned transcripts (cell_id != UNASSIGNED) are kept only if their cell
    is in selected_ids. Unassigned transcripts are kept if they fall inside
    the region bbox; if proportion < 1 they are additionally randomly thinned
    to match the cell downsampling rate.
    """
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


# ---------------------------------------------------------------------------
# Image cropping
# ---------------------------------------------------------------------------

def crop_bounds_pixels(region, pixel_size):
    """Convert a µm bounding box to pixel indices, flooring/ceiling to include all partial pixels."""
    x0 = max(0, int(np.floor(region.xmin / pixel_size)))
    y0 = max(0, int(np.floor(region.ymin / pixel_size)))
    x1 = max(x0 + 1, int(np.ceil(region.xmax / pixel_size)))
    y1 = max(y0 + 1, int(np.ceil(region.ymax / pixel_size)))
    return x0, y0, x1, y1


def crop_array_last_yx(arr, x0, y0, x1, y1):
    """Crop the last two axes (y, x) of an array, clamping to array bounds."""
    key = [slice(None)] * arr.ndim
    key[-2] = slice(y0, min(y1, arr.shape[-2]))
    key[-1] = slice(x0, min(x1, arr.shape[-1]))
    return arr[tuple(key)]


def scaled_crop_bounds(base_shape, level_shape, x0, y0, x1, y1):
    """Scale base-level pixel crop bounds to a lower pyramid level.

    OME-TIFF pyramids store each level at a reduced resolution. To crop the
    same physical region at every level, the base-level pixel bounds must be
    scaled by the ratio between base and level dimensions.
    """
    base_y, base_x = base_shape[-2], base_shape[-1]
    level_y, level_x = level_shape[-2], level_shape[-1]
    scale_y = base_y / level_y
    scale_x = base_x / level_x
    lx0 = max(0, int(np.floor(x0 / scale_x)))
    ly0 = max(0, int(np.floor(y0 / scale_y)))
    lx1 = min(level_x, max(lx0 + 1, int(np.ceil(x1 / scale_x))))
    ly1 = min(level_y, max(ly0 + 1, int(np.ceil(y1 / scale_y))))
    return lx0, ly0, lx1, ly1


def tile_crop_slices(crop_y0, crop_x0, crop_y1, crop_x1, src_y0, src_x0, src_y1, src_x1):
    """Compute the (dst, src) slice pair for a tile's overlap with a crop window.

    All bounds are in level pixel coordinates: the crop window (crop_*) and the
    tile's actual extent (src_*). Returns ((out_y, out_x), (in_y, in_x)) where
    out slices are relative to the crop origin and in slices to the tile origin,
    so `output[out_y, out_x] = tile[in_y, in_x]` copies the overlap. Returns None
    if the tile does not intersect the crop window.
    """
    oy0 = max(crop_y0, src_y0)
    ox0 = max(crop_x0, src_x0)
    oy1 = min(crop_y1, src_y1)
    ox1 = min(crop_x1, src_x1)
    if oy0 >= oy1 or ox0 >= ox1:
        return None
    dst = (slice(oy0 - crop_y0, oy1 - crop_y0), slice(ox0 - crop_x0, ox1 - crop_x0))
    src = (slice(oy0 - src_y0, oy1 - src_y0), slice(ox0 - src_x0, ox1 - src_x0))
    return dst, src


def crop_tiff_level(level, base_shape, x0, y0, x1, y1):
    """Crop one pyramid level of a Xenium morphology OME-TIFF (channels as pages, CYX layout)."""
    lx0, ly0, lx1, ly1 = scaled_crop_bounds(base_shape, level.shape, x0, y0, x1, y1)
    crop_height = ly1 - ly0
    crop_width = lx1 - lx0
    pages = list(level.pages)

    if not pages or not pages[0].is_tiled:
        return crop_array_last_yx(level.asarray(), lx0, ly0, lx1, ly1)

    # Allocate output array over all pages (one page per channel).
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

                # Copy only the overlap between this tile and the crop window.
                slices = tile_crop_slices(ly0, lx0, ly1, lx1, src_y0, src_x0, src_y1, src_x1)
                if slices is None:
                    continue
                (out_y, out_x), (in_y, in_x) = slices
                cropped[page_index, out_y, out_x] = tile[in_y, in_x]

    if len(level.shape) == 2:
        return cropped[0]
    return cropped


def _read_tiff_level_crop(level, y0, x0, y1, x1):
    """Crop one pyramid level of an H&E OME-TIFF (channels in last axis, YXS/RGB layout).

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

            slices = tile_crop_slices(y0, x0, y1, x1, src_y0, src_x0, src_y1, src_x1)
            if slices is None:
                continue
            (out_y, out_x), (in_y, in_x) = slices
            result[out_y, out_x] = tile[in_y, in_x]

    return result


def write_ome_pyramid(dst, crops, photometric, base_metadata):
    """Write pyramid-level crops to a BigTIFF OME-TIFF.

    crops[0] is the full-resolution base level; the rest are stored as
    reduced-resolution SubIFDs. Only the base level carries OME metadata.
    """
    import tifffile

    with tifffile.TiffWriter(dst, bigtiff=True, ome=True) as tif:
        tif.write(
            crops[0],
            subifds=max(0, len(crops) - 1),
            photometric=photometric,
            compression="deflate",
            metadata=base_metadata,
        )
        for crop in crops[1:]:
            tif.write(
                crop,
                subfiletype=1,
                photometric=photometric,
                compression="deflate",
                metadata=None,
            )


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
        # H&E is YXS (spatial axes first), so scale on shape[:2] to match
        # scaled_crop_bounds' last-two-axes (y, x) convention.
        base_yx = levels[0].shape[:2]

        crops = []
        for level in levels:
            lx0, ly0, lx1, ly1 = scaled_crop_bounds(base_yx, level.shape[:2], x0, y0, x1, y1)
            crops.append(_read_tiff_level_crop(level, ly0, lx0, ly1, lx1))

    write_ome_pyramid(dst, crops, photometric="rgb", base_metadata={"axes": "YXS"})


def crop_ome_tiff(src, dst, x0, y0, x1, y1):
    """Crop a Xenium morphology OME-TIFF to the given pixel bounds, preserving channel names."""
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
    write_ome_pyramid(dst, crops, photometric="minisblack", base_metadata=metadata)

    # tifffile auto-generates channel IDs as "Channel:<image_index>:<channel_index>"
    # (standard OME spec), but spatialdata_io's v4 xenium reader expects Xenium's
    # native format "Channel:<channel_index>" (no image prefix). Patch the OME-XML
    # comment in-place after writing so the output is readable without workarounds.
    ome_xml = tifffile.tiffcomment(str(dst))
    ome_xml_fixed = re.sub(r'Channel:(\d+):(\d+)', r'Channel:\2', ome_xml)
    if ome_xml_fixed != ome_xml:
        tifffile.tiffcomment(str(dst), ome_xml_fixed)


# ---------------------------------------------------------------------------
# H&E alignment
# ---------------------------------------------------------------------------

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
    All four corners of the bbox are transformed to handle non-axis-aligned
    rotations in the alignment, then we take the bounding box of the result.
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


# ---------------------------------------------------------------------------
# Zarr processing
# ---------------------------------------------------------------------------

def shift_zarr_dataset_if_spatial(key, data, region):
    """Rebase spatial coordinates in a zarr dataset array to the crop origin."""
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
    """Crop cells.zarr.zip to the selected cells and region pixel bounds.

    Returns (nucleus_label_id_remap, cell_label_id_remap): dicts mapping
    {old_label_id → new_label_id} for renumbering boundary parquets.

    Xenium native invariant: polygon at position p in polygon_sets/{0,1} has
    raster label p+1 in masks/{0,1}. After filtering to M kept positions, the
    new polygon at rank r has label r+1. The remap dict carries this mapping
    so nucleus_boundaries.parquet and cell_boundaries.parquet can be updated
    to match the new sequential labels that spatialdata-io 0.7.0 requires.

    Masks are read in row-chunks to avoid loading the full slide into memory.
    """
    store_in = open_zip_read_store(input_dir / "cells.zarr.zip", mode="r")
    z_in = zarr.open(store_in, mode="r")
    x0, y0, x1, y1 = crop_bounds_pixels(region, pixel_size)

    with tempfile.TemporaryDirectory() as tmpdir:
        tmpdir_path = Path(tmpdir)
        z_out = open_directory_group(tmpdir, mode="w")

        create_zarr_dataset(z_out, "cell_id", data=z_in["cell_id"][:][selected_indices])
        cell_summary = z_in["cell_summary"][:][selected_indices]
        summary_attrs = z_in["cell_summary"].attrs.asdict()
        # column_names may be stored under different attribute keys across Xenium versions.
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

        masks_out = z_out.create_group("masks")
        transform = np.asarray(z_in["masks/homogeneous_transform"][:], dtype=float)
        if transform.shape[0] >= 2 and transform.shape[1] >= 3:
            # Adjust the translation column of the transform to account for the
            # pixel crop offset (x0, y0) and the rebased coordinate origin.
            transform = transform.copy()
            transform[0, 2] = transform[0, 2] + x0 * transform[0, 0] - region.xmin
            transform[1, 2] = transform[1, 2] + y0 * transform[1, 1] - region.ymin
        create_zarr_dataset(masks_out, "homogeneous_transform", data=transform)

        old_to_new_idx = {old: new for new, old in enumerate(selected_indices)}
        remaps = {}

        for ps_name in ["0", "1"]:
            ps_in = z_in[f"polygon_sets/{ps_name}"]
            cell_index = ps_in["cell_index"][:]
            keep_mask = np.isin(cell_index, selected_indices)
            kept_positions = np.where(keep_mask)[0]
            n_ps_total = len(cell_index)

            # Build a lookup table (LUT) that maps old raster label → new raster label.
            # Position p in polygon_sets has raster label p+1 (Xenium native invariant).
            # After filtering, position new_rank has label new_rank+1.
            # LUT[old_pos+1] = new_rank+1; LUT[anything else] = 0 (erase pixel).
            lut = np.zeros(n_ps_total + 1, dtype=np.uint32)
            for new_rank, old_pos in enumerate(kept_positions):
                lut[int(old_pos) + 1] = new_rank + 1
            remaps[ps_name] = {
                int(old_pos) + 1: int(new_rank) + 1
                for new_rank, old_pos in enumerate(kept_positions)
            }

            # Crop and remap the mask in row-chunks to avoid loading the full
            # slide image into memory (can be multiple GB for large FOVs).
            print(f"    mask {ps_name} crop + remap...")
            mask_in = z_in[f"masks/{ps_name}"]
            crop_y0 = y0
            crop_y1 = min(y1, mask_in.shape[0])
            crop_x0 = x0
            crop_x1 = min(x1, mask_in.shape[1])
            out_h = max(0, crop_y1 - crop_y0)
            out_w = max(0, crop_x1 - crop_x0)

            chunk_rows = mask_in.chunks[0] if mask_in.chunks else 512
            mask_arr_out = masks_out.create_dataset(
                ps_name,
                shape=(out_h, out_w),
                dtype=np.uint32,
                chunks=(min(chunk_rows, out_h), out_w) if out_h > 0 else (1, 1),
            )
            for start in range(0, out_h, chunk_rows):
                end = min(start + chunk_rows, out_h)
                chunk = np.asarray(
                    mask_in[crop_y0 + start:crop_y0 + end, crop_x0:crop_x1]
                ).astype(np.uint32)
                # LUT maps kept labels to new sequential values; out-of-range → 0.
                chunk_clipped = np.minimum(chunk, n_ps_total)
                mask_arr_out[start:end, :] = lut[chunk_clipped]

            # Write filtered and rebased polygon_sets.
            print(f"    polygon_sets/{ps_name}...")
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
    # mask "0" ↔ polygon_sets "0" ↔ nucleus_boundaries; "1" ↔ cell_boundaries
    return remaps["0"], remaps["1"]


# ---------------------------------------------------------------------------
# Cell feature matrix
# ---------------------------------------------------------------------------

def subset_csc_columns(data, indices, indptr, col_indices):
    """Subset a CSC sparse matrix to the given columns (cells).

    Slices each selected column's value/row-index range out of data/indices and
    rebuilds the pointer array. Returns (data, indices, indptr) as numpy arrays;
    an empty selection yields correctly-typed empty arrays.
    """
    new_data, new_indices, new_indptr = [], [], [0]
    for col_idx in col_indices:
        s, e = indptr[col_idx], indptr[col_idx + 1]
        new_data.append(data[s:e])
        new_indices.append(indices[s:e])
        new_indptr.append(new_indptr[-1] + (e - s))
    out_data = np.concatenate(new_data) if new_data else np.array([], dtype=data.dtype)
    out_indices = (
        np.concatenate(new_indices) if new_indices else np.array([], dtype=indices.dtype)
    )
    out_indptr = np.array(new_indptr, dtype=indptr.dtype)
    return out_data, out_indices, out_indptr


def process_cell_feature_matrix_dir(input_dir, output_dir, selected_ids):
    """Subset cell_feature_matrix/ directory (barcodes, matrix, features)."""
    cfm_in = input_dir / "cell_feature_matrix"
    cfm_out = output_dir / "cell_feature_matrix"
    cfm_out.mkdir(exist_ok=True)

    # Barcodes — find selected column indices
    barcodes = pd.read_csv(
        cfm_in / "barcodes.tsv.gz", header=None, names=["barcode"]
    )
    barcode_list = barcodes["barcode"].tolist()
    selected_col_indices = [i for i, b in enumerate(barcode_list) if b in selected_ids]
    selected_barcodes = [barcode_list[i] for i in selected_col_indices]

    with gzip.open(cfm_out / "barcodes.tsv.gz", "wt") as f:
        f.write("\n".join(selected_barcodes) + "\n")

    # Features — copy unchanged
    shutil.copy(cfm_in / "features.tsv.gz", cfm_out / "features.tsv.gz")

    # Matrix — subset columns
    matrix = mmread(cfm_in / "matrix.mtx.gz")
    matrix_csc = csc_matrix(matrix)
    del matrix
    gc.collect()

    matrix_sub = matrix_csc[:, selected_col_indices]
    del matrix_csc
    gc.collect()

    with tempfile.NamedTemporaryFile(suffix=".mtx", delete=False) as tmp:
        tmp_path = Path(tmp.name)
    mmwrite(str(tmp_path), matrix_sub)
    del matrix_sub
    gc.collect()

    with open(tmp_path, "rb") as f_in:
        with gzip.open(cfm_out / "matrix.mtx.gz", "wb") as f_out:
            shutil.copyfileobj(f_in, f_out)
    tmp_path.unlink()

    return selected_col_indices


def process_cell_feature_matrix_h5(input_dir, output_dir, selected_ids):
    """Subset cell_feature_matrix.h5 (CSC sparse matrix)."""
    with h5py.File(input_dir / "cell_feature_matrix.h5", "r") as f_in:
        barcodes = f_in["matrix/barcodes"][:]
        barcode_list = [b.decode() for b in barcodes]
        selected_col_indices = [
            i for i, b in enumerate(barcode_list) if b in selected_ids
        ]

        data = f_in["matrix/data"][:]
        indices = f_in["matrix/indices"][:]
        indptr = f_in["matrix/indptr"][:]
        shape = f_in["matrix/shape"][:]

        new_data, new_indices, new_indptr = subset_csc_columns(
            data, indices, indptr, selected_col_indices
        )
        new_shape = np.array([shape[0], len(selected_col_indices)], dtype=shape.dtype)
        new_barcodes = np.array([barcodes[i] for i in selected_col_indices])

        with h5py.File(output_dir / "cell_feature_matrix.h5", "w") as f_out:
            f_out.create_dataset("matrix/barcodes", data=new_barcodes)
            f_out.create_dataset("matrix/data", data=new_data)
            f_out.create_dataset("matrix/indices", data=new_indices)
            f_out.create_dataset("matrix/indptr", data=new_indptr)
            f_out.create_dataset("matrix/shape", data=new_shape)
            for key in f_in["matrix/features"].keys():
                f_out.create_dataset(
                    f"matrix/features/{key}",
                    data=f_in[f"matrix/features/{key}"][:],
                )

    del data, indices, indptr, new_data, new_indices
    gc.collect()


def process_cell_feature_matrix_tar(output_dir):
    """Create cell_feature_matrix.tar.gz from the subset directory."""
    cfm_out = output_dir / "cell_feature_matrix"
    with tarfile.open(output_dir / "cell_feature_matrix.tar.gz", "w:gz") as tar:
        for f in sorted(cfm_out.iterdir()):
            tar.add(str(f), arcname=f"cell_feature_matrix/{f.name}")


# ---------------------------------------------------------------------------
# Analysis results
# ---------------------------------------------------------------------------

def process_analysis(input_dir, output_dir, selected_ids):
    """Subset analysis/ directory (clustering, pca, umap, diffexp)."""
    analysis_in = input_dir / "analysis"
    analysis_out = output_dir / "analysis"

    # Clustering
    for cluster_dir in sorted((analysis_in / "clustering").iterdir()):
        out_dir = analysis_out / "clustering" / cluster_dir.name
        out_dir.mkdir(parents=True, exist_ok=True)
        df = pd.read_csv(cluster_dir / "clusters.csv")
        df[df["Barcode"].isin(selected_ids)].to_csv(
            out_dir / "clusters.csv", index=False
        )

    # PCA
    for pca_dir in sorted((analysis_in / "pca").iterdir()):
        out_dir = analysis_out / "pca" / pca_dir.name
        out_dir.mkdir(parents=True, exist_ok=True)
        proj = pd.read_csv(pca_dir / "projection.csv")
        proj[proj["Barcode"].isin(selected_ids)].to_csv(
            out_dir / "projection.csv", index=False
        )
        for fname in [
            "components.csv",
            "dispersion.csv",
            "features_selected.csv",
            "variance.csv",
        ]:
            src = pca_dir / fname
            if src.exists():
                shutil.copy(src, out_dir / fname)

    # UMAP
    for umap_dir in sorted((analysis_in / "umap").iterdir()):
        out_dir = analysis_out / "umap" / umap_dir.name
        out_dir.mkdir(parents=True, exist_ok=True)
        proj = pd.read_csv(umap_dir / "projection.csv")
        proj[proj["Barcode"].isin(selected_ids)].to_csv(
            out_dir / "projection.csv", index=False
        )

    # Diffexp — copy unchanged
    for de_dir in sorted((analysis_in / "diffexp").iterdir()):
        out_dir = analysis_out / "diffexp" / de_dir.name
        out_dir.mkdir(parents=True, exist_ok=True)
        for f in de_dir.iterdir():
            shutil.copy(f, out_dir / f.name)


def process_analysis_zarr(input_dir, output_dir, selected_ids):
    """Subset analysis.zarr.zip (cluster groupings with remapped indices)."""
    # analysis.zarr indices refer to positions in the full cells.parquet table,
    # not to row numbers in clustering CSVs.
    cells = pd.read_parquet(input_dir / "cells.parquet", columns=["cell_id"])
    barcode_to_cell_idx = {cid: i for i, cid in enumerate(cells["cell_id"])}
    del cells
    gc.collect()

    # Preserve the filtered clustering CSV row order so the remapped zarr
    # indices match the subset analysis tables written by process_analysis().
    clusters = pd.read_csv(
        input_dir
        / "analysis"
        / "clustering"
        / "gene_expression_graphclust"
        / "clusters.csv"
    )
    selected_cluster_barcodes = [
        barcode for barcode in clusters["Barcode"].tolist() if barcode in selected_ids
    ]
    del clusters
    gc.collect()

    # Build full-cell-index → new analysis-row-position mapping.
    old_to_new = {
        barcode_to_cell_idx[barcode]: new_idx
        for new_idx, barcode in enumerate(selected_cluster_barcodes)
        if barcode in barcode_to_cell_idx
    }

    store_in = open_zip_read_store(input_dir / "analysis.zarr.zip", mode="r")
    z_in = zarr.open(store_in, mode="r")

    with tempfile.TemporaryDirectory() as tmpdir:
        tmpdir_path = Path(tmpdir)
        z_out = open_directory_group(tmpdir, mode="w")

        # Preserve attrs on cell_groups
        cg_out = z_out.create_group("cell_groups")
        cg_out.attrs.update(z_in["cell_groups"].attrs.asdict())

        for group_key in sorted(z_in["cell_groups"].keys(), key=int):
            grp = z_in["cell_groups"][group_key]
            indices = grp["indices"][:]
            indptr = grp["indptr"][:]

            new_indices_list = []
            new_indptr = [0]
            n_clusters = len(indptr) - 1

            for c in range(n_clusters):
                s, e = indptr[c], indptr[c + 1]
                cluster_indices = indices[s:e]
                filtered = sorted(
                    old_to_new[idx] for idx in cluster_indices if idx in old_to_new
                )
                new_indices_list.extend(filtered)
                new_indptr.append(len(new_indices_list))

            grp_out = cg_out.create_group(group_key)
            create_zarr_dataset(
                grp_out,
                "indices", data=np.array(new_indices_list, dtype=np.uint32)
            )
            create_zarr_dataset(
                grp_out,
                "indptr", data=np.array(new_indptr, dtype=np.uint32)
            )

        close_zarr_group(z_out)
        archive_directory_as_zip(tmpdir_path, output_dir / "analysis.zarr.zip")

    store_in.close()


def process_cfm_zarr(input_dir, output_dir, selected_ids):
    """Subset cell_feature_matrix.zarr.zip (CSC + CSR sparse, cell_id)."""
    # Map selected_ids to zarr row positions (same order as cells.parquet)
    cells = pd.read_parquet(input_dir / "cells.parquet", columns=["cell_id"])
    barcode_to_idx = {cid: i for i, cid in enumerate(cells["cell_id"])}
    selected_col_indices = sorted(
        barcode_to_idx[cid] for cid in selected_ids if cid in barcode_to_idx
    )
    del cells
    gc.collect()

    store_in = open_zip_read_store(input_dir / "cell_feature_matrix.zarr.zip", mode="r")
    z_in = zarr.open(store_in, mode="r")

    with tempfile.TemporaryDirectory() as tmpdir:
        tmpdir_path = Path(tmpdir)
        z_out = open_directory_group(tmpdir, mode="w")

        cf_out = z_out.create_group("cell_features")
        cf_out.attrs.update(z_in["cell_features"].attrs.asdict())

        # cell_id — subset rows
        cell_id = z_in["cell_features/cell_id"][:]
        create_zarr_dataset(cf_out, "cell_id", data=cell_id[selected_col_indices])
        del cell_id

        # --- CSC (column-oriented: one column per cell) ---
        csc_data = z_in["cell_features/csc/data"][:]
        csc_indices = z_in["cell_features/csc/indices"][:]
        csc_indptr = z_in["cell_features/csc/indptr"][:]

        new_csc_data, new_csc_indices, new_csc_indptr = subset_csc_columns(
            csc_data, csc_indices, csc_indptr, selected_col_indices
        )

        csc_out = cf_out.create_group("csc")
        create_zarr_dataset(csc_out, "data", data=new_csc_data)
        create_zarr_dataset(csc_out, "indices", data=new_csc_indices)
        create_zarr_dataset(csc_out, "indptr", data=new_csc_indptr)

        del csc_data, csc_indices, csc_indptr, new_csc_data, new_csc_indices
        gc.collect()

        # --- CSR (row-oriented: one row per feature) ---
        # Build from the subset CSC using scipy for correctness
        sub_csc_data = csc_out["data"][:]
        sub_csc_indices = csc_out["indices"][:]
        sub_csc_indptr = csc_out["indptr"][:]
        n_features = int(z_in["cell_features"].attrs["feature_ids"].__len__())
        n_cells_sub = len(selected_col_indices)

        sub_csc_mat = csc_matrix(
            (sub_csc_data, sub_csc_indices.astype(np.int32), sub_csc_indptr),
            shape=(n_features, n_cells_sub),
        )
        sub_csr_mat = sub_csc_mat.tocsr()

        create_zarr_dataset(
            cf_out,
            "data", data=sub_csr_mat.data.astype(np.uint32)
        )
        create_zarr_dataset(
            cf_out,
            "indices", data=sub_csr_mat.indices.astype(np.uint32)
        )
        create_zarr_dataset(
            cf_out,
            "indptr", data=sub_csr_mat.indptr.astype(np.uint32)
        )

        del sub_csc_mat, sub_csr_mat
        gc.collect()

        close_zarr_group(z_out)
        archive_directory_as_zip(tmpdir_path, output_dir / "cell_feature_matrix.zarr.zip")

    store_in.close()


def validate_output(output_dir):
    """Check internal consistency of all output files.

    Returns True if all checks pass, False otherwise.
    """
    print("  Running consistency checks...")
    passed = True

    def check(description, condition):
        nonlocal passed
        status = "PASS" if condition else "FAIL"
        if not condition:
            passed = False
        print(f"    [{status}] {description}")

    # --- Core cell count ---
    cells_table = pq.read_table(output_dir / "cells.parquet")
    n_cells = len(cells_table)
    cell_ids_set = set(cells_table.column("cell_id").to_pylist())
    del cells_table

    # barcodes.tsv.gz
    with gzip.open(
        output_dir / "cell_feature_matrix" / "barcodes.tsv.gz", "rt"
    ) as f:
        n_barcodes = sum(1 for line in f if line.strip())
    check(f"barcodes.tsv.gz lines ({n_barcodes}) == cells ({n_cells})",
          n_barcodes == n_cells)

    # cell_feature_matrix.h5
    with h5py.File(output_dir / "cell_feature_matrix.h5", "r") as h5:
        h5_n_barcodes = len(h5["matrix/barcodes"])
        h5_shape = h5["matrix/shape"][:]
        h5_n_features = h5_shape[0]
        h5_n_cells = h5_shape[1]
    check(f"h5 barcodes ({h5_n_barcodes}) == cells ({n_cells})",
          h5_n_barcodes == n_cells)
    check(f"h5 matrix columns ({h5_n_cells}) == cells ({n_cells})",
          h5_n_cells == n_cells)

    # cell_feature_matrix.zarr.zip
    store = open_zip_read_store(output_dir / "cell_feature_matrix.zarr.zip", mode="r")
    z = zarr.open(store, mode="r")
    zarr_cfm_n_cells = z["cell_features/cell_id"].shape[0]
    zarr_cfm_csc_indptr = z["cell_features/csc/indptr"].shape[0] - 1
    check(
        f"cfm zarr cell_id rows ({zarr_cfm_n_cells}) == cells ({n_cells})",
        zarr_cfm_n_cells == n_cells,
    )
    check(
        f"cfm zarr CSC indptr ({zarr_cfm_csc_indptr}) == cells ({n_cells})",
        zarr_cfm_csc_indptr == n_cells,
    )
    store.close()

    # cells.zarr.zip
    store = open_zip_read_store(output_dir / "cells.zarr.zip", mode="r")
    z = zarr.open(store, mode="r")
    zarr_cell_id_rows = z["cell_id"].shape[0]
    zarr_summary_rows = z["cell_summary"].shape[0]
    check(
        f"cells zarr cell_id rows ({zarr_cell_id_rows}) == cells ({n_cells})",
        zarr_cell_id_rows == n_cells,
    )
    check(
        f"cells zarr cell_summary rows ({zarr_summary_rows}) == cells ({n_cells})",
        zarr_summary_rows == n_cells,
    )
    store.close()

    # --- Feature count ---
    with gzip.open(
        output_dir / "cell_feature_matrix" / "features.tsv.gz", "rt"
    ) as f:
        n_features = sum(1 for line in f if line.strip())
    check(
        f"features.tsv.gz lines ({n_features}) == h5 features ({h5_n_features})",
        n_features == h5_n_features,
    )

    store = open_zip_read_store(output_dir / "cell_feature_matrix.zarr.zip", mode="r")
    z = zarr.open(store, mode="r")
    zarr_csr_indptr = z["cell_features/indptr"].shape[0] - 1
    zarr_feature_ids = list(z["cell_features"].attrs["feature_ids"])
    zarr_n_feature_ids = len(zarr_feature_ids)
    has_total_transcripts = bool(zarr_feature_ids) and zarr_feature_ids[-1] == "Total transcripts"
    check(
        f"cfm zarr CSR indptr ({zarr_csr_indptr}) == feature_ids attr ({zarr_n_feature_ids})",
        zarr_csr_indptr == zarr_n_feature_ids,
    )
    check(
        f"cfm zarr feature_ids attr ({zarr_n_feature_ids}) matches features ({n_features}) or features+1 with Total transcripts",
        zarr_n_feature_ids == n_features or (has_total_transcripts and zarr_n_feature_ids == n_features + 1),
    )
    store.close()

    # --- Analysis cell counts ---
    analysis_out = output_dir / "analysis"

    # Collect all clustering row counts
    cluster_counts = {}
    for cluster_dir in sorted((analysis_out / "clustering").iterdir()):
        df = pd.read_csv(cluster_dir / "clusters.csv")
        cluster_counts[cluster_dir.name] = len(df)
    unique_cluster_counts = set(cluster_counts.values())
    check(
        f"all clustering CSVs have same row count ({unique_cluster_counts})",
        len(unique_cluster_counts) == 1,
    )
    n_analysis_cells = list(cluster_counts.values())[0] if cluster_counts else 0

    # PCA projections
    for pca_dir in sorted((analysis_out / "pca").iterdir()):
        proj = pd.read_csv(pca_dir / "projection.csv")
        check(
            f"pca/{pca_dir.name} rows ({len(proj)}) == analysis cells ({n_analysis_cells})",
            len(proj) == n_analysis_cells,
        )

    # UMAP projections
    for umap_dir in sorted((analysis_out / "umap").iterdir()):
        proj = pd.read_csv(umap_dir / "projection.csv")
        check(
            f"umap/{umap_dir.name} rows ({len(proj)}) == analysis cells ({n_analysis_cells})",
            len(proj) == n_analysis_cells,
        )

    # analysis.zarr.zip — total indices per grouping should equal n_analysis_cells
    store = open_zip_read_store(output_dir / "analysis.zarr.zip", mode="r")
    z = zarr.open(store, mode="r")
    for group_key in sorted(z["cell_groups"].keys(), key=int):
        n_indices = z["cell_groups"][group_key]["indices"].shape[0]
        check(
            f"analysis zarr group {group_key} indices ({n_indices}) == analysis cells ({n_analysis_cells})",
            n_indices == n_analysis_cells,
        )
    store.close()

    # --- Boundary files ---
    cb = pq.read_table(
        output_dir / "cell_boundaries.parquet", columns=["cell_id"]
    )
    cb_ids = set(cb.column("cell_id").to_pylist())
    check(
        f"cell_boundaries cell_ids ({len(cb_ids)}) all in cells.parquet",
        cb_ids.issubset(cell_ids_set),
    )
    del cb

    nb = pq.read_table(
        output_dir / "nucleus_boundaries.parquet", columns=["cell_id"]
    )
    nb_ids = set(nb.column("cell_id").to_pylist())
    check(
        f"nucleus_boundaries cell_ids ({len(nb_ids)}) all in cells.parquet",
        nb_ids.issubset(cell_ids_set),
    )
    del nb

    # --- Transcripts ---
    parquet_file = pq.ParquetFile(output_dir / "transcripts.parquet")
    transcript_cell_ids = set()
    for batch in parquet_file.iter_batches(
        batch_size=1_000_000, columns=["cell_id"]
    ):
        table = pa.Table.from_batches([batch])
        ids = table.column("cell_id").to_pylist()
        transcript_cell_ids.update(ids)
    transcript_cell_ids.discard("UNASSIGNED")
    check(
        f"transcript assigned cell_ids ({len(transcript_cell_ids)}) all in cells.parquet",
        transcript_cell_ids.issubset(cell_ids_set),
    )

    gc.collect()
    return passed


# ---------------------------------------------------------------------------
# Orchestration
# ---------------------------------------------------------------------------

def copy_and_crop_images(input_dir, output_dir, region, pixel_size, he_image=None, he_alignment=None):
    """Copy metadata files and crop morphology/H&E images to the region."""
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


def run_region(input_dir, output_dir, region, proportion, grid_size, pixel_size, skip_validation, he_image=None, he_alignment=None):
    """Run all crop/subset steps for a single region and write the output directory."""
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

    # Zarr is processed before parquets so label remaps are available.
    # spatialdata-io 0.7.0 requires parquet label_id values to be sequential
    # 1..M matching polygon_sets positions. The remap dicts translate from
    # original label_ids to new sequential values.
    print("\nStep 2: Cropping cells.zarr.zip...")
    nucleus_label_id_remap, cell_label_id_remap = process_cells_zarr_region(
        input_dir, output_dir, selected_indices, n_total, region, pixel_size
    )

    print("\nStep 3: Subsetting and rebasing tabular cell files...")
    print("  cells...")
    subset_spatial_parquet_and_csv(input_dir, output_dir, "cells", selected_ids, region)
    print("  cell_boundaries...")
    subset_spatial_parquet_and_csv(
        input_dir, output_dir, "cell_boundaries", selected_ids, region,
        label_id_remap=cell_label_id_remap,
    )
    print("  nucleus_boundaries...")
    subset_spatial_parquet_and_csv(
        input_dir, output_dir, "nucleus_boundaries", selected_ids, region,
        label_id_remap=nucleus_label_id_remap,
    )

    print("\nStep 4: Subsetting and rebasing transcripts...")
    n_transcripts = process_transcripts_in_region(
        input_dir, output_dir, selected_ids, region, proportion
    )
    print(f"  Kept {n_transcripts:,} transcripts")

    print("\nStep 5: Subsetting cell feature matrix...")
    print("  cell_feature_matrix/ directory...")
    process_cell_feature_matrix_dir(input_dir, output_dir, selected_ids)
    print("  cell_feature_matrix.h5...")
    process_cell_feature_matrix_h5(input_dir, output_dir, selected_ids)
    print("  cell_feature_matrix.tar.gz...")
    process_cell_feature_matrix_tar(output_dir)

    print("\nStep 6: Subsetting analysis results...")
    process_analysis(input_dir, output_dir, selected_ids)

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
    pixel_size = args.pixel_size if args.pixel_size else XENIUM_PIXEL_SIZE_UM
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
