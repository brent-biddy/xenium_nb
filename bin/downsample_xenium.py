#!/usr/bin/env python3
"""
downsample_xenium.py - Downsample Xenium output files using spatial grid sampling.

Selects a proportion of cells from each spatial grid square, then regenerates
all Xenium output files for the selected cells.

Usage:
    python downsample_xenium.py <input_dir> [--proportion 0.05] [--grid_size 100.0]

Output directory is auto-named: <input_dir>_downsampled_<pct>pct
"""

import argparse
import gc
import gzip
import shutil
import sys
import tarfile
import tempfile
import time
from pathlib import Path

import xml.etree.ElementTree as ET

import h5py
import numpy as np
import pandas as pd
import pyarrow as pa
import pyarrow.parquet as pq
import tifffile
import zarr
from scipy.io import mmread, mmwrite
from scipy.sparse import csc_matrix


def open_zip_store(path, mode):
    zip_store_cls = getattr(zarr, "ZipStore", None)
    if zip_store_cls is None:
        zip_store_cls = zarr.storage.ZipStore
    return zip_store_cls(str(path), mode=mode)


def create_zarr_dataset(group, name, data=None, **kwargs):
    if data is not None:
        kwargs.setdefault("shape", np.shape(data))
        kwargs.setdefault("dtype", getattr(data, "dtype", np.asarray(data).dtype))
        return group.create_dataset(name, data=data, **kwargs)
    return group.create_dataset(name, **kwargs)


def parse_args():
    parser = argparse.ArgumentParser(
        description="Downsample Xenium output files using spatial grid sampling."
    )
    parser.add_argument("input_dir", type=Path, help="Path to Xenium output directory")
    parser.add_argument(
        "--proportion",
        type=float,
        default=0.05,
        help="Fraction of cells to keep (default: 0.05)",
    )
    parser.add_argument(
        "--grid_size",
        type=float,
        default=100.0,
        help="Side length of grid squares in micrometers (default: 100.0)",
    )
    parser.add_argument(
        "--image_level",
        type=int,
        default=5,
        help="Lowest pyramid level to keep in output images (0=full res, 7=smallest; default: 5)",
    )
    return parser.parse_args()


# ---------------------------------------------------------------------------
# Step 1: Cell selection
# ---------------------------------------------------------------------------

def select_cells(input_dir, proportion, grid_size):
    """Select cells using spatial grid sampling.

    Returns (selected_indices, selected_ids, n_total) where
    selected_indices is a sorted int array of 0-based row positions and
    selected_ids is a set of string cell_id values.
    """
    rng = np.random.default_rng(42)

    cells_table = pq.read_table(input_dir / "cells.parquet")
    n_total = len(cells_table)
    cell_ids = cells_table.column("cell_id").to_pylist()
    x = cells_table.column("x_centroid").to_numpy()
    y = cells_table.column("y_centroid").to_numpy()
    del cells_table
    gc.collect()

    grid_x = ((x - x.min()) / grid_size).astype(int)
    grid_y = ((y - y.min()) / grid_size).astype(int)

    # Group indices by grid square
    grid_dict = {}
    for i in range(len(x)):
        key = (int(grid_x[i]), int(grid_y[i]))
        grid_dict.setdefault(key, []).append(i)

    selected_mask = np.zeros(len(x), dtype=bool)
    for indices in grid_dict.values():
        n = max(1, round(len(indices) * proportion))
        chosen = rng.choice(indices, size=n, replace=False)
        selected_mask[chosen] = True

    selected_indices = np.where(selected_mask)[0]
    selected_ids = {cell_ids[i] for i in selected_indices}

    return selected_indices, selected_ids, n_total


# ---------------------------------------------------------------------------
# Step 2: Tabular cell files
# ---------------------------------------------------------------------------

def subset_parquet_and_csv(input_dir, output_dir, filename, selected_ids,
                           id_column="cell_id"):
    """Filter a parquet file by id_column, write parquet + csv.gz."""
    table = pq.read_table(input_dir / f"{filename}.parquet")
    mask = pa.compute.is_in(
        table.column(id_column), value_set=pa.array(list(selected_ids))
    )
    sub = table.filter(mask)
    pq.write_table(sub, output_dir / f"{filename}.parquet")
    sub.to_pandas().to_csv(
        output_dir / f"{filename}.csv.gz", index=False, compression="gzip"
    )
    n = len(sub)
    del table, sub, mask
    gc.collect()
    return n


# ---------------------------------------------------------------------------
# Step 3: Transcripts
# ---------------------------------------------------------------------------

def process_transcripts(input_dir, output_dir, selected_ids, proportion):
    """Subset transcripts with chunked reading.

    Assigned transcripts are kept if cell_id is in selected_ids.
    UNASSIGNED transcripts are randomly sampled at the target proportion.

    feature_name is written as dictionary-encoded (categorical) so that
    spatialdata-io/PointsModel can resolve categories without a full
    in-memory compute of the transcript table.
    """
    rng = np.random.default_rng(42)
    selected_ids_arrow = pa.array(list(selected_ids))

    # Pre-pass: collect all unique feature_name values for a consistent
    # dictionary across all row groups (required for dask to see known categories).
    print("    Collecting feature_name categories...")
    all_feature_names = set()
    for batch in pq.ParquetFile(input_dir / "transcripts.parquet").iter_batches(
        batch_size=2_000_000, columns=["feature_name"]
    ):
        all_feature_names.update(
            pa.Table.from_batches([batch]).column("feature_name").to_pylist()
        )
    all_feature_names.discard(None)
    feature_name_dict = pa.array(sorted(all_feature_names), type=pa.large_utf8())
    feature_to_idx = {v: i for i, v in enumerate(feature_name_dict.to_pylist())}
    dict_type = pa.dictionary(pa.int16(), pa.large_utf8())

    parquet_file = pq.ParquetFile(input_dir / "transcripts.parquet")
    writer = None
    output_schema = None
    csv_path = output_dir / "transcripts.csv"
    first_chunk = True
    total_kept = 0

    for batch in parquet_file.iter_batches(batch_size=1_000_000):
        table = pa.Table.from_batches([batch])
        cell_id_col = table.column("cell_id")

        # Assigned transcripts belonging to selected cells
        assigned_mask = pa.compute.is_in(cell_id_col, value_set=selected_ids_arrow)

        # UNASSIGNED transcripts sampled at target proportion
        unassigned_mask = pa.compute.equal(cell_id_col, "UNASSIGNED")
        n_unassigned = pa.compute.sum(unassigned_mask).as_py()

        if n_unassigned > 0:
            random_vals = rng.random(len(table))
            sample_mask = pa.array(random_vals < proportion)
            unassigned_keep = pa.compute.and_(unassigned_mask, sample_mask)
            combined_mask = pa.compute.or_(assigned_mask, unassigned_keep)
        else:
            combined_mask = assigned_mask

        chunk_sub = table.filter(combined_mask)

        if len(chunk_sub) > 0:
            # Encode feature_name as dictionary with the global dictionary so
            # that categories are known when dask reads the parquet back.
            fn_idx = chunk_sub.schema.get_field_index("feature_name")
            indices = pa.array(
                [feature_to_idx[v] for v in chunk_sub.column("feature_name").to_pylist()],
                type=pa.int16(),
            )
            dict_col = pa.DictionaryArray.from_arrays(indices, feature_name_dict)
            chunk_sub = chunk_sub.set_column(fn_idx, pa.field("feature_name", dict_type), dict_col)

            if writer is None:
                output_schema = chunk_sub.schema
                writer = pq.ParquetWriter(output_dir / "transcripts.parquet", output_schema)
            writer.write_table(chunk_sub)

            # Cast feature_name back to string for CSV output
            chunk_csv = chunk_sub.set_column(
                fn_idx, pa.field("feature_name", pa.large_utf8()),
                dict_col.cast(pa.large_utf8()),
            )
            chunk_csv.to_pandas().to_csv(
                csv_path, mode="a", index=False, header=first_chunk
            )
            first_chunk = False
            total_kept += len(chunk_sub)

        del table, chunk_sub
        gc.collect()

    if writer:
        writer.close()

    # Gzip the CSV
    with open(csv_path, "rb") as f_in:
        with gzip.open(output_dir / "transcripts.csv.gz", "wb") as f_out:
            shutil.copyfileobj(f_in, f_out)
    csv_path.unlink()

    return total_kept


# ---------------------------------------------------------------------------
# Step 4: Cell feature matrix (directory, h5, tar.gz)
# ---------------------------------------------------------------------------

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

        new_data, new_indices, new_indptr = [], [], [0]
        for col_idx in selected_col_indices:
            s, e = indptr[col_idx], indptr[col_idx + 1]
            new_data.append(data[s:e])
            new_indices.append(indices[s:e])
            new_indptr.append(new_indptr[-1] + (e - s))

        new_data = (
            np.concatenate(new_data)
            if new_data
            else np.array([], dtype=data.dtype)
        )
        new_indices = (
            np.concatenate(new_indices)
            if new_indices
            else np.array([], dtype=indices.dtype)
        )
        new_indptr = np.array(new_indptr, dtype=indptr.dtype)
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
# Step 5: Analysis results
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


# ---------------------------------------------------------------------------
# Step 6: cells.zarr.zip
# ---------------------------------------------------------------------------

def process_cells_zarr(input_dir, output_dir, selected_indices, n_total):
    """Subset cells.zarr.zip (cell_id, cell_summary, masks, polygon_sets)."""
    store_in = open_zip_store(input_dir / "cells.zarr.zip", mode="r")
    z_in = zarr.open(store_in, mode="r")

    store_out = open_zip_store(output_dir / "cells.zarr.zip", mode="w")
    z_out = zarr.open(store_out, mode="w")

    # cell_id and cell_summary — subset rows
    create_zarr_dataset(z_out, "cell_id", data=z_in["cell_id"][:][selected_indices])
    create_zarr_dataset(
        z_out,
        "cell_summary",
        data=z_in["cell_summary"][:][selected_indices],
    )
    z_out["cell_summary"].attrs.update(z_in["cell_summary"].attrs.asdict())

    # Build a fast lookup array for mask filtering
    # Mask pixel values are label_ids (1-based = row_index + 1)
    selected_label_ids = selected_indices + 1
    lookup = np.zeros(n_total + 1, dtype=bool)
    lookup[selected_label_ids] = True

    # Masks — process in row chunks
    masks_out = z_out.create_group("masks")
    create_zarr_dataset(
        masks_out,
        "homogeneous_transform", data=z_in["masks/homogeneous_transform"][:]
    )

    chunk_size = 500
    for mask_name in ["0", "1"]:
        print(f"    mask {mask_name}...")
        mask_in = z_in[f"masks/{mask_name}"]
        rows, cols = mask_in.shape
        mask_arr_out = masks_out.create_dataset(
            mask_name,
            shape=(rows, cols),
            dtype=mask_in.dtype,
            chunks=mask_in.chunks,
        )
        for start in range(0, rows, chunk_size):
            end = min(start + chunk_size, rows)
            chunk = mask_in[start:end, :]
            # Zero out pixels not belonging to selected cells
            nonzero = chunk > 0
            if nonzero.any():
                in_range = chunk <= n_total
                keep = np.zeros_like(chunk, dtype=bool)
                valid = nonzero & in_range
                keep[valid] = lookup[chunk[valid]]
                chunk[nonzero & ~keep] = 0
            mask_arr_out[start:end, :] = chunk

    # Polygon sets — filter by cell_index and remap to new consecutive indices
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
            create_zarr_dataset(ps_out, key, data=data)
        if ps_in.attrs:
            ps_out.attrs.update(ps_in.attrs.asdict())

    store_out.close()
    store_in.close()


# ---------------------------------------------------------------------------
# Step 7: analysis.zarr.zip
# ---------------------------------------------------------------------------

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

    store_in = open_zip_store(input_dir / "analysis.zarr.zip", mode="r")
    z_in = zarr.open(store_in, mode="r")

    store_out = open_zip_store(output_dir / "analysis.zarr.zip", mode="w")
    z_out = zarr.open(store_out, mode="w")

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

    store_out.close()
    store_in.close()


# ---------------------------------------------------------------------------
# Step 8: cell_feature_matrix.zarr.zip
# ---------------------------------------------------------------------------

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

    store_in = open_zip_store(input_dir / "cell_feature_matrix.zarr.zip", mode="r")
    z_in = zarr.open(store_in, mode="r")

    store_out = open_zip_store(output_dir / "cell_feature_matrix.zarr.zip", mode="w")
    z_out = zarr.open(store_out, mode="w")

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

    new_csc_data, new_csc_indices, new_csc_indptr = [], [], [0]
    for col_idx in selected_col_indices:
        s, e = csc_indptr[col_idx], csc_indptr[col_idx + 1]
        new_csc_data.append(csc_data[s:e])
        new_csc_indices.append(csc_indices[s:e])
        new_csc_indptr.append(new_csc_indptr[-1] + (e - s))

    csc_out = cf_out.create_group("csc")
    create_zarr_dataset(
        csc_out,
        "data",
        data=(
            np.concatenate(new_csc_data)
            if new_csc_data
            else np.array([], dtype=csc_data.dtype)
        ),
    )
    create_zarr_dataset(
        csc_out,
        "indices",
        data=(
            np.concatenate(new_csc_indices)
            if new_csc_indices
            else np.array([], dtype=csc_indices.dtype)
        ),
    )
    create_zarr_dataset(
        csc_out,
        "indptr", data=np.array(new_csc_indptr, dtype=csc_indptr.dtype)
    )

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

    store_out.close()
    store_in.close()


# ---------------------------------------------------------------------------
# Step 9: Copy unchanged files
# ---------------------------------------------------------------------------

def extract_pyramid_sublevels(input_dir, output_dir, image_level):
    """Write reduced-resolution pyramid OME-TIFFs for morphology images.

    Keeps pyramid levels >= image_level (e.g. image_level=5 keeps levels 5,6,7).
    Outputs OME-TIFF files with embedded OME-XML metadata.
    """
    write_opts = dict(tile=(256, 256), compression="deflate")

    _ome_ns = {"ome": "http://www.openmicroscopy.org/Schemas/OME/2016-06"}
    _scale_factor = 2 ** image_level  # each pyramid level is 2x downsampled

    def _parse_ome_metadata(tif_path):
        """Return (channel_names, physical_size_x, physical_size_y, size_unit) from OME-XML."""
        with tifffile.TiffFile(str(tif_path)) as tif:
            if tif.ome_metadata:
                root = ET.fromstring(tif.ome_metadata)
                names = [
                    ch.get("Name", "")
                    for ch in root.findall(".//ome:Channel", _ome_ns)
                ]
                pixels = root.find(".//ome:Pixels", _ome_ns)
                px = float(pixels.get("PhysicalSizeX", 1.0)) if pixels is not None else 1.0
                py = float(pixels.get("PhysicalSizeY", 1.0)) if pixels is not None else 1.0
                unit = pixels.get("PhysicalSizeXUnit", "µm") if pixels is not None else "µm"
                return names, px, py, unit
        return [], 1.0, 1.0, "µm"

    def _normalize_axes(series_axes, array_ndim):
        """Preserve source axes where possible while dropping non-data dims."""
        axes = series_axes.replace("S", "")
        if len(axes) == array_ndim:
            return axes
        if len(axes) > array_ndim:
            axes = axes[-array_ndim:]
        if len(axes) != array_ndim:
            raise ValueError(
                f"Cannot reconcile source axes '{series_axes}' with array ndim {array_ndim}"
            )
        return axes

    def _build_ome_metadata(axes, px, py, unit, channel_names=None):
        meta = {
            "axes": axes,
            "PhysicalSizeX": px,
            "PhysicalSizeXUnit": unit,
            "PhysicalSizeY": py,
            "PhysicalSizeYUnit": unit,
        }
        if "C" in axes and channel_names:
            c_len = len(channel_names)
            meta["Channel"] = {"Name": channel_names[:c_len]}
        return meta

    # --- morphology.ome.tif (ZYX pyramid via zarr) ---
    src = input_dir / "morphology.ome.tif"
    if src.exists():
        names, px, py, unit = _parse_ome_metadata(src)
        # Scale physical pixel size to match the new base level
        scaled_px = px * _scale_factor
        scaled_py = py * _scale_factor
        with tifffile.TiffFile(str(src)) as tif:
            source_axes = tif.series[0].axes
            store = tif.aszarr()
            z = zarr.open(store, mode="r")
            level_keys = sorted(z.keys(), key=int)
            use_keys = [k for k in level_keys if int(k) >= image_level]
            arrays = [z[k][:] for k in use_keys]
        shapes = " -> ".join(f"{a.shape[-2]}x{a.shape[-1]}" for a in arrays)
        print(f"    morphology.ome.tif: levels {use_keys[0]}-{use_keys[-1]} ({shapes})")
        print(f"    physical pixel size: {px} -> {scaled_px} {unit}/px")
        axes = _normalize_axes(source_axes, arrays[0].ndim)
        channel_names = names if "C" in axes else None
        meta = _build_ome_metadata(axes, scaled_px, scaled_py, unit, channel_names)
        with tifffile.TiffWriter(str(output_dir / "morphology.ome.tif")) as writer:
            writer.write(arrays[0], subifds=len(arrays) - 1, metadata=meta, **write_opts)
            for arr in arrays[1:]:
                writer.write(arr, subfiletype=1, **write_opts)
        del arrays
        gc.collect()

    # --- morphology_focus/ch*.ome.tif (single-channel SubIFD pyramids) ---
    focus_in = input_dir / "morphology_focus"
    focus_out = output_dir / "morphology_focus"
    focus_out.mkdir(exist_ok=True)
    if focus_in.exists():
        # Read channel names and physical size from the first focus file's OME metadata.
        # Each morphology_focus file is itself a CYX pyramid, so preserve that full
        # channel axis rather than flattening sub-pages to 2D planes.
        focus_files = sorted(focus_in.glob("*.ome.tif"))
        all_channel_names, px, py, unit = _parse_ome_metadata(focus_files[0]) if focus_files else ([], 1.0, 1.0, "µm")
        scaled_px = px * _scale_factor
        scaled_py = py * _scale_factor

        for src in focus_files:
            out_name = src.name

            with tifffile.TiffFile(str(src)) as tif:
                source_axes = tif.series[0].axes
                levels = tif.series[0].levels
                use_levels = levels[image_level:]
                arrays = [lvl.asarray() for lvl in use_levels]
            shapes = " -> ".join(f"{a.shape[-2]}x{a.shape[-1]}" for a in arrays)
            print(f"    {out_name}: levels {image_level}-{image_level + len(arrays) - 1} ({shapes})")
            axes = _normalize_axes(source_axes, arrays[0].ndim)
            channel_names = all_channel_names if "C" in axes else None
            with tifffile.TiffWriter(str(focus_out / out_name)) as writer:
                writer.write(
                    arrays[0],
                    subifds=len(arrays) - 1,
                    metadata=_build_ome_metadata(axes, scaled_px, scaled_py, unit, channel_names),
                    **write_opts,
                )
                for arr in arrays[1:]:
                    writer.write(arr, subfiletype=1, **write_opts)
            del arrays
            gc.collect()


def copy_unchanged(input_dir, output_dir, image_level):
    """Copy files that are not cell-indexed and need no modification.

    Morphology images are written as reduced-resolution pyramid TIFFs
    (levels >= image_level) instead of copying the full gigabyte-scale originals.
    """
    copy_files = [
        "experiment.xenium",
        "gene_panel.json",
        "metrics_summary.csv",
        "analysis_summary.html",
    ]
    for fname in copy_files:
        src = input_dir / fname
        if src.exists():
            print(f"    {fname} ({src.stat().st_size / 1e9:.1f} GB)...")
            shutil.copy2(str(src), str(output_dir / fname))

    print(f"    Extracting pyramid levels >= {image_level} from morphology images...")
    extract_pyramid_sublevels(input_dir, output_dir, image_level)


# ---------------------------------------------------------------------------
# Step 10: Validation
# ---------------------------------------------------------------------------

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
    store = open_zip_store(output_dir / "cell_feature_matrix.zarr.zip", mode="r")
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
    store = open_zip_store(output_dir / "cells.zarr.zip", mode="r")
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

    store = open_zip_store(output_dir / "cell_feature_matrix.zarr.zip", mode="r")
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
    store = open_zip_store(output_dir / "analysis.zarr.zip", mode="r")
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
# Main
# ---------------------------------------------------------------------------

def main():
    args = parse_args()
    input_dir = args.input_dir.resolve()
    proportion = args.proportion
    grid_size = args.grid_size
    image_level = args.image_level

    pct_str = f"{proportion * 100:g}"
    output_dir = input_dir.parent / f"{input_dir.name}_downsampled_{pct_str}pct"
    output_dir.mkdir(parents=True, exist_ok=True)

    print(f"Input:      {input_dir}")
    print(f"Output:     {output_dir}")
    print(f"Proportion: {proportion}")
    print(f"Grid size:  {grid_size} um")
    print()

    t0 = time.time()

    # Step 1 — Cell selection
    print("Step 1: Selecting cells...")
    selected_indices, selected_ids, n_total = select_cells(
        input_dir, proportion, grid_size
    )
    print(
        f"  Selected {len(selected_ids):,} / {n_total:,} cells "
        f"({len(selected_ids) / n_total * 100:.1f}%)"
    )

    # Step 2 — Tabular cell files
    print("\nStep 2: Subsetting tabular cell files...")
    for name in ["cells", "cell_boundaries", "nucleus_boundaries"]:
        print(f"  {name}...")
        subset_parquet_and_csv(input_dir, output_dir, name, selected_ids)

    # Step 3 — Transcripts
    print("\nStep 3: Subsetting transcripts (chunked)...")
    n_transcripts = process_transcripts(
        input_dir, output_dir, selected_ids, proportion
    )
    print(f"  Kept {n_transcripts:,} transcripts")

    # Step 4 — Cell feature matrix
    print("\nStep 4: Subsetting cell feature matrix...")
    print("  cell_feature_matrix/ directory...")
    process_cell_feature_matrix_dir(input_dir, output_dir, selected_ids)
    print("  cell_feature_matrix.h5...")
    process_cell_feature_matrix_h5(input_dir, output_dir, selected_ids)
    print("  cell_feature_matrix.tar.gz...")
    process_cell_feature_matrix_tar(output_dir)

    # Step 5 — Analysis
    print("\nStep 5: Subsetting analysis results...")
    process_analysis(input_dir, output_dir, selected_ids)

    # Step 6 — cells.zarr.zip
    print("\nStep 6: Subsetting cells.zarr.zip...")
    process_cells_zarr(input_dir, output_dir, selected_indices, n_total)

    # Step 7 — analysis.zarr.zip
    print("\nStep 7: Subsetting analysis.zarr.zip...")
    process_analysis_zarr(input_dir, output_dir, selected_ids)

    # Step 8 — cell_feature_matrix.zarr.zip
    print("\nStep 8: Subsetting cell_feature_matrix.zarr.zip...")
    process_cfm_zarr(input_dir, output_dir, selected_ids)

    # Step 9 — Copy unchanged files
    print("\nStep 9: Copying unchanged files...")
    copy_unchanged(input_dir, output_dir, image_level)

    # Step 10 — Validation
    print("\nStep 10: Validating output...")
    all_passed = validate_output(output_dir)

    elapsed = time.time() - t0
    print(f"\nDone in {elapsed / 60:.1f} minutes.")
    print(f"Output: {output_dir}")
    if not all_passed:
        print("WARNING: Some validation checks failed. Review output above.")
        sys.exit(1)


if __name__ == "__main__":
    main()
