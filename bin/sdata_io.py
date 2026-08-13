"""sdata_io.py - shared SpatialData write helpers for the bin/ scripts.

Exists for one reason: the cluster_sdata* steps change only the table, and
sdata.write() re-serialises the whole store to record that.
"""

import os
import shutil

import spatialdata


def write_table_only(sdata, input_path, output_path, table_key):
    """Write a clustered store that shares its images with the input store.

    The cluster_sdata* steps change exactly one element -- the table -- but
    `sdata.write()` re-serialises every element, so a whole-slide sample spends
    minutes copying morphology, H&E, focus and label images byte for byte to
    record ten new obs columns. Measured on roi_6 (638k cells), that write was
    418.7s of a 487.5s GPU run: 86% of the job, against 68.8s of actual compute.
    On the CPU step it was 709s.

    So: hardlink the input store to the output path, which costs a few hundred
    link() calls and no data movement (these are zarr v3 stores with large
    chunks -- 357 files for a 1.1 GB ROI), then replace only the table in place.

    The safety invariant is that the shared elements are never mutated. Deleting
    a hardlink does not touch the source, and writing the new table creates fresh
    inodes, so the input's table is unaffected -- verified on a real ROI. But
    anything that later modified an image chunk *in place* would corrupt the
    input store too. Nothing in this pipeline does; every downstream step reads
    images and writes new stores. Keep it that way.

    Note `write_element(overwrite=True)` cannot do this on its own: spatialdata
    refuses to overwrite a path it currently has open ("The target path of the
    write operation is in use"). Deleting the element from disk first is the
    supported route.
    """
    # Nextflow stages the input as a symlink into the work dir; hardlinks must be
    # made against the real file, and only work within one filesystem. The local
    # and oscer profiles both keep work/ and results/ under one root precisely so
    # that holds (see the hardlink publishing note in nextflow.config), but a
    # hand-run pointing across filesystems must still produce a correct store.
    source = os.path.realpath(input_path)
    try:
        shutil.copytree(source, output_path, copy_function=_link_data_copy_metadata)
    except (OSError, shutil.Error) as err:
        if not _is_cross_device(err):
            raise
        # EXDEV: fall back to a real copy. Correct, just as slow as before.
        print(f"Input is on a different filesystem than {output_path}; "
              f"copying instead of hardlinking.")
        shutil.rmtree(output_path, ignore_errors=True)
        shutil.copytree(source, output_path)

    written = spatialdata.read_zarr(output_path)
    written.delete_element_from_disk(table_key)
    written.tables[table_key] = sdata.tables[table_key]
    written.write_element(table_key)


def _link_data_copy_metadata(src, dst):
    """Hardlink chunk data, but copy zarr metadata files by content.

    Chunk files are the bulk and are only ever read, so sharing them is free and
    safe. Metadata is different: writing the new table rewrites the enclosing
    group's `zarr.json`, and on a shared inode that edit would reach back into
    the input store. Verified on a real ROI -- with a plain link of everything,
    `tables/zarr.json` was the one file left sharing an inode after the write.
    The content happened to be identical there, so nothing broke, but that is
    luck rather than a guarantee. These files are small and few (5 of 357 on a
    1.1 GB ROI), so copying them costs nothing and removes the whole class.
    """
    name = os.path.basename(src)
    if name == "zarr.json" or name.startswith(".z"):
        shutil.copy2(src, dst)
    else:
        os.link(src, dst)


def _is_cross_device(err):
    """True if err is (or wraps) an EXDEV cross-device link failure.

    shutil.copytree collects per-file failures into a shutil.Error whose args
    are (src, dst, message) triples with the message already stringified, so
    there is no errno to inspect on that path.
    """
    if isinstance(err, OSError):
        return err.errno == 18
    return "Invalid cross-device link" in str(err)
