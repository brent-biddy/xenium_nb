#!/bin/bash
#SBATCH --job-name=xeniumranger_snuc
#SBATCH --partition=normal
#SBATCH --cpus-per-task=32
#SBATCH --mem=128G
#SBATCH --time=04:00:00
#SBATCH --output=xeniumranger_snuc_%x_%j.out
#SBATCH --error=xeniumranger_snuc_%x_%j.err

# Standalone OSCER job. Not part of the Nextflow pipeline.
#
# Generates a single-nucleus (rather than single-cell) Xenium count matrix by
# re-importing segmentation with --expansion-distance=0, i.e. cells are not
# expanded past the nucleus boundary, so transcripts are only assigned to the
# nucleus itself.
#
# Usage:
#   sbatch xeniumranger_import_segmentation.sbatch <id> <xenium_bundle> <nuclei_zarr_zip> [expansion_distance]
#
# Example:
#   sbatch xeniumranger_import_segmentation.sbatch ROI1_A \
#       /path/to/xenium/files /path/to/cells.zarr.zip
#
# Output is written to ./<id>/ under the submission directory (xeniumranger's
# default behavior) — submit this from a scratch directory with enough space.
#
# --localcores/--localmem must be kept in sync with the #SBATCH --cpus-per-task
# and --mem directives above if you change the resource request.
#
# xeniumranger is not on PATH on OSCER. Set XENIUMRANGER_BIN below to the
# absolute path of the binary, or override at submit time without editing
# this file:
#   XENIUMRANGER_BIN=/path/to/xeniumranger sbatch xeniumranger_import_segmentation.sbatch ...

set -euo pipefail

: "${XENIUMRANGER_BIN:=/path/to/xeniumranger}"

if [ "$#" -lt 3 ]; then
    echo "Usage: sbatch $0 <id> <xenium_bundle> <nuclei_zarr_zip> [expansion_distance]" >&2
    exit 1
fi

if [ ! -x "${XENIUMRANGER_BIN}" ]; then
    echo "XENIUMRANGER_BIN is not an executable path: ${XENIUMRANGER_BIN}" >&2
    echo "Set it in this script's XENIUMRANGER_BIN default, or export XENIUMRANGER_BIN before sbatch." >&2
    exit 1
fi

id="$1"
xenium_bundle="$2"
nuclei="$3"
expansion_distance="${4:-0}"
localmem_gb=128

"${XENIUMRANGER_BIN}" import-segmentation \
    --id="${id}" \
    --xenium-bundle="${xenium_bundle}" \
    --nuclei="${nuclei}" \
    --expansion-distance="${expansion_distance}" \
    --localcores="${SLURM_CPUS_PER_TASK}" \
    --localmem="${localmem_gb}"
