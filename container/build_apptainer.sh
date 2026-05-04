#!/usr/bin/env bash
set -euo pipefail

script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
output_path="${1:-$script_dir/xenium_tools_squidpy_local.sif}"

mkdir -p "$(dirname "$output_path")"

cd "$script_dir"
apptainer build --force "$output_path" Apptainer.def
