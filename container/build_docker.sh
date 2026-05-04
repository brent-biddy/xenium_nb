#!/usr/bin/env bash
set -euo pipefail

script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
image_name="${1:-xenium_tools_squidpy:local}"

cd "$script_dir"
docker build -t "$image_name" .
