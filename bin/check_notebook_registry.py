#!/usr/bin/env python3
# Validates that every param listed in NotebookRegistry.groovy exists in the
# corresponding notebook's YAML front matter. Run from the repo root.

import re
import sys
from pathlib import Path

import yaml

REGISTRY_PATH = Path("lib/NotebookRegistry.groovy")


def parse_registry(text: str) -> dict[str, dict]:
    """Extract {notebook_id: {path, params}} from NotebookRegistry.groovy."""
    entries = {}
    # Match each top-level notebook block: id: [ path: "...", params: [...] ]
    for block in re.finditer(r'(\w+):\s*\[\s*path\s*:\s*"([^"]+)".*?params\s*:\s*\[([^\]]*)\]', text, re.DOTALL):
        notebook_id = block.group(1)
        path = block.group(2)
        params_raw = block.group(3)
        params = re.findall(r"'([^']+)'", params_raw)
        entries[notebook_id] = {"path": path, "params": params}
    return entries


def parse_notebook_params(qmd_path: Path) -> set[str]:
    """Extract declared params from a .qmd file's YAML front matter."""
    text = qmd_path.read_text()
    match = re.match(r"^---\n(.*?)\n---", text, re.DOTALL)
    if not match:
        return set()
    front_matter = yaml.safe_load(match.group(1))
    params = front_matter.get("params", {})
    if isinstance(params, dict):
        return set(params.keys())
    return set()


def main() -> int:
    registry_text = REGISTRY_PATH.read_text()
    # Run validation against both registry methods (create + analysis)
    entries = parse_registry(registry_text)

    errors = []
    for notebook_id, info in entries.items():
        # Resolve path relative to repo root; registry uses ${projectDir}/... prefix
        raw_path = info["path"].replace("${projectDir}/", "")
        qmd_path = Path(raw_path)
        if not qmd_path.exists():
            errors.append(f"[{notebook_id}] notebook not found: {qmd_path}")
            continue

        notebook_params = parse_notebook_params(qmd_path)
        for param in info["params"]:
            if param not in notebook_params:
                errors.append(
                    f"[{notebook_id}] registry param '{param}' not declared in {qmd_path}"
                )

    if errors:
        print("Notebook registry validation failed:\n")
        for err in errors:
            print(f"  {err}")
        return 1

    print(f"Registry OK — {len(entries)} notebook(s) validated.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
