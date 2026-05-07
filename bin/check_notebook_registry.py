#!/usr/bin/env python3
# Validates that every param listed in assets/notebook_registry.json exists
# in the corresponding notebook's YAML front matter. Run from the repo root.

import json
import re
import sys
from pathlib import Path

REGISTRY_PATH = Path("assets/notebook_registry.json")


def parse_notebook_params(qmd_path: Path) -> set[str]:
    """Extract declared params from the #| tags: [parameters] cell in a .qmd file."""
    text = qmd_path.read_text()
    match = re.search(
        r"```\{python\}.*?#\|\s*tags:\s*\[parameters\](.*?)```",
        text,
        re.DOTALL,
    )
    if not match:
        return set()
    return {m.group(1) for m in re.finditer(r"^(\w+)\s*=", match.group(1), re.MULTILINE)}


def main() -> int:
    registry = json.loads(REGISTRY_PATH.read_text())

    entries = registry

    errors = []
    for notebook_id, info in entries.items():
        qmd_path = Path(info["path"])
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
