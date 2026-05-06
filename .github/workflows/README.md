# CI Actions

## Validate notebook registry (`validate_registry.yml`)

Runs on push to `main` and on all PRs targeting `main`. Can also be triggered manually via `workflow_dispatch`.

Checks that every param listed in [`lib/NotebookRegistry.groovy`](../../lib/NotebookRegistry.groovy) exists in the corresponding notebook's YAML front matter. Fails if a registry param is missing from the notebook, which indicates the registry and notebook have drifted out of sync.

The validation script lives at [`bin/check_notebook_registry.py`](../../bin/check_notebook_registry.py).

## Stub run (`stub_run.yml`)

Runs on all PRs targeting `main`. Can also be triggered manually via `workflow_dispatch`.

Executes `create.nf` and `analyze.nf` in Nextflow stub mode (`-stub`), which skips actual process execution but verifies that the workflow parses correctly, all `include` statements resolve, and the process wiring is valid. Uses placeholder samplesheets from `assets/`.
