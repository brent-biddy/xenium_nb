# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## What this pipeline does

`xenium_nb` is a Nextflow pipeline for Xenium spatial transcriptomics analysis with two entry points:

- `create.nf` — converts raw Xenium outputs into reusable SpatialData zarr artifacts
- `analyze.nf` — renders Quarto analysis notebooks against artifact samplesheets

## Commands

### Run create workflow
Both `--samplesheet` and `--create` are required. Valid `--create` values: `downsample`, `sdata`, `follicle_sdata`, `all`.

```bash
nextflow run create.nf --samplesheet assets/samplesheet.csv --create downsample
nextflow run create.nf --samplesheet assets/samplesheet.csv --create all
nextflow run create.nf --samplesheet assets/samplesheet.csv --create sdata
nextflow run create.nf --samplesheet results/sample_sdata_samplesheet.csv --create follicle_sdata
```

The `downsample` mode requires the samplesheet to include `xmin,ymin,xmax,ymax` columns (µm coordinates) and an optional `region_name` column. The `region_name` defaults to the sample ID if omitted.

### Run analyze workflow
Both `--samplesheet` and `--analyze` are required. Valid `--analyze` values: `plot_follicle`, `all`.

```bash
nextflow run analyze.nf --samplesheet results/follicle_sdata_samplesheet.csv --analyze plot_follicle
```

### Profiles
Defined in `nextflow.config`:

| Profile | Executor | Container |
|---------|----------|-----------|
| (none)  | local, no container | requires activated conda env with Quarto + notebook deps |
| `local` | local, Apptainer | `babiddy755/xenium_nb:20260505-66addc7`, 2 CPUs, 8 GB |
| `oscer` | SLURM on OSCER HPC, Apptainer | same image, 8 CPUs, memory retries 32→64→96 GB |

The `local` profile defaults `samplesheet` and `cell_ids_file` to the test assets, so no extra flags are needed:

```bash
nextflow run create.nf --create all -profile local
nextflow run analyze.nf --analyze plot_follicle -profile local
```

### Stub run (CI-equivalent, no script/notebook execution)
```bash
nextflow run create.nf -stub --create downsample -profile local
nextflow run create.nf -stub --create all -profile local
nextflow run analyze.nf -stub --samplesheet assets/ci_analyze_samplesheet.csv --analyze plot_follicle
```

### Validate notebook registry
```bash
python bin/check_notebook_registry.py
```

### Config parse check
```bash
nextflow config .
```

## Architecture

### Two-stage data flow
1. `create.nf` reads a `sample,path` samplesheet and runs `CREATE_SDATA` and/or `CREATE_FOLLICLE_SDATA` (both plain Python scripts) to produce zarr artifacts. It writes `results/sample_sdata_samplesheet.csv` and `results/follicle_sdata_samplesheet.csv` as handoff inputs for the next stage.
2. `analyze.nf` reads an artifact samplesheet (produced by create) and runs analysis notebook processes, publishing reports under `results/<sample>/<notebook_id>/`.

### Create-stage scripts (`bin/`)
`CREATE_SDATA` runs `bin/create_sdata.py` and `CREATE_FOLLICLE_SDATA` runs `bin/create_follicle_sdata.py`. Both use `argparse` CLIs; parameters are passed directly from the pipeline rather than through a params YAML.

### Notebook registry (`assets/notebook_registry.json`)
Maps analysis notebook IDs to their `.qmd` path and the params they declare. This is the source of truth used by `modules/quarto_params.nf` at runtime and validated by `bin/check_notebook_registry.py` in CI. Every param listed in the registry must have a matching variable in the notebook's `#| tags: [parameters]` cell. Create-stage scripts are not registered here.

### Params YAML flow (`modules/quarto_params.nf`)
Used by `analyze.nf` only. `paramsFile()` writes `<outdir>/.quarto_params/<notebook>/params_<id>.yml` and returns the path for Nextflow staging. Writing to `outdir` (NFS) rather than `/tmp` is intentional — symlinks to head-node `/tmp` break on OSCER compute nodes.

### Process conventions
- Always use `script:` blocks, never `exec:` — processes must run through SLURM.
- Every process script sets `XDG_CACHE_HOME=$PWD/.cache` and `TMPDIR=$PWD/tmp` to avoid writing to a read-only compute-node `/tmp`.
- Keep named input variables; do not inline maps into process call arguments.
- Build command lines with optional arguments using a Groovy list + conditional append:

```groovy
def myArgs = ["--required_a ${val_a}", "--required_b ${val_b}"]
if (optional_c) myArgs << "--optional_c ${optional_c}"
"""
my_script.py ${myArgs.join(' ')}
"""
```

## Adding an analysis notebook

1. Create `notebooks/analyze/<name>.qmd` with a `#| tags: [parameters]` Python cell declaring all inputs.
2. Add an entry to `assets/notebook_registry.json` with the notebook ID, relative path, and `params` list matching the parameters cell exactly.
3. Wire a new process into `modules/analyze_notebooks.nf` and `analyze.nf`.
4. Run `python bin/check_notebook_registry.py` to verify.

## Adding a create-stage script

Create-stage steps use plain Python scripts, not notebooks.

1. Create `bin/<name>.py` with an `argparse` CLI (`parse_args()` function) declaring all inputs.
2. Wire a new process into `modules/create_notebooks.nf` and `create.nf`, passing args directly.
3. No registry entry is needed.

## CI

Two GitHub Actions run on PRs to `main`:
- **Validate notebook registry** — runs `python bin/check_notebook_registry.py`
- **Stub run** — runs both entry points with `-stub` to verify workflow wiring without executing scripts or notebooks

## Code style (`.nf` files)

- 4-space indentation
- Process names in `UPPER_SNAKE_CASE`; params, variables, CSV headers in `snake_case`
- Add file-level header comments, docstrings on helper functions, section markers, and WHY comments for non-obvious decisions
- Annotate channel shape at every `.set {}` call and after non-obvious transformations so the tuple structure is always visible without tracing back through the chain:

```groovy
.set { createSdataInputs } // tuple(sample, staged_path, he_image, he_alignment)
// createSdataRun.artifacts: tuple(sample, zarr)
```
