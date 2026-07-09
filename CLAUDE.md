# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## What this pipeline does

`xenium_nb` is a Nextflow pipeline for Xenium spatial transcriptomics analysis. All steps run through a single entry point, `main.nf`, selected with `--step`:

- `downsample_xenium_region` — crops a raw Xenium output directory to a bounding box region
- `create_sdata` — converts raw Xenium output into a sample-level SpatialData zarr artifact
- `create_follicle_sdata` — subsets a sample zarr into per-cell follicle zarrs
- `cluster_sdata` / `cluster_sdata_gpu` — QC, normalize, PCA, neighbors, UMAP, Leiden clustering (CPU vs. RAPIDS/GPU)
- `cluster_sdata_gpu_ooc` — same clustering pipeline, but streams the table's X matrix through Dask (rapids-singlecell out-of-core) so tables too large for VRAM (e.g. a merged cohort from `concat_sdata`) can still run on a single GPU
- `concat_sdata` — merges multiple sample zarrs into one
- `downsample_sdata` — subsamples cells from a SpatialData zarr
- `plot_follicle` — renders the `plot_follicle.qmd` Quarto notebook per follicle zarr

## Commands

### Run a step
`--samplesheet` is always required; columns vary by step (see `main.nf`'s header comment for the full table). Some steps take extra flags.

```bash
nextflow run main.nf --step downsample_xenium_region --samplesheet assets/samplesheet.csv
nextflow run main.nf --step create_sdata --samplesheet assets/downsampled_region_samplesheet.csv
nextflow run main.nf --step create_follicle_sdata --samplesheet my_sample_zarrs.csv --cell_ids_file assets/stage_quality_area_all_rois.csv
nextflow run main.nf --step cluster_sdata --samplesheet my_sample_zarrs.csv
nextflow run main.nf --step cluster_sdata_gpu --samplesheet my_sample_zarrs.csv
nextflow run main.nf --step cluster_sdata_gpu_ooc --samplesheet my_sample_zarrs.csv --chunk_size 20000 --n_top_genes 2000
nextflow run main.nf --step concat_sdata --samplesheet assets/concat_sdata_samplesheet.csv
nextflow run main.nf --step downsample_sdata --samplesheet my_sample_zarrs.csv --fraction 0.1
nextflow run main.nf --step plot_follicle --samplesheet assets/ci_analyze_samplesheet.csv
```

`my_sample_zarrs.csv` above is a stand-in for a `sample,path` CSV pointing at a prior step's output zarrs (e.g. `results/<sample>/create_sdata/<sample>.zarr`). Some producing steps now publish a ready-to-use handoff samplesheet next to their outputs (`<outdir>/create_sdata_samplesheet.csv`, `<outdir>/cluster_sdata_samplesheet.csv`) that you can point the next step's `--samplesheet` straight at; for steps without one yet, hand-build the CSV.

`downsample_xenium_region` requires the samplesheet to include `xmin,ymin,xmax,ymax` columns (µm coordinates) and an optional `region_name` column, which defaults to the sample ID if omitted. `downsample_sdata` requires `--fraction` or `--n_cells`.

### Profiles
Defined in `nextflow.config`:

| Profile | Executor | Container |
|---------|----------|-----------|
| (none)  | local, no container | requires activated conda env with Quarto + notebook deps |
| `local` | local, Apptainer | `babiddy755/python_spatial:1.2.0`, 8 CPUs, 16 GB |
| `oscer` | SLURM on OSCER HPC, Apptainer | same image, 16 CPUs, memory retries 48→96→144 GB (heavier for `CONCAT_SDATA`/`CLUSTER_SDATA`); GPU steps use the `sooner_gpu_test` partition with `--gres=gpu:1 --nv` |

**Run directories.** The `local` and `oscer` profiles set their own `workDir` and `outdir` so nothing lands in the repo (runs are typically launched from the repo root). Each run gets one self-contained directory, `<out_root>/<run_id>/{work,results}`, so a whole run is a single unit to size (`du -sh`) or prune (`rm -rf`). The shared Apptainer cache is a sibling of the run dirs (`<out_root>/apptainer_cache`), never nested under a `run_id`, so it survives across runs:

- `local` → `~/xenium_nb_out/<run_id>/{work,results}`
- `oscer` → `/scratch/$USER/xenium_nb_out/<run_id>/{work,results}`

Keeping `work` and `results` under the same root also keeps them on one filesystem, which the modules' hardlink publishing (`mode: 'link'`) relies on to avoid a second full copy of each zarr.

Because `run_id` defaults to a fresh timestamp, `-resume` across separate launches only works if you pin the id with `--run_id <name>` (or recover the prior timestamp from the run dir name / `.nextflow.log` and pass it back). `-resume` must also be run from the same launch directory, since its cache lives in `.nextflow/` there.

The `local` profile defaults `samplesheet` and `cell_ids_file` to the test assets, and also points `cluster_sdata_gpu` / `cluster_sdata_gpu_ooc` at the local RAPIDS container with WSL2 GPU passthrough settings:

```bash
nextflow run main.nf --step cluster_sdata_gpu -profile local
nextflow run main.nf --step cluster_sdata_gpu_ooc -profile local
```

`cluster_sdata_gpu_ooc` additionally needs `dask` and `zarr` in the container — both are present in `babiddy755/python_spatial:1.2.0` as rapids-singlecell/spatialdata dependencies (verified).

### Stub run (CI-equivalent, no script/notebook execution)
```bash
nextflow run main.nf --step create_sdata -stub --samplesheet assets/samplesheet.csv
nextflow run main.nf --step create_follicle_sdata -stub --samplesheet assets/ci_analyze_samplesheet.csv
nextflow run main.nf --step cluster_sdata -stub --samplesheet assets/ci_analyze_samplesheet.csv
nextflow run main.nf --step cluster_sdata_gpu -stub --samplesheet assets/ci_analyze_samplesheet.csv
nextflow run main.nf --step cluster_sdata_gpu_ooc -stub --samplesheet assets/ci_analyze_samplesheet.csv
nextflow run main.nf --step concat_sdata -stub --samplesheet assets/ci_analyze_samplesheet.csv
nextflow run main.nf --step downsample_sdata -stub --samplesheet assets/ci_analyze_samplesheet.csv --fraction 0.1
nextflow run main.nf --step downsample_xenium_region -stub --samplesheet assets/samplesheet.csv
nextflow run main.nf --step plot_follicle -stub --samplesheet assets/ci_analyze_samplesheet.csv
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

### Single entry point, one workflow per step
`main.nf` dispatches on `--step` to one of eight named workflows, each of which reads a samplesheet, builds a channel of tuples, and pipes it into a single process. There is no chaining between steps inside Nextflow — to run steps in sequence, point the next step's `--samplesheet` at a CSV listing the prior step's published output paths (e.g. `results/<sample>/create_sdata/<sample>.zarr`). As a convenience (not a control-flow link), every zarr-producing step publishes a handoff samplesheet into `outdir` (`<step>_samplesheet.csv`) that you can feed directly to the next step. Each process emits a `samplesheet_row` output whose published path comes from a per-module helper that also drives that module's `publishDir` — so the convention is single-sourced in the module and `main.nf` just `.map { it.text }` + `collectFile`s the rows (the `.text` read makes `collectFile`'s `sort` deterministic). The row fragment is kept out of the publish dir via `publishDir`'s `saveAs`. Schemas: the single-per-sample producers (`create_sdata`, `cluster_sdata`, `cluster_sdata_gpu`, `cluster_sdata_gpu_ooc`, `downsample_sdata`) and `concat_sdata` emit `sample,path`; `create_follicle_sdata` emits `sample,cell,path` (one row per per-cell zarr) for `plot_follicle`. `concat_sdata` and `create_follicle_sdata` also stage `.zarr` inputs into the work dir, so their row generation globs `*.zarr` and excludes the staged inputs.

### Create/cluster/downsample scripts (`bin/`)
Every step except `plot_follicle` runs a plain Python script with an `argparse` CLI (`bin/<step>.py`), invoked directly from its module's `script:` block — no params YAML involved.

### Notebook registry (`assets/notebook_registry.json`)
Maps analysis notebook IDs (currently just `plot_follicle`) to their `.qmd` path and the params they declare. This is the source of truth used by `modules/quarto_params.nf` at runtime and validated by `bin/check_notebook_registry.py` in CI. Every param listed in the registry must have a matching variable in the notebook's `#| tags: [parameters]` cell. The Python scripts under `bin/` are not registered here.

### Params YAML flow (`modules/quarto_params.nf`)
Used by the `plot_follicle` step only. `paramsFile()` writes `<outdir>/.quarto_params/<notebook>/params_<id>.yml` and returns the path for Nextflow staging. Writing to `outdir` (NFS) rather than `/tmp` is intentional — symlinks to head-node `/tmp` break on OSCER compute nodes.

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
3. Wire a new process into `modules/<name>.nf` and add a matching `--step` branch in `main.nf`.
4. Run `python bin/check_notebook_registry.py` to verify.

## Adding a create/cluster/downsample-stage script

These steps use plain Python scripts, not notebooks.

1. Create `bin/<name>.py` with an `argparse` CLI (`parse_args()` function) declaring all inputs.
2. Wire a new process into `modules/<name>.nf` and add a matching `--step` branch in `main.nf`, passing args directly.
3. No registry entry is needed.

## CI

Two GitHub Actions run on PRs to `main`:
- **Validate notebook registry** — runs `python bin/check_notebook_registry.py`
- **Stub run** — runs every `main.nf --step` with `-stub` to verify workflow wiring without executing scripts or notebooks

## Code style (`.nf` files)

- 4-space indentation
- Process names in `UPPER_SNAKE_CASE`; params, variables, CSV headers in `snake_case`
- Add file-level header comments, docstrings on helper functions, section markers, and WHY comments for non-obvious decisions
- Annotate channel shape at every `.set {}` call and after non-obvious transformations so the tuple structure is always visible without tracing back through the chain:

```groovy
.set { createSdataInputs } // tuple(sample, staged_path, he_image, he_alignment)
// createSdataRun.artifacts: tuple(sample, zarr)
```
