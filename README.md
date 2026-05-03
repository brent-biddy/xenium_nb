# xenium_nb

A Nextflow pipeline for analysing 10x Xenium spatial transcriptomics data using Quarto notebooks. Each notebook runs as an independent Nextflow process, producing a rendered HTML report alongside any output files.

Notebook parameters are staged into each task work directory as `params.json` and loaded explicitly by the notebook code.

---

## Overview

The pipeline runs in three passes, each driven by its own samplesheet:

| Pass | Notebook | Input | Output |
|------|----------|-------|--------|
| 1 | `01_create_spatialdata` | Raw Xenium output directory | SpatialData zarr store |
| 2 | `02_subset_follicle` | Base results directory | Per-cell follicle zarr stores |
| 3 | `03_plot_follicle` | Base results directory | Rendered tissue image plots |

---

## Requirements

- [Nextflow](https://www.nextflow.io/) ≥ 23.0
- [Quarto](https://quarto.org/) ≥ 1.4
- Python packages: `spatialdata`, `spatialdata-io`, `spatialdata-plot`, `session-info`, `pyyaml`, `nbformat`

---

## Repository structure

```
xenium_nb/
├── main.nf                    # Pipeline entry point
├── nextflow.config            # Parameters and profiles
├── conf/
│   └── base.config            # Resource defaults
├── modules/
│   └── run_notebook.nf        # RUN_NOTEBOOK process
├── notebooks/
│   ├── 01_create_spatialdata.qmd
│   ├── 02_subset_follicle.qmd
│   └── 03_plot_follicle.qmd
├── bin/
│   ├── timer.py               # Timing utilities for notebooks
│   └── make_follicle_samplesheet.py
└── assets/
    ├── samplesheet.csv        # Sample-level samplesheet (pass 1 & 2)
    ├── follicle_samplesheet.csv  # Cell-level samplesheet (pass 3, generated)
    └── stage_quality_area_all_rois.csv  # Cell ID reference file
```

---

## Samplesheets

### Sample-level (`assets/samplesheet.csv`)

Used for passes 1 and 2. `path` points to the raw Xenium output directory for pass 1, and to the base results directory for pass 2. Relative paths are resolved against the directory where you launch `nextflow run`.

```csv
sample,path
ROI1,/path/to/ROI1/xenium_output
ROI2,/path/to/ROI2/xenium_output
```

### Follicle-level (`assets/follicle_samplesheet.csv`)

Used for pass 3. Generated automatically by `bin/make_follicle_samplesheet.py`. Each row is one annotated follicle cell; `sample` is `<ROI>_<cell_id>`, `roi_id` and `cell_id` are emitted explicitly, and `path` is the base results directory. Relative paths are resolved against the directory where you launch `nextflow run`.

```csv
sample,roi_id,cell_id,path
ROI1_aaaaimck-1,ROI1,aaaaimck-1,results
ROI1_aaaalpdj-1,ROI1,aaaalpdj-1,results
ROI2_aaabfpcg-1,ROI2,aaabfpcg-1,results
```

### Cell ID reference file (`assets/stage_quality_area_all_rois.csv`)

Maps cell IDs to samples. Must contain `Donor.ROI` and `cell_id` columns. An optional `radius` column sets a per-cell bounding box radius (µm); missing values fall back to `params.radius`.

---

## Usage

### Pass 1 — Create SpatialData zarr stores

```bash
nextflow run main.nf \
    --samplesheet assets/samplesheet.csv \
    --notebooks "[${PWD}/notebooks/01_create_spatialdata.qmd]"
```

### Pass 2 — Subset follicles

Update `assets/samplesheet.csv` so `path` points to the base results directory (e.g. `results/`), then run:

```bash
nextflow run main.nf \
    --samplesheet assets/samplesheet.csv \
    --notebooks "[${PWD}/notebooks/02_subset_follicle.qmd]"
```

### Generate follicle samplesheet

```bash
python bin/make_follicle_samplesheet.py \
    --cell-ids assets/stage_quality_area_all_rois.csv \
    --outdir   results \
    --output   assets/follicle_samplesheet.csv
```

### Pass 3 — Plot follicles

```bash
nextflow run main.nf \
    --samplesheet assets/follicle_samplesheet.csv \
    --notebooks "[${PWD}/notebooks/03_plot_follicle.qmd]"
```

---

## Configuration

Key parameters (set in `nextflow.config` or passed via `--param value`):

| Parameter | Default | Description |
|-----------|---------|-------------|
| `samplesheet` | `null` | Path to samplesheet CSV |
| `outdir` | `results` | Output directory |
| `cell_ids_file` | `assets/stage_quality_area_all_rois.csv` | Cell ID reference file |
| `radius` | `250` | Default bounding box radius (µm) |
| `notebooks` | `[01_create_spatialdata.qmd]` | Notebooks to run. Override per pass via `--notebooks` (see Usage). |

### Profiles

| Profile | Description |
|---------|-------------|
| (default) | Local execution, no container |
| `oscer` | SLURM executor on OSCER HPC, Apptainer container, scratch-based work directory and image cache. Memory scales 32→64→96 GB across retries. |

Activate with `-profile oscer`:

```bash
nextflow run main.nf \
    --samplesheet assets/samplesheet.csv \
    --notebooks "[${PWD}/notebooks/01_create_spatialdata.qmd]" \
    -profile oscer
```

---

## Output structure

```
results/
├── pipeline_info/
│   ├── timeline.html
│   └── report.html
├── ROI1/
│   ├── 01_create_spatialdata/
│   │   ├── ROI1_01_create_spatialdata.html
│   │   └── output/
│   │       └── ROI1.zarr/
│   ├── 02_subset_follicle/
│   │   ├── ROI1_02_subset_follicle.html
│   │   └── output/
│   │       ├── aaaaimck-1.zarr/
│   │       └── aaaalpdj-1.zarr/
│   └── 03_plot_follicle/
│       ├── ROI1_aaaaimck-1_03_plot_follicle.html
│       └── ROI1_aaaalpdj-1_03_plot_follicle.html
└── ROI2/
    └── ...
```

---

## Adding notebooks

1. Create a new `.qmd` file in `notebooks/` with a YAML params block declaring `sample`, `path`, and any additional params needed.
2. Add the notebook path to `params.notebooks` in `nextflow.config` or pass it via `--notebooks` on the CLI.
3. Any params not declared in the notebook's front matter are automatically filtered out before rendering.
