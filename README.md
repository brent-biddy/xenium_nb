# xenium_nb

A Nextflow pipeline for Xenium spatial transcriptomics data with two separate entry points:

- `build.nf` builds reusable data artifacts
- `analyze.nf` runs analysis notebooks against artifact samplesheets

Both workflows use the same two-column samplesheet contract:

```csv
sample,path
```

Notebook parameters are staged into each task work directory as `params.json` and loaded explicitly by the notebook code.

---

## Requirements

- [Nextflow](https://www.nextflow.io/) ≥ 23.0
- [Quarto](https://quarto.org/) ≥ 1.4
- Python packages: `spatialdata`, `spatialdata-io`, `spatialdata-plot`, `session-info`, `pyyaml`, `nbformat`

---

## Repository structure

```
xenium_nb/
├── build.nf                   # Build workflow: raw Xenium -> sample and follicle artifacts
├── analyze.nf                 # Analysis workflow: artifact samplesheet -> notebook reports
├── nextflow.config            # Parameters and profiles
├── conf/
│   └── base.config            # Resource defaults
├── modules/
│   ├── create_spatialdata.nf  # Sample-level artifact producer
│   ├── subset_follicle.nf     # Follicle-level artifact producer
│   ├── run_notebook.nf        # Generic analysis notebook runner
│   └── write_samplesheet.nf   # Writes two-column artifact samplesheets
├── notebooks/
│   ├── create_spatialdata.qmd
│   ├── subset_follicle.qmd
│   └── plot_follicle.qmd
├── bin/
│   ├── timer.py               # Timing utilities for notebooks
│   └── make_follicle_samplesheet.py  # Legacy helper for manual/export workflows
└── assets/
    ├── samplesheet.csv        # Sample-level samplesheet
    └── stage_quality_area_all_rois.csv  # Cell ID reference file
```

---

## Samplesheets

### Build workflow input

Used by `build.nf`. `path` points to a raw Xenium output directory.

```csv
sample,path
ROI1,/path/to/ROI1/xenium_output
ROI2,/path/to/ROI2/xenium_output
```

### Analysis workflow input

Used by `analyze.nf`. `path` points to an already-built artifact, either sample-level or follicle-level depending on the notebook scope.

Sample artifact sheet:

```csv
sample,path
ROI1,results/ROI1/create_spatialdata/output/ROI1.zarr
ROI2,results/ROI2/create_spatialdata/output/ROI2.zarr
```

Follicle artifact sheet:

```csv
sample,path
ROI1_aaaaimck-1,results/ROI1/subset_follicle/output/aaaaimck-1.zarr
ROI1_aaaalpdj-1,results/ROI1/subset_follicle/output/aaaalpdj-1.zarr
```

### Cell ID reference file

`assets/stage_quality_area_all_rois.csv` is used by the build workflow to decide which follicles to subset. It must contain `Donor.ROI` and `cell_id`. An optional `radius` column sets a per-cell bounding box radius (µm); missing values fall back to `params.radius`.

---

## Usage

### Build artifacts

```bash
nextflow run build.nf \
    --samplesheet assets/samplesheet.csv
```

By default this runs both producer notebooks:

- `create_spatialdata.qmd`
- `subset_follicle.qmd`

and writes:

- sample zarrs under `results/<sample>/create_spatialdata/output/`
- follicle zarrs under `results/<sample>/subset_follicle/output/`
- `results/pipeline_info/sample_artifacts.csv`
- `results/pipeline_info/follicle_artifacts.csv`

### Build sample artifacts only

```bash
nextflow run build.nf \
    --samplesheet assets/samplesheet.csv \
    --run_subset_follicle false
```

This runs only `create_spatialdata.qmd` and writes `results/pipeline_info/sample_artifacts.csv`.

### Analyze follicle artifacts

```bash
nextflow run analyze.nf \
    --samplesheet results/pipeline_info/follicle_artifacts.csv \
    --notebooks plot_follicle
```

### Analyze sample artifacts

When you add sample-scoped analysis notebooks to the registry, point `analyze.nf` at the sample artifact sheet:

```bash
nextflow run analyze.nf \
    --samplesheet results/pipeline_info/sample_artifacts.csv \
    --notebooks your_sample_notebook_id
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
| `run_subset_follicle` | `true` | Whether `build.nf` should run `subset_follicle.qmd` after building sample-level zarrs |
| `producer_registry` | built-in map | The two producer notebooks used by `build.nf` |
| `analysis_notebook_registry` | built-in map | Notebook IDs, paths, and scopes used by `analyze.nf` |
| `notebooks` | `[]` | Analysis notebook IDs to run in `analyze.nf` |

### Analysis notebook IDs

The built-in analysis registry currently defines:

| ID | Scope | Notebook |
|----|-------|----------|
| `plot_follicle` | `follicle` | `notebooks/plot_follicle.qmd` |

### Profiles

| Profile | Description |
|---------|-------------|
| (default) | Local execution, no container |
| `oscer` | SLURM executor on OSCER HPC, Apptainer container, scratch-based work directory and image cache. Memory scales 32→64→96 GB across retries. |

Activate with `-profile oscer`:

```bash
nextflow run build.nf \
    --samplesheet assets/samplesheet.csv \
    -profile oscer
```

---

## Output structure

```
results/
├── pipeline_info/
│   ├── timeline.html
│   ├── report.html
│   ├── sample_artifacts.csv
│   └── follicle_artifacts.csv
├── ROI1/
│   ├── create_spatialdata/
│   │   ├── ROI1_create_spatialdata.html
│   │   └── output/
│   │       └── ROI1.zarr/
│   ├── subset_follicle/
│   │   ├── ROI1_subset_follicle.html
│   │   └── output/
│   │       ├── aaaaimck-1.zarr/
│   │       └── aaaalpdj-1.zarr/
│   └── plot_follicle/
│       ├── ROI1_aaaaimck-1_plot_follicle.html
│       └── ROI1_aaaalpdj-1_plot_follicle.html
└── ROI2/
    └── ...
```

Analysis outputs also publish under the parent sample directory, so follicle reports from `ROI1_aaaaimck-1`, `ROI1_aaaalpdj-1`, and similar artifacts all land under `results/ROI1/plot_follicle/`.

If `--run_subset_follicle false` is used, `subset_follicle/` outputs and `follicle_artifacts.csv` are not created.

---

## Adding notebooks

1. Create a new `.qmd` file in `notebooks/` with a YAML params block declaring at least `sample` and `path`.
2. If it is a build-stage producer, register it in `params.producer_registry` and wire it into `build.nf`.
3. If it is an analysis notebook, register it in `params.analysis_notebook_registry` with a unique ID and a scope of `sample` or `follicle`.
4. Run analysis notebooks with `nextflow run analyze.nf --samplesheet <artifact_sheet.csv> --notebooks <id1,id2>`.
5. Any params not declared in the notebook's front matter are automatically filtered out before rendering.
