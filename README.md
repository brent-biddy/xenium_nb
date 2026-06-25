# xenium_nb

A Nextflow pipeline for Xenium spatial transcriptomics data with two separate entry points:

- `create.nf` builds reusable data artifacts
- `analyze.nf` runs analysis notebooks against artifact samplesheets

Notebook-specific inputs, outputs, and samplesheet contracts are documented in [notebooks/README.md](notebooks/README.md).

---

## Requirements

- [Nextflow](https://www.nextflow.io/) ≥ 25.10.0
- For default non-container local runs: [Quarto](https://quarto.org/) ≥ 1.4 and the required Python notebook packages
- For containerized local runs: Apptainer

---

## Repository structure

```
xenium_nb/
├── create.nf                  # Create workflow: raw Xenium -> sample and follicle artifacts
├── analyze.nf                 # Analysis workflow: artifact samplesheet -> notebook reports
├── nextflow.config            # Parameters and profiles
├── modules/
│   ├── create_notebooks.nf    # Create-stage process definitions
│   ├── analyze_notebooks.nf   # Analysis notebook process definitions
│   └── quarto_params.nf       # Quarto params YAML helper (used by analyze.nf)
├── notebooks/
│   ├── README.md
│   └── analyze/
│       └── plot_follicle.qmd
├── bin/
│   ├── create_sdata.py              # Convert raw Xenium output to SpatialData zarr
│   ├── create_follicle_sdata.py     # Subset sample zarr into per-cell follicle zarrs
│   ├── downsample_xenium.py         # Downsample a full Xenium output directory
│   ├── downsample_xenium_region.py  # Crop a Xenium output to a bounding box region
│   ├── check_notebook_registry.py   # CI validator for notebook registry
│   └── timer.py                     # Timing utilities for scripts and notebooks
└── assets/
    ├── samplesheet.csv                    # Sample-level samplesheet
    ├── test_cell_ids.csv                  # Minimal cell ID file for test runs
    ├── stage_quality_area_all_rois.csv    # Full cell ID reference file
    └── notebook_registry.json             # Notebook metadata (paths and declared params)
```

---

## Usage

### Create artifacts

```bash
nextflow run create.nf \
    --samplesheet assets/samplesheet.csv \
    --create all
```

This runs both create-stage scripts:

- `bin/create_sdata.py` — raw Xenium → sample-level SpatialData zarr
- `bin/create_follicle_sdata.py` — sample zarr → per-cell follicle zarrs

and writes:

- sample zarrs under `results/<sample>/create_sdata/output/`
- follicle zarrs under `results/<sample>/follicle_sdata/output/`
- `results/sample_sdata_samplesheet.csv`
- `results/follicle_sdata_samplesheet.csv`

### Create sdata only

```bash
nextflow run create.nf \
    --samplesheet assets/samplesheet.csv \
    --create sdata
```

Writes `results/sample_sdata_samplesheet.csv`, which can be used as input to `--create follicle_sdata`.

### Create follicle sdata only

Run this after `create sdata` has produced `results/sample_sdata_samplesheet.csv`:

```bash
nextflow run create.nf \
    --samplesheet results/sample_sdata_samplesheet.csv \
    --create follicle_sdata
```

Writes `results/follicle_sdata_samplesheet.csv`.

### Analyze follicle artifacts

```bash
nextflow run analyze.nf \
    --samplesheet results/follicle_sdata_samplesheet.csv \
    --analyze plot_follicle
```

---

## Configuration

Key parameters (set in `nextflow.config` or passed via `--param value`):

| Parameter | Default | Description |
|-----------|---------|-------------|
| `samplesheet` | *(required)* | Path to samplesheet CSV |
| `outdir` | `results` | Output directory |
| `cell_ids_file` | `${projectDir}/assets/stage_quality_area_all_rois.csv` | Cell ID reference file path |
| `radius` | `100` | Default bounding box radius (µm) |
| `create` | *(required)* | Create workflow mode: `sdata`, `follicle_sdata`, or `all` |
| `analyze` | *(required)* | Analysis notebook selector: `all` or a notebook ID from `assets/notebook_registry.json` |

Analysis notebook IDs and how to add new notebooks are documented in [notebooks/README.md](notebooks/README.md).

### Profiles

| Profile | Description |
|---------|-------------|
| (default) | Local execution, no container (use an activated conda env that provides the notebook kernel). |
| `test` | Local execution with Apptainer, sized for a laptop / WSL2 box (2 CPUs, 8 GB). Defaults `samplesheet` and `cell_ids_file` to the test assets. |
| `oscer` | SLURM executor on OSCER HPC, Apptainer container, scratch-based work directory. Memory scales 32→64→96 GB across retries. |

Activate with `-profile oscer`:

```bash
nextflow run create.nf \
    --samplesheet assets/samplesheet.csv \
    --create all \
    -profile oscer
```

For local containerized runs:

```bash
nextflow run create.nf --create all -profile test
nextflow run analyze.nf --analyze plot_follicle -profile test
```

---

## Output structure

```
results/
├── pipeline_info/
│   ├── timeline.html
│   └── report.html
├── sample_sdata_samplesheet.csv
├── follicle_sdata_samplesheet.csv
├── ROI1_A/
│   ├── create_sdata/
│   │   └── output/
│   │       └── ROI1_A.zarr/
│   ├── follicle_sdata/
│   │   └── output/
│   │       ├── aaaaimck-1.zarr/
│   │       └── aaameida-1.zarr/
│   └── plot_follicle/
│       ├── aaaaimck-1_plot_follicle.pptx
│       ├── aaaaimck-1_plot_follicle.timing.tsv
│       ├── aaameida-1_plot_follicle.pptx
│       └── aaameida-1_plot_follicle.timing.tsv
└── ROI1_B/
    └── ...
```

If `--create sdata` is used, `follicle_sdata/` outputs and `results/follicle_sdata_samplesheet.csv` are not created.

The sample-stage sheet is the handoff input for `--create follicle_sdata`.
