# xenium_nb

A Nextflow pipeline for Xenium spatial transcriptomics data. Each module can be run individually or chained together via `create.nf` and `analyze.nf`.

---

## Requirements

- [Nextflow](https://www.nextflow.io/) ≥ 25.10.0
- For default non-container local runs: [Quarto](https://quarto.org/) ≥ 1.4 and the required Python notebook packages
- For containerized local runs: Apptainer

---

## Repository structure

```
xenium_nb/
├── create.nf                  # Chains create-stage modules into multi-step workflows
├── analyze.nf                 # Chains analysis modules into multi-step workflows
├── nextflow.config            # Parameters and profiles
├── modules/
│   ├── create_sdata.nf            # Raw Xenium → SpatialData zarr
│   ├── create_follicle_sdata.nf   # Sample zarr → per-cell follicle zarrs
│   ├── cluster_sdata.nf           # QC, normalize, PCA, UMAP, Leiden clustering
│   ├── concat_sdata.nf            # Merge multiple SpatialData zarrs into one
│   ├── downsample_xenium_region.nf # Crop raw Xenium output to a bounding box
│   ├── plot_follicle.nf           # Per-cell follicle plots (Quarto notebook)
│   └── quarto_params.nf           # Quarto params YAML helper
├── notebooks/
│   ├── README.md
│   └── analyze/
│       └── plot_follicle.qmd
├── bin/
│   ├── create_sdata.py              # Convert raw Xenium output to SpatialData zarr
│   ├── create_follicle_sdata.py     # Subset sample zarr into per-cell follicle zarrs
│   ├── cluster_sdata.py             # QC, normalize, cluster, write zarr
│   ├── concat_sdata.py              # Concatenate SpatialData zarrs
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

### Running modules individually

Each module can be run directly as an entry point. This is the typical way to run a single step:

```bash
nextflow run modules/create_sdata.nf \
    --samplesheet assets/samplesheet.csv

nextflow run modules/cluster_sdata.nf \
    --samplesheet results/sample_sdata_samplesheet.csv

nextflow run modules/concat_sdata.nf \
    --samplesheet results/sample_sdata_samplesheet.csv

nextflow run modules/plot_follicle.nf \
    --samplesheet results/follicle_sdata_samplesheet.csv
```

Samplesheet columns required by each module:

| Module | Samplesheet columns |
|--------|---------------------|
| `create_sdata` | `sample, path[, he_image, he_alignment]` |
| `create_follicle_sdata` | `sample, path` (+ `--cell_ids_file`) |
| `cluster_sdata` | `sample, path` |
| `concat_sdata` | `path` |
| `downsample_xenium_region` | `sample, path, xmin, ymin, xmax, ymax[, region_name, he_image, he_alignment]` |
| `plot_follicle` | `sample, cell, path` |

### Chaining with create.nf

`--create` accepts: `sdata`, `follicle_sdata`, `all`, `downsample`, `concat`.

```bash
# Raw Xenium → sample zarrs → follicle zarrs (full create chain)
nextflow run create.nf \
    --samplesheet assets/samplesheet.csv \
    --create all

# Sample zarrs only
nextflow run create.nf \
    --samplesheet assets/samplesheet.csv \
    --create sdata

# Follicle zarrs from existing sample zarrs
nextflow run create.nf \
    --samplesheet results/sample_sdata_samplesheet.csv \
    --create follicle_sdata

# Crop raw Xenium to a bounding box region
nextflow run create.nf \
    --samplesheet assets/samplesheet.csv \
    --create downsample

# Merge multiple sample zarrs into one
nextflow run create.nf \
    --samplesheet results/sample_sdata_samplesheet.csv \
    --create concat
```

`--create sdata` and `--create all` write handoff samplesheets:
- `results/sample_sdata_samplesheet.csv`
- `results/follicle_sdata_samplesheet.csv`

### Chaining with analyze.nf

`--analyze` accepts: `plot_follicle`, `cluster_sdata`, `all`.

```bash
nextflow run analyze.nf \
    --samplesheet results/follicle_sdata_samplesheet.csv \
    --analyze plot_follicle

nextflow run analyze.nf \
    --samplesheet results/sample_sdata_samplesheet.csv \
    --analyze cluster_sdata
```

---

## Configuration

Key parameters (set in `nextflow.config` or passed via `--param value`):

| Parameter | Default | Description |
|-----------|---------|-------------|
| `samplesheet` | *(required)* | Path to samplesheet CSV |
| `outdir` | `results` | Output directory |
| `cell_ids_file` | `assets/stage_quality_area_all_rois.csv` | Cell ID reference file for `create_follicle_sdata` |
| `radius` | `100` | Default bounding box radius (µm) around each cell centroid |
| `create` | *(required for create.nf)* | `sdata`, `follicle_sdata`, `all`, `downsample`, or `concat` |
| `analyze` | *(required for analyze.nf)* | `plot_follicle`, `cluster_sdata`, or `all` |

### Profiles

| Profile | Description |
|---------|-------------|
| (default) | Local execution, no container (use an activated conda env that provides the notebook kernel). |
| `local` | Local execution with Apptainer, sized for a laptop / WSL2 box (2 CPUs, 8 GB). Defaults `samplesheet` and `cell_ids_file` to the test assets. Also points `CLUSTER_SDATA_GPU` at the local RAPIDS container with WSL2 GPU passthrough settings. |
| `oscer` | SLURM executor on OSCER HPC, Apptainer container, scratch-based work directory. Memory scales 32→64→96 GB across retries. |

```bash
# Local profile (no --samplesheet needed)
nextflow run create.nf --create all -profile local
nextflow run analyze.nf --analyze plot_follicle -profile local

# OSCER HPC
nextflow run create.nf \
    --samplesheet assets/samplesheet.csv \
    --create all \
    -profile oscer
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
│   ├── cluster_sdata/
│   │   └── clustered.zarr/
│   └── plot_follicle/
│       ├── aaaaimck-1_plot_follicle.pptx
│       └── aaameida-1_plot_follicle.pptx
├── concat_sdata/
│   └── merged.zarr/
└── ROI1_B/
    └── ...
```
