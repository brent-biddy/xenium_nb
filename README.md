# xenium_nb

A Nextflow pipeline for Xenium spatial transcriptomics data. All steps run through a single entry point, `main.nf`, selected with `--step`.

---

## Requirements

- [Nextflow](https://www.nextflow.io/) ≥ 25.10.0
- For default non-container local runs: [Quarto](https://quarto.org/) ≥ 1.4 and the required Python notebook packages
- For containerized local runs: Apptainer

---

## Repository structure

```
xenium_nb/
├── main.nf                        # Single entry point; dispatches on --step
├── nextflow.config                # Parameters and profiles
├── modules/
│   ├── downsample_xenium_region.nf  # Crop raw Xenium output to a bounding box
│   ├── create_sdata.nf              # Raw Xenium → sample-level SpatialData zarr
│   ├── create_follicle_sdata.nf     # Sample zarr → per-cell follicle zarrs
│   ├── cluster_sdata.nf             # QC, normalize, PCA, UMAP, Leiden clustering (CPU)
│   ├── cluster_sdata_gpu.nf         # Same clustering pipeline, RAPIDS-accelerated (GPU)
│   ├── concat_sdata.nf              # Merge multiple SpatialData zarrs into one
│   ├── downsample_sdata.nf          # Subsample cells from a SpatialData zarr
│   ├── plot_follicle.nf             # Per-cell follicle plots (Quarto notebook)
│   └── quarto_params.nf             # Quarto params YAML helper (used by plot_follicle)
├── notebooks/
│   ├── README.md
│   └── analyze/
│       └── plot_follicle.qmd
├── bin/
│   ├── downsample_xenium_region.py  # Crop a Xenium output to a bounding box region
│   ├── create_sdata.py              # Convert raw Xenium output to SpatialData zarr
│   ├── create_follicle_sdata.py     # Subset sample zarr into per-cell follicle zarrs
│   ├── cluster_sdata.py             # QC, normalize, cluster, write zarr (CPU)
│   ├── cluster_sdata_gpu.py         # QC, normalize, cluster, write zarr (RAPIDS/GPU)
│   ├── concat_sdata.py              # Concatenate SpatialData zarrs
│   ├── downsample_sdata.py          # Subsample a SpatialData zarr
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

Every step is run the same way:

```bash
nextflow run main.nf --step <name> --samplesheet <path> [step-specific flags]
```

| Step | Samplesheet columns | Extra flags |
|------|----------------------|-------------|
| `downsample_xenium_region` | `sample, path, xmin, ymin, xmax, ymax[, region_name, he_image, he_alignment]` | |
| `create_sdata` | `sample, path[, he_image, he_alignment]` | |
| `create_follicle_sdata` | `sample, path` | `--cell_ids_file <path>` (required) |
| `cluster_sdata` | `sample, path` | |
| `cluster_sdata_gpu` | `sample, path` | |
| `concat_sdata` | `path` | |
| `downsample_sdata` | `sample, path` | `--fraction <float>` or `--n_cells <int>` (one required) |
| `plot_follicle` | `sample, cell, path` | |

Examples:

```bash
# Crop raw Xenium output to a bounding box region
nextflow run main.nf --step downsample_xenium_region \
    --samplesheet assets/samplesheet.csv

# Raw Xenium → sample-level SpatialData zarr
nextflow run main.nf --step create_sdata \
    --samplesheet assets/downsampled_region_samplesheet.csv

# Sample zarr → per-cell follicle zarrs
nextflow run main.nf --step create_follicle_sdata \
    --samplesheet results/sample_sdata_samplesheet.csv \
    --cell_ids_file assets/stage_quality_area_all_rois.csv

# Cluster a sample zarr (CPU)
nextflow run main.nf --step cluster_sdata \
    --samplesheet results/sample_sdata_samplesheet.csv

# Cluster a sample zarr (GPU, RAPIDS)
nextflow run main.nf --step cluster_sdata_gpu \
    --samplesheet results/sample_sdata_samplesheet.csv

# Merge multiple sample zarrs into one
nextflow run main.nf --step concat_sdata \
    --samplesheet assets/concat_sdata_samplesheet.csv

# Subsample a sample zarr
nextflow run main.nf --step downsample_sdata \
    --samplesheet results/sample_sdata_samplesheet.csv \
    --fraction 0.1

# Render per-cell follicle plots
nextflow run main.nf --step plot_follicle \
    --samplesheet results/follicle_sdata_samplesheet.csv
```

`create_sdata` writes zarrs under `results/<sample>/create_sdata/output/` — point a later step's `--samplesheet` at those paths (see [Output structure](#output-structure)) to chain steps together.

---

## Configuration

Key parameters (set in `nextflow.config` or passed via `--param value`):

| Parameter | Default | Description |
|-----------|---------|-------------|
| `step` | *(required)* | Pipeline step to run; see the table above |
| `samplesheet` | *(required)* | Path to samplesheet CSV (columns vary by step) |
| `outdir` | `results` | Output directory |
| `cell_ids_file` | `${launchDir}/assets/stage_quality_area_all_rois.csv` | Cell ID reference file for `create_follicle_sdata` |
| `radius` | `100` | Default bounding box radius (µm) around each cell centroid; overridable per-cell via a `radius` column in `cell_ids_file` |
| `fraction` | `0.1` | Fraction of cells to retain in `downsample_sdata` |
| `n_cells` | `null` | Absolute cell count to retain in `downsample_sdata` (alternative to `fraction`) |

### Profiles

| Profile | Description |
|---------|-------------|
| (default) | Local execution, no container (use an activated conda env that provides the notebook kernel). |
| `local` | Local execution with Apptainer, sized for a laptop / WSL2 box (2 CPUs, 8 GB). Defaults `samplesheet` and `cell_ids_file` to the test assets. Also points `CLUSTER_SDATA_GPU` at the local RAPIDS container with WSL2 GPU passthrough settings. |
| `oscer` | SLURM executor on OSCER HPC, Apptainer container, scratch-based work directory. Memory scales 32→64→96 GB across retries. |

```bash
# Local profile (no --samplesheet needed)
nextflow run main.nf --step cluster_sdata_gpu -profile local

# OSCER HPC
nextflow run main.nf \
    --step create_sdata \
    --samplesheet assets/samplesheet.csv \
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
├── concat_sdata/
│   └── merged.zarr/
├── ROI1_A/
│   ├── downsample_xenium_region/
│   │   └── ROI1_A/
│   ├── create_sdata/
│   │   └── output/
│   │       └── ROI1_A.zarr/
│   ├── follicle_sdata/
│   │   └── output/
│   │       ├── aaaaimck-1.zarr/
│   │       └── aaameida-1.zarr/
│   ├── cluster_sdata/
│   │   └── clustered.zarr/
│   ├── cluster_sdata_gpu/
│   │   └── clustered.zarr/
│   ├── downsample_sdata/
│   │   └── downsampled.zarr/
│   └── plot_follicle/
│       ├── aaaaimck-1_plot_follicle.pptx
│       └── aaameida-1_plot_follicle.pptx
└── ROI1_B/
    └── ...
```

Each step publishes under `results/<sample>/<step>/`, except `concat_sdata`, which merges multiple samples and publishes once under `results/concat_sdata/`.
