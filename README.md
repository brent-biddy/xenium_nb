# xenium_nb

A Nextflow pipeline for Xenium spatial transcriptomics data with two separate entry points:

- `create.nf` builds reusable data artifacts
- `analyze.nf` runs analysis notebooks against artifact samplesheets

The create workflow and sample-scoped analysis notebooks use a two-column samplesheet contract:

```csv
sample,path
```

Follicle-scoped analysis notebooks use an explicit three-column contract:

```csv
sample,cell,path
```

Notebook parameters are staged into each task work directory as `params.json` and loaded explicitly by the notebook code.

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
├── conf/
│   └── base.config            # Resource defaults
├── modules/
│   ├── create_spatialdata.nf  # Sample-level artifact producer
│   ├── subset_follicle.nf     # Follicle-level artifact producer
│   ├── run_notebook.nf        # Generic analysis notebook runner
│   └── write_samplesheet.nf   # Writes artifact samplesheets
├── notebooks/
│   ├── create_sdata.qmd
│   ├── create_follicle_sdata.qmd
│   └── plot_follicle.qmd
├── bin/
│   ├── timer.py               # Timing utilities for notebooks
│   ├── make_follicle_samplesheet.py  # Legacy helper for manual/export workflows
│   └── downsample_xenium.py   # Regenerates a smaller Xenium output for workflow testing
└── assets/
    ├── samplesheet.csv        # Sample-level samplesheet
    └── stage_quality_area_all_rois.csv  # Cell ID reference file
```

---

## Samplesheets

### Create workflow input

Used by `create.nf`. `path` points to a raw Xenium output directory.

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
ROI1,results/ROI1/create_sdata/output/ROI1.zarr
ROI2,results/ROI2/create_sdata/output/ROI2.zarr
```

Follicle artifact sheet:

```csv
sample,cell,path
ROI1,aaaaimck-1,results/ROI1/create_follicle_sdata/output/aaaaimck-1.zarr
ROI1,aaaalpdj-1,results/ROI1/create_follicle_sdata/output/aaaalpdj-1.zarr
```

### Cell ID reference file

`assets/stage_quality_area_all_rois.csv` is used by the create workflow to decide which follicles to subset. It must contain `Donor.ROI` and `cell_id`. An optional `radius` column sets a per-cell bounding box radius (µm); missing values fall back to `params.radius`.

`create.nf` can select this file by key. The built-in choices are:

- `--cell_ids_file full`
- `--cell_ids_file small`

### Downsampling Xenium test data

`bin/downsample_xenium.py` regenerates a smaller Xenium output directory by spatially subsampling cells and rebuilding the associated Xenium sidecar files, zarr archives, and reduced-resolution morphology OME-TIFF pyramids.

Example:

```bash
conda run -n squidpy python bin/downsample_xenium.py /path/to/xenium_output --proportion 0.05
```

The output directory is written alongside the input as `<input_dir>_downsampled_<pct>pct`.

Current follow-up items for `bin/downsample_xenium.py`:

- The script still emits many duplicate-entry warnings while writing zipped zarr outputs.
- The script currently relies on zarr compatibility helpers added for the local environment; this should be revisited and simplified later.
- It would be useful to add a deterministic option to force retention of specific cell IDs for workflow tests.

---

## Usage

### Create artifacts

```bash
nextflow run create.nf \
    --samplesheet assets/samplesheet.csv
```

By default this runs both producer notebooks:

- `create_sdata.qmd`
- `create_follicle_sdata.qmd`

and writes:

- sample zarrs under `results/<sample>/create_sdata/output/`
- follicle zarrs under `results/<sample>/create_follicle_sdata/output/`
- `results/pipeline_info/sample_analysis_inputs.csv`
- `results/pipeline_info/follicle_analysis_inputs.csv`

### Create sample artifacts only

```bash
nextflow run create.nf \
    --samplesheet assets/samplesheet.csv \
    --run_subset_follicle false
```

This runs only `create_sdata.qmd` and writes `results/pipeline_info/sample_analysis_inputs.csv`.

### Create with the small test follicle file

```bash
nextflow run create.nf \
    --samplesheet assets/samplesheet.csv \
    --cell_ids_file small
```

### Analyze follicle artifacts

```bash
nextflow run analyze.nf \
    --samplesheet results/pipeline_info/follicle_analysis_inputs.csv \
    --notebooks plot_follicle
```

### Analyze sample artifacts

When you add sample-scoped analysis notebooks to the registry, point `analyze.nf` at the sample artifact sheet:

```bash
nextflow run analyze.nf \
    --samplesheet results/pipeline_info/sample_analysis_inputs.csv \
    --notebooks your_sample_notebook_id
```

---

## Configuration

Key parameters (set in `nextflow.config` or passed via `--param value`):

| Parameter | Default | Description |
|-----------|---------|-------------|
| `samplesheet` | `null` | Path to samplesheet CSV |
| `outdir` | `results` | Output directory |
| `cell_ids_file` | `full` | Cell ID reference file key or direct path |
| `container_image` | `babiddy755/xenium_nb:latest` | Container reference pulled by the `local` and `oscer` profiles; may be a registry tag or a local `.sif` path |
| `cell_ids_registry` | built-in map | Named cell ID files available to `create.nf` |
| `radius` | `250` | Default bounding box radius (µm) |
| `run_subset_follicle` | `true` | Whether `create.nf` should run `create_follicle_sdata.qmd` after building sample-level zarrs |
| `producer_registry` | built-in map | The two producer notebooks used by `create.nf` |
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
| (default) | Local execution, no container (use an activated conda env that provides the notebook kernel). |
| `local` | Local execution with Apptainer, sized for a laptop / WSL2 box (16 GB per process). Override `--container_image` with a local `.sif` path after building `container/Apptainer.def`. |
| `oscer` | SLURM executor on OSCER HPC, Apptainer container, scratch-based work directory and image cache. Memory scales 32→64→96 GB across retries. |

Activate with `-profile oscer`:

```bash
nextflow run create.nf \
    --samplesheet assets/samplesheet.csv \
    -profile oscer
```

For local Apptainer validation:

```bash
HOME=/tmp/xenium_home conda run -n squidpy nextflow run analyze.nf \
    -profile local \
    --samplesheet /tmp/xenium_nb_test/follicle_analysis_inputs.csv \
    --notebooks plot_follicle \
    --outdir /home/babiddy/xenium_nb_results_fresh \
    -process.memory '16 GB'
```

To publish the same runtime for OSCER:

```bash
./container/build_docker.sh
docker tag xenium_tools_squidpy:local babiddy755/xenium_nb:<tag>
docker push babiddy755/xenium_nb:<tag>
nextflow run create.nf --samplesheet assets/samplesheet.csv -profile oscer --container_image babiddy755/xenium_nb:<tag>
```

---

## Output structure

```
results/
├── pipeline_info/
│   ├── timeline.html
│   ├── report.html
│   ├── sample_analysis_inputs.csv
│   └── follicle_analysis_inputs.csv
├── ROI1/
│   ├── create_sdata/
│   │   ├── ROI1_create_sdata.html
│   │   └── output/
│   │       └── ROI1.zarr/
│   ├── create_follicle_sdata/
│   │   ├── ROI1_create_follicle_sdata.html
│   │   └── output/
│   │       ├── aaaaimck-1.zarr/
│   │       └── aaaalpdj-1.zarr/
│   └── plot_follicle/
│       ├── ROI1_aaaaimck-1_plot_follicle.pptx
│       ├── ROI1_aaaaimck-1_plot_follicle.timing.tsv
│       ├── ROI1_aaaalpdj-1_plot_follicle.pptx
│       └── ROI1_aaaalpdj-1_plot_follicle.timing.tsv
└── ROI2/
    └── ...
```

Analysis outputs also publish under the parent sample directory, so follicle reports for cells like `aaaaimck-1` and `aaaalpdj-1` both land under `results/ROI1/plot_follicle/`.
For `plot_follicle`, the primary rendered artifact is a `.pptx` deck, with a companion timing TSV.

If `--run_subset_follicle false` is used, `create_follicle_sdata/` outputs and `follicle_analysis_inputs.csv` are not created.

---

## Adding notebooks

1. Create a new `.qmd` file in `notebooks/` with a YAML params block declaring at least `sample` and `path`, plus any scope-specific fields such as `cell` for follicle notebooks.
2. If it is a create-stage producer, register it in `params.producer_registry` and wire it into `create.nf`.
3. If it is an analysis notebook, register it in `params.analysis_notebook_registry` with a unique ID and a scope of `sample` or `follicle`.
4. Run analysis notebooks with `nextflow run analyze.nf --samplesheet <artifact_sheet.csv> --notebooks <id1,id2>`.
5. Any params not declared in the notebook's front matter are automatically filtered out before rendering.
