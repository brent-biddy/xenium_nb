# xenium_nb

A Nextflow pipeline for Xenium spatial transcriptomics data with two separate entry points:

- `create.nf` builds reusable data artifacts
- `analyze.nf` runs analysis notebooks against artifact samplesheets

Notebook-specific inputs, outputs, and samplesheet contracts are documented in [notebooks/README.md](notebooks/README.md).

---

## Requirements

- [Nextflow](https://www.nextflow.io/) ≥ 25.10.0
- For default non-container local runs: [Quarto](https://quarto.org/) ≥ 1.4 and the required Python notebook packages, including `papermill`
- For containerized local runs: Apptainer

---

## Repository structure

```
xenium_nb/
├── create.nf                  # Create workflow: raw Xenium -> sample and follicle artifacts
├── analyze.nf                 # Analysis workflow: artifact samplesheet -> notebook reports
├── lib/                       # (reserved for Nextflow helper classes)
├── nextflow.config            # Parameters and profiles
├── modules/
│   ├── create_notebooks.nf        # Create-stage notebook processes (sdata, follicle_sdata)
│   ├── analyze_notebooks.nf       # Analysis notebook processes
│   ├── run_notebook.nf            # Shared notebook runner
│   └── write_quarto_params.nf     # Renders params.yml for notebooks
├── notebooks/
│   ├── README.md
│   ├── create_sdata.qmd
│   ├── create_follicle_sdata.qmd
│   └── plot_follicle.qmd
├── bin/
│   └── timer.py                   # Timing utilities for notebooks
└── assets/
    ├── samplesheet.csv                    # Sample-level samplesheet
    ├── stage_quality_area_all_rois.csv    # Cell ID reference file
    └── notebook_registry.json             # Notebook metadata (paths and declared params)
```

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
- `results/sample_sdata_samplesheet.csv`
- `results/follicle_sdata_samplesheet.csv`

### Create sdata only

```bash
nextflow run create.nf \
    --samplesheet assets/samplesheet.csv \
    --create sdata
```

This runs only `create_sdata.qmd` and writes `results/sample_sdata_samplesheet.csv`.
That sheet can be used as input to `--create follicle_sdata`.

### Create follicle sdata only

Run this after `create sdata` has produced `results/sample_sdata_samplesheet.csv`:

```bash
nextflow run create.nf \
    --samplesheet results/sample_sdata_samplesheet.csv \
    --create follicle_sdata
```

This runs only `create_follicle_sdata.qmd` and writes `results/follicle_sdata_samplesheet.csv`.

### Create with the small test follicle file

```bash
nextflow run create.nf \
    --samplesheet assets/samplesheet.csv \
    --cell_ids_file assets/small_stage_quality_area_all_rois.csv
```

### Analyze follicle artifacts

```bash
nextflow run analyze.nf \
    --samplesheet results/follicle_sdata_samplesheet.csv \
    --analyze plot_follicle
```

### Analyze sample artifacts

When you add analysis notebooks that only require sample-level artifacts, point `analyze.nf` at the sample artifact sheet:

```bash
nextflow run analyze.nf \
    --samplesheet results/sample_sdata_samplesheet.csv \
    --analyze your_sample_notebook_id
```

---

## Configuration

Key parameters (set in `nextflow.config` or passed via `--param value`):

| Parameter | Default | Description |
|-----------|---------|-------------|
| `samplesheet` | `null` | Path to samplesheet CSV |
| `outdir` | `results` | Output directory |
| `cell_ids_file` | `${projectDir}/assets/stage_quality_area_all_rois.csv` | Cell ID reference file path |
| `radius` | `250` | Default bounding box radius (µm) |
| `create` | `all` | Create workflow mode: `sdata`, `follicle_sdata`, or `all` |
| `analyze` | `all` | Analysis notebook selector: `all` or a comma-separated list of notebook IDs from `assets/notebook_registry.json` |

Analysis notebook IDs and how to add new notebooks are documented in [notebooks/README.md](notebooks/README.md).

### Profiles

| Profile | Description |
|---------|-------------|
| (default) | Local execution, no container (use an activated conda env that provides the notebook kernel). |
| `test` | Local execution with Apptainer, sized for a laptop / WSL2 box (8 GB default). Build a local `.sif` with `container/build_apptainer.sh` and update `container` in the `test` profile. |
| `oscer` | SLURM executor on OSCER HPC, Apptainer container, scratch-based work directory and image cache. Memory scales 32→64→96 GB across retries. |

Activate with `-profile oscer`:

```bash
nextflow run create.nf \
    --samplesheet assets/samplesheet.csv \
    -profile oscer
```

For local Apptainer validation:

```bash
conda run -n squidpy nextflow run analyze.nf \
    -profile test \
    --samplesheet /path/to/follicle_sdata_samplesheet.csv \
    --outdir /path/to/outdir \
    -process.memory '16 GB'
```

To publish the same runtime for OSCER:

```bash
./container/build_docker.sh
docker tag xenium_tools_squidpy:local babiddy755/xenium_nb:<tag>
docker push babiddy755/xenium_nb:<tag>
# Update container in the oscer profile in nextflow.config, then:
nextflow run create.nf --samplesheet assets/samplesheet.csv -profile oscer
```

---

## Output structure

```
results/
├── pipeline_info/
│   ├── timeline.html
│   ├── report.html
├── sample_sdata_samplesheet.csv
├── follicle_sdata_samplesheet.csv
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

If `--create sdata` is used, `create_follicle_sdata/` outputs and `results/follicle_sdata_samplesheet.csv` are not created.

The sample-stage sheet is the handoff input for `--create follicle_sdata`.

