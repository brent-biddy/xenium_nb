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
├── lib/
│   ├── NotebookRegistry.groovy # Pipeline-internal notebook metadata catalog
│   └── QuartoParams.groovy    # Shared Quarto parameter filtering/merge helper
├── nextflow.config            # Parameters and profiles
├── modules/
│   ├── create_spatialdata.nf      # Sample-level artifact producer
│   ├── subset_follicle.nf         # Follicle-level artifact producer
│   ├── run_notebook.nf            # Shared notebook runner
│   ├── write_quarto_params.nf     # Renders params.yml for notebooks
│   └── write_samplesheet.nf       # Writes artifact samplesheets
├── notebooks/
│   ├── README.md
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

## Notebook Docs

Detailed notebook contracts and examples live in [notebooks/README.md](notebooks/README.md).

## Downsampling Xenium test data

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
- `results/create_sdata/sample_sdata_samplesheet.csv`
- `results/create_follicle_sdata/follicle_sdata_samplesheet.csv`

### Create sdata only

```bash
nextflow run create.nf \
    --samplesheet assets/samplesheet.csv \
    --create sdata
```

This runs only `create_sdata.qmd` and writes `results/create_sdata/sample_sdata_samplesheet.csv`.
That sheet can be used as input to `--create follicle_sdata`.

### Create follicle sdata only

Run this after `create sdata` has produced `results/create_sdata/sample_sdata_samplesheet.csv`:

```bash
nextflow run create.nf \
    --samplesheet results/create_sdata/sample_sdata_samplesheet.csv \
    --create follicle_sdata
```

This runs only `create_follicle_sdata.qmd` and writes `results/create_follicle_sdata/follicle_sdata_samplesheet.csv`.

### Create with the small test follicle file

```bash
nextflow run create.nf \
    --samplesheet assets/samplesheet.csv \
    --cell_ids_file assets/small_stage_quality_area_all_rois.csv
```

### Analyze follicle artifacts

```bash
nextflow run analyze.nf \
    --samplesheet results/create_follicle_sdata/follicle_sdata_samplesheet.csv \
    --analyze plot_follicle
```

### Analyze sample artifacts

When you add analysis notebooks that only require sample-level artifacts, point `analyze.nf` at the sample artifact sheet:

```bash
nextflow run analyze.nf \
    --samplesheet results/create_sdata/sample_sdata_samplesheet.csv \
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
| `container_image` | `babiddy755/xenium_nb:20260505-66addc7` | Container reference pulled by the `test` and `oscer` profiles; may be a registry tag or a local `.sif` path |
| `radius` | `250` | Default bounding box radius (µm) |
| `create` | `all` | Create workflow mode: `sdata`, `follicle_sdata`, or `all` |
| `analyze` | `all` | Analysis notebook selector: `all` or a comma-separated list of notebook IDs from `lib/NotebookRegistry.groovy` |

### Analysis notebook IDs

The built-in analysis registry currently defines:

| ID | Notebook | Registered params |
|----|----------|-------------------|
| `plot_follicle` | `notebooks/plot_follicle.qmd` | `sample`, `cell`, `path` |

Notebook metadata is defined in [`lib/NotebookRegistry.groovy`](lib/NotebookRegistry.groovy), not under `params`.

### Profiles

| Profile | Description |
|---------|-------------|
| (default) | Local execution, no container (use an activated conda env that provides the notebook kernel). |
| `test` | Local execution with Apptainer, sized for a laptop / WSL2 box (8 GB default). Override `--container_image` with a local `.sif` path after building `container/Apptainer.def`. |
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
    -profile test \
    --samplesheet /tmp/xenium_nb_test/follicle_sdata_samplesheet.csv \
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
├── create_sdata/
│   └── sample_sdata_samplesheet.csv
├── create_follicle_sdata/
│   └── follicle_sdata_samplesheet.csv
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

If `--create sdata` is used, `create_follicle_sdata/` outputs and `results/create_follicle_sdata/follicle_sdata_samplesheet.csv` are not created.

The sample-stage sheet is the handoff input for `--create follicle_sdata`.

---

## Adding notebooks

1. Create a new `.qmd` file in `notebooks/` with a Jupyter `parameters` cell declaring the notebook inputs it expects.
2. Register the notebook in `lib/NotebookRegistry.groovy` with a unique ID, a path, and an explicit `params` list naming the keys to pass through.
3. If it is a create-stage producer, wire it into `create.nf`.
4. For analysis notebooks, include every required row-level key in the registered `params` list (for example, include `cell` for follicle-level runs).
5. Run analysis notebooks with `nextflow run analyze.nf --samplesheet <artifact_sheet.csv> --analyze <id1,id2|all>`.
