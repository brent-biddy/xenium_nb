# Notebook Workflows

This directory contains the Quarto notebooks used by the Nextflow pipelines.
Notebook registry metadata lives in [`../assets/notebook_registry.json`](../assets/notebook_registry.json).

The create-stage steps (`CREATE_SDATA`, `CREATE_FOLLICLE_SDATA`) run plain Python scripts
in `bin/` rather than notebooks — see `bin/create_sdata.py` and `bin/create_follicle_sdata.py`.

## Notebooks

| Notebook | Purpose | Params | Main outputs |
|----------|---------|--------|--------------|
| `analyze/plot_follicle.qmd` | Render follicle zarrs into PowerPoint slides. | `sample`, `cell`, `path` | `<cell>_plot_follicle.pptx`, `<cell>_plot_follicle.timing.tsv` |

## Samplesheets

- Create workflow input: `sample,path` — optionally `he_image,he_alignment` for H&E alignment
- Analysis workflow input: `sample,path` for sample-level notebooks, `sample,cell,path` for notebooks that declare `cell`
- Follicle reference sheet: `Donor.ROI`, `cell_id`, optional `radius`

Examples:

```csv
sample,path
ROI1_A,/path/to/ROI1_A/xenium_output
```

```csv
sample,cell,path
ROI1_A,aaaaimck-1,results/ROI1_A/follicle_sdata/output/aaaaimck-1.zarr
```

## Running

```bash
nextflow run create.nf --samplesheet assets/samplesheet.csv --create sdata
nextflow run create.nf --samplesheet results/sample_sdata_samplesheet.csv --create follicle_sdata
nextflow run analyze.nf --samplesheet results/follicle_sdata_samplesheet.csv --analyze plot_follicle
```

## Analysis Notebook IDs

The built-in analysis registry currently defines:

| ID | Notebook | Registered params |
|----|----------|-------------------|
| `plot_follicle` | `notebooks/analyze/plot_follicle.qmd` | `sample`, `cell`, `path` |

Notebook metadata is defined in [`../assets/notebook_registry.json`](../assets/notebook_registry.json).

## Adding an analysis notebook

1. Create `notebooks/analyze/<name>.qmd` with a `#| tags: [parameters]` Python cell declaring all inputs.
2. Add an entry to [`../assets/notebook_registry.json`](../assets/notebook_registry.json) with a unique ID, the relative path, and a `params` list matching the parameters cell exactly.
3. Wire a new process into `modules/analyze_notebooks.nf` and `analyze.nf`.
4. Run `python bin/check_notebook_registry.py` to verify.

## Adding a create-stage script

Create-stage steps run plain Python scripts, not notebooks.

1. Create `bin/<name>.py` with an `argparse` CLI matching the inputs the process needs.
2. Wire a new process into `modules/create_notebooks.nf` and `create.nf`.
3. No registry entry is needed for scripts.
