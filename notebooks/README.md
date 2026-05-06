# Notebook Workflows

This directory contains the Quarto notebooks used by the Nextflow pipelines.
Notebook registry metadata lives in [`../lib/NotebookRegistry.groovy`](../lib/NotebookRegistry.groovy).

## Notebooks

| Notebook | Purpose | Params | Main outputs |
|----------|---------|--------|--------------|
| `create_sdata.qmd` | Convert a raw Xenium output into a sample-level zarr. | `sample`, `path`, `n_jobs` | `output/<sample>.zarr`, `<sample>_create_sdata.html`, `<sample>_create_sdata.timing.tsv` |
| `create_follicle_sdata.qmd` | Subset one sample-level zarr into one zarr per cell ID. | `sample`, `path`, `cell_ids_file`, `radius` | `output/<cell_id>.zarr`, `<sample>_create_follicle_sdata.html`, `<sample>_create_follicle_sdata.timing.tsv` |
| `plot_follicle.qmd` | Render follicle zarrs into PowerPoint slides. | `sample`, `cell`, `path` | `<sample>_<cell>_plot_follicle.pptx`, `<sample>_<cell>_plot_follicle.timing.tsv` |

## Samplesheets

- Create workflow input: `sample,path`
- Analysis workflow input: `sample,path` for sample-level notebooks, `sample,cell,path` for notebooks that declare `cell`
- Follicle reference sheet: `Donor.ROI`, `cell_id`, optional `radius`

Examples:

```csv
sample,path
ROI1,/path/to/ROI1/xenium_output
```

```csv
sample,cell,path
ROI1,aaaaimck-1,results/ROI1/create_follicle_sdata/output/aaaaimck-1.zarr
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
| `plot_follicle` | `notebooks/plot_follicle.qmd` | `sample`, `cell`, `path` |

Notebook metadata is defined in [`../lib/NotebookRegistry.groovy`](../lib/NotebookRegistry.groovy), not under `params`.

## Adding A Notebook

1. Create a new `.qmd` file in `notebooks/` with a Jupyter `parameters` cell declaring the notebook inputs it expects.
2. Register the notebook in [`../lib/NotebookRegistry.groovy`](../lib/NotebookRegistry.groovy) with a unique ID, a path, and an explicit `params` list naming the keys to pass through.
3. If it is a create-stage producer, wire it into `create.nf`.
4. For analysis notebooks, include every required row-level key in the registered `params` list (for example, include `cell` for follicle-level runs).
5. Run analysis notebooks with `nextflow run analyze.nf --samplesheet <artifact_sheet.csv> --analyze <id1,id2|all>`.
