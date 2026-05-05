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
nextflow run create.nf --samplesheet results/create_sdata/sample_sdata_samplesheet.csv --create follicle_sdata
nextflow run analyze.nf --samplesheet results/create_follicle_sdata/follicle_sdata_samplesheet.csv --analyze plot_follicle
```

## Adding A Notebook

1. Add a `.qmd` file with explicit Quarto `params` front matter.
2. Register it in [`../lib/NotebookRegistry.groovy`](../lib/NotebookRegistry.groovy) with the exact parameter list.
3. Wire it into `create.nf` or the analysis registry as appropriate.
