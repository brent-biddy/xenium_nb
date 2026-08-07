# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## What this pipeline does

`xenium_nb` is a Nextflow pipeline for Xenium spatial transcriptomics analysis. All steps run through a single entry point, `main.nf`, selected with `--step`:

- `downsample_xenium_region` — crops a raw Xenium output directory to a bounding box region
- `create_sdata` — converts raw Xenium output into a sample-level SpatialData zarr artifact. Annotates per-cell (`n_genes_by_transcripts`, `pct_transcripts_in_top_{5,10,20,50}_genes`, `control_counts`, `pct_control`, `transcripts_per_gene`, `nucleus_ratio`) and per-gene (`n_cells_by_transcripts`, `mean_transcripts`, `pct_dropout_by_transcripts`, `total_transcripts`) QC metrics but **filters nothing** — like `create_adata`, it produces the raw artifact and leaves thresholds to downstream analysis. Computing them here reads X once, at the only point it is already in memory, so a QC report can work from `obs`/`var` alone. Three details are Xenium-specific and easy to get wrong:
    - **`expr_type="transcripts"`, not the default `"counts"`.** The Xenium reader already puts a `total_counts` in `obs` from `cells.parquet`, and it is **not** a row sum of X — it is `transcript_counts` + the five control/codeword counters. The default `expr_type` would silently overwrite it with `X.sum(axis=1)`, destroying that meaning. Renaming the outputs sidesteps the collision. `cluster_sdata` and `cluster_sdata_gpu` used to clobber it this way and no longer recompute anything (see below); `cluster_sdata_gpu_ooc` still does
    - **Threshold on `transcript_counts`, never `total_counts`** — the latter includes negative controls, so a cell with 8 real transcripts and 4 deprecated codewords survives a `min_counts=10` cut. `transcript_counts` equals `X.sum(axis=1)` exactly. The redundant `obs["total_transcripts"]` that `calculate_qc_metrics` emits is dropped for this reason; per-gene `var["total_transcripts"]` is a different quantity and is kept
    - **No `qc_vars`.** The panel carries no `MT-` genes at all, so `pct_transcripts_mt` would read 0.0 for every cell — a metric that looks real and passes any threshold. Ribosomal/hemoglobin sets are only partially and arbitrarily represented on a targeted panel (and naive `^RP[SL]`/`^HB` patterns match non-members like `RPS6KA5`, `HBEGF`, `HBS1L`), so they are omitted rather than reported as if comparable to scRNA-seq. Negative-control % cannot use `qc_vars` either — controls are not features in X, only `obs` columns — so `control_counts` (the sum of the five control/codeword counters) and `pct_control` are written explicitly instead. `pct_control`'s denominator is `transcript_counts + control_counts`, not `obs["total_counts"]`: the two are identical on this object, but `total_counts` is the column `cluster_sdata_gpu_ooc` overwrites with a plain row sum, whereas `transcript_counts` always equals `X.sum(axis=1)`. It is NaN, not 0, where a cell has neither transcripts nor controls (2 cells in one test ROI) — 0 would read as a clean cell. `cluster_report` and `qc_report` both used to derive this number themselves; they now read it and keep their derivation only as a fallback for zarrs written before this column existed
    - **`percent_top` is `(5, 10, 20, 50)`, lower than `create_adata`'s `(10, 20, 50, 150)`** — a cell cannot detect more genes than the panel targets, and at a median 174 genes/cell, `top_150` pins 40% of cells at exactly 100%, inverting the metric so the emptiest cells score highest
- `create_adata` — converts a Cell Ranger `filtered_feature_bc_matrix` directory (scRNA-seq, not Xenium) into a sample-level `<sample>.h5ad`. Annotates per-cell/per-gene QC metrics (`total_counts`, `n_genes_by_counts`, `pct_counts_mt`, `percent_top`) but **filters nothing** — mirroring `create_sdata`, it produces the raw artifact and leaves thresholds to downstream analysis. Mito genes are auto-detected from the gene symbols with a case-insensitive `^[Mm][Tt]-` match (the same pattern used in the `oir-analysis` project), so human and mouse both work with no species flag. Output is a plain h5ad, so the `cluster_sdata*` steps (which read SpatialData zarrs) do **not** consume it — clustering scRNA-seq needs a new step
- `create_follicle_sdata` — subsets a sample zarr into per-cell follicle zarrs
- `cluster_sdata` / `cluster_sdata_gpu` — filter, normalize, PCA, neighbors, UMAP, Leiden clustering (CPU vs. RAPIDS/GPU). Neither **computes QC metrics** any more: `create_sdata` already annotated every per-cell and per-gene metric at the one point X was in memory, so both steps now only filter (`min_counts=10`, `min_cells=5`). Removing the redundant `calculate_qc_metrics` call is safe because `filter_cells`/`filter_genes` — scanpy's and rapids-singlecell's alike — derive their thresholds from X directly and never read `obs`; verified on a real ROI, the retained cells, retained genes, and all ten `leiden_res_*` columns are identical before and after. What it fixes is that the call ran with the default `expr_type="counts"` and so overwrote the Xenium-native `obs["total_counts"]`, and ran *before* `filter_genes`, leaving a value that described the full panel while the object shipped only the surviving genes (wrong for 802 of 21,719 cells on one ROI). Clustered zarrs therefore no longer carry `n_genes_by_counts`, `log1p_*`, or `pct_counts_in_top_*` — read `create_sdata`'s `n_genes_by_transcripts` and `pct_transcripts_in_top_*` instead. **Metrics on a clustered zarr are pre-filter by construction**, since they describe the raw artifact; the drift is negligible (mean transcripts/cell 277.74 vs 277.70) and it buys one definition of each metric everywhere
    - `cluster_sdata_gpu_ooc` is deliberately **not** part of this and still recomputes. Its filter reads `obs["total_counts"]`, which is correct *only because* its own `rsc.pp.calculate_qc_metrics` clobbers the column first — deleting the recompute without repointing that mask to `transcript_counts` turns accidental correctness into the exact bug described above. It also cannot simply inherit the gene-side metric: `anndata.concat` drops **all** `var` columns (obs survives whole), so the merged cohorts this step exists for arrive with an empty `var`, and cells-per-gene is a property of the merged object rather than of any one sample

All three `cluster_sdata*` steps **sweep Leiden resolutions** rather than committing to one. `--resolutions` takes a comma-separated list (`--resolutions 0.5,1.0,1.5`); the default sweep is `0.1` through `1.0` in steps of `0.1` (ten resolutions) and lives only in the Python scripts (`params.resolutions` is null in `nextflow.config`, and the modules omit the flag when it is unset, so the default is defined in one place). Each resolution gets its own obs column, `leiden_res_<r>` formatted to two decimals (`leiden_res_0.60`) — so two resolutions agreeing to two decimals (0.125, 0.126) collide on one column and the second silently wins; nothing checks for this. There is **no plain `leiden` column** — downstream code must select a resolution by name, and nothing records the sweep outside those columns. The sweep is cheap: neighbors/UMAP are computed once and only community detection re-runs per resolution (~0.5 s each on CPU for a ~100k-cell ROI, vs. ~19 s for the graph).
- `cluster_sdata_gpu_ooc` — same clustering pipeline as `cluster_sdata_gpu`, but streams the table's X matrix through Dask (rapids-singlecell out-of-core) so tables too large for VRAM (e.g. a merged cohort from `concat_sdata`) can still run on a single GPU. Optional `--n_top_genes` subsets to highly variable genes before PCA; it is **off by default** because it would cluster a different feature space than the other two steps, and a Xenium panel is already curated. Turn it on only if the materialized X does not fit. Note the chunked PCA differs from the in-memory one at ~1e-5, which is enough for Leiden to land ±1 cluster either way versus `cluster_sdata_gpu` — the embeddings themselves correlate at 1.000000
- `cluster_report` — renders `cluster_report.qmd` into one cohort-level pptx deck over every clustered zarr. Unlike the other notebook step it takes **no params**: every zarr is staged flat into the work dir and the notebook globs `*.zarr`, so staging is the whole input contract and there is no registry entry. Two details are easy to trip over: the zarrs must be staged under indexed names (`stageAs: 'sample*.zarr'`) because `cluster_sdata*` publishes every sample's output as `clustered.zarr` and a flat fan-in would be a name collision; and the qmd must pin `jupyter: python_spatial`, since Apptainer bind-mounts `$HOME` and a stray `~/.local` kernelspec otherwise shadows the container's, failing every import while still exiting 0 under `error: true`. Per sample it emits QC-on-UMAP, the resolution sweep on UMAP, the same sweep in tissue coordinates (both split across `CLUSTER_SLIDES` slides, so slide n of one holds the same resolutions as slide n of the other), and a full-slide tissue plot at `SPATIAL_RESOLUTION`. Marker size is derived from each axes' area and the cell count rather than fixed — one constant cannot serve both a five-panel sweep and a full-slide plot, nor a sparse strip and a dense section. Cluster colours come from tab10 below 11 clusters and tab20 above: tab20 is ten dark/light pairs, so the 2-cluster panels the low resolutions produce would be navy against pale blue
- `qc_report` — renders `qc_report.qmd` into one cohort-level pptx deck over every **raw `create_sdata` zarr**, to be read *before* clustering so the filtering thresholds are chosen from the data. Same staging-is-the-contract fan-in as `cluster_report` (indexed `stageAs`, glob `*.zarr`, no params, no registry entry, pinned `jupyter: python_spatial`), and the same `error: true` caveat — read the slides, not the exit status. Three things are specific to it:
    - **It never reads a value out of X.** Every metric it shows is already in `obs` from `create_sdata`, so the table is loaded with `ad.experimental.read_lazy` and X stays a Dask array whose graph is never computed. This is the payoff `create_sdata` was annotating QC metrics for. The one wrinkle: `read_lazy` returns `obs`/`var` as xarray `Dataset2D`, which anndata's own `strings_to_categoricals` chokes on (`unhashable type: 'DataArray'`) — and every `sc.pl`/`sq.pl` call sanitizes first, so both `to_memory()` lines in the loader are load-bearing. That is a [known anndata limitation](https://github.com/scverse/anndata/issues/981), not a quirk of this deck, and is expected to go away
    - **It must run on `create_sdata` output, not `cluster_sdata`'s.** Clustered zarrs have already had the `min_counts` cut applied, so the discarded tail is gone and the distributions look clean no matter how much was thrown away
    - **It derives nothing, and has no fallbacks.** Every obs column it reads is listed in `REQUIRED_OBS` and checked up front in the per-sample loop, so a zarr predating a column fails on the first sample with that column named rather than several slides later. A fallback derivation would be a second definition of the metric, letting a stale cohort render a deck whose numbers came from outside the pipeline — re-run `create_sdata` instead

    **The deck is currently an exploratory surface, not a settled report.** One list, `METRICS`, drives all three cell-level plot types, and every metric on it appears on every plot type it can — 40 slides per sample with the gene and transcript sections, ~60 MB for two ROIs. That is deliberate while the metrics are still being chosen: the panel is being read to decide which views are worth keeping, which cannot be decided from the ones that were left out. Expect it to be cut down, and treat the slide count and file size as temporary rather than as a target to optimise.

    Three sections per sample, **cells → genes → transcripts**, widest unit to narrowest, so each is a finer-grained account of the same sample rather than an unrelated appendix.

    `METRICS` entries are `(obs column, axis label, scale)`, and the **scale is what lets one list drive three plot types** — it picks the histogram's binning and axis and the tissue panel's colour norm, which are the same decision asked twice. Four kinds, and they are not interchangeable: `log` (counts and areas spanning decades; zeros cannot be drawn and are counted), `linear` (bounded ratios and percentages; zeros are **kept**, since `pct_control` is 0 for 87% of cells and the log panels' zero-dropping would empty the panel), `count` (small integers binned on half-offset integer edges — `nucleus_count` takes four distinct values, where 60 logspace bins are 56 empty bars), and `categorical` (no bins and no colourbar at all — a horizontal bar chart of category counts, and a discrete palette with a legend in tissue). Getting the norm wrong is not cosmetic: a log norm over `pct_control` (0–3%) or `nucleus_ratio` (0–1) spends the whole ramp below 1 and masks every cell at 0, so the panel reads as empty tissue.

    Per sample: per-cell QC histograms in 2x2 grids of `HIST_PANELS` (the metrics are small multiples of one measurement, so they stay gridded where the other two do not); transcripts against genes detected as one point per cell, once per metric — the joint view the margins cannot give, since a cell with many transcripts over few genes sits off the main arm while looking ordinary in either histogram; and the segmentation in tissue coordinates, once per metric. The scatter skips only the categorical metric, since `sc.pl.scatter` is the one call here whose colour handling this deck does not control. It **keeps the degenerate colours** — the two that are the scatter's own axes, and `transcripts_per_gene`, which on log axes is distance from the diagonal — because re-encoding position is not the same as showing nothing. There is no cohort section and the deck's only output is the pptx — it reports what the data looks like and stops, so choosing a threshold means reading the slides and setting the cut in `cluster_sdata*` yourself.

    **`codeword_category` is carried across sections as colour.** It lives on the transcripts, one value per molecule, but `feature_name` → `codeword_category` is strictly one-to-one (verified: no feature carries two), so a `groupby(...).first()` reindexed onto `var` is a **lookup, not a summary** — which is why it does not make the deck a second definition of anything. On the gene axis it resolves to predesigned vs custom panel probes, the question the histograms cannot answer: whether separately-designed custom probes land in the same abundance-detection cloud as the catalogue ones. Consequences worth knowing:
    - The molecule table is therefore loaded **before** the gene slides are built even though it is the last section to appear. If it fails the gene slides still render, falling back to a single colour — the geometry is theirs, only the colouring came from points.
    - `CATEGORY_INK` fixes one colour per category deck-wide, so a category keeps its colour between the gene scatter, the rank curve and the qv distributions, which are read against each other. Hue families carry the meaning: real genes blue/orange, every control and unassigned class red/purple/amber.
    - Reindexing leaves **`pd.NA`** for a gene with no transcript anywhere in the sample. `codeword_category` is pandas `string` dtype, so `pd.NA == "predesigned_gene"` evaluates to `pd.NA`, not `False`, and every comparison-built mask raises `boolean value of NA is ambiguous`. The lookup is normalised to object dtype with `None` **once, where it is built** — do not re-introduce the string dtype downstream. The unmapped state is real (absent even before cell assignment) and is named on the slides, not dropped.
    - The minority class is drawn **last and larger** on the scatter; at 100 in 5,101 plotted in index order it is buried and the slide shows nothing. The rank curve stays one colour and marks only the minority categories on it — marking 4,995 of 5,101 would redraw the line in a second colour and say nothing.

    The **gene section** reads `var`, which `create_sdata` annotated at the same point it annotated `obs`, so it derives nothing either. `GENE_METRICS` reuses the same `(column, label, scale)` shape and the same histogram function via its `level="var"` argument — the two populations are the same figure over a different frame, so writing them twice would only let them drift. `pct_dropout_by_transcripts` is the one linear entry: it is a bounded percentage running to 100 for the genes no cell detected, which a log axis drops entirely, and those are the genes worth finding. Beyond the histograms: a rank-abundance curve (rank linear, abundance log; zero-detection genes marked as a band at the foot rather than silently ending the curve early), an abundance-vs-detection scatter that separates "not expressed here" from "many transcripts in few cells", and a slide **naming** the undetected genes. That last one is the only list in the deck rather than a figure, because the useful output is which probes to look at and a name cannot be read off a distribution. On a test ROI it surfaces 30 of 5,101, and they are the expected absent-tissue markers for the section (Y-chromosome genes, eye, kidney, adrenal) — which is why the slide names them rather than calling them failures.

    The **transcript section** is the deck's **one deliberate exception to "derives nothing"**, and the notebook says so where the rule is stated. `qv`, codeword category, field of view and nucleus placement are properties of the molecule table, not summaries of X, so `create_sdata` never computed them and there is no pipeline definition for this deck to contradict — the rule exists to stop two definitions of one metric, and here there is only one. If any of these later earns a place in `create_sdata`'s `uns`, the deck should read it from there and the exception should shrink. Mechanically: `read_zarr` now takes `points` too (lazy Dask, free until asked), and `TX_COLUMNS` is materialized **once per sample** — never the whole frame, since coordinates and ids are most of its width. ~7M rows for a 21k-cell ROI takes about a second. This is the deck's only shared derived state, so it gets **its own `try`** rather than relying on the per-slide catch: if it fails, all four molecule slides are meaningless and are skipped together. A whole-slide sample is where it becomes the binding constraint, and the fix there is lazy Dask aggregation (`value_counts` and `groupby` both support it), not materialize-and-subsample. Four slides: the `qv` distribution with `QV_THRESHOLD` (20, Xenium's own default) drawn but **never applied** — the fraction below it is reported precisely because `transcript_counts` includes every transcript regardless of qv, and this is the only slide that says by how much; the **qv distribution per category**, each normalised to its own share so 162 negative-control probes are as readable as 5.4M genes — the control check the deck exists to make, since a negative control scoring like a real gene means the decoder cannot tell them apart and no per-cell metric would reveal it (on a test ROI: genes at median 40, `deprecated_codeword` 17, `negative_control_codeword` 10.5, `negative_control_probe` **0.0** with 96% below threshold), with the medians restated as a table under its **own explicit heading** because pandoc puts a table on a fresh slide regardless and it would otherwise arrive untitled; the codeword-category breakdown on a log axis (the categories span five decades, so on a linear axis every control bar is invisible, which is the opposite of what a control is for); per-FOV yield and median qv sharing a field order, the only view here that can see an *instrument* problem rather than a biological one, since a bad field's cells just look like poor cells in every per-cell distribution; and transcript placement (z, nucleus distance, nucleus overlap).

    Two details on the cell histograms are worth keeping straight. Cells fall off a **log** panel two ways and are counted separately: a **zero** cannot be drawn, and a value that was **never measured** is NaN — the reader leaves `nucleus_area` NaN, not 0, for a cell segmented without a nucleus (70 of 21,724 on one ROI, matching `nucleus_count == 0`), and `create_sdata`'s division carries that NaN into `nucleus_ratio`. Folding those together, or filtering non-finite values before counting, silently loses exactly the cells a QC deck exists to surface. Plotting an unmeasured cell at 0 would file it among genuinely cytoplasm-heavy ones.

    The tissue slide uses **`spatialdata_plot`'s `render_labels`**, not a scatter, so each cell is painted at its own segmented footprint. That is what retires the marker-size question — no constant travels between a 2,000 µm ROI and a whole slide, and deriving one from cell spacing still leaves a coverage factor that is pure taste. `render_labels` rather than `render_shapes` because the table annotates `cell_labels`, so it needs no schema change, while colouring `cell_boundaries` would need a second table keyed on `cell_id`; at this magnification a cell is ~2 pt across, so the polygons' extra fidelity is invisible anyway. Four traps here:
    - **`LogNorm` fails** (`spatialdata_plot` calls `norm.autoscale_None` on the 2D label array and matplotlib 3.11's generated `LogNorm` transforms it unconditionally), hence the `SafeLogNorm` shim — delete it when the libraries agree
    - spatialdata's `global` frame for a Xenium read is in **pixels**, so the axes are tick-less rather than mislabelled as µm, and the squidpy version's scale bar is gone with it
    - the slide is called **once per metric**, so each gets its own norm and colourbar; sharing a scale is not an option even in principle when genes top out near 2,400 against ~10,600 transcripts and 3% control
    - **never call `tight_layout` on this figure.** `spatialdata_plot` attaches its colourbar as a free-standing axes that `tight_layout` does not know about, so it expands the main axes from 77% to 92% of the figure width — straight under the colourbar — and the equal-aspect axes then centres itself in the oversized box, leaving a band of empty slide on its left. The panel is positioned by hand instead, fitted to the tissue's own aspect (read back from the drawn limits) and touching whichever edge of the slide binds first, with the colourbar re-seated beside it. These ROIs run both ways — one a 3:1 strip, the other a 1:2.2 block — so the binding edge differs per sample
- `concat_sdata` — merges multiple sample zarrs into one
- `downsample_sdata` — subsamples cells from a SpatialData zarr
- `plot_follicle` — renders the `plot_follicle.qmd` Quarto notebook per follicle zarr

## Commands

### Run a step
`--samplesheet` is always required; columns vary by step (see `main.nf`'s header comment for the full table). Some steps take extra flags.

```bash
nextflow run main.nf --step downsample_xenium_region --samplesheet assets/samplesheet.csv
nextflow run main.nf --step create_sdata --samplesheet assets/downsampled_region_samplesheet.csv
nextflow run main.nf --step create_adata --samplesheet assets/r21_adata_oscer_samplesheet.csv
nextflow run main.nf --step create_follicle_sdata --samplesheet my_sample_zarrs.csv --cell_ids_file assets/stage_quality_area_all_rois.csv
nextflow run main.nf --step cluster_sdata --samplesheet my_sample_zarrs.csv
nextflow run main.nf --step cluster_sdata --samplesheet my_sample_zarrs.csv --resolutions 0.5,1.0,1.5
nextflow run main.nf --step cluster_sdata_gpu --samplesheet my_sample_zarrs.csv
nextflow run main.nf --step cluster_sdata_gpu_ooc --samplesheet my_sample_zarrs.csv --chunk_size 20000
nextflow run main.nf --step cluster_report --samplesheet results/cluster_sdata_samplesheet.csv
nextflow run main.nf --step qc_report --samplesheet results/create_sdata_samplesheet.csv
nextflow run main.nf --step concat_sdata --samplesheet assets/concat_sdata_samplesheet.csv
nextflow run main.nf --step downsample_sdata --samplesheet my_sample_zarrs.csv --fraction 0.1
nextflow run main.nf --step plot_follicle --samplesheet assets/ci_analyze_samplesheet.csv
```

`my_sample_zarrs.csv` above is a stand-in for a `sample,path` CSV pointing at a prior step's output zarrs (e.g. `results/<sample>/create_sdata/<sample>.zarr`). Some producing steps now publish a ready-to-use handoff samplesheet next to their outputs (`<outdir>/create_sdata_samplesheet.csv`, `<outdir>/cluster_sdata_samplesheet.csv`) that you can point the next step's `--samplesheet` straight at; for steps without one yet, hand-build the CSV.

`downsample_xenium_region` requires the samplesheet to include `xmin,ymin,xmax,ymax` columns (µm coordinates) and an optional `region_name` column, which defaults to the sample ID if omitted. `downsample_sdata` requires `--fraction` or `--n_cells`.

### Profiles
Defined in `nextflow.config`:

| Profile | Executor | Container |
|---------|----------|-----------|
| (none)  | local, no container | requires activated conda env with Quarto + notebook deps |
| `local` | local, Apptainer | `babiddy755/python_spatial:1.2.0`, 8 CPUs, 16 GB |
| `oscer` | SLURM on OSCER HPC, Apptainer | same image, 16 CPUs, memory retries 48→96→144 GB (heavier for `CONCAT_SDATA`/`CLUSTER_SDATA`); GPU steps use the `sooner_gpu_test` partition with `--gres=gpu:1 --nv` |

**Run directories.** The `local` and `oscer` profiles set their own `workDir` and `outdir` so nothing lands in the repo (runs are typically launched from the repo root). Each run gets one self-contained directory, `<out_root>/<run_id>/{work,results}`, so a whole run is a single unit to size (`du -sh`) or prune (`rm -rf`). The shared Apptainer cache is a sibling of the run dirs (`<out_root>/apptainer_cache`), never nested under a `run_id`, so it survives across runs:

- `local` → `~/xenium_nb_out/<run_id>/{work,results}`
- `oscer` → `/scratch/$USER/xenium_nb_out/<run_id>/{work,results}`

Keeping `work` and `results` under the same root also keeps them on one filesystem, which the modules' hardlink publishing (`mode: 'link'`) relies on to avoid a second full copy of each zarr.

Both profiles set `cleanup = true`, so the work dir is deleted once the whole run **completes successfully** — leaving only the (hardlinked) results, i.e. each output stored exactly once. Nextflow scopes cleanup to success: a **failed** run keeps its work dir, so resume-after-failure (bump resources, skip the samples that already finished) still works. What cleanup forfeits is resume-*from*-success — reusing a finished run's cache on a later relaunch — which this single-process-per-step pipeline (no DAG to reuse) rarely needs.

Because `run_id` defaults to a fresh timestamp, `-resume` across separate launches only works if you pin the id with `--run_id <name>` (or recover the prior timestamp from the run dir name / `.nextflow.log` and pass it back). `-resume` must also be run from the same launch directory, since its cache lives in `.nextflow/` there. Note `cleanup = true` deletes that cache on success, so resume is available only after a failure.

The `local` profile defaults `samplesheet` and `cell_ids_file` to the test assets, and also points `cluster_sdata_gpu` / `cluster_sdata_gpu_ooc` at the local RAPIDS container with WSL2 GPU passthrough settings:

```bash
nextflow run main.nf --step cluster_sdata_gpu -profile local
nextflow run main.nf --step cluster_sdata_gpu_ooc -profile local
```

`cluster_sdata_gpu_ooc` additionally needs `dask` and `zarr` in the container — both are present in `babiddy755/python_spatial:1.2.0` as rapids-singlecell/spatialdata dependencies (verified).

### Stub run (CI-equivalent, no script/notebook execution)
```bash
nextflow run main.nf --step create_sdata -stub --samplesheet assets/samplesheet.csv
nextflow run main.nf --step create_adata -stub --samplesheet assets/ci_analyze_samplesheet.csv
nextflow run main.nf --step create_follicle_sdata -stub --samplesheet assets/ci_analyze_samplesheet.csv
nextflow run main.nf --step cluster_sdata -stub --samplesheet assets/ci_analyze_samplesheet.csv
nextflow run main.nf --step cluster_sdata_gpu -stub --samplesheet assets/ci_analyze_samplesheet.csv
nextflow run main.nf --step cluster_sdata_gpu_ooc -stub --samplesheet assets/ci_analyze_samplesheet.csv
nextflow run main.nf --step cluster_report -stub --samplesheet assets/ci_analyze_samplesheet.csv
nextflow run main.nf --step qc_report -stub --samplesheet assets/ci_analyze_samplesheet.csv
nextflow run main.nf --step concat_sdata -stub --samplesheet assets/ci_analyze_samplesheet.csv
nextflow run main.nf --step downsample_sdata -stub --samplesheet assets/ci_analyze_samplesheet.csv --fraction 0.1
nextflow run main.nf --step downsample_xenium_region -stub --samplesheet assets/samplesheet.csv
nextflow run main.nf --step plot_follicle -stub --samplesheet assets/ci_analyze_samplesheet.csv
```

### Validate notebook registry
```bash
python bin/check_notebook_registry.py
```

### Config parse check
```bash
nextflow config .
```

## Architecture

### Single entry point, one workflow per step
`main.nf` dispatches on `--step` to one of twelve named workflows, each of which reads a samplesheet, builds a channel of tuples, and pipes it into a single process. There is no chaining between steps inside Nextflow — to run steps in sequence, point the next step's `--samplesheet` at a CSV listing the prior step's published output paths (e.g. `results/<sample>/create_sdata/<sample>.zarr`). As a convenience (not a control-flow link), every artifact-producing step publishes a handoff samplesheet into `outdir` (`<step>_samplesheet.csv`) that you can feed directly to the next step. Each process emits a `samplesheet_row` output whose published path comes from a per-module helper that also drives that module's `publishDir` — so the convention is single-sourced in the module and `main.nf` just `.map { it.text }` + `collectFile`s the rows (the `.text` read makes `collectFile`'s `sort` deterministic). The row fragment is kept out of the publish dir via `publishDir`'s `saveAs`. Schemas: the single-per-sample producers (`create_sdata`, `create_adata`, `cluster_sdata`, `cluster_sdata_gpu`, `cluster_sdata_gpu_ooc`, `downsample_sdata`) and `concat_sdata` emit `sample,path`; `create_follicle_sdata` emits `sample,cell,path` (one row per per-cell zarr) for `plot_follicle`. `concat_sdata` and `create_follicle_sdata` also stage `.zarr` inputs into the work dir, so their row generation globs `*.zarr` and excludes the staged inputs.

### Create/cluster/downsample scripts (`bin/`)
Every step except `plot_follicle` runs a plain Python script with an `argparse` CLI (`bin/<step>.py`), invoked directly from its module's `script:` block — no params YAML involved.

### Notebook registry (`assets/notebook_registry.json`)
Maps analysis notebook IDs (currently just `plot_follicle`) to their `.qmd` path and the params they declare. This is the source of truth used by `modules/quarto_params.nf` at runtime and validated by `bin/check_notebook_registry.py` in CI. Every param listed in the registry must have a matching variable in the notebook's `#| tags: [parameters]` cell. The Python scripts under `bin/` are not registered here, and neither are `cluster_report.qmd` or `qc_report.qmd` — they declare no parameters, taking their input from what Nextflow stages beside them.

### Params YAML flow (`modules/quarto_params.nf`)
Used by the `plot_follicle` step only. `paramsFile()` writes `<outdir>/.quarto_params/<notebook>/params_<id>.yml` and returns the path for Nextflow staging. Writing to `outdir` (NFS) rather than `/tmp` is intentional — symlinks to head-node `/tmp` break on OSCER compute nodes.

### Process conventions
- Always use `script:` blocks, never `exec:` — processes must run through SLURM.
- Every process script sets `XDG_CACHE_HOME=$PWD/.cache` and `TMPDIR=$PWD/tmp` to avoid writing to a read-only compute-node `/tmp`.
- Keep named input variables; do not inline maps into process call arguments.
- Build command lines with optional arguments using a Groovy list + conditional append:

```groovy
def myArgs = ["--required_a ${val_a}", "--required_b ${val_b}"]
if (optional_c) myArgs << "--optional_c ${optional_c}"
"""
my_script.py ${myArgs.join(' ')}
"""
```

## Adding an analysis notebook

1. Create `notebooks/analyze/<name>.qmd` with a `#| tags: [parameters]` Python cell declaring all inputs.
2. Add an entry to `assets/notebook_registry.json` with the notebook ID, relative path, and `params` list matching the parameters cell exactly.
3. Wire a new process into `modules/<name>.nf` and add a matching `--step` branch in `main.nf`.
4. Run `python bin/check_notebook_registry.py` to verify.

## Adding a create/cluster/downsample-stage script

These steps use plain Python scripts, not notebooks.

1. Create `bin/<name>.py` with an `argparse` CLI (`parse_args()` function) declaring all inputs.
2. Wire a new process into `modules/<name>.nf` and add a matching `--step` branch in `main.nf`, passing args directly.
3. No registry entry is needed.

## CI

Two GitHub Actions run on PRs to `main`:
- **Validate notebook registry** — runs `python bin/check_notebook_registry.py`
- **Stub run** — runs every `main.nf --step` with `-stub` to verify workflow wiring without executing scripts or notebooks

## Code style (`.nf` files)

- 4-space indentation
- Process names in `UPPER_SNAKE_CASE`; params, variables, CSV headers in `snake_case`
- Add file-level header comments, docstrings on helper functions, section markers, and WHY comments for non-obvious decisions
- Annotate channel shape at every `.set {}` call and after non-obvious transformations so the tuple structure is always visible without tracing back through the chain:

```groovy
.set { createSdataInputs } // tuple(sample, staged_path, he_image, he_alignment)
// createSdataRun.artifacts: tuple(sample, zarr)
```
