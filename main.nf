#!/usr/bin/env nextflow

// Single entry point for all pipeline steps.
// Select a step with: nextflow run main.nf --step <name> --samplesheet <path>
//
// Steps:
//   downsample_xenium_region  samplesheet: sample, path, xmin, ymin, xmax, ymax[, region_name, he_image, he_alignment]
//   create_sdata              samplesheet: sample, path[, he_image, he_alignment]
//   create_adata              samplesheet: sample, path
//   create_follicle_sdata     samplesheet: sample, path  (+ --cell_ids_file)
//   cluster_sdata             samplesheet: sample, path[, min_counts, min_cells, max_counts_quantile]  (+ --resolutions)
//   cluster_sdata_gpu         samplesheet: sample, path[, min_counts, min_cells, max_counts_quantile]  (+ --resolutions)
//   cluster_sdata_gpu_ooc     samplesheet: sample, path[, min_counts, min_cells, max_counts_quantile]  (+ --chunk_size, --n_top_genes, --resolutions)
//   create_centroids          samplesheet: sample, path  (clustered zarrs; + --group_by)
//   cluster_report            samplesheet: sample, path  (clustered zarrs; one deck for the cohort)
//   qc_report                 samplesheet: sample, path  (raw create_sdata zarrs; one deck for the cohort)
//   sample_summary            samplesheet: sample, path, centroid_path  (create_centroids sheet; + --chosen_resolutions)
//   concat_sdata              samplesheet: path
//   downsample_sdata          samplesheet: sample, path  (+ --fraction or --n_cells)
//   plot_follicle             samplesheet: sample, cell, path

include { DOWNSAMPLE_XENIUM_REGION } from './modules/downsample_xenium_region'
include { CREATE_SDATA }             from './modules/create_sdata'
include { CREATE_ADATA }             from './modules/create_adata'
include { CREATE_FOLLICLE_SDATA }    from './modules/create_follicle_sdata'
include { CLUSTER_SDATA }            from './modules/cluster_sdata'
include { CLUSTER_SDATA_GPU }        from './modules/cluster_sdata_gpu'
include { CLUSTER_SDATA_GPU_OOC }    from './modules/cluster_sdata_gpu_ooc'
include { CREATE_CENTROIDS }         from './modules/create_centroids'
include { CLUSTER_REPORT }           from './modules/cluster_report'
include { QC_REPORT }                from './modules/qc_report'
include { SAMPLE_SUMMARY }           from './modules/sample_summary'
include { CONCAT_SDATA }             from './modules/concat_sdata'
include { DOWNSAMPLE_SDATA }         from './modules/downsample_sdata'
include { PLOT_FOLLICLE }            from './modules/plot_follicle'
include { paramsFile }               from './modules/quarto_params'

// ── Helpers ───────────────────────────────────────────────────────────────────

// Resolve one filtering threshold for one sample, most specific source winning:
// the samplesheet column, else the cohort-wide param, else '' meaning "pass no flag
// and let the script's own default stand".
//
// An empty or whitespace-only cell counts as absent, so a column filled in for only
// some samples falls back per row instead of passing an empty flag. "0" does NOT count
// as absent — zero is a meaningful threshold (keep every cell), which is also why the
// param arm tests != null rather than using ?:.
def resolveThreshold(rowValue, paramValue) {
    if (rowValue?.toString()?.trim()) return rowValue.toString().trim()
    return paramValue != null ? paramValue : ''
}

// ── Entry workflow ────────────────────────────────────────────────────────────

workflow {
    if (!params.step) error "Please provide --step <name>. Valid steps: downsample_xenium_region, create_sdata, create_adata, create_follicle_sdata, cluster_sdata, cluster_sdata_gpu, cluster_sdata_gpu_ooc, create_centroids, cluster_report, qc_report, sample_summary, concat_sdata, downsample_sdata, plot_follicle"

    if      (params.step == 'downsample_xenium_region')  downsample_xenium_region()
    else if (params.step == 'create_sdata')              create_sdata()
    else if (params.step == 'create_adata')              create_adata()
    else if (params.step == 'create_follicle_sdata')     create_follicle_sdata()
    else if (params.step == 'cluster_sdata')             cluster_sdata()
    else if (params.step == 'cluster_sdata_gpu')         cluster_sdata_gpu()
    else if (params.step == 'cluster_sdata_gpu_ooc')     cluster_sdata_gpu_ooc()
    else if (params.step == 'create_centroids')          create_centroids()
    else if (params.step == 'cluster_report')            cluster_report()
    else if (params.step == 'qc_report')                 qc_report()
    else if (params.step == 'sample_summary')            sample_summary()
    else if (params.step == 'concat_sdata')              concat_sdata()
    else if (params.step == 'downsample_sdata')          downsample_sdata()
    else if (params.step == 'plot_follicle')             plot_follicle()
    else error "Unknown --step '${params.step}'. Valid steps: downsample_xenium_region, create_sdata, create_adata, create_follicle_sdata, cluster_sdata, cluster_sdata_gpu, cluster_sdata_gpu_ooc, create_centroids, cluster_report, qc_report, sample_summary, concat_sdata, downsample_sdata, plot_follicle"
}

// ── downsample_xenium_region ──────────────────────────────────────────────────

workflow downsample_xenium_region {
    if (!params.samplesheet) error "Please provide --samplesheet"

    channel
        .fromPath(params.samplesheet)
        .splitCsv(header: true)        // Map(sample, path, xmin, ymin, xmax, ymax[, region_name, he_image, he_alignment])
        .map { row ->
            if (!row.sample) error "Samplesheet row missing 'sample': ${row}"
            if (!row.path)   error "Samplesheet row missing 'path': ${row}"
            // `file()`, not a bare path string: these are staged `path` inputs, so the
            // container can read them wherever they live. `[]` is the optional-input
            // sentinel — it stages nothing and renders falsy in the process script.
            def heImage    = row.he_image     ? file(row.he_image     as String, checkIfExists: true) : []
            def heAlign    = row.he_alignment ? file(row.he_alignment as String, checkIfExists: true) : []
            def regionName = row.region_name ?: row.sample
            tuple(row.sample, file(row.path), row.xmin, row.ymin, row.xmax, row.ymax, regionName, heImage, heAlign)
        }                              // tuple(sample, path, xmin, ymin, xmax, ymax, region_name, he_image, he_alignment)
        | DOWNSAMPLE_XENIUM_REGION
}

// ── create_sdata ──────────────────────────────────────────────────────────────

workflow create_sdata {
    if (!params.samplesheet) error "Please provide --samplesheet"

    channel
        .fromPath(params.samplesheet)
        .splitCsv(header: true)        // Map(sample, path[, he_image, he_alignment])
        .map { row ->
            if (!row.sample) error "Samplesheet row missing 'sample': ${row}"
            if (!row.path)   error "Samplesheet row missing 'path': ${row}"
            // `file()`, not a bare path string: these are staged `path` inputs, so the
            // container can read them wherever they live. `[]` is the optional-input
            // sentinel — it stages nothing and renders falsy in the process script.
            def heImage = row.he_image     ? file(row.he_image     as String, checkIfExists: true) : []
            def heAlign = row.he_alignment ? file(row.he_alignment as String, checkIfExists: true) : []
            tuple(row.sample, file(row.path), heImage, heAlign)
        }                              // tuple(sample, path, he_image, he_alignment)
        | CREATE_SDATA

    // Aggregate the per-sample rows the process emits into a ready-to-use handoff
    // samplesheet, so a downstream step (cluster_sdata, downsample_sdata,
    // concat_sdata, create_follicle_sdata) can be pointed straight at it instead of
    // hand-building a sample,path CSV. The published path lives in the module (its
    // publishDir and the emitted row share one helper), so main.nf stays agnostic.
    CREATE_SDATA.out.samplesheet_row
        .map { it.text }             // read row content so collectFile's sort is deterministic
        .collectFile(name: 'create_sdata_samplesheet.csv', storeDir: params.outdir,
                     seed: 'sample,path', newLine: true, sort: true)
}

// ── create_adata ──────────────────────────────────────────────────────────────

workflow create_adata {
    if (!params.samplesheet) error "Please provide --samplesheet"

    channel
        .fromPath(params.samplesheet)
        .splitCsv(header: true)      // Map(sample, path)
        .map { row ->
            if (!row.sample) error "Samplesheet row missing 'sample': ${row}"
            if (!row.path)   error "Samplesheet row missing 'path': ${row}"
            tuple(row.sample, file(row.path))
        }                            // tuple(sample, path)
        | CREATE_ADATA

    // Handoff samplesheet of the per-sample h5ads (see create_sdata for the general
    // rationale). No clustering step reads h5ad yet, so this currently serves ad hoc
    // downstream use rather than another --step.
    CREATE_ADATA.out.samplesheet_row
        .map { it.text }             // read row content so collectFile's sort is deterministic
        .collectFile(name: 'create_adata_samplesheet.csv', storeDir: params.outdir,
                     seed: 'sample,path', newLine: true, sort: true)
}

// ── create_follicle_sdata ─────────────────────────────────────────────────────

workflow create_follicle_sdata {
    if (!params.samplesheet)   error "Please provide --samplesheet"
    if (!params.cell_ids_file) error "Please provide --cell_ids_file"

    def cellIdsFile = file(params.cell_ids_file)

    def inputs = channel
        .fromPath(params.samplesheet)
        .splitCsv(header: true)      // Map(sample, path)
        .map { row ->
            if (!row.sample) error "Samplesheet row missing 'sample': ${row}"
            if (!row.path)   error "Samplesheet row missing 'path': ${row}"
            tuple(row.sample, file(row.path))
        }                            // tuple(sample, path)

    CREATE_FOLLICLE_SDATA(inputs, cellIdsFile, params.radius)

    // Handoff samplesheet of the per-cell follicle zarrs for plot_follicle. Uses
    // the sample,cell,path schema (see create_sdata for the general rationale).
    CREATE_FOLLICLE_SDATA.out.samplesheet_row
        .map { it.text }             // read row content so collectFile's sort is deterministic
        .collectFile(name: 'create_follicle_sdata_samplesheet.csv', storeDir: params.outdir,
                     seed: 'sample,cell,path', newLine: true, sort: true)
}

// ── cluster_sdata ─────────────────────────────────────────────────────────────

workflow cluster_sdata {
    if (!params.samplesheet) error "Please provide --samplesheet"

    def inputs = channel
        .fromPath(params.samplesheet)
        .splitCsv(header: true)      // Map(sample, path[, min_counts, min_cells, max_counts_quantile])
        .map { row ->
            if (!row.sample) error "Samplesheet row missing 'sample': ${row}"
            if (!row.path)   error "Samplesheet row missing 'path': ${row}"
            // Per-sample cut: the six samples of a cohort rarely share one threshold,
            // and putting it in the samplesheet records which sample got which value
            // alongside the sample itself rather than in a run's flags.
            tuple(row.sample, file(row.path),
                  resolveThreshold(row.min_counts, params.min_counts),
                  resolveThreshold(row.min_cells,  params.min_cells),
                  resolveThreshold(row.max_counts_quantile, params.max_counts_quantile))
        }                            // tuple(sample, path, min_counts, min_cells, max_counts_quantile)

    // Leiden resolution sweep. Null by default (see nextflow.config) so the list
    // lives only in the clustering script; a `val` process input cannot be null,
    // so pass an empty string and let the module's conditional append omit the flag.
    // Stays a plain val, not a tuple field: the sweep is cohort-wide by design, since
    // comparing resolutions across samples is the whole point of cluster_report.
    def resolutions = params.resolutions ?: ''

    CLUSTER_SDATA(inputs, resolutions)

    // Handoff samplesheet of the clustered zarrs (see create_sdata for rationale).
    CLUSTER_SDATA.out.samplesheet_row
        .map { it.text }             // read row content so collectFile's sort is deterministic
        .collectFile(name: 'cluster_sdata_samplesheet.csv', storeDir: params.outdir,
                     seed: 'sample,path', newLine: true, sort: true)
}

// ── cluster_sdata_gpu ─────────────────────────────────────────────────────────

workflow cluster_sdata_gpu {
    if (!params.samplesheet) error "Please provide --samplesheet"

    def inputs = channel
        .fromPath(params.samplesheet)
        .splitCsv(header: true)      // Map(sample, path[, min_counts, min_cells, max_counts_quantile])
        .map { row ->
            if (!row.sample) error "Samplesheet row missing 'sample': ${row}"
            if (!row.path)   error "Samplesheet row missing 'path': ${row}"
            // See cluster_sdata above for why the cut is per-sample.
            tuple(row.sample, file(row.path),
                  resolveThreshold(row.min_counts, params.min_counts),
                  resolveThreshold(row.min_cells,  params.min_cells),
                  resolveThreshold(row.max_counts_quantile, params.max_counts_quantile))
        }                            // tuple(sample, path, min_counts, min_cells, max_counts_quantile)

    // See cluster_sdata above for why an empty string stands in for "unset".
    def resolutions = params.resolutions ?: ''

    CLUSTER_SDATA_GPU(inputs, resolutions)

    // Handoff samplesheet of the clustered zarrs (see create_sdata for rationale).
    CLUSTER_SDATA_GPU.out.samplesheet_row
        .map { it.text }             // read row content so collectFile's sort is deterministic
        .collectFile(name: 'cluster_sdata_gpu_samplesheet.csv', storeDir: params.outdir,
                     seed: 'sample,path', newLine: true, sort: true)
}

// ── cluster_sdata_gpu_ooc ─────────────────────────────────────────────────────

workflow cluster_sdata_gpu_ooc {
    if (!params.samplesheet) error "Please provide --samplesheet"

    def inputs = channel
        .fromPath(params.samplesheet)
        .splitCsv(header: true)      // Map(sample, path[, min_counts, min_cells, max_counts_quantile])
        .map { row ->
            if (!row.sample) error "Samplesheet row missing 'sample': ${row}"
            if (!row.path)   error "Samplesheet row missing 'path': ${row}"
            // See cluster_sdata above for why the cut is per-sample.
            tuple(row.sample, file(row.path),
                  resolveThreshold(row.min_counts, params.min_counts),
                  resolveThreshold(row.min_cells,  params.min_cells),
                  resolveThreshold(row.max_counts_quantile, params.max_counts_quantile))
        }                            // tuple(sample, path, min_counts, min_cells, max_counts_quantile)

    // HVG selection is off by default (params.n_top_genes = null) so this step
    // matches cluster_sdata/cluster_sdata_gpu. A `val` process input cannot be
    // null, so pass an empty string — the module's conditional append then omits
    // the flag and the script falls back to its own default of no filtering.
    def nTopGenes = params.n_top_genes ?: ''

    // See cluster_sdata above for why an empty string stands in for "unset".
    def resolutions = params.resolutions ?: ''

    CLUSTER_SDATA_GPU_OOC(inputs, params.chunk_size, nTopGenes, resolutions)

    // Handoff samplesheet of the clustered zarrs (see create_sdata for rationale).
    CLUSTER_SDATA_GPU_OOC.out.samplesheet_row
        .map { it.text }             // read row content so collectFile's sort is deterministic
        .collectFile(name: 'cluster_sdata_gpu_ooc_samplesheet.csv', storeDir: params.outdir,
                     seed: 'sample,path', newLine: true, sort: true)
}

// ── create_centroids ──────────────────────────────────────────────────────────

workflow create_centroids {
    if (!params.samplesheet) error "Please provide --samplesheet"

    // Point --samplesheet at a cluster_sdata* handoff sheet — this step reads the
    // clustered zarrs (layers["counts"] and the leiden_res_* columns), not
    // create_sdata's raw ones, which carry no clusters to group by.
    channel
        .fromPath(params.samplesheet)
        .splitCsv(header: true)      // Map(sample, path)
        .map { row ->
            if (!row.sample) error "Samplesheet row missing 'sample': ${row}"
            if (!row.path)   error "Samplesheet row missing 'path': ${row}"
            // row.path rides along as a val as well as a staged path: inside the task the
            // staged name means nothing outside that work dir, and the handoff row has to
            // forward a location sample_summary can still resolve. That deck needs both
            // artifacts — the centroids for expression, the clustered zarr for obs, obsm
            // and the spatial elements — so the row it collects carries both.
            tuple(row.sample, file(row.path), row.path)
        }                            // tuple(sample, path, cluster_path)
        | CREATE_CENTROIDS

    // Handoff samplesheet of the clustered zarrs AND their centroid stores (see
    // create_sdata for rationale). Three columns rather than two, because the decks
    // downstream of this step read both artifacts.
    CREATE_CENTROIDS.out.samplesheet_row
        .map { it.text }             // read row content so collectFile's sort is deterministic
        .collectFile(name: 'create_centroids_samplesheet.csv', storeDir: params.outdir,
                     seed: 'sample,path,centroid_path', newLine: true, sort: true)
}

// ── cluster_report ────────────────────────────────────────────────────────────

workflow cluster_report {
    if (!params.samplesheet) error "Please provide --samplesheet"

    // Point --samplesheet at a cluster_sdata* handoff sheet — this report reads the
    // clustered zarrs (X_umap and the leiden_res_* columns), not create_sdata's raw ones.
    channel
        .fromPath(params.samplesheet)
        .splitCsv(header: true)      // Map(sample, path)
        .map { row ->
            if (!row.path) error "Samplesheet row missing 'path': ${row}"
            file(row.path)
        }                            // path(zarr)
        // sort so the staged order — and therefore the report — is reproducible
        // regardless of the order tasks happen to finish upstream.
        .toSortedList()              // one list of every zarr: fans in to a single task
        .set { clusterReportZarrs }  // val(list of zarr paths)

    CLUSTER_REPORT(
        clusterReportZarrs,
        file("${projectDir}/notebooks/analyze/cluster_report.qmd"),
        file("${projectDir}/resources/ouhsc_ppt_template.pptx"),
    )
}

// ── qc_report ─────────────────────────────────────────────────────────────────

workflow qc_report {
    if (!params.samplesheet) error "Please provide --samplesheet"

    // Point --samplesheet at create_sdata's handoff sheet — this report reads the RAW
    // zarrs, before any filtering. Run against cluster_sdata* output it would still
    // render, but every distribution would already have had the min_counts cut applied,
    // which is the cut the deck exists to choose.
    channel
        .fromPath(params.samplesheet)
        .splitCsv(header: true)      // Map(sample, path)
        .map { row ->
            if (!row.path) error "Samplesheet row missing 'path': ${row}"
            file(row.path)
        }                            // path(zarr)
        // sort so the staged order — and therefore the report — is reproducible
        // regardless of the order tasks happen to finish upstream.
        .toSortedList()              // one list of every zarr: fans in to a single task
        .set { qcReportZarrs }       // val(list of zarr paths)

    QC_REPORT(
        qcReportZarrs,
        file("${projectDir}/notebooks/analyze/qc_report.qmd"),
        file("${projectDir}/resources/ouhsc_ppt_template.pptx"),
    )
}

// ── sample_summary ────────────────────────────────────────────────────────────

workflow sample_summary {
    if (!params.samplesheet)        error "Please provide --samplesheet"
    // Defaulted in nextflow.config to a header-only asset, so this only fires if it is
    // explicitly unset. Its rows override the notebook's per-sample resolution; samples
    // without one fall back to the notebook's default.
    if (!params.chosen_resolutions) error "Please provide --chosen_resolutions"
    // Likewise defaulted in nextflow.config, to the asset built by
    // scripts/convert_ovary_reference.py. Tissue-specific, so it is a path rather than
    // anything derived — a different tissue needs a different marker file.
    if (!params.follicle_markers)   error "Please provide --follicle_markers"
    if (!params.reference_major)    error "Please provide --reference_major"
    // The immune subtype assets, both defaulted in nextflow.config to the same
    // converter's output. Tissue-specific in the same way the two above are.
    if (!params.immune_markers)     error "Please provide --immune_markers"
    if (!params.reference_immune)   error "Please provide --reference_immune"
    // Defaulted to a header-only asset, so this only fires if explicitly unset. Its rows
    // override the argmax per cell type; a sample with no rows takes the defaults.
    if (!params.cell_type_annotations) error "Please provide --cell_type_annotations"
    // Defaulted in nextflow.config to the curated oocyte asset. Tissue-specific in the
    // same way the markers are — it drives the per-stage zoom slides, and a sample with
    // no curated cells in it simply gets none.
    if (!params.curated_oocytes)    error "Please provide --curated_oocytes"

    // Point --samplesheet at a CREATE_CENTROIDS handoff sheet, which carries both
    // locations per sample: `path` is the clustered zarr (read for obs, obsm and the
    // spatial elements) and `centroid_path` the centroid store (read for all expression
    // evidence). A cluster_sdata* sheet no longer suffices — it has no centroid column,
    // and this deck stopped deriving centroids from the counts matrix.
    def sampleSummaryRows = channel
        .fromPath(params.samplesheet)
        .splitCsv(header: true)          // Map(sample, path, centroid_path)
        .map { row ->
            if (!row.path) error "Samplesheet row missing 'path': ${row}"
            if (!row.centroid_path) error(
                "Samplesheet row missing 'centroid_path': ${row}. sample_summary reads " +
                "create_centroids' stores — point --samplesheet at " +
                "<outdir>/create_centroids_samplesheet.csv, not cluster_sdata's sheet.")
            tuple(file(row.path), file(row.centroid_path))
        }                                // tuple(zarr, centroid_h5ad)

    // Split into two fan-in lists rather than one list of pairs: the process takes them
    // as separate `path` inputs so each gets its own stageAs rule. Both are sorted, so
    // the staged order — and therefore the report — is reproducible regardless of the
    // order tasks happen to finish upstream.
    sampleSummaryRows.map { zarr, _centroid -> zarr }
        .toSortedList()
        .set { sampleSummaryZarrs }      // val(list of zarr paths)
    sampleSummaryRows.map { _zarr, centroid -> centroid }
        .toSortedList()
        .set { sampleSummaryCentroids }  // val(list of centroid h5ad paths)

    SAMPLE_SUMMARY(
        sampleSummaryZarrs,
        sampleSummaryCentroids,
        file("${projectDir}/notebooks/analyze/sample_summary.qmd"),
        file("${projectDir}/resources/ouhsc_ppt_template.pptx"),
        file(params.chosen_resolutions),
        file(params.follicle_markers),
        file(params.reference_major),
        file(params.immune_markers),
        file(params.reference_immune),
        file(params.cell_type_annotations),
        file(params.curated_oocytes),
    )
}

// ── concat_sdata ──────────────────────────────────────────────────────────────

workflow concat_sdata {
    if (!params.samplesheet) error "Please provide --samplesheet"

    channel
        .fromPath(params.samplesheet)
        .splitCsv(header: true)  // Map(path, ...)
        .map { row ->
            if (!row.path) error "Samplesheet row missing 'path': ${row}"
            file(row.path)
        }                        // path
        .collect()               // List<path>
        | CONCAT_SDATA

    // Handoff samplesheet for the merged zarr (see create_sdata for rationale).
    CONCAT_SDATA.out.samplesheet_row
        .map { it.text }             // read row content so collectFile's sort is deterministic
        .collectFile(name: 'concat_sdata_samplesheet.csv', storeDir: params.outdir,
                     seed: 'sample,path', newLine: true, sort: true)
}

// ── downsample_sdata ──────────────────────────────────────────────────────────

workflow downsample_sdata {
    if (!params.samplesheet) error "Please provide --samplesheet"
    if (!params.fraction && !params.n_cells) error "Please provide --fraction or --n_cells"

    channel
        .fromPath(params.samplesheet)
        .splitCsv(header: true)      // Map(sample, path)
        .map { row ->
            if (!row.sample) error "Samplesheet row missing 'sample': ${row}"
            if (!row.path)   error "Samplesheet row missing 'path': ${row}"
            tuple(row.sample, file(row.path))
        }                            // tuple(sample, path)
        | DOWNSAMPLE_SDATA

    // Handoff samplesheet of the downsampled zarrs (see create_sdata for rationale).
    DOWNSAMPLE_SDATA.out.samplesheet_row
        .map { it.text }             // read row content so collectFile's sort is deterministic
        .collectFile(name: 'downsample_sdata_samplesheet.csv', storeDir: params.outdir,
                     seed: 'sample,path', newLine: true, sort: true)
}

// ── plot_follicle ─────────────────────────────────────────────────────────────

workflow plot_follicle {
    if (!params.samplesheet) error "Please provide --samplesheet"

    def plotFollicleNotebook = file("${projectDir}/notebooks/analyze/plot_follicle.qmd")
    def timerScript          = file("${projectDir}/bin/timer.py")

    def plotInputs = channel
        .fromPath(params.samplesheet)
        .splitCsv(header: true)        // Map(sample, cell, path)
        .map { row ->
            if (!row.sample) error "Samplesheet row missing 'sample': ${row}"
            if (!row.path)   error "Samplesheet row missing 'path': ${row}"
            def follicleId = "${row.sample}_${row.cell}"
            tuple(follicleId, row.sample, file(row.path), paramsFile(follicleId, plotFollicleNotebook, row))
        }                              // tuple(follicle_id, sample, path, params_yml)

    PLOT_FOLLICLE(plotInputs, plotFollicleNotebook, timerScript)
}
