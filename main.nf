#!/usr/bin/env nextflow

// Single entry point for all pipeline steps.
// Select a step with: nextflow run main.nf --step <name> --samplesheet <path>
//
// Steps:
//   downsample_xenium_region  samplesheet: sample, path, xmin, ymin, xmax, ymax[, region_name, he_image, he_alignment]
//   create_sdata              samplesheet: sample, path[, he_image, he_alignment]
//   create_follicle_sdata     samplesheet: sample, path  (+ --cell_ids_file)
//   cluster_sdata             samplesheet: sample, path
//   cluster_sdata_gpu         samplesheet: sample, path
//   cluster_sdata_gpu_ooc     samplesheet: sample, path  (+ --chunk_size, --n_top_genes)
//   concat_sdata              samplesheet: path
//   downsample_sdata          samplesheet: sample, path  (+ --fraction or --n_cells)
//   plot_follicle             samplesheet: sample, cell, path

include { DOWNSAMPLE_XENIUM_REGION } from './modules/downsample_xenium_region'
include { CREATE_SDATA }             from './modules/create_sdata'
include { CREATE_FOLLICLE_SDATA }    from './modules/create_follicle_sdata'
include { CLUSTER_SDATA }            from './modules/cluster_sdata'
include { CLUSTER_SDATA_GPU }        from './modules/cluster_sdata_gpu'
include { CLUSTER_SDATA_GPU_OOC }    from './modules/cluster_sdata_gpu_ooc'
include { CONCAT_SDATA }             from './modules/concat_sdata'
include { DOWNSAMPLE_SDATA }         from './modules/downsample_sdata'
include { PLOT_FOLLICLE }            from './modules/plot_follicle'
include { paramsFile }               from './modules/quarto_params'

// ── Entry workflow ────────────────────────────────────────────────────────────

workflow {
    if (!params.step) error "Please provide --step <name>. Valid steps: downsample_xenium_region, create_sdata, create_follicle_sdata, cluster_sdata, cluster_sdata_gpu, cluster_sdata_gpu_ooc, concat_sdata, downsample_sdata, plot_follicle"

    if      (params.step == 'downsample_xenium_region')  downsample_xenium_region()
    else if (params.step == 'create_sdata')              create_sdata()
    else if (params.step == 'create_follicle_sdata')     create_follicle_sdata()
    else if (params.step == 'cluster_sdata')             cluster_sdata()
    else if (params.step == 'cluster_sdata_gpu')         cluster_sdata_gpu()
    else if (params.step == 'cluster_sdata_gpu_ooc')     cluster_sdata_gpu_ooc()
    else if (params.step == 'concat_sdata')              concat_sdata()
    else if (params.step == 'downsample_sdata')          downsample_sdata()
    else if (params.step == 'plot_follicle')             plot_follicle()
    else error "Unknown --step '${params.step}'. Valid steps: downsample_xenium_region, create_sdata, create_follicle_sdata, cluster_sdata, cluster_sdata_gpu, cluster_sdata_gpu_ooc, concat_sdata, downsample_sdata, plot_follicle"
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
            def heImage    = row.he_image     ? new File(row.he_image     as String).absolutePath : ""
            def heAlign    = row.he_alignment ? new File(row.he_alignment as String).absolutePath : ""
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
            def heImage = row.he_image     ? new File(row.he_image     as String).absolutePath : ""
            def heAlign = row.he_alignment ? new File(row.he_alignment as String).absolutePath : ""
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

    channel
        .fromPath(params.samplesheet)
        .splitCsv(header: true)      // Map(sample, path)
        .map { row ->
            if (!row.sample) error "Samplesheet row missing 'sample': ${row}"
            if (!row.path)   error "Samplesheet row missing 'path': ${row}"
            tuple(row.sample, file(row.path))
        }                            // tuple(sample, path)
        | CLUSTER_SDATA

    // Handoff samplesheet of the clustered zarrs (see create_sdata for rationale).
    CLUSTER_SDATA.out.samplesheet_row
        .map { it.text }             // read row content so collectFile's sort is deterministic
        .collectFile(name: 'cluster_sdata_samplesheet.csv', storeDir: params.outdir,
                     seed: 'sample,path', newLine: true, sort: true)
}

// ── cluster_sdata_gpu ─────────────────────────────────────────────────────────

workflow cluster_sdata_gpu {
    if (!params.samplesheet) error "Please provide --samplesheet"

    channel
        .fromPath(params.samplesheet)
        .splitCsv(header: true)      // Map(sample, path)
        .map { row ->
            if (!row.sample) error "Samplesheet row missing 'sample': ${row}"
            if (!row.path)   error "Samplesheet row missing 'path': ${row}"
            tuple(row.sample, file(row.path))
        }                            // tuple(sample, path)
        | CLUSTER_SDATA_GPU

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
        .splitCsv(header: true)      // Map(sample, path)
        .map { row ->
            if (!row.sample) error "Samplesheet row missing 'sample': ${row}"
            if (!row.path)   error "Samplesheet row missing 'path': ${row}"
            tuple(row.sample, file(row.path))
        }                            // tuple(sample, path)

    // HVG selection is off by default (params.n_top_genes = null) so this step
    // matches cluster_sdata/cluster_sdata_gpu. A `val` process input cannot be
    // null, so pass an empty string — the module's conditional append then omits
    // the flag and the script falls back to its own default of no filtering.
    def nTopGenes = params.n_top_genes ?: ''

    CLUSTER_SDATA_GPU_OOC(inputs, params.chunk_size, nTopGenes)

    // Handoff samplesheet of the clustered zarrs (see create_sdata for rationale).
    CLUSTER_SDATA_GPU_OOC.out.samplesheet_row
        .map { it.text }             // read row content so collectFile's sort is deterministic
        .collectFile(name: 'cluster_sdata_gpu_ooc_samplesheet.csv', storeDir: params.outdir,
                     seed: 'sample,path', newLine: true, sort: true)
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
