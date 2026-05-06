#!/usr/bin/env nextflow

// Renders analysis Quarto notebooks against pre-built artifacts (typically
// produced by create.nf). The samplesheet must carry whatever per-row params
// the selected notebooks declare (e.g. `cell` for plot_follicle).
//
// --analyze accepts 'all' or a notebook ID.

nextflow.enable.dsl = 2

include { WRITE_QUARTO_PARAMS as PLOT_FOLLICLE_PARAMS } from './modules/write_quarto_params'
include { PLOT_FOLLICLE } from './modules/analyze_notebooks'

workflow {
    def analyzeMode = (params.analyze ?: 'all').toLowerCase()

    if (!params.samplesheet) error "Please provide --samplesheet"
    if (!(analyzeMode in ['plot_follicle', 'all'])) {
        error "Invalid --analyze '${analyzeMode}'. Valid values are: plot_follicle, all"
    }

    def timerScript = Channel.fromPath("${projectDir}/bin/timer.py")

    // ---- samplesheet ----
    Channel
        .fromPath(params.samplesheet)
        .splitCsv(header: true)
        .map { row ->
            if (!row.sample) error "Analysis samplesheet row missing 'sample': ${row}"
            if (!row.path)   error "Analysis samplesheet row missing 'path': ${row}"
            tuple(row.sample.toString(), file(row.path), row)
        }
        .set { rows }

    // ---- plot_follicle: per-cell follicle plots ----
    if (analyzeMode == 'plot_follicle' || analyzeMode == 'all') {
        def notebook = Channel.fromPath("${projectDir}/notebooks/plot_follicle.qmd")

        rows
            .map { sample, artifactPath, rowParams ->
                def cell = rowParams.cell.toString()
                def sampleId = "${sample}_${cell}"
                tuple(sampleId, artifactPath, rowParams, ['sample', 'cell'])
            }
            .set { plotParamsInputs }
        PLOT_FOLLICLE_PARAMS(plotParamsInputs) | set { plotParams }

        rows
            .map { sample, artifactPath, rowParams ->
                def cell = rowParams.cell.toString()
                def sampleId = "${sample}_${cell}"
                tuple(sampleId, sample, cell, artifactPath)
            }
            .join(plotParams.params_file)
            .set { plotInputs }
        PLOT_FOLLICLE(plotInputs, notebook, timerScript)
    }
}
