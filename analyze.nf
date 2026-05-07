#!/usr/bin/env nextflow

// Renders analysis Quarto notebooks against pre-built artifacts (typically
// produced by create.nf). The samplesheet must carry whatever per-row params
// the selected notebooks declare (e.g. `cell` for plot_follicle).
//
// --analyze accepts 'all' or a notebook ID.

nextflow.enable.dsl = 2

include { paramsFile } from './modules/quarto_params'
include { PLOT_FOLLICLE } from './modules/analyze_notebooks'

workflow {
    if (!params.samplesheet) error "Please provide --samplesheet"
    if (!params.analyze)     error "Please provide --analyze (plot_follicle, all)"

    def analyzeMode = params.analyze.toLowerCase()

    if (!(analyzeMode in ['plot_follicle', 'all'])) {
        error "Invalid --analyze '${analyzeMode}'. Valid values are: plot_follicle, all"
    }

    def timerScript = file("${projectDir}/bin/timer.py")

    // ---- samplesheet ----
    Channel
        .fromPath(params.samplesheet)
        .splitCsv(header: true)
        .map { row ->
            if (!row.sample) error "Analysis samplesheet row missing 'sample': ${row}"
            if (!row.path)   error "Analysis samplesheet row missing 'path': ${row}"
            tuple(row.sample, file(row.path), row)
        }
        .set { rowsList } // tuple(sample, staged_path, row_map)

    // ---- plot_follicle: per-cell follicle plots ----
    if (analyzeMode == 'plot_follicle' || analyzeMode == 'all') {
        def plotFollicleNotebook = file("${projectDir}/notebooks/analyze_plot_follicle.qmd")

        rowsList
            .map { sample, stagedPath, rowMap ->
                def follicleId = "${sample}_${rowMap.cell}"
                tuple(follicleId, sample, stagedPath, paramsFile(follicleId, plotFollicleNotebook, rowMap))
            }
            .set { plotInputs } // tuple(follicle_id, sample, staged_path, params_yml)
        PLOT_FOLLICLE(plotInputs, plotFollicleNotebook, timerScript)
    }
}
