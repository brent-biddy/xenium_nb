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
    def analyzeMode = (params.analyze ?: 'all').toLowerCase()

    if (!params.samplesheet) error "Please provide --samplesheet"
    if (!(analyzeMode in ['plot_follicle', 'all'])) {
        error "Invalid --analyze '${analyzeMode}'. Valid values are: plot_follicle, all"
    }

    def timerScript = file("${projectDir}/bin/timer.py")
    def analysisRegistry = new groovy.json.JsonSlurper()
        .parse(new File("${projectDir}/assets/notebook_registry.json"))
        .analysis

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
        def notebook = file("${projectDir}/notebooks/plot_follicle.qmd")

        rowsList
            .map { sample, stagedPath, rowMap ->
                def cell = rowMap.cell
                def sampleId = "${sample}_${cell}"
                tuple(sampleId, sample, cell, stagedPath, paramsFile(sampleId, analysisRegistry.plot_follicle.params, rowMap))
            }
            .set { plotInputs } // tuple(sample_cell_id, sample, cell, staged_path, params_yml)
        PLOT_FOLLICLE(plotInputs, notebook, timerScript)
    }
}
