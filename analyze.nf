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
            tuple(row.sample.toString(), file(row.path), row)
        }
        .collect(flat: false)
        .set { rowsList } // List<tuple(sample, staged_path, row_map)>

    // ---- plot_follicle: per-cell follicle plots ----
    if (analyzeMode == 'plot_follicle' || analyzeMode == 'all') {
        def notebook = file("${projectDir}/notebooks/plot_follicle.qmd")

        rowsList
            .flatMap { rows ->
                rows.collect { row ->
                    def sample = row[0]
                    def artifactPath = row[1]
                    def rowParams = row[2]
                    def cell = rowParams.cell.toString()
                    def sampleId = "${sample}_${cell}"
                    tuple(sampleId, artifactPath, rowParams, analysisRegistry.plot_follicle.params)
                }
            }
            .set { plotParamsInputs } // tuple(sample_cell_id, staged_path, row_map, declared_params)
        PLOT_FOLLICLE_PARAMS(plotParamsInputs) | set { plotParams } // tuple(sample_cell_id, params_yml)

        rowsList
            .flatMap { rows ->
                rows.collect { row ->
                    def sample = row[0]
                    def artifactPath = row[1]
                    def rowParams = row[2]
                    def cell = rowParams.cell.toString()
                    def sampleId = "${sample}_${cell}"
                    tuple(sampleId, sample, cell, artifactPath)
                }
            }
            .join(plotParams.params_file)
            .set { plotInputs } // tuple(sample_cell_id, sample, cell, staged_path, params_yml)
        PLOT_FOLLICLE(plotInputs, notebook, timerScript)
    }
}
