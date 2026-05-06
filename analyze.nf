#!/usr/bin/env nextflow

// Renders analysis Quarto notebooks against pre-built artifacts (typically
// produced by create.nf). The samplesheet must carry whatever per-row params
// the selected notebooks declare (e.g. `cell` for plot_follicle).
//
// --analyze accepts 'all', a comma-separated list of notebook IDs, or a
// Groovy list.

nextflow.enable.dsl = 2

include { WRITE_QUARTO_PARAMS as WRITE_PLOT_FOLLICLE_PARAMS } from './modules/write_quarto_params'
include { PLOT_FOLLICLE } from './modules/analyze_notebooks'

// Resolves --analyze input ('all' / comma-separated / List) into a list of
// known notebook IDs, validated against the registry.
def resolveNotebookIds(selection, registry) {
    def raw = selection == null ? 'all' : selection
    def ids = raw instanceof List
        ? raw.collect { it.toString().trim() }.findAll { it }
        : raw.toString().trim().equalsIgnoreCase('all')
            ? (registry.keySet() as List).sort()
            : raw.toString().split(',').collect { it.trim() }.findAll { it }
    if (!ids) {
        error "Please provide at least one analysis notebook ID via --analyze, or use 'all'"
    }
    def unknown = ids.findAll { !registry.containsKey(it) }.unique()
    if (unknown) {
        error "Unknown analysis notebook IDs: ${unknown.join(', ')}. Known IDs: ${registry.keySet().sort().join(', ')}"
    }
    ids
}

workflow {
    if (!params.samplesheet) error "Please provide --samplesheet"
    def registry = [
        plot_follicle: [
            path            : file("${projectDir}/notebooks/plot_follicle.qmd"),
            required_columns: ['sample', 'path', 'cell'],
        ],
    ]
    def notebookIds = resolveNotebookIds(params.analyze, registry)
    log.info "Running analysis notebooks: ${notebookIds.join(', ')}"

    // Union of selected notebook requirements becomes the samplesheet contract.
    def requiredColumns = notebookIds
        .collectMany { registry[it].required_columns ?: [] }
        .collect { it.toString() } as Set

    def rows = Channel
        .fromPath(params.samplesheet)
        .splitCsv(header: true)
        .map { row ->
            def columns = row.keySet().findAll { it != null && it != '' } as Set
            if (!columns.containsAll(requiredColumns)) {
                error "Analysis samplesheet must contain at least these columns: ${requiredColumns.join(',')}. Found: ${columns.join(',')}"
            }
            requiredColumns.each { col ->
                if (!row[col]) error "Samplesheet row missing '${col}': ${row}"
            }
            tuple(row.sample.toString(), file(row.path), row)
        }

    def timerScript = file("${projectDir}/bin/timer.py")
    if ('plot_follicle' in notebookIds) {
        def notebook = registry.plot_follicle.path
        def plotParamsInputs = rows.map { sample, artifactPath, rowParams ->
            def cell = rowParams.cell.toString()
            def sampleId = "${sample}_${cell}"
            tuple(sampleId, artifactPath, rowParams + [sample: sample, cell: cell], ['sample', 'cell', 'path'])
        }
        def plotRunInputs = rows.map { sample, artifactPath, rowParams ->
            def cell = rowParams.cell.toString()
            def sampleId = "${sample}_${cell}"
            tuple(sampleId, sample, cell, artifactPath)
        }

        def plotParams = WRITE_PLOT_FOLLICLE_PARAMS(plotParamsInputs)
        PLOT_FOLLICLE(
            plotRunInputs.join(plotParams.params_file),
            Channel.value(notebook),
            Channel.value(timerScript),
        )
    }
}
