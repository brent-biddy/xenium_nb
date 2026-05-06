#!/usr/bin/env nextflow

// Renders analysis Quarto notebooks against pre-built artifacts (typically
// produced by create.nf). The samplesheet must carry whatever per-row params
// the selected notebooks declare (e.g. `cell` for plot_follicle).
//
// --analyze accepts 'all', a comma-separated list of notebook IDs, or a
// Groovy list. IDs are resolved against NotebookRegistry.analysis().

nextflow.enable.dsl = 2

include { WRITE_QUARTO_PARAMS } from './modules/write_quarto_params'
include { RUN_NOTEBOOK } from './modules/run_notebook'

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
    def registry = NotebookRegistry.analysis(projectDir.toString())
    def notebookIds = resolveNotebookIds(params.analyze, registry)
    log.info "Running analysis notebooks: ${notebookIds.join(', ')}"

    // Params that RUN_NOTEBOOK injects from pipeline config rather than the
    // samplesheet, so they don't need to appear as samplesheet columns.
    def PIPELINE_PARAM_KEYS = ['cell_ids_file', 'radius', 'n_jobs'] as Set

    // Union of declared notebook params, minus those Nextflow injects, becomes
    // the set of samplesheet columns each row must populate.
    def notebookParamKeys = notebookIds
        .collectMany { registry[it].params ?: [] }
        .collect { it.toString() } as Set
    def requiredColumns = ((['sample', 'path'] as Set) + (notebookParamKeys - PIPELINE_PARAM_KEYS)) as Set

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
    // TODO: split notebook staging so analysis no longer carries cell_ids_file.
    def cellIdsFile = file(params.cell_ids_file)
    def notebookSpecs = notebookIds.collect { id ->
        [path: file(registry[id].path), params: registry[id].params ?: []]
    }
    def notebookChannel = Channel.fromList(notebookSpecs)

    rows
        .combine(notebookChannel)
        .map { sample, artifactPath, rowParams, spec ->
            // Notebooks with a 'cell' param fan out per-cell, so disambiguate
            // their published filenames with the cell ID.
            def publishSample = spec.params.contains('cell')
                ? "${sample}_${rowParams.cell}"
                : sample
            def publishDir = "${params.outdir}/${sample}/${spec.path.baseName}"
            tuple(
                spec.path.toString(),
                spec.path.baseName,
                timerScript.toString(),
                artifactPath.toString(),
                sample,
                publishDir,
                publishSample,
                rowParams,
                cellIdsFile,
                spec.params
            )
        }
        | WRITE_QUARTO_PARAMS
        | RUN_NOTEBOOK
}
