#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { RUN_NOTEBOOK } from './modules/run_notebook'

workflow {
    if (!params.samplesheet) error "Please provide --samplesheet"

    def notebookIds = (
        params.notebooks == null
            ? []
            : params.notebooks instanceof List
                ? params.notebooks.collect { it.toString().trim() }.findAll { it }
                : params.notebooks.toString().split(',').collect { it.trim() }.findAll { it }
    )
    if (!notebookIds) {
        error "Please provide at least one notebook ID via --notebooks"
    }

    def registry = params.analysis_notebook_registry ?: [:]
    def unknownNotebookIds = notebookIds.findAll { !registry.containsKey(it) }.unique()
    if (unknownNotebookIds) {
        error "Unknown notebook IDs: ${unknownNotebookIds.join(', ')}. Known IDs: ${registry.keySet().sort().join(', ')}"
    }

    def scopes = notebookIds.collect { registry[it].scope }.unique()
    if (scopes.size() != 1) {
        error "All selected notebooks must share one scope. Requested scopes: ${scopes.join(', ')}"
    }
    def scope = scopes[0]
    def requiredColumns = scope == 'follicle'
        ? ['sample', 'cell', 'path'] as Set
        : ['sample', 'path'] as Set

    def rows = Channel
        .fromPath(params.samplesheet)
        .splitCsv(header: true)
        .map { row ->
            def columns = row.keySet().findAll { it != null && it != '' } as Set
            if (!columns.containsAll(requiredColumns)) {
                error "Analysis samplesheet must contain at least these columns: ${requiredColumns.join(',')}. Found: ${columns.join(',')}"
            }
            if (!row.sample) error "Samplesheet row missing 'sample': ${row}"
            if (scope == 'follicle' && !row.cell) error "Samplesheet row missing 'cell': ${row}"
            if (!row.path)   error "Samplesheet row missing 'path': ${row}"
            def rowMap = new LinkedHashMap(row)
            def sample = rowMap.sample.toString()
            tuple(sample, file(rowMap.path), rowMap)
        }

    def timerScript = file("${projectDir}/bin/timer.py")
    def notebookSpecs = notebookIds.collect { id ->
        [
            id   : id,
            scope: registry[id].scope,
            path : file(registry[id].path),
            params: registry[id].params ?: [],
        ]
    }
    def notebookChannel = Channel.fromList(notebookSpecs)

    rows
        .combine(notebookChannel)
        .map { sample, artifactPath, rowParams, spec ->
            def publishSample = spec.scope == 'follicle'
                ? "${sample}_${rowParams.cell}"
                : sample
            def publishDir = "${params.outdir}/${sample}/${spec.path.baseName}"
            tuple(spec.path, timerScript, artifactPath, publishSample, publishDir, rowParams, spec.params)
        }
        | RUN_NOTEBOOK
}
