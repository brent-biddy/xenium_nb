#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { RUN_NOTEBOOK } from './modules/run_notebook'

workflow {
    if (!params.samplesheet) error "Please provide --samplesheet"
    def expectedColumns = ['sample', 'path'] as Set

    def normalizeNotebookSelection = { value ->
        if (value == null) {
            return []
        }
        if (value instanceof List) {
            return value.collect { it.toString().trim() }.findAll { it }
        }
        value
            .toString()
            .split(',')
            .collect { it.trim() }
            .findAll { it }
    }

    def notebookIds = normalizeNotebookSelection(params.notebooks)
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
    def rows = Channel
        .fromPath(params.samplesheet)
        .splitCsv(header: true)
        .map { row ->
            def columns = row.keySet().findAll { it != null && it != '' } as Set
            if (columns != expectedColumns) {
                error "Analysis samplesheet must contain exactly these columns: sample,path. Found: ${columns.join(',')}"
            }
            if (!row.sample) error "Samplesheet row missing 'sample': ${row}"
            if (!row.path)   error "Samplesheet row missing 'path': ${row}"
            def rowMap = new LinkedHashMap(row)
            def sample = rowMap.sample.toString()
            tuple(sample, file(rowMap.path), rowMap)
        }

    def timerScript = file("${projectDir}/bin/timer.py")
    def notebookSpecs = notebookIds.collect { id ->
        [id: id, scope: registry[id].scope, path: file(registry[id].path)]
    }
    def notebookChannel = Channel.of(*notebookSpecs)

    rows
        .combine(notebookChannel)
        .map { sample, artifactPath, rowParams, spec ->
            def parentSample = artifactPath.parent.parent.parent.name
            def publishDir = "${params.outdir}/${parentSample}/${spec.path.baseName}"
            def outputName = "${sample}_${spec.path.baseName}.html"
            tuple(spec.path, timerScript, artifactPath, sample, publishDir, outputName, rowParams)
        }
        | RUN_NOTEBOOK
}
