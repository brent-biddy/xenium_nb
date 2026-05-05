#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { RUN_NOTEBOOK } from './modules/run_notebook'

workflow {
    if (!params.samplesheet) error "Please provide --samplesheet"
    def registry = NotebookRegistry.analysis(projectDir.toString())
    def pipelineParamKeys = ['cell_ids_file', 'radius', 'n_jobs'] as Set

    def analyzeSelection = params.analyze == null ? 'all' : params.analyze
    def notebookIds = (
        analyzeSelection instanceof List
            ? analyzeSelection.collect { it.toString().trim() }.findAll { it }
            : analyzeSelection.toString().trim().equalsIgnoreCase('all')
                ? (registry.keySet() as List).sort()
                : analyzeSelection.toString().split(',').collect { it.trim() }.findAll { it }
    )
    if (!notebookIds) {
        error "Please provide at least one analysis notebook ID via --analyze, or use 'all'"
    }
    def unknownNotebookIds = notebookIds.findAll { !registry.containsKey(it) }.unique()
    if (unknownNotebookIds) {
        error "Unknown analysis notebook IDs: ${unknownNotebookIds.join(', ')}. Known IDs: ${registry.keySet().sort().join(', ')}"
    }
    log.info "Running analysis notebooks: ${notebookIds.join(', ')}"

    def notebookParamKeys = notebookIds
        .collectMany { registry[it].params ?: [] }
        .collect { it.toString() } as Set
    def requiredColumns = ((['sample', 'path'] as Set) + (notebookParamKeys - pipelineParamKeys)) as Set

    def rows = Channel
        .fromPath(params.samplesheet)
        .splitCsv(header: true)
        .map { row ->
            def columns = row.keySet().findAll { it != null && it != '' } as Set
            if (!columns.containsAll(requiredColumns)) {
                error "Analysis samplesheet must contain at least these columns: ${requiredColumns.join(',')}. Found: ${columns.join(',')}"
            }
            if (!row.sample) error "Samplesheet row missing 'sample': ${row}"
            if (!row.path)   error "Samplesheet row missing 'path': ${row}"
            requiredColumns.each { col ->
                if (!row[col]) {
                    error "Samplesheet row missing '${col}': ${row}"
                }
            }
            def rowMap = new LinkedHashMap(row)
            def sample = rowMap.sample.toString()
            tuple(sample, file(rowMap.path), rowMap)
        }

    def timerScript = file("${projectDir}/bin/timer.py")
    def notebookSpecs = notebookIds.collect { id ->
        [
            id    : id,
            path  : file(registry[id].path),
            params: registry[id].params ?: [],
        ]
    }
    def notebookChannel = Channel.fromList(notebookSpecs)

    rows
        .combine(notebookChannel)
        .map { sample, artifactPath, rowParams, spec ->
            def publishSample = spec.params.contains('cell')
                ? "${sample}_${rowParams.cell}"
                : sample
            def publishDir = "${params.outdir}/${sample}/${spec.path.baseName}"
            tuple(spec.path, timerScript, artifactPath, publishSample, publishDir, rowParams, spec.params)
        }
        | RUN_NOTEBOOK
}
