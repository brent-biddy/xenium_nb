#!/usr/bin/env nextflow

// Builds Xenium SpatialData artifacts and the samplesheets that downstream
// workflows (analyze.nf) consume. Modes:
//   sdata           - render create_sdata.qmd per ROI from raw Xenium output
//   follicle_sdata  - render create_follicle_sdata.qmd per cell ID from existing sample zarrs
//   all             - run both, chaining sdata outputs into the follicle step

nextflow.enable.dsl = 2

include { RUN_CREATE_NOTEBOOK as CREATE_SDATA } from './modules/run_create_notebook'
include { RUN_CREATE_NOTEBOOK as CREATE_FOLLICLE_SDATA } from './modules/run_create_notebook'
include { WRITE_SAMPLESHEET as WRITE_SAMPLE_ANALYSIS_INPUTS } from './modules/write_samplesheet'
include { WRITE_SAMPLESHEET as WRITE_FOLLICLE_ANALYSIS_INPUTS } from './modules/write_samplesheet'

// Reads a (sample, path, ...) CSV and emits (sample, file, rowMap) per row.
// Extra columns are preserved in rowMap so notebooks can read per-sample params.
def parseSamplesheet(sheetPath, label) {
    Channel
        .fromPath(sheetPath)
        .splitCsv(header: true)
        .map { row ->
            if (!row.sample) error "${label} row missing 'sample': ${row}"
            if (!row.path)   error "${label} row missing 'path': ${row}"
            tuple(row.sample.toString(), file(row.path), row)
        }
}

// Collects per-sample row maps into a single JSON-encoded payload tuple
// shaped for WRITE_SAMPLESHEET.
def buildSamplesheetInput(rowsChannel, outputName, publishDir) {
    rowsChannel
        .collect()
        .map { rows -> groovy.json.JsonOutput.toJson(rows) }
        .map { rowsJson -> tuple(outputName, rowsJson, publishDir) }
}

workflow {
    if (!params.samplesheet) error "Please provide --samplesheet"
    def producerRegistry = NotebookRegistry.producer(projectDir.toString())
    def cellIdsFilePath = file(params.cell_ids_file)
    def timerScript = file("${projectDir}/bin/timer.py")
    def createMode = params.create.toLowerCase()
    if (!(createMode in ['sdata', 'follicle_sdata', 'all'])) {
        error "Invalid create '${createMode}'. Valid values are: sdata, follicle_sdata, all"
    }

    def sampleArtifacts = null

    // ---- sdata: raw Xenium -> per-sample SpatialData zarr ----
    if (createMode in ['sdata', 'all']) {
        def sampleRows = parseSamplesheet(params.samplesheet, 'Create samplesheet')

        def createNotebook = file(producerRegistry.create_sdata.path)
        def createInputs = sampleRows.map { sample, inputPath, rowParams ->
            def publishDir = "${params.outdir}/${sample}/${createNotebook.baseName}"
            tuple(createNotebook, timerScript, inputPath, cellIdsFilePath, sample, publishDir, rowParams, producerRegistry.create_sdata.params)
        }

        sampleArtifacts = CREATE_SDATA(createInputs).artifacts

        def sampleArtifactRows = sampleArtifacts.map { sample, sampleZarr, rowParams ->
            [
                sample: sample,
                path  : "${params.outdir}/${sample}/${createNotebook.baseName}/output/${sample}.zarr",
            ]
        }

        WRITE_SAMPLE_ANALYSIS_INPUTS(buildSamplesheetInput(
            sampleArtifactRows,
            'sample_sdata_samplesheet.csv',
            "${params.outdir}/${createNotebook.baseName}",
        ))
    }

    // Skip CREATE_SDATA: caller's samplesheet already points at existing sample zarrs.
    if (createMode == 'follicle_sdata') {
        sampleArtifacts = parseSamplesheet(params.samplesheet, 'Create follicle samplesheet')
    }

    // ---- follicle_sdata: per-sample SpatialData -> per-cell-ID subset zarrs ----
    if (createMode in ['follicle_sdata', 'all']) {
        def subsetNotebook = file(producerRegistry.create_follicle_sdata.path)
        def subsetInputs = sampleArtifacts.map { sample, sampleZarr, rowParams ->
            def publishDir = "${params.outdir}/${sample}/${subsetNotebook.baseName}"
            tuple(subsetNotebook, timerScript, sampleZarr, cellIdsFilePath, sample, publishDir, rowParams, producerRegistry.create_follicle_sdata.params)
        }

        def follicleArtifacts = CREATE_FOLLICLE_SDATA(subsetInputs).artifacts

        def follicleArtifactRows = follicleArtifacts.flatMap { sample, zarrPaths, _rowParams ->
            // Nextflow emits a single Path for one match and a List<Path> for many; normalize.
            def zarrs = zarrPaths instanceof List ? zarrPaths : [zarrPaths]
            zarrs.collect { zarr ->
                [
                    sample: sample,
                    cell  : zarr.baseName,
                    path  : "${params.outdir}/${sample}/${subsetNotebook.baseName}/output/${zarr.baseName}.zarr",
                ]
            }
        }

        WRITE_FOLLICLE_ANALYSIS_INPUTS(buildSamplesheetInput(
            follicleArtifactRows,
            'follicle_sdata_samplesheet.csv',
            "${params.outdir}/${subsetNotebook.baseName}",
        ))
    }
}
