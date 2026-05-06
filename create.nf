#!/usr/bin/env nextflow

// Builds Xenium SpatialData artifacts and the samplesheets that downstream
// workflows (analyze.nf) consume. Modes:
//   sdata           - render create_sdata.qmd per ROI from raw Xenium output
//   follicle_sdata  - render create_follicle_sdata.qmd per cell ID from existing sample zarrs
//   all             - run both, chaining sdata outputs into the follicle step

nextflow.enable.dsl = 2

include { WRITE_QUARTO_PARAMS as SDATA_QUARTO_PARAMS } from './modules/write_quarto_params'
include { WRITE_QUARTO_PARAMS as FOLLICLE_SDATA_QUARTO_PARAMS } from './modules/write_quarto_params'
include { RUN_NOTEBOOK as CREATE_SDATA } from './modules/run_notebook'
include { RUN_NOTEBOOK as CREATE_FOLLICLE_SDATA } from './modules/run_notebook'
include { WRITE_SAMPLESHEET as WRITE_SDATA_SAMPLESHEET } from './modules/write_samplesheet'
include { WRITE_SAMPLESHEET as WRITE_FOLLICLE_SAMPLESHEET } from './modules/write_samplesheet'

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
    def createRegistry = NotebookRegistry.create(projectDir.toString())
    def timerScript = file("${projectDir}/bin/timer.py")
    def createMode = params.create.toLowerCase()
    if (!(createMode in ['sdata', 'follicle_sdata', 'all'])) {
        error "Invalid create '${createMode}'. Valid values are: sdata, follicle_sdata, all"
    }

    def sampleArtifacts = null

    // ---- sdata: raw Xenium -> per-sample SpatialData zarr ----
    if (createMode in ['sdata', 'all']) {
        def sampleRows = parseSamplesheet(params.samplesheet, 'Create samplesheet')

        def createNotebook = file(createRegistry.create_sdata.path)
        def cellIdsFile = file(params.cell_ids_file)
        def createInputs = sampleRows.map { sample, inputPath, rowParams ->
            tuple(
                createNotebook.toString(),
                createNotebook.baseName,
                timerScript.toString(),
                inputPath.toString(),
                sample,
                "${params.outdir}/${sample}/${createNotebook.baseName}",
                sample,
                rowParams,
                cellIdsFile,
                createRegistry.create_sdata.params
            )
        }

        sampleArtifacts = CREATE_SDATA(SDATA_QUARTO_PARAMS(createInputs).notebook_inputs).artifacts

        def sampleArtifactRows = sampleArtifacts.map { sample, sampleZarr, _rowParams ->
            def imageScaleFactor = _rowParams.image_scale_factor ?: 1.0
            [
                sample            : sample,
                path              : "${params.outdir}/${sample}/${createNotebook.baseName}/output/${sampleZarr.name}",
                image_scale_factor: imageScaleFactor,
            ]
        }

        WRITE_SDATA_SAMPLESHEET(buildSamplesheetInput(
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
        def follicleNotebook = file(createRegistry.create_follicle_sdata.path)
        def cellIdsFile = file(params.cell_ids_file)
        def follicleInputs = sampleArtifacts.map { sample, sampleZarr, rowParams ->
            tuple(
                follicleNotebook.toString(),
                follicleNotebook.baseName,
                timerScript.toString(),
                sampleZarr.toString(),
                sample,
                "${params.outdir}/${sample}/${follicleNotebook.baseName}",
                sample,
                rowParams,
                cellIdsFile,
                createRegistry.create_follicle_sdata.params
            )
        }

        def follicleArtifacts = CREATE_FOLLICLE_SDATA(FOLLICLE_SDATA_QUARTO_PARAMS(follicleInputs).notebook_inputs).artifacts

        def follicleArtifactRows = follicleArtifacts.flatMap { sample, zarrPaths, _rowParams ->
            // Nextflow emits a single Path for one match and a List<Path> for many; normalize.
            def zarrs = zarrPaths instanceof List ? zarrPaths : [zarrPaths]
            def imageScaleFactor = _rowParams.image_scale_factor ?: 1.0
            zarrs.collect { zarr ->
                [
                    sample            : sample,
                    cell              : zarr.baseName,
                    path              : "${params.outdir}/${sample}/${follicleNotebook.baseName}/output/${zarr.name}",
                    image_scale_factor: imageScaleFactor,
                ]
            }
        }

        WRITE_FOLLICLE_SAMPLESHEET(buildSamplesheetInput(
            follicleArtifactRows,
            'follicle_sdata_samplesheet.csv',
            "${params.outdir}/${follicleNotebook.baseName}",
        ))
    }
}
