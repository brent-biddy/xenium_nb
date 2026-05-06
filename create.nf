#!/usr/bin/env nextflow

// Builds Xenium SpatialData artifacts and the samplesheets that downstream
// workflows (analyze.nf) consume. Modes:
//   sdata           - render create_sdata.qmd per ROI from raw Xenium output
//   follicle_sdata  - render create_follicle_sdata.qmd per cell ID from existing sample zarrs
//   all             - run both, chaining sdata outputs into the follicle step

nextflow.enable.dsl = 2

include { WRITE_QUARTO_PARAMS as WRITE_CREATE_SDATA_PARAMS } from './modules/write_quarto_params'
include { WRITE_QUARTO_PARAMS as WRITE_CREATE_FOLLICLE_SDATA_PARAMS } from './modules/write_quarto_params'
include { CREATE_SDATA; CREATE_FOLLICLE_SDATA } from './modules/create_notebooks'
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
    def timerScript = file("${projectDir}/bin/timer.py")
    def createMode = params.create.toLowerCase()
    if (!(createMode in ['sdata', 'follicle_sdata', 'all'])) {
        error "Invalid create '${createMode}'. Valid values are: sdata, follicle_sdata, all"
    }

    def createSdataArtifacts = null
    def follicleSourceArtifacts = null
    def cellIdsFile = file(params.cell_ids_file)
    def createNotebook = file("${projectDir}/notebooks/create_sdata.qmd")
    def follicleNotebook = file("${projectDir}/notebooks/create_follicle_sdata.qmd")

    // ---- sdata: raw Xenium -> per-sample SpatialData zarr ----
    if (createMode in ['sdata', 'all']) {
        def sampleRows = parseSamplesheet(params.samplesheet, 'Create samplesheet')

        def createParamsInputs = sampleRows.map { sample, inputPath, rowParams ->
            tuple(sample, inputPath, rowParams + [sample: sample], ['sample', 'path', 'n_jobs'])
        }
        def createRunInputs = sampleRows.map { sample, inputPath, rowParams ->
            tuple(sample, inputPath)
        }
        def createRowParams = sampleRows.map { sample, inputPath, rowParams ->
            tuple(sample, rowParams + [sample: sample])
        }

        def createSdataParams = WRITE_CREATE_SDATA_PARAMS(createParamsInputs)
        def createSdataRun = CREATE_SDATA(
            createRunInputs.join(createSdataParams.params_file),
            Channel.value(createNotebook),
            Channel.value(timerScript),
        )
        createSdataArtifacts = createSdataRun.artifacts
        follicleSourceArtifacts = createSdataArtifacts.join(createRowParams)

        def sampleArtifactRows = follicleSourceArtifacts.map { sample, sampleZarr, rowParams ->
            def imageScaleFactor = rowParams.image_scale_factor ?: 1.0
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
        follicleSourceArtifacts = parseSamplesheet(params.samplesheet, 'Create follicle samplesheet')
    }

    // ---- follicle_sdata: per-sample SpatialData -> per-cell-ID subset zarrs ----
    if (createMode in ['follicle_sdata', 'all']) {
        def follicleParamsInputs = follicleSourceArtifacts.map { sample, sampleZarr, rowParams ->
            tuple(sample, sampleZarr, rowParams + [sample: sample], ['sample', 'path', 'cell_ids_file', 'radius', 'image_scale_factor'])
        }
        def follicleRunInputs = follicleSourceArtifacts.map { sample, sampleZarr, rowParams ->
            tuple(sample, sampleZarr)
        }
        def follicleRowParams = follicleSourceArtifacts.map { sample, sampleZarr, rowParams ->
            tuple(sample, rowParams + [sample: sample])
        }

        def follicleParams = WRITE_CREATE_FOLLICLE_SDATA_PARAMS(follicleParamsInputs)
        def follicleRun = CREATE_FOLLICLE_SDATA(
            follicleRunInputs.join(follicleParams.params_file),
            Channel.value(cellIdsFile),
            Channel.value(follicleNotebook),
            Channel.value(timerScript),
        )
        def follicleArtifacts = follicleRun.artifacts

        def follicleArtifactRows = follicleArtifacts.join(follicleRowParams).flatMap { sample, zarrPaths, rowParams ->
            // Nextflow emits a single Path for one match and a List<Path> for many; normalize.
            def zarrs = zarrPaths instanceof List ? zarrPaths : [zarrPaths]
            def imageScaleFactor = rowParams.image_scale_factor ?: 1.0
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
