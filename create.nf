#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { CREATE_SPATIALDATA } from './modules/create_spatialdata'
include { SUBSET_FOLLICLE } from './modules/subset_follicle'
include { WRITE_SAMPLESHEET as WRITE_SAMPLE_ANALYSIS_INPUTS } from './modules/write_samplesheet'
include { WRITE_SAMPLESHEET as WRITE_FOLLICLE_ANALYSIS_INPUTS } from './modules/write_samplesheet'

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

workflow {
    if (!params.samplesheet) error "Please provide --samplesheet"
    def producerRegistry = NotebookRegistry.producer(projectDir.toString())
    def cellIdsFilePath = file(params.cell_ids_file)
    def createMode = params.create.toLowerCase()
    if (!(createMode in ['sdata', 'follicle_sdata', 'all'])) {
        error "Invalid create '${createMode}'. Valid values are: sdata, follicle_sdata, all"
    }

    def sampleArtifacts = null
    if (createMode in ['sdata', 'all']) {
        def sampleRows = parseSamplesheet(params.samplesheet, 'Create samplesheet')

        def timerScript = file("${projectDir}/bin/timer.py")
        def createNotebook = file(producerRegistry.create_sdata.path)
        def createNotebookParams = producerRegistry.create_sdata.params
        def createInputs = sampleRows.map { sample, inputPath, rowParams ->
            def publishDir = "${params.outdir}/${sample}/${createNotebook.baseName}"
            tuple(createNotebook, timerScript, inputPath, sample, publishDir, rowParams, createNotebookParams)
        }

        def createSpatialdata = CREATE_SPATIALDATA(createInputs)
        sampleArtifacts = createSpatialdata.artifacts

        def sampleArtifactRows = sampleArtifacts
            .map { sample, sampleZarr, rowParams ->
                [
                    sample: sample,
                    path  : "${params.outdir}/${sample}/${createNotebook.baseName}/output/${sample}.zarr",
                ]
            }
            .collect()
            .map { rows -> groovy.json.JsonOutput.toJson(rows) }
            .map { rowsJson ->
                tuple(
                    'sample_sdata_samplesheet.csv',
                    rowsJson,
                    "${params.outdir}/${createNotebook.baseName}"
                )
            }

        WRITE_SAMPLE_ANALYSIS_INPUTS(sampleArtifactRows)
    }

    if (createMode == 'follicle_sdata') {
        sampleArtifacts = parseSamplesheet(params.samplesheet, 'Create follicle samplesheet')
    }

    if (createMode in ['follicle_sdata', 'all']) {
        def timerScript = file("${projectDir}/bin/timer.py")
        def subsetNotebook = file(producerRegistry.create_follicle_sdata.path)
        def subsetNotebookParams = producerRegistry.create_follicle_sdata.params
        def subsetInputs = sampleArtifacts.map { sample, sampleZarr, rowParams ->
            def publishDir = "${params.outdir}/${sample}/${subsetNotebook.baseName}"
            tuple(sample, sampleZarr, rowParams, cellIdsFilePath, subsetNotebook, timerScript, publishDir, subsetNotebookParams)
        }

        def subsetFollicle = SUBSET_FOLLICLE(subsetInputs)
        def follicleArtifacts = subsetFollicle.artifacts

        def follicleArtifactRows = follicleArtifacts
            .flatMap { sample, zarrPaths ->
                def zarrs = zarrPaths instanceof List ? zarrPaths : [zarrPaths]
                zarrs.collect { zarr ->
                    def cellId = zarr.baseName
                    [
                        sample: sample,
                        cell  : cellId,
                        path  : "${params.outdir}/${sample}/${subsetNotebook.baseName}/output/${cellId}.zarr",
                    ]
                }
            }
            .collect()
            .map { rows -> groovy.json.JsonOutput.toJson(rows) }
            .map { rowsJson ->
                tuple(
                    'follicle_sdata_samplesheet.csv',
                    rowsJson,
                    "${params.outdir}/${subsetNotebook.baseName}"
                )
            }

        WRITE_FOLLICLE_ANALYSIS_INPUTS(follicleArtifactRows)
    }
}
