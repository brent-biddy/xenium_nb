#!/usr/bin/env nextflow

import groovy.json.JsonOutput

nextflow.enable.dsl = 2

include { CREATE_SPATIALDATA } from './modules/create_spatialdata'
include { SUBSET_FOLLICLE } from './modules/subset_follicle'
include { WRITE_SAMPLESHEET as WRITE_SAMPLE_ANALYSIS_INPUTS } from './modules/write_samplesheet'
include { WRITE_SAMPLESHEET as WRITE_FOLLICLE_ANALYSIS_INPUTS } from './modules/write_samplesheet'

workflow {
    if (!params.samplesheet) error "Please provide --samplesheet"
    def expectedColumns = ['sample', 'path'] as Set
    def cellIdsRegistry = params.cell_ids_registry ?: [:]
    def cellIdsFileValue = params.cell_ids_file?.toString()
    def cellIdsFilePath = file(cellIdsRegistry.get(cellIdsFileValue, cellIdsFileValue))

    def sampleRows = Channel
        .fromPath(params.samplesheet)
        .splitCsv(header: true)
        .map { row ->
            def columns = row.keySet().findAll { it != null && it != '' } as Set
            if (columns != expectedColumns) {
                error "Build samplesheet must contain exactly these columns: sample,path. Found: ${columns.join(',')}"
            }
            if (!row.sample) error "Samplesheet row missing 'sample': ${row}"
            if (!row.path)   error "Samplesheet row missing 'path': ${row}"
            def rowMap = new LinkedHashMap(row)
            def sample = rowMap.sample.toString()
            tuple(sample, file(rowMap.path), rowMap)
        }

    def timerScript = file("${projectDir}/bin/timer.py")
    def createNotebook = file(params.producer_registry.create_spatialdata.path)
    def subsetNotebook = file(params.producer_registry.subset_follicle.path)
    def runSubsetFollicle = params.run_subset_follicle as boolean

    def createInputs = sampleRows.map { sample, inputPath, rowParams ->
        def publishDir = "${params.outdir}/${sample}/${createNotebook.baseName}"
        tuple(createNotebook, timerScript, inputPath, sample, publishDir, rowParams)
    }

    def createSpatialdata = CREATE_SPATIALDATA(createInputs)
    def sampleArtifacts = createSpatialdata.artifacts

    def sampleArtifactRows = sampleArtifacts
        .map { sample, sampleZarr, rowParams ->
            [sample, "${params.outdir}/${sample}/${createNotebook.baseName}/output/${sample}.zarr"]
        }
        .collect()
        .map { rows -> JsonOutput.toJson(rows) }
        .map { rowsJson -> tuple('sample_analysis_inputs.csv', rowsJson) }

    WRITE_SAMPLE_ANALYSIS_INPUTS(sampleArtifactRows)

    if (runSubsetFollicle) {
        def subsetInputs = sampleArtifacts.map { sample, sampleZarr, rowParams ->
            def publishDir = "${params.outdir}/${sample}/${subsetNotebook.baseName}"
            tuple(sample, sampleZarr, rowParams, cellIdsFilePath, subsetNotebook, timerScript, publishDir)
        }

        def subsetFollicle = SUBSET_FOLLICLE(subsetInputs)
        def follicleArtifactGroups = subsetFollicle.artifacts

        def follicleArtifactRows = follicleArtifactGroups
            .flatMap { sample, outputDir ->
                def zarrs = outputDir.toFile().listFiles()?.findAll { it.name.endsWith('.zarr') } ?: []
                zarrs.collect { zarr ->
                    def cellId = zarr.baseName
                    ["${sample}_${cellId}", "${params.outdir}/${sample}/${subsetNotebook.baseName}/output/${cellId}.zarr"]
                }
            }
            .collect()
            .map { rows -> JsonOutput.toJson(rows) }
            .map { rowsJson -> tuple('follicle_analysis_inputs.csv', rowsJson) }

        WRITE_FOLLICLE_ANALYSIS_INPUTS(follicleArtifactRows)
    }
}
