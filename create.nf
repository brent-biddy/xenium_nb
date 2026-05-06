#!/usr/bin/env nextflow

// Builds Xenium SpatialData artifacts and the samplesheets that downstream
// workflows (analyze.nf) consume. Modes:
//   sdata           - render create_sdata.qmd per ROI from raw Xenium output
//   follicle_sdata  - render create_follicle_sdata.qmd per cell ID from existing sample zarrs
//   all             - run both, chaining sdata outputs into the follicle step

nextflow.enable.dsl = 2

include { WRITE_QUARTO_PARAMS as SDATA_PARAMS } from './modules/write_quarto_params'
include { WRITE_QUARTO_PARAMS as FOLLICLE_SDATA_PARAMS } from './modules/write_quarto_params'
include { CREATE_SDATA; CREATE_FOLLICLE_SDATA } from './modules/create_notebooks'
include { WRITE_SAMPLESHEET as SDATA_SAMPLESHEET } from './modules/write_samplesheet'
include { WRITE_SAMPLESHEET as FOLLICLE_SAMPLESHEET } from './modules/write_samplesheet'

workflow {
    def createMode = params.create.toLowerCase()

    if (!params.samplesheet) error "Please provide --samplesheet"
    if (!(createMode in ['sdata', 'follicle_sdata', 'all'])) {
        error "Invalid create '${createMode}'. Valid values are: sdata, follicle_sdata, all"
    }

    def follicleSourceArtifacts = null
    def timerScript = file("${projectDir}/bin/timer.py")
    def cellIdsFile = file(params.cell_ids_file)
    def createNotebook = file("${projectDir}/notebooks/create_sdata.qmd")
    def follicleNotebook = file("${projectDir}/notebooks/create_follicle_sdata.qmd")

    Channel
        .fromPath(params.samplesheet)
        .splitCsv(header: true)
        .map { row ->
            if (!row.sample) error "Samplesheet row missing 'sample': ${row}"
            if (!row.path)   error "Samplesheet row missing 'path': ${row}"
            tuple(row.sample.toString(), file(row.path), row)
        }
        .collect(flat: false)
        .set { sampleRowsList }

    // ---- sdata: raw Xenium -> per-sample SpatialData zarr ----
    if (createMode == 'sdata' || createMode == 'all') {

        sampleRowsList
            .flatMap { rows ->
                rows.collect { row ->
                    tuple(row[0], row[1], row[2], ['sample'])
                }
            }
            .set { sdataParamsInputs }

        SDATA_PARAMS(sdataParamsInputs) | set { createSdataParams }

        sampleRowsList
            .flatMap { rows ->
                rows.collect { row -> tuple(row[0], row[1]) }
            }
            .join(createSdataParams.params_file)
            .set { createSdataInputs }
            
        CREATE_SDATA(createSdataInputs, createNotebook, timerScript) | set { createSdataRun }

        follicleSourceArtifacts = createSdataRun.artifacts
            .join(
                sampleRowsList.flatMap { rows -> rows }
            )
            .map { sample, zarr, inputPath, rowParams ->
                tuple(sample, zarr, rowParams)
            }

        follicleSourceArtifacts
            .map { sample, sampleZarr, rowParams ->
                def imageScaleFactor = rowParams.image_scale_factor ?: 1.0
                [
                    sample            : sample,
                    path              : "${params.outdir}/${sample}/create_sdata/output/${sampleZarr.name}",
                    image_scale_factor: imageScaleFactor,
                ]
            }
            .set { sampleArtifactRows }
            
        sampleArtifactRows
            .collect()
            .map { rows -> tuple('sample_sdata_samplesheet.csv', rows) }
            .set { sdataSamplesheetInput }

        SDATA_SAMPLESHEET(sdataSamplesheetInput)
    }

    // Skip CREATE_SDATA: caller's samplesheet already points at existing sample zarrs.
    if (createMode == 'follicle_sdata') {
        follicleSourceArtifacts = sampleRowsList.flatMap { rows -> rows }
    }

    // ---- follicle_sdata: per-sample SpatialData -> per-cell-ID subset zarrs ----
    if (createMode == 'follicle_sdata' || createMode == 'all') {

        follicleSourceArtifacts
            .collect(flat: false)
            .set { follicleSourceArtifactRows }

        follicleSourceArtifactRows
            .flatMap { rows ->
                rows.collect { row ->
                    tuple(row[0], row[1], row[2], ['sample', 'cell_ids_file', 'radius', 'image_scale_factor'])
                }
            }
            .set { follicleParamsInputs }

        FOLLICLE_SDATA_PARAMS(follicleParamsInputs) | set { follicleParams }

        follicleSourceArtifactRows
            .flatMap { rows ->
                rows.collect { row -> tuple(row[0], row[1]) }
            }
            .join(follicleParams.params_file)
            .set { follicleInputs }

        CREATE_FOLLICLE_SDATA(follicleInputs, cellIdsFile, follicleNotebook, timerScript) | set { follicleRun }

        follicleSourceArtifactRows
            .flatMap { rows ->
                rows.collect { row -> tuple(row[0], row[2] + [sample: row[0]]) }
            }
            .set { follicleRowParams }

        follicleRun.artifacts
            .join(follicleRowParams)
            .flatMap { sample, zarrPaths, rowParams ->
                // Nextflow emits a single Path for one match and a List<Path> for many; normalize.
                def zarrs = zarrPaths instanceof List ? zarrPaths : [zarrPaths]
                def imageScaleFactor = rowParams.image_scale_factor ?: 1.0
                zarrs.collect { zarr ->
                    [
                        sample            : sample,
                        cell              : zarr.baseName,
                        path              : "${params.outdir}/${sample}/create_follicle_sdata/output/${zarr.name}",
                        image_scale_factor: imageScaleFactor,
                    ]
                }
            }
            .set { follicleArtifactRows }

        follicleArtifactRows
            .collect()
            .map { rows -> tuple('follicle_sdata_samplesheet.csv', rows) }
            .set { follicleSheetInput }

        FOLLICLE_SAMPLESHEET(follicleSheetInput)
    }
}
