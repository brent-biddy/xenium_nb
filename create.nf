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
    def createRegistry = new groovy.json.JsonSlurper()
        .parse(new File("${projectDir}/assets/notebook_registry.json"))
        .create

    Channel
        .fromPath(params.samplesheet)
        .splitCsv(header: true)
        .map { row ->
            if (!row.sample) error "Samplesheet row missing 'sample': ${row}"
            if (!row.path)   error "Samplesheet row missing 'path': ${row}"
            tuple(row.sample.toString(), file(row.path), row)
        }
        .collect(flat: false)
        .set { sampleRowsList } // List<tuple(sample, staged_path, row_map)>

    // ---- sdata: raw Xenium -> per-sample SpatialData zarr ----
    if (createMode == 'sdata' || createMode == 'all') {

        sampleRowsList
            .flatMap { rows ->
                rows.collect { row ->
                    tuple(row[0], row[2], createRegistry.create_sdata.params)
                }
            }
            .set { sdataParamsInputs } // tuple(sample, row_map, declared_params)

        SDATA_PARAMS(sdataParamsInputs) | set { createSdataParams } // tuple(sample, params_yml)

        sampleRowsList
            .flatMap { rows ->
                rows.collect { row -> tuple(row[0], row[1]) }
            }
            .join(createSdataParams.params_file)
            .set { createSdataInputs } // tuple(sample, staged_path, params_yml)

        CREATE_SDATA(createSdataInputs, createNotebook, timerScript) | set { createSdataRun }
        // createSdataRun.artifacts: tuple(sample, zarr)

        follicleSourceArtifacts = createSdataRun.artifacts
            .join(
                sampleRowsList.flatMap { rows -> rows }
            )
            .map { sample, zarr, inputPath, rowParams ->
                tuple(sample, zarr, rowParams)
            }
        // follicleSourceArtifacts: tuple(sample, zarr, row_map)

        follicleSourceArtifacts
            .map { sample, sampleZarr, rowParams ->
                [
                    sample: sample,
                    path  : "${params.outdir}/${sample}/create_sdata/output/${sampleZarr.name}",
                ]
            }
            .set { sampleArtifactRows } // Map(sample, path) per item

        Channel.of('sample,path')
            .concat(sampleArtifactRows.map { row -> "${row.sample},${row.path}" })
            .collectFile(name: 'sample_sdata_samplesheet.csv', newLine: true, storeDir: params.outdir, sort: false)
    }

    // Skip CREATE_SDATA: caller's samplesheet already points at existing sample zarrs.
    if (createMode == 'follicle_sdata') {
        follicleSourceArtifacts = sampleRowsList.flatMap { rows -> rows }
        // follicleSourceArtifacts: tuple(sample, staged_path, row_map)
    }

    // ---- follicle_sdata: per-sample SpatialData -> per-cell-ID subset zarrs ----
    if (createMode == 'follicle_sdata' || createMode == 'all') {

        follicleSourceArtifacts
            .collect(flat: false)
            .set { follicleSourceArtifactRows } // List<tuple(sample, zarr, row_map)>

        follicleSourceArtifactRows
            .flatMap { rows ->
                rows.collect { row ->
                    // Override path with the zarr so the notebook receives the correct input path.
                    def rowMap = row[2] + [path: row[1].toString()]
                    tuple(row[0], rowMap, createRegistry.create_follicle_sdata.params)
                }
            }
            .set { follicleParamsInputs } // tuple(sample, row_map, declared_params)

        FOLLICLE_SDATA_PARAMS(follicleParamsInputs) | set { follicleParams } // tuple(sample, params_yml)

        follicleSourceArtifactRows
            .flatMap { rows ->
                rows.collect { row -> tuple(row[0], row[1]) }
            }
            .join(follicleParams.params_file)
            .set { follicleInputs } // tuple(sample, zarr, params_yml)

        CREATE_FOLLICLE_SDATA(follicleInputs, cellIdsFile, follicleNotebook, timerScript) | set { follicleRun }
        // follicleRun.artifacts: tuple(sample, List<zarr>)

        follicleSourceArtifactRows
            .flatMap { rows ->
                rows.collect { row -> tuple(row[0], row[2] + [sample: row[0]]) }
            }
            .set { follicleRowParams } // tuple(sample, row_map)

        follicleRun.artifacts
            .join(follicleRowParams)
            .flatMap { sample, zarrPaths, rowParams ->
                // Nextflow emits a single Path for one match and a List<Path> for many; normalize.
                def zarrs = zarrPaths instanceof List ? zarrPaths : [zarrPaths]
                zarrs.collect { zarr ->
                    [
                        sample: sample,
                        cell  : zarr.baseName,
                        path  : "${params.outdir}/${sample}/create_follicle_sdata/output/${zarr.name}",
                    ]
                }
            }
            .set { follicleArtifactRows } // Map(sample, cell, path) per item

        Channel.of('sample,cell,path')
            .concat(follicleArtifactRows.map { row -> "${row.sample},${row.cell},${row.path}" })
            .collectFile(name: 'follicle_sdata_samplesheet.csv', newLine: true, storeDir: params.outdir, sort: false)
    }
}
