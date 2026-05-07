#!/usr/bin/env nextflow

// Builds Xenium SpatialData artifacts and the samplesheets that downstream
// workflows (analyze.nf) consume. Modes:
//   sdata           - render create_sdata.qmd per ROI from raw Xenium output
//   follicle_sdata  - render create_follicle_sdata.qmd per cell ID from existing sample zarrs
//   all             - run both, chaining sdata outputs into the follicle step

nextflow.enable.dsl = 2

include { paramsFile } from './modules/quarto_params'
include { CREATE_SDATA; CREATE_FOLLICLE_SDATA } from './modules/create_notebooks'

workflow {
    if (!params.samplesheet) error "Please provide --samplesheet"
    if (!params.create)      error "Please provide --create (sdata, follicle_sdata, all)"

    def createMode = params.create.toLowerCase()

    if (!(createMode in ['sdata', 'follicle_sdata', 'all'])) {
        error "Invalid --create '${createMode}'. Valid values are: sdata, follicle_sdata, all"
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
            tuple(row.sample, file(row.path), row)
        }
        .set { sampleRowsList } // tuple(sample, staged_path, row_map)

    // ---- sdata: raw Xenium -> per-sample SpatialData zarr ----
    if (createMode == 'sdata' || createMode == 'all') {

        sampleRowsList
            .map { sample, stagedPath, rowMap ->
                tuple(sample, stagedPath, paramsFile(sample, 'create_sdata', createRegistry.create_sdata.params, rowMap, params.outdir))
            }
            .set { createSdataInputs } // tuple(sample, staged_path, params_yml)

        CREATE_SDATA(createSdataInputs, createNotebook, timerScript) | set { createSdataRun }
        // createSdataRun.artifacts: tuple(sample, zarr)

        follicleSourceArtifacts = createSdataRun.artifacts
            .join(sampleRowsList) // tuple(sample, zarr, staged_path, row_map)
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
        follicleSourceArtifacts = sampleRowsList
        // follicleSourceArtifacts: tuple(sample, staged_path, row_map)
    }

    // ---- follicle_sdata: per-sample SpatialData -> per-cell-ID subset zarrs ----
    if (createMode == 'follicle_sdata' || createMode == 'all') {

        follicleSourceArtifacts
            .map { sample, stagedPath, rowMap ->
                // Override path with the zarr so the notebook receives the correct input path.
                tuple(sample, stagedPath, paramsFile(sample, 'create_follicle_sdata', createRegistry.create_follicle_sdata.params, rowMap + [path: stagedPath.toString()], params.outdir))
            }
            .set { follicleInputs } // tuple(sample, staged_path, params_yml)

        CREATE_FOLLICLE_SDATA(follicleInputs, cellIdsFile, follicleNotebook, timerScript) | set { follicleRun }
        // follicleRun.artifacts: tuple(sample, List<zarr>)

        follicleRun.artifacts
            .flatMap { sample, zarrPaths ->
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
