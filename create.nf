#!/usr/bin/env nextflow

// Builds Xenium SpatialData artifacts and the samplesheets that downstream
// workflows (analyze.nf) consume. Modes:
//   downsample      - crop raw Xenium output to a bounding box region
//   sdata           - run create_sdata.py per ROI from raw Xenium output
//   follicle_sdata  - run create_follicle_sdata.py per cell ID from existing sample zarrs
//   all             - run sdata + follicle_sdata, chaining outputs
//   concat          - concatenate existing SpatialData zarrs into one (requires --output_name)

nextflow.enable.dsl = 2

include { DOWNSAMPLE_XENIUM_REGION; CREATE_SDATA; CREATE_FOLLICLE_SDATA; CONCAT_SDATA } from './modules/create_notebooks'

workflow {
    if (!params.samplesheet) error "Please provide --samplesheet"
    if (!params.create)      error "Please provide --create (sdata, follicle_sdata, all)"

    def createMode = params.create.toLowerCase()

    if (!(createMode in ['downsample', 'sdata', 'follicle_sdata', 'all', 'concat'])) {
        error "Invalid --create '${createMode}'. Valid values are: downsample, sdata, follicle_sdata, all, concat"
    }

    def follicleSourceArtifacts = null
    def cellIdsFile = file(params.cell_ids_file)

    Channel
        .fromPath(params.samplesheet)
        .splitCsv(header: true)
        .map { row ->
            if (!row.sample) error "Samplesheet row missing 'sample': ${row}"
            if (!row.path)   error "Samplesheet row missing 'path': ${row}"
            tuple(row.sample, file(row.path), row)
        }
        .set { sampleRowsList } // tuple(sample, staged_path, row_map)

    // ---- downsample: crop raw Xenium output to a bounding box region ----
    if (createMode == 'downsample') {

        sampleRowsList
            .map { sample, stagedPath, rowMap ->
                def heImage = rowMap.he_image     ? new File(rowMap.he_image     as String).absolutePath : ""
                def heAlign = rowMap.he_alignment ? new File(rowMap.he_alignment as String).absolutePath : ""
                def regionName = rowMap.region_name ?: sample
                tuple(sample, stagedPath, rowMap.xmin, rowMap.ymin, rowMap.xmax, rowMap.ymax, regionName, heImage, heAlign)
            }
            .set { downsampleInputs } // tuple(sample, staged_path, xmin, ymin, xmax, ymax, region_name, he_image, he_alignment)

        DOWNSAMPLE_XENIUM_REGION(downsampleInputs)
    }

    // ---- sdata: raw Xenium -> per-sample SpatialData zarr ----
    if (createMode == 'sdata' || createMode == 'all') {

        sampleRowsList
            .map { sample, stagedPath, rowMap ->
                def heImage = rowMap.he_image    ? new File(rowMap.he_image    as String).absolutePath : ""
                def heAlign = rowMap.he_alignment ? new File(rowMap.he_alignment as String).absolutePath : ""
                tuple(sample, stagedPath, heImage, heAlign)
            }
            .set { createSdataInputs } // tuple(sample, staged_path, he_image, he_alignment)

        CREATE_SDATA(createSdataInputs) | set { createSdataRun }
        // createSdataRun.artifacts: tuple(sample, zarr)

        follicleSourceArtifacts = createSdataRun.artifacts
            .join(sampleRowsList) // tuple(sample, zarr, staged_path, row_map)
            .map { sample, zarr, inputPath, rowParams ->
                tuple(sample, zarr)
            }
        // follicleSourceArtifacts: tuple(sample, zarr)

        follicleSourceArtifacts
            .map { sample, sampleZarr ->
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
            .map { sample, stagedPath, rowMap -> tuple(sample, stagedPath) }
        // follicleSourceArtifacts: tuple(sample, staged_path)
    }

    // ---- follicle_sdata: per-sample SpatialData -> per-cell-ID subset zarrs ----
    if (createMode == 'follicle_sdata' || createMode == 'all') {

        CREATE_FOLLICLE_SDATA(follicleSourceArtifacts, cellIdsFile, params.radius) | set { follicleRun }
        // follicleRun.artifacts: tuple(sample, List<zarr>)

        follicleRun.artifacts
            .flatMap { sample, zarrPaths ->
                // Nextflow emits a single Path for one match and a List<Path> for many; normalize.
                def zarrs = zarrPaths instanceof List ? zarrPaths : [zarrPaths]
                zarrs.collect { zarr ->
                    [
                        sample: sample,
                        cell  : zarr.baseName,
                        path  : "${params.outdir}/${sample}/follicle_sdata/output/${zarr.name}",
                    ]
                }
            }
            .set { follicleArtifactRows } // Map(sample, cell, path) per item

        Channel.of('sample,cell,path')
            .concat(follicleArtifactRows.map { row -> "${row.sample},${row.cell},${row.path}" })
            .collectFile(name: 'follicle_sdata_samplesheet.csv', newLine: true, storeDir: params.outdir, sort: false)
    }

    // ---- concat: merge all input SpatialData zarrs into a single zarr ----
    if (createMode == 'concat') {

        sampleRowsList
            .map { sample, stagedPath, rowMap -> tuple(sample, stagedPath) }
            .collect()
            .map { pairs -> tuple(pairs.collect { it[0] }.join('_'), pairs.collect { it[1] }) }
            .set { concatInputs } // tuple(sample_id, List<Path>)

        CONCAT_SDATA(concatInputs)
    }
}
