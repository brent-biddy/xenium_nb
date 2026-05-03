#!/usr/bin/env nextflow

import groovy.json.JsonOutput

nextflow.enable.dsl = 2

include { RUN_NOTEBOOK } from './modules/run_notebook'

workflow {
    // Validate required params
    if (!params.samplesheet)                              error "Please provide --samplesheet"
    if (!params.notebooks || params.notebooks.size() < 1) error "Please provide --notebooks"

    // Parse samplesheet — all columns are forwarded as notebook params.
    // Pipeline-level params (cell_ids_file, radius, n_jobs) are merged in
    // by the process only when a notebook declares them.
    //
    // Required columns: sample_id, data_path
    // Optional columns: roi_id (used to group outputs); falls back to sample_id
    def samples = Channel
        .fromPath(params.samplesheet)
        .splitCsv(header: true)
        .map { row ->
            if (!row.sample_id) error "Samplesheet row missing 'sample_id': ${row}"
            if (!row.data_path) error "Samplesheet row missing 'data_path': ${row}"
            def roi_id = row.roi_id ?: row.sample_id
            tuple(row.sample_id, roi_id, JsonOutput.toJson(row))
        }

    // Resolve notebook paths from config
    def notebooks = Channel.fromList(params.notebooks).map { file(it) }

    // Stage timer.py into each process work directory so notebooks can
    // import it without any path manipulation
    def timer_script = file("${projectDir}/bin/timer.py")

    // Cross-product: every sample x every notebook runs as a separate process.
    // publish_dir groups outputs by roi_id so follicle-level runs land back
    // under their parent ROI directory.
    samples
        .combine(notebooks)
        .map { sample_id, roi_id, json, nb ->
            def publish_dir = "${params.outdir}/${roi_id}/${nb.baseName}"
            def output_name = "${sample_id}_${nb.baseName}.html"
            tuple(nb, timer_script, sample_id, publish_dir, output_name, json)
        }
        | RUN_NOTEBOOK
}
