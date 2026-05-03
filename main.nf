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
    // Required columns: sample, path
    // Optional columns: roi_id (used to group outputs); falls back to sample
    def samples = Channel
        .fromPath(params.samplesheet)
        .splitCsv(header: true)
        .map { row ->
            if (!row.sample) error "Samplesheet row missing 'sample': ${row}"
            if (!row.path)   error "Samplesheet row missing 'path': ${row}"
            // The 'path' column is staged into each work dir. The notebook
            // sees the staged basename in its CWD, so rewrite the row before
            // serializing it to params.json.
            def input_path = file(row.path)
            row.path = input_path.getName()
            def roi_id = row.roi_id ?: row.sample
            tuple(row.sample, roi_id, input_path, JsonOutput.toJson(row))
        }

    // Resolve notebook paths from config. Accept either a single string or a list.
    def notebook_list = params.notebooks instanceof List ? params.notebooks : [params.notebooks]
    def notebooks = Channel.of(*notebook_list).map { file(it) }

    // Stage timer.py into each process work directory so notebooks can
    // import it without any path manipulation
    def timer_script = file("${projectDir}/bin/timer.py")

    // Cross-product: every sample x every notebook runs as a separate process.
    // publish_dir groups outputs by roi_id so follicle-level runs land back
    // under their parent ROI directory.
    samples
        .combine(notebooks)
        .map { sample, roi_id, input_path, json, nb ->
            def publish_dir = "${params.outdir}/${roi_id}/${nb.baseName}"
            def output_name = "${sample}_${nb.baseName}.html"
            tuple(nb, timer_script, input_path, sample, publish_dir, output_name, json)
        }
        | RUN_NOTEBOOK
}
