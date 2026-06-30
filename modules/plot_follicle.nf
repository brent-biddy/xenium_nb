include { paramsFile } from './quarto_params'

process PLOT_FOLLICLE {
    tag "${follicle_id}"

    publishDir { "${params.outdir}/${sample}/plot_follicle" },
        mode: 'copy',
        saveAs: { fn -> fn.startsWith('output/') ? fn : "${follicle_id}_${fn}" }

    input:
    tuple val(follicle_id), val(sample), path(input_path), path('params.yml')
    path notebook
    path 'timer.py'

    output:
    path "plot_follicle.*", emit: reports
    path "output/**", optional: true, hidden: true, emit: output_tree

    script:
    """
    # Redirect cache and temp dirs into the writable work dir so quarto/deno
    # don't try to write to a read-only /tmp on HPC compute nodes.
    export XDG_CACHE_HOME="\$PWD/.cache"
    export TMPDIR="\$PWD/tmp"
    mkdir -p "\$XDG_CACHE_HOME" "\$TMPDIR"

    quarto render ${notebook} --execute-params params.yml --output-dir .
    """

    stub:
    """
    touch plot_follicle.pptx
    touch plot_follicle.timing.tsv
    """
}

workflow {
    if (!params.samplesheet) error "Please provide --samplesheet"

    def plotFollicleNotebook = file("${projectDir}/notebooks/analyze/plot_follicle.qmd")
    def timerScript          = file("${projectDir}/bin/timer.py")

    Channel
        .fromPath(params.samplesheet)
        .splitCsv(header: true)
        .map { row ->
            if (!row.sample) error "Samplesheet row missing 'sample': ${row}"
            if (!row.path)   error "Samplesheet row missing 'path': ${row}"
            def follicleId = "${row.sample}_${row.cell}"
            tuple(follicleId, row.sample, file(row.path), paramsFile(follicleId, plotFollicleNotebook, row))
        }
        .set { plotInputs } // tuple(follicle_id, sample, staged_path, params_yml)

    PLOT_FOLLICLE(plotInputs, plotFollicleNotebook, timerScript)
}
