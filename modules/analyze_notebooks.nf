// Quarto notebook processes for analyze.nf. Each process renders one notebook
// against a pre-built SpatialData artifact and publishes the resulting reports.

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
