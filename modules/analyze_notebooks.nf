// Quarto notebook processes for analyze.nf. Each process renders one notebook
// against a pre-built SpatialData artifact and publishes the resulting reports.

process PLOT_FOLLICLE {
    tag "${sample}:${cell}"

    publishDir { "${params.outdir}/${sample}/plot_follicle" },
        mode: 'copy',
        saveAs: { fn -> fn.startsWith('output/') ? fn : "${sample_id}_${fn}" }

    input:
    tuple val(sample_id), val(sample), val(cell), path(input_path), val(params_b64)
    path notebook
    path 'timer.py'

    output:
    path "plot_follicle.*", emit: reports
    path "output/**", optional: true, hidden: true, emit: output_tree

    script:
    """
    python3 -c "import base64; open('params.yml', 'w').write(base64.b64decode('${params_b64}').decode())"

    # Redirect cache and temp dirs into the writable work dir so quarto/deno
    # don't try to write to a read-only /tmp on HPC compute nodes.
    export XDG_CACHE_HOME="\$PWD/.cache"
    export TMPDIR="\$PWD/tmp"
    mkdir -p "\$XDG_CACHE_HOME" "\$TMPDIR"

    quarto render ${notebook} --execute-params params.yml --output-dir .
    """

    stub:
    """
    touch params.yml
    touch plot_follicle.pptx
    touch plot_follicle.timing.tsv
    """
}
