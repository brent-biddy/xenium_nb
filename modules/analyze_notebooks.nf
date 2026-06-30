// Quarto notebook processes for analyze.nf. Each process renders one notebook
// against a pre-built SpatialData artifact and publishes the resulting reports.

process CLUSTER_SDATA {
    tag "${sample}"

    publishDir { "${params.outdir}/${sample}/cluster_sdata" },
        mode: 'copy'

    input:
    tuple val(sample), path(input_path)

    output:
    tuple val(sample), path("clustered.zarr"), emit: zarr
    path "cluster_sdata_timing.tsv", emit: timing

    script:
    """
    export XDG_CACHE_HOME="\$PWD/.cache"
    export TMPDIR="\$PWD/tmp"
    mkdir -p "\$XDG_CACHE_HOME" "\$TMPDIR"

    cluster_sdata.py --sample ${sample} --path ${input_path}
    """

    stub:
    """
    mkdir -p clustered.zarr
    touch clustered.zarr/.zgroup
    touch cluster_sdata_timing.tsv
    """
}

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
