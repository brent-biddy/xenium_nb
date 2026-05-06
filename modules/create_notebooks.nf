// Quarto notebook processes for create.nf. Each process renders one notebook
// to build a SpatialData artifact and publishes the result.

process CREATE_SDATA {
    tag "${sample}"

    publishDir { "${params.outdir}/${sample}/create_sdata" }, mode: 'copy'

    input:
    tuple val(sample), path(input_path), path('params.yml')
    path notebook
    path 'timer.py'

    output:
    tuple val(sample), path('output/*.zarr'), emit: artifacts
    path "create_sdata.*", emit: reports
    path "output/**", optional: true, hidden: true, emit: output_tree

    script:
    """
    # Redirect cache and temp dirs into the writable work dir so quarto/deno
    # don't try to write to a read-only /tmp on HPC compute nodes.
    export XDG_CACHE_HOME="\$PWD/.cache"
    export TMPDIR="\$PWD/tmp"
    mkdir -p "\$XDG_CACHE_HOME" "\$TMPDIR"

    quarto render ${notebook} --execute-params params.yml -P n_jobs:${task.cpus} --output-dir .
    """

    stub:
    """
    mkdir -p output/${sample}.zarr
    touch output/${sample}.zarr/.zgroup
    touch output/${sample}.zarr/.zattrs
    touch output/${sample}.zarr/.zmetadata
    touch create_sdata.html
    touch create_sdata.timing.tsv
    """
}

process CREATE_FOLLICLE_SDATA {
    tag "${sample}"

    publishDir { "${params.outdir}/${sample}/create_follicle_sdata" }, mode: 'copy'

    input:
    tuple val(sample), path(input_path), path('params.yml')
    path cell_ids_file
    path notebook
    path 'timer.py'

    output:
    tuple val(sample), path('output/*.zarr'), emit: artifacts
    path "create_follicle_sdata.*", emit: reports
    path "output/**", optional: true, hidden: true, emit: output_tree

    script:
    """
    # Redirect cache and temp dirs into the writable work dir so quarto/deno
    # don't try to write to a read-only /tmp on HPC compute nodes.
    export XDG_CACHE_HOME="\$PWD/.cache"
    export TMPDIR="\$PWD/tmp"
    mkdir -p "\$XDG_CACHE_HOME" "\$TMPDIR"

    quarto render ${notebook} --execute-params params.yml -P n_jobs:${task.cpus} --output-dir .
    """

    stub:
    """
    mkdir -p output/${sample}.zarr
    touch output/${sample}.zarr/.zgroup
    touch output/${sample}.zarr/.zattrs
    touch output/${sample}.zarr/.zmetadata
    touch create_follicle_sdata.html
    touch create_follicle_sdata.timing.tsv
    """
}
