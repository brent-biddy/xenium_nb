// Processes for create.nf. CREATE_SDATA runs a plain Python script; CREATE_FOLLICLE_SDATA
// still renders a Quarto notebook so the per-cell spatial plots are published as HTML.

process CREATE_SDATA {
    tag "${sample}"

    publishDir { "${params.outdir}/${sample}/create_sdata" },
        mode: 'copy'

    input:
    tuple val(sample), path(input_path), val(he_image), val(he_alignment)

    output:
    tuple val(sample), path('output/*.zarr'), emit: artifacts
    path "output/**", optional: true, hidden: true, emit: output_tree

    script:
    def heImageArg = he_image    ? "--he_image ${he_image}"       : ""
    def heAlignArg = he_alignment ? "--he_alignment ${he_alignment}" : ""
    """
    export XDG_CACHE_HOME="\$PWD/.cache"
    export TMPDIR="\$PWD/tmp"
    mkdir -p "\$XDG_CACHE_HOME" "\$TMPDIR"

    create_sdata.py \
        --sample ${sample} \
        --path ${input_path} \
        --n_jobs ${task.cpus} \
        ${heImageArg} \
        ${heAlignArg}
    """

    stub:
    """
    mkdir -p output/${sample}.zarr
    touch output/${sample}.zarr/.zgroup
    touch output/${sample}.zarr/.zattrs
    touch output/${sample}.zarr/.zmetadata
    """
}

process CREATE_FOLLICLE_SDATA {
    tag "${sample}"

    publishDir { "${params.outdir}/${sample}/follicle_sdata" },
        mode: 'copy',
        saveAs: { fn -> fn.startsWith('output/') ? fn : "${sample}_${fn}" }

    input:
    tuple val(sample), path(input_path), path('params.yml')
    path cell_ids_file
    path notebook
    path 'timer.py'

    output:
    tuple val(sample), path('output/*.zarr'), emit: artifacts
    path "follicle_sdata.*", emit: reports
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
    touch follicle_sdata.html
    touch follicle_sdata.timing.tsv
    """
}
