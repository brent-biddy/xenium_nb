// Processes for create.nf. Both CREATE_SDATA and CREATE_FOLLICLE_SDATA run
// plain Python scripts — no Quarto/Deno overhead for pure ETL steps.

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
    def sdataArgs = ["--sample ${sample}", "--path ${input_path}", "--n_jobs ${task.cpus}"]
    if (he_image)     sdataArgs << "--he_image ${he_image}"
    if (he_alignment) sdataArgs << "--he_alignment ${he_alignment}"
    """
    export XDG_CACHE_HOME="\$PWD/.cache"
    export TMPDIR="\$PWD/tmp"
    mkdir -p "\$XDG_CACHE_HOME" "\$TMPDIR"

    create_sdata.py ${sdataArgs.join(' ')}
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
        mode: 'copy'

    input:
    tuple val(sample), path(input_path)
    path cell_ids_file
    val radius

    output:
    tuple val(sample), path('output/*.zarr'), emit: artifacts
    path "output/**", optional: true, hidden: true, emit: output_tree

    script:
    def follicleArgs = [
        "--sample ${sample}",
        "--path ${input_path}",
        "--cell_ids_file ${cell_ids_file}",
        "--radius ${radius}",
    ]
    """
    export XDG_CACHE_HOME="\$PWD/.cache"
    export TMPDIR="\$PWD/tmp"
    mkdir -p "\$XDG_CACHE_HOME" "\$TMPDIR"

    create_follicle_sdata.py ${follicleArgs.join(' ')}
    """

    stub:
    """
    mkdir -p output/${sample}.zarr
    touch output/${sample}.zarr/.zgroup
    touch output/${sample}.zarr/.zattrs
    touch output/${sample}.zarr/.zmetadata
    """
}
