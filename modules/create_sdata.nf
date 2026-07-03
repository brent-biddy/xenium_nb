process CREATE_SDATA {
    tag "${sample}"

    publishDir { "${params.outdir}/${sample}/create_sdata" },
        mode: 'copy'

    input:
    tuple val(sample), path(input_path), val(he_image), val(he_alignment)

    output:
    tuple val(sample), path("${sample}.zarr"), emit: artifacts
    path "${sample}_timing.tsv", emit: timing
    path "${sample}_session_info.txt", emit: session_info

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
    mkdir -p ${sample}.zarr
    touch ${sample}.zarr/.zgroup
    touch ${sample}.zarr/.zattrs
    touch ${sample}.zarr/.zmetadata
    touch ${sample}_timing.tsv
    touch ${sample}_session_info.txt
    """
}
