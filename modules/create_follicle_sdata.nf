process CREATE_FOLLICLE_SDATA {
    tag "${sample}"

    publishDir { "${params.outdir}/${sample}/follicle_sdata" },
        mode: 'copy'

    input:
    tuple val(sample), path(input_path)
    path cell_ids_file
    val radius

    output:
    tuple val(sample), path('*.zarr'), emit: artifacts
    path "${sample}_timing.tsv", emit: timing
    path "${sample}_session_info.txt", emit: session_info

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
    mkdir -p ${sample}.zarr
    touch ${sample}.zarr/.zgroup
    touch ${sample}.zarr/.zattrs
    touch ${sample}.zarr/.zmetadata
    touch ${sample}_timing.tsv
    touch ${sample}_session_info.txt
    """
}
