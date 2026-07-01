process DOWNSAMPLE_SDATA {
    tag "${sample}"

    publishDir { "${params.outdir}/${sample}/downsample_sdata" },
        mode: 'copy'

    input:
    tuple val(sample), path(input_path)

    output:
    tuple val(sample), path("downsampled.zarr"), emit: zarr
    path "downsample_sdata_timing.tsv", emit: timing

    script:
    def downsampleArgs = ["--sample ${sample}", "--path ${input_path}"]
    if (params.fraction) downsampleArgs << "--fraction ${params.fraction}"
    if (params.n_cells)  downsampleArgs << "--n_cells ${params.n_cells}"
    """
    export XDG_CACHE_HOME="\$PWD/.cache"
    export TMPDIR="\$PWD/tmp"
    mkdir -p "\$XDG_CACHE_HOME" "\$TMPDIR"

    downsample_sdata.py ${downsampleArgs.join(' ')}
    """

    stub:
    """
    mkdir -p downsampled.zarr
    touch downsampled.zarr/.zgroup
    touch downsample_sdata_timing.tsv
    """
}
