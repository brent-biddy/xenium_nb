process CLUSTER_SDATA_GPU {
    tag "${sample}"

    // --nv passes through the host NVIDIA driver and CUDA libs into the container.
    containerOptions '--nv'

    publishDir { "${params.outdir}/${sample}/cluster_sdata_gpu" },
        mode: 'copy'

    input:
    tuple val(sample), path(input_path)

    output:
    tuple val(sample), path("clustered.zarr"), emit: zarr
    path "cluster_sdata_gpu_timing.tsv", emit: timing

    script:
    """
    export XDG_CACHE_HOME="\$PWD/.cache"
    export TMPDIR="\$PWD/tmp"
    mkdir -p "\$XDG_CACHE_HOME" "\$TMPDIR"

    cluster_sdata_gpu.py --sample ${sample} --path ${input_path}
    """

    stub:
    """
    mkdir -p clustered.zarr
    touch clustered.zarr/.zgroup
    touch cluster_sdata_gpu_timing.tsv
    """
}
