process CLUSTER_SDATA_GPU_OOC {
    tag "${sample}"

    // --nv passes through the host NVIDIA driver and CUDA libs into the container.
    containerOptions '--nv'

    publishDir { "${params.outdir}/${sample}/cluster_sdata_gpu_ooc" },
        mode: 'copy'

    input:
    tuple val(sample), path(input_path)
    val chunk_size
    val n_top_genes

    output:
    tuple val(sample), path("clustered.zarr"), emit: zarr
    path "cluster_sdata_gpu_ooc_timing.tsv", emit: timing

    script:
    def clusterArgs = ["--sample ${sample}", "--path ${input_path}"]
    if (chunk_size)  clusterArgs << "--chunk-size ${chunk_size}"
    if (n_top_genes) clusterArgs << "--n-top-genes ${n_top_genes}"
    """
    export XDG_CACHE_HOME="\$PWD/.cache"
    export TMPDIR="\$PWD/tmp"
    mkdir -p "\$XDG_CACHE_HOME" "\$TMPDIR"

    cluster_sdata_gpu_ooc.py ${clusterArgs.join(' ')}
    """

    stub:
    """
    mkdir -p clustered.zarr
    touch clustered.zarr/.zgroup
    touch cluster_sdata_gpu_ooc_timing.tsv
    """
}
