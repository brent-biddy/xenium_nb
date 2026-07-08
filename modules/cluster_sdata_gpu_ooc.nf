// Published output directory for this step's per-sample artifacts. Single-sourced
// here so the publishDir directive, the emitted zarr, and the handoff samplesheet
// row all reference the same location and cannot drift apart.
def clusterSdataGpuOocPublishDir(sample) {
    "${params.outdir}/${sample}/cluster_sdata_gpu_ooc"
}

process CLUSTER_SDATA_GPU_OOC {
    tag "${sample}"

    // --nv passes through the host NVIDIA driver and CUDA libs into the container.
    containerOptions '--nv'

    // saveAs drops the per-sample row fragment from the published dir; it is only
    // needed on the channel for main.nf to collectFile into the aggregate sheet.
    publishDir { clusterSdataGpuOocPublishDir(sample) },
        mode: 'copy',
        saveAs: { fn -> fn.endsWith('.samplesheet_row.csv') ? null : fn }

    input:
    tuple val(sample), path(input_path)
    val chunk_size
    val n_top_genes

    output:
    tuple val(sample), path("clustered.zarr"), emit: zarr
    path "cluster_sdata_gpu_ooc_timing.tsv", emit: timing
    // One `sample,path` line pointing at the published zarr; main.nf collectFiles
    // these into a ready-to-use handoff samplesheet. No trailing newline — the
    // collectFile(newLine: true) call adds the separator.
    path "${sample}.samplesheet_row.csv", emit: samplesheet_row

    script:
    def clusterArgs = ["--sample ${sample}", "--path ${input_path}"]
    if (chunk_size)  clusterArgs << "--chunk-size ${chunk_size}"
    if (n_top_genes) clusterArgs << "--n-top-genes ${n_top_genes}"
    """
    export XDG_CACHE_HOME="\$PWD/.cache"
    export TMPDIR="\$PWD/tmp"
    mkdir -p "\$XDG_CACHE_HOME" "\$TMPDIR"

    cluster_sdata_gpu_ooc.py ${clusterArgs.join(' ')}

    printf '%s' '${sample},${clusterSdataGpuOocPublishDir(sample)}/clustered.zarr' > ${sample}.samplesheet_row.csv
    """

    stub:
    """
    mkdir -p clustered.zarr
    touch clustered.zarr/.zgroup
    touch cluster_sdata_gpu_ooc_timing.tsv

    printf '%s' '${sample},${clusterSdataGpuOocPublishDir(sample)}/clustered.zarr' > ${sample}.samplesheet_row.csv
    """
}
