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
    // Hardlink (not copy) into results: workDir and outdir share the scratch
    // filesystem, so linking avoids a second full copy of the large zarr.
    publishDir { clusterSdataGpuOocPublishDir(sample) },
        mode: 'link',
        saveAs: { fn -> fn.endsWith('.samplesheet_row.csv') ? null : fn }

    input:
    tuple val(sample), path(input_path)
    val chunk_size
    val n_top_genes
    val resolutions
    val min_counts
    val min_cells

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
    // The script takes a space-separated nargs list; the param is comma-separated
    // so it can be given as a single --resolutions value on the Nextflow CLI.
    // toString() first: a single value (--resolutions 1.0) arrives as a Number.
    if (resolutions) clusterArgs << "--resolutions ${resolutions.toString().tokenize(',').join(' ')}"
    // Omitted when unset so the filtering cut, like the resolution sweep, keeps its
    // single definition in the Python script rather than being restated in the config.
    if (min_counts != '') clusterArgs << "--min_counts ${min_counts}"
    if (min_cells != '') clusterArgs << "--min_cells ${min_cells}"
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
