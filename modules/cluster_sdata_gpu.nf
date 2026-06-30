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

workflow {
    if (!params.samplesheet) error "Please provide --samplesheet"

    Channel
        .fromPath(params.samplesheet)
        .splitCsv(header: true)      // Map(sample, path)
        .map { row ->
            if (!row.sample) error "Samplesheet row missing 'sample': ${row}"
            if (!row.path)   error "Samplesheet row missing 'path': ${row}"
            tuple(row.sample, file(row.path))
        }                            // tuple(sample, path)
        | CLUSTER_SDATA_GPU
}
