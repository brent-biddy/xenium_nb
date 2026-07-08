// Published output directory for this step's per-sample artifacts. Single-sourced
// here so the publishDir directive, the emitted zarr, and the handoff samplesheet
// row all reference the same location and cannot drift apart.
def clusterSdataPublishDir(sample) {
    "${params.outdir}/${sample}/cluster_sdata"
}

process CLUSTER_SDATA {
    tag "${sample}"

    publishDir { clusterSdataPublishDir(sample) },
        mode: 'copy'

    input:
    tuple val(sample), path(input_path)

    output:
    tuple val(sample), path("clustered.zarr"), emit: zarr
    path "cluster_sdata_timing.tsv", emit: timing
    // One `sample,path` line pointing at the published zarr; main.nf collectFiles
    // these into a ready-to-use handoff samplesheet. No trailing newline — the
    // collectFile(newLine: true) call adds the separator.
    path "${sample}.samplesheet_row.csv", emit: samplesheet_row

    script:
    """
    export XDG_CACHE_HOME="\$PWD/.cache"
    export TMPDIR="\$PWD/tmp"
    mkdir -p "\$XDG_CACHE_HOME" "\$TMPDIR"

    cluster_sdata.py --sample ${sample} --path ${input_path}

    printf '%s' '${sample},${clusterSdataPublishDir(sample)}/clustered.zarr' > ${sample}.samplesheet_row.csv
    """

    stub:
    """
    mkdir -p clustered.zarr
    touch clustered.zarr/.zgroup
    touch cluster_sdata_timing.tsv

    printf '%s' '${sample},${clusterSdataPublishDir(sample)}/clustered.zarr' > ${sample}.samplesheet_row.csv
    """
}
