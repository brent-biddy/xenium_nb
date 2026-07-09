// Published output directory for this step's per-sample artifacts. Single-sourced
// here so the publishDir directive, the emitted zarr, and the handoff samplesheet
// row all reference the same location and cannot drift apart.
def createSdataPublishDir(sample) {
    "${params.outdir}/${sample}/create_sdata"
}

process CREATE_SDATA {
    tag "${sample}"

    // saveAs drops the per-sample row fragment from the published dir; it is only
    // needed on the channel for main.nf to collectFile into the aggregate sheet.
    // Hardlink (not copy) into results: workDir and outdir share the scratch
    // filesystem, so linking avoids a second full copy of the large zarr.
    publishDir { createSdataPublishDir(sample) },
        mode: 'link',
        saveAs: { fn -> fn.endsWith('.samplesheet_row.csv') ? null : fn }

    input:
    tuple val(sample), path(input_path), val(he_image), val(he_alignment)

    output:
    tuple val(sample), path("${sample}.zarr"), emit: artifacts
    path "${sample}_timing.tsv", emit: timing
    path "${sample}_session_info.txt", emit: session_info
    // One `sample,path` line pointing at the published zarr; main.nf collectFiles
    // these into a ready-to-use handoff samplesheet. No trailing newline — the
    // collectFile(newLine: true) call adds the separator.
    path "${sample}.samplesheet_row.csv", emit: samplesheet_row

    script:
    def sdataArgs = ["--sample ${sample}", "--path ${input_path}", "--n_jobs ${task.cpus}"]
    if (he_image)     sdataArgs << "--he_image ${he_image}"
    if (he_alignment) sdataArgs << "--he_alignment ${he_alignment}"
    """
    export XDG_CACHE_HOME="\$PWD/.cache"
    export TMPDIR="\$PWD/tmp"
    mkdir -p "\$XDG_CACHE_HOME" "\$TMPDIR"

    create_sdata.py ${sdataArgs.join(' ')}

    printf '%s' '${sample},${createSdataPublishDir(sample)}/${sample}.zarr' > ${sample}.samplesheet_row.csv
    """

    stub:
    """
    mkdir -p ${sample}.zarr
    touch ${sample}.zarr/.zgroup
    touch ${sample}.zarr/.zattrs
    touch ${sample}.zarr/.zmetadata
    touch ${sample}_timing.tsv
    touch ${sample}_session_info.txt

    printf '%s' '${sample},${createSdataPublishDir(sample)}/${sample}.zarr' > ${sample}.samplesheet_row.csv
    """
}
