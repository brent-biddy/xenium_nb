// Published output directory for this step's per-sample artifacts. Single-sourced
// here so the publishDir directive, the emitted h5ad, and the handoff samplesheet
// row all reference the same location and cannot drift apart.
def createAdataPublishDir(sample) {
    "${params.outdir}/${sample}/create_adata"
}

process CREATE_ADATA {
    tag "${sample}"

    // saveAs drops the per-sample row fragment from the published dir; it is only
    // needed on the channel for main.nf to collectFile into the aggregate sheet.
    // Hardlink (not copy) into results: workDir and outdir share the scratch
    // filesystem, so linking avoids a second full copy of the matrix.
    publishDir { createAdataPublishDir(sample) },
        mode: 'link',
        saveAs: { fn -> fn.endsWith('.samplesheet_row.csv') ? null : fn }

    input:
    tuple val(sample), path(input_path)

    output:
    tuple val(sample), path("${sample}.h5ad"), emit: artifacts
    path "${sample}_timing.tsv", emit: timing
    path "${sample}_session_info.txt", emit: session_info
    // One `sample,path` line pointing at the published h5ad; main.nf collectFiles
    // these into a ready-to-use handoff samplesheet. No trailing newline — the
    // collectFile(newLine: true) call adds the separator.
    path "${sample}.samplesheet_row.csv", emit: samplesheet_row

    script:
    """
    export XDG_CACHE_HOME="\$PWD/.cache"
    export TMPDIR="\$PWD/tmp"
    mkdir -p "\$XDG_CACHE_HOME" "\$TMPDIR"

    create_adata.py --sample ${sample} --path ${input_path}

    printf '%s' '${sample},${createAdataPublishDir(sample)}/${sample}.h5ad' > ${sample}.samplesheet_row.csv
    """

    stub:
    """
    touch ${sample}.h5ad
    touch ${sample}_timing.tsv
    touch ${sample}_session_info.txt

    printf '%s' '${sample},${createAdataPublishDir(sample)}/${sample}.h5ad' > ${sample}.samplesheet_row.csv
    """
}
