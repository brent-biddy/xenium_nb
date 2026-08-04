// Cell type annotation: score each clustered zarr's cells against a marker YAML and
// label every Leiden cluster with its best-matching cell type.

// Published output directory for this step's per-sample artifacts. Single-sourced
// here so the publishDir directive, the emitted zarr, and the handoff samplesheet
// row all reference the same location and cannot drift apart.
def annotateSdataPublishDir(sample) {
    "${params.outdir}/${sample}/annotate_sdata"
}

process ANNOTATE_SDATA {
    tag "${sample}"

    // saveAs drops the per-sample row fragment from the published dir; it is only
    // needed on the channel for main.nf to collectFile into the aggregate sheet.
    // Hardlink (not copy) into results: workDir and outdir share the scratch
    // filesystem, so linking avoids a second full copy of the large zarr.
    publishDir { annotateSdataPublishDir(sample) },
        mode: 'link',
        saveAs: { fn -> fn.endsWith('.samplesheet_row.csv') ? null : fn }

    input:
    tuple val(sample), path(input_path)
    // The marker YAML, staged into the work dir rather than read from its repo
    // path: the task runs inside a container that need not have the launch dir
    // bound, and on OSCER the compute node may not see it at all.
    path markers
    val resolutions
    val exclude_nonspecific

    output:
    tuple val(sample), path("annotated.zarr"), emit: zarr
    // Cluster x cell type score table. Published because it is how a surprising
    // label gets traced back to the numbers that produced it.
    path "annotate_sdata_scores.tsv", emit: scores
    path "annotate_sdata_timing.tsv", emit: timing
    // One `sample,path` line pointing at the published zarr; main.nf collectFiles
    // these into a ready-to-use handoff samplesheet. No trailing newline — the
    // collectFile(newLine: true) call adds the separator.
    path "${sample}.samplesheet_row.csv", emit: samplesheet_row

    script:
    def annotateArgs = ["--sample ${sample}", "--path ${input_path}", "--markers ${markers}"]
    // The script takes a space-separated nargs list; the param is comma-separated
    // so it can be given as a single --resolutions value on the Nextflow CLI.
    // toString() first: a single value (--resolutions 1.0) arrives as a Number.
    if (resolutions) annotateArgs << "--resolutions ${resolutions.toString().tokenize(',').join(' ')}"
    if (exclude_nonspecific) annotateArgs << "--exclude_nonspecific"
    """
    export XDG_CACHE_HOME="\$PWD/.cache"
    export TMPDIR="\$PWD/tmp"
    mkdir -p "\$XDG_CACHE_HOME" "\$TMPDIR"

    annotate_sdata.py ${annotateArgs.join(' ')}

    printf '%s' '${sample},${annotateSdataPublishDir(sample)}/annotated.zarr' > ${sample}.samplesheet_row.csv
    """

    stub:
    """
    mkdir -p annotated.zarr
    touch annotated.zarr/.zgroup
    touch annotate_sdata_scores.tsv
    touch annotate_sdata_timing.tsv

    printf '%s' '${sample},${annotateSdataPublishDir(sample)}/annotated.zarr' > ${sample}.samplesheet_row.csv
    """
}
