// Published output directory for this step's per-sample artifacts. Single-sourced
// here so the publishDir directive, the emitted zarr, and the handoff samplesheet
// row all reference the same location and cannot drift apart.
def downsampleSdataPublishDir(sample) {
    "${params.outdir}/${sample}/downsample_sdata"
}

process DOWNSAMPLE_SDATA {
    tag "${sample}"

    // saveAs drops the per-sample row fragment from the published dir; it is only
    // needed on the channel for main.nf to collectFile into the aggregate sheet.
    // Hardlink (not copy) into results: workDir and outdir share the scratch
    // filesystem, so linking avoids a second full copy of the large zarr.
    publishDir { downsampleSdataPublishDir(sample) },
        mode: 'link',
        saveAs: { fn -> fn.endsWith('.samplesheet_row.csv') ? null : fn }

    input:
    tuple val(sample), path(input_path)

    output:
    tuple val(sample), path("downsampled.zarr"), emit: zarr
    path "downsample_sdata_timing.tsv", emit: timing
    // One `sample,path` line pointing at the published zarr; main.nf collectFiles
    // these into a ready-to-use handoff samplesheet. No trailing newline — the
    // collectFile(newLine: true) call adds the separator.
    path "${sample}.samplesheet_row.csv", emit: samplesheet_row

    script:
    def downsampleArgs = ["--sample ${sample}", "--path ${input_path}"]
    // Mutually exclusive downstream, and `fraction` carries a non-null default, so
    // an explicit --n_cells has to suppress it or the script gets handed both.
    if (params.n_cells) downsampleArgs << "--n_cells ${params.n_cells}"
    else if (params.fraction) downsampleArgs << "--fraction ${params.fraction}"
    if (params.stratify_by) downsampleArgs << "--stratify_by ${params.stratify_by}"
    // Compared against null, not truthiness: Groovy treats 0 as false, and 0 is a
    // meaningful value for both (floor disabled / seed zero).
    if (params.min_cells_per_cluster != null) downsampleArgs << "--min_cells_per_cluster ${params.min_cells_per_cluster}"
    if (params.seed != null) downsampleArgs << "--seed ${params.seed}"
    """
    export XDG_CACHE_HOME="\$PWD/.cache"
    export TMPDIR="\$PWD/tmp"
    mkdir -p "\$XDG_CACHE_HOME" "\$TMPDIR"

    downsample_sdata.py ${downsampleArgs.join(' ')}

    printf '%s' '${sample},${downsampleSdataPublishDir(sample)}/downsampled.zarr' > ${sample}.samplesheet_row.csv
    """

    stub:
    """
    mkdir -p downsampled.zarr
    touch downsampled.zarr/.zgroup
    touch downsample_sdata_timing.tsv

    printf '%s' '${sample},${downsampleSdataPublishDir(sample)}/downsampled.zarr' > ${sample}.samplesheet_row.csv
    """
}
