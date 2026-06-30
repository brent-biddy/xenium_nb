process DOWNSAMPLE_SDATA {
    tag "${sample}"

    publishDir { "${params.outdir}/${sample}/downsample_sdata" },
        mode: 'copy'

    input:
    tuple val(sample), path(input_path)

    output:
    tuple val(sample), path("downsampled.zarr"), emit: zarr
    path "downsample_sdata_timing.tsv", emit: timing

    script:
    def downsampleArgs = ["--sample ${sample}", "--path ${input_path}"]
    if (params.fraction) downsampleArgs << "--fraction ${params.fraction}"
    if (params.n_cells)  downsampleArgs << "--n_cells ${params.n_cells}"
    """
    export XDG_CACHE_HOME="\$PWD/.cache"
    export TMPDIR="\$PWD/tmp"
    mkdir -p "\$XDG_CACHE_HOME" "\$TMPDIR"

    downsample_sdata.py ${downsampleArgs.join(' ')}
    """

    stub:
    """
    mkdir -p downsampled.zarr
    touch downsampled.zarr/.zgroup
    touch downsample_sdata_timing.tsv
    """
}

workflow {
    if (!params.samplesheet) error "Please provide --samplesheet"
    if (!params.fraction && !params.n_cells) error "Please provide --fraction or --n_cells"

    Channel
        .fromPath(params.samplesheet)
        .splitCsv(header: true)      // Map(sample, path)
        .map { row ->
            if (!row.sample) error "Samplesheet row missing 'sample': ${row}"
            if (!row.path)   error "Samplesheet row missing 'path': ${row}"
            tuple(row.sample, file(row.path))
        }                            // tuple(sample, path)
        | DOWNSAMPLE_SDATA
}
