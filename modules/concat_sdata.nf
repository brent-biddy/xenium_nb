process CONCAT_SDATA {

    publishDir { "${params.outdir}/concat_sdata" },
        mode: 'copy'

    input:
    path input_paths

    output:
    path 'merged.zarr', emit: artifacts
    path "concat_sdata_timing.tsv", emit: timing
    path "concat_sdata_session_info.txt", emit: session_info

    script:
    """
    export XDG_CACHE_HOME="\$PWD/.cache"
    export TMPDIR="\$PWD/tmp"
    mkdir -p "\$XDG_CACHE_HOME" "\$TMPDIR"

    concat_sdata.py --paths ${input_paths}
    """

    stub:
    """
    mkdir -p merged.zarr
    touch merged.zarr/.zgroup
    touch merged.zarr/.zattrs
    touch merged.zarr/.zmetadata
    touch concat_sdata_timing.tsv
    touch concat_sdata_session_info.txt
    """
}

workflow {
    if (!params.samplesheet) error "Please provide --samplesheet"

    Channel
        .fromPath(params.samplesheet)
        .splitCsv(header: true)  // Map(path, ...)
        .map { row ->
            if (!row.path) error "Samplesheet row missing 'path': ${row}"
            file(row.path)
        }                        // path
        .collect()               // List<path>
        | CONCAT_SDATA
}
