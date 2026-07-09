// Published output directory for the merged artifact. Single-sourced here so the
// publishDir directive and the handoff samplesheet row can't drift apart.
def concatSdataPublishDir() {
    "${params.outdir}/concat_sdata"
}

process CONCAT_SDATA {

    // saveAs drops the row fragment from the published dir; it is only needed on
    // the channel for main.nf to collectFile into the aggregate sheet.
    publishDir { concatSdataPublishDir() },
        mode: 'copy',
        saveAs: { fn -> fn.endsWith('.samplesheet_row.csv') ? null : fn }

    input:
    path input_paths

    output:
    path '*.zarr', emit: artifacts
    path "concat_sdata_timing.tsv", emit: timing
    path "concat_sdata_session_info.txt", emit: session_info
    // One `sample,path` line for the merged zarr (sample = its basename). The name
    // is derived from the merged sample keys inside concat_sdata.py, so it is only
    // known at runtime — identify it as the *.zarr that is not a staged input.
    path "concat_sdata.samplesheet_row.csv", emit: samplesheet_row

    script:
    """
    export XDG_CACHE_HOME="\$PWD/.cache"
    export TMPDIR="\$PWD/tmp"
    mkdir -p "\$XDG_CACHE_HOME" "\$TMPDIR"

    concat_sdata.py --paths ${input_paths}

    # The staged input zarrs also match *.zarr; the merged output is the one that
    # is not among the inputs.
    merged=""
    for z in \$(ls -d *.zarr); do
        case " ${input_paths} " in
            *" \$z "*) ;;
            *) merged="\$z" ;;
        esac
    done
    printf '%s,%s' "\${merged%.zarr}" "${concatSdataPublishDir()}/\$merged" > concat_sdata.samplesheet_row.csv
    """

    stub:
    """
    mkdir -p merged.zarr
    touch merged.zarr/.zgroup
    touch merged.zarr/.zattrs
    touch merged.zarr/.zmetadata
    touch concat_sdata_timing.tsv
    touch concat_sdata_session_info.txt

    merged=""
    for z in \$(ls -d *.zarr); do
        case " ${input_paths} " in
            *" \$z "*) ;;
            *) merged="\$z" ;;
        esac
    done
    printf '%s,%s' "\${merged%.zarr}" "${concatSdataPublishDir()}/\$merged" > concat_sdata.samplesheet_row.csv
    """
}
