// Published output directory for this step's per-sample follicle zarrs. Single-
// sourced here so the publishDir directive and the handoff samplesheet rows can't
// drift apart.
def follicleSdataPublishDir(sample) {
    "${params.outdir}/${sample}/follicle_sdata"
}

process CREATE_FOLLICLE_SDATA {
    tag "${sample}"

    // saveAs drops the row fragment from the published dir; it is only needed on
    // the channel for main.nf to collectFile into the aggregate sheet.
    // Hardlink (not copy) into results: workDir and outdir share the scratch
    // filesystem, so linking avoids a second full copy of the large zarrs.
    publishDir { follicleSdataPublishDir(sample) },
        mode: 'link',
        saveAs: { fn -> fn.endsWith('.samplesheet_row.csv') ? null : fn }

    input:
    tuple val(sample), path(input_path)
    path cell_ids_file
    val radius

    output:
    tuple val(sample), path('*.zarr'), emit: artifacts
    path "${sample}_timing.tsv", emit: timing
    path "${sample}_session_info.txt", emit: session_info
    // One `sample,cell,path` line per produced per-cell zarr (cell = its basename),
    // for the plot_follicle step. main.nf collectFiles these into the aggregate
    // sheet. The staged input zarr also matches *.zarr, so exclude it.
    path "${sample}.samplesheet_row.csv", emit: samplesheet_row

    script:
    def follicleArgs = [
        "--sample ${sample}",
        "--path ${input_path}",
        "--cell_ids_file ${cell_ids_file}",
        "--radius ${radius}",
    ]
    """
    export XDG_CACHE_HOME="\$PWD/.cache"
    export TMPDIR="\$PWD/tmp"
    mkdir -p "\$XDG_CACHE_HOME" "\$TMPDIR"

    create_follicle_sdata.py ${follicleArgs.join(' ')}

    in_name=\$(basename ${input_path})
    first=1
    for z in \$(ls -d *.zarr | sort); do
        [ "\$z" = "\$in_name" ] && continue
        if [ \$first -eq 1 ]; then first=0; else printf '\\n'; fi
        printf '%s,%s,%s' "${sample}" "\${z%.zarr}" "${follicleSdataPublishDir(sample)}/\$z"
    done > ${sample}.samplesheet_row.csv
    """

    stub:
    """
    mkdir -p ${sample}.zarr
    touch ${sample}.zarr/.zgroup
    touch ${sample}.zarr/.zattrs
    touch ${sample}.zarr/.zmetadata
    touch ${sample}_timing.tsv
    touch ${sample}_session_info.txt

    in_name=\$(basename ${input_path})
    first=1
    for z in \$(ls -d *.zarr | sort); do
        [ "\$z" = "\$in_name" ] && continue
        if [ \$first -eq 1 ]; then first=0; else printf '\\n'; fi
        printf '%s,%s,%s' "${sample}" "\${z%.zarr}" "${follicleSdataPublishDir(sample)}/\$z"
    done > ${sample}.samplesheet_row.csv
    """
}
