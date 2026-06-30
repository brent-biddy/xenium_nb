process CREATE_FOLLICLE_SDATA {
    tag "${sample}"

    publishDir { "${params.outdir}/${sample}/follicle_sdata" },
        mode: 'copy'

    input:
    tuple val(sample), path(input_path)
    path cell_ids_file
    val radius

    output:
    tuple val(sample), path('output/*.zarr'), emit: artifacts
    path "output/**", optional: true, hidden: true, emit: output_tree

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
    """

    stub:
    """
    mkdir -p output/${sample}.zarr
    touch output/${sample}.zarr/.zgroup
    touch output/${sample}.zarr/.zattrs
    touch output/${sample}.zarr/.zmetadata
    """
}

workflow {
    if (!params.samplesheet)   error "Please provide --samplesheet"
    if (!params.cell_ids_file) error "Please provide --cell_ids_file"

    def cellIdsFile = file(params.cell_ids_file)

    Channel
        .fromPath(params.samplesheet)
        .splitCsv(header: true)      // Map(sample, path)
        .map { row ->
            if (!row.sample) error "Samplesheet row missing 'sample': ${row}"
            if (!row.path)   error "Samplesheet row missing 'path': ${row}"
            tuple(row.sample, file(row.path))
        }
        .set { inputs }              // tuple(sample, path)

    CREATE_FOLLICLE_SDATA(inputs, cellIdsFile, params.radius)
}
