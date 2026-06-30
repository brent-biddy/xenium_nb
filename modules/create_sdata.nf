process CREATE_SDATA {
    tag "${sample}"

    publishDir { "${params.outdir}/${sample}/create_sdata" },
        mode: 'copy'

    input:
    tuple val(sample), path(input_path), val(he_image), val(he_alignment)

    output:
    tuple val(sample), path('output/*.zarr'), emit: artifacts
    path "output/**", optional: true, hidden: true, emit: output_tree

    script:
    def sdataArgs = ["--sample ${sample}", "--path ${input_path}", "--n_jobs ${task.cpus}"]
    if (he_image)     sdataArgs << "--he_image ${he_image}"
    if (he_alignment) sdataArgs << "--he_alignment ${he_alignment}"
    """
    export XDG_CACHE_HOME="\$PWD/.cache"
    export TMPDIR="\$PWD/tmp"
    mkdir -p "\$XDG_CACHE_HOME" "\$TMPDIR"

    create_sdata.py ${sdataArgs.join(' ')}
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
    if (!params.samplesheet) error "Please provide --samplesheet"

    Channel
        .fromPath(params.samplesheet)
        .splitCsv(header: true)        // Map(sample, path[, he_image, he_alignment])
        .map { row ->
            if (!row.sample) error "Samplesheet row missing 'sample': ${row}"
            if (!row.path)   error "Samplesheet row missing 'path': ${row}"
            def heImage = row.he_image     ? new File(row.he_image     as String).absolutePath : ""
            def heAlign = row.he_alignment ? new File(row.he_alignment as String).absolutePath : ""
            tuple(row.sample, file(row.path), heImage, heAlign)
        }                              // tuple(sample, path, he_image, he_alignment)
        | CREATE_SDATA
}
