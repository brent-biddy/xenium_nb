process DOWNSAMPLE_XENIUM_REGION {
    tag "${sample}/${region_name}"

    publishDir { "${params.outdir}/${sample}/downsample_xenium_region" },
        mode: 'copy',
        saveAs: { it.replaceFirst('output/', '') }

    input:
    tuple val(sample), path(input_path), val(xmin), val(ymin), val(xmax), val(ymax), val(region_name), val(he_image), val(he_alignment)

    output:
    tuple val(sample), path('output/*'), emit: artifacts

    script:
    def downsampleArgs = [
        "${input_path}",
        "--bbox ${xmin} ${ymin} ${xmax} ${ymax}",
        "--region_name ${region_name}",
        "--output_dir output",
    ]
    if (he_image)     downsampleArgs << "--he_image ${he_image}"
    if (he_alignment) downsampleArgs << "--he_alignment ${he_alignment}"
    """
    export XDG_CACHE_HOME="\$PWD/.cache"
    export TMPDIR="\$PWD/tmp"
    mkdir -p "\$XDG_CACHE_HOME" "\$TMPDIR"

    downsample_xenium_region.py ${downsampleArgs.join(' ')}
    """

    stub:
    """
    mkdir -p output/${region_name}
    touch output/${region_name}/experiment.xenium
    """
}

workflow {
    if (!params.samplesheet) error "Please provide --samplesheet"

    Channel
        .fromPath(params.samplesheet)
        .splitCsv(header: true)        // Map(sample, path, xmin, ymin, xmax, ymax[, region_name, he_image, he_alignment])
        .map { row ->
            if (!row.sample) error "Samplesheet row missing 'sample': ${row}"
            if (!row.path)   error "Samplesheet row missing 'path': ${row}"
            def heImage = row.he_image     ? new File(row.he_image     as String).absolutePath : ""
            def heAlign = row.he_alignment ? new File(row.he_alignment as String).absolutePath : ""
            def regionName = row.region_name ?: row.sample
            tuple(row.sample, file(row.path), row.xmin, row.ymin, row.xmax, row.ymax, regionName, heImage, heAlign)
        }                              // tuple(sample, path, xmin, ymin, xmax, ymax, region_name, he_image, he_alignment)
        | DOWNSAMPLE_XENIUM_REGION
}
