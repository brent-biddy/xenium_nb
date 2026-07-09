process DOWNSAMPLE_XENIUM_REGION {
    tag "${sample}/${region_name}"

    // Hardlink (not copy) into results: workDir and outdir share the scratch
    // filesystem, so linking avoids a second full copy of the cropped region.
    publishDir { "${params.outdir}/${sample}/downsample_xenium_region" },
        mode: 'link'

    input:
    tuple val(sample), path(input_path), val(xmin), val(ymin), val(xmax), val(ymax), val(region_name), val(he_image), val(he_alignment)

    output:
    tuple val(sample), path("${region_name}/*"), emit: artifacts

    script:
    def downsampleArgs = [
        "${input_path}",
        "--bbox ${xmin} ${ymin} ${xmax} ${ymax}",
        "--region_name ${region_name}",
        "--output_dir .",
        "--threads ${task.cpus}",
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
    mkdir -p ${region_name}
    touch ${region_name}/experiment.xenium
    """
}
