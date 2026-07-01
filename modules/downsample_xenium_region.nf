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
    mkdir -p output/${region_name}
    touch output/${region_name}/experiment.xenium
    """
}
