// Processes for create.nf. All processes run plain Python scripts —
// no Quarto/Deno overhead for pure ETL steps.

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

process CONCAT_SDATA {

    publishDir { "${params.outdir}/concat_sdata" },
        mode: 'copy'

    input:
    path input_paths

    output:
    path 'output/merged.zarr', emit: artifacts
    path "output/**", optional: true, hidden: true, emit: output_tree

    script:
    """
    export XDG_CACHE_HOME="\$PWD/.cache"
    export TMPDIR="\$PWD/tmp"
    mkdir -p "\$XDG_CACHE_HOME" "\$TMPDIR"

    concat_sdata.py --paths ${input_paths}
    """

    stub:
    """
    mkdir -p output/merged.zarr
    touch output/merged.zarr/.zgroup
    touch output/merged.zarr/.zattrs
    touch output/merged.zarr/.zmetadata
    """
}

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
