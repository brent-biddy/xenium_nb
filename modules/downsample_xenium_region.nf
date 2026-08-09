process DOWNSAMPLE_XENIUM_REGION {
    tag "${sample}/${region_name}"

    // Hardlink (not copy) into results: workDir and outdir share the scratch
    // filesystem, so linking avoids a second full copy of the cropped region.
    publishDir { "${params.outdir}/${sample}/downsample_xenium_region" },
        mode: 'link'

    // he_image/he_alignment are `path`, not `val`, and are optional — main.nf passes
    // `[]` when the samplesheet omits them, which stages nothing and leaves the
    // variable falsy. They must be staged: an unstaged absolute path is a host path
    // the container cannot see. Apptainer's default binds ($HOME, /tmp, cwd) hide
    // this locally, but an H&E anywhere else — e.g. OSCER /scratch, which is bound
    // only at the work dir Nextflow mounts — is simply missing inside the container.
    // Staging is safe against the output glob: inputs land at the work dir root,
    // while the published artifacts are everything under ${region_name}/.
    input:
    tuple val(sample), path(input_path), val(xmin), val(ymin), val(xmax), val(ymax), val(region_name), path(he_image), path(he_alignment)

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
