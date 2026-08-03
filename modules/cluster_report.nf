// Cohort-level clustering report: one pptx deck over every clustered zarr, rather than
// a per-sample artifact. No publishDir helper or samplesheet row here — this step is a
// terminal fan-in and produces nothing another step consumes.

process CLUSTER_REPORT {
    tag "CLUSTER_REPORT"

    // No sample in the path — this is one fan-in task over the whole cohort.
    // 'copy' not 'link': the report is small and is the thing you scp off the
    // cluster, so it should survive the work dir being cleaned.
    publishDir "${params.outdir}/cluster_report", mode: 'copy'

    input:
    // Every sample's clustered zarr, staged flat into the work dir. The notebook globs
    // *.zarr from its own directory, so staging IS the input contract — there are no
    // params to pass and nothing to keep in sync with the notebook, which is why this
    // notebook has no entry in assets/notebook_registry.json.
    //
    // stageAs with an index is required, not cosmetic: cluster_sdata* publish every
    // sample's output as `clustered.zarr`, so staging them under their own names is a
    // hard "input file name collision" error the moment there are two samples. The
    // staged name carries no meaning as a result — the notebook takes each sample's id
    // from obs["sample"], not from the path.
    path zarrs, stageAs: 'sample*.zarr'
    path notebook
    // Quarto resolves reference-doc relative to the qmd's own directory, so the
    // template has to land beside the staged notebook, not at its repo path.
    path template

    output:
    path "cluster_report.pptx", emit: report

    script:
    """
    export XDG_CACHE_HOME="\$PWD/.cache"
    export TMPDIR="\$PWD/tmp"
    mkdir -p "\$XDG_CACHE_HOME" "\$TMPDIR"

    quarto render ${notebook} --output-dir .
    """

    stub:
    """
    touch cluster_report.pptx
    """
}
