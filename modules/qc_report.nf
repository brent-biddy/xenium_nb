// Cohort-level QC report: one pptx deck over every raw create_sdata zarr, read
// BEFORE clustering so the filtering thresholds clustering will use are chosen from
// the data rather than inherited. Like cluster_report this is a terminal fan-in over
// the whole cohort, so there is no publishDir helper and no `sample,path` handoff row.
//
// The deck is the only output. It reports what the data looks like and stops there —
// it neither infers a threshold nor emits one, so acting on it means reading the slides
// and setting the cut in cluster_sdata* yourself.

process QC_REPORT {
    tag "QC_REPORT"

    // No sample in the path — this is one fan-in task over the whole cohort.
    // 'copy' not 'link': the deck is small and is the thing you scp off the cluster,
    // so it should survive the work dir being cleaned.
    publishDir "${params.outdir}/qc_report", mode: 'copy'

    input:
    // Every sample's raw zarr, staged flat into the work dir. The notebook globs
    // *.zarr from its own directory, so staging IS the input contract — there are no
    // params to pass and nothing to keep in sync with the notebook, which is why this
    // notebook has no entry in assets/notebook_registry.json.
    //
    // stageAs with an index for the same reason cluster_report needs it: create_sdata
    // publishes as <sample>.zarr, which is unique today, but concat_sdata and
    // downsample_sdata both publish fixed names, and a flat fan-in of those would be a
    // hard "input file name collision". Indexing makes the staged name meaningless
    // either way — the notebook takes each sample's id from obs["sample"], not the path.
    path zarrs, stageAs: 'sample*.zarr'
    path notebook
    // Quarto resolves reference-doc relative to the qmd's own directory, so the
    // template has to land beside the staged notebook, not at its repo path.
    path template

    output:
    path "qc_report.pptx", emit: report

    script:
    """
    export XDG_CACHE_HOME="\$PWD/.cache"
    export TMPDIR="\$PWD/tmp"
    mkdir -p "\$XDG_CACHE_HOME" "\$TMPDIR"

    quarto render ${notebook} --output-dir .
    """

    stub:
    """
    touch qc_report.pptx
    """
}
