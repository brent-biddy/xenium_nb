// Cohort-level sample summary: one pptx deck over every clustered zarr, holding each
// sample's per-cell QC distributions split by cluster. No publishDir helper or
// samplesheet row here — this step is a terminal fan-in and produces nothing another
// step consumes.

process SAMPLE_SUMMARY {
    tag "SAMPLE_SUMMARY"

    // No sample in the path — this is one fan-in task over the whole cohort.
    // 'copy' not 'link': the deck and its CSV are small and are the thing you scp off
    // the cluster, so they should survive the work dir being cleaned.
    publishDir "${params.outdir}/sample_summary", mode: 'copy'

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
    // The curated per-sample resolution table, staged flat beside the notebook. The
    // cluster count is reported at the resolution chosen for each sample, so this is a
    // required input, not a param with a default.
    //
    // stageAs pins the name the notebook reads: without it the file stages under
    // whatever basename --chosen_resolutions happened to point at, and the notebook's
    // fixed-name read would fail for every CSV not already called chosen_resolutions.csv.
    path resolutions, stageAs: 'chosen_resolutions.csv'

    output:
    path "sample_summary.pptx", emit: report

    script:
    // Redirect caches and temp files into the task dir: an OSCER compute node's /tmp is
    // read-only, so anything defaulting there fails.
    """
    export XDG_CACHE_HOME="\$PWD/.cache"
    export TMPDIR="\$PWD/tmp"
    mkdir -p "\$XDG_CACHE_HOME" "\$TMPDIR"

    quarto render ${notebook} --output-dir .
    """

    stub:
    """
    touch sample_summary.pptx
    """
}
