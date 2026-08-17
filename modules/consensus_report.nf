// Cohort-wide consensus clustering: one pptx deck grouping every sample's clusters into
// consensus clusters, plus the two sheets that let a grouping be curated and committed.
//
// A SIBLING of SAMPLE_SUMMARY, not a step downstream of it — both read create_centroids'
// stores directly, so neither blocks the other. The split exists because consensus
// grouping is the part that gets iterated on: inside sample_summary, changing the linkage
// or the cut meant re-rendering every sample's fourteen per-sample slides to look at one
// heatmap. This deck reads centroids only, so it renders in seconds.
//
// No publishDir helper or samplesheet row: like the other decks this is a terminal
// fan-in. Its two CSVs are edited by hand and committed into assets/, which is what makes
// a grouping authoritative — nothing reads them back out of a run directory.

process CONSENSUS_REPORT {
    tag "CONSENSUS_REPORT"

    // No sample in the path — one fan-in task over the whole cohort. 'copy' not 'link':
    // the deck and its sheets are small and are what you scp off the cluster, so they
    // should survive the work dir being cleaned.
    publishDir "${params.outdir}/consensus_report", mode: 'copy'

    input:
    // Every sample's centroid store, staged flat. The notebook globs *_centroids.h5ad
    // from its own directory, so staging IS the input contract — no params to pass and
    // nothing to keep in sync, which is why this notebook has no registry entry.
    //
    // No indexed stageAs: create_centroids names its output <sample>_centroids.h5ad, so
    // two samples cannot collide the way two `clustered.zarr` do.
    path centroids
    path notebook
    // Quarto resolves reference-doc relative to the qmd's own directory, so the template
    // has to land beside the staged notebook, not at its repo path.
    path template
    // The curated per-sample resolution table. The cohort matrix is built at the
    // resolution chosen for each sample, so this decides which rows of each store are
    // read — required, not a param with a default.
    path resolutions, stageAs: 'chosen_resolutions.csv'
    // The two evidence assets, the same ones sample_summary uses and staged under the
    // same fixed names. A consensus cluster is named by the rule that names a per-sample
    // cluster, so both decks read one definition of each family.
    path markers, stageAs: 'ovary_follicle_markers.yaml'
    path reference, stageAs: 'ovary_reference_major.csv.gz'
    // The immune subtype centroids, for the one panel that asks what the group called
    // Immune is made of. Not a subtyping section — see the notebook for why the
    // compartment does not split into subtypes at any cut.
    path immune_reference, stageAs: 'ovary_reference_immune.csv.gz'
    // The committed consensus grouping, read by fixed name. OPTIONAL IN SUBSTANCE: it
    // does not exist as a real grouping until one has been promoted from
    // consensus_scaffold.csv and committed, and this deck is what you render to produce
    // that scaffold in the first place. Nextflow has no optional path input, so the
    // shipped asset is HEADER-ONLY — the same convention cell_type_annotations.csv uses,
    // and simpler than staging a sentinel file under a different name. The notebook reads
    // zero rows as "no grouping exists yet" and proposes one.
    path grouping, stageAs: 'consensus_clusters.csv'

    output:
    path "consensus_report.pptx", emit: report
    // One row per (sample, cluster) carrying the PROPOSED grouping. Published so it can
    // be edited by hand against the deck, then promoted to assets/consensus_clusters.csv
    // and committed — the commit is what makes a grouping authoritative, which is why
    // nothing reads this back from a run dir.
    path "consensus_scaffold.csv", emit: scaffold
    // One row per consensus cluster with its cell type call and the evidence behind it.
    // Always written, so this declared output exists whether or not a grouping is staged.
    path "consensus_annotations.csv", emit: annotations

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
    // Every declared output, or the stub run fails on the missing ones — which makes the
    // documented wiring check useless. Keep in step with the output block above.
    """
    touch consensus_report.pptx
    touch consensus_scaffold.csv
    touch consensus_annotations.csv
    """
}
