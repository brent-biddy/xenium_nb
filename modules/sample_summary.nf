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
    // Every sample's centroid store from create_centroids, staged flat beside the zarrs.
    // This is where the deck's expression evidence comes from — it no longer opens a
    // counts matrix at all, so the zarrs above are read only for obs, obsm and the
    // spatial elements the tissue and zoom slides draw.
    //
    // No indexed stageAs here, unlike the zarrs: create_centroids names its output
    // <sample>_centroids.h5ad, so two samples cannot collide the way two `clustered.zarr`
    // do. The notebook still keys them on obs["sample"] rather than the file name, so a
    // renamed file changes nothing.
    path centroids
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
    // The follicle marker sets, likewise staged under the fixed name the notebook reads.
    // Required rather than optional: absent, the marker section would vanish silently,
    // and it covers exactly the cell types the major-type reference cannot name.
    path markers, stageAs: 'ovary_follicle_markers.yaml'
    // The major cell-type reference the cluster centroids are correlated against,
    // likewise under the fixed name the notebook reads.
    path reference, stageAs: 'ovary_reference_major.csv.gz'
    // The hand-made cell type calls. Optional in substance — the notebook falls back to
    // the argmax per cell type when the file has no row for a (sample, cell type) — but
    // staged unconditionally, since Nextflow has no optional path input and the shipped
    // asset is a header-only sheet.
    path cell_types, stageAs: 'cell_type_annotations.csv'
    // The hand-curated oocyte list, likewise under the fixed name the notebook reads. It
    // drives the per-stage zoom slides; a sample with no curated cells simply gets none,
    // which the notebook says on a slide rather than passing over in silence.
    path oocytes, stageAs: 'curated_oocytes.csv'

    output:
    path "sample_summary.pptx", emit: report
    // The draft cell type map: argmax per cell type, with both cluster numberings and
    // the score behind each call. Edited by hand and fed back via --cell_type_annotations.
    path "cluster_annotations.csv", emit: annotations
    // Cell type composition per sample, the same numbers as the deck's tables — so a
    // cohort's composition can be read without opening the pptx.
    path "cell_type_composition.csv", emit: composition

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
    touch cluster_annotations.csv
    touch cell_type_composition.csv
    """
}
