// Reduce each clustered zarr to per-group sums, so the report decks never open a counts
// matrix: ~1.1 GB per clustered zarr against ~6.5 MB of centroids on a test ROI, and
// sample_summary was re-reading all of it every render to recompute the same numbers.
// Per-sample fan-out, the same shape as CLUSTER_SDATA — one zarr in, one h5ad out,
// publishing its own handoff row.
//
// Deliberately its own step rather than folded into CLUSTER_SDATA, which already holds
// the matrix in memory: Nextflow hashes the task script, so folding it in would make any
// change to the centroid recipe re-run PCA, UMAP and the whole Leiden sweep for the
// cohort. See bin/create_centroids.py for what the artifact carries and why it stores
// sums rather than means.
//
// --group_by picks the obs column the cells are summed over. Unset (the default) means
// the leiden sweep, which is the store every deck reads; set, it names a column some
// upstream step wrote — per-cell cell type being the case it exists for. The two land in
// the same published dir under different names, so a grouping run never displaces the
// sweep and both stay available. See the script's docstring for when --group_by is NOT
// the answer: a grouping that only coarsens the sweep's clusters is a row-wise add in
// the deck.

// The artifact stem for a run, and with it every output name. Single-sourced because the
// output block, the samplesheet row, and the stub must all agree on it, and a --group_by
// run changes all three at once.
def centroidStem(sample) {
    params.group_by ? "${sample}_centroids_${params.group_by}" : "${sample}_centroids"
}

// Published output directory for this step's per-sample artifacts. Single-sourced here
// so the publishDir directive, the emitted h5ad, and the handoff samplesheet row all
// reference the same location and cannot drift apart.
def createCentroidsPublishDir(sample) {
    "${params.outdir}/${sample}/create_centroids"
}

process CREATE_CENTROIDS {
    tag "${sample}"

    // saveAs drops the per-sample row fragment from the published dir; it is only
    // needed on the channel for main.nf to collectFile into the aggregate sheet.
    // Hardlink (not copy) into results: workDir and outdir share the scratch
    // filesystem, so linking avoids a second copy.
    publishDir { createCentroidsPublishDir(sample) },
        mode: 'link',
        saveAs: { fn -> fn.endsWith('.samplesheet_row.csv') ? null : fn }

    input:
    // No stageAs needed here, unlike SAMPLE_SUMMARY: this is a per-sample fan-out, so
    // only one zarr is ever staged into a task dir and there is nothing to collide with.
    //
    // cluster_path rides along as a val because inside the task `input_path` resolves to
    // a staged symlink that means nothing outside this work dir, and the handoff sheet
    // has to forward a location the next step can still resolve. sample_summary needs
    // BOTH artifacts — the centroids for expression, and the clustered zarr for obs,
    // obsm and the spatial elements its tissue slides draw — so the row carries both.
    tuple val(sample), path(input_path), val(cluster_path)

    output:
    // One row per group. In the default (sweep) mode that is every (resolution, cluster)
    // across the whole sweep, not just the resolution in assets/chosen_resolutions.csv —
    // that stays a report-level knob, so revising a sample's chosen resolution re-renders
    // a deck rather than re-running this step. X is summed CP10K and layers["counts"]
    // summed raw counts; both are sums so a union of clusters is a row-wise add, and
    // n_cells in obs makes the mean-of-normalized centroid one division away.
    tuple val(sample), path("${centroidStem(sample)}.h5ad"), emit: centroids
    path "${centroidStem(sample)}_timing.tsv", emit: timing
    path "${centroidStem(sample)}_session_info.txt", emit: session_info
    // One `sample,path,centroid_path` line; main.nf collectFiles these into a
    // ready-to-use handoff samplesheet. No trailing newline — the
    // collectFile(newLine: true) call adds the separator.
    path "${sample}.samplesheet_row.csv", emit: samplesheet_row

    script:
    def centroidArgs = ["--sample ${sample}", "--path ${input_path}"]
    // Omitted when unset so the default grouping keeps its single definition in the
    // Python script rather than being restated in the config.
    if (params.group_by) centroidArgs << "--group_by ${params.group_by}"
    def published = "${createCentroidsPublishDir(sample)}/${centroidStem(sample)}.h5ad"
    """
    export XDG_CACHE_HOME="\$PWD/.cache"
    export TMPDIR="\$PWD/tmp"
    mkdir -p "\$XDG_CACHE_HOME" "\$TMPDIR"

    create_centroids.py ${centroidArgs.join(' ')}

    printf '%s' '${sample},${cluster_path},${published}' > ${sample}.samplesheet_row.csv
    """

    stub:
    // Every declared output, or the stub run fails on the missing ones — which makes the
    // documented wiring check useless. Keep in step with the output block above.
    """
    touch ${centroidStem(sample)}.h5ad
    touch ${centroidStem(sample)}_timing.tsv
    touch ${centroidStem(sample)}_session_info.txt

    printf '%s' '${sample},${cluster_path},${createCentroidsPublishDir(sample)}/${centroidStem(sample)}.h5ad' > ${sample}.samplesheet_row.csv
    """
}
