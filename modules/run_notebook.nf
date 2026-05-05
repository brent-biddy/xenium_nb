process RUN_NOTEBOOK {
    tag "${sample}:${notebook.baseName}"

    // Publish the rendered notebook under the sample-scoped analysis directory.
    // Other outputs (zarr stores etc.) are published under output/ as written.
    publishDir({ publish_dir }), mode: 'copy', saveAs: { fn ->
        fn.startsWith("${notebook.baseName}.") ? "${sample}_${fn}" : fn
    }

    input:
    tuple path(notebook),
          path('timer.py'),       // staged into work dir so notebooks can `from timer import timer`
          path(artifact_path),    // upstream artifact staged so notebooks read it from CWD
          val(sample),
          val(publish_dir),       // resolved by the caller: <outdir>/<sample>/<notebook_basename>
          val(row_params),        // row params map; path is rewritten to the staged artifact basename
          val(declared_params)

    output:
    path "${notebook.baseName}.*"
    // hidden: true is required so dotfiles inside zarr stores (.zgroup,
    // .zattrs, .zmetadata) are collected — without it the published zarr
    // is missing root metadata and downstream readers fail.
    path "output/**", optional: true, hidden: true

    script:
    def paramsYaml = QuartoParams.renderParamsYaml(declared_params, new LinkedHashMap(row_params) + [path: artifact_path.getName()], [
        cell_ids_file: params.cell_ids_file,
        radius       : params.radius,
        n_jobs       : task.cpus,
    ])
    """
    # Redirect cache and temp dirs into the writable work dir so quarto/deno
    # don't try to write to a read-only /tmp on HPC compute nodes.
    export XDG_CACHE_HOME="\$PWD/.cache"
    export TMPDIR="\$PWD/tmp"
    mkdir -p "\$XDG_CACHE_HOME" "\$TMPDIR"

    cat > params.yml << 'PARAMS_EOF'
${paramsYaml}
PARAMS_EOF

    quarto render ${notebook} --execute-params params.yml --output-dir .
    """

    stub:
    """
    mkdir -p output
    touch ${notebook.baseName}.pptx
    touch ${notebook.baseName}.timing.tsv
    touch output/.keep
    """
}
