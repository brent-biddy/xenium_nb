process CREATE_SPATIALDATA {
    tag "${sample}:${notebook.baseName}"

    publishDir({ publish_dir }), mode: 'copy', saveAs: { fn ->
        fn.endsWith('.html') ? "${sample}_${notebook.baseName}.html" : fn
    }

    input:
    tuple path(notebook),
          path('timer.py'),
          path(input_path),
          val(sample),
          val(publish_dir),
          val(row_params),
          val(declared_params)

    output:
    tuple val(sample),
          path("output/${sample}.zarr"),
          val(row_params), emit: artifacts
    path "*.html", emit: reports

    script:
    def paramsYaml = QuartoParams.renderParamsYaml(declared_params, new LinkedHashMap(row_params) + [path: input_path.getName()], [
        cell_ids_file: params.cell_ids_file,
        radius       : params.radius,
        n_jobs       : task.cpus,
    ])
    """
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
    mkdir -p output/${sample}.zarr
    touch output/${sample}.zarr/.zgroup
    touch output/${sample}.zarr/.zattrs
    touch output/${sample}.zarr/.zmetadata
    touch ${sample}_${notebook.baseName}.html
    """
}
