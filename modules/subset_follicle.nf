process SUBSET_FOLLICLE {
    tag "${sample}:${notebook.baseName}"

    publishDir({ publish_dir }), mode: 'copy', saveAs: { fn ->
        fn.endsWith('.html') ? "${sample}_${notebook.baseName}.html" : fn
    }

    input:
    tuple val(sample),
          path(sample_zarr),
          val(row_params),
          path(cell_ids_file),
          path(notebook),
          path('timer.py'),
          val(publish_dir),
          val(declared_params)

    output:
    tuple val(sample),
          path("output/*.zarr"), emit: artifacts
    path "*.html", emit: reports

    script:
    def paramsYaml = QuartoParams.renderParamsYaml(declared_params, new LinkedHashMap(row_params) + [path: sample_zarr.getName()], [
        cell_ids_file: cell_ids_file.getName(),
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
