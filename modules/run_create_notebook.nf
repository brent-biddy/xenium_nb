// Renders a create-stage Quarto notebook (e.g. create_sdata, create_follicle_sdata)
// for one sample and emits the resulting zarr artifacts plus the rendered HTML.
// The output glob `output/*.zarr` accommodates both single-zarr producers and
// multi-zarr producers (Nextflow emits a Path for one match, a List<Path> for many).
process RUN_CREATE_NOTEBOOK {
    tag "${sample}:${notebook.baseName}"

    publishDir({ publish_dir }), mode: 'copy', saveAs: { fn ->
        fn.endsWith('.html') ? "${sample}_${notebook.baseName}.html" : fn
    }

    input:
    tuple path(notebook),
          path('timer.py'),
          path(input_path),
          path(cell_ids_file),
          val(sample),
          val(publish_dir),
          val(row_params),
          val(declared_params)

    output:
    tuple val(sample), path("output/*.zarr"), val(row_params), emit: artifacts
    path "*.html", emit: reports

    script:
    def paramsYaml = QuartoParams.renderParamsYaml(declared_params, row_params + [path: input_path.getName()], [
        cell_ids_file: cell_ids_file.getName(),
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
    mkdir -p output/${sample}.zarr
    touch output/${sample}.zarr/.zgroup
    touch output/${sample}.zarr/.zattrs
    touch output/${sample}.zarr/.zmetadata
    touch ${sample}_${notebook.baseName}.html
    """
}
