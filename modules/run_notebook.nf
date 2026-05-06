// Generic Quarto notebook runner. Notebook-specific params are prepared in a
// separate WRITE_QUARTO_PARAMS step so this process only executes Quarto and
// publishes the resulting reports/artifacts.
process RUN_NOTEBOOK {
    tag "${sample}:${publish_name}"

    publishDir({ publish_dir }), mode: 'copy', saveAs: { fn ->
        fn.startsWith('output/') ? fn : "${publish_name}_${fn}"
    }

    input:
    tuple path(notebook),
          val(notebook_base),
          path('timer.py'),
          path(input_path),
          val(sample),
          val(publish_dir),
          val(publish_name),
          val(row_params),
          path(cell_ids_file),
          path('params.yml')

    output:
    tuple val(sample), path('output/*.zarr'), val(row_params), optional: true, emit: artifacts
    path "${notebook_base}.*", emit: reports
    path "output/**", optional: true, hidden: true, emit: output_tree

    script:
    """
    # Redirect cache and temp dirs into the writable work dir so quarto/deno
    # don't try to write to a read-only /tmp on HPC compute nodes.
    export XDG_CACHE_HOME="\$PWD/.cache"
    export TMPDIR="\$PWD/tmp"
    mkdir -p "\$XDG_CACHE_HOME" "\$TMPDIR"

    quarto render ${notebook} --execute-params params.yml --output-dir .
    """

    stub:
    """
    mkdir -p output/${sample}.zarr
    touch output/${sample}.zarr/.zgroup
    touch output/${sample}.zarr/.zattrs
    touch output/${sample}.zarr/.zmetadata
    touch ${notebook_base}.html
    touch ${notebook_base}.pptx
    touch ${notebook_base}.timing.tsv
    """
}
