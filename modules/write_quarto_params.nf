// Writes a Quarto params.yml file from notebook-specific and pipeline-specific
// params. This keeps YAML generation outside the notebook execution process.
process WRITE_QUARTO_PARAMS {
    tag "${sample}:${publish_name}"

    input:
    tuple val(notebook),
          val(notebook_base),
          val(timer_script),
          val(input_path),
          val(sample),
          val(publish_dir),
          val(publish_name),
          val(row_params),
          val(cell_ids_file),
          val(declared_params)

    output:
    tuple val(notebook),
          val(notebook_base),
          val(timer_script),
          val(input_path),
          val(sample),
          val(publish_dir),
          val(publish_name),
          val(row_params),
          val(cell_ids_file),
          path('params.yml'), emit: notebook_inputs

    script:
    def paramsYaml = QuartoParams.renderParamsYaml(
        declared_params,
        row_params + [path: new File(input_path.toString()).name],
        [
            cell_ids_file: new File(cell_ids_file.toString()).name,
            radius       : params.radius,
            n_jobs       : task.cpus,
        ]
    )
    """
    cat > params.yml << 'PARAMS_EOF'
${paramsYaml}
PARAMS_EOF
    """

    stub:
    """
    touch params.yml
    """
}
