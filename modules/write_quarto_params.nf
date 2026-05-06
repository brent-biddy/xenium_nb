// Writes a Quarto params.yml file from notebook-specific and pipeline-specific
// params. This keeps YAML generation outside the notebook execution process.
process WRITE_QUARTO_PARAMS {
    tag "${sample_id}"

    input:
    tuple val(sample_id), val(input_path), val(row_params), val(declared_params)

    output:
    tuple val(sample_id), path('params.yml'), emit: params_file

    script:
    def paramsYaml = QuartoParams.renderParamsYaml(
        declared_params,
        row_params + [path: new File(input_path.toString()).name],
        [
            cell_ids_file: new File(params.cell_ids_file.toString()).name,
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
