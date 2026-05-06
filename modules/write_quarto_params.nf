// Writes a Quarto params.yml file from notebook-specific and pipeline-specific
// params. This keeps YAML generation outside the notebook execution process.

// Builds a YAML string from declared per-notebook params and row data.
// path is always included; cell_ids_file and radius are read from pipeline
// params when declared. n_jobs is passed directly by the execution process.
def renderParamsYaml(Collection declaredParams, String inputPath, Map rowParams) {
    def declared = declaredParams as Set
    def yaml = new org.yaml.snakeyaml.Yaml()
    def merged = new LinkedHashMap()

    merged['path'] = inputPath

    if (declared.contains('cell_ids_file'))
        merged['cell_ids_file'] = new File(params.cell_ids_file as String).absolutePath
    if (declared.contains('radius'))
        merged['radius'] = params.radius

    rowParams.each { key, value ->
        if (declared.contains(key) && !merged.containsKey(key)) merged[key] = value
    }

    return yaml.dump(merged)
}

process WRITE_QUARTO_PARAMS {
    tag "${sample_id}"

    input:
    tuple val(sample_id), val(input_path), val(row_params), val(declared_params)

    output:
    tuple val(sample_id), path('params.yml'), emit: params_file

    script:
    def paramsB64 = renderParamsYaml(declared_params, input_path as String, row_params)
        .bytes.encodeBase64().toString()
    """
    python3 -c "import base64; open('params.yml', 'w').write(base64.b64decode('${paramsB64}').decode())"
    """

    stub:
    """
    touch params.yml
    """
}
