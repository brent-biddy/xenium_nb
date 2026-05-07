// Builds a base64-encoded Quarto params YAML string from declared per-notebook
// params and samplesheet row data. Consumed directly by notebook processes.
def renderParamsYaml(Collection declaredParams, Map rowParams) {
    def declared = declaredParams as Set
    def yaml = new org.yaml.snakeyaml.Yaml()
    def merged = new LinkedHashMap()

    merged['path'] = new File(rowParams.path as String).absolutePath

    if (declared.contains('cell_ids_file'))
        merged['cell_ids_file'] = new File(params.cell_ids_file as String).absolutePath
    if (declared.contains('radius'))
        merged['radius'] = params.radius

    rowParams.each { key, value ->
        if (declared.contains(key) && !merged.containsKey(key)) merged[key] = value
    }

    return yaml.dump(merged)
}
