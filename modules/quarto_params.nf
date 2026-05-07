// Builds a Quarto params YAML string from declared per-notebook params and
// samplesheet row data. Samplesheet columns take precedence; cell_ids_file
// and radius fall back to pipeline params when not present in the row.
def renderParamsYaml(Collection declaredParams, Map rowParams) {
    def declared = declaredParams as Set
    def yaml = new org.yaml.snakeyaml.Yaml()
    def merged = new LinkedHashMap()

    merged['path'] = new File(rowParams.path as String).absolutePath

    rowParams.each { key, value ->
        if (declared.contains(key) && !merged.containsKey(key)) merged[key] = value
    }

    if (declared.contains('cell_ids_file') && !merged.containsKey('cell_ids_file'))
        merged['cell_ids_file'] = new File(params.cell_ids_file as String).absolutePath
    if (declared.contains('radius') && !merged.containsKey('radius'))
        merged['radius'] = params.radius

    return yaml.dump(merged)
}
