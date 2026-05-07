// Writes a Quarto params YAML file to <outdir>/.params/ and returns its Path
// for channel staging. Writing to outdir (on shared NFS on HPC) avoids the
// symlink-to-local-/tmp problem that breaks staged files on compute nodes.
def paramsFile(String id, Collection declaredParams, Map rowParams, def outdir) {
    def paramsDir = new File("${outdir}/.quarto_params")
    paramsDir.mkdirs()
    def paramsFile = new File(paramsDir, "params_${id}.yml")
    paramsFile.text = renderParamsYaml(declaredParams, rowParams)
    return paramsFile.toPath()
}

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
