// Writes a Quarto params YAML file to <outdir>/.quarto_params/<notebook>/
// and returns its Path for channel staging. Writing to outdir (on shared NFS
// on HPC) avoids the symlink-to-local-/tmp problem that breaks staged files
// on compute nodes.
def paramsFile(String id, Path notebook, Map rowParams) {
    def notebookName = notebook.fileName.toString().replaceAll('\\.qmd$', '')
    def registry = new groovy.json.JsonSlurper()
        .parse(new File("${projectDir}/assets/notebook_registry.json"))
    def entry = registry[notebookName]
    if (!entry) throw new IllegalArgumentException("Notebook '${notebookName}' not found in registry")
    def paramsDir = new File("${params.outdir}/.quarto_params/${notebookName}")
    paramsDir.mkdirs()
    def paramsFile = new File(paramsDir, "params_${id}.yml")
    paramsFile.text = renderParamsYaml(entry.params, rowParams)
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
