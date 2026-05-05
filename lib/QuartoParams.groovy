class QuartoParams {
    private static final org.yaml.snakeyaml.Yaml YAML = new org.yaml.snakeyaml.Yaml()

    static String renderParamsYaml(Collection declaredParams, Map rowParams = [:], Map pipelineParams = [:]) {
        def declared = declaredParams as Set
        def merged = new LinkedHashMap()

        rowParams.each { key, value ->
            if (declared.contains(key)) {
                merged[key] = value
            }
        }

        pipelineParams.each { key, value ->
            if (declared.contains(key) && !merged.containsKey(key)) {
                merged[key] = value
            }
        }

        return YAML.dump(merged)
    }
}
