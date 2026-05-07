class NotebookRegistry {
    static Map create(String projectDir) {
        [
            create_sdata: [
                path  : "${projectDir}/notebooks/create_sdata.qmd",
                params: ['sample', 'path', 'he_image', 'he_alignment'],
            ],
            create_follicle_sdata: [
                path  : "${projectDir}/notebooks/create_follicle_sdata.qmd",
                params: ['sample', 'path', 'cell_ids_file', 'radius'],
            ],
        ]
    }

    static Map analysis(String projectDir) {
        [
            plot_follicle: [
                path  : "${projectDir}/notebooks/plot_follicle.qmd",
                params: ['sample', 'cell', 'path'],
            ],
        ]
    }
}
