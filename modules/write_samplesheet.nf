// Writes a CSV samplesheet from a list of row maps. Downstream workflows use
// the published CSV to locate artifacts built by create.nf.

process WRITE_SAMPLESHEET {
    tag "${output_name}"
    publishDir params.outdir, mode: 'copy'

    input:
    tuple val(output_name), val(rows)

    output:
    path output_name

    exec:
    def keys = rows[0].keySet().toList()
    task.workDir.resolve(output_name).text =
        ([keys.join(',')] + rows.collect { row -> keys.collect { k -> row[k] ?: '' }.join(',') }).join('\n') + '\n'
}
