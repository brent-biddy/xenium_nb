// Writes a CSV samplesheet from a list of row maps. Downstream workflows use
// the published CSV to locate artifacts built by create.nf.

process WRITE_SAMPLESHEET {
    tag "${output_name}"
    publishDir params.outdir, mode: 'copy'

    input:
    tuple val(output_name), val(rows)

    output:
    path output_name

    script:
    def keys = rows[0].keySet().toList()
    def csv = ([keys.join(',')] + rows.collect { row -> keys.collect { k -> row[k] ?: '' }.join(',') }).join('\n')
    """
    cat > ${output_name} <<'CSVEOF'
${csv}
CSVEOF
    """
}
