// Writes a CSV samplesheet from a JSON-encoded list of row maps. Downstream
// workflows use the published CSV to locate artifacts built by create.nf.

process WRITE_SAMPLESHEET {
    tag "${output_name}"
    publishDir params.outdir, mode: 'copy'

    input:
    tuple val(output_name), val(rows_json)

    output:
    path output_name

    script:
    """
    python3 << 'PYEOF'
import csv
import json

rows = json.loads('''${rows_json}''')
with open('${output_name}', 'w', newline='') as fh:
    writer = csv.DictWriter(fh, fieldnames=list(rows[0].keys()))
    writer.writeheader()
    writer.writerows(rows)
PYEOF
    """

    stub:
    """
    touch ${output_name}
    """
}
