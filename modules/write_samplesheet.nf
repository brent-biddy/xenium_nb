process WRITE_SAMPLESHEET {
    tag "${output_name}"

    publishDir "${params.outdir}/pipeline_info", mode: 'copy'

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
    writer = csv.writer(fh)
    writer.writerow(['sample', 'path'])
    writer.writerows(rows)
PYEOF
    """
}
