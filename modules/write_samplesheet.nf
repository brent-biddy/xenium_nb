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
if rows and isinstance(rows[0], str):
    if len(rows) % 2 != 0:
        raise ValueError(f"Expected an even number of flat row values, got {len(rows)}")
    rows = [rows[i:i + 2] for i in range(0, len(rows), 2)]

with open('${output_name}', 'w', newline='') as fh:
    if rows and isinstance(rows[0], dict):
        headers = list(rows[0].keys())
        writer = csv.DictWriter(fh, fieldnames=headers)
        writer.writeheader()
        writer.writerows(rows)
    else:
        writer = csv.writer(fh)
        writer.writerow(['sample', 'path'])
        writer.writerows(rows)
PYEOF
    """
}
