process WRITE_SAMPLESHEET {
    tag "${output_name}"
    publishDir({ publish_dir }), mode: 'copy'

    input:
    tuple val(output_name), val(rows_json), val(publish_dir)

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
