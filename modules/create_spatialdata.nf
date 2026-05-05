process CREATE_SPATIALDATA {
    tag "${sample}:${notebook.baseName}"

    publishDir({ publish_dir }), mode: 'copy', saveAs: { fn ->
        fn.endsWith('.html') ? "${sample}_${notebook.baseName}.html" : fn
    }

    input:
    tuple path(notebook),
          path('timer.py'),
          path(input_path),
          val(sample),
          val(publish_dir),
          val(row_params)

    output:
    tuple val(sample),
          path("output/${sample}.zarr"),
          val(row_params), emit: artifacts
    path "*.html", emit: reports

    script:
    def stagedRowJson = groovy.json.JsonOutput.toJson(new LinkedHashMap(row_params) + [path: input_path.getName()])
    def pipeline_params_json = groovy.json.JsonOutput.toJson([
        cell_ids_file: params.cell_ids_file,
        radius       : params.radius,
        n_jobs       : task.cpus,
    ])
    """
    export XDG_CACHE_HOME="\$PWD/.cache"
    export TMPDIR="\$PWD/tmp"
    mkdir -p "\$XDG_CACHE_HOME" "\$TMPDIR"

    cat > row_params.json << 'ROW_EOF'
${stagedRowJson}
ROW_EOF

    cat > pipeline_params.json << 'PIPELINE_EOF'
${pipeline_params_json}
PIPELINE_EOF

    python3 << 'PYEOF'
import json, re, yaml

with open('${notebook}') as f:
    content = f.read()

match = re.match(r'^---\\r?\\n(.*?)\\r?\\n---', content, re.DOTALL)
if not match:
    raise ValueError("Notebook ${notebook} has no YAML front matter")
declared = set((yaml.safe_load(match.group(1)) or {}).get('params', {}).keys())

with open('row_params.json') as f:
    full = json.load(f)

with open('pipeline_params.json') as f:
    for k, v in json.load(f).items():
        if k in declared and k not in full:
            full[k] = v

with open('params.json', 'w') as f:
    json.dump({k: v for k, v in full.items() if k in declared}, f)
PYEOF

    python3 << 'PYEOF'
import importlib.util
import sys

missing = [pkg for pkg in ("nbformat",) if importlib.util.find_spec(pkg) is None]
if missing:
    sys.stderr.write(
        "Missing Python packages required for Quarto Jupyter execution: "
        + ", ".join(missing)
        + "\\n"
    )
    sys.exit(1)
PYEOF

    quarto render ${notebook} --output-dir .
    """

    stub:
    """
    mkdir -p output/${sample}.zarr
    touch output/${sample}.zarr/.zgroup
    touch output/${sample}.zarr/.zattrs
    touch output/${sample}.zarr/.zmetadata
    touch ${sample}_${notebook.baseName}.html
    """
}
