process RUN_NOTEBOOK {
    tag "${sample}:${notebook.baseName}"

    // Publish HTML under <roi_id>/<notebook>/ with the sample-scoped filename.
    // Other outputs (zarr stores etc.) are published under output/ as written.
    publishDir "${publish_dir}", mode: 'copy', saveAs: { fn ->
        fn.endsWith('.html') ? output_name : fn
    }

    input:
    tuple path(notebook),
          path('timer.py'),       // staged into work dir so notebooks can `from timer import timer`
          path(input_path),       // samplesheet 'path' column, staged so notebooks read it from CWD
          val(sample),
          val(publish_dir),       // resolved in main.nf: <outdir>/<roi_id>/<notebook_basename>
          val(output_name),       // resolved in main.nf: <sample>_<notebook_basename>.html
          val(params_json)        // params JSON built from the samplesheet row (path = staged basename)

    output:
    path "*.html"
    // hidden: true is required so dotfiles inside zarr stores (.zgroup,
    // .zattrs, .zmetadata) are collected — without it the published zarr
    // is missing root metadata and downstream readers fail.
    path "output/**", optional: true, hidden: true

    script:
    // Pipeline-level params are written only if the notebook declares them.
    // n_jobs comes from task.cpus so it tracks whatever the executor allocated.
    def pipeline_params_json = groovy.json.JsonOutput.toJson([
        cell_ids_file: params.cell_ids_file,
        radius       : params.radius,
        n_jobs       : task.cpus,
    ])
    """
    # Redirect cache and temp dirs into the writable work dir so quarto/deno
    # don't try to write to a read-only /tmp on HPC compute nodes.
    export XDG_CACHE_HOME="\$PWD/.cache"
    export TMPDIR="\$PWD/tmp"
    mkdir -p "\$XDG_CACHE_HOME" "\$TMPDIR"

    cat > row_params.json << 'ROW_EOF'
${params_json}
ROW_EOF

    cat > pipeline_params.json << 'PIPELINE_EOF'
${pipeline_params_json}
PIPELINE_EOF

    python3 << 'PYEOF'
import json, re, yaml

with open('${notebook}') as f:
    content = f.read()

# Match YAML front matter, tolerating both Unix and Windows line endings.
match = re.match(r'^---\\r?\\n(.*?)\\r?\\n---', content, re.DOTALL)
if not match:
    raise ValueError("Notebook ${notebook} has no YAML front matter")
declared = set((yaml.safe_load(match.group(1)) or {}).get('params', {}).keys())

with open('row_params.json') as f:
    full = json.load(f)

# Merge in pipeline-level params only when the notebook declares them.
# Row params win over pipeline defaults if both are present.
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
}
