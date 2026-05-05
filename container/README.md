# Container Build

This directory captures the runtime used by the Xenium notebooks.
The primary local workflow is:

1. build a local `.sif` with Apptainer
2. validate the pipeline locally with `-profile test`
3. build and push the matching Docker image to a registry for OSCER

## Files

- `Apptainer.def`: native Apptainer recipe for building a local `.sif`
- `build_apptainer.sh`: wrapper that builds the local `.sif` from the correct directory
- `build_docker.sh`: wrapper that builds the Docker image from the same curated environment
- `environment.container.yml`: curated Conda environment for notebook execution
- `environment.squidpy.yml`: full exported local environment kept as a reference snapshot
- `Dockerfile`: installs Quarto and recreates the curated runtime in a `micromamba` image

## Build A Local SIF

```bash
./container/build_apptainer.sh
```

If `/tmp` is too small on your machine, point Apptainer at a larger writable temp area:

```bash
mkdir -p /home/babiddy/xenium_nb/container/.apptainer-tmp
APPTAINER_TMPDIR=/home/babiddy/xenium_nb/container/.apptainer-tmp \
TMPDIR=/home/babiddy/xenium_nb/container/.apptainer-tmp \
  ./container/build_apptainer.sh
```

You can also pass an explicit output path:

```bash
./container/build_apptainer.sh /absolute/path/to/xenium_tools_squidpy_local.sif
```

Use the SIF directly with local Apptainer-backed Nextflow runs:

```bash
nextflow run create.nf \
  -profile test \
  --container_image /absolute/path/to/container/xenium_tools_squidpy_local.sif \
  --samplesheet assets/samplesheet.csv
```

Validate the SIF directly:

```bash
apptainer exec /absolute/path/to/xenium_tools_squidpy_local.sif \
  python -c "import spatialdata, spatialdata_io, spatialdata_plot, scanpy, squidpy, nbclient, nbformat, papermill, session_info, yaml; print('ok')"
```

## Build And Push A Registry Image

Build the Docker image from the same curated environment:

```bash
./container/build_docker.sh
```

Validate it:

```bash
docker run --rm xenium_tools_squidpy:local \
  python -c "import spatialdata, spatialdata_io, spatialdata_plot, scanpy, squidpy, nbclient, nbformat, papermill, session_info, yaml; print('ok')"
```

You can also pass an explicit image tag:

```bash
./container/build_docker.sh babiddy755/xenium_nb:<tag>
```

Tag and push it to Docker Hub:

```bash
docker tag xenium_tools_squidpy:local babiddy755/xenium_nb:<tag>
docker push babiddy755/xenium_nb:<tag>
```

Then point OSCER at that tag with `--container_image` or update the default in `nextflow.config`.
