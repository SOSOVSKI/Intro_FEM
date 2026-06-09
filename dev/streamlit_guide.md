# Streamlit local and container deployment guide

This guide documents how to run the `IntroFem` Streamlit app locally, how to package it in Docker so the runtime behaves like a cloud deployment, and how to promote the same container image to a cloud service.

## Short answer: is conda needed for the Streamlit app?

No — **conda is not required for the Streamlit app itself**.

Why:

- The runnable app is `streamlit_app.py`.
- The app imports `symbolic_fem_workbench`, `apps.presets`, SymPy, and Streamlit modules.
- The symbolic workflow helpers used by the app live in `src/symbolic_fem_workbench/workflow.py` and are pure Python/SymPy.
- The repository's FEniCSx dependency only appears in course material and notebooks for the FEniCS chapters, not in the Streamlit UI execution path.
- A local smoke test succeeded without conda by importing `streamlit`, `symbolic_fem_workbench`, and `streamlit_app`, and by starting the app through `uv run streamlit run streamlit_app.py`.

Use conda only if you also need the broader course environment for FEniCSx-backed notebooks/chapters.

## What you need for the app

- Python `3.11`
- `uv` (recommended for local development in this repo)
- Docker Desktop or another Docker engine for containerized local deployment

No `.env` file is currently required for the Streamlit app.

## Local development run (no Docker)

From the repository root:

```bash
uv sync --dev
uv run streamlit run streamlit_app.py
```

Open:

```text
http://localhost:8501
```

### Useful variations

Run headless on a different port:

```bash
uv run streamlit run streamlit_app.py --server.headless true --server.port 8502
```

Smoke-test imports only:

```bash
uv run python -c "import streamlit, symbolic_fem_workbench, streamlit_app; print('imports-ok')"
```

## Why the Docker path is closer to cloud hosting

The Docker image in this repo is designed to behave like a small production service:

- the package is installed into the image from a built wheel
- the app runs as a non-root user
- the server binds to `0.0.0.0`
- the image honors the standard `PORT` environment variable
- the image exposes a health endpoint through Streamlit's `/_stcore/health`

That means the same image you run locally is suitable for pushing to a registry and deploying to most container platforms.

## Files added for container deployment

- `Dockerfile` — two-stage build for the Streamlit app
- `.dockerignore` — reduces build context size and avoids copying local/editor artefacts

## Build the Docker image locally

From the repository root:

```bash
docker build -t introfem-streamlit:local .
```

If you want to mirror a typical Linux cloud target from Apple Silicon, build for `linux/amd64` explicitly:

```bash
docker build --platform linux/amd64 -t introfem-streamlit:local .
```

## Run the container locally

Basic local run:

```bash
docker run --rm -p 8501:8501 introfem-streamlit:local
```

Then open:

```text
http://localhost:8501
```

### Run with an explicit cloud-style `PORT`

```bash
docker run --rm -e PORT=8501 -p 8501:8501 introfem-streamlit:local
```

### Run detached

```bash
docker run -d --name introfem-streamlit -e PORT=8501 -p 8501:8501 introfem-streamlit:local
```

Stop it later:

```bash
docker stop introfem-streamlit
```

## Validate the container locally

### Browser check

Open the app in a browser and make sure the main pages load:

- `Guided Workflow`
- `Load Preset`
- `Example Demos`
- `Module Surface`

### Health endpoint check

```bash
curl http://localhost:8501/_stcore/health
```

Expected response is plain text similar to:

```text
ok
```

### Logs

If the container is running detached:

```bash
docker logs -f introfem-streamlit
```

## Troubleshooting

### `ModuleNotFoundError: symbolic_fem_workbench`

For non-Docker local runs, sync/install the project in the repo environment first:

```bash
uv sync --dev
```

Also make sure you launch from the repository root.

### Port already in use

Use another local port mapping:

```bash
docker run --rm -p 8502:8501 introfem-streamlit:local
```

Then browse to `http://localhost:8502`.

### Docker image works locally but fails in cloud

Check these first:

- the platform expects the app to bind to `0.0.0.0`
- the platform sets a custom `PORT`
- the service target port is configured to match the container port
- the registry image was built for the cloud CPU architecture you selected

## Deployment instructions

The recommended promotion flow is:

1. build the image locally
2. validate it locally with `docker run`
3. tag it for your registry
4. push it
5. deploy the same image to your cloud container service

### Tag the image

Replace the placeholders below with your registry and tag:

```bash
docker tag introfem-streamlit:local <registry>/<namespace>/introfem-streamlit:<tag>
```

### Push the image

```bash
docker push <registry>/<namespace>/introfem-streamlit:<tag>
```

### Deploy on a generic container platform

Use these settings on the cloud side:

- **Image**: `<registry>/<namespace>/introfem-streamlit:<tag>`
- **Container port**: `8501`
- **Environment variable**: `PORT=8501` if the platform does not inject it automatically
- **Startup command override**: not needed
- **Health check path**: `/_stcore/health`

### Good cloud targets for this image

This image should work on any service that runs OCI/Docker containers, for example:

- Azure Container Apps
- Azure App Service for Containers
- AWS App Runner / ECS / Fargate
- Google Cloud Run
- Render
- Fly.io
- Railway

## Recommended pre-cloud checklist

Before pushing to the cloud, confirm all of the following locally:

- `uv run streamlit run streamlit_app.py` works without conda
- `docker build -t introfem-streamlit:local .` succeeds
- `docker run --rm -p 8501:8501 introfem-streamlit:local` starts the app
- `http://localhost:8501/_stcore/health` returns `ok`
- the main Streamlit pages render without import/runtime errors

## Summary

For this repository:

- **Conda is not needed for the Streamlit app**
- **Conda is only needed for the FEniCSx course material path**
- **Docker is now the preferred local deployment path before cloud rollout** because it gives you an immutable, cloud-like runtime
