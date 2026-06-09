FROM python:3.11-slim AS builder

ENV PYTHONDONTWRITEBYTECODE=1 \
    PYTHONUNBUFFERED=1 \
    PIP_NO_CACHE_DIR=1

WORKDIR /app

COPY pyproject.toml README.md ./
COPY src ./src

RUN python -m pip install --upgrade pip build && \
    python -m build --wheel --outdir /dist


FROM python:3.11-slim AS runtime

ENV PYTHONDONTWRITEBYTECODE=1 \
    PYTHONUNBUFFERED=1 \
    PIP_NO_CACHE_DIR=1 \
    STREAMLIT_BROWSER_GATHER_USAGE_STATS=false \
    STREAMLIT_SERVER_HEADLESS=true \
    PORT=8501

WORKDIR /app

RUN useradd --create-home --shell /bin/bash appuser

COPY --from=builder /dist/*.whl /tmp/dist/
RUN python -m pip install --upgrade pip && \
    python -m pip install /tmp/dist/*.whl && \
    rm -rf /tmp/dist

COPY apps ./apps
COPY streamlit_app.py ./

USER appuser

EXPOSE 8501

HEALTHCHECK --interval=30s --timeout=5s --start-period=20s --retries=3 \
  CMD python -c "import os, urllib.request; port = os.environ.get('PORT', '8501'); urllib.request.urlopen(f'http://127.0.0.1:{port}/_stcore/health', timeout=3)"

CMD ["sh", "-c", "python -m streamlit run streamlit_app.py --server.address=0.0.0.0 --server.port=${PORT:-8501}"]
