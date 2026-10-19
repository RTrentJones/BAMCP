# ============================================================================
# BAMCP Production Dockerfile
# Multi-stage build: compile pysam from source, then copy into slim runtime.
# ============================================================================

# ---------------------------------------------------------------------------
# Stage 1a — Frontend build (Vite bundles viewer into single HTML)
# ---------------------------------------------------------------------------
# Base images are pinned by multi-arch index digest (amd64 + arm64) so a rebuild of the same
# commit gets the same base; Dependabot (docker ecosystem) bumps the digests.
FROM node:26-slim@sha256:193fe51b64e77981119c98c2002c9e32a70e2f006fb4d25068ce0558998917f0 AS frontend

WORKDIR /frontend
COPY src/bamcp/static/package.json src/bamcp/static/package-lock.json* ./
RUN npm ci --ignore-scripts
COPY src/bamcp/static/ ./
RUN npm run build

# ---------------------------------------------------------------------------
# Stage 1b — Python build
# ---------------------------------------------------------------------------
FROM python:3.14-slim@sha256:a2b82f3c48559aa0a8446d9af49826b6e2b2016f4cd2afabfe6013ec53729170 AS builder

RUN apt-get update && apt-get install -y --no-install-recommends \
        gcc \
        libc6-dev \
        zlib1g-dev \
        libbz2-dev \
        liblzma-dev \
        libcurl4-openssl-dev \
        libssl-dev \
        make \
    && rm -rf /var/lib/apt/lists/*

WORKDIR /build

# uv is used only to export the lockfile; it never reaches the runtime stage. Copied straight from
# its digest-pinned image (not a separate FROM stage, so every FROM stays a slim base). Bump the
# tag + digest by hand alongside the uv version used to write uv.lock.
COPY --from=ghcr.io/astral-sh/uv:0.11.32@sha256:df4cae8f3a96d175e2e5f992e597550000edbe78fdc2594d5cd8de1a217f504c /uv /usr/local/bin/uv

COPY pyproject.toml uv.lock README.md ./
COPY src/ src/

# Copy bundled viewer into source tree so hatchling includes it in the wheel
COPY --from=frontend /frontend/dist/viewer.html src/bamcp/static/dist/viewer.html

# Install the dependencies FROM uv.lock (exact versions + hashes) — the same set CI tests and
# pip-audit checks — rather than letting pip resolve the newest versions within pyproject's
# ranges. `--frozen` fails the build if the lock is stale instead of silently re-resolving.
# Include [telemetry] so the OTel SDK ships with the image — JSONL is built-in, but the
# OpenTelemetry exporter requires the extras. [eval] is intentionally excluded: it pulls in
# anthropic/openai/playwright which are dev-only. The project itself is then installed with
# --no-deps so pip cannot pull anything the lock did not pin.
RUN uv export --frozen --no-dev --extra telemetry --no-emit-project -o /tmp/requirements.txt \
    && pip install --no-cache-dir --prefix=/install --require-hashes -r /tmp/requirements.txt \
    && pip install --no-cache-dir --prefix=/install --no-deps .

# ---------------------------------------------------------------------------
# Stage 2 — Runtime
# ---------------------------------------------------------------------------
FROM python:3.14-slim@sha256:a2b82f3c48559aa0a8446d9af49826b6e2b2016f4cd2afabfe6013ec53729170 AS runtime

# pysam needs the htslib shared libs at runtime
RUN apt-get update && apt-get install -y --no-install-recommends \
        zlib1g \
        libbz2-1.0 \
        liblzma5 \
        libcurl4 \
    && rm -rf /var/lib/apt/lists/*

# Copy installed packages from builder
COPY --from=builder /install /usr/local

# Copy source (needed for src layout discovery)
WORKDIR /app
COPY src/ src/
COPY pyproject.toml ./
COPY docker/ docker/

# Create a non-root user with a proper writable home directory.
# Libraries (pysam, cffi, etc.) need $HOME/.cache for runtime caches.
RUN groupadd -r bamcp && useradd -r -g bamcp -m -d /home/bamcp bamcp \
    && mkdir -p /data && chown bamcp:bamcp /data \
    && chmod +x /app/docker/entrypoint.sh

ENV HOME=/home/bamcp
USER bamcp

# Default env vars (can be overridden). Calling thresholds (BAMCP_MIN_VAF/_DEPTH/_MAPQ/_BASEQ/
# _ALT_COUNT) are deliberately NOT set here: the code defaults in constants.py are the single
# source of truth, and an image-level override silently disabled the MAPQ filter in prod.
ENV BAMCP_REFERENCE=""
ENV BAMCP_MAX_READS="10000"
ENV BAMCP_DEFAULT_WINDOW="500"
ENV BAMCP_TRANSPORT="sse"
ENV BAMCP_HOST="0.0.0.0"
ENV BAMCP_PORT="8000"
ENV BAMCP_AUTH_ENABLED=""

# Telemetry — opt-in at deploy time by setting these on the container instance.
# JSONL works with the base install; OTel uses the bundled opentelemetry SDK
# and also requires OTEL_EXPORTER_OTLP_ENDPOINT to be set.
ENV BAMCP_TELEMETRY_ENABLED=""
ENV BAMCP_TELEMETRY_PATH=""
ENV BAMCP_TELEMETRY_OTEL=""

# Volume for BAM/CRAM data
VOLUME ["/data"]

# HTTP transport port
EXPOSE 8000

# Health check using the dedicated script
HEALTHCHECK --interval=30s --timeout=5s --retries=3 \
    CMD python /app/docker/healthcheck.py || exit 1

# Deploy identity: the commit this image was built from, served at /__version so Greenlight's
# verify gate can prove the running container is the build CI just shipped. Declared last so a
# per-commit value never invalidates the cached layers above. Unset (local builds) → sha: null.
ARG GREENLIGHT_SHA=""
ENV GREENLIGHT_SHA=${GREENLIGHT_SHA}

ENTRYPOINT ["/app/docker/entrypoint.sh"]
CMD ["python", "-m", "bamcp"]
