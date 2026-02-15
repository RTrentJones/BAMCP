# ============================================================================
# BAMCP Production Dockerfile
# Multi-stage build: compile pysam from source, then copy into slim runtime.
# ============================================================================

# ---------------------------------------------------------------------------
# Stage 1 — Build
# ---------------------------------------------------------------------------
FROM python:3.12-slim AS builder

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

COPY pyproject.toml README.md ./
COPY src/ src/

# Install the package and its dependencies into a prefix we can copy later.
RUN pip install --no-cache-dir --prefix=/install .

# ---------------------------------------------------------------------------
# Stage 2 — Runtime
# ---------------------------------------------------------------------------
FROM python:3.12-slim AS runtime

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

# Default env vars (can be overridden)
ENV BAMCP_REFERENCE=""
ENV BAMCP_MAX_READS="10000"
ENV BAMCP_DEFAULT_WINDOW="500"
ENV BAMCP_MIN_VAF="0.1"
ENV BAMCP_MIN_DEPTH="10"
ENV BAMCP_MIN_MAPQ="0"
ENV BAMCP_TRANSPORT="sse"
ENV BAMCP_HOST="0.0.0.0"
ENV BAMCP_PORT="8000"
ENV BAMCP_AUTH_ENABLED=""

# Volume for BAM/CRAM data
VOLUME ["/data"]

# HTTP transport port
EXPOSE 8000

# Health check using the dedicated script
HEALTHCHECK --interval=30s --timeout=5s --retries=3 \
    CMD python /app/docker/healthcheck.py || exit 1

ENTRYPOINT ["/app/docker/entrypoint.sh"]
CMD ["python", "-m", "bamcp"]
