"""Shared constants for BAMCP runtime defaults and thresholds.

This module is the single source of truth for default values that are consumed
across configuration loading, tool behavior, and cache management.
"""

from __future__ import annotations

from pathlib import Path

# Paths and networking defaults
DEFAULT_HOST = "0.0.0.0"  # noqa: S104
DEFAULT_PORT = 8000
DEFAULT_ISSUER_URL = "http://localhost:8000"
DEFAULT_RESOURCE_SERVER_URL = "http://localhost:8000"
DEFAULT_TRANSPORT = "stdio"
# Project root is 3 levels up from this file: src/bamcp/constants.py -> BAMCP/
_PROJECT_ROOT = Path(__file__).parent.parent.parent
DEFAULT_CACHE_DIR = _PROJECT_ROOT / ".cache"

# Genomics operation defaults
DEFAULT_MAX_READS = 10_000
DEFAULT_WINDOW_SIZE = 500
DEFAULT_MIN_VAF = 0.1
DEFAULT_MIN_DEPTH = 3
DEFAULT_MIN_MAPQ = 0
DEFAULT_CONTIG = "chr1"

# Auth and integrations
DEFAULT_TOKEN_EXPIRY_SECONDS = 3_600
DEFAULT_GNOMAD_DATASET = "gnomad_r4"
DEFAULT_GENOME_BUILD = "GRCh38"

# Cache behavior
DEFAULT_CACHE_TTL_SECONDS = 86_400  # 24 hours
CACHE_SESSION_ID_LENGTH = 8
REMOTE_FILE_SCHEMES = ("http://", "https://")

# Timeout for BAM/CRAM parsing operations (seconds)
BAM_PARSE_TIMEOUT_SECONDS = 30.0

# UI constants
VIEWER_RESOURCE_URI = "ui://bamcp/viewer"

# Variant quality heuristics
LOW_CONFIDENCE_MIN_DEPTH = 10
LOW_CONFIDENCE_MIN_VAF = 0.10
LOW_CONFIDENCE_MIN_MEAN_QUALITY = 20
LOW_CONFIDENCE_MAX_STRAND_BIAS = 0.9

QUALITY_HISTOGRAM_BINS = [0, 10, 20, 30, 40]
POSITION_HISTOGRAM_BINS = [0, 25, 50, 75, 100, 150]
MAPQ_HISTOGRAM_BINS = [0, 10, 20, 30, 40, 50, 60]

# Artifact detection thresholds
ARTIFACT_STRAND_BIAS_THRESHOLD = 0.8
ARTIFACT_NEAR_END_FRACTION_THRESHOLD = 0.5
ARTIFACT_LOW_MAPQ_FRACTION_THRESHOLD = 0.3
ARTIFACT_HOMOPOLYMER_LENGTH_THRESHOLD = 4

# Display mode configurations for IGV-like read visualization
DISPLAY_MODES = {
    "squished": {"read_height": 6, "read_gap": 1, "show_labels": False},
    "compact": {"read_height": 12, "read_gap": 2, "show_labels": False},
    "expanded": {"read_height": 24, "read_gap": 4, "show_labels": True},
}
DEFAULT_DISPLAY_MODE = "compact"
