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
# XDG-compliant cache directory: ~/.cache/bamcp
DEFAULT_CACHE_DIR = Path.home() / ".cache" / "bamcp"

# Genomics operation defaults
DEFAULT_MAX_READS = 10_000
DEFAULT_WINDOW_SIZE = 500
DEFAULT_MIN_VAF = 0.02
DEFAULT_MIN_DEPTH = 2
DEFAULT_MIN_MAPQ = 0
DEFAULT_MIN_BASEQ = 0

# Candidate indel detection runs a pysam pileup pass (costlier than count_coverage).
# Skip it above this region size so large visualize_region requests stay responsive;
# SNV/coverage detection (C-backed count_coverage) still runs at any size.
INDEL_DETECTION_MAX_REGION = 100_000

# Readless coverage/variant queries are computed and cached per fixed-width genomic
# tile so panning the viewport reuses work instead of recomputing each region.
REGION_TILE_SIZE = 4_096

# Version of the RegionData payload contract sent to the viewer (_meta.ui/init). Bump on any
# breaking shape change so the viewer (and a schema-compat test) can detect a mismatch instead
# of silently mis-decoding. A first step toward a generated, typed server<->viewer contract.
PAYLOAD_SCHEMA_VERSION = 1

# Cap on distinct regions held in a server's in-memory region cache (LRU-evicted).
# Bounds memory: each entry can hold a full RegionData payload.
REGION_CACHE_MAX_ENTRIES = 128
# ...and by approximate in-memory bytes (see tools._approx_region_bytes): a full-mode 1 Mb entry
# can weigh tens of MB, so the count cap alone could exceed a 2 GB container. An entry larger
# than a quarter of this budget is served uncached.
REGION_CACHE_MAX_BYTES = 512 * 1024 * 1024

# Max concurrent pysam parses. All blocking BAM/CRAM/VCF work runs in a dedicated thread pool
# of this size, so a burst of requests (or parses that hang past their timeout) is confined to
# this pool and cannot starve the event loop's default executor.
PARSE_MAX_CONCURRENCY = 8

# Auth and integrations
DEFAULT_TOKEN_EXPIRY_SECONDS = 3_600
# When auth is enabled with no explicit BAMCP_REQUIRED_SCOPES, require this scope so tokens
# are never scope-less. BAMCP is a read-only evidence viewer, hence a single read scope.
DEFAULT_REQUIRED_SCOPES = ("bamcp:read",)
DEFAULT_GNOMAD_DATASET = "gnomad_r4"
DEFAULT_GENOME_BUILD = "GRCh38"

# Cache behavior
DEFAULT_CACHE_TTL_SECONDS = 86_400  # 24 hours
CACHE_SESSION_ID_LENGTH = 8
REMOTE_FILE_SCHEMES = ("http://", "https://")
DEFAULT_MAX_REMOTE_INDEX_BYTES = 64 * 1024 * 1024
# Total on-disk budget for one session's cached remote indexes. Each distinct URL (and a new
# query string is a new URL) can add up to DEFAULT_MAX_REMOTE_INDEX_BYTES, so without a cap an
# open server's cache could be grown without bound; oldest files are evicted past this.
DEFAULT_MAX_INDEX_CACHE_BYTES = 512 * 1024 * 1024

# Timeout for BAM/CRAM parsing operations (seconds)
BAM_PARSE_TIMEOUT_SECONDS = 30.0

# Whole-contig variant scanning
SCAN_VARIANTS_TIMEOUT_SECONDS = 120.0
SCAN_VARIANTS_CHUNK_SIZE = 50_000
SCAN_VARIANTS_MAX_REGION = 250_000_000
SCAN_VARIANTS_MAX_RESULTS = 500

# UI constants
VIEWER_RESOURCE_URI = "ui://bamcp/viewer"

# Variant quality heuristics
LOW_CONFIDENCE_MIN_DEPTH = 10
LOW_CONFIDENCE_MIN_VAF = 0.10
LOW_CONFIDENCE_MIN_MEAN_QUALITY = 20
LOW_CONFIDENCE_MAX_STRAND_BIAS = 0.9

# Thresholds for a confident call, used by the ACMG quality assessment (classify_variant).
# Also the "high" tier of compute_confidence and the curation recommendation/balance cutoffs.
HIGH_CONFIDENCE_MIN_DEPTH = 20
HIGH_CONFIDENCE_MIN_MEAN_QUALITY = 25
HIGH_CONFIDENCE_MIN_VAF = 0.2
HIGH_CONFIDENCE_MAX_STRAND_BIAS = 0.7

QUALITY_HISTOGRAM_BINS = [0, 10, 20, 30, 40]
POSITION_HISTOGRAM_BINS = [0, 25, 50, 75, 100, 150]
MAPQ_HISTOGRAM_BINS = [0, 10, 20, 30, 40, 50, 60]

# Artifact detection thresholds
ARTIFACT_STRAND_BIAS_THRESHOLD = 0.8
# Read-position bias: a supporting base is "near an end" when it lies within this many
# ALIGNED (non-soft-clipped) bases of the nearest aligned read end. An absolute ~10bp
# cutoff (not the display histogram's 0-50bp bins, which cover 2/3 of a 150bp read and
# all of a 100bp read) keeps uniformly positioned support on clean 2x150 data unflagged.
ARTIFACT_NEAR_END_DISTANCE_BP = 10
ARTIFACT_NEAR_END_FRACTION_THRESHOLD = 0.5
# Short reads put most positions near an end even under uniform placement, so the
# observed near-end fraction must also exceed the uniform expectation by this margin.
ARTIFACT_NEAR_END_EXCESS_OVER_EXPECTED = 0.25
ARTIFACT_LOW_MAPQ_FRACTION_THRESHOLD = 0.3
ARTIFACT_HOMOPOLYMER_LENGTH_THRESHOLD = 4
