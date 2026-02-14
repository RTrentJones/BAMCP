"""Core BAM/CRAM data pipeline modules."""

from .cache import BAMIndexCache
from .parsers import AlignedRead, RegionData, SoftClip, fetch_region
from .reference import detect_genome_build, get_public_reference_url
from .serialization import serialize_region_data
from .tools import (
    get_cache,
    handle_get_coverage,
    handle_get_region_summary,
    handle_get_variants,
    handle_jump_to,
    handle_list_contigs,
    handle_lookup_clinvar,
    handle_lookup_gnomad,
    handle_visualize_region,
)
from .validation import (
    validate_path,
    validate_region,
    validate_remote_url,
    validate_variant_input,
)

__all__ = [
    "AlignedRead",
    "BAMIndexCache",
    "RegionData",
    "SoftClip",
    "detect_genome_build",
    "fetch_region",
    "get_cache",
    "get_public_reference_url",
    "handle_get_coverage",
    "handle_get_region_summary",
    "handle_get_variants",
    "handle_jump_to",
    "handle_list_contigs",
    "handle_lookup_clinvar",
    "handle_lookup_gnomad",
    "handle_visualize_region",
    "serialize_region_data",
    "validate_path",
    "validate_region",
    "validate_remote_url",
    "validate_variant_input",
]
