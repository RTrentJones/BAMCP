"""RegionData serialization for BAMCP tool responses.

Converts parsed BAM/CRAM region data into JSON-compatible dicts with
variant evidence enhancement, confidence scoring, and optional compaction.
"""

from __future__ import annotations

from typing import Any

from ..analysis.evidence import (
    build_mismatch_index,
    compute_artifact_risk,
    compute_confidence,
    compute_variant_evidence,
)
from .parsers import AlignedRead, RegionData


def _serialize_read(r: AlignedRead, include_sequence: bool) -> dict:
    """Serialize a read, omitting null paired-end fields to reduce payload."""
    d: dict[str, Any] = {
        "name": r.name,
        "cigar": r.cigar,
        "position": r.position,
        "end_position": r.end_position,
        "mapping_quality": r.mapping_quality,
        "is_reverse": r.is_reverse,
        "mismatches": r.mismatches,
    }
    if include_sequence:
        d["sequence"] = r.sequence
    # Only include paired-end fields if the read is actually paired
    if r.is_paired:
        d["mate_position"] = r.mate_position
        d["mate_contig"] = r.mate_contig
        d["insert_size"] = r.insert_size
        d["is_proper_pair"] = r.is_proper_pair
        d["is_read1"] = r.is_read1
        d["is_paired"] = True
    # Include soft clips if present
    if r.soft_clips:
        d["soft_clips"] = [
            {
                "position": sc.position,
                "length": sc.length,
                "sequence": sc.sequence,
                "side": sc.side,
            }
            for sc in r.soft_clips
        ]
    return d


def serialize_region_data(data: RegionData, compact: bool | None = None) -> dict:
    """Serialize RegionData to a JSON-compatible dict.

    Args:
        data: The region data to serialize.
        compact: If True, omit sequences to reduce payload size.
                 If None (default), auto-detect based on region size:
                 include sequences for regions <= 500bp (base-level view possible).
    """
    # Auto-detect compact mode based on region size
    # Base rendering requires scale >= 10, which means ~80bp at 800px width
    # Use 500bp threshold for safety margin
    if compact is None:
        region_size = data.end - data.start
        compact = region_size > 500

    # Build mismatch index once for O(1) variant evidence lookups
    mismatch_index = build_mismatch_index(data.reads)

    # Compute variant evidence using index (O(k) per variant instead of O(n*m) total)
    variant_evidence = {}
    enhanced_variants = []

    for variant in data.variants:
        key = f"{variant['position']}:{variant['ref']}>{variant['alt']}"
        evidence = compute_variant_evidence(mismatch_index, variant)

        # Compute artifact risk
        artifact_risk = compute_artifact_risk(
            variant, evidence, data.reference_sequence, data.start
        )
        evidence["artifact_risk"] = artifact_risk
        variant_evidence[key] = evidence

        # Enhance variant with evidence data for table display
        enhanced = dict(variant)
        enhanced["strand_forward"] = evidence["forward_count"]
        enhanced["strand_reverse"] = evidence["reverse_count"]
        enhanced["mean_quality"] = evidence["mean_quality"]
        enhanced["artifact_risk"] = artifact_risk

        confidence = compute_confidence(variant, evidence, artifact_risk)
        enhanced["confidence"] = confidence
        # Keep is_low_confidence for backwards compatibility
        enhanced["is_low_confidence"] = confidence == "low"

        enhanced_variants.append(enhanced)

    # Serialize reads - compact mode omits sequences for smaller payload
    reads_data = [_serialize_read(r, include_sequence=not compact) for r in data.reads]

    return {
        "contig": data.contig,
        "start": data.start,
        "end": data.end,
        "reads": reads_data,
        "coverage": data.coverage,
        "variants": enhanced_variants,
        "reference_sequence": data.reference_sequence,
        "variant_evidence": variant_evidence,
        "compact": compact,
    }
