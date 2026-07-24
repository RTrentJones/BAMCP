"""RegionData serialization for BAMCP tool responses.

Converts parsed BAM/CRAM region data into JSON-compatible dicts with
variant evidence enhancement, confidence scoring, and optional compaction.
"""

from __future__ import annotations

from typing import Any

from ..analysis.evidence import enhance_variants_with_evidence
from ..constants import PAYLOAD_SCHEMA_VERSION
from .parsers import AlignedRead, RegionData


def _serialize_reads_columnar(reads: list[AlignedRead], include_sequence: bool) -> dict[str, Any]:
    """Serialize reads as parallel arrays (columnar) to shrink the payload.

    A per-read object repeats every JSON key for every read; columnar layout stores
    each field once as a position-aligned array. The viewer decodes this back into
    individual read objects on ingestion, so nothing downstream of that decode
    changes. Optional groups (sequence/qualities, paired-end, soft clips) are
    emitted only when at least one read needs them, and are position-aligned with
    ``None`` where a given read doesn't carry that field.
    """
    cols: dict[str, Any] = {
        "count": len(reads),
        "name": [r.name for r in reads],
        "cigar": [r.cigar for r in reads],
        "position": [r.position for r in reads],
        "end_position": [r.end_position for r in reads],
        "mapping_quality": [r.mapping_quality for r in reads],
        "is_reverse": [r.is_reverse for r in reads],
        "mismatches": [r.mismatches for r in reads],
    }

    if include_sequence:
        cols["sequence"] = [r.sequence for r in reads]
        # Base qualities ride along with the sequence — only useful at base-level
        # zoom, and they enable the DeepVariant-style channels.
        cols["qualities"] = [list(r.qualities) if r.qualities else None for r in reads]

    if any(r.is_paired for r in reads):
        cols["is_paired"] = [r.is_paired for r in reads]
        cols["mate_position"] = [r.mate_position if r.is_paired else None for r in reads]
        cols["mate_contig"] = [r.mate_contig if r.is_paired else None for r in reads]
        cols["insert_size"] = [r.insert_size if r.is_paired else None for r in reads]
        cols["is_proper_pair"] = [r.is_proper_pair if r.is_paired else None for r in reads]
        cols["is_read1"] = [r.is_read1 if r.is_paired else None for r in reads]

    if any(r.soft_clips for r in reads):
        cols["soft_clips"] = [
            [
                {
                    "position": sc.position,
                    "length": sc.length,
                    "sequence": sc.sequence,
                    "side": sc.side,
                }
                for sc in r.soft_clips
            ]
            for r in reads
        ]

    return cols


def build_enhanced_variants(data: RegionData) -> tuple[list[dict], dict]:
    """Evidence-enhanced variants + the per-variant evidence map for a region.

    Shared by :func:`serialize_region_data` (viewer payload) and the evidence-fusion path
    (``classify_variant``), so the strand/quality/artifact/confidence signals reach both the
    viewer and the ACMG reasoning scaffold. Positions stay 0-based here (the internal frame).

    Returns ``(enhanced_variants, variant_evidence)``.
    """
    return enhance_variants_with_evidence(
        data.variants, data.reads, data.reference_sequence, data.start
    )


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

    # Compute read-level evidence for every variant (shared with get_variants).
    enhanced_variants, variant_evidence = build_enhanced_variants(data)

    # Serialize reads columnar - compact mode omits sequences for smaller payload
    reads_data = _serialize_reads_columnar(data.reads, include_sequence=not compact)

    return {
        "schema_version": PAYLOAD_SCHEMA_VERSION,
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
