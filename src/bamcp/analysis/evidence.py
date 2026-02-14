"""Variant evidence computation and artifact risk assessment.

Pure functions for analyzing variant support quality: mismatch indexing,
strand/position/quality evidence extraction, artifact risk scoring,
and confidence level computation.
"""

from __future__ import annotations

import re

from ..constants import (
    ARTIFACT_HOMOPOLYMER_LENGTH_THRESHOLD,
    ARTIFACT_LOW_MAPQ_FRACTION_THRESHOLD,
    ARTIFACT_NEAR_END_FRACTION_THRESHOLD,
    ARTIFACT_STRAND_BIAS_THRESHOLD,
    LOW_CONFIDENCE_MAX_STRAND_BIAS,
    LOW_CONFIDENCE_MIN_DEPTH,
    LOW_CONFIDENCE_MIN_MEAN_QUALITY,
    LOW_CONFIDENCE_MIN_VAF,
    MAPQ_HISTOGRAM_BINS,
    POSITION_HISTOGRAM_BINS,
    QUALITY_HISTOGRAM_BINS,
)
from ..core.parsers import AlignedRead


def get_query_position(read: AlignedRead, ref_pos: int) -> int | None:
    """Get the query (read) position corresponding to a reference position.

    Uses CIGAR string to map reference coordinates to read coordinates.
    Returns None if the reference position is not aligned to the read.
    """
    cigar_ops = re.findall(r"(\d+)([MIDNSHP=X])", read.cigar)
    query_pos = 0
    ref_cursor = read.position

    for length_str, op in cigar_ops:
        length = int(length_str)

        if op in ("M", "=", "X"):
            # Match/mismatch: consumes both query and reference
            if ref_cursor <= ref_pos < ref_cursor + length:
                return query_pos + (ref_pos - ref_cursor)
            query_pos += length
            ref_cursor += length
        elif op == "I":
            # Insertion: consumes query only
            query_pos += length
        elif op in ("D", "N"):
            # Deletion/skip: consumes reference only
            if ref_cursor <= ref_pos < ref_cursor + length:
                return None  # Position is in a deletion
            ref_cursor += length
        elif op == "S":
            # Soft clip: consumes query only
            query_pos += length
        elif op == "H":
            # Hard clip: consumes nothing
            pass

    return None


def bin_values(values: list[int], bins: list[int]) -> list[int]:
    """Bin values into histogram buckets.

    Args:
        values: List of values to bin.
        bins: List of bin boundaries (e.g., [0, 10, 20, 30, 40]).
              Values >= bins[i] and < bins[i+1] go into bin i.
              The last bin captures all values >= bins[-1].

    Returns:
        List of counts per bin, length = len(bins).
    """
    counts = [0] * len(bins)
    for v in values:
        for i in range(len(bins) - 1, -1, -1):
            if v >= bins[i]:
                counts[i] += 1
                break
    return counts


def detect_homopolymer(seq: str, pos: int) -> int:
    """Detect homopolymer run length at a given position.

    Args:
        seq: Reference sequence.
        pos: Position within the sequence (0-indexed).

    Returns:
        Length of homopolymer run containing the position.
    """
    if not seq or pos < 0 or pos >= len(seq):
        return 0

    base = seq[pos].upper()
    if base not in "ACGT":
        return 0

    length = 1

    # Extend left
    i = pos - 1
    while i >= 0 and seq[i].upper() == base:
        length += 1
        i -= 1

    # Extend right
    i = pos + 1
    while i < len(seq) and seq[i].upper() == base:
        length += 1
        i += 1

    return length


def build_mismatch_index(
    reads: list[AlignedRead],
) -> dict[tuple[int, str], list[tuple[AlignedRead, dict]]]:
    """Build an index of mismatches by (position, alt) for O(1) variant lookup.

    Returns:
        Dict mapping (position, alt) to list of (read, mismatch) tuples.
    """
    index: dict[tuple[int, str], list[tuple[AlignedRead, dict]]] = {}
    for read in reads:
        for mm in read.mismatches:
            key = (mm["pos"], mm["alt"])
            if key not in index:
                index[key] = []
            index[key].append((read, mm))
    return index


def compute_variant_evidence(
    index: dict[tuple[int, str], list[tuple[AlignedRead, dict]]],
    variant: dict,
) -> dict:
    """Compute variant evidence using pre-built mismatch index.

    This is O(k) where k is the number of reads supporting the variant,
    rather than O(n*m) when iterating all reads for each variant.

    Returns evidence including histogram distributions for curation.
    """
    pos = variant["position"]
    alt = variant["alt"]

    matches = index.get((pos, alt), [])

    forward_count = 0
    reverse_count = 0
    qualities: list[int] = []
    positions_in_read: list[int] = []
    mapq_values: list[int] = []

    for read, mm in matches:
        if read.is_reverse:
            reverse_count += 1
        else:
            forward_count += 1

        # Get quality at mismatch position
        query_pos = get_query_position(read, mm["pos"])
        if query_pos is not None and query_pos < len(read.qualities):
            qualities.append(read.qualities[query_pos])

            # Compute position relative to nearest read end
            read_length = len(read.qualities)
            dist_from_end = min(query_pos, read_length - query_pos - 1)
            positions_in_read.append(dist_from_end)

        # Collect MAPQ for all reads
        mapq_values.append(read.mapping_quality)

    # Strand bias: 0 = balanced, 1 = all one strand
    total = forward_count + reverse_count
    strand_bias = abs(forward_count - reverse_count) / max(total, 1)

    # Compute histogram distributions
    quality_histogram = bin_values(qualities, QUALITY_HISTOGRAM_BINS)
    position_histogram = bin_values(positions_in_read, POSITION_HISTOGRAM_BINS)
    mapq_histogram = bin_values(mapq_values, MAPQ_HISTOGRAM_BINS)

    return {
        "forward_count": forward_count,
        "reverse_count": reverse_count,
        "strand_bias": round(strand_bias, 3),
        "mean_quality": round(sum(qualities) / len(qualities), 1) if qualities else 0,
        "median_quality": sorted(qualities)[len(qualities) // 2] if qualities else 0,
        "quality_histogram": quality_histogram,
        "position_histogram": position_histogram,
        "mapq_histogram": mapq_histogram,
    }


def compute_artifact_risk(
    variant: dict,
    evidence: dict,
    reference_sequence: str | None,
    region_start: int,
) -> dict:
    """Compute artifact risk indicators for variant curation.

    Args:
        variant: Variant dict with position, ref, alt, depth, vaf.
        evidence: Evidence dict from compute_variant_evidence.
        reference_sequence: Reference sequence string (may be None).
        region_start: Start position of the region.

    Returns:
        Dict with risks list, risk_score, and artifact_likelihood.
    """
    risks: list[dict] = []
    risk_score = 0.0

    # 1. Position-in-read bias (ReadPosRankSum equivalent)
    pos_hist = evidence.get("position_histogram", [])
    if pos_hist:
        total = sum(pos_hist)
        if total > 0:
            # First two bins (0-25, 25-50) represent near-end positions
            near_end = sum(pos_hist[:2])
            near_end_fraction = near_end / total
            if near_end_fraction > ARTIFACT_NEAR_END_FRACTION_THRESHOLD:
                risks.append(
                    {
                        "type": "read_position_bias",
                        "severity": "medium",
                        "description": f"{near_end_fraction:.0%} of bases near read ends",
                        "value": round(near_end_fraction, 3),
                    }
                )
                risk_score += 0.3

    # 2. Strand bias
    strand_bias = evidence.get("strand_bias", 0)
    if strand_bias > ARTIFACT_STRAND_BIAS_THRESHOLD:
        severity = "high" if strand_bias > 0.95 else "medium"
        risks.append(
            {
                "type": "strand_bias",
                "severity": severity,
                "description": f"Strand bias: {strand_bias:.0%}",
                "value": round(strand_bias, 3),
            }
        )
        risk_score += 0.4 if severity == "high" else 0.2

    # 3. Low MAPQ fraction
    mapq_hist = evidence.get("mapq_histogram", [])
    if mapq_hist:
        total = sum(mapq_hist)
        if total > 0:
            # First three bins (0-10, 10-20, 20-30) represent low MAPQ
            low_mapq = sum(mapq_hist[:3])
            low_mapq_fraction = low_mapq / total
            if low_mapq_fraction > ARTIFACT_LOW_MAPQ_FRACTION_THRESHOLD:
                risks.append(
                    {
                        "type": "low_mapq",
                        "severity": "medium",
                        "description": f"{low_mapq_fraction:.0%} reads with MAPQ < 30",
                        "value": round(low_mapq_fraction, 3),
                    }
                )
                risk_score += 0.25

    # 4. Homopolymer context
    if reference_sequence:
        pos_in_region = variant["position"] - region_start
        if 0 <= pos_in_region < len(reference_sequence):
            hp_len = detect_homopolymer(reference_sequence, pos_in_region)
            if hp_len >= ARTIFACT_HOMOPOLYMER_LENGTH_THRESHOLD:
                severity = "high" if hp_len >= 6 else "medium"
                risks.append(
                    {
                        "type": "homopolymer",
                        "severity": severity,
                        "description": f"In {hp_len}bp homopolymer run",
                        "value": hp_len,
                    }
                )
                risk_score += 0.3 if severity == "high" else 0.15

    # 5. Low depth
    depth = variant.get("depth", 0)
    if depth < LOW_CONFIDENCE_MIN_DEPTH:
        severity = "high" if depth < 5 else "medium"
        risks.append(
            {
                "type": "low_depth",
                "severity": severity,
                "description": f"Low coverage depth ({depth}x)",
                "value": depth,
            }
        )
        risk_score += 0.3 if severity == "high" else 0.15

    # Cap risk score at 1.0
    risk_score = min(risk_score, 1.0)

    # Determine artifact likelihood
    if risk_score >= 0.6:
        artifact_likelihood = "high"
    elif risk_score >= 0.3:
        artifact_likelihood = "medium"
    else:
        artifact_likelihood = "low"

    return {
        "risks": risks,
        "risk_score": round(risk_score, 2),
        "artifact_likelihood": artifact_likelihood,
    }


def compute_confidence(
    variant: dict,
    evidence: dict,
    artifact_risk: dict,
) -> str:
    """Compute confidence level for a variant call.

    Evaluates VAF, depth, quality, and strand bias against thresholds,
    then downgrades if artifact risk is elevated.

    Args:
        variant: Variant dict with vaf and depth.
        evidence: Evidence dict with mean_quality and strand_bias.
        artifact_risk: Artifact risk dict with artifact_likelihood.

    Returns:
        Confidence level: "high", "medium", or "low".
    """
    vaf = variant["vaf"]
    depth = variant["depth"]
    quality = evidence["mean_quality"]
    strand_bias = evidence["strand_bias"]

    if vaf >= 0.2 and depth >= 20 and quality >= 25 and strand_bias <= 0.7:
        confidence = "high"
    elif (
        vaf >= LOW_CONFIDENCE_MIN_VAF
        and depth >= LOW_CONFIDENCE_MIN_DEPTH
        and quality >= LOW_CONFIDENCE_MIN_MEAN_QUALITY
        and strand_bias <= LOW_CONFIDENCE_MAX_STRAND_BIAS
    ):
        confidence = "medium"
    else:
        confidence = "low"

    # Downgrade confidence if artifact risk is high
    if artifact_risk["artifact_likelihood"] == "high":
        confidence = "low"
    elif artifact_risk["artifact_likelihood"] == "medium" and confidence == "high":
        confidence = "medium"

    return confidence
