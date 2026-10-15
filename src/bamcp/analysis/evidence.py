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
    ARTIFACT_NEAR_END_DISTANCE_BP,
    ARTIFACT_NEAR_END_EXCESS_OVER_EXPECTED,
    ARTIFACT_NEAR_END_FRACTION_THRESHOLD,
    ARTIFACT_STRAND_BIAS_THRESHOLD,
    HIGH_CONFIDENCE_MAX_STRAND_BIAS,
    HIGH_CONFIDENCE_MIN_DEPTH,
    HIGH_CONFIDENCE_MIN_MEAN_QUALITY,
    HIGH_CONFIDENCE_MIN_VAF,
    LOW_CONFIDENCE_MAX_STRAND_BIAS,
    LOW_CONFIDENCE_MIN_DEPTH,
    LOW_CONFIDENCE_MIN_MEAN_QUALITY,
    LOW_CONFIDENCE_MIN_VAF,
    MAPQ_HISTOGRAM_BINS,
    POSITION_HISTOGRAM_BINS,
    QUALITY_HISTOGRAM_BINS,
)
from ..core.parsers import AlignedRead

# Mean alt-base quality below this is a low_base_quality concern (Phred 20 = 1% error).
ARTIFACT_LOW_BASEQ_THRESHOLD = LOW_CONFIDENCE_MIN_MEAN_QUALITY


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


def homopolymer_context_length(variant: dict, seq: str, pos_in_region: int) -> int:
    """Homopolymer run length relevant to a variant at ``pos_in_region`` of ``seq``.

    SNV: the run containing the base. Indel: the longest run touching the anchor, the
    deleted/inserted span, or the base just after it (``pos+1 .. pos+len+1``).
    """
    if not is_indel(variant):
        return detect_homopolymer(seq, pos_in_region)
    span = max(len(str(variant.get("ref", ""))), len(str(variant.get("alt", ""))))
    stop = min(len(seq), pos_in_region + span + 1)
    return max(
        (detect_homopolymer(seq, i) for i in range(pos_in_region, stop)),
        default=0,
    )


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


_CIGAR_RE = re.compile(r"(\d+)([MIDNSHP=X])")


def _aligned_query_span(read: AlignedRead) -> tuple[int, int]:
    """Return ``(first, last)`` query indices of the ALIGNED (non-soft-clipped) part.

    Soft-clipped bases are in the query but were not aligned, so a base 3bp from the
    end of the aligned block is "near the end" even when 20 clipped bases follow it.
    """
    ops = [(int(n), op) for n, op in _CIGAR_RE.findall(read.cigar)]
    qlen = sum(n for n, op in ops if op in "MIS=X")
    if qlen == 0:
        qlen = len(read.qualities) or len(read.sequence)
    lead = 0
    for n, op in ops:
        if op == "H":
            continue
        if op == "S":
            lead += n
            continue
        break
    trail = 0
    for n, op in reversed(ops):
        if op == "H":
            continue
        if op == "S":
            trail += n
            continue
        break
    return lead, qlen - trail - 1


def summarize_read_positions(
    distances: list[int], aligned_lengths: list[int] | None = None
) -> dict:
    """Near-end summary for supporting bases.

    ``distances`` are each supporting base's distance (bp) to the nearest aligned read
    end. A base is "near the end" when that distance is below
    ``ARTIFACT_NEAR_END_DISTANCE_BP``. When the aligned read lengths are known, the
    fraction expected under UNIFORM placement (``min(L, 2c) / L`` per read, averaged)
    is reported too, so short reads — where most positions are near an end anyway —
    are not mistaken for a position bias.
    """
    cutoff = ARTIFACT_NEAR_END_DISTANCE_BP
    n = len(distances)
    near = sum(1 for d in distances if d < cutoff)
    expected: float | None = None
    if aligned_lengths:
        per_read = [min(length, 2 * cutoff) / length for length in aligned_lengths if length > 0]
        if per_read:
            expected = round(sum(per_read) / len(per_read), 3)
    return {
        "near_end_cutoff_bp": cutoff,
        "near_end_count": near,
        "position_count": n,
        "near_end_fraction": round(near / n, 3) if n else 0.0,
        "near_end_expected_fraction": expected,
    }


def is_read_position_biased(evidence: dict) -> bool:
    """True when supporting bases cluster near aligned read ends beyond chance.

    Requires both an absolute majority near an end (> ``ARTIFACT_NEAR_END_FRACTION_THRESHOLD``)
    and, when the uniform expectation is known, an excess over it of at least
    ``ARTIFACT_NEAR_END_EXCESS_OVER_EXPECTED``. Uniform positions never flag.
    """
    if not evidence.get("position_count"):
        return False
    frac = float(evidence.get("near_end_fraction") or 0.0)
    if frac <= ARTIFACT_NEAR_END_FRACTION_THRESHOLD:
        return False
    expected = evidence.get("near_end_expected_fraction")
    return expected is None or frac >= expected + ARTIFACT_NEAR_END_EXCESS_OVER_EXPECTED


def _not_computed_evidence(reason: str) -> dict:
    """Evidence record for a variant whose supporting reads could not be located.

    Numeric fields stay numeric (the viewer renders them) but ``read_evidence_computed``
    is False and ``read_evidence_note`` says why; analysis consumers must treat strand,
    quality and position as UNKNOWN rather than as balanced/clean zeros.
    """
    return {
        "read_evidence_computed": False,
        "read_evidence_note": reason,
        "forward_count": 0,
        "reverse_count": 0,
        "strand_bias": 0.0,
        "mean_quality": 0,
        "median_quality": 0,
        "quality_histogram": bin_values([], QUALITY_HISTOGRAM_BINS),
        "position_histogram": bin_values([], POSITION_HISTOGRAM_BINS),
        "mapq_histogram": bin_values([], MAPQ_HISTOGRAM_BINS),
        **summarize_read_positions([]),
    }


def _aggregate_read_evidence(supporting: list[tuple[AlignedRead, int | None, int | None]]) -> dict:
    """Aggregate ``(read, query_position, base_quality)`` support into an evidence dict.

    ``query_position`` may be None when the site isn't aligned in a read; ``base_quality``
    is the quality attributed to the supporting event (the base itself for an SNV, the
    inserted/flanking bases for an indel), or None if unknown.
    """
    forward_count = 0
    reverse_count = 0
    qualities: list[int] = []
    positions_in_read: list[int] = []
    aligned_lengths: list[int] = []
    mapq_values: list[int] = []

    for read, query_pos, base_quality in supporting:
        if read.is_reverse:
            reverse_count += 1
        else:
            forward_count += 1

        if base_quality is not None:
            qualities.append(base_quality)
        if query_pos is not None:
            first, last = _aligned_query_span(read)
            if first <= query_pos <= last:
                # Distance to the nearest ALIGNED read end (soft clips excluded).
                positions_in_read.append(min(query_pos - first, last - query_pos))
                aligned_lengths.append(last - first + 1)

        mapq_values.append(read.mapping_quality)

    # Strand bias: 0 = balanced, 1 = all one strand
    total = forward_count + reverse_count
    strand_bias = abs(forward_count - reverse_count) / max(total, 1)

    return {
        "read_evidence_computed": True,
        "forward_count": forward_count,
        "reverse_count": reverse_count,
        "strand_bias": round(strand_bias, 3),
        "mean_quality": round(sum(qualities) / len(qualities), 1) if qualities else 0,
        "median_quality": sorted(qualities)[len(qualities) // 2] if qualities else 0,
        "quality_histogram": bin_values(qualities, QUALITY_HISTOGRAM_BINS),
        "position_histogram": bin_values(positions_in_read, POSITION_HISTOGRAM_BINS),
        "mapq_histogram": bin_values(mapq_values, MAPQ_HISTOGRAM_BINS),
        **summarize_read_positions(positions_in_read, aligned_lengths),
    }


def _indel_key(variant: dict) -> tuple[int, str, int | str] | None:
    """Map a VCF-style anchored indel to ``(anchor_pos, kind, length|inserted_seq)``.

    Returns None for SNVs and for complex/unanchored alleles BAMCP can't key from CIGAR.
    """
    ref = str(variant.get("ref", "")).upper()
    alt = str(variant.get("alt", "")).upper()
    if len(ref) == len(alt) or not ref or not alt or ref[0] != alt[0]:
        return None
    if len(ref) > len(alt):
        if len(alt) != 1:
            return None
        return (variant["position"], "del", len(ref) - 1)
    if len(ref) != 1:
        return None
    return (variant["position"], "ins", alt[1:])


def build_indel_index(
    reads: list[AlignedRead],
) -> dict[tuple[int, str, int | str], list[tuple[AlignedRead, int, int | None]]]:
    """Index CIGAR insertions/deletions by ``(anchor_pos, kind, length|seq)``.

    ``read.mismatches`` carries SNVs only, so indel candidates need their own index. The
    key matches :func:`_indel_key` (VCF-style anchor = last reference base before the
    event, 0-based — the same convention ``detect_indels`` reports via pileup). Each entry
    is ``(read, query_pos, event_quality)``: for an insertion the query position of its
    first inserted base and the minimum inserted-base quality; for a deletion the anchor
    base's query position and the minimum quality of the two flanking bases.
    """
    index: dict[tuple[int, str, int | str], list[tuple[AlignedRead, int, int | None]]] = {}
    for read in reads:
        quals = read.qualities
        seq = read.sequence.upper()
        qpos = 0
        rpos = read.position
        prev_aligned = False  # whether the previous op ended on an aligned base
        for n_str, op in _CIGAR_RE.findall(read.cigar):
            n = int(n_str)
            if op in ("M", "=", "X"):
                qpos += n
                rpos += n
                prev_aligned = True
            elif op == "I":
                if prev_aligned and qpos + n <= len(seq):
                    ins = seq[qpos : qpos + n]
                    q = min(quals[qpos : qpos + n]) if len(quals) >= qpos + n else None
                    index.setdefault((rpos - 1, "ins", ins), []).append((read, qpos, q))
                qpos += n
            elif op == "D":
                if prev_aligned:
                    flank = [quals[i] for i in (qpos - 1, qpos) if 0 <= i < len(quals)]
                    q = min(flank) if flank else None
                    index.setdefault((rpos - 1, "del", n), []).append((read, qpos - 1, q))
                rpos += n
            elif op == "N":
                rpos += n
                prev_aligned = False
            elif op == "S":
                qpos += n
                prev_aligned = False
            # H/P consume nothing
    return index


def compute_variant_evidence(
    index: dict[tuple[int, str], list[tuple[AlignedRead, dict]]],
    variant: dict,
) -> dict:
    """Compute variant evidence using a pre-built mismatch index.

    This is O(k) where k is the number of reads supporting the variant, rather
    than O(n*m) when iterating all reads for each variant. Returns evidence
    including histogram distributions for curation.
    """
    matches = index.get((variant["position"], variant["alt"]), [])
    supporting: list[tuple[AlignedRead, int | None, int | None]] = []
    for read, mm in matches:
        qpos = get_query_position(read, mm["pos"])
        qual = read.qualities[qpos] if qpos is not None and qpos < len(read.qualities) else None
        supporting.append((read, qpos, qual))
    if not supporting and variant.get("alt_count", 0):
        return _not_computed_evidence(
            "No alt-supporting reads were found in the fetched read set (it may be "
            "capped by max_reads); strand/quality/position evidence was not computed."
        )
    return _aggregate_read_evidence(supporting)


def compute_indel_evidence(
    indel_index: dict[tuple[int, str, int | str], list[tuple[AlignedRead, int, int | None]]],
    variant: dict,
) -> dict:
    """Read-level evidence for an anchored indel, from CIGAR-derived support.

    Falls back to a not-computed record (explicit, never zero-filled "clean") when the
    allele can't be keyed or no supporting reads are in the fetched read set.
    """
    key = _indel_key(variant)
    if key is None:
        return _not_computed_evidence(
            "Complex/unanchored indel allele: read-level support was not computed."
        )
    matches = indel_index.get(key, [])
    if not matches:
        return _not_computed_evidence(
            "No indel-supporting reads were found in the fetched read set (it may be "
            "capped by max_reads); strand/quality/position evidence was not computed."
        )
    return _aggregate_read_evidence(list(matches))


def is_indel(variant: dict) -> bool:
    """True for a length-changing allele (insertion/deletion)."""
    return len(str(variant.get("ref", ""))) != len(str(variant.get("alt", "")))


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

    # 0. Read evidence not computed (e.g. supporting reads not in the fetched set). The
    # strand/quality/position zeros are placeholders, NOT a clean bill of health.
    computed = evidence.get("read_evidence_computed", True)
    if not computed:
        risks.append(
            {
                "type": "read_evidence_unavailable",
                "severity": "medium",
                "description": evidence.get("read_evidence_note")
                or "Read-level evidence was not computed for this variant",
                "value": None,
            }
        )
        risk_score += 0.3

    # 1. Position-in-read bias: supporting bases clustered within a few bp of an aligned
    # read end, beyond what uniform placement would give (see is_read_position_biased).
    if computed and is_read_position_biased(evidence):
        near_end_fraction = float(evidence.get("near_end_fraction") or 0.0)
        cutoff = evidence.get("near_end_cutoff_bp", ARTIFACT_NEAR_END_DISTANCE_BP)
        risks.append(
            {
                "type": "read_position_bias",
                "severity": "medium",
                "description": (
                    f"{near_end_fraction:.0%} of supporting bases within {cutoff}bp "
                    "of an aligned read end"
                ),
                "value": round(near_end_fraction, 3),
            }
        )
        risk_score += 0.3

    # 2. Strand bias
    strand_bias = evidence.get("strand_bias", 0)
    if computed and strand_bias > ARTIFACT_STRAND_BIAS_THRESHOLD:
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

    # 2b. Low base quality: the alt bases themselves are low-confidence calls. Only when
    # qualities were actually observed (an empty histogram means "unknown", not "Q0").
    quality_observed = sum(evidence.get("quality_histogram") or []) > 0
    mean_quality = float(evidence.get("mean_quality") or 0)
    low_baseq = computed and quality_observed and mean_quality < ARTIFACT_LOW_BASEQ_THRESHOLD
    if low_baseq:
        risks.append(
            {
                "type": "low_base_quality",
                "severity": "medium",
                "description": (
                    f"Low mean alt base quality (Q{mean_quality:.1f} < "
                    f"Q{ARTIFACT_LOW_BASEQ_THRESHOLD})"
                ),
                "value": round(mean_quality, 1),
            }
        )
        risk_score += 0.3

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

    # 4. Homopolymer context. For an indel the VCF anchor base precedes the event, so the
    # anchor alone under-reads the run (CAAAAAAG: anchor C → 1, real A-run → 6); take the
    # longest run over the anchor plus the affected/flanking bases instead.
    if reference_sequence:
        pos_in_region = variant["position"] - region_start
        if 0 <= pos_in_region < len(reference_sequence):
            hp_len = homopolymer_context_length(variant, reference_sequence, pos_in_region)
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

    # Low-quality alt bases that all sit on one strand is the classic sequencing-error
    # signature (e.g. Q8 support, 11F/0R): never less than HIGH, whatever else scored.
    if low_baseq and any(r["type"] == "strand_bias" for r in risks):
        risk_score = max(risk_score, 0.6)

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

    if not evidence.get("read_evidence_computed", True):
        # Strand/quality are unknown, not clean — never let placeholders earn confidence.
        return "low"

    if (
        vaf >= HIGH_CONFIDENCE_MIN_VAF
        and depth >= HIGH_CONFIDENCE_MIN_DEPTH
        and quality >= HIGH_CONFIDENCE_MIN_MEAN_QUALITY
        and strand_bias <= HIGH_CONFIDENCE_MAX_STRAND_BIAS
    ):
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


# Transient per-read arrays the annotator attaches for aggregation; not serialized.
_VCF_RAW_FIELDS = ("_alt_qualities", "_alt_positions", "_alt_mapqs")


def _vcf_evidence(
    variant: dict, reference_sequence: str | None, region_start: int
) -> tuple[dict, dict]:
    """Build a full evidence record for a VCF SNV from the annotator's per-read metrics.

    :func:`bamcp.core.parsers.annotate_vcf_snv_support` pileups the site and attaches
    strand counts, read depth, and raw quality/position/MAPQ arrays for the
    alt-supporting reads. Aggregate them into the same evidence shape as BAMCP
    candidates — strand bias, quality/position/MAPQ histograms, mean quality — then
    score artifact risk and confidence, so VCF-primary SNVs get the full read-level
    curation contract rather than zero-filled quality.
    """
    fwd = int(variant.get("strand_forward", 0))
    rev = int(variant.get("strand_reverse", 0))
    depth = int(variant.get("read_depth", 0))
    vaf = float(variant.get("read_support_vaf", 0.0))
    quals = variant.get("_alt_qualities", [])
    positions = variant.get("_alt_positions", [])
    mapqs = variant.get("_alt_mapqs", [])
    strand_bias = round(abs(fwd - rev) / max(fwd + rev, 1), 3)

    evidence = {
        "forward_count": fwd,
        "reverse_count": rev,
        "strand_bias": strand_bias,
        "mean_quality": round(sum(quals) / len(quals), 1) if quals else 0,
        "median_quality": sorted(quals)[len(quals) // 2] if quals else 0,
        "quality_histogram": bin_values(quals, QUALITY_HISTOGRAM_BINS),
        "position_histogram": bin_values(positions, POSITION_HISTOGRAM_BINS),
        "mapq_histogram": bin_values(mapqs, MAPQ_HISTOGRAM_BINS),
        # The annotator reports distance to the nearest read end without the read length,
        # so the uniform expectation is unknown; the absolute cutoff still keeps uniformly
        # placed support on >=40bp reads below the majority threshold.
        **summarize_read_positions(list(positions)),
    }
    artifact_risk = compute_artifact_risk(
        {**variant, "depth": depth}, evidence, reference_sequence, region_start
    )
    evidence["artifact_risk"] = artifact_risk
    # Score confidence from BAMCP's measured read support (depth/VAF/quality/strand).
    confidence = compute_confidence(
        {**variant, "vaf": vaf, "depth": depth}, evidence, artifact_risk
    )

    enhanced = {k: val for k, val in variant.items() if k not in _VCF_RAW_FIELDS}
    enhanced["strand_forward"] = fwd
    enhanced["strand_reverse"] = rev
    enhanced["mean_quality"] = evidence["mean_quality"]
    enhanced["strand_bias"] = strand_bias
    enhanced["artifact_risk"] = artifact_risk
    enhanced["confidence"] = confidence
    enhanced["is_low_confidence"] = confidence == "low"
    return enhanced, evidence


def enhance_vcf_variants(
    variants: list[dict], reference_sequence: str | None, region_start: int
) -> list[dict]:
    """Attach the read-support evidence record to VCF SNVs only, leaving others as-is.

    Used by ``get_variants`` (which is readless, so BAMCP candidates must stay raw
    rather than be zero-scored from an empty read list) to give VCF-primary SNVs the
    same ``confidence``/``artifact_risk`` contract the viewer path produces. VCF
    indels/SVs and BAMCP candidates pass through unchanged.
    """
    out: list[dict] = []
    for v in variants:
        if v.get("source") == "vcf" and len(v["ref"]) == 1 == len(v["alt"]):
            enhanced, _ = _vcf_evidence(v, reference_sequence, region_start)
            out.append(enhanced)
        else:
            out.append(v)
    return out


def enhance_variants_with_evidence(
    variants: list[dict],
    reads: list[AlignedRead],
    reference_sequence: str | None,
    region_start: int,
) -> tuple[list[dict], dict]:
    """Attach BAMCP read-level evidence to BAMCP-candidate variants.

    For each BAMCP candidate, computes strand/quality/position evidence from the
    reads (via a mismatch index), an artifact-risk assessment, and a confidence
    level, returning the enhanced variants plus a keyed evidence map for the viewer.

    VCF SNVs are scored by :func:`_vcf_evidence` from the per-read metrics that
    :func:`bamcp.core.parsers.annotate_vcf_snv_support` pileups at each site (immune
    to the ``max_reads`` cap that bounds ``reads`` here); VCF indels/SVs, which BAMCP
    can't corroborate from single-base support, pass through unscored.
    """
    mismatch_index = build_mismatch_index(reads)
    indel_index: dict | None = None
    variant_evidence: dict = {}
    enhanced_variants: list[dict] = []

    for variant in variants:
        key = f"{variant['position']}:{variant['ref']}>{variant['alt']}"

        if variant.get("source") == "vcf":
            if len(variant["ref"]) == 1 == len(variant["alt"]):
                # Annotated SNV: build an evidence record from the count_coverage
                # support so the viewer has a confidence + evidence entry to filter on.
                enhanced, evidence = _vcf_evidence(variant, reference_sequence, region_start)
                variant_evidence[key] = evidence
                enhanced_variants.append(enhanced)
            else:
                # Indel/SV: BAMCP can't count single-base support, so pass the trusted
                # VCF record through unscored rather than zero-scoring it to low
                # confidence (which the viewer's high-confidence filter would hide).
                enhanced_variants.append(dict(variant))
            continue

        if is_indel(variant):
            if indel_index is None:
                indel_index = build_indel_index(reads)
            evidence = compute_indel_evidence(indel_index, variant)
        else:
            evidence = compute_variant_evidence(mismatch_index, variant)
        artifact_risk = compute_artifact_risk(variant, evidence, reference_sequence, region_start)
        evidence["artifact_risk"] = artifact_risk
        variant_evidence[key] = evidence

        enhanced = dict(variant)
        computed = evidence.get("read_evidence_computed", True)
        # Unknown read evidence is reported as null, never as a clean 0/0 / 0.0.
        enhanced["strand_forward"] = evidence["forward_count"] if computed else None
        enhanced["strand_reverse"] = evidence["reverse_count"] if computed else None
        enhanced["mean_quality"] = evidence["mean_quality"] if computed else None
        enhanced["read_evidence_computed"] = computed
        enhanced["artifact_risk"] = artifact_risk

        confidence = compute_confidence(variant, evidence, artifact_risk)
        enhanced["confidence"] = confidence
        # Keep is_low_confidence for backwards compatibility
        enhanced["is_low_confidence"] = confidence == "low"

        enhanced_variants.append(enhanced)

    return enhanced_variants, variant_evidence


# ---------------------------------------------------------------------------
# Plain-language interpretation
# ---------------------------------------------------------------------------

# Germline zygosity bands on VAF. Deliberately wide: allelic imbalance, mapping bias and
# sampling at modest depth all move a real het well off 0.5.
_HOM_ALT_MIN_VAF = 0.85
_HET_MIN_VAF = 0.30
_HET_MAX_VAF = 0.70


def zygosity_hint(vaf: float | None) -> str | None:
    """A germline zygosity hint from VAF — a hint, not a genotype call."""
    if not isinstance(vaf, (int, float)):
        return None
    if vaf >= _HOM_ALT_MIN_VAF:
        return "consistent with homozygous alt"
    if _HET_MIN_VAF <= vaf <= _HET_MAX_VAF:
        return "consistent with heterozygous"
    if vaf > _HET_MAX_VAF:
        return "between het and hom-alt (possible allelic imbalance)"
    return "below the germline het range (subclonal/mosaic, or an artifact)"


def variant_concerns(variant: dict) -> list[str]:
    """Artifact-risk types flagged for a variant, as a short list (empty when none)."""
    risk = variant.get("artifact_risk")
    if not isinstance(risk, dict):
        return []
    return [str(r.get("type")) for r in risk.get("risks") or [] if r.get("type")]


def interpret_variant(variant: dict, *, one_based: bool = True) -> str:
    """One plain-language sentence for a candidate/VCF variant.

    Says what was seen (alt reads / depth / VAF), a zygosity hint, the artifact risk with
    its concerns, and the confidence — so a newcomer can read a call without a glossary.
    ``position`` is printed as given when ``one_based`` (tool-boundary dicts), else +1.
    """
    pos = variant.get("position", 0) + (0 if one_based else 1)
    label = (
        f"{variant.get('contig', '?')}:{pos} {variant.get('ref', '?')}>{variant.get('alt', '?')}"
    )
    if variant.get("variant_kind") == "sv":
        return f"{label}: structural variant from the VCF; BAMCP does not assess SV evidence."
    vaf = variant.get("vaf")
    depth = variant.get("depth")
    alt = variant.get("alt_count")
    parts: list[str] = []
    if isinstance(vaf, (int, float)) and depth:
        seen = f"{alt}/{depth} reads" if alt is not None else f"depth {depth}"
        parts.append(f"{seen} (VAF {vaf:.0%})")
        hint = zygosity_hint(vaf)
        if hint:
            parts.append(hint)
    risk = variant.get("artifact_risk")
    if isinstance(risk, dict):
        likelihood = risk.get("artifact_likelihood", "unknown")
        concerns = [str(r.get("description") or r.get("type")) for r in risk.get("risks") or []]
        risk_text = f"artifact risk {likelihood}"
        if concerns:
            risk_text += " (" + "; ".join(concerns) + ")"
        parts.append(risk_text)
    else:
        parts.append("read-level artifact evidence not assessed")
    if variant.get("confidence"):
        parts.append(f"confidence {variant['confidence']}")
    source = "from your VCF" if variant.get("source") == "vcf" else "BAMCP candidate"
    return f"{label} ({source}): " + "; ".join(parts) + "."


def annotate_interpretation(variant: dict, *, one_based: bool = True) -> dict:
    """Return a copy of ``variant`` with ``concerns`` and ``interpretation`` attached."""
    out = dict(variant)
    out["concerns"] = variant_concerns(variant)
    out["interpretation"] = interpret_variant(variant, one_based=one_based)
    return out
