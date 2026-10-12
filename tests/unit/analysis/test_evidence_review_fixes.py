"""Regression tests for the read-level evidence review fixes.

Covers: read-position bias measured on aligned read ends (clean 2x150 data must not
flag), CIGAR-derived indel evidence (or an explicit not-computed record — never
zero-filled "clean" evidence), indel homopolymer context, and the confidence
thresholds living in constants.
"""

from __future__ import annotations

import json

import pytest

from bamcp.analysis.curation import (
    compute_rubric_scores,
    format_curation_summary,
    generate_curator_recommendations,
    handle_get_variant_curation_summary,
)
from bamcp.analysis.evidence import (
    _vcf_evidence,
    build_indel_index,
    build_mismatch_index,
    compute_artifact_risk,
    compute_confidence,
    compute_indel_evidence,
    compute_variant_evidence,
    enhance_variants_with_evidence,
    homopolymer_context_length,
    is_read_position_biased,
    summarize_read_positions,
)
from bamcp.config import BAMCPConfig
from bamcp.constants import (
    HIGH_CONFIDENCE_MAX_STRAND_BIAS,
    HIGH_CONFIDENCE_MIN_DEPTH,
    HIGH_CONFIDENCE_MIN_MEAN_QUALITY,
    HIGH_CONFIDENCE_MIN_VAF,
)
from bamcp.core.parsers import AlignedRead, RegionData

pytestmark = pytest.mark.unit

VARIANT_POS = 1_000  # 0-based reference position of the test SNV


def _snv_read(i: int, qpos: int, length: int, cigar: str | None = None) -> AlignedRead:
    """A read carrying an A>T mismatch at VARIANT_POS, placed so it sits at ``qpos``."""
    lead_clip = 0
    if cigar and "S" in cigar.split("M")[0]:
        lead_clip = int(cigar.split("S")[0])
    return AlignedRead(
        name=f"r{i}",
        sequence="A" * length,
        qualities=[35] * length,
        cigar=cigar or f"{length}M",
        position=VARIANT_POS - (qpos - lead_clip),
        end_position=VARIANT_POS - (qpos - lead_clip) + length,
        mapping_quality=60,
        is_reverse=bool(i % 2),
        mismatches=[{"pos": VARIANT_POS, "ref": "A", "alt": "T"}],
    )


def _snv_variant(depth: int = 60) -> dict:
    return {
        "position": VARIANT_POS,
        "ref": "A",
        "alt": "T",
        "vaf": 0.5,
        "depth": depth,
        "alt_count": depth // 2,
        "source": "bamcp",
    }


def _risk_types(reads: list[AlignedRead]) -> tuple[dict, list[str]]:
    variant = _snv_variant()
    evidence = compute_variant_evidence(build_mismatch_index(reads), variant)
    risk = compute_artifact_risk(variant, evidence, None, 0)
    return evidence, [r["type"] for r in risk["risks"]]


# -- 1. read-position bias -------------------------------------------------------


@pytest.mark.parametrize("length", [150, 100, 50])
def test_uniform_positions_do_not_flag_position_bias(length):
    """Clean data (support spread evenly along the read) must not read as an artifact."""
    reads = [_snv_read(i, qpos, length) for i, qpos in enumerate(range(length))]
    evidence, types = _risk_types(reads)
    assert "read_position_bias" not in types
    assert evidence["near_end_fraction"] <= 0.5


def test_clean_2x150_reaches_high_confidence():
    """The review's headline bug: mid-read 2x150 support was capped at medium."""
    reads = [_snv_read(i, 30 + i, 150) for i in range(30)]  # distances 30..59 from an end
    variant = _snv_variant()
    evidence = compute_variant_evidence(build_mismatch_index(reads), variant)
    risk = compute_artifact_risk(variant, evidence, None, 0)
    assert risk["artifact_likelihood"] == "low"
    assert compute_confidence(variant, evidence, risk) == "high"


def test_support_clustered_at_read_ends_flags():
    reads = [_snv_read(i, i % 6, 150) for i in range(20)]
    evidence, types = _risk_types(reads)
    assert "read_position_bias" in types
    assert evidence["near_end_fraction"] == 1.0


def test_distance_measured_from_aligned_end_not_soft_clip():
    """A base 3bp into the aligned block is near the end, even behind a 20bp soft clip."""
    reads = [_snv_read(i, 23, 150, cigar="20S130M") for i in range(10)]
    evidence, types = _risk_types(reads)
    assert evidence["near_end_count"] == 10
    assert "read_position_bias" in types


def test_short_reads_need_excess_over_uniform_expectation():
    """On 20bp reads every position is near an end; uniform support must not flag."""
    reads = [_snv_read(i, i % 20, 20) for i in range(40)]
    evidence, types = _risk_types(reads)
    assert evidence["near_end_expected_fraction"] == 1.0
    assert "read_position_bias" not in types


def test_is_read_position_biased_requires_positions():
    assert not is_read_position_biased({})
    assert not is_read_position_biased(summarize_read_positions([]))


def test_vcf_path_uses_absolute_cutoff():
    base = {
        "position": VARIANT_POS,
        "ref": "A",
        "alt": "T",
        "source": "vcf",
        "strand_forward": 10,
        "strand_reverse": 10,
        "read_depth": 40,
        "read_support_vaf": 0.5,
        "_alt_qualities": [35] * 20,
        "_alt_mapqs": [60] * 20,
    }
    uniform = {**base, "_alt_positions": list(range(0, 75, 4))[:20]}
    enhanced, _ = _vcf_evidence(uniform, None, 0)
    assert "read_position_bias" not in [r["type"] for r in enhanced["artifact_risk"]["risks"]]

    clustered = {**base, "_alt_positions": [i % 4 for i in range(20)]}
    enhanced, _ = _vcf_evidence(clustered, None, 0)
    assert "read_position_bias" in [r["type"] for r in enhanced["artifact_risk"]["risks"]]


# -- 4. indel read evidence ------------------------------------------------------


def _indel_read(i: int, cigar: str, start: int = 100, length: int = 50) -> AlignedRead:
    return AlignedRead(
        name=f"ind{i}",
        sequence="ACGTACGTAC" * (length // 10),
        qualities=[30] * length,
        cigar=cigar,
        position=start,
        end_position=start + length,
        mapping_quality=60,
        is_reverse=bool(i % 2),
    )


def test_build_indel_index_keys_match_vcf_anchor():
    # 20M3D30M at 100: anchor = last aligned ref base = 119; 20M2I28M: anchor 119, ins
    # seq = query[20:22].
    reads = [_indel_read(0, "20M3D30M"), _indel_read(1, "20M2I28M")]
    index = build_indel_index(reads)
    assert (119, "del", 3) in index
    ins_seq = reads[1].sequence[20:22]
    assert (119, "ins", ins_seq) in index
    # Insertion query position = first inserted base; deletion = anchor base.
    assert index[(119, "ins", ins_seq)][0][1] == 20
    assert index[(119, "del", 3)][0][1] == 19


def test_indel_evidence_computed_from_cigar():
    reads = [_indel_read(i, "20M3D30M") for i in range(10)]
    variant = {
        "position": 119,
        "ref": "ACGT",
        "alt": "A",
        "vaf": 0.5,
        "depth": 20,
        "alt_count": 10,
        "variant_kind": "indel",
        "source": "bamcp",
    }
    enhanced, evmap = enhance_variants_with_evidence([variant], reads, None, 0)
    ev = evmap["119:ACGT>A"]
    assert ev["read_evidence_computed"] is True
    assert ev["forward_count"] == 5 and ev["reverse_count"] == 5
    assert ev["mean_quality"] == 30
    assert enhanced[0]["strand_forward"] == 5


def _orphan_indel() -> dict:
    return {
        "position": 119,
        "ref": "A",
        "alt": "ATT",
        "vaf": 0.5,
        "depth": 30,
        "alt_count": 15,
        "variant_kind": "indel",
        "source": "bamcp",
    }


def test_indel_without_supporting_reads_is_not_computed_not_clean():
    """No supporting reads located → null scores + explicit concern, never 'balanced'."""
    enhanced, evmap = enhance_variants_with_evidence([_orphan_indel()], [], None, 0)
    v = enhanced[0]
    assert v["read_evidence_computed"] is False
    assert v["strand_forward"] is None and v["strand_reverse"] is None
    assert v["mean_quality"] is None
    assert v["confidence"] == "low"
    types = [r["type"] for r in v["artifact_risk"]["risks"]]
    assert "read_evidence_unavailable" in types
    assert "strand_bias" not in types

    evidence = evmap["119:A>ATT"]
    scores = compute_rubric_scores(v, evidence, v["artifact_risk"])
    assert scores["strand_balance"] is None
    assert scores["position_quality"] is None
    assert scores["mapq_quality"] is None

    recs = generate_curator_recommendations(v, evidence, v["artifact_risk"])
    joined = " ".join(recs)
    assert "NOT COMPUTED" in joined
    assert "No automated quality flags raised" not in joined


def test_complex_indel_allele_not_computed():
    ev = compute_indel_evidence({}, {"position": 5, "ref": "AC", "alt": "GTT"})
    assert ev["read_evidence_computed"] is False


def test_snv_missing_from_capped_read_set_not_computed():
    variant = _snv_variant()
    ev = compute_variant_evidence(build_mismatch_index([]), variant)
    assert ev["read_evidence_computed"] is False
    risk = compute_artifact_risk(variant, ev, None, 0)
    assert compute_confidence(variant, ev, risk) == "low"


async def test_curation_handler_indel_not_computed_text(monkeypatch):
    data = RegionData(
        contig="chr1",
        start=69,
        end=170,
        reads=[],
        coverage=[30] * 101,
        variants=[{**_orphan_indel(), "contig": "chr1"}],
    )

    async def fake_fetch(*_a, **_k):
        return data

    monkeypatch.setattr("bamcp.core.tools._fetch_region_with_timeout", fake_fetch)
    monkeypatch.setattr("bamcp.analysis.curation.validate_path", lambda *_a: None)
    args = {"file_path": "x.bam", "chrom": "chr1", "pos": 120, "ref": "A", "alt": "ATT"}
    result = await handle_get_variant_curation_summary(args, BAMCPConfig())
    text = result["content"][0]["text"]
    assert "NOT COMPUTED" in text
    assert "balanced" not in text.lower().replace("not evidence of balance", "")

    rubric = await handle_get_variant_curation_summary({**args, "format": "rubric"}, BAMCPConfig())
    payload = json.loads(rubric["content"][0]["text"])
    assert payload["read_evidence_computed"] is False
    assert payload["strand_analysis"]["strand_bias"] is None
    assert payload["scores"]["strand_balance"] is None


def test_format_summary_renders_computed_evidence():
    summary = {
        "variant": {"location": "chr1:1", "change": "A>T", "vaf": 0.5, "depth": 30},
        "quality_metrics": {"mean_quality": 30.0, "quality_distribution": [0, 0, 0, 5, 5]},
        "strand_analysis": {"forward": 5, "reverse": 5, "strand_bias": 0.0, "is_balanced": True},
        "position_in_read": {"near_end_fraction": 0.1},
        "artifact_assessment": {"artifact_likelihood": "low", "risk_score": 0.0, "risks": []},
        "confidence": "high",
        "recommendations": ["x"],
    }
    text = format_curation_summary(summary)
    assert "Strand Bias: 0.00 (balanced)" in text
    assert "Near-end fraction: 10.0%" in text


# -- 5. indel homopolymer context ------------------------------------------------


def test_deletion_homopolymer_measured_past_anchor():
    seq = "GTCAAAAAAGTC"  # anchor C at index 2, then a 6bp A-run
    deletion = {"position": 2, "ref": "CA", "alt": "C"}
    assert homopolymer_context_length(deletion, seq, 2) == 6
    risk = compute_artifact_risk(
        {**deletion, "depth": 50}, {"strand_bias": 0.0}, reference_sequence=seq, region_start=0
    )
    hp = [r for r in risk["risks"] if r["type"] == "homopolymer"]
    assert hp and hp[0]["value"] == 6 and hp[0]["severity"] == "high"


def test_insertion_homopolymer_measured_past_anchor():
    seq = "GTCAAAAAGTC"
    insertion = {"position": 2, "ref": "C", "alt": "CA"}
    assert homopolymer_context_length(insertion, seq, 2) == 5


def test_snv_homopolymer_unchanged():
    assert homopolymer_context_length({"ref": "C", "alt": "T"}, "GTCAAAAAAGTC", 2) == 1


# -- 10. confidence thresholds come from constants --------------------------------


def test_compute_confidence_uses_high_confidence_constants():
    evidence = {
        "mean_quality": HIGH_CONFIDENCE_MIN_MEAN_QUALITY,
        "strand_bias": HIGH_CONFIDENCE_MAX_STRAND_BIAS,
    }
    clean = {"artifact_likelihood": "low"}
    at_bar = {"vaf": HIGH_CONFIDENCE_MIN_VAF, "depth": HIGH_CONFIDENCE_MIN_DEPTH}
    assert compute_confidence(at_bar, evidence, clean) == "high"
    below = {"vaf": HIGH_CONFIDENCE_MIN_VAF, "depth": HIGH_CONFIDENCE_MIN_DEPTH - 1}
    assert compute_confidence(below, evidence, clean) != "high"
