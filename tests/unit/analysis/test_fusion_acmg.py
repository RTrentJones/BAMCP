"""Unit tests for evidence fusion + the classify_variant / ACMG scaffold."""

from __future__ import annotations

import pytest

from bamcp.analysis.acmg import build_acmg_scaffold, handle_classify_variant
from bamcp.analysis.fusion import assemble_evidence, clinical_context
from bamcp.config import BAMCPConfig


@pytest.fixture
def offline_cfg(comprehensive_ref_fasta_path):
    """Config with external lookups disabled so fusion stays offline/deterministic."""
    return BAMCPConfig(
        reference=comprehensive_ref_fasta_path,
        clinvar_enabled=False,
        gnomad_enabled=False,
    )


@pytest.mark.unit
async def test_assemble_evidence_detects_fixture_variant(offline_cfg, comprehensive_bam_path):
    # chr1:1051 A>T is the planted strand-bias site (0-based 1050).
    ev = await assemble_evidence("chr1", 1051, "A", "T", comprehensive_bam_path, offline_cfg)
    assert ev["observation"]["detected"] is True
    assert ev["observation"]["confidence"] in ("low", "medium", "high")
    # External sources disabled -> blocks are None, not fabricated.
    assert ev["clinvar"] is None
    assert ev["population_frequency"] is None


@pytest.mark.unit
async def test_assemble_evidence_absent_variant_not_detected(offline_cfg, comprehensive_bam_path):
    # A position with no planted variant should report not-detected, no crash.
    ev = await assemble_evidence("chr1", 1500, "A", "G", comprehensive_bam_path, offline_cfg)
    assert ev["observation"]["detected"] is False


@pytest.mark.unit
async def test_assemble_evidence_missing_contig_degrades(offline_cfg, comprehensive_bam_path):
    # chr17 is not in the fixture BAM -> observation degrades to not-detected.
    ev = await assemble_evidence(
        "chr17", 7674220, "C", "T", comprehensive_bam_path, offline_cfg, gene="TP53"
    )
    assert ev["observation"]["detected"] is False
    assert ev["gene_context"]["gene"] == "TP53"


@pytest.mark.unit
def test_build_scaffold_does_not_state_a_classification():
    evidence = {
        "variant": {"location": "chr1:100", "change": "A>T"},
        "observation": {"detected": True, "depth": 30, "confidence": "high"},
        "clinvar": None,
        "population_frequency": {"global_af": 0.2, "max_pop_af": 0.25, "max_pop": "afr"},
    }
    scaffold = build_acmg_scaffold(evidence)
    assert "criteria_to_evaluate" in scaffold
    assert scaffold["classification_options"][0] == "Pathogenic"
    # High-frequency variant: BA1 must be offered as a criterion to evaluate,
    # but the scaffold must NOT pre-decide the classification.
    codes = {c["code"] for c in scaffold["criteria_to_evaluate"]}
    assert any("BA1" in c for c in codes)
    assert "final_classification" not in scaffold  # only in response_format template
    assert scaffold["response_format"]["final_classification"]


@pytest.mark.unit
async def test_handle_classify_variant_offline(offline_cfg, comprehensive_bam_path):
    result = await handle_classify_variant(
        {
            "file_path": comprehensive_bam_path,
            "chrom": "chr1",
            "pos": 1051,
            "ref": "A",
            "alt": "T",
        },
        offline_cfg,
    )
    text = result["content"][0]["text"]
    assert "ACMG classification scaffold" in text
    assert "INTENDED USE" in text
    sc = result["structuredContent"]
    assert "evidence" in sc and "scaffold" in sc


@pytest.mark.unit
async def test_handle_classify_variant_rejects_bad_input(offline_cfg, comprehensive_bam_path):
    import json

    result = await handle_classify_variant(
        {
            "file_path": comprehensive_bam_path,
            "chrom": "chr1",
            "pos": -5,
            "ref": "A",
            "alt": "T",
        },
        offline_cfg,
    )
    payload = json.loads(result["content"][0]["text"])
    assert "error" in payload


@pytest.mark.unit
async def test_clinical_context_disabled_returns_empty_blocks(offline_cfg):
    ctx = await clinical_context(offline_cfg, "chr1", 100, "A", "T")
    assert ctx == {"clinvar": None, "population_frequency": None}


@pytest.mark.unit
async def test_curation_include_clinical_bridge(offline_cfg, comprehensive_bam_path):
    """get_variant_curation_summary(include_clinical=True) attaches a context block."""
    from bamcp.analysis.curation import handle_get_variant_curation_summary

    result = await handle_get_variant_curation_summary(
        {
            "file_path": comprehensive_bam_path,
            "chrom": "chr1",
            "pos": 1051,
            "ref": "A",
            "alt": "T",
            "include_clinical": True,
        },
        offline_cfg,
    )
    text = result["content"][0]["text"]
    assert "CLINICAL CONTEXT" in text
    # External lookups disabled -> "no record found", never fabricated data.
    assert "no record found" in text


class _RaisingClient:
    """A lookup client whose every call errors — stands in for a network failure."""

    async def lookup(self, *_a, **_k):
        raise RuntimeError("network down")


# -- Phase 2: error-vs-absence -------------------------------------------------


@pytest.mark.unit
def test_frequency_criteria_unavailable_does_not_apply_pm2_on_absence():
    """An ERRORED gnomAD lookup must not become PM2 'absence supports rarity' evidence."""
    from bamcp.analysis.acmg import _frequency_criteria

    crit = _frequency_criteria({"status": "unavailable"})
    text = " ".join(c["evaluate_with"] for c in crit)
    assert "UNAVAILABLE" in text
    assert "Do NOT apply PM2" in text
    # Must NOT emit the genuine-absence phrasing that would license rarity.
    assert "Absence supports rarity" not in text


@pytest.mark.unit
def test_frequency_criteria_none_is_genuine_absence():
    """A genuine not-found (None) still yields the normal PM2-on-absence line."""
    from bamcp.analysis.acmg import _frequency_criteria

    crit = _frequency_criteria(None)
    assert "No gnomAD entry was found" in crit[0]["evaluate_with"]


@pytest.mark.unit
def test_clinical_criteria_unavailable_is_distinct_from_not_found():
    from bamcp.analysis.acmg import _clinical_criteria

    unavail = _clinical_criteria({"status": "unavailable"})[0]["evaluate_with"]
    assert "UNAVAILABLE" in unavail
    not_found = _clinical_criteria(None)[0]["evaluate_with"]
    assert "genuinely not" in not_found
    assert "UNAVAILABLE" not in not_found


@pytest.mark.unit
async def test_assemble_evidence_marks_failed_lookups_unavailable(
    comprehensive_ref_fasta_path, comprehensive_bam_path, monkeypatch
):
    """A failed ClinVar/gnomAD fetch → {'status': 'unavailable'}, never None (absence)."""
    monkeypatch.setattr("bamcp.core.tools.get_clinvar_client", lambda config: _RaisingClient())
    monkeypatch.setattr("bamcp.core.tools.get_gnomad_client", lambda config: _RaisingClient())
    cfg = BAMCPConfig(
        reference=comprehensive_ref_fasta_path, clinvar_enabled=True, gnomad_enabled=True
    )

    ev = await assemble_evidence("chr1", 1051, "A", "T", comprehensive_bam_path, cfg)
    assert ev["clinvar"] == {"status": "unavailable"}
    assert ev["population_frequency"] == {"status": "unavailable"}

    # And the scaffold built over it refuses to treat the failure as absence.
    scaffold = build_acmg_scaffold(ev)
    blob = str(scaffold)
    assert "UNAVAILABLE" in blob


@pytest.mark.unit
async def test_clinical_context_marks_failed_lookups_unavailable(
    comprehensive_ref_fasta_path, monkeypatch
):
    monkeypatch.setattr("bamcp.core.tools.get_clinvar_client", lambda config: _RaisingClient())
    monkeypatch.setattr("bamcp.core.tools.get_gnomad_client", lambda config: _RaisingClient())
    cfg = BAMCPConfig(
        reference=comprehensive_ref_fasta_path, clinvar_enabled=True, gnomad_enabled=True
    )

    ctx = await clinical_context(cfg, "chr1", 100, "A", "T")
    assert ctx["clinvar"] == {"status": "unavailable"}
    assert ctx["population_frequency"] == {"status": "unavailable"}
