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
    # External sources disabled -> an explicit "disabled" block: neither fabricated
    # data nor None (which would read as a genuine not-found / rarity evidence).
    assert ev["clinvar"] == {"status": "disabled"}
    assert ev["population_frequency"] == {"status": "disabled"}


@pytest.mark.unit
async def test_assemble_evidence_absent_variant_not_detected(offline_cfg, comprehensive_bam_path):
    # A position with no planted variant should report not-detected, no crash. REF is the
    # real reference base there (C): a wrong REF is rejected, see the test below.
    ev = await assemble_evidence("chr1", 1500, "C", "G", comprehensive_bam_path, offline_cfg)
    assert ev["observation"]["detected"] is False


@pytest.mark.unit
async def test_assemble_evidence_rejects_ref_mismatch(offline_cfg, comprehensive_bam_path):
    # A REF that isn't the reference base would make the site probe count ordinary reads as
    # alt support (and the ClinVar/gnomAD lookups would query the wrong allele): caller error.
    from bamcp.core.errors import ToolInputError

    with pytest.raises(ToolInputError, match="does not match the reference"):
        await assemble_evidence("chr1", 1500, "A", "G", comprehensive_bam_path, offline_cfg)


@pytest.mark.unit
async def test_assemble_evidence_missing_contig_degrades(offline_cfg, comprehensive_bam_path):
    # chr17 is not in the fixture BAM -> observation is UNAVAILABLE (a failed read),
    # never "not detected" (which would read as evidence of absence).
    ev = await assemble_evidence(
        "chr17", 7674220, "C", "T", comprehensive_bam_path, offline_cfg, gene="TP53"
    )
    assert ev["observation"]["status"] == "unavailable"
    assert ev["observation"]["detected"] is None
    assert ev["observation"]["reason"]
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
    from bamcp.core.errors import ToolInputError

    # Invalid input is a tool error (isError), never a success payload carrying "error".
    with pytest.raises(ToolInputError, match="Position must be positive"):
        await handle_classify_variant(
            {
                "file_path": comprehensive_bam_path,
                "chrom": "chr1",
                "pos": -5,
                "ref": "A",
                "alt": "T",
            },
            offline_cfg,
        )


@pytest.mark.unit
async def test_clinical_context_disabled_returns_disabled_blocks(offline_cfg):
    ctx = await clinical_context(offline_cfg, "chr1", 100, "A", "T")
    assert ctx == {
        "clinvar": {"status": "disabled"},
        "population_frequency": {"status": "disabled"},
    }


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
    # External lookups disabled -> said so; never fabricated data, never "no record found".
    assert "ClinVar: disabled" in text
    assert "gnomAD: disabled" in text
    assert "no record found" not in text


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


# -- Review fixes: disabled sources, failed BAM reads, labels, PM2 QC -------------


@pytest.mark.unit
def test_disabled_sources_render_like_unavailable_not_absence():
    from bamcp.analysis.acmg import _clinical_criteria, _frequency_criteria

    freq = " ".join(c["evaluate_with"] for c in _frequency_criteria({"status": "disabled"}))
    assert "DISABLED" in freq
    assert "Do NOT apply PM2" in freq
    assert "Absence supports rarity" not in freq
    clin = _clinical_criteria({"status": "disabled"})[0]["evaluate_with"]
    assert "DISABLED" in clin
    assert "genuinely not" not in clin


@pytest.mark.unit
async def test_classify_text_with_disabled_sources(offline_cfg, comprehensive_bam_path):
    result = await handle_classify_variant(
        {"file_path": comprehensive_bam_path, "chrom": "chr1", "pos": 1051, "ref": "A", "alt": "T"},
        offline_cfg,
    )
    text = result["content"][0]["text"]
    assert "No ClinVar record found" not in text
    assert "No gnomAD record found" not in text
    assert text.count("DISABLED on this server") >= 2


@pytest.mark.unit
def test_clinvar_criteria_labelled_pp5_bp6_never_ps1():
    from bamcp.analysis.acmg import _clinical_criteria

    for block in (None, {"status": "unavailable"}, {"clinical_significance": "Pathogenic"}):
        crit = _clinical_criteria(block)
        assert all("PS1" not in c["code"] for c in crit)
        assert crit[0]["code"] == "PP5/BP6"
        assert "ClinGen SVI" in crit[0]["evaluate_with"]


@pytest.mark.unit
@pytest.mark.parametrize(
    "pop",
    [
        {"global_af": 0.00001, "max_pop_af": 0.00002, "filters": ["AC0"], "an": 150_000},
        {"global_af": 0.00001, "max_pop_af": 0.00002, "filters": [], "an": 812},
    ],
)
def test_pm2_refuses_filtered_or_low_an_sites(pop):
    from bamcp.analysis.acmg import _frequency_criteria

    pm2 = next(c for c in _frequency_criteria(pop) if c["code"] == "PM2_supporting")
    assert "NOT evidence of rarity" in pm2["evaluate_with"]
    assert "Do NOT" in pm2["evaluate_with"]


@pytest.mark.unit
def test_pm2_tolerates_missing_qc_fields():
    from bamcp.analysis.acmg import _frequency_criteria

    pm2 = next(
        c
        for c in _frequency_criteria({"global_af": 0.00001, "max_pop_af": 0.00002})
        if c["code"] == "PM2_supporting"
    )
    assert "NOT evidence of rarity" not in pm2["evaluate_with"]
    assert "Applies only when AF is below" in pm2["evaluate_with"]


@pytest.mark.unit
def test_population_block_surfaces_filters_and_an():
    from types import SimpleNamespace

    from bamcp.analysis.fusion import _population_block

    r = SimpleNamespace(
        global_af=0.001,
        populations=[SimpleNamespace(id="afr", af=0.002)],
        homozygote_count=0,
        ac=3,
        an=812,
        filters=["RF"],
        source="genome",
    )
    block = _population_block(r)
    assert block["an"] == 812 and block["filters"] == ["RF"]
    # max_pop_af comes only from the client's gated field, never a raw max over
    # `populations` (a lone allele in a tiny subgroup must not read as BA1).
    assert block["max_pop_af"] is None and block["max_pop"] is None
    gated = _population_block(
        SimpleNamespace(
            global_af=0.001,
            populations=[SimpleNamespace(id="adygei_XY", af=0.071)],
            max_pop_af=0.002,
            max_pop="afr",
            max_pop_an=41000,
            frequency_reliable=False,
            warnings=["low coverage"],
        )
    )
    assert gated["max_pop_af"] == 0.002 and gated["max_pop"] == "afr"
    assert gated["max_pop_an"] == 41000
    assert gated["frequency_reliable"] is False and gated["warnings"] == ["low coverage"]
    # A client surfacing fewer fields still yields a block.
    sparse = _population_block(SimpleNamespace(global_af=0.1))
    assert sparse["filters"] == [] and sparse["an"] is None


@pytest.mark.unit
async def test_assemble_evidence_failed_read_is_unavailable(offline_cfg, monkeypatch):
    async def boom(*_a, **_k):
        raise TimeoutError("Timed out after 30s reading chr1:1-100")

    monkeypatch.setattr("bamcp.core.tools._fetch_region_with_timeout", boom)
    ev = await assemble_evidence("chr1", 1051, "A", "T", "x.bam", offline_cfg)
    obs = ev["observation"]
    assert obs["status"] == "unavailable"
    assert obs["detected"] is None
    assert "timeout" in obs["reason"]

    scaffold = build_acmg_scaffold(ev)
    from bamcp.analysis.acmg import _format_classify_text

    text = _format_classify_text(ev, scaffold)
    assert "NOT observed" not in text
    assert "UNAVAILABLE" in text


@pytest.mark.unit
async def test_assemble_evidence_filtered_variant_echoes_filters(
    comprehensive_ref_fasta_path, comprehensive_bam_path
):
    """A real variant dropped by min_vaf must not read as a reference call: the requested
    site is evaluated anyway and reported as below_filters with the filter it fails."""
    cfg = BAMCPConfig(
        reference=comprehensive_ref_fasta_path,
        clinvar_enabled=False,
        gnomad_enabled=False,
        min_vaf=0.95,
    )
    ev = await assemble_evidence("chr1", 1051, "A", "T", comprehensive_bam_path, cfg)
    obs = ev["observation"]
    assert obs["status"] == "below_filters"
    assert obs["detected"] is False
    assert obs["passes_filters"] is False
    assert obs["failed_filters"] == ["min_vaf"]
    assert obs["applied_filters"]["min_vaf"] == 0.95
    assert obs["alt_reads_seen"] == 8
    assert obs["alt_count"] == 8 and obs["artifact_likelihood"] is not None

    from bamcp.analysis.acmg import _format_classify_text

    text = _format_classify_text(ev, build_acmg_scaffold(ev))
    assert "min_vaf=0.95" in text
    assert "NOT CALLED under the applied filters (fails min_vaf)" in text


@pytest.mark.unit
async def test_assemble_evidence_observed_variant_echoes_filters(
    offline_cfg, comprehensive_bam_path
):
    ev = await assemble_evidence("chr1", 1051, "A", "T", comprehensive_bam_path, offline_cfg)
    assert ev["observation"]["status"] == "observed"
    assert "min_vaf" in ev["observation"]["applied_filters"]


@pytest.mark.unit
async def test_assemble_evidence_lowercase_alleles(offline_cfg, comprehensive_bam_path):
    ev = await assemble_evidence("chr1", 1051, "a", "t", comprehensive_bam_path, offline_cfg)
    assert ev["observation"]["detected"] is True
    assert ev["variant"]["change"] == "A>T"


@pytest.mark.unit
async def test_assemble_evidence_overlaps_lookups_with_bam_read(
    comprehensive_ref_fasta_path, monkeypatch
):
    """The BAM read waits on the ClinVar lookup: only true concurrency finishes."""
    import asyncio

    started = asyncio.Event()

    class _SignalClient:
        async def lookup(self, *_a, **_k):
            started.set()
            return None

    async def fetch(*_a, **_k):
        await asyncio.wait_for(started.wait(), timeout=2)
        raise KeyError("chr1")

    monkeypatch.setattr("bamcp.core.tools.get_clinvar_client", lambda config: _SignalClient())
    monkeypatch.setattr("bamcp.core.tools._fetch_region_with_timeout", fetch)
    cfg = BAMCPConfig(
        reference=comprehensive_ref_fasta_path, clinvar_enabled=True, gnomad_enabled=False
    )
    ev = await assemble_evidence("chr1", 1051, "A", "T", "x.bam", cfg)
    assert ev["clinvar"] is None  # genuine not-found from the stub
    assert ev["observation"]["status"] == "unavailable"


@pytest.mark.unit
async def test_assemble_evidence_cancels_lookups_on_unexpected_error(
    comprehensive_ref_fasta_path, monkeypatch
):
    import asyncio

    cancelled = asyncio.Event()

    class _SlowClient:
        async def lookup(self, *_a, **_k):
            try:
                await asyncio.sleep(10)
            except asyncio.CancelledError:
                cancelled.set()
                raise

    async def fetch(*_a, **_k):
        await asyncio.sleep(0)
        raise RuntimeError("pysam exploded")

    monkeypatch.setattr("bamcp.core.tools.get_clinvar_client", lambda config: _SlowClient())
    monkeypatch.setattr("bamcp.core.tools._fetch_region_with_timeout", fetch)
    cfg = BAMCPConfig(
        reference=comprehensive_ref_fasta_path, clinvar_enabled=True, gnomad_enabled=False
    )
    with pytest.raises(RuntimeError):
        await assemble_evidence("chr1", 1051, "A", "T", "x.bam", cfg)
    await asyncio.wait_for(cancelled.wait(), timeout=1)


# -- Curation handler: error handling, filters, case ------------------------------


def _curation_args(path: str, **kw) -> dict:
    return {"file_path": path, "chrom": "chr1", "pos": 1051, "ref": "A", "alt": "T", **kw}


@pytest.mark.unit
@pytest.mark.parametrize(
    ("exc", "expected"),
    [
        (TimeoutError("slow"), "Timeout"),
        (KeyError("chrZ"), "contig not present"),
        (ValueError("bad region"), "contig not present"),
        (OSError("unreadable"), "Could not read the alignment file"),
    ],
)
async def test_curation_fetch_errors_become_error_results(offline_cfg, monkeypatch, exc, expected):
    from bamcp.analysis.curation import handle_get_variant_curation_summary

    async def boom(*_a, **_k):
        raise exc

    monkeypatch.setattr("bamcp.core.tools._fetch_region_with_timeout", boom)
    monkeypatch.setattr("bamcp.analysis.curation.validate_path", lambda *_a: None)
    result = await handle_get_variant_curation_summary(_curation_args("x.bam"), offline_cfg)
    assert expected in result["content"][0]["text"]


@pytest.mark.unit
async def test_curation_not_found_echoes_applied_filters(offline_cfg, comprehensive_bam_path):
    import json

    from bamcp.analysis.curation import handle_get_variant_curation_summary

    result = await handle_get_variant_curation_summary(
        _curation_args(comprehensive_bam_path, pos=1500, ref="A", alt="G", format="rubric"),
        offline_cfg,
    )
    payload = json.loads(result["content"][0]["text"])
    assert "error" in payload
    assert payload["applied_filters"]["min_vaf"] == offline_cfg.min_vaf

    text_result = await handle_get_variant_curation_summary(
        _curation_args(comprehensive_bam_path, pos=1500, ref="A", alt="G"), offline_cfg
    )
    assert "Applied filters:" in text_result["content"][0]["text"]


@pytest.mark.unit
async def test_curation_accepts_lowercase_alleles(offline_cfg, comprehensive_bam_path):
    from bamcp.analysis.curation import handle_get_variant_curation_summary

    result = await handle_get_variant_curation_summary(
        _curation_args(comprehensive_bam_path, ref="a", alt="t"), offline_cfg
    )
    assert result["content"][0]["text"].startswith("Variant Curation Summary: chr1:1051 A>T")


@pytest.mark.unit
def test_unreliable_frequency_blocks_pm2():
    from bamcp.analysis.acmg import _frequency_criteria

    pm2 = next(
        c
        for c in _frequency_criteria(
            {"global_af": 0.00001, "frequency_reliable": False, "warnings": ["AN below 2000"]}
        )
        if c["code"] == "PM2_supporting"
    )
    assert "NOT evidence of rarity" in pm2["evaluate_with"]
    assert "AN below 2000" in pm2["evaluate_with"]


@pytest.mark.unit
async def test_bam_build_routes_lookups_to_the_bams_assembly(monkeypatch):
    from bamcp.analysis import fusion
    from bamcp.config import BAMCPConfig
    from bamcp.core import tools

    async def grch37(file_path, config):
        return {"build": "GRCh37"}

    async def unknown(file_path, config):
        return {"build": "unknown"}

    async def broken(file_path, config):
        raise OSError("unreadable")

    cfg = BAMCPConfig()
    monkeypatch.setattr(tools, "_detect_bam_build", grch37)
    assert await fusion.bam_build("x.bam", cfg) == ("GRCh37", "detected-from-bam")
    monkeypatch.setattr(tools, "_detect_bam_build", unknown)
    build, source = await fusion.bam_build("x.bam", cfg)
    assert build is None and source.startswith("server-default")
    monkeypatch.setattr(tools, "_detect_bam_build", broken)
    build, source = await fusion.bam_build("x.bam", cfg)
    assert build is None and "unreadable" in source

    seen = {}

    class _Client:
        async def lookup(self, chrom, pos, ref, alt, build=None):
            seen["build"] = build
            return None

    monkeypatch.setattr(tools, "get_clinvar_client", lambda config: _Client())
    await fusion._lookup_clinvar_safe(cfg, "17", 7579472, "G", "C", "GRCh37")
    assert seen["build"] == "GRCh37"


@pytest.mark.unit
def test_gnomad_caution_is_direction_neutral_and_not_repeated():
    from bamcp.analysis.acmg import _benign_frequency_caveat, _format_classify_text

    pop = {
        "global_af": 0.207,
        "max_pop_af": 0.235,
        "filters": ["joint:discrepant_frequencies"],
        "frequency_reliable": False,
        "warnings": ["site fails gnomAD QC filters in joint (['joint:discrepant_frequencies'])"],
    }
    ev = {
        "variant": {"location": "chr17:43099914", "change": "G>A"},
        "observation": {"status": "observed", "detected": True},
        "clinvar": None,
        "population_frequency": pop,
    }
    text = _format_classify_text(ev, build_acmg_scaffold(ev))
    caution = next(line for line in text.splitlines() if "CAUTION" in line)
    # A common (AF 0.2) site must not be cautioned about "rarity".
    assert "NOT evidence of rarity" not in caution
    assert "unreliable in either direction" in caution
    # The failed-QC reason appears once, not restated by gnomAD's own warning.
    assert caution.count("discrepant_frequencies") == 1
    assert _benign_frequency_caveat(pop).count("discrepant_frequencies") == 1
