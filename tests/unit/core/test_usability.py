"""Newcomer-usability contracts: auto-reference, contig-style mapping, compact list_contigs,
noise filters + inline per-variant evidence, and a compact viewer text summary."""

from __future__ import annotations

import copy
import json

import pysam
import pytest

from bamcp.analysis.acmg import _frequency_criteria
from bamcp.analysis.evidence import (
    annotate_interpretation,
    compute_artifact_risk,
    interpret_variant,
    zygosity_hint,
)
from bamcp.config import BAMCPConfig
from bamcp.core import reference as reference_mod
from bamcp.core import tools as tools_mod
from bamcp.core.parsers import RegionData, contig_not_found_message, scan_contig
from bamcp.core.reference import alternate_contig_names, map_contig_to_style
from bamcp.core.tools import (
    VIEWER_SUMMARY_MAX_CHARS,
    _viewer_summary,
    handle_get_region_summary,
    handle_get_variants,
    handle_jump_to,
    handle_list_contigs,
    handle_search_gene,
    handle_visualize_region,
)
from bamcp.core.validation import validate_remote_url

REGION = "chr1:1000-2900"


def _payload(result: dict) -> dict:
    return json.loads(result["content"][0]["text"])


def _write_header_bam(path, contigs: list[tuple[str, int]]) -> str:
    header = {"HD": {"VN": "1.6"}, "SQ": [{"SN": n, "LN": ln} for n, ln in contigs]}
    with pysam.AlignmentFile(str(path), "wb", header=header):
        pass
    return str(path)


@pytest.fixture
def auto_registry(monkeypatch, comprehensive_ref_fasta_path):
    """Point the GRCh38 chr-style registry entry at the local comprehensive FASTA, and make
    the comprehensive BAM's chr1 (10,000bp) read as GRCh38."""
    builds = copy.deepcopy(reference_mod.GENOME_BUILDS)
    builds["GRCh38"]["chr1_length"] = 10_000
    builds["GRCh38"]["fasta_urls"][reference_mod.CONTIG_STYLE_CHR] = comprehensive_ref_fasta_path
    monkeypatch.setattr(reference_mod, "GENOME_BUILDS", builds)
    return comprehensive_ref_fasta_path


# -- Auto-reference -----------------------------------------------------------


class TestAutoReference:
    @pytest.mark.unit
    @pytest.mark.asyncio
    async def test_get_variants_uses_auto_reference(self, comprehensive_bam_path, auto_registry):
        cfg = BAMCPConfig(allow_remote_files=True)
        p = _payload(
            await handle_get_variants({"file_path": comprehensive_bam_path, "region": REGION}, cfg)
        )
        assert p["applied_filters"]["reference_source"] == "auto"
        assert p["applied_filters"]["snv_detection"] == "enabled"
        assert p["count"] > 0

    @pytest.mark.unit
    @pytest.mark.asyncio
    async def test_viewer_echoes_auto_reference_and_refetch_validates(
        self, comprehensive_bam_path, auto_registry
    ):
        cfg = BAMCPConfig(allow_remote_files=True)
        r = await handle_visualize_region(
            {"file_path": comprehensive_bam_path, "region": "chr1:1040-1060"}, cfg
        )
        payload = r["_meta"]["ui/init"]
        assert payload["applied_filters"]["reference_source"] == "auto"
        # Echoed so the viewer's pan/zoom refetches keep it ...
        assert payload["reference"] == auto_registry
        # ... and a refetch passing it back as caller input is accepted.
        r2 = await handle_visualize_region(
            {
                "file_path": comprehensive_bam_path,
                "region": "chr1:1060-1080",
                "reference": payload["reference"],
            },
            cfg,
        )
        assert r2["_meta"]["ui/init"]["applied_filters"]["reference_source"] == "explicit"

    @pytest.mark.unit
    @pytest.mark.asyncio
    async def test_no_auto_reference_without_remote_files(
        self, comprehensive_bam_path, auto_registry
    ):
        p = _payload(
            await handle_get_variants(
                {"file_path": comprehensive_bam_path, "region": REGION}, BAMCPConfig()
            )
        )
        assert p["applied_filters"]["reference_source"] == "none"
        assert p["applied_filters"]["snv_detection"].startswith("disabled")

    @pytest.mark.unit
    @pytest.mark.asyncio
    async def test_auto_reference_respects_host_allow_list(
        self, comprehensive_bam_path, monkeypatch
    ):
        builds = copy.deepcopy(reference_mod.GENOME_BUILDS)
        builds["GRCh38"]["chr1_length"] = 10_000
        builds["GRCh38"]["fasta_urls"]["chr"] = "https://not-allowed.example/ref.fa"
        monkeypatch.setattr(reference_mod, "GENOME_BUILDS", builds)
        cfg = BAMCPConfig(allow_remote_files=True, allowed_remote_hosts=["ftp.ensembl.org"])
        p = _payload(
            await handle_get_variants({"file_path": comprehensive_bam_path, "region": REGION}, cfg)
        )
        assert p["applied_filters"]["reference_source"] == "none"
        assert "not permitted" in p["applied_filters"]["reference_note"]

    @pytest.mark.unit
    @pytest.mark.asyncio
    async def test_config_and_explicit_sources(
        self, comprehensive_bam_path, comprehensive_ref_fasta_path
    ):
        cfg = BAMCPConfig(reference=comprehensive_ref_fasta_path)
        p = _payload(
            await handle_get_variants({"file_path": comprehensive_bam_path, "region": REGION}, cfg)
        )
        assert p["applied_filters"]["reference_source"] == "config"
        p = _payload(
            await handle_get_variants(
                {
                    "file_path": comprehensive_bam_path,
                    "region": REGION,
                    "reference": comprehensive_ref_fasta_path,
                },
                BAMCPConfig(),
            )
        )
        assert p["applied_filters"]["reference_source"] == "explicit"

    @pytest.mark.unit
    @pytest.mark.asyncio
    async def test_header_cache_is_bounded(self, tmp_path, monkeypatch):
        monkeypatch.setattr(tools_mod, "_HEADER_CACHE_MAX", 2)
        cfg = BAMCPConfig()
        for i in range(4):
            bam = _write_header_bam(tmp_path / f"h{i}.bam", [("chr1", 100)])
            await tools_mod._bam_contigs(bam, cfg)
        assert len(tools_mod.get_services(cfg).header_cache) == 2


# -- reference_source echo on every variant + viewer response -----------------


@pytest.mark.unit
@pytest.mark.asyncio
async def test_reference_source_echoed_everywhere(
    comprehensive_bam_path, comprehensive_ref_fasta_path
):
    cfg = BAMCPConfig(reference=comprehensive_ref_fasta_path)
    args = {"file_path": comprehensive_bam_path, "region": "chr1:1040-1060"}
    viz = await handle_visualize_region(args, cfg)
    assert viz["_meta"]["ui/init"]["applied_filters"]["reference_source"] == "config"
    # A config reference is never echoed as caller input.
    assert "reference" not in viz["_meta"]["ui/init"]
    jump = await handle_jump_to(
        {"file_path": comprehensive_bam_path, "contig": "chr1", "position": 1051, "window": 40},
        cfg,
    )
    assert jump["_meta"]["ui/init"]["applied_filters"]["reference_source"] == "config"
    summary = await handle_get_region_summary(args, cfg)
    assert "reference: config" in summary["content"][0]["text"]


# -- Contig naming style (P0-2) -----------------------------------------------


class TestContigStyle:
    @pytest.mark.unit
    def test_alternate_names(self):
        assert alternate_contig_names("chr13") == ["13"]
        assert alternate_contig_names("13") == ["chr13"]
        assert "MT" in alternate_contig_names("chrM")
        assert map_contig_to_style("chr13", {"13", "14"}) == "13"
        assert map_contig_to_style("13", {"chr13"}) == "chr13"
        assert map_contig_to_style("chrM", {"MT"}) == "MT"
        assert map_contig_to_style("chr13", {"chr14"}) is None

    @pytest.mark.unit
    def test_not_found_message_suggests_alternate(self):
        msg = contig_not_found_message("chr13", ["13", "14"])
        assert "'13'" in msg and "not found" in msg
        assert "list_contigs" in contig_not_found_message("chrZ", ["1"])

    @pytest.mark.unit
    @pytest.mark.asyncio
    async def test_wrong_style_region_error_names_alternate(self, tmp_path):
        bam = _write_header_bam(tmp_path / "nochr.bam", [("1", 5000)])
        with pytest.raises(ValueError, match="retry with '1'"):
            await handle_get_variants({"file_path": bam, "region": "chr1:100-200"}, BAMCPConfig())

    @pytest.mark.unit
    @pytest.mark.asyncio
    async def test_search_gene_maps_to_nochr_bam(self, tmp_path, monkeypatch):
        from bamcp.clients.genes import GeneInfo

        bam = _write_header_bam(tmp_path / "g37.bam", [("1", 249250621), ("13", 115169878)])

        class _Client:
            async def search(self, symbol, build=None):
                return GeneInfo("BRCA2", "BRCA2", "chr13", 32889611, 32973805, "+", build=build)

        monkeypatch.setattr(tools_mod, "get_gene_client", lambda config: _Client())
        p = _payload(await handle_search_gene({"symbol": "BRCA2", "file_path": bam}, BAMCPConfig()))
        assert p["region"] == "13:32889611-32973805"
        assert p["contig_mapped_from"] == "chr13"
        assert p["bam_build"] == "GRCh37"

    @pytest.mark.unit
    @pytest.mark.asyncio
    async def test_search_gene_warns_when_contig_absent(self, tmp_path, monkeypatch):
        from bamcp.clients.genes import GeneInfo

        bam = _write_header_bam(tmp_path / "g.bam", [("1", 249250621)])

        class _Client:
            async def search(self, symbol, build=None):
                return GeneInfo("BRCA2", "BRCA2", "chr13", 1, 2, "+", build=build)

        monkeypatch.setattr(tools_mod, "get_gene_client", lambda config: _Client())
        p = _payload(await handle_search_gene({"symbol": "BRCA2", "file_path": bam}, BAMCPConfig()))
        assert "not in the BAM header" in p["warning"]


# -- list_contigs (P0-4) -------------------------------------------------------


def _hg002_like_contigs() -> list[tuple[str, int]]:
    primary = [(f"chr{i}", 248956422 if i == 1 else 1_000_000) for i in range(1, 23)]
    primary += [("chrX", 156040895), ("chrY", 57227415), ("chrM", 16569)]
    other: list[tuple[str, int]] = []
    for i in range(2580 - len(primary) - 1):
        kind = i % 4
        if kind == 0:
            other.append((f"chr1_KI27076{i}v1_alt", 1000 + i))
        elif kind == 1:
            other.append((f"chrUn_JTFH0100{i}v1_decoy", 1000 + i))
        elif kind == 2:
            other.append((f"HLA-A*01:01:01:{i:02d}", 3000 + i))
        else:
            other.append((f"chrUn_KN70{i}v1", 2000 + i))
    return primary + other + [("chrEBV", 171823)]


class TestListContigs:
    @pytest.mark.unit
    @pytest.mark.asyncio
    async def test_summary_first_primary_only_and_small(self, tmp_path):
        bam = _write_header_bam(tmp_path / "hg002.bam", _hg002_like_contigs())
        result = await handle_list_contigs({"file_path": bam}, BAMCPConfig())
        text = result["content"][0]["text"]
        assert len(text.encode()) < 5000
        p = json.loads(text)
        keys = list(p)
        assert keys[0] == "genome_build" and keys.index("contigs") > keys.index("contig_counts")
        assert p["genome_build"]["build"] == "GRCh38"
        assert p["contig_style"] == "chr"
        assert p["suggested_reference_url"]
        assert p["contigs_shown"] == "primary"
        assert len(p["contigs"]) == 25
        assert p["contig_counts"]["total"] == 2580
        assert p["contig_counts"]["primary"] == 25
        assert p["contig_counts"]["alt"] > 0 and p["contig_counts"]["decoy"] > 0
        assert p["contigs_omitted"] == 2555
        assert "include_all" in p["contigs_note"]
        assert p["coordinate_system"] == "1-based-inclusive"

    @pytest.mark.unit
    @pytest.mark.asyncio
    async def test_include_all(self, tmp_path):
        bam = _write_header_bam(tmp_path / "hg002.bam", _hg002_like_contigs())
        p = _payload(
            await handle_list_contigs({"file_path": bam, "include_all": True}, BAMCPConfig())
        )
        assert p["contigs_shown"] == "all"
        assert len(p["contigs"]) == 2580
        assert p["contigs_omitted"] == 0

    @pytest.mark.unit
    @pytest.mark.asyncio
    async def test_non_human_falls_back_to_listing(self, tmp_path):
        bam = _write_header_bam(tmp_path / "x.bam", [("contigA", 10), ("contigB", 20)])
        p = _payload(await handle_list_contigs({"file_path": bam}, BAMCPConfig()))
        assert p["contigs_shown"] == "all"
        assert [c["name"] for c in p["contigs"]] == ["contigA", "contigB"]


# -- Noise filters + inline evidence (P1-5) -------------------------------------


class TestMinAltCount:
    @pytest.mark.unit
    @pytest.mark.asyncio
    async def test_filters_and_echoes(self, small_bam_path, ref_fasta_path):
        cfg = BAMCPConfig(reference=ref_fasta_path, min_depth=1, min_vaf=0.05)
        args = {"file_path": small_bam_path, "region": "chr1:490-600"}
        loose = _payload(await handle_get_variants({**args, "min_alt_count": 1}, cfg))
        strict = _payload(await handle_get_variants(args, cfg))
        assert loose["applied_filters"]["min_alt_count"] == 1
        assert strict["applied_filters"]["min_alt_count"] == 3
        assert loose["count"] > strict["count"]
        assert all(v["alt_count"] >= 3 for v in strict["variants"])

    @pytest.mark.unit
    def test_defaults(self):
        cfg = BAMCPConfig()
        assert (cfg.min_vaf, cfg.min_mapq, cfg.min_baseq, cfg.min_alt_count) == (0.10, 20, 13, 3)
        with pytest.raises(ValueError, match="min_alt_count"):
            BAMCPConfig(min_alt_count=0)

    @pytest.mark.unit
    def test_env(self, monkeypatch):
        monkeypatch.setenv("BAMCP_MIN_ALT_COUNT", "5")
        assert BAMCPConfig.from_env().min_alt_count == 5

    @pytest.mark.unit
    def test_scan_contig_applies_alt_count_before_cap(
        self, comprehensive_bam_path, comprehensive_ref_fasta_path
    ):
        loose = scan_contig(
            comprehensive_bam_path, "chr1", comprehensive_ref_fasta_path, min_vaf=0.01, min_depth=1
        )
        strict = scan_contig(
            comprehensive_bam_path,
            "chr1",
            comprehensive_ref_fasta_path,
            min_vaf=0.01,
            min_depth=1,
            min_alt_count=5,
        )
        assert strict.total_candidates < loose.total_candidates
        assert all(v["alt_count"] >= 5 for v in strict.variants)

    @pytest.mark.unit
    @pytest.mark.asyncio
    async def test_viewer_applies_config_alt_count(self, small_bam_path, ref_fasta_path):
        cfg = BAMCPConfig(reference=ref_fasta_path, min_depth=1, min_vaf=0.05)
        r = await handle_visualize_region(
            {"file_path": small_bam_path, "region": "chr1:490-600"}, cfg
        )
        payload = r["_meta"]["ui/init"]
        assert payload["applied_filters"]["min_alt_count"] == 3
        assert all(v["alt_count"] >= 3 for v in payload["variants"])


class TestInlineEvidence:
    @pytest.mark.unit
    @pytest.mark.asyncio
    async def test_every_variant_has_evidence_and_interpretation(
        self, comprehensive_bam_path, comprehensive_ref_fasta_path
    ):
        cfg = BAMCPConfig(reference=comprehensive_ref_fasta_path)
        p = _payload(
            await handle_get_variants({"file_path": comprehensive_bam_path, "region": REGION}, cfg)
        )
        assert p["count"] > 0
        for v in p["variants"]:
            assert v["artifact_risk"]["artifact_likelihood"] in ("low", "medium", "high")
            assert v["confidence"] in ("low", "medium", "high")
            assert isinstance(v["concerns"], list)
            assert f":{v['position']} " in v["interpretation"]
            assert "_alt_qualities" not in v
        pse = p["applied_filters"]["per_site_evidence"]
        assert pse["capped"] is False and pse["scored_sites"] > 0
        # The planted strand-bias site is flagged inline.
        sb = next(v for v in p["variants"] if v["position"] == 1051)
        assert "strand_bias" in sb["concerns"]

    @pytest.mark.unit
    @pytest.mark.asyncio
    async def test_evidence_cap_is_flagged(
        self, comprehensive_bam_path, comprehensive_ref_fasta_path, monkeypatch
    ):
        monkeypatch.setattr(tools_mod, "GET_VARIANTS_EVIDENCE_MAX_SITES", 1)
        cfg = BAMCPConfig(reference=comprehensive_ref_fasta_path)
        p = _payload(
            await handle_get_variants({"file_path": comprehensive_bam_path, "region": REGION}, cfg)
        )
        pse = p["applied_filters"]["per_site_evidence"]
        assert pse["capped"] is True and pse["scored_sites"] == 1
        unscored = [v for v in p["variants"] if v.get("read_evidence_computed") is False]
        assert unscored and all(v["confidence"] == "low" for v in unscored)
        assert all("read_evidence_unavailable" in v["concerns"] for v in unscored)

    @pytest.mark.unit
    def test_zygosity_hints(self):
        assert zygosity_hint(0.5) == "consistent with heterozygous"
        assert zygosity_hint(0.97) == "consistent with homozygous alt"
        assert "allelic imbalance" in zygosity_hint(0.78)
        assert "below the germline het range" in zygosity_hint(0.12)
        assert zygosity_hint(None) is None

    @pytest.mark.unit
    def test_interpretation_text(self):
        v = {
            "contig": "chr17",
            "position": 43099914,
            "ref": "G",
            "alt": "A",
            "vaf": 0.52,
            "depth": 60,
            "alt_count": 31,
            "confidence": "high",
            "artifact_risk": {"artifact_likelihood": "low", "risks": []},
        }
        text = interpret_variant(v)
        assert text.startswith("chr17:43099914 G>A (BAMCP candidate)")
        assert "31/60 reads" in text and "heterozygous" in text and "confidence high" in text
        # 0-based payload dicts print 1-based.
        assert "chr17:43099915" in interpret_variant(v, one_based=False)
        assert annotate_interpretation(v)["concerns"] == []
        sv = interpret_variant({"contig": "chr1", "position": 5, "variant_kind": "sv"})
        assert "structural variant" in sv
        assert "not assessed" in interpret_variant({"contig": "c", "position": 1, "vaf": 0.5})


# -- low_base_quality (P1-6) ---------------------------------------------------


def _evidence(strand_bias: float, mean_q: float, quals_seen: bool = True) -> dict:
    return {
        "read_evidence_computed": True,
        "strand_bias": strand_bias,
        "mean_quality": mean_q,
        "quality_histogram": [3, 1, 0, 0, 0] if quals_seen else [0, 0, 0, 0, 0],
        "mapq_histogram": [],
        "near_end_fraction": 0.0,
    }


class TestLowBaseQuality:
    VARIANT = {"position": 10, "ref": "A", "alt": "T", "depth": 50, "vaf": 0.3}

    @pytest.mark.unit
    def test_low_bq_alone_is_a_concern(self):
        risk = compute_artifact_risk(self.VARIANT, _evidence(0.1, 8.3), None, 0)
        assert [r["type"] for r in risk["risks"]] == ["low_base_quality"]
        assert risk["artifact_likelihood"] == "medium"

    @pytest.mark.unit
    def test_low_bq_with_strand_bias_is_high(self):
        # Medium strand bias (0.2) + low BQ (0.3) = 0.5 would read medium; the combination
        # is the sequencing-error signature and is forced to HIGH.
        risk = compute_artifact_risk(self.VARIANT, _evidence(0.85, 12.0), None, 0)
        types = {r["type"] for r in risk["risks"]}
        assert types == {"low_base_quality", "strand_bias"}
        assert risk["artifact_likelihood"] == "high"
        # The prod example: Q8.3 support, 11F/0R.
        risk = compute_artifact_risk(self.VARIANT, _evidence(1.0, 8.3), None, 0)
        assert risk["artifact_likelihood"] == "high"

    @pytest.mark.unit
    def test_no_concern_without_observed_qualities_or_when_good(self):
        risk = compute_artifact_risk(self.VARIANT, _evidence(0.1, 0, quals_seen=False), None, 0)
        assert risk["risks"] == []
        risk = compute_artifact_risk(self.VARIANT, _evidence(0.1, 35.0), None, 0)
        assert risk["risks"] == []


# -- BA1/BS1 reliability (P2-11) ---------------------------------------------


class TestBenignFrequencyReliability:
    @pytest.mark.unit
    def test_ba1_bs1_suppressed_when_unreliable(self):
        crit = {
            c["code"]: c["evaluate_with"]
            for c in _frequency_criteria(
                {
                    "global_af": 0.2,
                    "max_pop_af": 0.3,
                    "frequency_reliable": False,
                    "warnings": ["low coverage"],
                }
            )
        }
        assert "Do NOT apply BA1/BS1" in crit["BA1"]
        assert "Do NOT apply BA1/BS1" in crit["BS1"]
        assert "low coverage" in crit["BA1"]

    @pytest.mark.unit
    def test_ba1_bs1_normal_when_reliable(self):
        crit = {
            c["code"]: c["evaluate_with"]
            for c in _frequency_criteria({"global_af": 0.2, "max_pop_af": 0.3})
        }
        assert "Do NOT" not in crit["BA1"] and "Do NOT" not in crit["BS1"]

    @pytest.mark.unit
    def test_failed_filters_also_suppress(self):
        crit = {
            c["code"]: c["evaluate_with"]
            for c in _frequency_criteria({"global_af": 0.2, "max_pop_af": 0.3, "filters": ["AC0"]})
        }
        assert "FAILED gnomAD QC" in crit["BA1"]


# -- Disallowed host error (P2-10) ---------------------------------------------


@pytest.mark.unit
def test_disallowed_host_error_lists_allowed_hosts():
    cfg = BAMCPConfig(
        allow_remote_files=True, allowed_remote_hosts=["ftp.ensembl.org", "ftp.ncbi.nlm.nih.gov"]
    )
    with pytest.raises(ValueError) as exc:
        validate_remote_url("https://evil.example/x.bam", cfg)
    msg = str(exc.value)
    assert "evil.example" in msg
    assert "ftp.ensembl.org" in msg and "ftp.ncbi.nlm.nih.gov" in msg


# -- Compact viewer text (P1-8) ------------------------------------------------


class TestViewerTextSize:
    @pytest.mark.unit
    @pytest.mark.asyncio
    async def test_visualize_text_is_compact(
        self, comprehensive_bam_path, comprehensive_ref_fasta_path
    ):
        cfg = BAMCPConfig(reference=comprehensive_ref_fasta_path)
        r = await handle_visualize_region(
            {"file_path": comprehensive_bam_path, "region": "chr1:1000-1300"}, cfg
        )
        text = r["content"][0]["text"]
        assert len(text) < VIEWER_SUMMARY_MAX_CHARS
        # Reads stay in the payload only.
        assert r["_meta"]["ui/init"]["reads"]["count"] > 0
        assert "sequence" not in text and "cigar" not in text
        assert "chr1:1051 A>T" in text  # the planted variant is summarized

    @pytest.mark.unit
    def test_summary_bound_with_many_long_variants(self):
        data = RegionData(
            contig="chr1", start=0, end=100, reads=[], coverage=[10] * 100, variants=[]
        )
        variants = [
            {"interpretation": "x" * 500, "confidence": "low", "alt_count": i} for i in range(200)
        ]
        filters = {
            "min_vaf": 0.1,
            "min_depth": 2,
            "min_mapq": 20,
            "min_baseq": 13,
            "reference_source": "auto",
        }
        text = _viewer_summary("Region chr1:1-100", data, variants, filters)
        assert len(text) <= VIEWER_SUMMARY_MAX_CHARS
        assert "195 more" in text
