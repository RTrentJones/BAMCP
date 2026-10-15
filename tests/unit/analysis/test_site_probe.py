"""The requested site is evaluated even when it fails the calling filters.

The demo artifact (HG002 chr17:43096151 T>G) is low-base-quality support on one strand: the
default min_baseq=13 drops most of its alt bases, so get_variants doesn't call it — and the
curation/classify tools used to answer "not found", hiding exactly the artifact a curator asked
about. They now probe the site (min_baseq=0, no VAF/alt floor), score its evidence, and report
passes_filters=false with the filters it fails.
"""

from __future__ import annotations

import json

import pysam
import pytest

from bamcp.analysis.curation import handle_get_variant_curation_summary
from bamcp.analysis.fusion import assemble_evidence
from bamcp.config import BAMCPConfig
from bamcp.core.tools import handle_get_variants

SEQ = "ACGT" * 250  # 1000 bp, no homopolymers
SITE0 = 500  # 0-based; 1-based 501
REF_BASE = SEQ[SITE0]
ALT_BASE = "C" if REF_BASE != "C" else "G"


def _write(tmp_path, alt_reads: int, alt_qual: int) -> tuple[str, str]:
    ref = str(tmp_path / "ref.fa")
    with open(ref, "w") as fh:
        fh.write(">chr1\n" + SEQ + "\n")
    pysam.faidx(ref)
    header = {"HD": {"VN": "1.6", "SO": "coordinate"}, "SQ": [{"SN": "chr1", "LN": len(SEQ)}]}
    reads = []
    for i in range(20):  # reference reads, both strands, Q40
        reads.append((470 + i % 20, 16 if i % 2 else 0, None))
    for i in range(alt_reads):  # alt reads: forward strand only, low-BQ alt base
        reads.append((472 + i * 2, 0, alt_qual))
    bam = str(tmp_path / "probe.bam")
    with pysam.AlignmentFile(bam, "wb", header=header) as out:
        for n, (start, flag, aq) in enumerate(sorted(reads, key=lambda r: r[0])):
            seq = list(SEQ[start : start + 60])
            quals = [40] * 60
            if aq is not None:
                seq[SITE0 - start] = ALT_BASE
                quals[SITE0 - start] = aq
            a = pysam.AlignedSegment()
            a.query_name = f"r{n}"
            a.query_sequence = "".join(seq)
            a.flag = flag
            a.reference_id = 0
            a.reference_start = start
            a.mapping_quality = 60
            a.cigartuples = [(0, 60)]
            a.query_qualities = pysam.qualitystring_to_array("".join(chr(q + 33) for q in quals))
            out.write(a)
    pysam.index(bam)
    return bam, ref


def _cfg(ref: str) -> BAMCPConfig:
    # The shipped defaults: min_baseq=13, min_alt_count=3, min_vaf=0.1.
    return BAMCPConfig(reference=ref, clinvar_enabled=False, gnomad_enabled=False)


def _args(bam: str, **kw) -> dict:
    return {
        "file_path": bam,
        "chrom": "chr1",
        "pos": SITE0 + 1,
        "ref": REF_BASE,
        "alt": ALT_BASE,
        **kw,
    }


@pytest.mark.unit
async def test_low_bq_one_strand_site_is_not_called_by_get_variants(tmp_path):
    bam, ref = _write(tmp_path, alt_reads=8, alt_qual=6)
    out = await handle_get_variants({"file_path": bam, "region": "chr1:480-520"}, _cfg(ref))
    assert json.loads(out["content"][0]["text"])["variants"] == []


@pytest.mark.unit
async def test_curation_scores_a_site_that_fails_the_filters(tmp_path):
    bam, ref = _write(tmp_path, alt_reads=8, alt_qual=6)
    out = await handle_get_variant_curation_summary(_args(bam, format="rubric"), _cfg(ref))
    doc = json.loads(out["content"][0]["text"])
    assert "error" not in doc, doc
    risk = doc["artifact_assessment"]
    assert risk["artifact_likelihood"] == "high"
    concerns = {r["type"] for r in risk["risks"]}
    assert {"low_base_quality", "strand_bias"} <= concerns
    check = doc["filter_check"]
    assert check["passes_filters"] is False
    assert check["failed_filters"][0] == "min_baseq"
    assert "min_alt_count" in check["failed_filters"]
    assert check["unfiltered"]["alt_count"] == 8
    assert check["at_filters"]["alt_count"] == 0
    assert check["low_baseq_alt_reads_removed"] == 8
    assert doc["variant"]["alt_count"] == 8
    assert doc["applied_filters"]["min_alt_count"] == 3
    assert doc["applied_filters"]["reference_source"] == "config"
    assert doc["confidence"] == "low"

    text = (await handle_get_variant_curation_summary(_args(bam), _cfg(ref)))["content"][0]["text"]
    assert "ARTIFACT RISK: HIGH" in text
    assert "FAILS the applied filters: min_baseq" in text


@pytest.mark.unit
async def test_curation_passing_site_reports_passes_filters(tmp_path):
    bam, ref = _write(tmp_path, alt_reads=8, alt_qual=40)
    out = await handle_get_variant_curation_summary(_args(bam, format="json"), _cfg(ref))
    doc = json.loads(out["content"][0]["text"])  # format="json" is an alias of rubric
    assert doc["filter_check"]["passes_filters"] is True
    assert doc["filter_check"]["failed_filters"] == []


@pytest.mark.unit
async def test_curation_no_support_is_still_not_found(tmp_path):
    bam, ref = _write(tmp_path, alt_reads=0, alt_qual=6)
    out = await handle_get_variant_curation_summary(_args(bam, format="rubric"), _cfg(ref))
    doc = json.loads(out["content"][0]["text"])
    assert "not found" in doc["error"]
    assert doc["applied_filters"]["min_alt_count"] == 3
    assert doc["applied_filters"]["reference_source"] == "config"


@pytest.mark.unit
async def test_classify_observation_reports_below_filters_with_evidence(tmp_path):
    bam, ref = _write(tmp_path, alt_reads=8, alt_qual=6)
    ev = await assemble_evidence("chr1", SITE0 + 1, REF_BASE, ALT_BASE, bam, _cfg(ref))
    obs = ev["observation"]
    assert obs["status"] == "below_filters"
    assert obs["passes_filters"] is False
    assert obs["failed_filters"][0] == "min_baseq"
    assert obs["artifact_likelihood"] == "high"
    assert {"low_base_quality", "strand_bias"} <= set(obs["artifact_concerns"])
    assert obs["alt_reads_seen"] == 8
