"""Per-call resource bounds: multi-sample/huge VCFs, per-read payload, get_variants output.

The prod server is an open demo on a 2 GB container, so every one of these inputs (a
1000 Genomes cohort VCF, long noisy reads, a min_vaf=0 scan of 1 Mb) must stay bounded —
and every truncation must be stated in applied_filters so absence is never silent.
"""

from __future__ import annotations

import json

import pysam
import pytest

import bamcp.core.parsers as parsers
import bamcp.core.tools as tm
from bamcp.config import BAMCPConfig
from bamcp.core.parsers import RegionData, fetch_region, load_vcf_variants, vcf_truncated
from bamcp.core.serialization import serialize_region_data

REF_LEN = 2000


def _write_vcf(tmp_path, n_samples: int, positions: list[int], name: str = "calls.vcf.gz"):
    h = pysam.VariantHeader()
    h.add_line(f"##contig=<ID=chr1,length={REF_LEN}>")
    h.add_line('##INFO=<ID=AF,Number=A,Type=Float,Description="AF">')
    h.add_line('##INFO=<ID=DP,Number=1,Type=Integer,Description="Depth">')
    h.add_line('##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">')
    h.add_line('##FORMAT=<ID=AD,Number=R,Type=Integer,Description="Allele depths">')
    names = [f"S{i:04d}" for i in range(n_samples)]
    for s in names:
        h.add_sample(s)
    path = str(tmp_path / name)
    with pysam.VariantFile(path, "wz", header=h) as out:
        for pos in positions:
            rec = out.new_record(contig="chr1", start=pos - 1, alleles=("A", "G"))
            rec.info["AF"] = (0.5,)
            rec.info["DP"] = 20
            for s in names:
                rec.samples[s]["GT"] = (0, 1)
                rec.samples[s]["AD"] = (10, 10)
            out.write(rec)
    pysam.tabix_index(path, preset="vcf", force=True)
    return path, names


@pytest.fixture
def bam_and_ref(tmp_path):
    """A reference plus a BAM of 20 reads, each carrying 20 mismatches (every 5th base)."""
    # Deterministic, non-repetitive sequence (no homopolymer runs to trip the evidence rules).
    ref = "".join("ACGT"[(i * 7 + i // 4) % 4] for i in range(REF_LEN))
    ref_path = str(tmp_path / "ref.fa")
    with open(ref_path, "w") as f:
        f.write(">chr1\n" + ref + "\n")
    pysam.faidx(ref_path)
    header = {"HD": {"VN": "1.6", "SO": "coordinate"}, "SQ": [{"SN": "chr1", "LN": REF_LEN}]}
    bam_path = str(tmp_path / "reads.bam")
    swap = {"A": "C", "C": "G", "G": "T", "T": "A"}
    with pysam.AlignmentFile(bam_path, "wb", header=header) as out:
        for i in range(20):
            start = 100 + i
            seq = list(ref[start : start + 100])
            for j in range(0, 100, 5):
                seq[j] = swap[seq[j]]
            a = pysam.AlignedSegment()
            a.query_name = f"r{i}"
            a.query_sequence = "".join(seq)
            a.flag = 0
            a.reference_id = 0
            a.reference_start = start
            a.mapping_quality = 60
            a.cigartuples = [(0, 100)]
            a.query_qualities = pysam.qualitystring_to_array("I" * 100)
            out.write(a)
    pysam.index(bam_path)
    return bam_path, ref_path


# -- 1. VCF bounds ------------------------------------------------------------------------


@pytest.mark.unit
def test_few_samples_keep_full_detail(tmp_path):
    path, names = _write_vcf(tmp_path, 3, [100])
    (v,) = load_vcf_variants(path, "chr1:1-2000")
    assert v["sample_names"] == names
    assert v["samples"]["S0000"] == {"GT": [0, 1], "AD": [10, 10]}
    assert v["sample_count"] == 3
    assert v["samples_truncated"] is False


@pytest.mark.unit
def test_cohort_vcf_keeps_gt_of_first_samples_only(tmp_path):
    cap = parsers.VCF_MAX_SAMPLES_DETAIL
    path, names = _write_vcf(tmp_path, cap + 10, [100, 200])
    variants = load_vcf_variants(path, "chr1:1-2000")
    for v in variants:
        assert v["sample_count"] == cap + 10
        assert v["samples_truncated"] is True
        assert v["sample_names"] == names[:cap]
        assert set(v["samples"]) == set(names[:cap])
        # GT only (the shape the viewer/_format_genotypes read), not every FORMAT field.
        assert all(s == {"GT": [0, 1]} for s in v["samples"].values())
    assert tm._format_genotypes(variants[0]["samples"])[0] == "S0000=0/1"


@pytest.mark.unit
def test_vcf_record_cap(tmp_path, monkeypatch):
    path, _ = _write_vcf(tmp_path, 1, [100, 200, 300, 400, 500])
    monkeypatch.setattr(parsers, "VCF_MAX_RECORDS", 3)
    variants = load_vcf_variants(path, "chr1:1-2000")
    assert [v["position"] for v in variants] == [99, 199, 299]
    assert vcf_truncated(variants)
    # At or under the cap: not truncated.
    assert not vcf_truncated(load_vcf_variants(path, "chr1:1-350"))
    # A plain list (e.g. a stub) reads as not truncated.
    assert not vcf_truncated([])


@pytest.mark.unit
@pytest.mark.asyncio
@pytest.mark.parametrize("variant_source", ["vcf", "auto"])
async def test_vcf_truncation_reaches_applied_filters(
    tmp_path, bam_and_ref, monkeypatch, variant_source
):
    bam_path, ref_path = bam_and_ref
    vcf_path, _ = _write_vcf(tmp_path, 1, [150, 160, 170, 180])
    monkeypatch.setattr(parsers, "VCF_MAX_RECORDS", 2)
    monkeypatch.setattr(tm, "VCF_MAX_RECORDS", 2)
    config = BAMCPConfig(reference=ref_path)
    args = {
        "file_path": bam_path,
        "region": "chr1:101-300",
        "vcf_path": vcf_path,
        "variant_source": variant_source,
    }

    payload = json.loads((await tm.handle_get_variants(args, config))["content"][0]["text"])
    assert payload["applied_filters"]["vcf_truncated"] is True
    assert payload["applied_filters"]["vcf_max_records"] == 2
    vcf_sites = [v for v in payload["variants"] if v["source"] == "vcf"]
    assert len(vcf_sites) == 2

    viewer = await tm.handle_visualize_region(args, config)
    assert viewer["_meta"]["ui/init"]["applied_filters"]["vcf_truncated"] is True
    assert "VCF TRUNCATED" in viewer["content"][0]["text"]

    summary = await tm.handle_get_region_summary(args, config)
    assert "VCF TRUNCATED at 2 records" in summary["content"][0]["text"]


@pytest.mark.unit
@pytest.mark.asyncio
async def test_untruncated_vcf_states_nothing(tmp_path, bam_and_ref):
    bam_path, ref_path = bam_and_ref
    vcf_path, _ = _write_vcf(tmp_path, 1, [150])
    payload = json.loads(
        (
            await tm.handle_get_variants(
                {"file_path": bam_path, "region": "chr1:101-300", "vcf_path": vcf_path},
                BAMCPConfig(reference=ref_path),
            )
        )["content"][0]["text"]
    )
    assert "vcf_truncated" not in payload["applied_filters"]


@pytest.mark.unit
def test_region_bytes_count_sample_entries():
    def region(n_samples: int) -> RegionData:
        samples = {f"S{i}": {"GT": [0, 1], "AD": [1, 1]} for i in range(n_samples)}
        # Two ALTs of one record share one samples dict: counted once.
        variants = [
            {"position": 1, "ref": "A", "alt": alt, "samples": samples} for alt in ("C", "G")
        ]
        return RegionData("chr1", 0, 10, [], [], variants)

    base = tm._approx_region_bytes(region(0))
    many = tm._approx_region_bytes(region(1000))
    assert many - base >= 1000 * 300
    assert many - base < 2 * 1000 * 400  # shared dict not double-counted


# -- 2. Per-read payload ------------------------------------------------------------------


@pytest.mark.unit
def test_normal_reads_keep_every_mismatch(bam_and_ref):
    bam_path, ref_path = bam_and_ref
    data = fetch_region(bam_path, "chr1:1-2000", ref_path)
    assert len(data.reads) == 20
    assert all(len(r.mismatches) == 20 and not r.mismatches_truncated for r in data.reads)
    assert not data.reads_truncated
    assert "mismatches_truncated" not in serialize_region_data(data)["reads"]


@pytest.mark.unit
def test_mismatches_capped_per_read(bam_and_ref, monkeypatch):
    bam_path, ref_path = bam_and_ref
    monkeypatch.setattr(parsers, "READ_MAX_MISMATCHES", 5)
    data = fetch_region(bam_path, "chr1:1-2000", ref_path)
    for r in data.reads:
        assert len(r.mismatches) == 5
        assert r.mismatches_truncated
        # The first mismatches by reference position are the ones kept.
        assert [m["pos"] for m in r.mismatches] == [r.position + 5 * k for k in range(5)]
    cols = serialize_region_data(data)["reads"]
    assert cols["mismatches_truncated"] == [True] * 20


@pytest.mark.unit
def test_read_payload_budget_sets_reads_truncated(bam_and_ref, monkeypatch):
    bam_path, ref_path = bam_and_ref
    # Budget below one read's estimate: the first read is kept, then extraction stops.
    monkeypatch.setattr(parsers, "READ_PAYLOAD_MAX_BYTES", 100)
    data = fetch_region(bam_path, "chr1:1-2000", ref_path, min_vaf=0.01, min_depth=1)
    assert len(data.reads) == 1
    assert data.reads_truncated
    # Coverage still comes from the uncapped pass.
    assert max(data.coverage) == 20


# -- 3. get_variants / get_region_summary output cap --------------------------------------


@pytest.mark.unit
def test_cap_variants_priority_and_order():
    vs = [
        {"position": 5, "ref": "A", "alt": "C", "source": "bamcp", "alt_count": 3, "vaf": 0.1},
        {"position": 1, "ref": "A", "alt": "C", "source": "bamcp", "alt_count": 9, "vaf": 0.5},
        {"position": 9, "ref": "A", "alt": "C", "source": "vcf", "alt_count": 0, "vaf": 0.0},
        {"position": 3, "ref": "A", "alt": "C", "source": "bamcp", "alt_count": 9, "vaf": 0.2},
    ]
    assert tm._cap_variants(vs, 10) is vs
    # VCF first, then most alt reads (VAF breaks the 9/9 tie); kept in input order.
    assert [v["position"] for v in tm._cap_variants(vs, 2)] == [1, 9]
    assert [v["position"] for v in tm._cap_variants(vs, 3)] == [1, 9, 3]


@pytest.mark.unit
@pytest.mark.asyncio
async def test_get_variants_output_capped(bam_and_ref, monkeypatch):
    bam_path, ref_path = bam_and_ref
    config = BAMCPConfig(reference=ref_path)
    args = {
        "file_path": bam_path,
        "region": "chr1:101-300",
        "min_vaf": 0,
        "min_depth": 1,
        "min_alt_count": 1,
    }
    full = json.loads((await tm.handle_get_variants(args, config))["content"][0]["text"])
    total = full["count"]
    assert total > 5
    assert full["applied_filters"]["truncated"] is False
    assert full["applied_filters"]["total_candidates"] == total
    assert full["applied_filters"]["max_results"] == tm.GET_VARIANTS_MAX_RESULTS

    monkeypatch.setattr(tm, "GET_VARIANTS_MAX_RESULTS", 5)
    tm._services_registry.clear()
    capped = json.loads((await tm.handle_get_variants(args, config))["content"][0]["text"])
    filters = capped["applied_filters"]
    assert capped["count"] == 5
    assert filters["truncated"] is True
    assert filters["total_candidates"] == total
    assert filters["max_results"] == 5
    assert "truncation_note" in filters
    positions = [v["position"] for v in capped["variants"]]
    assert positions == sorted(positions)
    # The kept ones are the most alt-supported candidates.
    top = sorted((v["alt_count"] for v in full["variants"]), reverse=True)[:5]
    assert min(v["alt_count"] for v in capped["variants"]) >= min(top)


@pytest.mark.unit
@pytest.mark.asyncio
async def test_region_summary_listing_capped(bam_and_ref, monkeypatch):
    bam_path, ref_path = bam_and_ref
    args = {"file_path": bam_path, "region": "chr1:101-300", "min_alt_count": 1}
    config = BAMCPConfig(reference=ref_path, min_vaf=0.0, min_depth=1)
    monkeypatch.setattr(tm, "GET_VARIANTS_MAX_RESULTS", 3)
    text = (await tm.handle_get_region_summary(args, config))["content"][0]["text"]
    assert "Variant listing TRUNCATED: 3 of" in text
    assert sum(1 for line in text.splitlines() if line.startswith("  chr1:")) == 3


@pytest.mark.unit
@pytest.mark.asyncio
@pytest.mark.parametrize("tool", ["visualize_region", "jump_to"])
async def test_viewer_variants_capped(bam_and_ref, monkeypatch, tool):
    """The viewer payload's variant list is capped like get_variants, and says so."""
    bam_path, ref_path = bam_and_ref
    config = BAMCPConfig(reference=ref_path, min_alt_count=1)
    if tool == "visualize_region":
        handler, args = tm.handle_visualize_region, {"region": "chr1:101-300"}
    else:
        handler, args = tm.handle_jump_to, {"contig": "chr1", "position": 200, "window": 200}
    args.update(file_path=bam_path, min_vaf=0, min_depth=1)

    full = (await handler(args, config))["_meta"]["ui/init"]
    total = len(full["variants"])
    assert total > 3
    assert full["applied_filters"]["variants_truncated"] is False
    assert full["applied_filters"]["total_candidates"] == total

    monkeypatch.setattr(tm, "GET_VARIANTS_MAX_RESULTS", 3)
    tm._services_registry.clear()
    result = await handler(args, config)
    payload = result["_meta"]["ui/init"]
    filters = payload["applied_filters"]
    assert len(payload["variants"]) == 3
    assert filters["variants_truncated"] is True
    assert filters["total_candidates"] == total
    assert filters["max_results"] == 3
    assert "truncation_note" in filters
    assert f"variant list TRUNCATED to 3 of {total}" in result["content"][0]["text"]
