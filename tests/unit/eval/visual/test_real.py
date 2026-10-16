"""Real-data tier: truth classification, candidate scanning, blinding, and the run set.

Hermetic: a tiny local "GIAB" (reference, BAM with a het SNV, a noisy non-variant site and
duplicate reads, a bgzipped truth VCF, BEDs) stands in for the public files.
"""

from __future__ import annotations

import gzip
import json
import random
from pathlib import Path

import pysam
import pytest

from bamcp.eval.visual import real
from bamcp.eval.visual.grade import grade_run_set
from bamcp.eval.visual.leaklint import lint_run_set

pytestmark = pytest.mark.unit

CHROM = "chr20"
LEN = 12_000
HET = 3000  # 0-based truth het SNV
NOISY = 6000  # 0-based: 4 of ~40 reads carry a non-ref base, no truth variant
DUP_SITE = 9000  # 0-based: alt carried only by duplicate-flagged reads


def _ref() -> str:
    rng = random.Random(5)  # noqa: S311
    seq = []
    while len(seq) < LEN:
        b = rng.choice("ACGT")
        if len(seq) >= 2 and seq[-1] == seq[-2] == b:
            continue
        seq.append(b)
    return "".join(seq)


def _other(b: str) -> str:
    return {"A": "C", "C": "G", "G": "T", "T": "A"}[b]


@pytest.fixture(scope="module")
def giab(tmp_path_factory: pytest.TempPathFactory) -> dict[str, Path]:
    d = tmp_path_factory.mktemp("giab")
    ref = _ref()
    fa = d / "ref.fa"
    fa.write_text(f">{CHROM}\n{ref}\n")
    pysam.faidx(str(fa))

    header = {"HD": {"VN": "1.6", "SO": "coordinate"}, "SQ": [{"SN": CHROM, "LN": LEN}]}
    rows = []
    rng = random.Random(9)  # noqa: S311
    n = 0
    for site, kind in ((HET, "het"), (NOISY, "noisy"), (DUP_SITE, "dup")):
        for i in range(40):
            start = site - 140 + rng.randint(0, 130)
            seq = list(ref[start : start + 150])
            flag = 16 if i % 2 else 0
            if (kind == "het" and i % 2 == 0) or (kind == "noisy" and i < 4):
                seq[site - start] = _other(ref[site])
            if kind == "dup" and i < 6:
                seq[site - start] = _other(ref[site])
                flag |= 1024
            rows.append((start, f"read{n}", flag, "".join(seq)))
            n += 1
    rows.sort()
    bam = d / "sample.bam"
    with pysam.AlignmentFile(str(bam), "wb", header=header) as out:
        for start, name, flag, seq in rows:
            a = pysam.AlignedSegment(out.header)
            a.query_name = name
            a.query_sequence = seq
            a.flag = flag
            a.reference_id = 0
            a.reference_start = start
            a.mapping_quality = 60
            a.cigartuples = [(0, 150)]
            a.query_qualities = pysam.qualitystring_to_array("I" * 150)
            out.write(a)
    pysam.index(str(bam))

    vcf_txt = d / "truth.vcf"
    vcf_txt.write_text(
        "##fileformat=VCFv4.2\n"
        f"##contig=<ID={CHROM},length={LEN}>\n"
        '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n'
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tHG002\n"
        f"{CHROM}\t{HET + 1}\t.\t{ref[HET]}\t{_other(ref[HET])}\t50\tPASS\t.\tGT\t0|1\n"
    )
    vcf = d / "truth.vcf.gz"
    pysam.tabix_compress(str(vcf_txt), str(vcf), force=True)
    pysam.tabix_index(str(vcf), preset="vcf", force=True)

    conf = d / "conf.bed"
    conf.write_text(f"{CHROM}\t1000\t11000\n")
    rep = d / "rep.bed.gz"
    with gzip.open(rep, "wt") as fh:
        fh.write(f"{CHROM}\t200\t300\n")
    return {"ref": fa, "bam": bam, "vcf": vcf, "conf": conf, "rep": rep}


# --- truth classification ------------------------------------------------------------


@pytest.mark.parametrize(
    ("ref", "alt", "gt", "expect"),
    [
        ("A", "G", (0, 1), ("SNV", "G", "het", 0)),
        ("A", "G", (1, 1), ("SNV", "G", "hom_alt", 0)),
        ("ACT", "A", (0, 1), ("del", "CT", "het", 1)),
        ("A", "AGG", (1, 1), ("ins", "GG", "hom_alt", 0)),
    ],
)
def test_classify_record(ref, alt, gt, expect):
    got = real.classify_record(ref, alt, gt)
    lab = got["label"]
    assert (lab["variant_type"], lab["alt_allele"], lab["zygosity"], got["offset"]) == expect
    assert lab["artifact_verdict"] == "real" and lab["artifact_type"] == "none"


@pytest.mark.parametrize(
    ("ref", "alt", "gt"),
    [
        ("A", "G", (1, 2)),  # multi-allelic genotype
        ("A", "G", (0, None)),  # missing call
        ("AC", "GT", (0, 1)),  # MNV / complex
        ("A" + "C" * 11, "A", (0, 1)),  # deletion over MAX_INDEL
        ("A", "G", (0, 0)),  # hom-ref
    ],
)
def test_classify_record_rejects(ref, alt, gt):
    assert real.classify_record(ref, alt, gt) is None


def test_intervals_merge_and_lookup():
    iv = real.Intervals([(10, 20), (15, 30), (40, 50)])
    assert len(iv) == 2
    assert iv.contains(10) and iv.contains(29) and not iv.contains(30)
    assert iv.overlaps(28, 45) and not iv.overlaps(30, 40)
    assert iv.covers(12, 30) and not iv.covers(25, 45)


def test_load_bed_plain_and_gz(giab):
    assert real.load_bed(giab["conf"], CHROM).covers(1000, 11000)
    assert real.load_bed(giab["rep"], CHROM).contains(250)
    assert len(real.load_bed(giab["conf"], "chr1")) == 0


# --- scanning --------------------------------------------------------------------------


def test_scan_window_finds_each_kind(giab):
    conf = real.load_bed(giab["conf"], CHROM)
    rep = real.load_bed(giab["rep"], CHROM)
    with (
        pysam.AlignmentFile(str(giab["bam"])) as bam,
        pysam.FastaFile(str(giab["ref"])) as fa,
        pysam.VariantFile(str(giab["vcf"])) as vcf,
    ):
        pools = real.scan_window(bam, fa, vcf, CHROM, 1000, 11000, conf, rep)
    het = pools["real_snv_het"]
    assert [s.site0 for s in het] == [HET]
    assert het[0].label["zygosity"] == "het" and het[0].vcf_pos == HET + 1
    noisy = {s.site0: s for s in pools["real_noisy_negative"]}
    assert NOISY in noisy and noisy[NOISY].alt_reads == 4
    assert noisy[NOISY].label["artifact_type"] == "unknown"
    # Duplicate-flagged alt reads are not evidence: that site is not a noisy candidate,
    # and the truth het's neighbourhood never yields negatives.
    assert DUP_SITE not in noisy
    assert all(abs(s - HET) > real.NEG_CLEAR_BP for s in noisy)
    clean = pools["real_clean_negative"]
    assert clean and all(s.label["artifact_type"] == "none" for s in clean)


def test_select_sites_draws_quotas(giab):
    quotas = {"real_snv_het": 1, "real_noisy_negative": 1, "real_clean_negative": 1}
    sites = real.select_sites(
        3,
        chrom=CHROM,
        bam_path=str(giab["bam"]),
        ref_path=str(giab["ref"]),
        vcf_path=str(giab["vcf"]),
        confident_bed=giab["conf"],
        repeat_bed=giab["rep"],
        quotas=quotas,
        max_windows=40,
        pool_factor=1,
        log=lambda *_: None,
    )
    assert sorted(s.case_type for s in sites) == sorted(quotas)
    again = real.select_sites(
        3,
        chrom=CHROM,
        bam_path=str(giab["bam"]),
        ref_path=str(giab["ref"]),
        vcf_path=str(giab["vcf"]),
        confident_bed=giab["conf"],
        repeat_bed=giab["rep"],
        quotas=quotas,
        max_windows=40,
        pool_factor=1,
        log=lambda *_: None,
    )
    assert [s.site0 for s in sites] == [s.site0 for s in again]  # deterministic


def test_select_sites_fails_loudly_when_a_stratum_is_short(giab):
    with pytest.raises(ValueError, match="real_snv_hom"):
        real.select_sites(
            3,
            chrom=CHROM,
            bam_path=str(giab["bam"]),
            ref_path=str(giab["ref"]),
            vcf_path=str(giab["vcf"]),
            confident_bed=giab["conf"],
            repeat_bed=giab["rep"],
            quotas={"real_snv_hom": 1},
            max_windows=5,
            log=lambda *_: None,
        )


# --- blinding + run set ------------------------------------------------------------------


def _sites(giab) -> list[real.RealSite]:
    het_label = real.classify_record("A", "C", (0, 1))["label"]
    with pysam.FastaFile(str(giab["ref"])) as fa:
        r = fa.fetch(CHROM, HET, HET + 1)
    het_label = {**het_label, "alt_allele": _other(r)}
    return [
        real.RealSite("real_snv_het", CHROM, HET, het_label, False, vcf_pos=HET + 1),
        real.RealSite(
            "real_noisy_negative",
            CHROM,
            NOISY,
            {**real.NEGATIVE_LABEL, "artifact_type": "unknown"},
            False,
            alt_reads=4,
            depth=40,
        ),
        real.RealSite(
            "real_clean_negative",
            CHROM,
            DUP_SITE,
            {**real.NEGATIVE_LABEL, "artifact_type": "none"},
            False,
        ),
    ]


def _counts(bam_path: str, contig: str, pos0: int) -> tuple[int, int]:
    with pysam.AlignmentFile(bam_path) as bam:
        for col in bam.pileup(contig, pos0, pos0 + 1, truncate=True, stepper="nofilter"):
            bases = [
                p.alignment.query_sequence[p.query_position]
                for p in col.pileups
                if p.query_position is not None
            ]
            return len(bases), len(set(bases))
    return 0, 0


def test_write_real_fixtures_blinds_and_preserves_evidence(giab, tmp_path):
    sites = _sites(giab)
    paths, items, mapping = real.write_real_fixtures(
        sites, tmp_path, 3, bam_path=str(giab["bam"]), ref_path=str(giab["ref"])
    )
    by_type = {it.case_type: it for it in items}
    het = by_type["real_snv_het"]
    assert het.contig.startswith("chr") and het.contig != CHROM
    depth, n_alleles = _counts(str(paths["bam"]), het.contig, het.site)
    assert depth >= 30 and n_alleles == 2
    with pysam.FastaFile(str(paths["reference"])) as new, pysam.FastaFile(str(giab["ref"])) as old:
        assert new.fetch(het.contig, het.site - 20, het.site + 20) == old.fetch(
            CHROM, HET - 20, HET + 20
        )
    with pysam.AlignmentFile(str(paths["bam"])) as bam:
        reads = list(bam.fetch(until_eof=True))
    assert not any(r.query_name.startswith("read") for r in reads)  # names hashed
    assert not any(r.is_duplicate for r in reads)
    assert all(r.get_tags() == [] for r in reads)
    # The duplicate-only alt signal is gone: the "clean" site shows one allele.
    clean = by_type["real_clean_negative"]
    assert _counts(str(paths["bam"]), clean.contig, clean.site)[1] == 1
    assert {m["source"] for m in mapping} == {f"{CHROM}:{s.site0 + 1}" for s in sites}


def test_build_real_run_set_lints_and_grades(giab, tmp_path):
    sites_path = tmp_path / "sites.jsonl"
    digest = real.sites_to_jsonl(_sites(giab), sites_path)
    assert len(digest) == 64
    assert [s.site0 for s in real.sites_from_jsonl(sites_path)] == [HET, NOISY, DUP_SITE]
    out = tmp_path / "run"
    res = real.build_real_run_set(
        out,
        sites_path,
        seed=3,
        conditions=["Traw", "Tjson"],
        models=["m"],
        png=False,
        bam_path=str(giab["bam"]),
        ref_path=str(giab["ref"]),
        log=lambda *_: None,
    )
    assert res.n_items == 3 and res.n_runs == 6
    report = lint_run_set(out)
    assert report.ok, report.problems
    meta = json.loads((out / "build.json").read_text())
    assert meta["tier"] == "real" and meta["sites_sha256"] == digest
    caption = next((out / "stim").iterdir()) / "caption.md"
    assert real.SAMPLE_LINE in caption.read_text()
    assert (out / "labels" / "mapping.jsonl").exists()

    labels = {
        json.loads(x)["item"]: json.loads(x)
        for x in (out / "labels" / "labels.jsonl").read_text().splitlines()
    }
    for run in (json.loads(x) for x in (out / "runs.jsonl").read_text().splitlines()):
        lab = labels[run["item"]]
        ans = {k: lab[k] for k in real.NEGATIVE_LABEL} | {
            "variant_type": lab["variant_type"],
            "alt_allele": lab["alt_allele"],
            "zygosity": lab["zygosity"],
            "artifact_type": "low_baseq" if lab["artifact_type"] == "unknown" else "none",
            "confidence": 0.9,
            "rationale": "oracle",
        }
        (out / run["answer_path"]).write_text(json.dumps(ans))
    summary = grade_run_set(out, charts=False)
    g = summary["groups"]["m|Traw"]
    assert g["primary_joint_accuracy"] == 1.0
    # The negative's artifact_type is "unknown" and is not scored, so an arbitrary type
    # answer there costs nothing.
    assert g["field_accuracy"]["artifact_type"] == 1.0
    assert "unknown" not in g["artifact_type_confusion"]
    # ...and the per-run export leaves that cell empty rather than marking it wrong.
    import csv

    rows = list(csv.DictReader((out / "results" / "results.csv").open()))
    noisy = [r for r in rows if r["case_type"] == "real_noisy_negative"]
    assert noisy and all(r["artifact_type_correct"] == "" for r in noisy)
    assert all(r["artifact_type_correct"] == "1" for r in rows if r not in noisy)


def _vcf_with(tmp_path: Path, giab, records: list[tuple[int, str, str]]) -> Path:
    txt = tmp_path / "t.vcf"
    body = "".join(f"{CHROM}\t{p}\t.\t{r}\t{a}\t50\tPASS\t.\tGT\t0/1\n" for p, r, a in records)
    txt.write_text(
        "##fileformat=VCFv4.2\n"
        f"##contig=<ID={CHROM},length={LEN}>\n"
        '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n'
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tHG002\n" + body
    )
    out = tmp_path / "t.vcf.gz"
    pysam.tabix_compress(str(txt), str(out), force=True)
    pysam.tabix_index(str(out), preset="vcf", force=True)
    return out


def test_truth_deletion_span_excludes_negatives_inside_it(giab, tmp_path):
    """A long truth deletion anchored > NEG_CLEAR_BP upstream still covers the noisy site:
    reads there reflect the true variant, so the site must not become a 'negative'."""
    with pysam.FastaFile(str(giab["ref"])) as fa:
        anchor0 = NOISY - 25
        ref = fa.fetch(CHROM, anchor0, anchor0 + 41)  # deletes NOISY-24 .. NOISY+15
    vcf = _vcf_with(tmp_path, giab, [(anchor0 + 1, ref, ref[0])])
    conf = real.load_bed(giab["conf"], CHROM)
    rep = real.load_bed(giab["rep"], CHROM)
    with (
        pysam.AlignmentFile(str(giab["bam"])) as bam,
        pysam.FastaFile(str(giab["ref"])) as fa,
        pysam.VariantFile(str(vcf)) as v,
    ):
        pools = real.scan_window(bam, fa, v, CHROM, 1000, 11000, conf, rep)
    assert NOISY not in {s.site0 for s in pools["real_noisy_negative"]}
    lo, hi = anchor0 - real.NEG_CLEAR_BP, anchor0 + 41 + real.NEG_CLEAR_BP
    for k in ("real_noisy_negative", "real_clean_negative"):
        assert not [s.site0 for s in pools[k] if lo <= s.site0 < hi], k
    # ...while sites just outside the cleared span are still eligible.
    assert any(s.site0 == lo - 1 for s in pools["real_clean_negative"])


def test_noisy_negatives_must_split_across_repeat_context(giab):
    """No noisy candidate sits in the repeat BED, so a half-in-repeats quota cannot be met
    and the selection fails instead of filling from outside repeats."""
    with pytest.raises(ValueError, match="inside repeats"):
        real.select_sites(
            3,
            chrom=CHROM,
            bam_path=str(giab["bam"]),
            ref_path=str(giab["ref"]),
            vcf_path=str(giab["vcf"]),
            confident_bed=giab["conf"],
            repeat_bed=giab["rep"],
            quotas={"real_noisy_negative": 2},
            max_windows=40,
            pool_factor=1,
            log=lambda *_: None,
        )


def test_snv_next_to_a_repeat_is_not_in_it(giab, tmp_path):
    """An SNV's repeat context is its own base: a repeat starting one base later does not
    move it into real_snv_in_repeat."""
    rep_path = tmp_path / "adjacent.bed.gz"
    with gzip.open(rep_path, "wt") as fh:
        fh.write(f"{CHROM}\t{HET + 1}\t{HET + 20}\n")
    conf = real.load_bed(giab["conf"], CHROM)
    rep = real.load_bed(rep_path, CHROM)
    with (
        pysam.AlignmentFile(str(giab["bam"])) as bam,
        pysam.FastaFile(str(giab["ref"])) as fa,
        pysam.VariantFile(str(giab["vcf"])) as vcf,
    ):
        pools = real.scan_window(bam, fa, vcf, CHROM, 1000, 11000, conf, rep)
    assert [s.site0 for s in pools["real_snv_het"]] == [HET]
    assert pools["real_snv_in_repeat"] == []


def _sites_at(case: str, in_repeat: bool, n: int, start: int = 0) -> list[real.RealSite]:
    """``n`` sites MIN_SPACING apart."""
    return [
        real.RealSite(case, CHROM, start + i * real.MIN_SPACING, {}, in_repeat) for i in range(n)
    ]


def test_pools_full_needs_both_noisy_subpools():
    quotas = {"real_noisy_negative": 4, "real_clean_negative": 1}
    clean = _sites_at("real_clean_negative", False, 1, start=10**6)
    lopsided = {"real_noisy_negative": _sites_at("n", False, 8), "real_clean_negative": clean}
    assert not real._pools_full(lopsided, quotas, 2)  # 8 >= 2x4 overall, but none in repeats
    split = {
        "real_noisy_negative": _sites_at("n", True, 4) + _sites_at("n", False, 4, start=10**5),
        "real_clean_negative": clean,
    }
    assert not real._pools_full(split, quotas, 2)  # clean pool short of 2x1
    split["real_clean_negative"] = _sites_at("real_clean_negative", False, 2, start=10**6)
    assert real._pools_full(split, quotas, 2)


def test_pools_full_counts_adjacent_candidates_once():
    """Ten adjacent clean positions from one locus can supply only one spaced site."""
    adjacent = [real.RealSite("real_clean_negative", CHROM, 5000 + i, {}, False) for i in range(10)]
    assert real._spaced(adjacent) == 1
    assert not real._pools_full({"real_clean_negative": adjacent}, {"real_clean_negative": 2}, 1)


def test_selftest_oracle_passes_on_a_real_run_set(giab, tmp_path):
    from bamcp.eval.visual.selftest import run_selftests

    sites_path = tmp_path / "sites.jsonl"
    real.sites_to_jsonl(_sites(giab), sites_path)
    out = tmp_path / "run"
    real.build_real_run_set(
        out,
        sites_path,
        seed=3,
        conditions=["Traw", "Tjson"],
        models=["m"],
        png=False,
        bam_path=str(giab["bam"]),
        ref_path=str(giab["ref"]),
        log=lambda *_: None,
    )
    results = {r.name: r for r in run_selftests(out, charts=False)}
    assert results["oracle = 100%"].ok, results["oracle = 100%"].detail


def test_pick_honours_zero_and_small_quotas():
    """An odd noisy quota asks for 0 sites from one side: none may be taken (or reserved)."""
    rng = random.Random(1)  # noqa: S311
    pool = [real.RealSite("n", CHROM, i * 5000, {}, True) for i in range(5)]
    taken: list[real.RealSite] = []
    assert real._pick(list(pool), 0, rng, taken) == [] and taken == []
    got = real._pick(list(pool), 2, rng, taken)
    assert len(got) == 2 and taken == got


@pytest.mark.parametrize(("flag", "expect"), [([], None), (["--png"], True), (["--no-png"], False)])
def test_build_real_cli_forwards_png_mode(monkeypatch, tmp_path, flag, expect):
    """A registered real-tier build can require the PNG step, like the synthetic `build`."""
    from types import SimpleNamespace

    from bamcp.eval.visual import cli

    seen: dict = {}

    def fake_build(out, sites, **kw):
        seen.update(kw)
        return SimpleNamespace(n_stimuli=0, n_items=0, n_runs=0, out=tmp_path)

    monkeypatch.setattr(real, "build_real_run_set", fake_build)
    monkeypatch.setattr(cli, "_lint", lambda *_: True)
    argv = ["build-real", "--out", str(tmp_path), "--sites", "s.jsonl", "--seed", "3", *flag]
    assert cli.main(argv) == 0
    assert seen["png"] is expect


def test_same_pos_truth_records_block_each_other(giab, tmp_path):
    """Two truth records at one POS are not 'isolated' just because they share a start."""
    with pysam.FastaFile(str(giab["ref"])) as fa:
        r = fa.fetch(CHROM, HET, HET + 2)
    vcf = _vcf_with(
        tmp_path, giab, [(HET + 1, r[0], _other(r[0])), (HET + 1, r, r[0])]
    )  # an SNV and a 1 bp deletion anchored at the same base
    conf = real.load_bed(giab["conf"], CHROM)
    rep = real.load_bed(giab["rep"], CHROM)
    with (
        pysam.AlignmentFile(str(giab["bam"])) as bam,
        pysam.FastaFile(str(giab["ref"])) as fa,
        pysam.VariantFile(str(vcf)) as v,
    ):
        pools = real.scan_window(bam, fa, v, CHROM, 1000, 11000, conf, rep)
    assert pools["real_snv_het"] == [] and pools["real_indel_het"] == []


def test_rejected_window_draws_do_not_count_as_scans(giab):
    """With room for only a few non-overlapping windows, max_windows still means scans."""
    lines: list[str] = []
    with pytest.raises(ValueError):  # the fixture has no hom SNV, so the pools never fill
        real.select_sites(
            3,
            chrom=CHROM,
            bam_path=str(giab["bam"]),
            ref_path=str(giab["ref"]),
            vcf_path=str(giab["vcf"]),
            confident_bed=giab["conf"],
            repeat_bed=giab["rep"],
            quotas={"real_snv_hom": 1},
            max_windows=3,
            pool_factor=1,
            log=lines.append,
        )
    summary = lines[-1]
    assert summary.startswith("scanned 3 windows in ") and "draws" in summary
