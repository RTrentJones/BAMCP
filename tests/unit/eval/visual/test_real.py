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


def test_nearby_sites_never_share_reads(giab, tmp_path):
    """Two sites closer than their read windows plus a read length: each read is copied into
    one item only, and the second item records how many it gave up."""
    het_label = real.classify_record("A", "C", (0, 1))["label"]
    sites = [
        real.RealSite("real_snv_het", CHROM, HET, het_label, False),
        real.RealSite(
            "real_clean_negative",
            CHROM,
            HET + 500,
            {**real.NEGATIVE_LABEL, "artifact_type": "none"},
            False,
        ),
    ]
    paths, items, mapping = real.write_real_fixtures(
        sites, tmp_path, 3, bam_path=str(giab["bam"]), ref_path=str(giab["ref"])
    )
    with pysam.AlignmentFile(str(giab["bam"])) as src:
        near = {
            r.query_name
            for r in src.fetch(CHROM, HET - real.READ_HALF_SPAN, HET + 500 + real.READ_HALF_SPAN)
            if not r.is_duplicate
        }
    with pysam.AlignmentFile(str(paths["bam"])) as out:
        written = sum(1 for _ in out.fetch(until_eof=True))
    assert written == len(near)  # every source read once, none twice
    assert mapping[0]["shared_reads_dropped"] == 0
    assert mapping[1]["shared_reads_dropped"] > 0


# =========================================================================================
# Round 4: features, paired fixtures, generators, lazy checks, build
# =========================================================================================

from collections import Counter  # noqa: E402

from bamcp.eval.visual import strata_r4  # noqa: E402

from . import r4_world as W  # noqa: E402


@pytest.fixture(scope="module")
def world(tmp_path_factory: pytest.TempPathFactory) -> W.World:
    return W.build_world(tmp_path_factory.mktemp("r4world"))


@pytest.fixture(scope="module")
def res(world):
    return W.resources(world)


def test_allele_features_snv(world, res):
    alt = W.other(world.ref[W.FP])
    f = real._allele_features(res.bam, res.fasta, W.CHROM, W.FP, allele=alt)
    assert f["depth"] >= 30 and 0.08 < f["vaf"] < 0.2
    assert f["alt_fwd"] + f["alt_rev"] == f["alt_reads"]
    assert f["alt_mean_bq"] == 12 and f["ref_mean_bq"] == 35  # planted low-BQ alt reads
    assert f["alt_mapq_lt20_frac"] == 0 and f["alt_proper_pair_frac"] == 1
    assert f["alt_mate_mapq_median"] == 60  # no MQ tag: taken from the mate itself
    assert f["alt_softclip_frac"] == 0 and f["alt_mismatch_burden"] == 0
    assert f["linked_mismatch_max_frac"] == 0 and f["local_nonref_sites_100bp"] == 0
    assert "stutter_frac" not in f
    cm = real._allele_features(
        res.bam, res.fasta, W.CHROM, W.CMRG_SNV, allele=W.other(world.ref[W.CMRG_SNV])
    )
    assert cm["alt_mapq_lt20_frac"] == 1.0 and 0.4 < cm["vaf"] < 0.6
    # The PSV site sits 500 bp away: not within the 100 bp local window...
    assert cm["local_nonref_sites_100bp"] == 0
    # ...while next to the truth deletion the indel columns are local non-reference sites.
    adj = real._allele_features(
        res.bam, res.fasta, W.CHROM, W.ADJ_NOISE, allele=W.other(world.ref[W.ADJ_NOISE])
    )
    assert adj["local_nonref_sites_100bp"] >= 2


def test_allele_features_length_mode(res):
    run = (W.HP_POS, W.HP_POS + 12)
    f = real._allele_features(
        res.bam, res.fasta, W.CHROM, W.HP_POS, run=run, delta=-1, unit="A", truth_deltas=(0, -1)
    )
    assert f["run_length"] == 12 and 0.3 < f["vaf"] < 0.6
    assert f["stutter_frac"] > 0.05  # the planted -2 reads
    assert f["alt_mean_bq"] == 35  # anchor base, never a deleted one
    assert f["linked_mismatch_max_frac"] == 0  # its own deletion is not a linked mismatch


def test_scan_cache_persists(tmp_path):
    p = tmp_path / "c.jsonl"
    c = real.ScanCache(p)
    calls = []
    assert c.memo("k", lambda: calls.append(1) or {"a": 1}) == {"a": 1}
    assert c.memo("k", lambda: calls.append(1) or {"a": 2}) == {"a": 1} and calls == [1]
    again = real.ScanCache(p)
    assert "k" in again and again.get("k") == {"a": 1}
    assert real.ScanCache(None).get("x") is None


def test_pair_fields_mates_elsewhere_and_unmapped():
    from types import SimpleNamespace

    def src(**kw):
        base = dict(
            flag=0x1 | 0x2 | 0x20 | 0x40 | 0x400,
            is_paired=True,
            mate_is_unmapped=False,
            next_reference_name="chr20",
            next_reference_start=10_500,
            template_length=400,
        )
        return SimpleNamespace(**{**base, **kw})

    def out():
        return real._OutRead(1, 1500, "n", 0, "A", None, None, 60)

    rec = out()
    real._pair_fields(rec, src(), "chr20", 9000, -8000, 1, 3)
    assert rec.flag == 0x1 | 0x2 | 0x20 | 0x40  # duplicate bit dropped
    assert (rec.mate_tid, rec.mate_start, rec.tlen) == (1, 2500, 400)
    rec = out()
    real._pair_fields(rec, src(next_reference_name="chr7"), "chr20", 9000, -8000, 1, 3)
    assert (rec.mate_tid, rec.mate_start, rec.tlen) == (3, 0, 0)  # sentinel, no coordinate
    rec = out()
    real._pair_fields(rec, src(next_reference_start=20_000), "chr20", 9000, -8000, 1, 3)
    assert rec.mate_tid == 3
    rec = out()
    real._pair_fields(rec, src(mate_is_unmapped=True), "chr20", 9000, -8000, 1, 3)
    assert (rec.mate_tid, rec.mate_start) == (1, 1500)
    rec = out()
    real._pair_fields(rec, src(is_paired=False, flag=0x10), "chr20", 9000, -8000, 1, 3)
    assert rec.flag == 0x10 and rec.mate_tid is None


def test_write_real_fixtures_keeps_pairs(world, tmp_path):
    sites = [
        real.RealSite("mosaic_snv_pos", W.CHROM, W.MOSAIC, {"variant_present": "yes"}, False),
        real.RealSite("anchor_het_snv", W.CHROM, W.ANCHOR_HET, {"variant_present": "yes"}, False),
    ]
    paths, items, mapping = real.write_real_fixtures(
        sites,
        tmp_path,
        3,
        bam_path=str(world.paths["bam"]),
        ref_path=str(world.paths["ref"]),
        keep_pairs=True,
    )
    with pysam.AlignmentFile(str(paths["bam"])) as bam:
        assert real.MATE_SENTINEL in bam.references
        assert bam.get_reference_length(real.MATE_SENTINEL) == real.MATE_SENTINEL_LEN
        reads = list(bam.fetch(until_eof=True))
    assert reads and all(r.is_paired for r in reads)
    assert all(r.flag & ~real.PAIR_FLAG_MASK == 0 for r in reads)
    in_slot = [r for r in reads if r.next_reference_name != real.MATE_SENTINEL]
    assert in_slot  # every planted mate lies in its slot (mates elsewhere: see below)
    # A re-indexed mate keeps its distance from the read.
    assert all(
        abs(r.next_reference_start - r.reference_start) == W.INSERT - W.READ for r in in_slot
    )
    with_mq = [r for r in reads if r.has_tag("MQ")]
    assert with_mq and all(r.get_tag("MQ") == 60 for r in with_mq)
    assert all({t for t, _ in r.get_tags()} <= {"MQ"} for r in reads)
    with pysam.FastaFile(str(paths["reference"])) as fa:
        assert set(fa.fetch(real.MATE_SENTINEL)) == {"N"}
    # Default mode is unchanged (single-end, no tags, no sentinel).
    paths2, _, _ = real.write_real_fixtures(
        sites,
        tmp_path / "se",
        3,
        bam_path=str(world.paths["bam"]),
        ref_path=str(world.paths["ref"]),
    )
    with pysam.AlignmentFile(str(paths2["bam"])) as bam:
        assert real.MATE_SENTINEL not in bam.references
        assert all(r.flag in (0, 16) and not r.get_tags() for r in bam.fetch(until_eof=True))


# --- generators ---------------------------------------------------------------------------


def _by(rows, stratum):
    return [r for r in rows if r["stratum"] == stratum]


def test_iter_mosaic_sites(res):
    rows = list(real.iter_mosaic_sites(res, [W.CHROM]))
    assert [r["site0"] for r in rows] == [W.MOSAIC]
    r = rows[0]
    assert r["label"]["zygosity"] == "low_vaf" and r["truth_source"] == "mosaic_v1.1"
    assert r["exclude"] is None and r["alt_reads"] >= 4 and r["pair"] == "mosaic"


def test_iter_mosaic_caller_fps(world, res):
    rows = list(
        real.iter_mosaic_caller_fps(res, [W.CHROM], sources={"c": str(world.paths["caller"])})
    )
    by = {r["site0"]: r for r in rows}
    assert by[W.FP]["exclude"] is None and by[W.FP]["context"]["fp_sources"] == ["c"]
    assert by[W.FP_NEAR_TRUTH]["exclude"] == "v5.0q_within_30bp"
    assert by[W.MOSAIC]["exclude"] == "mosaic_v1.1_within_30bp"


def test_iter_tr_sites(res):
    pos = list(real.iter_tr_sites(res, [W.CHROM], polarity="pos", seed=1))
    assert [r["site0"] for r in pos] == [W.HP_POS]  # first deleted base
    p = pos[0]
    assert p["context"]["run"] == [W.HP_POS, W.HP_POS + 12] and p["context"]["delta"] == -1
    assert strata_r4.eligible(p, "default")[0]
    neg = list(real.iter_tr_sites(res, [W.CHROM], polarity="neg", seed=1))
    by = {tuple(r["context"]["run"]): r for r in neg}
    assert by[(W.HP_POS, W.HP_POS + 12)]["exclude"] == "truth_in_run_30bp"
    n = by[(W.HP_NEG, W.HP_NEG + 12)]
    assert n["exclude"] is None and n["site0"] == W.HP_NEG and n["context"]["delta"] == -1
    assert n["label"] == {} and strata_r4.eligible(n, "default")[0]
    assert strata_r4.offset_in_run(n) == strata_r4.offset_in_run(p) == 0.0


def test_iter_cmrg_and_psv(world, res):
    cm = list(real.iter_cmrg_sites(res, [W.CHROM]))
    assert [r["site0"] for r in cm] == [W.CMRG_SNV]
    assert cm[0]["truth_source"] == "cmrg_v1.00+v5.0q" and cm[0]["exclude"] is None
    assert strata_r4.cmrg_difficult(cm[0])  # all reads MAPQ 15
    psv = list(
        real.iter_psv_candidates(
            res, [W.CHROM], seed=1, dragen_vcf=str(world.paths["dragen"]), scan_windows=10
        )
    )
    by = {r["site0"]: r for r in psv}
    assert by[W.PSV]["exclude"] is None and by[W.PSV]["context"]["psv_source"] == "dragen_fp"
    # A CMRG truth record is never a false positive (the scan sees it; isolation drops it).
    assert by[W.CMRG_SNV]["exclude"] == "cmrg_v1.00_within_50bp"
    assert all(20000 <= s < 23000 for s in by)
    terr = real.psv_territory(res, W.CHROM)
    assert terr == [(20000, 23000)]


def test_iter_indel_adjacent(world, res):
    pos = list(real.iter_indel_adjacent(res, [W.CHROM], polarity="pos", seed=1))
    assert [r["site0"] for r in pos if not r["exclude"]] == [W.ADJ_SNV]
    assert pos[0]["context"]["gap"] == 6 and pos[0]["truth_source"] == "v4.2.1+v5.0q"
    neg = list(real.iter_indel_adjacent(res, [W.CHROM], polarity="neg", seed=1))
    by = {r["site0"]: r for r in neg}
    assert by[W.ADJ_NOISE]["exclude"] is None
    # The truth SNV is a pileup signal too, but a truth record at the POI excludes it here.
    assert by[W.ADJ_SNV]["exclude"] == "truth_within_3bp"


def test_adjacent_gap_and_anchoring():
    assert real._gap(10, 12, 15) == 2 and real._gap(15, 12, 15) == 1 and real._gap(13, 12, 15) == 0
    seq = "ACGT" * 10 + "AAAAAA" + "CGTACGTA" * 5

    def fetch(c, s, e):
        return seq[max(0, s) : e]

    # A POI inside the A run: a truth span 2 bp past the run's end is reachable.
    assert not real.anchoring_safe(fetch, "c", 42, [(47, 48)])
    assert real.anchoring_safe(fetch, "c", 42, [(52, 53)])


def test_iter_anchor_sites(res):
    rows = list(real.iter_anchor_sites(res, [W.CHROM], seed=1, limit=3))
    n = Counter(r["stratum"] for r in rows)
    assert n["anchor_het_snv"] == 1 and n["anchor_hom_snv"] == 1 and n["anchor_clean_neg"] == 3
    for r in _by(rows, "anchor_clean_neg"):
        assert r["alt_reads"] == 0 and r["depth"] >= 20


def _inject(ts, *recs, drop=None):
    """Add truth records to (and optionally drop ``drop(rec)`` ones from) a loaded set."""
    keep = [r for r in ts.records(W.CHROM) if drop is None or not drop(r)]
    ts._index(W.CHROM, keep + list(recs))


def _rec(pos0, source, ref="A", alt="C"):
    from bamcp.eval.visual.truth import TruthRecord

    return TruthRecord(W.CHROM, pos0 + 1, ref, (alt,), (0, 1), "PASS", source)


def test_tr_positive_isolation_ignores_only_its_own_v421_record(world):
    """A different v4.2.1 record in the run +/- 30 bp breaks isolation even where v4.2.1
    does not cover the variant (it used to pass as "the one v4.2.1 record")."""
    from bamcp.eval.visual.truth import BedSet, Intervals

    r = W.resources(world)
    r.v421.bed = BedSet({W.CHROM: Intervals([(23000, 59000)])})  # no v4.2.1 coverage here
    _inject(r.v421, drop=lambda x: x.start0 < W.HP_POS + 20 and x.end0 > W.HP_POS - 20)
    clean = list(real.iter_tr_sites(r, [W.CHROM], polarity="pos", seed=1))
    assert [x["exclude"] for x in clean] == [None]
    _inject(r.v421, _rec(W.HP_POS + 12 + 20, "v4.2.1"))
    got = list(real.iter_tr_sites(r, [W.CHROM], polarity="pos", seed=1))
    assert [x["exclude"] for x in got] == ["truth_in_run_30bp"]


def test_isolation_and_clearance_count_every_truth_set(world):
    """Design 3.0: no other truth record of *any* set near the POI (CMRG / mosaic too)."""
    r = W.resources(world)
    _inject(r.cmrg, _rec(W.HP_NEG + 12 + 10, "cmrg_v1.00"), _rec(W.FP + 5, "cmrg_v1.00"))
    _inject(r.cmrg, _rec(W.ANCHOR_HET + 5, "cmrg_v1.00"))
    _inject(r.mosaic, _rec(W.PSV + 10, "mosaic_v1.1"))
    neg = {
        tuple(x["context"]["run"]): x
        for x in real.iter_tr_sites(r, [W.CHROM], polarity="neg", seed=1)
    }
    assert neg[(W.HP_NEG, W.HP_NEG + 12)]["exclude"] == "truth_in_run_30bp"
    fps = {
        x["site0"]: x
        for x in real.iter_mosaic_caller_fps(
            r, [W.CHROM], sources={"caller": str(world.paths["caller"])}
        )
    }
    assert fps[W.FP]["exclude"] == "cmrg_v1.00_within_30bp"
    psv = {
        x["site0"]: x
        for x in real.iter_psv_candidates(
            r, [W.CHROM], seed=1, dragen_vcf=str(world.paths["dragen"]), scan_windows=0
        )
    }
    assert psv[W.PSV]["exclude"] == "mosaic_v1.1_within_50bp"
    anchors = list(real.iter_anchor_sites(r, [W.CHROM], seed=1, limit=3))
    assert not _by(anchors, "anchor_het_snv")


# --- pools, checks, draw --------------------------------------------------------------------


def test_collect_pools_and_label_checks(world, res):
    rows = real.collect_pools_r4(
        res,
        chroms=[W.CHROM],
        truth_chroms=[W.CHROM],
        seed=1,
        fp_sources={"c": str(world.paths["caller"])},
        dragen_vcf=str(world.paths["dragen"]),
        limits={"psv_windows": 5, "anchor": 2},
    )
    assert all(r["cluster"] and r["partition"] in ("pilot", "main") for r in rows)
    psv = next(r for r in rows if r["stratum"] == "segdup_psv_neg" and r["site0"] == W.PSV)
    assert psv["cluster"] == "sd:chr20:20000"
    seen: list = []
    chk = real.LabelChecks(res, gnomad_fn=lambda *a: seen.append(a) or 0.001)
    get = {(r["stratum"], r["site0"]): r for r in rows}
    assert chk(get[("mosaic_snv_pos", W.MOSAIC)]) is None
    assert get[("mosaic_snv_pos", W.MOSAIC)]["checks"]["hifi"]["pass"]
    fp = get[("mosaic_caller_fp_neg", W.FP)]
    assert chk(fp) is None and fp["checks"]["gnomad_af"] == 0.001 and len(seen) == 1
    assert fp["checks"]["trio"]["recurrent"] is False  # recorded, not required
    assert chk(psv) is None and psv["checks"]["trio"]["recurrent"]
    stut = get[("tr_stutter_neg", W.HP_NEG)]
    assert chk(stut) is None
    # Common in gnomAD: excluded.
    keep = res.cache
    try:
        res.cache = real.ScanCache(None)
        common = real.LabelChecks(res, gnomad_fn=lambda *a: 0.2)
        assert common({**fp, "checks": {}}) == "gnomad_af_ge_1pct"
        assert common.counts == {"gnomad_af_ge_1pct": 1}

        def boom(*a):
            raise OSError("down")

        res.cache = real.ScanCache(None)
        down = real.LabelChecks(res, gnomad_fn=boom)
        assert down({**fp, "checks": {}}) == "gnomad_unavailable"  # never read as "rare"
    finally:
        res.cache = keep
    rep = real.pools_report_r4(rows)
    assert rep["mosaic_snv_pos"]["candidates"] == 1


def test_label_checks_without_handles(world):
    r0 = W.resources(world, hifi=False, parents=False)
    chk = real.LabelChecks(r0, gnomad=False)
    row = {
        "stratum": "segdup_psv_neg",
        "polarity": "neg",
        "chrom": W.CHROM,
        "site0": W.PSV,
        "context": {"candidate_allele": "A"},
        "label": {},
    }
    assert chk(row) is None
    assert row["checks"] == {"hifi": "not_run", "trio": "not_run"}


def test_manifest_round_trip(tmp_path):
    rows = [{"chrom": "chr1", "site0": 5, "stratum": "x", "checks": {"s": {1, 2}}}]
    digest = real.manifest_to_jsonl(rows, tmp_path / "m.jsonl")
    assert len(digest) == 64
    assert real.manifest_from_jsonl(tmp_path / "m.jsonl")[0]["checks"]["s"] == [1, 2]


def test_build_r4_run_set_writes_r4_labels(world, tmp_path, monkeypatch):
    from types import SimpleNamespace

    from bamcp.eval.visual import stimuli

    rows = [
        {
            **real._new_row(
                "mosaic_snv_pos",
                W.CHROM,
                W.MOSAIC,
                truth_source="mosaic_v1.1",
                label={
                    "variant_present": "yes",
                    "artifact_verdict": "real",
                    "variant_type": "SNV",
                    "alt_allele": "G",
                    "zygosity": "low_vaf",
                    "artifact_type": "none",
                    "phase": "na",
                },
            ),
            "cluster": "bin:chr20:0",
            "partition": "pilot",
            "knob": "default",
        },
        {
            **real._new_row(
                "tr_stutter_neg",
                W.CHROM,
                W.HP_NEG,
                truth_source="v5.0q",
                context={"run": [W.HP_NEG, W.HP_NEG + 12]},
            ),
            "label": strata_r4.label_for("tr_stutter_neg", {}),
            "cluster": "bin:chr20:0",
            "partition": "pilot",
            "knob": "alt",
        },
    ]
    man = tmp_path / "m.jsonl"
    real.manifest_to_jsonl(rows, man)
    seen: dict = {}

    def fake(out, **kw):
        seen.update(kw)
        (Path(out) / "labels").mkdir(parents=True, exist_ok=True)
        return SimpleNamespace(n_items=len(kw["items"]))

    monkeypatch.setattr(stimuli, "build_items_run_set", fake)
    res = real.build_r4_run_set(
        tmp_path / "run",
        man,
        seed=3,
        models=["m"],
        bam_path=str(world.paths["bam"]),
        ref_path=str(world.paths["ref"]),
        log=lambda *_: None,
    )
    assert res.n_items == 2 and seen["r4"] is True and seen["sample_kind"] == "sequenced"
    assert seen["sample_line"] == real.SAMPLE_LINE_R4 and "mosaic" in real.SAMPLE_LINE_R4
    assert list(seen["conditions"]) == list(real.R4_CONDITIONS)
    labs = {lab["case_type"]: lab for lab in seen["label_records"]}
    assert labs["tr_stutter_neg"]["artifact_type"] == "homopolymer_indel"
    for lab in labs.values():
        assert lab["stratum"] == lab["case_type"]
        assert set(real.R4_LABEL_KEYS) <= set(lab)
    assert labs["tr_stutter_neg"]["knob"] == "alt" and labs["mosaic_snv_pos"]["pair"] == "mosaic"
    with pysam.AlignmentFile(seen["bam"]) as bam:
        assert real.MATE_SENTINEL in bam.references  # pairs kept
    mapping = (tmp_path / "run" / "labels" / "mapping.jsonl").read_text().splitlines()
    assert len(mapping) == 2 and '"stratum"' in mapping[0]


def test_r4_cli_dispatch(monkeypatch, tmp_path):
    from types import SimpleNamespace

    from bamcp.eval.visual import cli

    calls: dict = {}
    rows = [{"stratum": "mosaic_snv_pos", "exclude": None}]
    monkeypatch.setattr(cli, "_r4_rows", lambda args, checks: (rows, SimpleNamespace()))

    def fake_pools(r, overrides=None):
        calls["pools_overrides"] = overrides
        return {"mosaic_snv_pos": {"candidates": 1, "eligible": {}}}

    monkeypatch.setattr(real, "pools_report_r4", fake_pools)
    monkeypatch.setattr(real, "LabelChecks", lambda res, gnomad: SimpleNamespace(counts={}))

    def fake_select(r, **kw):
        calls["select"] = kw
        return [{"chrom": "chr1", "site0": 1}], {"short": {"x": 1}}

    monkeypatch.setattr(real, "select_r4", fake_select)
    out = tmp_path / "pools.json"
    assert cli.main(["pools-r4", "--out", str(out)]) == 0 and out.exists()
    man = tmp_path / "m.jsonl"
    argv = [
        "select-r4",
        "--split",
        "pilot",
        "--seed",
        "5",
        "--out",
        str(man),
        "--knob",
        "mosaic_snv_pos=alt",
        "--contamination",
        "0.02",
        "--no-gnomad",
    ]
    assert cli.main(argv) == 0
    assert calls["select"]["split"] == "pilot" and calls["select"]["seed"] == 5
    assert calls["select"]["knob_overrides"] == {"mosaic_snv_pos": "alt"}
    assert calls["select"]["param_overrides"]["mosaic_caller_fp_neg"]["vaf_floor_override"] == 0.08
    assert man.exists() and (tmp_path / "m.jsonl.report.json").exists()
    assert cli.main([*argv, "--strict"]) == 1  # a short stratum fails a strict draw

    seen: dict = {}

    def fake_build(out, manifest, **kw):
        seen.update(kw)
        return SimpleNamespace(n_stimuli=0, n_items=0, n_runs=0, out=tmp_path)

    monkeypatch.setattr(real, "build_r4_run_set", fake_build)
    monkeypatch.setattr(cli, "_lint", lambda *_: True)
    assert cli.main(["build-r4", "--out", str(tmp_path), "--manifest", "m", "--seed", "1"]) == 0
    assert seen["conditions"] == list(real.R4_CONDITIONS)


def test_r4_resources_retry_reopens_handles(res):
    """A dropped remote read is retried after re-opening every handle."""
    import time as _time

    opened: list[str] = []
    keep = dict(res.openers), res.bam
    res.openers = {"bam": lambda: opened.append("bam") or keep[1]}
    calls = {"n": 0}

    def flaky():
        calls["n"] += 1
        if calls["n"] < 3:
            raise OSError("truncated file")
        return "ok"

    sleeps: list = []
    orig = _time.sleep
    _time.sleep = sleeps.append  # type: ignore[assignment]
    try:
        assert res.retry(flaky) == "ok" and opened == ["bam", "bam"]
        with pytest.raises(OSError):
            res.retry(lambda: (_ for _ in ()).throw(OSError("down")), tries=2)
    finally:
        _time.sleep = orig  # type: ignore[assignment]
        res.openers = keep[0]


def test_r4_resources_fetch_blocks_match_fasta(world, res):
    for s, e in ((0, 10), (16_380, 16_390), (5, 40_000), (59_990, 60_000), (100, 100)):
        assert res.fetch(W.CHROM, s, e) == world.ref[s:e]
    assert len(res._blocks) <= real.FETCH_BLOCKS_MAX


def test_sv_spans_cover_the_duplicated_copy_of_an_insertion(tmp_path):
    """A 500 bp insertion (tandem dup) makes the next 500 reference bases variable-copy;
    the anchor +/- 50 bp alone missed them."""
    vcf = W._vcf(
        tmp_path / "sv.vcf.gz",
        [
            (10001, "A", "A" + "C" * 500, "PASS", "0/1"),
            (20001, "G" + "T" * 100, "G", "PASS", "0/1"),
        ],
    )
    bs = real._sv_spans(vcf, [W.CHROM])
    assert bs.contains(W.CHROM, 10000 + 400) and bs.contains(W.CHROM, 10000 - 400)
    assert not bs.contains(W.CHROM, 10000 + 600)
    assert bs.contains(W.CHROM, 20000 + 100) and not bs.contains(W.CHROM, 20000 + 200)


def test_anchor_spans_reach_the_whole_poi_run(world):
    """A truth record 100 bp from the indel can still sit inside a long run holding the
    POI; the anchoring check used to see only records within 40 bp of the indel."""
    r = W.resources(world)
    far = _rec(W.ADJ_DEL + 100, "v5.0q")
    _inject(r.v5q, far)
    _inject(r.cmrg, _rec(W.ADJ_DEL - 150, "cmrg_v1.00"))
    spans = real._anchor_spans(r, W.CHROM, W.ADJ_DEL, W.ADJ_DEL + 4)
    assert (far.start0, far.end0) in spans
    assert (W.ADJ_DEL - 150, W.ADJ_DEL - 149) in spans
    # And anchoring_safe then rejects a POI whose run reaches that record.
    seq = "C" * 5 + "A" * 120 + "G" * 5
    fetch = lambda chrom, a, b: seq[max(0, a) : b]  # noqa: E731
    assert real.anchoring_safe(fetch, W.CHROM, 10, [(110, 111)]) is False
