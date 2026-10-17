"""Round-4 truth handling: FILTER-aware truth sets, normalisation, concordance, HiFi and
trio checks, clusters and the pilot/main partition. Hermetic (tiny local files only)."""

from __future__ import annotations

import gzip
from collections import Counter
from pathlib import Path

import pysam
import pytest

from bamcp.eval.visual import truth as T

from . import r4_world as W

pytestmark = pytest.mark.unit


@pytest.fixture(scope="module")
def world(tmp_path_factory: pytest.TempPathFactory) -> W.World:
    return W.build_world(tmp_path_factory.mktemp("r4world"))


def _fetch(ref: str):
    return lambda chrom, s, e: ref[max(0, s) : e]


# --- records ---------------------------------------------------------------------------


@pytest.mark.parametrize(
    ("ref", "alt", "kind"),
    [
        ("A", "G", "SNV"),
        ("AC", "A", "del"),
        ("A", "ACG", "ins"),
        ("AC", "GT", "mnp"),
        ("ACG", "TA", "complex"),
    ],
)
def test_variant_kind(ref, alt, kind):
    assert T.variant_kind(ref, alt) == kind


def test_record_spans_and_zygosity():
    r = T.TruthRecord("chr1", 101, "ACGT", ("A",), (0, 1), "PASS", "x")
    assert (r.start0, r.end0) == (100, 104)  # full reference span
    assert r.zygosity == "het" and r.is_pass and r.kind == "del"
    assert r.overlaps(104, 105, bp=1) and not r.overlaps(105, 106, bp=1)
    assert T.TruthRecord("c", 1, "A", ("G",), (1, 1), "HET1", "x").zygosity == "hom_alt"
    assert T.TruthRecord("c", 1, "A", ("G", "T"), (1, 2), "PASS", "x").kind == "multi"
    assert T.TruthRecord("c", 1, "A", ("G",), (None, 1), ".", "x").zygosity == "unknown"


# --- truth sets --------------------------------------------------------------------------


def test_truthset_filter_aware(world):
    p = world.paths
    v5 = T.TruthSet.load("v5.0q", p["v5q"], p["v5q_bed"], chroms=[W.CHROM])
    pos = {r.start0 for r in v5.positives(W.CHROM)}
    assert 40000 not in pos  # HET1 is not a positive...
    assert v5.records_near(W.CHROM, 39990, 39991, 10)  # ...but blocks isolation
    assert W.MOSAIC not in pos and W.ANCHOR_HET in pos
    assert v5.covers(W.CHROM, 1000, 59000) and not v5.covers(W.CHROM, 500, 1500)
    assert v5.asserts_reference(W.CHROM, 50000, 50001, 30)
    assert not v5.asserts_reference(W.CHROM, 40000, 40001, 0)


def test_truthset_region_queries_match_loaded(world):
    """A chromosome not preloaded is answered by tabix region queries, identically."""
    p = world.paths
    full = T.TruthSet.load("v5", p["v5q"], p["v5q_bed"], chroms=[W.CHROM])
    lazy = T.TruthSet.load("v5", p["v5q"], p["v5q_bed"], chroms=[])
    assert W.CHROM not in lazy.by_chrom
    for site in (W.FP_NEAR_TRUTH, 40000, W.ADJ_DEL + 2, 50000):
        a = [(r.pos, r.ref) for r in full.records_near(W.CHROM, site, site + 1, 30)]
        b = [(r.pos, r.ref) for r in lazy.records_near(W.CHROM, site, site + 1, 30)]
        assert a == b
    assert lazy.records_near("chrZ", 1, 2, 30) == []
    # positives() loads the chromosome on demand.
    assert {r.pos for r in lazy.positives(W.CHROM)} == {r.pos for r in full.positives(W.CHROM)}


def test_truthset_unindexed_reads_whole_file(world):
    cm = T.TruthSet.load("cmrg", world.paths["cmrg"], world.paths["cmrg_bed"], chroms=[W.CHROM])
    assert [r.start0 for r in cm.positives(W.CHROM)] == [W.CMRG_SNV]
    assert cm.find(W.CHROM, W.CMRG_SNV + 1, world.ref[W.CMRG_SNV], W.other(world.ref[W.CMRG_SNV]))


def test_records_near_uses_full_reference_span(world):
    v4 = T.TruthSet.load("v4", world.paths["v421"], world.paths["v421_bed"], chroms=[W.CHROM])
    # The 3 bp deletion anchored at ADJ_DEL spans ADJ_DEL..ADJ_DEL+3: a site 2 bp past its
    # last deleted base is within 2 bp, though 5 bp from the anchor.
    end = W.ADJ_DEL + 4

    def dels(s0: int) -> list:
        return [r for r in v4.records_near(W.CHROM, s0, s0 + 1, 2) if r.kind == "del"]

    assert dels(end + 1) and not dels(end + 2)


def test_bedset_and_load_bed_multi(tmp_path):
    bed = tmp_path / "b.bed.gz"
    with gzip.open(bed, "wt") as fh:
        fh.write("chr1\t10\t20\nchr1\t15\t30\nchr2\t5\t6\ntrack x\n")
    bs = T.BedSet.load(bed)
    assert bs.covers("chr1", 10, 30) and bs.contains("chr2", 5) and not bs.contains("chr3", 5)
    assert bs.total_bp() == 21
    assert set(T.load_bed_multi(bed, ["chr2"])) == {"chr2"}


# --- normalisation + concordance ------------------------------------------------------------


def test_left_normalize_shifts_indels_left_in_repeats():
    ref = "GCAAAAT"  # 1-based: G1 C2 A3 A4 A5 A6 T7
    f = _fetch(ref)
    # Deleting the last A, written right-shifted (POS 6: "AT" -> "T" is not valid VCF, so
    # use "AA" -> "A" at POS 5): normalises to the anchor before the run.
    assert T.left_normalize(f, "c", 5, "AA", "A") == (2, "CA", "C")
    # Insertion of an A written at the run end.
    assert T.left_normalize(f, "c", 6, "A", "AA") == (2, "C", "CA")
    # Shared leading bases are trimmed for an MNP-style record.
    assert T.left_normalize(f, "c", 2, "CAA", "CAT") == (4, "A", "T")
    assert T.left_normalize(f, "c", 3, "A", "A") == (3, "A", "A")


def test_concordant_is_representation_and_zygosity_aware(tmp_path):
    ref = "G" + "CAAAAT" + "G" * 200
    fa = tmp_path / "r.fa"
    fa.write_text(f">c\n{ref}\n")
    f = _fetch(ref)

    def ts(recs, bed=((0, 200),)):
        bs = T.BedSet({"c": T.Intervals(list(bed))})
        return T.TruthSet("o", recs, bs)

    # 1-based: G1 C2 A3 A4 A5 A6 T7
    a = T.TruthRecord("c", 2, "CA", ("C",), (0, 1), "PASS", "a")  # left-aligned deletion
    b = T.TruthRecord("c", 5, "AA", ("A",), (0, 1), "PASS", "b")  # same, right-shifted
    hom = T.TruthRecord("c", 5, "AA", ("A",), (1, 1), "PASS", "b")
    assert T.concordant(a, ts([b]), f) is True
    assert T.concordant(a, ts([hom]), f) is False  # zygosity differs
    assert T.concordant(a, ts([]), f) is False  # the other set asserts reference
    assert T.concordant(a, ts([b], bed=((100, 200),)), f) is None  # not covered: no claim
    assert T.reference_concordant(ts([]), "c", 50, 51, 5) is True
    assert T.reference_concordant(ts([b]), "c", 3, 4, 5) is False
    assert T.reference_concordant(ts([], bed=((100, 200),)), "c", 4, 5, 5) is None


@pytest.mark.parametrize(
    ("seq", "i", "expect"),
    [
        ("GCAAAAAT", 3, (2, 7, "A")),
        ("GACACACAT", 4, (1, 8, "AC")),
        ("TTAGCAGCAGCAGG", 6, (2, 13, "AGC")),
        ("ACGT", 1, (1, 2, "C")),
    ],
)
def test_repeat_run(seq, i, expect):
    assert T.repeat_run(seq, i) == expect


def test_is_tr_run():
    assert T.is_tr_run(10, "A", 10) and not T.is_tr_run(9, "A", 10)
    assert T.is_tr_run(20, "AC", 10) and not T.is_tr_run(19, "AC", 10)
    assert not T.is_tr_run(30, "ACGT", 10)


# --- read-level ---------------------------------------------------------------------------


def test_spanning_hist_counts_run_length_alleles(world):
    with pysam.AlignmentFile(str(world.paths["bam"])) as bam:
        h = T.spanning_hist(bam, W.CHROM, W.HP_POS, W.HP_POS + 12)
        hn = T.spanning_hist(bam, W.CHROM, W.HP_NEG, W.HP_NEG + 12)
    n = sum(h.values())
    assert n > 20 and set(h) <= {0, -1, -2}
    assert 0.35 < h[-1] / n < 0.65 and h[-2] > 0
    assert 0.10 < T.unit_fraction(hn, 1) < 0.25
    with pysam.AlignmentFile(str(world.paths["bam"])) as bam:
        c = T.base_counts(bam, W.CHROM, W.MOSAIC)
    alt = W.other(world.ref[W.MOSAIC])
    assert 0.05 < c[alt] / sum(c.values()) < 0.25


def test_read_run_delta_requires_spanning():
    class R:
        reference_start, reference_end = 100, 200
        cigartuples = [(0, 50), (2, 2), (0, 48)]

    assert T.read_run_delta(R(), 148, 160) == -2
    assert T.read_run_delta(R(), 98, 110) is None  # does not reach SPAN_FLANK before the run

    class Ins(R):
        cigartuples = [(0, 50), (1, 1), (0, 50)]

    assert T.read_run_delta(Ins(), 150, 155) == 1
    assert T.read_run_delta(Ins(), 160, 170) == 0


# --- orthogonal (HiFi) check ---------------------------------------------------------------


def _row(stratum, site0, *, polarity, label=None, context=None, features=None):
    return {
        "stratum": stratum,
        "chrom": W.CHROM,
        "site0": site0,
        "polarity": polarity,
        "label": label or {},
        "context": context or {},
        "features": features or {},
    }


def test_orthogonal_check_snv(world):
    alt = W.other(world.ref[W.MOSAIC])
    with pysam.AlignmentFile(str(world.paths["hifi"])) as hifi:
        pos = T.orthogonal_check(
            _row(
                "mosaic_snv_pos",
                W.MOSAIC,
                polarity="pos",
                label={"zygosity": "low_vaf", "alt_allele": alt},
            ),
            hifi,
        )
        assert pos["pass"] and pos["alt_reads"] >= 3
        # The same allele read as a het needs >= 25%: fails.
        het = T.orthogonal_check(
            _row("x", W.MOSAIC, polarity="pos", label={"zygosity": "het", "alt_allele": alt}),
            hifi,
        )
        assert not het["pass"] and het["reason"] == "hifi_alt_absent"
        fp_alt = W.other(world.ref[W.FP])
        neg = T.orthogonal_check(
            _row(
                "mosaic_caller_fp_neg", W.FP, polarity="neg", context={"candidate_allele": fp_alt}
            ),
            hifi,
        )
        assert neg["pass"] and neg["alt_reads"] == 0
        # A "negative" at a real het fails: HiFi carries the allele.
        bad = T.orthogonal_check(
            _row(
                "segdup_psv_neg",
                W.ANCHOR_HET,
                polarity="neg",
                context={"candidate_allele": W.other(world.ref[W.ANCHOR_HET])},
            ),
            hifi,
        )
        assert not bad["pass"] and bad["reason"] == "hifi_carries_allele"


def test_orthogonal_check_repeat_lengths(world):
    run = [W.HP_POS, W.HP_POS + 12]
    with pysam.AlignmentFile(str(world.paths["hifi"])) as hifi:
        het = T.orthogonal_check(
            _row(
                "tr_indel_pos",
                W.HP_POS,
                polarity="pos",
                label={"zygosity": "het"},
                context={"run": run, "unit": "A", "delta": -1},
            ),
            hifi,
        )
        assert het["pass"], het
        hom = T.orthogonal_check(
            _row(
                "tr_indel_pos",
                W.HP_POS,
                polarity="pos",
                label={"zygosity": "hom_alt"},
                context={"run": run, "unit": "A", "delta": -1},
            ),
            hifi,
        )
        assert not hom["pass"]  # half the HiFi reads are reference length
        nrun = [W.HP_NEG, W.HP_NEG + 12]
        neg = T.orthogonal_check(
            _row(
                "tr_stutter_neg",
                W.HP_NEG,
                polarity="neg",
                context={"run": nrun, "unit": "T", "delta": -1},
                features={"stutter_frac": 0.16},
            ),
            hifi,
        )
        assert neg["pass"] and neg["unit_frac"] == 0.0
        # Stutter-free Illumina (0%) cannot pass the "<= half the Illumina fraction" rule
        # when HiFi also shows the allele: use the Illumina BAM as the "HiFi".
    with pysam.AlignmentFile(str(world.paths["bam"])) as ill:
        dirty = T.orthogonal_check(
            _row(
                "tr_stutter_neg",
                W.HP_NEG,
                polarity="neg",
                context={"run": nrun, "unit": "T", "delta": -1},
                features={"stutter_frac": 0.16},
            ),
            ill,
        )
    assert not dirty["pass"]


@pytest.mark.parametrize(
    ("hist", "ok"),
    [
        ({0: 59, -1: 1}, True),  # 1 read, 1.7%
        ({0: 29, -1: 1}, False),  # 1 read but 3.3% > 2%
        ({0: 98, -1: 2}, False),  # 2% but 2 reads > 1
    ],
)
def test_orthogonal_check_stutter_negative_needs_both_read_and_fraction_caps(monkeypatch, hist, ok):
    """Design 3.0: a negative's candidate allele on <= 1 HiFi read *and* <= 2%."""
    monkeypatch.setattr(T, "spanning_hist", lambda *a, **k: Counter(hist))
    row = _row(
        "tr_stutter_neg",
        100,
        polarity="neg",
        context={"run": [100, 112], "unit": "A", "delta": -1},
        features={"stutter_frac": 0.5},  # the ratio rule passes in every case
    )
    assert T.orthogonal_check(row, None)["pass"] is ok


def test_orthogonal_check_low_depth(tmp_path, world):
    header = {"HD": {"VN": "1.6"}, "SQ": [{"SN": W.CHROM, "LN": W.LEN}]}
    empty = tmp_path / "e.bam"
    with pysam.AlignmentFile(str(empty), "wb", header=header):
        pass
    pysam.index(str(empty))
    with pysam.AlignmentFile(str(empty)) as e:
        got = T.orthogonal_check(
            _row("x", 100, polarity="neg", context={"candidate_allele": "A"}), e
        )
    assert not got["pass"] and got["reason"] == "hifi_low_depth"


# --- trio --------------------------------------------------------------------------------


def test_trio_recurrence(world):
    psv_alt = W.other(world.ref[W.PSV])
    with (
        pysam.AlignmentFile(str(world.paths["HG003"])) as f,
        pysam.AlignmentFile(str(world.paths["HG004"])) as m,
    ):
        parents = {"HG003": f, "HG004": m}
        got = T.trio_recurrence(
            _row("segdup_psv_neg", W.PSV, polarity="neg", context={"candidate_allele": psv_alt}),
            parents,
        )
        assert got["recurrent"] and got["parent_truth"] == "unchecked"
        assert got["parents"]["HG003"]["frac"] >= 0.05 and got["parents"]["HG004"]["frac"] == 0
        tr = T.trio_recurrence(
            _row(
                "tr_stutter_neg",
                W.HP_NEG,
                polarity="neg",
                context={"run": [W.HP_NEG, W.HP_NEG + 12], "unit": "T", "delta": -1},
            ),
            parents,
        )
        assert tr["recurrent"]
        # A parent whose own truth carries a variant there is not applicable.
        rec = T.TruthRecord(W.CHROM, W.PSV + 1, "A", ("C",), (0, 1), "PASS", "HG003")
        pt = {"HG003": T.TruthSet("HG003", [rec]), "HG004": T.TruthSet("HG004", [])}
        got = T.trio_recurrence(
            _row("segdup_psv_neg", W.PSV, polarity="neg", context={"candidate_allele": psv_alt}),
            parents,
            pt,
        )
        assert not got["recurrent"] and got["parents"]["HG003"]["applicable"] is False
        assert got["parent_truth"] == "checked"


# --- clusters + partition --------------------------------------------------------------------


def test_cluster_index_components(tmp_path):
    sd = tmp_path / "sd.txt.gz"
    rows = [
        ("chr1", 100, 200, "chr2", 500, 600),
        ("chr2", 550, 700, "chr3", 10, 50),  # overlaps chr2:500-600 -> same component
        ("chr5", 1000, 2000, "chr5", 9000, 9500),
    ]
    with gzip.open(sd, "wt") as fh:
        for c1, s1, e1, c2, s2, e2 in rows:
            fh.write(f"585\t{c1}\t{s1}\t{e1}\tx\t0\t+\t{c2}\t{s2}\t{e2}\n")
    idx = T.cluster_index(sd)
    a = idx.cluster_id("chr1", 150)
    assert a == idx.cluster_id("chr2", 650) == idx.cluster_id("chr3", 20) == "sd:chr1:100"
    assert idx.cluster_id("chr5", 9100) == idx.cluster_id("chr5", 1500) != a
    assert idx.cluster_id("chr1", 2_500_000) == "bin:chr1:2"
    assert idx.sizes[a] == 3


def test_cluster_index_min_identity(tmp_path):
    sd = tmp_path / "sd.txt.gz"
    pad = "\t".join(["0"] * 16)
    with gzip.open(sd, "wt") as fh:
        for c1, s1, c2, s2, ident in (
            ("chr1", 100, "chr2", 500, 0.995),
            ("chr2", 550, "chr3", 10, 0.96),
        ):
            fh.write(
                f"585\t{c1}\t{s1}\t{s1 + 100}\tx\t0\t+\t{c2}\t{s2}\t{s2 + 100}\t{pad}\t{ident}\n"
            )
    loose, strict = T.cluster_index(sd), T.cluster_index(sd, min_identity=0.98)
    assert loose.cluster_id("chr1", 150) == loose.cluster_id("chr3", 20)
    assert strict.cluster_id("chr1", 150) == strict.cluster_id("chr2", 550)
    assert strict.cluster_id("chr3", 20) == "bin:chr3:0"  # its only pairing is below 0.98


def test_partition_is_seeded_by_cluster_and_about_a_quarter():
    assert T.partition("sd:chr1:100") == T.partition("sd:chr1:100")
    parts = Counter(T.partition(f"bin:chr{c}:{i}") for c in range(1, 23) for i in range(60))
    frac = parts["pilot"] / sum(parts.values())
    assert 0.20 < frac < 0.30
    assert Counter(T.partition(f"bin:chr1:{i}", salt="other") for i in range(400)) != Counter(
        T.partition(f"bin:chr1:{i}") for i in range(400)
    )


# --- cache + contamination -------------------------------------------------------------------


def test_cache_dir_and_local_passthrough(tmp_path, monkeypatch):
    monkeypatch.setenv("BAMCP_EVAL_CACHE", str(tmp_path / "c"))
    assert T.cache_dir() == tmp_path / "c" and (tmp_path / "c").is_dir()
    assert T.cached(tmp_path / "x.bed") == tmp_path / "x.bed"
    assert T.is_remote("https://h/x") and not T.is_remote("/x")
    monkeypatch.delenv("BAMCP_EVAL_CACHE")
    monkeypatch.setenv("XDG_CACHE_HOME", str(tmp_path / "xdg"))
    assert T.cache_dir() == tmp_path / "xdg" / "bamcp-eval"


def test_cached_downloads_once(tmp_path, monkeypatch):
    monkeypatch.setenv("BAMCP_EVAL_CACHE", str(tmp_path / "c"))
    calls: list[str] = []

    class Resp:
        def __enter__(self):
            return self

        def __exit__(self, *a):
            return None

        def read(self, n=-1):
            if calls[-1] == "done":
                return b""
            calls.append("done")
            return b"payload"

    def fake_open(url, timeout=0):
        calls.append(url)
        return Resp()

    monkeypatch.setattr(T.urllib.request, "urlopen", fake_open)
    p1 = T.cached("https://example.org/a/file.txt")
    p2 = T.cached("https://example.org/a/file.txt")
    assert p1 == p2 and p1.read_bytes() == b"payload"
    assert calls.count("https://example.org/a/file.txt") == 1


def test_estimate_contamination_clean_sample(world):
    v4 = T.TruthSet.load("v4", world.paths["v421"], world.paths["v421_bed"], chroms=[W.CHROM])
    with pysam.AlignmentFile(str(world.paths["bam"])) as bam:
        got = T.estimate_contamination(bam, v4, [W.CHROM], n_sites=10)
    assert got["sites"] == 1 and got["contamination_est"] == 0.0 and got["ref_fraction"] == 0.0


def test_open_bam_local_and_vcf(world):
    with T.open_bam(str(world.paths["bam"])) as b:
        assert b.references == (W.CHROM,)
    with T.open_vcf(world.paths["v5q"]) as v:
        assert W.CHROM in v.header.contigs


def test_paths_are_not_in_the_repo(tmp_path, monkeypatch):
    """Downloads go to the user cache dir, never into the source tree (the default
    location is checked under a temporary HOME, so the test never writes to the real one)."""
    monkeypatch.delenv("BAMCP_EVAL_CACHE", raising=False)
    monkeypatch.delenv("XDG_CACHE_HOME", raising=False)
    monkeypatch.setenv("HOME", str(tmp_path / "home"))
    repo = Path(__file__).resolve().parents[4]
    got = T.cache_dir().resolve()
    assert repo not in got.parents and got == (tmp_path / "home" / ".cache" / "bamcp-eval")


def test_cached_rejects_a_short_download(tmp_path, monkeypatch):
    """A proxy that drops the connection must not leave a silently truncated cache file."""
    import urllib.error

    monkeypatch.setenv("BAMCP_EVAL_CACHE", str(tmp_path / "c"))
    calls: list[str] = []

    class Short:
        headers = {"Content-Length": "100"}

        def __enter__(self):
            return self

        def __exit__(self, *a):
            return None

        def __init__(self):
            self.sent = False

        def read(self, n=-1):
            if self.sent:
                return b""
            self.sent = True
            return b"x" * 10

    def short(url, timeout=0):
        calls.append(url)
        return Short()

    monkeypatch.setattr(T.urllib.request, "urlopen", short)
    with pytest.raises(OSError, match="short download"):
        T.cached("https://example.org/big.vcf.gz")
    assert len(calls) == T.DOWNLOAD_RETRIES
    assert not list((tmp_path / "c").iterdir())  # no partial file left behind

    def missing(url, timeout=0):
        calls.append(url)
        raise urllib.error.HTTPError(url, 404, "nf", {}, None)  # type: ignore[arg-type]

    calls.clear()
    monkeypatch.setattr(T.urllib.request, "urlopen", missing)
    with pytest.raises(OSError):
        T.cached("https://example.org/none.tbi")
    assert len(calls) == 1  # a 404 is not retried
