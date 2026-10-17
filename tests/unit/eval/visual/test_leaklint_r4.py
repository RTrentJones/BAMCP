"""Round-4 leak-lint: section-10 vocabulary, the sentinel mate contig, r4 strata names, and
the selection-side KS / anchoring lints."""

from __future__ import annotations

import itertools
import json
import random
from pathlib import Path

import pysam
import pytest

from bamcp.eval.visual import leaklint, real, stimuli
from bamcp.eval.visual.leaklint import (
    CAPTION_RE,
    R4_CASE_TYPES,
    anchor_window,
    ks_two_sample,
    lint_anchor_collisions,
    lint_poi_offsets,
    lint_run_set,
    lint_sentinel_contig,
    poi_offset,
    r4_case_types,
    r4_vocab_problems,
)

from .test_stimuli_r4 import build_r4, write_pair_fixture

pytestmark = pytest.mark.unit


@pytest.fixture(scope="module")
def pairs(tmp_path_factory: pytest.TempPathFactory) -> dict[str, Path]:
    return write_pair_fixture(tmp_path_factory.mktemp("pairs_lint"))


# ---------------------------------------------------------------------------
# Vocabulary
# ---------------------------------------------------------------------------


@pytest.mark.parametrize("name", ["traw.md", "tcur.md", "tfull.md", "tjson.json"])
@pytest.mark.parametrize(
    "word", ["GIAB", "CMRG", "v4.2.1", "v5.0q", "benchmark", "segdups", "HiFi", "PacBio", "trio"]
)
def test_truth_set_names_are_flagged(name, word):
    problems = r4_vocab_problems(name, f"see the {word} call here\n")
    assert problems and "truth-set vocabulary" in problems[0]


def test_trio_needs_a_word_boundary():
    assert r4_vocab_problems("traw.md", "patrion triode\n") == []


@pytest.mark.parametrize("name", ["traw.md", "tcur.md", "tfull.md"])
def test_mosaic_only_in_caption_and_prompt(name):
    assert r4_vocab_problems(name, "a Mosaic site\n")
    assert r4_vocab_problems("caption.md", stimuli.SAMPLE_LINE_R4) == []
    assert r4_vocab_problems("prompt.md", "cells (mosaic)") == []


def test_product_phrases_are_allowed_in_tool_text_only():
    text = "VAF 8%: below the germline het range (subclonal/mosaic, or an artifact); mosaicism\n"
    assert r4_vocab_problems("tcur.md", text) == []
    assert r4_vocab_problems("tfull.md", text) == []
    assert r4_vocab_problems("traw.md", text)
    # the allow-list strips the phrase, not any other mention
    assert r4_vocab_problems("tcur.md", text + "a mosaic call\n")


def test_sentinel_coordinate_in_text_is_flagged():
    assert r4_vocab_problems("traw.md", "mate at chrU:17\n")
    assert r4_vocab_problems("traw.md", "mate aligned to another contig: 2\n") == []


def test_caption_template_accepts_the_r4_sample_line():
    cap = (
        "# Site\n\n- Region: chr1:1901-2100 (1-based, inclusive)\n"
        f"- Position(s) of interest: chr1:2001\n- Sample: {stimuli.SAMPLE_LINE_R4}\n"
    )
    assert CAPTION_RE.match(cap)
    assert not CAPTION_RE.match(cap.replace("(mosaic)", "(mosaic, low VAF)"))


# ---------------------------------------------------------------------------
# Strata names
# ---------------------------------------------------------------------------


def test_r4_case_types_match_strata_r4_when_present():
    strata = pytest.importorskip("bamcp.eval.visual.strata_r4")
    assert tuple(strata.STRATA_R4) == R4_CASE_TYPES


def test_r4_case_types_include_the_local_tuple():
    assert set(R4_CASE_TYPES) <= set(r4_case_types())
    assert len(R4_CASE_TYPES) == 11


# ---------------------------------------------------------------------------
# Run-set lint (round-4 build)
# ---------------------------------------------------------------------------


@pytest.fixture
def r4_run_set(pairs, tmp_path, monkeypatch) -> Path:
    return build_r4(pairs, tmp_path / "rs", monkeypatch)


def _stim(run_set: Path, condition: str) -> Path:
    runs = [json.loads(x) for x in (run_set / "runs.jsonl").read_text().splitlines()]
    return run_set / next(r["stim_dir"] for r in runs if r["condition"] == condition)


def test_r4_run_set_lints_clean(r4_run_set):
    report = lint_run_set(r4_run_set)
    assert report.ok, report.problems
    assert report.files_checked > 0


@pytest.mark.parametrize(
    ("file", "text", "expect"),
    [
        ("traw.md", "HiFi reads agree\n", "truth-set vocabulary"),
        ("traw.md", "a mosaic site\n", "outside caption.md"),
        ("tcur.md", "inside a segdup\n", "truth-set vocabulary"),
        ("tcur.md", "tr_stutter_neg\n", "case-type name"),
        ("traw.md", "mate at chrU:400\n", "sentinel contig"),
    ],
)
def test_r4_lint_catches_planted_leaks(r4_run_set, file, text, expect):
    d = _stim(r4_run_set, "Traw-r4+Tcur")
    with (d / file).open("a") as f:
        f.write(text)
    report = lint_run_set(r4_run_set)
    assert any(expect in p for p in report.problems), report.problems


def test_r4_prompt_must_be_the_r4_prompt(r4_run_set):
    d = _stim(r4_run_set, "Traw-r4")
    (d / "prompt.md").write_text(stimuli.prompt_text(["traw.md"], "sequenced"))  # round 3
    report = lint_run_set(r4_run_set)
    assert any("differs from the uniform prompt" in p for p in report.problems)


# ---------------------------------------------------------------------------
# Sentinel contig
# ---------------------------------------------------------------------------


def _sentinel_bam(
    d: Path, *, length: int = 1000, mate_pos=(0, 0), on_sentinel: bool = False, seq: str = "N"
) -> tuple[str, str]:
    fa = d / "ref.fa"
    fa.write_text(f">chr1\n{'ACGT' * 500}\n>chrU\n{seq * length}\n")
    pysam.faidx(str(fa))
    header = {"HD": {"VN": "1.6"}, "SQ": [{"SN": "chr1", "LN": 2000}, {"SN": "chrU", "LN": length}]}
    bam = d / "s.bam"
    with pysam.AlignmentFile(str(bam), "wb", header=header) as out:
        for k, mp in enumerate(mate_pos):
            a = pysam.AlignedSegment(out.header)
            a.query_name, a.query_sequence, a.flag = f"q{k}", "ACGTACGTAC", 0x1 | 0x40
            a.reference_id = 1 if on_sentinel and k == 0 else 0
            a.reference_start, a.mapping_quality = 100 + k, 60
            a.cigartuples = [(0, 10)]
            a.next_reference_id, a.next_reference_start = 1, mp
            out.write(a)
    return str(bam), str(fa)


def test_sentinel_contig_clean(tmp_path, pairs):
    assert lint_sentinel_contig(*_sentinel_bam(tmp_path)) == []
    assert lint_sentinel_contig(str(pairs["bam"]), str(pairs["reference"])) == []
    assert lint_sentinel_contig(None) == []


def test_sentinel_contig_absent_is_fine(tmp_path):
    fa = tmp_path / "r.fa"
    fa.write_text(">chr1\nACGT\n")
    header = {"HD": {"VN": "1.6"}, "SQ": [{"SN": "chr1", "LN": 4}]}
    bam = tmp_path / "x.bam"
    with pysam.AlignmentFile(str(bam), "wb", header=header):
        pass
    assert lint_sentinel_contig(str(bam), str(fa)) == []


@pytest.mark.parametrize(
    ("kw", "expect"),
    [
        ({"length": 5000}, "length 5000"),
        ({"mate_pos": (3, 7)}, "distinct mate positions"),
        ({"mate_pos": (2000, 2000)}, "outside the contig"),
        ({"on_sentinel": True}, "aligned to the sentinel"),
        ({"seq": "A"}, "not all N"),
    ],
)
def test_sentinel_contig_problems(tmp_path, kw, expect):
    problems = lint_sentinel_contig(*_sentinel_bam(tmp_path, **kw))
    assert any(expect in p for p in problems), problems


# ---------------------------------------------------------------------------
# KS (strata 3 vs 4 POI offset)
# ---------------------------------------------------------------------------


def _brute_force_p(a: list[float], b: list[float]) -> float:
    """Exact permutation p-value of the KS statistic (all relabellings of the pooled data)."""
    d_obs, _ = ks_two_sample(a, b)
    pool = a + b
    n = len(a)
    hits = total = 0
    for idx in itertools.combinations(range(len(pool)), n):
        s = set(idx)
        x = [pool[i] for i in idx]
        y = [pool[i] for i in range(len(pool)) if i not in s]
        total += 1
        hits += ks_two_sample(x, y)[0] >= d_obs - 1e-12
    return hits / total


def test_ks_separated_samples():
    d, p = ks_two_sample([1, 2, 3], [4, 5, 6])
    assert d == 1.0
    assert p == pytest.approx(2 / 20)  # 2 of C(6, 3) orderings are fully separated


def test_ks_identical_samples():
    assert ks_two_sample([1, 2, 3, 4], [1, 2, 3, 4]) == (0.0, 1.0)


@pytest.mark.parametrize("seed", [1, 2, 3])
def test_ks_exact_matches_permutation_without_ties(seed):
    rng = random.Random(seed)  # noqa: S311
    a = [rng.random() for _ in range(6)]
    b = [rng.random() + 0.3 for _ in range(5)]
    assert ks_two_sample(a, b)[1] == pytest.approx(_brute_force_p(a, b), abs=1e-9)


def test_ks_asymptotic_close_to_exact_for_large_samples():
    rng = random.Random(7)  # noqa: S311
    a = [rng.gauss(0, 1) for _ in range(300)]
    b = [rng.gauss(0.15, 1) for _ in range(300)]
    d, p_exact = ks_two_sample(a, b)
    assert leaklint._ks_asymptotic_p(300, 300, d) == pytest.approx(p_exact, abs=0.02)
    big = ks_two_sample(a * 2, b * 2)  # n + m > 1000: asymptotic branch
    assert 0.0 <= big[1] <= 1.0


def test_ks_rejects_empty():
    with pytest.raises(ValueError):
        ks_two_sample([], [1.0])


@pytest.mark.parametrize(
    ("a", "b"),
    [
        ([0, 0, 0, -1, -1, 0], [-1, -1, -1, 0, -1]),
        ([0, 1, 1, 2, 0, -1], [2, 2, 1, 3, 2, 0]),
        ([0, 0, 0, 0, 0, 0, -1], [0, -1, -1, -1, 0, -1]),
    ],
)
def test_ks_exact_matches_permutation_with_ties(a, b):
    # Integer POI offsets are heavily tied; the p-value must be the permutation one.
    assert ks_two_sample(a, b)[1] == pytest.approx(_brute_force_p(a, b), abs=1e-9)


def test_lint_poi_offsets_ties_do_not_inflate_p():
    # 15 vs 15 offsets in {-1, 0} (insertion anchor vs first deleted base): the permutation
    # p is about 0.14, so the draw must be redone; the continuous-data null gave 0.39.
    pos = [0] * 10 + [-1] * 5
    neg = [0] * 5 + [-1] * 10
    res = lint_poi_offsets(pos, neg)
    assert res.p == pytest.approx(0.1431, abs=0.002)
    assert not res.ok and res.problems


def test_lint_poi_offsets():
    same = lint_poi_offsets([0, 1, 2, 3, 0, 1, 2, 3], [0, 1, 2, 3, 1, 2, 0, 3])
    assert same.ok and not same.problems and same.p > 0.2
    apart = lint_poi_offsets([0] * 15, [5] * 15)
    assert not apart.ok and apart.d == 1.0 and "redraw by seed" in apart.problems[0]


def test_poi_offset():
    assert poi_offset(10, (10, 20)) == 0
    assert poi_offset(9, (10, 20)) == -1  # left-normalised insertion anchor before the run
    assert poi_offset(19, (10, 20)) == 9
    with pytest.raises(ValueError):
        poi_offset(20, (10, 20))


# ---------------------------------------------------------------------------
# Anchoring (strata 7-8)
# ---------------------------------------------------------------------------


def test_anchor_window():
    assert anchor_window(100, None) == (98, 103)
    assert anchor_window(100, (95, 110)) == (95, 110)


def test_lint_anchor_collisions():
    rows = [
        {"id": "clear", "site0": 100, "run": None, "other_truth": [[106, 109]]},
        {"id": "near", "site0": 100, "run": None, "other_truth": [[102, 105]]},
        {"chrom": "chr20", "site0": 200, "run": [195, 215], "other_truth": [[212, 214]]},
    ]
    problems = lint_anchor_collisions(rows)
    assert len(problems) == 2
    assert problems[0].startswith("near:") and problems[1].startswith("chr20:201:")


def test_mate_position_outside_its_contig_is_flagged(tmp_path):
    # A real (un-re-indexed) mate coordinate on a fixture contig, with or without chrU.
    bam, fa = _sentinel_bam(tmp_path)
    with pysam.AlignmentFile(bam, "rb") as src:
        header = src.header.to_dict()
    leaky = tmp_path / "leaky.bam"
    with pysam.AlignmentFile(str(leaky), "wb", header=header) as out:
        a = pysam.AlignedSegment(out.header)
        a.query_name, a.query_sequence, a.flag = "q9", "ACGTACGTAC", 0x1 | 0x40
        a.reference_id, a.reference_start, a.mapping_quality = 0, 100, 60
        a.cigartuples = [(0, 10)]
        a.next_reference_id, a.next_reference_start = 0, 31_000_000
        out.write(a)
    problems = lint_sentinel_contig(str(leaky), fa)
    assert any("mate position outside its contig" in p for p in problems), problems


def test_r4_caption_must_use_the_r4_sample_line(r4_run_set):
    d = _stim(r4_run_set, "Traw-r4")
    cap = d / "caption.md"
    cap.write_text(cap.read_text().replace(stimuli.SAMPLE_LINE_R4, real.SAMPLE_LINE))
    report = lint_run_set(r4_run_set)
    assert any("SAMPLE_LINE_R4" in p for p in report.problems), report.problems
