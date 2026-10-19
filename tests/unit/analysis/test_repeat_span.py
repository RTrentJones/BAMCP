"""Evidence-based homopolymer risk for repeat indels (eval round 4).

A one-unit indel in a long homopolymer is judged from the reads that span the whole run:
the length-only rule flagged every true repeat indel, and per-column counts make a het
insertion look like a ~30% minority because reads ending inside the run cannot show it.
"""

from __future__ import annotations

import pytest

from bamcp.analysis.curation import format_curation_summary
from bamcp.analysis.evidence import (
    compute_artifact_risk,
    interpret_variant,
    repeat_span_evidence,
)
from bamcp.core.parsers import AlignedRead

pytestmark = pytest.mark.unit

# Reference: 20 flank bases, a 12bp T run, 20 flank bases; the region starts at 1000.
START = 1000
FLANK = "ACGTCAGGACATGCAGTCAG"
REF = FLANK + "T" * 12 + "GACCTAGCATGCAGGTCAAC"
RUN0 = START + len(FLANK)  # first base of the run
RUN_END = RUN0 + 12
ANCHOR = RUN0 - 1  # VCF anchor of an indel in the run (0-based)
A = REF[len(FLANK) - 1]  # the anchor base
INS = {"position": ANCHOR, "ref": A, "alt": A + "T"}
DEL = {"position": ANCHOR, "ref": A + "T", "alt": A}


def _read(
    i: int,
    delta: int,
    *,
    start: int = START,
    end: int = START + len(REF),
    inserted: str = "T",
) -> AlignedRead:
    """A read spanning ``[start, end)`` that carries a ``delta``-bp length change in the run
    (an insertion is made of ``inserted`` bases at the run's left edge)."""
    if delta > 0:
        left = RUN0 - start
        cigar = f"{left}M{delta}I{end - RUN0}M"
    elif delta < 0:
        left = RUN0 - start
        cigar = f"{left}M{-delta}D{end - RUN0 + delta}M"
    else:
        cigar = f"{end - start}M"
    qlen = (end - start) + max(delta, 0)
    sequence = "A" * qlen
    if delta > 0:
        left = RUN0 - start
        sequence = sequence[:left] + inserted * delta + sequence[left + delta :]
    return AlignedRead(
        name=f"r{i}",
        sequence=sequence,
        qualities=[30] * qlen,
        cigar=cigar,
        position=start,
        end_position=end,
        mapping_quality=60,
        is_reverse=bool(i % 2),
    )


def _reads(**counts: int) -> list[AlignedRead]:
    """``ref=``, ``ins=``, ``dels=``, ``ins2=`` spanning reads plus ``short=`` reads that end
    inside the run (which must not count)."""
    out: list[AlignedRead] = []
    for kind, n in counts.items():
        for _ in range(n):
            i = len(out)
            if kind == "short":
                out.append(_read(i, 0, end=RUN0 + 5))
            else:
                out.append(_read(i, {"ref": 0, "ins": 1, "ins2": 2, "dels": -1}[kind]))
    return out


def _risk(variant: dict, reads: list[AlignedRead] | None) -> dict:
    v = {**variant, "depth": 60, "vaf": 0.3}
    evidence = {"strand_bias": 0.0, "quality_histogram": [], "mean_quality": 30}
    return compute_artifact_risk(v, evidence, REF, START, reads=reads)


def test_spanning_evidence_counts_only_reads_that_cover_the_run():
    span = repeat_span_evidence(INS, _reads(ref=10, ins=10, dels=2, short=20), REF, START)
    assert span is not None
    assert (span["run_start"], span["run_length"], span["run_base"]) == (RUN0, 12, "T")
    assert span["spanning_reads"] == 22  # the 20 reads ending inside the run are excluded
    assert span["allele_reads"] == 10 and span["allele_fraction"] == pytest.approx(
        10 / 22, abs=1e-3
    )
    assert span["other_length_fraction"] == pytest.approx(2 / 22, abs=1e-3)
    assert span["length_counts"] == {"-1": 2, "0": 10, "1": 10}


def test_spanning_evidence_deletion_and_snv():
    span = repeat_span_evidence(DEL, _reads(ref=8, dels=6), REF, START)
    assert span is not None and span["allele_reads"] == 6
    assert repeat_span_evidence({"position": 1005, "ref": "A", "alt": "G"}, [], REF, START) is None
    assert repeat_span_evidence(INS, [], None, START) is None


def test_well_supported_repeat_indel_is_not_flagged():
    # Het insertion: half the spanning reads, little other-length noise.
    risk = _risk(INS, _reads(ref=10, ins=10, dels=1, short=25))
    assert not [r for r in risk["risks"] if r["type"] == "homopolymer"]
    assert risk["repeat_span"]["allele_fraction"] == pytest.approx(10 / 21, abs=1e-3)


def test_rare_length_is_flagged_as_slippage():
    # Stutter: 1 of 20 spanning reads (5%) -> HIGH.
    risk = _risk(DEL, _reads(ref=19, dels=1))
    hp = [r for r in risk["risks"] if r["type"] == "homopolymer"]
    assert hp and hp[0]["severity"] == "high" and "1/20" in hp[0]["description"]
    # 3/20 (15%) -> MEDIUM.
    hp = [r for r in _risk(DEL, _reads(ref=17, dels=3))["risks"] if r["type"] == "homopolymer"]
    assert hp and hp[0]["severity"] == "medium"


def test_length_no_commoner_than_other_lengths_is_flagged():
    # 25% at the called length but 20% at other lengths: slippage spread, not a clear allele.
    risk = _risk(INS, _reads(ref=11, ins=5, ins2=2, dels=2))
    assert [r for r in risk["risks"] if r["type"] == "homopolymer"]


def test_without_reads_or_too_few_spanning_reads_the_length_rule_applies():
    for reads in (None, _reads(ref=3, ins=3, short=30)):
        hp = [r for r in _risk(INS, reads)["risks"] if r["type"] == "homopolymer"]
        assert hp and hp[0]["description"] == "In 12bp homopolymer run"


def test_interpretation_reads_zygosity_from_spanning_reads():
    variant = {
        "contig": "chr1",
        **INS,
        "alt_count": 18,
        "depth": 60,
        "vaf": 0.30,
        "artifact_risk": _risk(INS, _reads(ref=10, ins=10, short=25)),
    }
    text = interpret_variant(variant, one_based=False)
    assert "10/20 of the reads spanning the 12bp run (50%)" in text
    assert "consistent with heterozygous" in text and "below the germline" not in text


def test_curation_text_has_a_repeat_run_block():
    risk = _risk(INS, _reads(ref=10, ins=10, dels=1, short=25))
    summary = {
        "variant": {
            "location": "chr1:1020",
            "change": "C>CT",
            "vaf": 0.3,
            "alt_count": 10,
            "depth": 46,
        },
        "quality_metrics": {"mean_quality": None},
        "strand_analysis": {},
        "position_in_read": {},
        "artifact_assessment": risk,
        "confidence": "medium",
        "recommendations": [],
    }
    text = format_curation_summary(summary)
    assert "REPEAT RUN" in text
    assert "21 reads span the whole run" in text
    assert "This length on 10/21 spanning reads (47.6%)" in text
    assert "Lengths among spanning reads: -1: 1, ref: 10, +1: 10" in text
    assert "From spanning reads: consistent with heterozygous" in text


def test_recommendations_name_the_homopolymer_concern():
    from bamcp.analysis.curation import generate_curator_recommendations

    variant = {**DEL, "depth": 60, "vaf": 0.05}
    evidence = {"strand_bias": 0.0, "mapq_histogram": []}
    recs = generate_curator_recommendations(variant, evidence, _risk(DEL, _reads(ref=19, dels=1)))
    assert any(r.startswith("HOMOPOLYMER SLIPPAGE") for r in recs)
    recs = generate_curator_recommendations(variant, evidence, _risk(INS, None))
    assert any(r.startswith("HOMOPOLYMER CONTEXT") for r in recs)
    clean = _risk(INS, _reads(ref=10, ins=10))
    recs = generate_curator_recommendations({**INS, "depth": 60}, evidence, clean)
    assert not any(r.startswith("HOMOPOLYMER") for r in recs)


def test_alleles_that_are_not_a_pure_run_length_change_get_no_span_evidence():
    # Deleting the run plus the next base (here 12 T + G) is not a run-length change.
    nxt = REF[len(FLANK) + 12]
    crossing = {"position": ANCHOR, "ref": A + "T" * 12 + nxt, "alt": A}
    assert repeat_span_evidence(crossing, _reads(ref=10), REF, START) is None
    wrong_base = {"position": ANCHOR, "ref": A, "alt": A + "G"}
    assert repeat_span_evidence(wrong_base, _reads(ref=10), REF, START) is None
    # ...and they keep the run-length rule.
    hp = [r for r in _risk(crossing, _reads(ref=10))["risks"] if r["type"] == "homopolymer"]
    assert hp and "repeat_span" not in _risk(crossing, _reads(ref=10))


def test_run_at_the_window_edge_gets_no_span_evidence():
    # Same run, but the reference window ends 3 bases into it: its extent is unknown.
    clipped = REF[: len(FLANK) + 3]
    assert repeat_span_evidence(INS, _reads(ref=10, ins=10), clipped, START) is None
    starts_in_run = REF[len(FLANK) + 2 :]
    shifted = {**INS, "position": RUN0 + 1}
    assert repeat_span_evidence(shifted, [], starts_in_run, RUN0 + 2) is None


def test_an_insertion_of_another_base_is_not_the_called_allele():
    reads = _reads(ref=1) + [_read(100 + i, 1, inserted="G") for i in range(19)] + _reads(ins=1)
    span = repeat_span_evidence(INS, reads, REF, START)
    assert span is not None
    assert span["allele_reads"] == 1 and span["spanning_reads"] == 21
    assert span["length_counts"] == {"0": 1, "1": 1, "other": 19}
    hp = [r for r in _risk(INS, reads)["risks"] if r["type"] == "homopolymer"]
    assert hp and hp[0]["severity"] == "high"


def test_an_indel_outside_the_run_gets_no_span_evidence():
    from bamcp.analysis.evidence import homopolymer_run_bounds

    # A T insertion one base before the run's anchor: the run-length lookup reaches the T
    # run, but the insertion (between two flank bases) does not change it.
    b = REF[len(FLANK) - 2]
    early = {"position": ANCHOR - 1, "ref": b, "alt": b + "T"}
    assert homopolymer_run_bounds(early, REF, ANCHOR - 1 - START) == (RUN0 - START, RUN_END - START)
    assert repeat_span_evidence(early, _reads(ref=10), REF, START) is None
    # A deletion reaching past the run's end is not inside it.
    long_del = {"position": RUN_END - 3, "ref": "T" * 3 + "G", "alt": "T"}
    assert repeat_span_evidence(long_del, _reads(ref=10), REF, START) is None


def test_a_deletion_that_also_removes_a_flanking_base_is_not_the_called_allele():
    # Reads deleting the anchor base plus one T: same in-run length as the called -1T, but
    # a different allele.
    crossing = []
    for i in range(19):
        crossing.append(
            AlignedRead(
                name=f"x{i}",
                sequence="A" * (len(REF) - 2),
                qualities=[30] * (len(REF) - 2),
                cigar=f"{RUN0 - START - 1}M2D{len(REF) - (RUN0 - START) - 1}M",
                position=START,
                end_position=START + len(REF),
                mapping_quality=60,
                is_reverse=False,
            )
        )
    span = repeat_span_evidence(DEL, crossing + _reads(dels=1), REF, START)
    assert span is not None and span["allele_reads"] == 1
    assert span["length_counts"].get("other") == 19


def test_too_few_spanning_reads_keep_the_vaf_hint():
    risk = _risk(INS, _reads(ins=1, short=30))
    assert risk["repeat_span"]["spanning_reads"] == 1
    variant = {
        "contig": "chr1",
        **INS,
        "alt_count": 50,
        "depth": 100,
        "vaf": 0.5,
        "artifact_risk": risk,
    }
    text = interpret_variant(variant, one_based=False)
    assert "consistent with heterozygous" in text and "spanning" not in text
    summary = {
        "variant": {
            "location": "chr1:1020",
            "change": "C>CT",
            "vaf": 0.5,
            "alt_count": 50,
            "depth": 100,
        },
        "quality_metrics": {"mean_quality": None},
        "strand_analysis": {},
        "position_in_read": {},
        "artifact_assessment": risk,
        "confidence": "medium",
        "recommendations": [],
    }
    out = format_curation_summary(summary)
    assert "Too few spanning reads" in out and "From spanning reads" not in out
