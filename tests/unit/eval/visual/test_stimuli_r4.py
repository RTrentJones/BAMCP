"""Round-4 stimuli: ``traw_r4_text`` (mate fields, 5' distance, linkage table), ``tcur_text``
(J5) and the r4 conditions.

Hermetic: one hand-built paired-end BAM whose expected numbers are worked out below.
"""

from __future__ import annotations

import hashlib
import json
import random
from pathlib import Path
from typing import Any

import pysam
import pytest

from bamcp.eval.visual import stimuli
from bamcp.eval.visual.fixtures import Item

pytestmark = pytest.mark.unit

CHROM = "chr1"
LEN = 4000
P = 2000  # 0-based position of interest
START = P - 50  # every read starts at START + (i % 5): all cover P +/- 15


def _ref() -> str:
    rng = random.Random(11)  # noqa: S311
    seq: list[str] = []
    while len(seq) < LEN:
        b = rng.choice("ACGT")
        if len(seq) >= 2 and seq[-1] == seq[-2] == b:
            continue
        seq.append(b)
    return "".join(seq)


REF = _ref()


def _other(b: str) -> str:
    return {"A": "C", "C": "G", "G": "T", "T": "A"}[b]


ALT = _other(REF[P])
LINKED = _other(REF[P + 5])


def _make(
    start: int,
    *,
    mism: dict[int, str] | None = None,
    deletion: tuple[int, int] | None = None,
    insertion: tuple[int, str] | None = None,
    end_before: int | None = None,
    clip: int = 0,
    length: int = 100,
) -> tuple[str, list[tuple[int, int]]]:
    """(sequence, CIGAR) of a read aligned from ``start``: mismatches, one deletion, one
    insertion (after the given base), and an optional right soft clip of ``clip`` bases after
    the last aligned base ``end_before - 1``."""
    mism = mism or {}
    seq: list[str] = []
    cig: list[tuple[int, int]] = []

    def push(op: int, n: int) -> None:
        if cig and cig[-1][0] == op:
            cig[-1] = (op, cig[-1][1] + n)
        else:
            cig.append((op, n))

    pos = start
    aligned = length - clip
    while len(seq) < aligned:
        if end_before is not None and pos >= end_before:
            break
        if deletion and pos == deletion[0]:
            push(2, deletion[1])
            pos += deletion[1]
            continue
        seq.append(mism.get(pos, REF[pos]))
        push(0, 1)
        if insertion and pos == insertion[0]:
            seq.extend(insertion[1])
            push(1, len(insertion[1]))
        pos += 1
    if clip:
        seq.extend("A" * clip)
        push(4, clip)
    return "".join(seq), cig


# name, sequence/CIGAR spec, flag, mapq, MQ tag (None = absent), mate (contig, pos) or None
def _reads() -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    paired = 0x1 | 0x40
    proper = paired | 0x2
    alt_specs = [
        # a0, a1: proper pair, MQ 60, linked mismatch at P+5
        ({"mism": {P: ALT, P + 5: LINKED}}, proper, 60, (CHROM, P + 150)),
        ({"mism": {P: ALT, P + 5: LINKED}}, proper, 60, (CHROM, P - 200)),
        # a2, a3: mate on the sentinel contig, MQ 5, linked mismatch; a3 soft-clipped
        ({"mism": {P: ALT, P + 5: LINKED}}, paired, 5, ("chrU", 0)),
        (
            {"mism": {P: ALT, P + 5: LINKED}, "end_before": P + 8, "clip": 20},
            paired,
            5,
            ("chrU", 0),
        ),
        # a4: proper, 1 bp deletion at P-6; a5: mate unmapped, no MQ tag, same deletion
        ({"mism": {P: ALT}, "deletion": (P - 6, 1)}, proper, 60, (CHROM, P + 120)),
        ({"mism": {P: ALT}, "deletion": (P - 6, 1)}, paired | 0x8, None, None),
    ]
    ref_specs = [
        ({"mism": {P + 5: LINKED}}, proper, 30, (CHROM, P + 130)),  # r0
        ({}, paired | 0x8, None, None),  # r1 mate unmapped
        ({"mism": {P - 10: _other(REF[P - 10])}}, proper, 60, (CHROM, P + 100)),  # r2 singleton
        # r3, r4: insertion of TT after P+12; r5..r7 plain
        *[({"insertion": (P + 12, "TT")}, proper, 60, (CHROM, P + 110)) for _ in range(2)],
        *[({}, proper, 60, (CHROM, P + 110)) for _ in range(3)],
        ({}, 0, None, None),  # r8 unpaired
        ({}, 0, None, None),  # r9 unpaired
    ]
    i = 0
    for prefix, specs in (("a", alt_specs), ("r", ref_specs)):
        for k, (spec, flag, mq, mate) in enumerate(specs):
            start = START + (i % 5)
            seq, cig = _make(start, **spec)
            f = flag
            if i % 2:
                f |= 0x10
            elif f & 0x1 and not f & 0x8:
                f |= 0x20
            rows.append(
                {
                    "name": f"{prefix}{k}",
                    "start": start,
                    "seq": seq,
                    "cigar": cig,
                    "flag": f,
                    "mapq": 60,
                    "mq": mq,
                    "mate": mate,
                }
            )
            i += 1
    return rows


def write_pair_fixture(d: Path) -> dict[str, Path]:
    """Reference (``chr1`` + all-N ``chrU``) and the paired BAM described above."""
    fa = d / "ref.fa"
    fa.write_text(f">{CHROM}\n{REF}\n>chrU\n{'N' * 1000}\n")
    pysam.faidx(str(fa))
    header = {
        "HD": {"VN": "1.6", "SO": "coordinate"},
        "SQ": [{"SN": CHROM, "LN": LEN}, {"SN": "chrU", "LN": 1000}],
    }
    bam = d / "pairs.bam"
    rows = sorted(_reads(), key=lambda r: (r["start"], r["name"]))
    with pysam.AlignmentFile(str(bam), "wb", header=header) as out:
        for r in rows:
            a = pysam.AlignedSegment(out.header)
            a.query_name = r["name"]
            a.query_sequence = r["seq"]
            a.flag = r["flag"]
            a.reference_id = 0
            a.reference_start = r["start"]
            a.mapping_quality = r["mapq"]
            a.cigartuples = r["cigar"]
            a.query_qualities = pysam.qualitystring_to_array("I" * len(r["seq"]))
            if r["mate"] is not None:
                a.next_reference_id = out.get_tid(r["mate"][0])
                a.next_reference_start = r["mate"][1]
            elif r["flag"] & 0x8:
                a.next_reference_id = 0
                a.next_reference_start = r["start"]
            if r["mq"] is not None:
                a.set_tag("MQ", r["mq"])
            out.write(a)
    pysam.index(str(bam))
    return {"bam": bam, "reference": fa}


@pytest.fixture(scope="module")
def pairs(tmp_path_factory: pytest.TempPathFactory) -> dict[str, Path]:
    return write_pair_fixture(tmp_path_factory.mktemp("pairs"))


def _item(sites: list[int], item_id: str = "iabc") -> Item:
    return Item(item_id, "anchor_het_snv", 0, CHROM, sites, 16)


def _sha(text: str) -> str:
    return hashlib.sha256(text.encode()).hexdigest()


# sha256 of round-3 ``traw_text`` on this fixture, computed with the round-3 generator (base
# commit 51c8616) before the r4 refactor: the frozen Traw panel must stay byte-identical.
ROUND3_TRAW_SHA = {
    "snv": "58cc6f7919973a5378b3d2f36013e6aa844b6dd2543ed4f1540104732a42f4cb",
    "two_sites": "60ab508f96f5deb744f8bf8e939bf2e4005b2074d74b3c90fbf328428ddea99c",
    "del": "6229a3391663514d9f3f6184837eb4e91b4633828367a2232e67d083e8e0b5b6",
    "ins": "61d5f1c0e3458a4f4e01cc049782498ca552ac5dcd848fc6e00135de5356964d",
}
ROUND3_SITES = {"snv": [P], "two_sites": [P, P + 5], "del": [P - 6], "ins": [P + 12]}
# sha256 of the round-3 prompt for (files, sample_kind), same provenance.
ROUND3_PROMPT_SHA = {
    ("traw.md",): {
        "simulated": "295775806abf2bfb4923dbd90e35e4a644f5c1a54962a15230cc7ff02a1972ad",
        "sequenced": "5d20b250eba78daab81894ff574c7925cb2f67d9f2f9a4ad4d2dfa564953ff96",
    },
    ("traw.md", "view_wide.png", "view_zoom.png"): {
        "simulated": "6a6de2065c6863ffaaaf47b78dda48268310f990792d299ce0763a34572226c3",
        "sequenced": "5643fc634d58e8e3fc3b1886ff647f825272ae633f0598bf8c031da8f0b3ac66",
    },
}


def _paths(pairs: dict[str, Path]) -> tuple[str, str]:
    return str(pairs["bam"]), str(pairs["reference"])


@pytest.mark.parametrize("case", sorted(ROUND3_TRAW_SHA))
def test_round3_traw_is_byte_identical(pairs, case):
    text = stimuli.traw_text(*_paths(pairs), _item(ROUND3_SITES[case]))
    assert _sha(text) == ROUND3_TRAW_SHA[case]


@pytest.mark.parametrize("files", sorted(ROUND3_PROMPT_SHA))
@pytest.mark.parametrize("kind", ["simulated", "sequenced"])
def test_round3_prompt_is_byte_identical(files, kind):
    assert _sha(stimuli.prompt_text(list(files), kind)) == ROUND3_PROMPT_SHA[files][kind]
    assert stimuli.prompt_text(list(files), kind, r4=False) == stimuli.prompt_text(
        list(files), kind
    )


R4_ONLY_LINE_PREFIXES = (
    "Mate fields come from",
    "- distance (bp) from this position to the 5' end",
    "- reads with soft-clipped bases at their 5' end",
    "- read pairing:",
    "- mate mapping quality",
)


def _strip_r4(text: str) -> str:
    """Round-4 text minus its additions (r4-only lines and the linkage sections)."""
    out: list[str] = []
    skipping = False
    for line in text.split("\n"):
        if line.startswith("## Co-occurrence"):
            # the blank line before the first section stands in for the one this skips
            # before the section that follows the linkage tables
            skipping = True
            continue
        if skipping and line.startswith("## ") and not line.startswith("## Co-occurrence"):
            skipping = False
        if skipping or line.startswith(R4_ONLY_LINE_PREFIXES):
            continue
        out.append(line)
    return "\n".join(out).rstrip("\n") + "\n"


@pytest.mark.parametrize("case", sorted(ROUND3_SITES))
def test_r4_panel_is_round3_panel_plus_additions(pairs, case):
    item = _item(ROUND3_SITES[case])
    r3 = stimuli.traw_text(*_paths(pairs), item)
    r4 = stimuli.traw_r4_text(*_paths(pairs), item)
    assert r4 != r3
    assert _strip_r4(r4) == r3


def _section(text: str, header: str) -> list[str]:
    lines = text.split("\n")
    i = lines.index(header)
    out = []
    for line in lines[i + 1 :]:
        if line.startswith("#"):
            break
        out.append(line)
    return out


def test_r4_mate_fields_and_5prime_distance(pairs):
    """Hand counts. Alt G (6 reads): a0, a1, a4 proper with MQ 60; a2, a3 mate on chrU with MQ
    5; a5 mate unmapped, no MQ. a3 is reverse and soft-clipped on the right = its 5' end, 7 bp
    from the site. Ref C (10 reads): r0 MQ 30, r1 mate unmapped (no MQ), r2-r7 MQ 60, r8/r9
    unpaired (no flag 0x1, no MQ)."""
    text = stimuli.traw_r4_text(*_paths(pairs), _item([P]))
    ref, alt = REF[P], ALT
    g = _section(text, f"### allele {alt} (6 reads)")
    assert (
        "- read pairing: 6 of 6 reads are paired; in a proper pair: 3/6 (0.500); mate "
        "unmapped: 1; mate aligned to another contig (not in this region): 2"
    ) in g
    assert (
        "- mate mapping quality (MQ tag): 0: 0 | 1-9: 2 | 10-29: 0 | 30-59: 0 | 60+: 3 | "
        "no MQ tag: 1" in g
    )
    assert (
        "- distance (bp) from this position to the 5' end of the aligned part of the read: "
        "0-4: 0 | 5-9: 1 | 10-24: 0 | 25+: 5"
    ) in g
    assert "- reads with soft-clipped bases at their 5' end: 1" in g
    assert "- reads with soft-clipped bases at their 3' end: 0" in g
    c = _section(text, f"### allele {ref} (10 reads)")
    assert (
        "- read pairing: 8 of 10 reads are paired; in a proper pair: 7/8 (0.875); mate "
        "unmapped: 1; mate aligned to another contig (not in this region): 0"
    ) in c
    assert (
        "- mate mapping quality (MQ tag): 0: 0 | 1-9: 0 | 10-29: 0 | 30-59: 1 | 60+: 6 | "
        "no MQ tag: 3" in c
    )


def test_r4_linkage_table_hand_counts(pairs):
    """Within +/-15 bp of the site: the 1 bp deletion at P-6 is in a4, a5 (2 of 6 alt, 0 of
    10 ref); the P+5 mismatch is the linked band a0-a3 plus r0 (4/6 vs 1/10); a3's soft clip
    after P+7 (1/6); the TT insertion after P+12 in r3, r4 (a3 ends at P+7, so 0/5 alt
    reads cover it, 2/10 ref); r2's P-10 mismatch is a single-read event."""
    text = stimuli.traw_r4_text(*_paths(pairs), _item([P]))
    ref, alt = REF[P], ALT
    rows = _section(text, f"## Co-occurrence with other non-reference events near {CHROM}:{P + 1}")
    assert f"Most frequent non-reference allele at {CHROM}:{P + 1}: `{alt}`." in rows[1]
    table = [r for r in rows if r.startswith("| ") and r[2].isdigit()]
    assert table == [
        f"| {P - 5} | deletion of 1 bp starting here | 2 | 2/6 (0.333) | 0/10 (0.000) |",
        f"| {P + 6} | base {REF[P + 5]}>{LINKED} | 5 | 4/6 (0.667) | 1/10 (0.100) |",
        f"| {P + 8} | soft-clipped bases after this base | 1 | 1/6 (0.167) | 0/10 (0.000) |",
        f"| {P + 13} | insertion of TT after this base | 2 | 0/5 (0.000) | 2/10 (0.200) |",
    ]
    assert f"among reads with `{ref}` (reference)" in "\n".join(rows)
    assert "Other events within 15 bp, each in a single read without that allele: 1." in rows


def test_r4_linkage_excludes_the_site_own_event(pairs):
    """At the deletion site the deletion itself is the allele, not a co-occurring event; the
    alt G column at P (in both deletion reads) is."""
    text = stimuli.traw_r4_text(*_paths(pairs), _item([P - 6]))
    rows = _section(text, f"## Co-occurrence with other non-reference events near {CHROM}:{P - 5}")
    assert "allele at chr1:1995: `-`." in rows[1]
    table = [r for r in rows if r.startswith("| ") and r[2].isdigit()]
    assert not any("deletion" in r for r in table)
    assert f"| {P + 1} | base {REF[P]}>{ALT} | 6 | 2/2 (1.000) | 4/14 (0.286) |" in table


def test_r4_linkage_always_printed_at_a_clean_site(pairs):
    text = stimuli.traw_r4_text(*_paths(pairs), _item([P + 40]))
    assert f"## Co-occurrence with other non-reference events near {CHROM}:{P + 41}" in text
    assert f"No read carries a non-reference allele at {CHROM}:{P + 41}." in text
    assert "| - | (no other non-reference event) | 0 | - | - |" in text


def test_r4_two_sites_get_a_table_each_and_keep_combinations(pairs):
    text = stimuli.traw_r4_text(*_paths(pairs), _item([P, P + 5]))
    assert text.count("## Co-occurrence with other non-reference events") == 2
    assert "## Allele combinations in reads covering all of" in text


def test_r4_unpaired_reads(tmp_path):
    """A BAM without any pairing information (round-3 fixtures) still renders."""
    fa = tmp_path / "ref.fa"
    fa.write_text(f">{CHROM}\n{REF}\n")
    pysam.faidx(str(fa))
    header = {"HD": {"VN": "1.6", "SO": "coordinate"}, "SQ": [{"SN": CHROM, "LN": LEN}]}
    bam = tmp_path / "se.bam"
    with pysam.AlignmentFile(str(bam), "wb", header=header) as out:
        for k in range(4):
            seq, cig = _make(START + k, mism={P: ALT} if k < 2 else None)
            a = pysam.AlignedSegment(out.header)
            a.query_name, a.query_sequence, a.flag = f"s{k}", seq, 16 * (k % 2)
            a.reference_id, a.reference_start, a.mapping_quality = 0, START + k, 60
            a.cigartuples = cig
            a.query_qualities = pysam.qualitystring_to_array("I" * len(seq))
            out.write(a)
    pysam.index(str(bam))
    text = stimuli.traw_r4_text(str(bam), str(fa), _item([P]))
    assert text.count("- read pairing: 0 of 2 reads are flagged as paired (no mate fields)") == 2
    assert "- mate mapping quality" not in text
    assert f"| {P + 1} |" not in "\n".join(
        _section(text, f"## Co-occurrence with other non-reference events near {CHROM}:{P + 1}")
    )


def test_r4_traw_has_no_label_vocabulary(pairs):
    from bamcp.eval.visual.leaklint import LEAK_VOCAB, r4_vocab_problems

    for sites in ROUND3_SITES.values():
        text = stimuli.traw_r4_text(*_paths(pairs), _item(sites))
        assert not LEAK_VOCAB.search(text)
        assert r4_vocab_problems("traw.md", text) == []


# ---------------------------------------------------------------------------
# Tcur (J5)
# ---------------------------------------------------------------------------


@pytest.mark.parametrize(
    ("site", "expect"),
    [
        (P, {"pos": P + 1, "ref": REF[P], "alt": ALT, "reads": 6}),
        # deletion of REF[P-6]: VCF anchor is the base before it
        (P - 6, {"pos": P - 6, "ref": REF[P - 7 : P - 5], "alt": REF[P - 7], "reads": 2}),
        (P + 12, {"pos": P + 13, "ref": REF[P + 12], "alt": REF[P + 12] + "TT", "reads": 2}),
        (P - 5, None),  # the base after the deletion: every read is reference there
    ],
)
def test_top_nonref_allele(pairs, site, expect):
    got = stimuli.top_nonref_allele(*_paths(pairs), CHROM, site)
    assert got == (None if expect is None else {"chrom": CHROM, **expect})


def test_top_nonref_allele_inside_a_deletion(pairs):
    """A deletion spanning (not starting at) the site is still that site's allele."""
    got = stimuli.top_nonref_allele(*_paths(pairs), CHROM, P - 6)
    assert got is not None and got["ref"] == REF[P - 7 : P - 5]


async def test_tcur_always_has_a_curation_block_at_the_poi(pairs):
    bam, ref = _paths(pairs)
    item = _item([P, P + 12, P - 6, P + 40])
    text = await stimuli.tcur_text(bam, ref, item)
    s0, e0 = stimuli.window(item, stimuli.WIDE_BP)
    assert text.startswith(f"## get_region_summary(region={CHROM}:{s0 + 1}-{e0})")
    blocks = text.split("## get_variant_curation_summary(")[1:]
    assert len(blocks) == 4
    # passes the calling filters: the shipped evidence path
    assert blocks[0].startswith(f"chrom={CHROM}, pos={P + 1}, ref={REF[P]}, alt={ALT})")
    assert "Passes the applied filters" in blocks[0]
    # 2 alt reads < min_alt_count: probed (probe_site semantics), filter_check reported
    for b in blocks[1:3]:
        assert "FAILS the applied filters: min_alt_count" in b
        assert "Applied filters: min_vaf=" in b
    assert blocks[1].startswith(f"chrom={CHROM}, pos={P + 13}, ref={REF[P + 12]}, alt=")
    assert blocks[2].startswith(f"chrom={CHROM}, pos={P - 6}, ref={REF[P - 7 : P - 5]}, ")
    assert blocks[3].startswith(f"chrom={CHROM}, pos={P + 41})")
    assert "not run: every read covering" in blocks[3]


async def test_tcur_reuses_the_tfull_machinery(pairs):
    """Same region summary as Tfull; at a called site, the same curation text."""
    bam, ref = _paths(pairs)
    item = _item([P])
    assert await stimuli.tcur_text(bam, ref, item) == await stimuli.tfull_text(bam, ref, item)


async def test_tcur_runs_where_tfull_says_not_called(pairs):
    bam, ref = _paths(pairs)
    item = _item([P - 10])  # r2's single-read mismatch
    assert "(not called" in await stimuli.tfull_text(bam, ref, item)
    text = await stimuli.tcur_text(bam, ref, item)
    assert f"pos={P - 9}, ref={REF[P - 10]}, alt=" in text
    assert "FAILS the applied filters" in text and "Applied filters: min_vaf=" in text


# ---------------------------------------------------------------------------
# Conditions, prompt, build
# ---------------------------------------------------------------------------


def test_r4_conditions():
    assert stimuli.condition_files("Traw-r4") == ["traw.md"]
    assert stimuli.condition_files("Traw-r4+Tcur") == ["traw.md", "tcur.md"]
    assert stimuli.condition_files("V+Traw-r4") == ["traw.md", "view_wide.png", "view_zoom.png"]
    assert "V+Traw-r4" in stimuli.VISUAL_CONDITIONS
    assert not set(stimuli.CONDITIONS_R4) & set(stimuli.CONDITIONS)
    assert "tcur.md" in stimuli.FILE_DESCRIPTIONS and "tcur.md" in stimuli.FILE_ORDER


def test_r4_prompt_swaps_the_anchoring_sentence():
    r4 = stimuli.prompt_text(["traw.md", "tcur.md"], "sequenced", r4=True)
    assert stimuli.ANCHOR_SENTENCE_R4 in r4
    assert stimuli.ANCHOR_SENTENCE not in r4
    assert stimuli.FILE_DESCRIPTIONS_R4["traw.md"] in r4
    assert stimuli.FILE_DESCRIPTIONS["tcur.md"] in r4
    r3 = stimuli.prompt_text(["traw.md"], "sequenced")
    assert stimuli.ANCHOR_SENTENCE in r3 and stimuli.ANCHOR_SENTENCE_R4 not in r3


def test_sample_line_r4_is_the_frozen_caption_line():
    assert stimuli.SAMPLE_LINE_R4 == (
        "sequenced human sample (Illumina short reads, ~60x whole-genome); the sample may "
        "carry germline variants or variants present in only a fraction of cells (mosaic); "
        "coordinates are re-indexed"
    )


def _label(item: Item) -> dict[str, Any]:
    return {
        "item": item.item_id,
        "case_type": item.case_type,
        "positions": item.positions_1based(),
        "variant_present": "yes",
        "stratum": item.case_type,
        "pair": "anchor",
        "polarity": "pos",
        "cluster": "bin:chr20:12",
        "truth_source": "v4.2.1+v5.0q",
        "knob": "default",
    }


def build_r4(pairs: dict[str, Path], out: Path, monkeypatch, conditions=None) -> Path:
    import bamcp.eval.renderer as renderer_mod

    from .conftest import FakeViewerRenderer

    async def _available() -> None:
        return None

    monkeypatch.setattr(renderer_mod, "ViewerRenderer", FakeViewerRenderer)
    monkeypatch.setattr(stimuli, "chromium_available", _available)
    monkeypatch.setattr(stimuli, "viewer_bundle_problem", lambda: None)
    items = [
        Item("i0000aaaa", "anchor_het_snv", 0, CHROM, [P], 16),
        Item("i0000bbbb", "tr_indel_pos", 1, CHROM, [P + 12], 16),
    ]
    stimuli.build_items_run_set(
        out,
        bam=str(pairs["bam"]),
        ref=str(pairs["reference"]),
        items=items,
        label_records=[_label(it) for it in items],
        seed=4,
        conditions=list(conditions or stimuli.CONDITIONS_R4),
        models=["m"],
        meta={"tier": "real"},
        sample_line=stimuli.SAMPLE_LINE_R4,
        sample_kind="sequenced",
        r4=True,
        log=lambda *_: None,
    )
    return out


def test_build_items_run_set_r4(pairs, tmp_path, monkeypatch):
    from bamcp.eval.visual.leaklint import lint_run_set

    out = build_r4(pairs, tmp_path / "rs", monkeypatch)
    runs = [json.loads(x) for x in (out / "runs.jsonl").read_text().splitlines()]
    assert sorted({r["condition"] for r in runs}) == sorted(stimuli.CONDITIONS_R4)
    meta = json.loads((out / "build.json").read_text())
    assert meta["r4"] is True and meta["conditions"] == list(stimuli.CONDITIONS_R4)
    labels = [json.loads(x) for x in (out / "labels" / "labels.jsonl").read_text().splitlines()]
    assert labels[0]["pair"] == "anchor" and labels[0]["truth_source"] == "v4.2.1+v5.0q"
    for r in runs:
        d = out / r["stim_dir"]
        names = sorted(p.name for p in d.iterdir())
        assert names == sorted(
            ["caption.md", "prompt.md", *stimuli.condition_files(r["condition"])]
        )
        if "traw.md" in names:
            assert "## Co-occurrence with other non-reference" in (d / "traw.md").read_text()
        else:
            assert (d / "mpileup.txt").read_text().startswith("# samtools mpileup")
        assert stimuli.SAMPLE_LINE_R4 in (d / "caption.md").read_text()
        assert stimuli.ANCHOR_SENTENCE_R4 in (d / "prompt.md").read_text()
        if r["condition"].endswith("+Tcur"):
            assert "## get_variant_curation_summary(" in (d / "tcur.md").read_text()
    report = lint_run_set(out)
    assert report.ok, report.problems


@pytest.mark.parametrize(
    ("conditions", "r4"), [(["Traw"], True), (["Traw-r4"], False), (["Tcur"], True)]
)
def test_r4_flag_and_conditions_must_agree(pairs, tmp_path, conditions, r4):
    with pytest.raises(ValueError, match="unknown condition"):
        stimuli.build_items_run_set(
            tmp_path,
            bam=str(pairs["bam"]),
            ref=str(pairs["reference"]),
            items=[],
            label_records=[],
            seed=1,
            conditions=conditions,
            r4=r4,
        )


def test_mpileup_text_covers_the_site_window_unfiltered(pairs):
    item = _item([P])
    text = stimuli.mpileup_text(*_paths(pairs), item, half=5)
    header, _, body = text.partition("\n\n")
    assert "-B -Q 0 -q 0 -s --output-BP" in header
    lines = [ln.split("\t") for ln in body.splitlines()]
    assert [int(f[1]) for f in lines] == list(range(P + 1 - 5, P + 1 + 6))  # 1-based
    site = next(f for f in lines if int(f[1]) == P + 1)
    assert len(site) == 8  # + mapping qualities (-s) and read positions (--output-BP)
    assert int(site[3]) == len(site[5]) == len(site[6])  # one quality per read


def test_tmp_conditions_list_mpileup_first():
    assert stimuli.condition_files("Tmp") == ["mpileup.txt"]
    assert stimuli.condition_files("Tmp+Tcur") == ["mpileup.txt", "tcur.md"]
    assert stimuli.condition_files("V+Tmp") == ["mpileup.txt", "view_wide.png", "view_zoom.png"]
    assert "V+Tmp" in stimuli.VISUAL_CONDITIONS
