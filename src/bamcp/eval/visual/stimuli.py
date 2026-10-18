"""Stimulus builder: one opaque directory per (item, condition).

Layout of a run set (``<out>``)::

    stim/<sha8>/            caption.md, prompt.md + the condition's evidence files
    labels/labels.jsonl     ground truth — never inside stim/, never shown to a model
    runs.jsonl              run manifest: run_id, item, condition, model, repeat, files, answer_path
    answers/                where model answers go (answers/<run_id>.json)
    build.json              seed, conditions, item count, whether PNGs were rendered

Conditions (pre-registered):

- ``Traw``     site pileup numbers only (no verdicts / risk / confidence / recommendations)
- ``Tjson``    the viewer's render payload, stripped of read names / variants / evidence
- ``Tfull``    the shipped ``get_region_summary`` + ``get_variant_curation_summary`` text
- ``V``        two blinded viewer PNGs (~200 bp and ~60 bp)
- ``V+Traw``   V plus Traw
- ``Tfull+V``  Tfull plus V

Round-4 conditions (``docs/eval/visual-design-r4.md`` section 4; built with ``r4=True``):

- ``Traw-r4``       the Traw panel plus mate fields, 5' distance and a linkage table (B)
- ``Traw-r4+Tcur``  B plus BAMCP's text tools at the POI's top non-reference allele (C1, J5)
- ``V+Traw-r4``     B plus the two viewer PNGs (C2)
"""

from __future__ import annotations

import asyncio
import json
import shutil
from collections import Counter, defaultdict
from dataclasses import dataclass
from pathlib import Path
from typing import Any

from .answers import ANSWER_SCHEMA_TEXT, FIELD_GUIDE
from .cases import CASE_BY_NAME
from .fixtures import Item, Plan, build_plan, ensure_fixtures, opaque_id

CONDITIONS = ("Traw", "Tjson", "Tfull", "V", "V+Traw", "Tfull+V")
PILOT_CONDITIONS = ("Traw", "Tjson", "V", "V+Traw")
# "Tmp*": round-4 raw-input arms -- the baseline is the samtools mpileup text a user without
# BAMCP would paste; BAMCP's text tools (tcur.md) or viewer are added on top.
CONDITIONS_R4 = (
    "Traw-r4",
    "Traw-r4+Tcur",
    "V+Traw-r4",
    "Tmp",
    "Tmp+Tcur",
    "V+Tmp",
)
ALL_CONDITIONS = CONDITIONS + CONDITIONS_R4
VISUAL_CONDITIONS = frozenset({"V", "V+Traw", "Tfull+V", "V+Traw-r4", "V+Tmp"})

WIDE_BP = 200
ZOOM_BP = 60
VIEW_WIDTH = 1024

FILE_ORDER = (
    "traw.md",
    "mpileup.txt",
    "tjson.json",
    "tfull.md",
    "tcur.md",
    "view_wide.png",
    "view_zoom.png",
)
MPILEUP_HALF = 30  # mpileup.txt covers the position(s) of interest +/- this many bp
FILE_DESCRIPTIONS = {
    "traw.md": (
        "Pileup statistics computed from every aligned read overlapping the position(s) "
        "of interest: per-allele read counts by alignment orientation, mapping-quality, "
        "base-quality and read-position histograms, and the reference context."
    ),
    "tjson.json": (
        "The genome viewer's data payload for the region: reads (start/end, CIGAR, "
        "sequence, base qualities, mapping quality, orientation, mismatches, soft clips), "
        "the reference sequence and per-base coverage. Positions in this file are "
        "0-based, half-open (add 1 to compare with the caption)."
    ),
    "tfull.md": (
        "Text returned by the BAMCP MCP tools get_region_summary and "
        "get_variant_curation_summary for this region."
    ),
    "tcur.md": (
        "Text returned by the BAMCP MCP tools at their default settings: get_region_summary "
        "for the region, and get_variant_curation_summary for the most frequent "
        "non-reference allele at the position(s) of interest."
    ),
    "mpileup.txt": (
        "Output of samtools mpileup for the reads around the position(s) of interest, with "
        "no base- or mapping-quality filter; the file's header explains the columns."
    ),
    "view_wide.png": "Screenshot of the read pileup viewer, ~200 bp centred on the site.",
    "view_zoom.png": "Screenshot of the read pileup viewer, ~60 bp centred on the site.",
}

# Round 4 describes traw.md's added fields; every other file keeps its round-3 description.
FILE_DESCRIPTIONS_R4 = {
    **FILE_DESCRIPTIONS,
    "traw.md": (
        "Pileup statistics computed from every aligned read overlapping the position(s) "
        "of interest: per-allele read counts by alignment orientation, mapping-quality, "
        "base-quality and read-position (3' and 5' end) histograms, mate-pair fields, a "
        "table of how often other non-reference events within 15 bp occur in the same "
        "reads, and the reference context."
    ),
}

VIEWER_KEY = """\
Viewer key (screenshots; describes the BAMCP viewer as rendered at build time): top to
bottom — coordinate ruler (1-based), reference sequence (coloured base tiles: A green,
C blue, G orange, T red), per-base coverage track (the dashed line marks 10x), then one row
per read, reads packed into rows. Every read is a light-grey bar whatever its strand; the
strand is shown only by a coloured stripe along the read's bottom edge and its coloured
pointed end: blue = aligned forward (points right), red = aligned in reverse (points left).
Reads with lower mapping quality are drawn paler (faded toward white); reads with mapping
quality 0 are white boxes with a dashed grey outline. At the wide zoom, read bases that
differ from the reference are vertical ticks in the read base's colour (A green, C blue,
G orange, T red, N dark grey) and matching bases are not marked. At base-level zoom every
read base is a lettered tile: bases that match the reference are grey tiles with dark
letters, bases that differ are tiles in the base's colour with white letters, so a
differing base stands out from the grey. Base quality: a differing base with base quality
below 20 is drawn hollow (a pale centre inside a border of the base's colour; at base zoom
its letter is also in the base's colour), while base quality 20 or above is drawn solid.
Tick and tile colours are not faded by mapping quality. A deletion is a gap in the read
bridged by a thin dark line; an insertion is a purple vertical mark between two bases.
Soft-clipped bases are not shown."""


@dataclass
class BuildResult:
    out: Path
    n_items: int
    n_stimuli: int
    n_runs: int
    pngs: bool
    skipped_conditions: list[str]


# ---------------------------------------------------------------------------
# Regions
# ---------------------------------------------------------------------------


def window(item: Item, width: int) -> tuple[int, int]:
    """0-based half-open window of ``width`` bp centred on the item."""
    start = item.center - width // 2
    return start, start + width


def region_1based(contig: str, start0: int, end0: int) -> str:
    return f"{contig}:{start0 + 1}-{end0}"


# ---------------------------------------------------------------------------
# Traw — raw pileup numbers
# ---------------------------------------------------------------------------

_MAPQ_BINS = ((0, 0, "0"), (1, 9, "1-9"), (10, 29, "10-29"), (30, 59, "30-59"), (60, 255, "60+"))
_BQ_BINS = ((0, 9, "<10"), (10, 19, "10-19"), (20, 29, "20-29"), (30, 39, "30-39"), (40, 99, "40+"))
_DIST_BINS = ((0, 4, "0-4"), (5, 9, "5-9"), (10, 24, "10-24"), (25, 10_000, "25+"))


def _bin(value: int, bins: tuple[tuple[int, int, str], ...]) -> str:
    for lo, hi, name in bins:
        if lo <= value <= hi:
            return name
    return bins[-1][2]


@dataclass(frozen=True)
class _Mate:
    """Pairing fields of one read (FLAG bits and the ``MQ`` tag)."""

    paired: bool
    proper: bool
    mate_unmapped: bool
    elsewhere: bool  # mate mapped to another contig (the fixtures' sentinel contig included)
    mq: int | None


@dataclass
class _Obs:
    allele: str
    reverse: bool
    mapq: int
    bq: int | None
    dist3: int | None
    clipped3: bool
    span: tuple[int, int, bool]
    dist5: int | None = None
    clipped5: bool = False
    mate: _Mate | None = None


def _mate_info(read: Any) -> _Mate:
    paired = bool(read.is_paired)
    unmapped = paired and bool(read.mate_is_unmapped)
    elsewhere = (
        paired
        and not unmapped
        and read.next_reference_id >= 0
        and read.next_reference_id != read.reference_id
    )
    mq = int(read.get_tag("MQ")) if paired and read.has_tag("MQ") else None
    return _Mate(paired, paired and bool(read.is_proper_pair), unmapped, elsewhere, mq)


def _read_observations(read: Any, positions: set[int]) -> dict[int, _Obs]:
    """Allele of one read at each of ``positions`` (0-based)."""
    seq = read.query_sequence or ""
    quals = read.query_qualities
    rev = bool(read.is_reverse)
    q_first, q_last = read.query_alignment_start, read.query_alignment_end - 1
    cig = read.cigartuples or []
    clipped3 = bool(cig) and ((cig[0][0] == 4) if rev else (cig[-1][0] == 4))
    clipped5 = bool(cig) and ((cig[-1][0] == 4) if rev else (cig[0][0] == 4))
    mate = _mate_info(read)
    span = (read.reference_start, read.reference_end, rev)
    out: dict[int, _Obs] = {}
    last_ref: int | None = None
    pending_ins: dict[int, str] = {}
    pairs = read.get_aligned_pairs()
    for qpos, rpos in pairs:
        if rpos is None:
            if qpos is not None and last_ref is not None and q_first <= qpos <= q_last:
                pending_ins[last_ref] = pending_ins.get(last_ref, "") + seq[qpos]
            continue
        last_ref = rpos
        if rpos not in positions:
            continue
        if qpos is None:
            out[rpos] = _Obs(
                "-", rev, read.mapping_quality, None, None, clipped3, span, None, clipped5, mate
            )
        else:
            dist3 = (qpos - q_first) if rev else (q_last - qpos)
            dist5 = (q_last - qpos) if rev else (qpos - q_first)
            bq = int(quals[qpos]) if quals is not None else None
            out[rpos] = _Obs(
                seq[qpos],
                rev,
                read.mapping_quality,
                bq,
                dist3,
                clipped3,
                span,
                dist5,
                clipped5,
                mate,
            )
    for rpos, ins in pending_ins.items():
        if rpos in out and out[rpos].allele != "-":
            out[rpos].allele = f"{out[rpos].allele}+{ins}"
    return out


# A non-reference event in one read: (0-based position, kind, detail). Deletions sit at their
# first deleted base, insertions at the base they follow, soft clips at the aligned base they
# touch ("before" = clipped bases precede it, "after" = they follow it).
_Event = tuple[int, str, str]
_EVENT_ORDER = {"mismatch": 0, "deletion": 1, "insertion": 2, "soft clip": 3}
_LINK_BP = 15  # linkage table: events within this many bp of the site


def _read_events(read: Any, lo: int, hi: int, ref_at: Any) -> set[_Event]:
    """Every non-reference event of ``read`` at a position in [lo, hi)."""
    seq = read.query_sequence or ""
    events: set[_Event] = set()
    cig = read.cigartuples or []
    if cig and cig[0][0] == 4 and lo <= read.reference_start < hi:
        events.add((read.reference_start, "soft clip", "before"))
    if cig and cig[-1][0] == 4 and lo <= read.reference_end - 1 < hi:
        events.add((read.reference_end - 1, "soft clip", "after"))
    pos = read.reference_start
    qpos = 0
    for op, n in cig:
        if op in (0, 7, 8):  # M / = / X
            for k in range(n):
                p = pos + k
                if lo <= p < hi:
                    base = seq[qpos + k].upper()
                    if base in "ACGT" and base != ref_at(p).upper():
                        events.add((p, "mismatch", f"{ref_at(p).upper()}>{base}"))
            pos += n
            qpos += n
        elif op == 1:  # I
            if lo <= pos - 1 < hi:
                events.add((pos - 1, "insertion", seq[qpos : qpos + n].upper()))
            qpos += n
        elif op in (2, 3):  # D / N
            if op == 2 and lo <= pos < hi:
                events.add((pos, "deletion", str(n)))
            pos += n
        elif op == 4:  # S
            qpos += n
    return events


def _event_label(ev: _Event) -> str:
    _, kind, detail = ev
    if kind == "mismatch":
        return f"base {detail}"
    if kind == "deletion":
        return f"deletion of {detail} bp starting here"
    if kind == "insertion":
        return f"insertion of {detail} after this base"
    return f"soft-clipped bases {detail} this base"


def _overlaps_site(ev: _Event, s: int) -> bool:
    """True for the site's own column: a base, insertion or deletion at ``s``, or a deletion
    spanning it. Soft-clip boundaries are never the site's own event."""
    p, kind, detail = ev
    if kind == "soft clip":
        return False
    if kind == "deletion":
        return p <= s < p + int(detail)
    return p == s


def _hist(values: list[str], bins: tuple[tuple[int, int, str], ...]) -> str:
    c = Counter(values)
    return " | ".join(f"{name}: {c.get(name, 0)}" for _, _, name in bins)


def _frac(k: int, n: int) -> str:
    return f"{k}/{n} ({k / n:.3f})" if n else "0/0"


def _mate_lines(g: list[_Obs]) -> list[str]:
    mates = [o.mate for o in g if o.mate is not None and o.mate.paired]
    if not mates:
        return [f"- read pairing: 0 of {len(g)} reads are flagged as paired (no mate fields)"]
    no_mq = sum(1 for o in g if o.mate is None or o.mate.mq is None)
    mq = _hist([_bin(m.mq, _MAPQ_BINS) for m in mates if m.mq is not None], _MAPQ_BINS)
    return [
        f"- read pairing: {len(mates)} of {len(g)} reads are paired; in a proper pair: "
        f"{_frac(sum(m.proper for m in mates), len(mates))}; mate unmapped: "
        f"{sum(m.mate_unmapped for m in mates)}; mate aligned to another contig "
        f"(not in this region): {sum(m.elsewhere for m in mates)}",
        f"- mate mapping quality (MQ tag): {mq} | no MQ tag: {no_mq}",
    ]


def _linkage_lines(
    contig: str,
    s: int,
    ref_base: str,
    reads: list[tuple[dict[int, _Obs], set[_Event], int, int]],
) -> list[str]:
    """The always-printed linkage table for site ``s``: how often each other non-reference
    event within the context window occurs in reads with the site's most frequent
    non-reference allele vs reads with the reference base, among the reads covering it."""
    at_site = Counter(obs[s].allele for obs, _, _, _ in reads if s in obs)
    nonref = sorted(
        ((n, a) for a, n in at_site.items() if a != ref_base), key=lambda t: (-t[0], t[1])
    )
    top = nonref[0][1] if nonref else None

    def others(evs: set[_Event]) -> list[_Event]:
        return [e for e in evs if abs(e[0] - s) <= _LINK_BP and not _overlaps_site(e, s)]

    carriers: Counter[_Event] = Counter()
    in_top: Counter[_Event] = Counter()
    for obs, evs, _, _ in reads:
        carriers.update(others(evs))
        if top is not None and s in obs and obs[s].allele == top:
            in_top.update(others(evs))
    shown = sorted(
        (e for e in carriers if carriers[e] >= 2 or in_top[e] >= 1),
        key=lambda e: (e[0], _EVENT_ORDER[e[1]], e[2]),
    )
    rest = [e for e in carriers if e not in shown]
    label = f"`{top}`" if top is not None else "(none)"
    lines = [
        "",
        f"## Co-occurrence with other non-reference events near {contig}:{s + 1}",
        "",
        f"Most frequent non-reference allele at {contig}:{s + 1}: {label}. For each other "
        f"non-reference event within {_LINK_BP} bp, the reads carrying it among the reads "
        "with that allele and among the reads with the reference base, counting only reads "
        "whose alignment covers the event's position. Listed: every event in two or more "
        "reads, and every event in a read with that allele.",
        "",
    ]
    if top is None:
        lines.append(f"No read carries a non-reference allele at {contig}:{s + 1}.")
        lines.append("")
    lines.append(
        f"| position | event | reads with event | among reads with {label} | "
        f"among reads with `{ref_base}` (reference) |"
    )
    lines.append("|---|---|---|---|---|")
    for e in shown:
        p = e[0]

        def count(allele: str | None, ev: _Event = e, p: int = p) -> str:
            if allele is None:
                return "n/a"
            cov = [
                evs
                for obs, evs, a, b in reads
                if s in obs and obs[s].allele == allele and a <= p < b
            ]
            return _frac(sum(1 for evs in cov if ev in evs), len(cov))

        lines.append(
            f"| {p + 1} | {_event_label(e)} | {carriers[e]} | {count(top)} | {count(ref_base)} |"
        )
    if not shown:
        lines.append("| - | (no other non-reference event) | 0 | - | - |")
    lines.append("")
    lines.append(
        f"Other events within {_LINK_BP} bp, each in a single read without that allele: "
        f"{len(rest)}."
    )
    return lines


def traw_text(bam_path: str, ref_path: str, item: Item) -> str:
    """The round-3 ``Traw`` panel (frozen: byte-identical to the round-3 generator)."""
    return _traw(bam_path, ref_path, item, r4=False)


def traw_r4_text(bam_path: str, ref_path: str, item: Item) -> str:
    """The round-4 ``Traw-r4`` panel (B): the round-3 panel plus, per allele, the 5'-end
    distance and 5' soft clips and the mate fields (proper pair, mate unmapped, mate on
    another contig, ``MQ`` histogram), and per site an always-printed linkage table."""
    return _traw(bam_path, ref_path, item, r4=True)


def _traw(bam_path: str, ref_path: str, item: Item, *, r4: bool) -> str:
    import pysam

    sites = item.sites
    context = 15
    table_pos: set[int] = set()
    for s in sites:
        table_pos.update(range(s - context, s + context + 1))
    with pysam.FastaFile(ref_path) as fa:
        lo, hi = min(table_pos), max(table_pos) + 1
        ref_seq = fa.fetch(item.contig, lo, hi)

    def ref_at(p: int) -> str:
        return ref_seq[p - lo]

    per_pos: dict[int, list[_Obs]] = defaultdict(list)
    per_read: list[dict[int, _Obs]] = []
    # r4 linkage: (observations, events in the window, aligned start, aligned end) per read
    linked: list[tuple[dict[int, _Obs], set[_Event], int, int]] = []
    with pysam.AlignmentFile(bam_path, "rb") as bam:
        for read in bam.fetch(item.contig, lo, hi):
            if read.is_unmapped or read.is_secondary or read.is_supplementary:
                continue
            seen = _read_observations(read, table_pos)
            for p, o in seen.items():
                per_pos[p].append(o)
            per_read.append(seen)
            if r4:
                evs = _read_events(read, lo, hi, ref_at)
                linked.append((seen, evs, read.reference_start, read.reference_end or 0))

    lines = ["# Pileup statistics", ""]
    lines.append("All primary alignments are counted (no mapping- or base-quality filter).")
    lines.append(
        "Allele notation: a base; `-` = the base is deleted in the read; "
        "`X+SEQ` = base X followed by inserted bases SEQ."
    )
    if r4:
        lines.append(
            "Mate fields come from each read's FLAG and `MQ` (mate mapping quality) tag; a "
            "read with no pairing flag is counted as unpaired."
        )
    lines.append("")
    lines.append("## Reference context")
    for s in sites:
        left = "".join(ref_at(p) for p in range(s - context, s))
        right = "".join(ref_at(p) for p in range(s + 1, s + context + 1))
        lines.append(
            f"- {item.contig}:{s - context + 1}-{s + context + 1}: "
            f"`{left}[{ref_at(s)}]{right}` (brackets = {item.contig}:{s + 1})"
        )
    lines.append("")
    lines.append("## Per-position allele counts")
    lines.append("")
    lines.append("| position | ref | depth | A | C | G | T | - | with insertion after |")
    lines.append("|---|---|---|---|---|---|---|---|---|")
    for p in sorted(table_pos):
        obs = per_pos.get(p, [])
        c = Counter(o.allele[0] if o.allele != "-" else "-" for o in obs)
        n_ins = sum(1 for o in obs if "+" in o.allele)
        lines.append(
            f"| {p + 1} | {ref_at(p)} | {len(obs)} | {c['A']} | {c['C']} | {c['G']} | "
            f"{c['T']} | {c['-']} | {n_ins} |"
        )

    detail: list[int] = []
    for p in sorted(table_pos):
        near = any(abs(p - s) <= 2 for s in sites)
        alleles = Counter(o.allele for o in per_pos.get(p, []))
        nonref = sum(n for a, n in alleles.items() if a != ref_at(p))
        if p in sites or (near and nonref >= 2):
            detail.append(p)

    for p in detail:
        obs = per_pos.get(p, [])
        lines.append("")
        lines.append(f"## Alleles at {item.contig}:{p + 1} (reference {ref_at(p)})")
        lines.append("")
        groups: dict[str, list[_Obs]] = defaultdict(list)
        for o in obs:
            groups[o.allele].append(o)
        order = sorted(groups, key=lambda a: (-len(groups[a]), a))
        lines.append("| allele | reads | fraction | forward | reverse | distinct alignments |")
        lines.append("|---|---|---|---|---|---|")
        for a in order:
            g = groups[a]
            fwd = sum(1 for o in g if not o.reverse)
            frac = len(g) / len(obs) if obs else 0.0
            lines.append(
                f"| {a} | {len(g)} | {frac:.3f} | {fwd} | {len(g) - fwd} | "
                f"{len({o.span for o in g})} |"
            )
        lines.append("")
        lines.append("Distinct alignments = distinct (start, end, orientation) combinations.")
        for a in order:
            g = groups[a]
            lines.append("")
            lines.append(f"### allele {a} ({len(g)} reads)")
            lines.append(
                f"- mapping quality: {_hist([_bin(o.mapq, _MAPQ_BINS) for o in g], _MAPQ_BINS)}"
            )
            bqs = [_bin(o.bq, _BQ_BINS) for o in g if o.bq is not None]
            if bqs:
                lines.append(f"- base quality at this position: {_hist(bqs, _BQ_BINS)}")
            dists = [_bin(o.dist3, _DIST_BINS) for o in g if o.dist3 is not None]
            if dists:
                lines.append(
                    "- distance (bp) from this position to the 3' end of the aligned part "
                    f"of the read: {_hist(dists, _DIST_BINS)}"
                )
            if r4:
                d5 = [_bin(o.dist5, _DIST_BINS) for o in g if o.dist5 is not None]
                if d5:
                    lines.append(
                        "- distance (bp) from this position to the 5' end of the aligned part "
                        f"of the read: {_hist(d5, _DIST_BINS)}"
                    )
            lines.append(
                f"- reads with soft-clipped bases at their 3' end: {sum(o.clipped3 for o in g)}"
            )
            if r4:
                lines.append(
                    f"- reads with soft-clipped bases at their 5' end: {sum(o.clipped5 for o in g)}"
                )
                lines.extend(_mate_lines(g))

    if r4:
        for s in sites:
            lines.extend(_linkage_lines(item.contig, s, ref_at(s), linked))

    if len(sites) > 1:
        combos: Counter[tuple[str, ...]] = Counter()
        for obs_map in per_read:
            if all(s in obs_map for s in sites):
                combos[tuple(obs_map[s].allele for s in sites)] += 1
        names = ", ".join(f"{item.contig}:{s + 1}" for s in sites)
        lines.append("")
        lines.append(f"## Allele combinations in reads covering all of {names}")
        lines.append("")
        lines.append("| alleles (in position order) | reads |")
        lines.append("|---|---|")
        for combo, n in sorted(combos.items(), key=lambda kv: (-kv[1], kv[0])):
            lines.append(f"| {' / '.join(combo)} | {n} |")
    return "\n".join(lines) + "\n"


# ---------------------------------------------------------------------------
# Tjson / V payloads
# ---------------------------------------------------------------------------


def sanitize_payload(payload: dict[str, Any]) -> dict[str, Any]:
    """Blind a ``ui/init`` payload via the renderer's single implementation
    (:func:`bamcp.eval.renderer.blind_payload`: drawable keys only, ``variants: []``,
    hashed read names) and pin the coordinate system for the Tjson reader."""
    from ..renderer import blind_payload

    out = blind_payload(payload)
    out["coordinate_system"] = "0-based-half-open"
    return out


def build_payload(bam_path: str, ref_path: str, contig: str, start0: int, end0: int) -> dict:
    from ...core.parsers import fetch_region
    from ...core.serialization import serialize_region_data

    data = fetch_region(
        bam_path,
        region_1based(contig, start0, end0),
        ref_path,
        max_reads=10_000,
        min_mapq=0,
        min_baseq=0,
        detect=False,
    )
    return sanitize_payload(serialize_region_data(data, compact=False))


# ---------------------------------------------------------------------------
# Tfull — shipped tool text
# ---------------------------------------------------------------------------


async def _region_summary_part(bam_path: str, item: Item, config: Any) -> tuple[str, list[str]]:
    """``get_region_summary`` over the ``WIDE_BP`` window at shipped defaults."""
    from ...core.tools import handle_get_region_summary

    start0, end0 = window(item, WIDE_BP)
    region = region_1based(item.contig, start0, end0)
    summary = await handle_get_region_summary({"file_path": bam_path, "region": region}, config)
    return region, [f"## get_region_summary(region={region})", "", summary["content"][0]["text"]]


async def _curation_part(bam_path: str, args: dict[str, Any], config: Any) -> list[str]:
    from ...analysis.curation import handle_get_variant_curation_summary

    res = await handle_get_variant_curation_summary({"file_path": bam_path, **args}, config)
    return [
        f"## get_variant_curation_summary(chrom={args['chrom']}, pos={args['pos']}, "
        f"ref={args['ref']}, alt={args['alt']})",
        "",
        res["content"][0]["text"],
    ]


async def tfull_text(bam_path: str, ref_path: str, item: Item) -> str:
    from ...config import BAMCPConfig
    from ...core.tools import _fetch_region_with_timeout

    config = BAMCPConfig(reference=ref_path)
    region, parts = await _region_summary_part(bam_path, item, config)

    data = await _fetch_region_with_timeout(bam_path, region, ref_path, config)
    for s in item.sites:
        near = [v for v in data.variants if abs(int(v["position"]) - s) <= 2]
        parts.append("")
        if not near:
            parts.append("## get_variant_curation_summary")
            parts.append("")
            parts.append(
                f"(not called: no candidate variant was reported within 2 bp of "
                f"{item.contig}:{s + 1})"
            )
            continue
        v = max(near, key=lambda v: (v.get("alt_count", 0), v.get("vaf", 0)))
        args = {
            "chrom": item.contig,
            "pos": int(v["position"]) + 1,
            "ref": v["ref"],
            "alt": v["alt"],
        }
        parts.extend(await _curation_part(bam_path, args, config))
    return "\n".join(parts) + "\n"


# ---------------------------------------------------------------------------
# Tcur — shipped tool text, curation always at the POI (J5)
# ---------------------------------------------------------------------------


def top_nonref_allele(bam_path: str, ref_path: str, contig: str, site0: int) -> dict | None:
    """The most frequent non-reference allele at ``site0`` (0-based), VCF-style.

    Counted over the same reads as ``traw.md`` (every primary alignment, no quality filter).
    Alleles are distinct events: an SNV base, an insertion after the site (whatever the
    read's base there), or one deletion (start, length) spanning the site. Ties go to an
    event anchored at the site, then to the smaller key. Returns ``{"chrom", "pos"`` (1-based,
    the VCF anchor), ``"ref", "alt", "reads"}``, or None when every read shows the reference.
    """
    import pysam

    counts: Counter[tuple[str, int, str]] = Counter()
    with pysam.FastaFile(ref_path) as fa, pysam.AlignmentFile(bam_path, "rb") as bam:
        site_ref = fa.fetch(contig, site0, site0 + 1).upper()
        for read in bam.fetch(contig, site0, site0 + 1):
            if read.is_unmapped or read.is_secondary or read.is_supplementary:
                continue
            obs = _read_observations(read, {site0}).get(site0)
            if obs is None:
                continue
            allele = obs.allele.upper()
            if allele == "-":
                pos = read.reference_start
                for op, n in read.cigartuples or []:
                    if op == 2 and pos <= site0 < pos + n:
                        counts[("del", pos, str(n))] += 1
                        break
                    if op in (0, 2, 3, 7, 8):
                        pos += n
            elif "+" in allele:
                counts[("ins", site0, allele.split("+", 1)[1])] += 1
            elif allele in ("A", "C", "G", "T") and allele != site_ref:
                counts[("snv", site0, allele)] += 1
        if not counts:
            return None
        (kind, pos, detail), n = min(
            counts.items(), key=lambda kv: (-kv[1], kv[0][1] != site0, kv[0])
        )
        if kind == "snv":
            out = {"pos": site0 + 1, "ref": site_ref, "alt": detail}
        elif kind == "ins":
            out = {"pos": site0 + 1, "ref": site_ref, "alt": site_ref + detail}
        else:  # VCF anchor = the base before the deletion
            ref = fa.fetch(contig, pos - 1, pos + int(detail)).upper()
            out = {"pos": pos, "ref": ref, "alt": ref[0]}
    return {"chrom": contig, **out, "reads": n}


async def tcur_text(bam_path: str, ref_path: str, item: Item) -> str:
    """``Tcur`` (J5): ``get_region_summary`` over the ``WIDE_BP`` window, then
    ``get_variant_curation_summary`` at every POI's most frequent non-reference allele.

    Both at shipped defaults. The curation handler itself falls back to ``probe_site``
    (min_baseq=0, no VAF / alt-read floor) when the allele fails the calling filters, and its
    text reports ``filter_check`` and ``applied_filters``; so, unlike ``Tfull``, every site
    gets a curation block even when nothing was called there.
    """
    from ...config import BAMCPConfig

    config = BAMCPConfig(reference=ref_path)
    _, parts = await _region_summary_part(bam_path, item, config)
    for s in item.sites:
        parts.append("")
        top = top_nonref_allele(bam_path, ref_path, item.contig, s)
        if top is None:
            parts.append(f"## get_variant_curation_summary(chrom={item.contig}, pos={s + 1})")
            parts.append("")
            parts.append(
                f"(not run: every read covering {item.contig}:{s + 1} shows the reference "
                "base, so there is no non-reference allele to summarise)"
            )
            continue
        args = {k: top[k] for k in ("chrom", "pos", "ref", "alt")}
        parts.extend(await _curation_part(bam_path, args, config))
    return "\n".join(parts) + "\n"


# ---------------------------------------------------------------------------
# V — blinded PNGs
# ---------------------------------------------------------------------------


async def render_pngs(jobs: list[tuple[dict[str, Any], Path]]) -> None:
    """Render each (payload, out_path) blinded (``capture(..., blind=True)``: pinned display
    settings, ``#viewer`` clip, overlays hidden) in expanded mode, with the viewport
    height auto-fitted to the read depth so no rows are cut."""
    from ..renderer import ViewerRenderer

    async with ViewerRenderer() as renderer:
        for payload, out_path in jobs:
            png = await renderer.capture(
                payload, "expanded", blind=True, viewport=(VIEW_WIDTH, None)
            )
            out_path.write_bytes(png)


def _viewer_hash() -> str:
    from ...viewer_build import viewer_source_hash

    return viewer_source_hash()


def viewer_bundle_problem() -> str | None:
    """Why ``dist/viewer.html`` must not be used for stimuli, or None when it is current.

    PNG stimuli are only meaningful for the viewer they claim to show: a bundle built from
    older sources renders the old encoding while ``VIEWER_KEY`` describes the new one."""
    from ..render_cli import _check_rebuild_needed

    return _check_rebuild_needed()


async def chromium_available() -> str | None:
    """None when a headless Chromium can be launched, else the reason it cannot."""
    bundle = Path(__file__).resolve().parents[2] / "static" / "dist" / "viewer.html"
    if not bundle.exists():
        return "viewer bundle missing — run `make build-viewer`"
    try:
        from ..renderer import ViewerRenderer

        async with ViewerRenderer():
            return None
    except Exception as e:  # noqa: BLE001 — any launch failure means "no PNGs"
        return f"{type(e).__name__}: {str(e).splitlines()[0] if str(e) else ''}"


# ---------------------------------------------------------------------------
# Text pieces
# ---------------------------------------------------------------------------


SYNTHETIC_SAMPLE_LINE = "simulated diploid sample on a synthetic GRCh38-like reference"
# Round-4 caption sample line (design section 4; frozen before the pilot).
SAMPLE_LINE_R4 = (
    "sequenced human sample (Illumina short reads, ~60x whole-genome); the sample may carry "
    "germline variants or variants present in only a fraction of cells (mosaic); "
    "coordinates are re-indexed"
)
ANCHOR_SENTENCE = (
    "An insertion or deletion may be anchored at, or within 2 bp of, the position of "
    "interest (e.g. left-aligned in a repeat)."
)
# Round 4 (design section 4): replaces ANCHOR_SENTENCE when ``r4=True``.
ANCHOR_SENTENCE_R4 = (
    "An insertion or deletion may be anchored anywhere within the repeat run that contains "
    "the position of interest, or within 2 bp of it."
)


def caption_text(item: Item, sample_line: str = SYNTHETIC_SAMPLE_LINE) -> str:
    start0, end0 = window(item, WIDE_BP)
    poi = ", ".join(f"{item.contig}:{s + 1}" for s in item.sites)
    return (
        "# Site\n\n"
        f"- Region: {region_1based(item.contig, start0, end0)} (1-based, inclusive)\n"
        f"- Position(s) of interest: {poi}\n"
        f"- Sample: {sample_line}\n"
    )


def prompt_text(files: list[str], sample_kind: str = "simulated", *, r4: bool = False) -> str:
    """The uniform task prompt for a file set. ``r4`` switches to the round-4 ``traw.md``
    description and indel anchoring sentence; everything else is the round-3 text."""
    descriptions = FILE_DESCRIPTIONS_R4 if r4 else FILE_DESCRIPTIONS
    anchor = ANCHOR_SENTENCE_R4 if r4 else ANCHOR_SENTENCE
    listed = "\n".join(f"- `{f}` — {descriptions[f]}" for f in files)
    key = f"\n\n{VIEWER_KEY}" if any(f.endswith(".png") for f in files) else ""
    return (
        "# Task\n\n"
        f"You are reviewing sequencing evidence for one site in a {sample_kind} sample. Use ONLY "
        "the files in this directory listed below; do not open any other file. Decide "
        "whether the sample truly carries a variant at the position(s) of interest and "
        "whether any non-reference signal there is a technical artifact.\n\n"
        "## Files\n\n"
        f"- `caption.md` — the region and position(s) of interest.\n{listed}{key}\n\n"
        "## Answer\n\n"
        "Reply with ONLY one JSON object (no prose, no code fence) matching:\n\n"
        f"```\n{ANSWER_SCHEMA_TEXT}\n```\n\n"
        f"{FIELD_GUIDE}\n\n"
        f"Coordinates are 1-based. {anchor}\n"
    )


MPILEUP_ARGS = ("-B", "-Q", "0", "-q", "0", "-s", "--output-BP")
MPILEUP_HEADER = """\
# samtools mpileup -f <reference> -r {region} -B -Q 0 -q 0 -s --output-BP <bam>
#
# One line per reference position: contig, position (1-based), reference base, read depth,
# read bases, base qualities, mapping qualities, and each read's position within the read.
# Read bases: `.` / `,` = matches the reference on the forward / reverse strand; ACGTN /
# acgtn = a mismatch on the forward / reverse strand; `+3ACG` / `-2tt` = an insertion /
# deletion of those bases starting after this position; `*` = a base deleted in this read;
# `^` = start of a read (the next character is its mapping quality + 33); `$` = end of a
# read. Base and mapping qualities are Phred + 33 ASCII, one character per read, in the
# same order as the read bases. Default read filters: unmapped, secondary, QC-fail and
# duplicate reads are skipped.
"""


def mpileup_text(bam_path: str, ref_path: str, item: Item, half: int = MPILEUP_HALF) -> str:
    """samtools mpileup over the position(s) of interest +/- ``half`` bp (round 4's
    raw-input baseline: what a user without BAMCP would paste)."""
    import pysam

    lo = max(0, min(item.sites) - half)
    hi = max(item.sites) + half + 1
    region = region_1based(item.contig, lo, hi)
    body = pysam.mpileup("-f", ref_path, "-r", region, *MPILEUP_ARGS, bam_path)
    text = body if isinstance(body, str) else "".join(body)
    return MPILEUP_HEADER.format(region=region) + "\n" + text


def condition_files(condition: str) -> list[str]:
    pngs = ["view_wide.png", "view_zoom.png"]
    files = {
        "Traw": ["traw.md"],
        "Tjson": ["tjson.json"],
        "Tfull": ["tfull.md"],
        "V": pngs,
        "V+Traw": [*pngs, "traw.md"],
        "Tfull+V": ["tfull.md", *pngs],
        "Traw-r4": ["traw.md"],
        "Traw-r4+Tcur": ["traw.md", "tcur.md"],
        "V+Traw-r4": ["traw.md", *pngs],
        "Tmp": ["mpileup.txt"],
        "Tmp+Tcur": ["mpileup.txt", "tcur.md"],
        "V+Tmp": ["mpileup.txt", *pngs],
    }[condition]
    return sorted(files, key=FILE_ORDER.index)


def stim_id(seed: int, item_id: str, condition: str) -> str:
    return opaque_id(seed, "stim", item_id, condition)


def run_id(stim: str, model: str, repeat: int) -> str:
    return f"{stim}-{opaque_id(stim, model, repeat, n=6)}"


def select_items(plan: Plan, n_items: int | None) -> list[Item]:
    """Stratified round-robin over case types (instance 0 of every type first, ...)."""
    if n_items is None or n_items >= len(plan.items):
        return list(plan.items)
    by_inst = sorted(plan.items, key=lambda it: (it.instance, it.case_type))
    return by_inst[:n_items]


def label_record(item: Item) -> dict[str, Any]:
    spec = CASE_BY_NAME[item.case_type]
    return {
        "item": item.item_id,
        "case_type": item.case_type,
        "instance": item.instance,
        "contig": item.contig,
        "positions": item.positions_1based(),
        **item.label,
        "lenient_verdicts": [item.label["artifact_verdict"], *spec.lenient_verdicts],
    }


# ---------------------------------------------------------------------------
# Build
# ---------------------------------------------------------------------------


def build_run_set(
    out: str | Path,
    *,
    seed: int = 1,
    conditions: tuple[str, ...] | list[str] = CONDITIONS,
    models: tuple[str, ...] | list[str] = ("model",),
    repeats: int = 1,
    n_items: int | None = None,
    instances: int = 3,
    png: bool | None = None,
    fixture_dir: str | Path | None = None,
    log: Any = print,
) -> BuildResult:
    """Build fixtures (if needed), stimuli, labels and the runs manifest.

    ``png``: True = PNGs required (fail without Chromium); False = skip visual
    conditions; None = render when Chromium is available, else skip with a message.
    """
    unknown = [c for c in conditions if c not in CONDITIONS]
    if unknown:
        raise ValueError(f"unknown condition(s): {unknown}; choose from {CONDITIONS}")
    paths = ensure_fixtures(seed, fixture_dir, instances)
    plan = build_plan(seed, instances)
    return build_items_run_set(
        out,
        bam=str(paths["bam"]),
        ref=str(paths["reference"]),
        items=select_items(plan, n_items),
        label_records=[label_record(it) for it in plan.items],
        seed=seed,
        conditions=conditions,
        models=models,
        repeats=repeats,
        png=png,
        log=log,
        meta={"tier": "synthetic", "instances": instances},
    )


def build_items_run_set(
    out: str | Path,
    *,
    bam: str,
    ref: str,
    items: list[Item],
    label_records: list[dict[str, Any]],
    seed: int,
    conditions: tuple[str, ...] | list[str] = CONDITIONS,
    models: tuple[str, ...] | list[str] = ("model",),
    repeats: int = 1,
    png: bool | None = None,
    log: Any = print,
    meta: dict[str, Any] | None = None,
    sample_line: str = SYNTHETIC_SAMPLE_LINE,
    sample_kind: str = "simulated",
    r4: bool = False,
) -> BuildResult:
    """Stimuli, labels and the runs manifest for ``items`` read from ``bam``/``ref``.

    Shared by the synthetic tier (:func:`build_run_set`) and the real-data tier
    (:mod:`.real`); ``label_records`` is written verbatim to ``labels/labels.jsonl``.
    ``r4`` builds round-4 conditions only (:data:`CONDITIONS_R4`): ``traw.md`` from
    :func:`traw_r4_text`, ``tcur.md`` from :func:`tcur_text`, and the round-4 prompt; it is
    recorded in ``build.json`` so the leak-lint can rebuild the uniform prompt.
    """
    allowed = CONDITIONS_R4 if r4 else CONDITIONS
    unknown = [c for c in conditions if c not in allowed]
    if unknown:
        raise ValueError(f"unknown condition(s): {unknown}; choose from {allowed}")
    out = Path(out)
    skipped: list[str] = []
    want_png = any(c in VISUAL_CONDITIONS for c in conditions)
    if want_png and png is not True:
        reason = "disabled (--no-png)" if png is False else asyncio.run(chromium_available())
        if reason:
            skipped = [c for c in conditions if c in VISUAL_CONDITIONS]
            log(
                f"SKIP PNG step — Chromium unavailable ({reason}); not building "
                f"{', '.join(skipped)}. Set PLAYWRIGHT_CHROMIUM_EXECUTABLE to render them."
            )
            want_png = False
    if want_png:
        problem = viewer_bundle_problem()
        if problem:
            raise RuntimeError(f"refusing to render stimuli from this viewer bundle: {problem}")
    active = [c for c in conditions if c not in skipped]

    for sub in ("stim", "labels", "answers"):
        shutil.rmtree(out / sub, ignore_errors=True)
    (out / "stim").mkdir(parents=True)
    (out / "labels").mkdir(parents=True)
    (out / "answers").mkdir(parents=True)

    needs = {f for c in active for f in condition_files(c)}
    evidence: dict[str, dict[str, Any]] = {it.item_id: {} for it in items}
    png_jobs: list[tuple[dict[str, Any], Path]] = []
    tmp_png = out / "labels" / "_png"  # staging area outside stim/
    for it in items:
        ev = evidence[it.item_id]
        if "traw.md" in needs:
            ev["traw.md"] = (traw_r4_text if r4 else traw_text)(bam, ref, it)
        if "mpileup.txt" in needs:
            ev["mpileup.txt"] = mpileup_text(bam, ref, it)
        if "tjson.json" in needs or want_png:
            s0, e0 = window(it, WIDE_BP)
            wide = build_payload(bam, ref, it.contig, s0, e0)
            ev["tjson.json"] = json.dumps(wide, separators=(",", ":"))
            if want_png:
                z0, z1 = window(it, ZOOM_BP)
                zoom = build_payload(bam, ref, it.contig, z0, z1)
                tmp_png.mkdir(parents=True, exist_ok=True)
                png_jobs.append((wide, tmp_png / f"{it.item_id}_wide.png"))
                png_jobs.append((zoom, tmp_png / f"{it.item_id}_zoom.png"))

    if "tfull.md" in needs:

        async def _all_tfull() -> list[str]:
            return [await tfull_text(bam, ref, it) for it in items]

        for it, text in zip(items, asyncio.run(_all_tfull()), strict=True):
            evidence[it.item_id]["tfull.md"] = text

    if "tcur.md" in needs:

        async def _all_tcur() -> list[str]:
            return [await tcur_text(bam, ref, it) for it in items]

        for it, text in zip(items, asyncio.run(_all_tcur()), strict=True):
            evidence[it.item_id]["tcur.md"] = text

    if png_jobs:
        log(f"rendering {len(png_jobs)} PNGs ...")
        asyncio.run(render_pngs(png_jobs))
        for it in items:
            evidence[it.item_id]["view_wide.png"] = (
                tmp_png / f"{it.item_id}_wide.png"
            ).read_bytes()
            evidence[it.item_id]["view_zoom.png"] = (
                tmp_png / f"{it.item_id}_zoom.png"
            ).read_bytes()
        shutil.rmtree(tmp_png, ignore_errors=True)

    with (out / "labels" / "labels.jsonl").open("w") as f:
        for rec in label_records:
            f.write(json.dumps(rec, sort_keys=True) + "\n")

    n_stim = 0
    runs: list[dict[str, Any]] = []
    for it in items:
        cap = caption_text(it, sample_line)
        for cond in active:
            sid = stim_id(seed, it.item_id, cond)
            d = out / "stim" / sid
            d.mkdir()
            files = condition_files(cond)
            (d / "caption.md").write_text(cap)
            (d / "prompt.md").write_text(prompt_text(files, sample_kind, r4=r4))
            for name in files:
                content = evidence[it.item_id][name]
                if isinstance(content, bytes):
                    (d / name).write_bytes(content)
                else:
                    (d / name).write_text(content)
            n_stim += 1
            rel = [f"stim/{sid}/{n}" for n in ("prompt.md", "caption.md", *files)]
            for model in models:
                for rep in range(repeats):
                    rid = run_id(sid, model, rep)
                    runs.append(
                        {
                            "run_id": rid,
                            "item": it.item_id,
                            "condition": cond,
                            "model": model,
                            "repeat": rep,
                            "stim_dir": f"stim/{sid}",
                            "files": rel,
                            "answer_path": f"answers/{rid}.json",
                        }
                    )
    # Shuffle-free but item-interleaved order: sorted by run_id (opaque) so the manifest
    # order carries no case-type information.
    runs.sort(key=lambda r: r["run_id"])
    with (out / "runs.jsonl").open("w") as f:
        for r in runs:
            f.write(json.dumps(r, sort_keys=True) + "\n")
    build_meta = {
        **(meta or {}),
        "seed": seed,
        "conditions": active,
        "skipped_conditions": skipped,
        "models": list(models),
        "repeats": repeats,
        "n_items": len(items),
        "pngs": bool(png_jobs),
        "viewer_sources": _viewer_hash() if png_jobs else None,
        "bam": bam,
        "reference": ref,
    }
    if r4:
        build_meta["r4"] = True
    (out / "build.json").write_text(json.dumps(build_meta, indent=2) + "\n")
    return BuildResult(out, len(items), n_stim, len(runs), bool(png_jobs), skipped)
