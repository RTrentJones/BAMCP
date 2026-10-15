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
VISUAL_CONDITIONS = frozenset({"V", "V+Traw", "Tfull+V"})

WIDE_BP = 200
ZOOM_BP = 60
VIEW_WIDTH = 1024

FILE_ORDER = ("traw.md", "tjson.json", "tfull.md", "view_wide.png", "view_zoom.png")
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
    "view_wide.png": "Screenshot of the read pileup viewer, ~200 bp centred on the site.",
    "view_zoom.png": "Screenshot of the read pileup viewer, ~60 bp centred on the site.",
}

VIEWER_KEY = """\
Viewer key (screenshots; describes the BAMCP viewer as rendered at build time): top to
bottom — coordinate ruler (1-based), reference sequence (coloured base tiles), per-base
coverage track (the dashed line marks 10x), then one row per read, reads packed into rows.
Reads aligned forward are blue, reads aligned in reverse are red; the pointed end shows the
alignment direction. Reads with lower mapping quality are drawn paler; reads with mapping
quality 0 are drawn as outlined, nearly transparent boxes. At the wide zoom, read bases that
differ from the reference are drawn as coloured ticks (A green, C blue, G orange, T red)
and matching bases are not marked. At base-level zoom every read base is drawn as a coloured
letter tile (same colours) with a blue/red bar under the read for its orientation, so
compare the rows with the reference row to spot differences. A deletion is a gap in the read
bridged by a thin black line; an insertion is a purple vertical mark between two bases.
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


@dataclass
class _Obs:
    allele: str
    reverse: bool
    mapq: int
    bq: int | None
    dist3: int | None
    clipped3: bool
    span: tuple[int, int, bool]


def _read_observations(read: Any, positions: set[int]) -> dict[int, _Obs]:
    """Allele of one read at each of ``positions`` (0-based)."""
    seq = read.query_sequence or ""
    quals = read.query_qualities
    rev = bool(read.is_reverse)
    q_first, q_last = read.query_alignment_start, read.query_alignment_end - 1
    cig = read.cigartuples or []
    clipped3 = bool(cig) and ((cig[0][0] == 4) if rev else (cig[-1][0] == 4))
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
            out[rpos] = _Obs("-", rev, read.mapping_quality, None, None, clipped3, span)
        else:
            dist3 = (qpos - q_first) if rev else (q_last - qpos)
            bq = int(quals[qpos]) if quals is not None else None
            out[rpos] = _Obs(seq[qpos], rev, read.mapping_quality, bq, dist3, clipped3, span)
    for rpos, ins in pending_ins.items():
        if rpos in out and out[rpos].allele != "-":
            out[rpos].allele = f"{out[rpos].allele}+{ins}"
    return out


def _hist(values: list[str], bins: tuple[tuple[int, int, str], ...]) -> str:
    c = Counter(values)
    return " | ".join(f"{name}: {c.get(name, 0)}" for _, _, name in bins)


def traw_text(bam_path: str, ref_path: str, item: Item) -> str:
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
    with pysam.AlignmentFile(bam_path, "rb") as bam:
        for read in bam.fetch(item.contig, lo, hi):
            if read.is_unmapped or read.is_secondary or read.is_supplementary:
                continue
            seen = _read_observations(read, table_pos)
            for p, o in seen.items():
                per_pos[p].append(o)
            per_read.append(seen)

    lines = ["# Pileup statistics", ""]
    lines.append("All primary alignments are counted (no mapping- or base-quality filter).")
    lines.append(
        "Allele notation: a base; `-` = the base is deleted in the read; "
        "`X+SEQ` = base X followed by inserted bases SEQ."
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
            lines.append(
                f"- reads with soft-clipped bases at their 3' end: {sum(o.clipped3 for o in g)}"
            )

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


async def tfull_text(bam_path: str, ref_path: str, item: Item) -> str:
    from ...analysis.curation import handle_get_variant_curation_summary
    from ...config import BAMCPConfig
    from ...core.tools import _fetch_region_with_timeout, handle_get_region_summary

    config = BAMCPConfig(reference=ref_path)
    start0, end0 = window(item, WIDE_BP)
    region = region_1based(item.contig, start0, end0)
    summary = await handle_get_region_summary({"file_path": bam_path, "region": region}, config)
    parts = [f"## get_region_summary(region={region})", "", summary["content"][0]["text"]]

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
            "file_path": bam_path,
            "chrom": item.contig,
            "pos": int(v["position"]) + 1,
            "ref": v["ref"],
            "alt": v["alt"],
        }
        res = await handle_get_variant_curation_summary(args, config)
        parts.append(
            f"## get_variant_curation_summary(chrom={item.contig}, pos={args['pos']}, "
            f"ref={args['ref']}, alt={args['alt']})"
        )
        parts.append("")
        parts.append(res["content"][0]["text"])
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


def caption_text(item: Item) -> str:
    start0, end0 = window(item, WIDE_BP)
    poi = ", ".join(f"{item.contig}:{s + 1}" for s in item.sites)
    return (
        "# Site\n\n"
        f"- Region: {region_1based(item.contig, start0, end0)} (1-based, inclusive)\n"
        f"- Position(s) of interest: {poi}\n"
        "- Sample: simulated diploid sample on a synthetic GRCh38-like reference\n"
    )


def prompt_text(files: list[str]) -> str:
    listed = "\n".join(f"- `{f}` — {FILE_DESCRIPTIONS[f]}" for f in files)
    key = f"\n\n{VIEWER_KEY}" if any(f.endswith(".png") for f in files) else ""
    return (
        "# Task\n\n"
        "You are reviewing sequencing evidence for one site in a simulated sample. Use ONLY "
        "the files in this directory listed below; do not open any other file. Decide "
        "whether the sample truly carries a variant at the position(s) of interest and "
        "whether any non-reference signal there is a technical artifact.\n\n"
        "## Files\n\n"
        f"- `caption.md` — the region and position(s) of interest.\n{listed}{key}\n\n"
        "## Answer\n\n"
        "Reply with ONLY one JSON object (no prose, no code fence) matching:\n\n"
        f"```\n{ANSWER_SCHEMA_TEXT}\n```\n\n"
        f"{FIELD_GUIDE}\n\n"
        "Coordinates are 1-based. An insertion or deletion may be anchored at, or within "
        "2 bp of, the position of interest (e.g. left-aligned in a repeat).\n"
    )


def condition_files(condition: str) -> list[str]:
    pngs = ["view_wide.png", "view_zoom.png"]
    files = {
        "Traw": ["traw.md"],
        "Tjson": ["tjson.json"],
        "Tfull": ["tfull.md"],
        "V": pngs,
        "V+Traw": [*pngs, "traw.md"],
        "Tfull+V": ["tfull.md", *pngs],
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
    out = Path(out)
    paths = ensure_fixtures(seed, fixture_dir, instances)
    bam, ref = str(paths["bam"]), str(paths["reference"])
    plan = build_plan(seed, instances)
    items = select_items(plan, n_items)

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
            ev["traw.md"] = traw_text(bam, ref, it)
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
        for it in plan.items:
            f.write(json.dumps(label_record(it), sort_keys=True) + "\n")

    n_stim = 0
    runs: list[dict[str, Any]] = []
    for it in items:
        cap = caption_text(it)
        for cond in active:
            sid = stim_id(seed, it.item_id, cond)
            d = out / "stim" / sid
            d.mkdir()
            files = condition_files(cond)
            (d / "caption.md").write_text(cap)
            (d / "prompt.md").write_text(prompt_text(files))
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
    meta = {
        "seed": seed,
        "instances": instances,
        "conditions": active,
        "skipped_conditions": skipped,
        "models": list(models),
        "repeats": repeats,
        "n_items": len(items),
        "pngs": bool(png_jobs),
        "bam": bam,
        "reference": ref,
    }
    (out / "build.json").write_text(json.dumps(meta, indent=2) + "\n")
    return BuildResult(out, len(items), n_stim, len(runs), bool(png_jobs), skipped)
