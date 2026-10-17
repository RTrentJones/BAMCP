"""Leak-lint: fail a run set whose stimuli could reveal labels.

Checks, per pre-registration:

1. ``stim/`` children are opaque 8-hex directory names; each holds only the allowed files.
2. Evidence files (everything except ``prompt.md``, whose shared task/schema text names the
   answer vocabulary on purpose) contain no label vocabulary. ``tfull.md`` is exempt from the
   vocabulary rule — it is the shipped tool's own wording, which is what Tfull measures —
   but not from the other checks.
3. No stimulus file (PNG text chunks included) contains a read name from the BAM, an item
   id, or a case-type name.
4. Every caption matches the single uniform template.
5. Nothing label-bearing (``labels.jsonl``) sits under ``stim/``.

Round 4 (``docs/eval/visual-design-r4.md`` section 10) adds:

6. No evidence file names a truth set or stratification (GIAB, CMRG, v4.2.1, v5.0q,
   benchmark, segdup, HiFi, PacBio, trio); ``mosaic`` appears only in ``caption.md`` and
   ``prompt.md`` (the shipped tool text's own fixed phrases are stripped first, see
   :data:`PRODUCT_PHRASES`). ``tcur.md`` is, like ``tfull.md``, exempt from rule 2 only.
7. The fixtures' "mate elsewhere" sentinel contig carries no real coordinate.

Selection-side helpers (run on the manifest, not on stimuli): :func:`ks_two_sample` /
:func:`lint_poi_offsets` (strata 3 vs 4 POI offset within the repeat run) and
:func:`lint_anchor_collisions` (strata 7-8: the anchoring sentence never maps the POI onto
another truth record).
"""

from __future__ import annotations

import json
import math
import re
import zlib
from collections.abc import Iterable, Sequence
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

from .cases import CASE_TYPES
from .real import REAL_CASE_TYPES
from .real import SAMPLE_LINE as REAL_SAMPLE_LINE
from .stimuli import FILE_ORDER, SAMPLE_LINE_R4, SYNTHETIC_SAMPLE_LINE, prompt_text

# Round-4 strata. ``strata_r4.STRATA_R4`` is owned by the selection code; this local copy keeps
# the lint usable without it, and a test asserts the two match whenever both exist.
R4_CASE_TYPES = (
    "mosaic_snv_pos",
    "mosaic_caller_fp_neg",
    "tr_indel_pos",
    "tr_stutter_neg",
    "segdup_cmrg_pos",
    "segdup_psv_neg",
    "indel_adjacent_misalign_neg",
    "indel_adjacent_snv_pos",
    "anchor_het_snv",
    "anchor_hom_snv",
    "anchor_clean_neg",
)
SENTINEL_CONTIG = "chrU"
SENTINEL_MAX_LEN = 1000


def r4_case_types() -> tuple[str, ...]:
    """``R4_CASE_TYPES`` plus anything ``strata_r4.STRATA_R4`` adds (imported lazily)."""
    try:
        from .strata_r4 import STRATA_R4  # type: ignore[import-not-found,unused-ignore]
    except ImportError:
        return R4_CASE_TYPES
    return tuple(dict.fromkeys([*R4_CASE_TYPES, *STRATA_R4]))


ALLOWED_FILES = frozenset(
    {
        "caption.md",
        "prompt.md",
        "traw.md",
        "tjson.json",
        "tfull.md",
        "tcur.md",
        "view_wide.png",
        "view_zoom.png",
    }
)
VOCAB_EXEMPT = frozenset({"prompt.md", "tfull.md", "tcur.md"})
_OPAQUE_DIR = re.compile(r"^[0-9a-f]{8}$")
LEAK_VOCAB = re.compile(
    r"strand|bias|homopolymer|artifact|artefact|mapq_low|low_mapq|low_baseq|low_bq|dup|"
    r"pcr|risk|imbalanced|alt_fwd|alt_rev|_ref_|paralog|verdict|concern|recommend|"
    r"confidence|\bhet\b|heterozyg|\bhom\b|homozyg|hom_alt|\bcis\b|\btrans\b|subclonal|"
    r"low_vaf|\bnoise\b|\berror|misalign|read_end|negative|scatter|\breal\b|\bmnv\b|"
    r"multi_allelic|phas",
    re.IGNORECASE,
)
# Section 10 of the round-4 design: truth-set / stratification names, in every evidence file.
TRUTH_VOCAB = re.compile(
    r"giab|cmrg|v4\.2\.1|v5\.0q|benchmark|segdup|hifi|pacbio|\btrio\b", re.IGNORECASE
)
# ``mosaic`` only in caption.md / prompt.md (the round-4 sample line names it on purpose).
MOSAIC_VOCAB = re.compile(r"\bmosaic", re.IGNORECASE)
MOSAIC_ALLOWED = frozenset({"caption.md", "prompt.md"})
# Fixed wording of the shipped tools (``analysis/evidence.zygosity_hint`` and the J6 low-VAF
# note), stripped from tool text before the mosaic check: it depends on VAF alone, never on
# the label, and Tfull/Tcur measure the tool as shipped.
PRODUCT_PHRASES = ("subclonal/mosaic", "mosaicism")
TOOL_TEXT_FILES = frozenset({"tfull.md", "tcur.md"})
CAPTION_RE = re.compile(
    r"^# Site\n\n"
    r"- Region: chr\d+:\d+-\d+ \(1-based, inclusive\)\n"
    r"- Position\(s\) of interest: chr\d+:\d+(, chr\d+:\d+)?\n"
    r"- Sample: ("
    + "|".join(re.escape(x) for x in (SYNTHETIC_SAMPLE_LINE, REAL_SAMPLE_LINE, SAMPLE_LINE_R4))
    + r")\n$"
)
_SENTINEL_COORD = re.compile(re.escape(SENTINEL_CONTIG) + r"\s*:\s*\d")


@dataclass
class LintReport:
    problems: list[str] = field(default_factory=list)
    files_checked: int = 0

    @property
    def ok(self) -> bool:
        return not self.problems


def png_text_chunks(data: bytes) -> str:
    """Concatenated tEXt/zTXt/iTXt payloads of a PNG (metadata a model could read)."""
    out: list[str] = []
    i = 8
    while i + 8 <= len(data):
        length = int.from_bytes(data[i : i + 4], "big")
        ctype = data[i + 4 : i + 8]
        body = data[i + 8 : i + 8 + length]
        if ctype == b"tEXt":
            out.append(body.decode("latin-1"))
        elif ctype == b"zTXt":
            key, _, rest = body.partition(b"\0")
            try:
                out.append(key.decode("latin-1") + " " + zlib.decompress(rest[1:]).decode())
            except zlib.error:
                out.append(key.decode("latin-1"))
        elif ctype == b"iTXt":
            out.append(body.decode("utf-8", "replace"))
        i += 12 + length
    return "\n".join(out)


def _bam_read_names(bam_path: str | None) -> set[str]:
    if not bam_path or not Path(bam_path).exists():
        return set()
    import pysam

    with pysam.AlignmentFile(bam_path, "rb") as bam:
        return {r.query_name for r in bam.fetch(until_eof=True) if r.query_name}


def lint_run_set(run_dir: str | Path, bam_path: str | None = None) -> LintReport:
    run_dir = Path(run_dir)
    report = LintReport()
    stim = run_dir / "stim"
    if not stim.is_dir():
        report.problems.append(f"{stim}: missing stim/ directory")
        return report
    if bam_path is None:
        meta = run_dir / "build.json"
        if meta.exists():
            bam_path = json.loads(meta.read_text()).get("bam")
    names = _bam_read_names(bam_path)
    labels_path = run_dir / "labels" / "labels.jsonl"
    ids: set[str] = set()
    if labels_path.exists():
        ids = {json.loads(line)["item"] for line in labels_path.read_text().splitlines() if line}
    case_names = {c.name for c in CASE_TYPES} | set(REAL_CASE_TYPES) | set(r4_case_types())
    meta_path = run_dir / "build.json"
    build = json.loads(meta_path.read_text()) if meta_path.exists() else {}
    r4 = bool(build.get("r4"))
    sample_kind = "sequenced" if build.get("tier") == "real" or r4 else "simulated"
    report.problems.extend(lint_sentinel_contig(bam_path, build.get("reference")))

    for path in sorted(stim.rglob("*")):
        rel = path.relative_to(run_dir)
        if path.is_dir():
            if path.parent != stim:
                report.problems.append(f"{rel}: nested directory under a stimulus dir")
            elif not _OPAQUE_DIR.match(path.name):
                report.problems.append(f"{rel}: stimulus directory name is not opaque (8 hex)")
            continue
        if path.parent.parent != stim:
            report.problems.append(f"{rel}: file outside a stimulus directory")
            continue
        if path.name not in ALLOWED_FILES:
            report.problems.append(f"{rel}: unexpected file name")
            continue
        report.files_checked += 1
        raw = path.read_bytes()
        text = png_text_chunks(raw) if path.suffix == ".png" else raw.decode("utf-8", "replace")
        if path.name not in VOCAB_EXEMPT:
            m = LEAK_VOCAB.search(text)
            if m:
                report.problems.append(f"{rel}: label vocabulary {m.group(0)!r}")
        report.problems.extend(f"{rel}: {p}" for p in r4_vocab_problems(path.name, text))
        if path.name == "caption.md" and not CAPTION_RE.match(text):
            report.problems.append(f"{rel}: caption does not match the uniform template")
        elif path.name == "caption.md" and r4 and f"- Sample: {SAMPLE_LINE_R4}\n" not in text:
            # Round 4 freezes one caption sample line (design section 4).
            report.problems.append(f"{rel}: round-4 caption must use SAMPLE_LINE_R4")
        if path.name == "prompt.md":
            # The shared prompt names the answer vocabulary on purpose; instead of a
            # vocabulary scan it must be exactly the uniform text for its file set.
            present = sorted(p.name for p in path.parent.iterdir() if p.name in FILE_ORDER)
            present.sort(key=FILE_ORDER.index)
            if text != prompt_text(present, sample_kind, r4=r4):
                report.problems.append(f"{rel}: prompt differs from the uniform prompt")
            continue
        for token in _tokens(text):
            if token in names:
                report.problems.append(f"{rel}: contains a BAM read name {token!r}")
                break
            if token in ids:
                report.problems.append(f"{rel}: contains an item id {token!r}")
                break
            if token in case_names:
                report.problems.append(f"{rel}: contains a case-type name {token!r}")
                break
    return report


def r4_vocab_problems(name: str, text: str) -> list[str]:
    """Section-10 vocabulary problems of one stimulus file (``name`` = its file name)."""
    if name in ("caption.md", "prompt.md"):
        # Both are checked against their uniform templates; only the truth names apply.
        m = TRUTH_VOCAB.search(text)
        return [f"truth-set vocabulary {m.group(0)!r}"] if m else []
    out: list[str] = []
    m = TRUTH_VOCAB.search(text)
    if m:
        out.append(f"truth-set vocabulary {m.group(0)!r}")
    scan = text
    if name in TOOL_TEXT_FILES:
        for phrase in PRODUCT_PHRASES:
            scan = re.sub(re.escape(phrase), " ", scan, flags=re.IGNORECASE)
    m = MOSAIC_VOCAB.search(scan)
    if m and name not in MOSAIC_ALLOWED:
        out.append(f"{m.group(0)!r} outside caption.md / prompt.md")
    if _SENTINEL_COORD.search(text):
        out.append(f"names a coordinate on the sentinel contig {SENTINEL_CONTIG}")
    return out


def lint_sentinel_contig(bam_path: str | None, ref_path: str | None = None) -> list[str]:
    """The "mate elsewhere" sentinel contig must carry no real coordinate.

    Every mapped mate's position must lie inside its (fixture) contig, else it is an
    un-re-indexed real coordinate. When the BAM declares ``SENTINEL_CONTIG``: it is at most
    ``SENTINEL_MAX_LEN`` bp, no read is aligned to it, every mate placed on it sits inside
    it at one shared position (a varying position could encode the mate's real one), and
    its reference sequence (when ``ref_path`` has it) is all ``N``.
    """
    if not bam_path or not Path(bam_path).exists():
        return []
    import pysam

    problems: list[str] = []
    with pysam.AlignmentFile(bam_path, "rb") as bam:
        has_sentinel = SENTINEL_CONTIG in bam.references
        tid = bam.get_tid(SENTINEL_CONTIG) if has_sentinel else -2
        length = bam.get_reference_length(SENTINEL_CONTIG) if has_sentinel else 0
        if length > SENTINEL_MAX_LEN:
            problems.append(f"{SENTINEL_CONTIG}: length {length} > {SENTINEL_MAX_LEN}")
        lengths = list(bam.lengths)
        placed = 0
        mate_pos: set[int] = set()
        off_contig = 0
        for r in bam.fetch(until_eof=True):
            if r.reference_id == tid:
                placed += 1
            if r.is_paired and has_sentinel and r.next_reference_id == tid:
                mate_pos.add(int(r.next_reference_start))
            elif (
                r.is_paired
                and not r.mate_is_unmapped
                and 0 <= r.next_reference_id < len(lengths)
                and not 0 <= r.next_reference_start < lengths[r.next_reference_id]
            ):
                # A mate position past its fixture contig is an un-re-indexed real coordinate.
                off_contig += 1
    if off_contig:
        problems.append(
            f"{off_contig} read(s) with a mate position outside its contig (not re-indexed)"
        )
    if not has_sentinel:
        return problems
    if placed:
        problems.append(f"{SENTINEL_CONTIG}: {placed} read(s) aligned to the sentinel contig")
    outside = sorted(p for p in mate_pos if not 0 <= p < length)
    if outside:
        problems.append(f"{SENTINEL_CONTIG}: mate position(s) outside the contig: {outside[:5]}")
    if len(mate_pos) > 1:
        problems.append(
            f"{SENTINEL_CONTIG}: {len(mate_pos)} distinct mate positions (must be one constant)"
        )
    if ref_path and Path(ref_path).exists():
        with pysam.FastaFile(ref_path) as fa:
            if SENTINEL_CONTIG in fa.references:
                seq = fa.fetch(SENTINEL_CONTIG).upper()
                if set(seq) - {"N"}:
                    problems.append(f"{SENTINEL_CONTIG}: reference sequence is not all N")
    return problems


# ---------------------------------------------------------------------------
# Selection-side lints (site manifest, before any stimulus is built)
# ---------------------------------------------------------------------------


_KS_EXACT_MAX_NM = 250_000  # exact DP is O(n*m)


def _ks_exact_p(n: int, m: int, d: float, checkpoints: set[int] | None = None) -> float:
    """P(D >= d) for the two-sample KS statistic under the permutation null, exact.

    Walks the monotone lattice paths from (0, 0) to (n, m), each interleaving of the two
    samples equally likely, and sums the probability of the paths that stay strictly inside
    the band |i/n - j/m| < d. ``checkpoints`` holds the pooled counts ``i + j`` at which the
    ECDFs are compared, i.e. the end of each block of tied values: with ties ``D`` is only
    observed there, so constraining every lattice point (the continuous-data null) would
    overstate the p-value. ``None`` means no ties (every point is a checkpoint).
    """
    # |i/n - j/m| < d  <=>  |i*m - j*n| < d*n*m; the epsilon keeps the observed D itself
    # outside the band (P(D >= d) includes D == d).
    bound = d * n * m - 1e-9 * n * m
    prev = [0.0] * (m + 1)
    for i in range(n + 1):
        row = [0.0] * (m + 1)
        for j in range(m + 1):
            k = i + j
            if k == 0:
                val = 1.0
            else:
                left = n + m - k + 1  # elements not yet placed before this step
                val = 0.0
                if i > 0:
                    val += prev[j] * (n - i + 1) / left
                if j > 0:
                    val += row[j - 1] * (m - j + 1) / left
            if (checkpoints is None or k in checkpoints) and abs(i * m - j * n) >= bound:
                val = 0.0
            row[j] = val
        prev = row
    return min(1.0, max(0.0, 1.0 - prev[m]))


def _ks_asymptotic_p(n: int, m: int, d: float) -> float:
    """Kolmogorov limiting distribution with Stephens' small-sample correction."""
    en = math.sqrt(n * m / (n + m))
    lam = (en + 0.12 + 0.11 / en) * d
    if lam < 1e-3:
        return 1.0
    s = sum((1.0 if k % 2 else -1.0) * math.exp(-2.0 * k * k * lam * lam) for k in range(1, 101))
    return min(1.0, max(0.0, 2 * s))


def ks_two_sample(a: Sequence[float], b: Sequence[float]) -> tuple[float, float]:
    """Two-sample Kolmogorov-Smirnov test (two-sided): ``(D, p)``.

    Pure Python. Exact permutation p-value when ``n * m <= 250_000``, ties included: the
    ECDFs are compared after every distinct value, and the null distribution is computed for
    the same tie blocks (POI offsets are small integers, mostly -1 or 0, so ties are the rule;
    the continuous-data null would inflate p and let a differing draw pass the lint). Larger
    samples use the asymptotic Kolmogorov distribution, which ignores ties and is therefore
    lenient (p too high) when they are present.
    """
    n, m = len(a), len(b)
    if not n or not m:
        raise ValueError("ks_two_sample needs two non-empty samples")
    xs, ys = sorted(a), sorted(b)
    d = 0.0
    i = j = 0
    checkpoints: set[int] = set()
    for v in sorted(set(xs) | set(ys)):
        while i < n and xs[i] <= v:
            i += 1
        while j < m and ys[j] <= v:
            j += 1
        checkpoints.add(i + j)
        d = max(d, abs(i / n - j / m))
    if d <= 0.0:
        return 0.0, 1.0
    if n * m <= _KS_EXACT_MAX_NM:
        tied = len(checkpoints) < n + m
        p = _ks_exact_p(n, m, d, checkpoints if tied else None)
    else:
        p = _ks_asymptotic_p(n, m, d)
    return d, p


def poi_offset(site0: int, run: tuple[int, int]) -> int:
    """Offset (bp) of the POI from the start of its repeat run (``run`` 0-based half-open).

    -1 is allowed: a left-normalised insertion of the run's unit is anchored on the base
    just before the run."""
    start, end = run
    if not start - 1 <= site0 < end:
        raise ValueError(f"POI {site0} is not inside the run [{start}, {end}) or just before it")
    return site0 - start


@dataclass
class OffsetLint:
    ok: bool
    d: float
    p: float
    n_pos: int
    n_neg: int
    problems: list[str] = field(default_factory=list)


def lint_poi_offsets(
    pos_offsets: Sequence[float], neg_offsets: Sequence[float], *, min_p: float = 0.2
) -> OffsetLint:
    """Section 3.4: the POI-offset-within-run distributions of ``tr_indel_pos`` (``pos``)
    and ``tr_stutter_neg`` (``neg``) must not differ: two-sample KS p > ``min_p``, else the
    draw is redone by seed (data-only)."""
    d, p = ks_two_sample(pos_offsets, neg_offsets)
    res = OffsetLint(p > min_p, d, p, len(pos_offsets), len(neg_offsets))
    if not res.ok:
        res.problems.append(
            f"POI offset within run differs between strata 3 and 4: KS D={d:.3f}, "
            f"p={p:.3f} <= {min_p} (n={len(pos_offsets)}/{len(neg_offsets)}); redraw by seed"
        )
    return res


def anchor_window(site0: int, run: tuple[int, int] | None, slop: int = 2) -> tuple[int, int]:
    """Half-open span the round-4 anchoring sentence lets an indel at the POI occupy: the
    repeat run containing the POI (if any) and ``slop`` bp either side of the POI."""
    lo, hi = site0 - slop, site0 + slop + 1
    if run is not None:
        lo, hi = min(lo, run[0]), max(hi, run[1])
    return lo, hi


def lint_anchor_collisions(rows: Iterable[dict[str, Any]]) -> list[str]:
    """Strata 7-8 (section 3.7): fail any item where the anchoring sentence could map its
    POI onto a truth record other than its own.

    Each row: ``site0``, optional ``run`` ([start, end) of the repeat run containing the
    POI), ``other_truth`` (list of [start, end) reference spans of every *other* truth
    record nearby, any set and FILTER) and an identifier in ``id`` or ``chrom``.
    """
    problems: list[str] = []
    for r in rows:
        run = tuple(r["run"]) if r.get("run") else None
        lo, hi = anchor_window(int(r["site0"]), run)  # type: ignore[arg-type]
        hits = [tuple(t) for t in r.get("other_truth", []) if t[0] < hi and lo < t[1]]
        if hits:
            name = r.get("id") or f"{r.get('chrom', '?')}:{int(r['site0']) + 1}"
            problems.append(
                f"{name}: anchoring window [{lo}, {hi}) overlaps truth record(s) {hits[:3]}"
            )
    return problems


_TOKEN = re.compile(r"[A-Za-z0-9_]{4,}")


def _tokens(text: str) -> set[str]:
    return set(_TOKEN.findall(text))
