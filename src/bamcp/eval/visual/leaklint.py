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
"""

from __future__ import annotations

import json
import re
import zlib
from dataclasses import dataclass, field
from pathlib import Path

from .cases import CASE_TYPES
from .real import REAL_CASE_TYPES
from .stimuli import FILE_ORDER, prompt_text

ALLOWED_FILES = frozenset(
    {
        "caption.md",
        "prompt.md",
        "traw.md",
        "tjson.json",
        "tfull.md",
        "view_wide.png",
        "view_zoom.png",
    }
)
VOCAB_EXEMPT = frozenset({"prompt.md", "tfull.md"})
_OPAQUE_DIR = re.compile(r"^[0-9a-f]{8}$")
LEAK_VOCAB = re.compile(
    r"strand|bias|homopolymer|artifact|artefact|mapq_low|low_mapq|low_baseq|low_bq|dup|"
    r"pcr|risk|imbalanced|alt_fwd|alt_rev|_ref_|paralog|verdict|concern|recommend|"
    r"confidence|\bhet\b|heterozyg|\bhom\b|homozyg|hom_alt|\bcis\b|\btrans\b|subclonal|"
    r"low_vaf|\bnoise\b|\berror|misalign|read_end|negative|scatter|\breal\b|\bmnv\b|"
    r"multi_allelic|phas",
    re.IGNORECASE,
)
CAPTION_RE = re.compile(
    r"^# Site\n\n"
    r"- Region: chr\d+:\d+-\d+ \(1-based, inclusive\)\n"
    r"- Position\(s\) of interest: chr\d+:\d+(, chr\d+:\d+)?\n"
    r"- Sample: (simulated diploid sample on a synthetic GRCh38-like reference"
    r"|sequenced human sample \(Illumina short reads\); coordinates are re-indexed)\n$"
)


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
    case_names = {c.name for c in CASE_TYPES} | set(REAL_CASE_TYPES)
    meta = run_dir / "build.json"
    sample_kind = (
        "sequenced"
        if meta.exists() and json.loads(meta.read_text()).get("tier") == "real"
        else "simulated"
    )

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
        if path.name == "caption.md" and not CAPTION_RE.match(text):
            report.problems.append(f"{rel}: caption does not match the uniform template")
        if path.name == "prompt.md":
            # The shared prompt names the answer vocabulary on purpose; instead of a
            # vocabulary scan it must be exactly the uniform text for its file set.
            present = sorted(p.name for p in path.parent.iterdir() if p.name in FILE_ORDER)
            present.sort(key=FILE_ORDER.index)
            if text != prompt_text(present, sample_kind):
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


_TOKEN = re.compile(r"[A-Za-z0-9_]{4,}")


def _tokens(text: str) -> set[str]:
    return set(_TOKEN.findall(text))
