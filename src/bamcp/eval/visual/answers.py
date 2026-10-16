"""Answer schema, validation and normalization for the visual-vs-text eval.

An answer that fails validation is scored as wrong on every field (pre-registered).
"""

from __future__ import annotations

import json
import re
from dataclasses import dataclass
from typing import Any

ENUMS: dict[str, tuple[str, ...]] = {
    "variant_present": ("yes", "no"),
    "variant_type": ("SNV", "ins", "del", "MNV", "none"),
    "zygosity": ("het", "hom_alt", "low_vaf", "none"),
    "artifact_verdict": ("real", "artifact", "uncertain"),
    "artifact_type": (
        "none",
        "strand_bias",
        "read_end",
        "low_mapq",
        "low_baseq",
        "homopolymer_indel",
        "pcr_duplicate",
        "misalignment",
        "low_depth",
    ),
    "phase": ("cis", "trans", "na"),
}
FIELDS = (
    "variant_present",
    "variant_type",
    "alt_allele",
    "zygosity",
    "artifact_verdict",
    "artifact_type",
    "phase",
    "confidence",
    "rationale",
)
SCORED_FIELDS = (
    "variant_present",
    "variant_type",
    "alt_allele",
    "zygosity",
    "artifact_verdict",
    "artifact_type",
    "phase",
)
MAX_RATIONALE_WORDS = 60  # prompt asks for <=40; 50% tolerance before it counts as invalid
_NONE_ALLELES = {"", "NONE", "NA", "N/A", "-", ".", "NULL"}

ANSWER_SCHEMA_TEXT = """\
{
  "variant_present": "yes" | "no",
  "variant_type": "SNV" | "ins" | "del" | "MNV" | "none",
  "alt_allele": "<string>",
  "zygosity": "het" | "hom_alt" | "low_vaf" | "none",
  "artifact_verdict": "real" | "artifact" | "uncertain",
  "artifact_type": "none" | "strand_bias" | "read_end" | "low_mapq" | "low_baseq" |
                   "homopolymer_indel" | "pcr_duplicate" | "misalignment" | "low_depth",
  "phase": "cis" | "trans" | "na",
  "confidence": <number between 0 and 1>,
  "rationale": "<at most 40 words>"
}"""

FIELD_GUIDE = """\
- variant_present: "yes" if the sample truly carries a variant at the position(s) of interest.
- variant_type: type of that variant; "none" when variant_present is "no".
- alt_allele: SNV: the alternate base (e.g. "T"); MNV: the alternate bases in order (e.g. "TG");
  deletion: the deleted reference base(s); insertion: the inserted bases. Several alleles or
  several positions: comma-separated (e.g. "A,G"). Empty string when variant_present is "no".
- zygosity: "het" (about half the reads), "hom_alt" (essentially all reads), "low_vaf"
  (a real variant clearly below heterozygous fraction), "none" when no variant.
- artifact_verdict: "real" if you conclude a true variant is present; "artifact" if you
  conclude no true variant is present (any non-reference signal is noise or a technical
  artifact); "uncertain" if you cannot decide.
- artifact_type: the technical mechanism if you judge the signal to be an artifact, else "none".
- phase: only when two positions of interest are listed: "cis" if the two variants sit on the
  same reads/haplotype, "trans" if on different ones; otherwise "na".
- confidence: your probability (0-1) that your variant_present and artifact_verdict are correct.
- rationale: at most 40 words."""


@dataclass
class Validation:
    ok: bool
    answer: dict[str, Any] | None
    error: str | None = None


def _strip_fence(text: str) -> str:
    t = text.strip()
    m = re.fullmatch(r"```(?:json)?\s*(.*?)\s*```", t, flags=re.DOTALL)
    return m.group(1) if m else t


def parse_answer(text: str) -> Validation:
    """Parse + validate raw answer text (a single optional code fence is tolerated)."""
    try:
        obj = json.loads(_strip_fence(text))
    except (json.JSONDecodeError, ValueError) as e:
        return Validation(False, None, f"not JSON: {e}")
    return validate_answer(obj)


def validate_answer(obj: Any) -> Validation:
    if not isinstance(obj, dict):
        return Validation(False, None, "answer is not a JSON object")
    missing = [f for f in FIELDS if f not in obj]
    if missing:
        return Validation(False, None, f"missing fields: {', '.join(missing)}")
    for name, allowed in ENUMS.items():
        if obj[name] not in allowed:
            return Validation(False, None, f"{name}={obj[name]!r} not in {allowed}")
    if obj["alt_allele"] is not None and not isinstance(obj["alt_allele"], str):
        return Validation(False, None, "alt_allele must be a string")
    conf = obj["confidence"]
    if isinstance(conf, bool) or not isinstance(conf, (int, float)) or not 0 <= conf <= 1:
        return Validation(False, None, "confidence must be a number in [0, 1]")
    if not isinstance(obj["rationale"], str):
        return Validation(False, None, "rationale must be a string")
    if len(obj["rationale"].split()) > MAX_RATIONALE_WORDS:
        return Validation(False, None, f"rationale longer than {MAX_RATIONALE_WORDS} words")
    return Validation(True, obj)


def normalize_allele(value: Any) -> str:
    """Canonical comparison form: uppercase, no +/- markers, sorted comma list."""
    if value is None:
        return ""
    s = str(value).upper().replace(" ", "")
    if s in _NONE_ALLELES:
        return ""
    parts = [p.lstrip("+-") for p in re.split(r"[,/]", s) if p and p not in _NONE_ALLELES]
    return ",".join(sorted(parts))


def field_correct(name: str, answer: dict[str, Any], label: dict[str, Any]) -> bool:
    if name == "alt_allele":
        return normalize_allele(answer.get(name)) == normalize_allele(label.get(name))
    return answer.get(name) == label.get(name)
