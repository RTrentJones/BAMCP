"""Round-4 strata: frozen quotas, knobs, per-stratum predicates, spacing, matching, labels.

Everything here is data-only and frozen before the pilot (design section 3 and the knob
grid of section 7). Candidates are manifest rows (dicts, see ``CONTRACT.md``): the
generators in :mod:`.real` fill ``features`` / ``context`` from the 60x Illumina reads, and
this module decides, without ever looking at a model answer, whether a row is eligible for
its stratum at a knob level, whether it is recoverable, which matching bin it falls in,
and which rows a seeded draw takes.

Pairs (the primary endpoint averages balanced accuracy over them):

=========  =============================  ==================================
pair       positive                       negative
=========  =============================  ==================================
mosaic     ``mosaic_snv_pos``             ``mosaic_caller_fp_neg``
tr         ``tr_indel_pos``               ``tr_stutter_neg``
segdup     ``segdup_cmrg_pos``            ``segdup_psv_neg``
indel_adj  ``indel_adjacent_snv_pos``     ``indel_adjacent_misalign_neg``
anchor     ``anchor_het_snv``, ``anchor_hom_snv``  ``anchor_clean_neg``
=========  =============================  ==================================
"""

from __future__ import annotations

import hashlib
import math
import random
from collections import Counter
from collections.abc import Callable, Iterable
from typing import Any

from .real import NEGATIVE_LABEL

STRATA_R4 = (
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
PAIRS: dict[str, tuple[str, str]] = {
    "mosaic_snv_pos": ("mosaic", "pos"),
    "mosaic_caller_fp_neg": ("mosaic", "neg"),
    "tr_indel_pos": ("tr", "pos"),
    "tr_stutter_neg": ("tr", "neg"),
    "segdup_cmrg_pos": ("segdup", "pos"),
    "segdup_psv_neg": ("segdup", "neg"),
    "indel_adjacent_misalign_neg": ("indel_adj", "neg"),
    "indel_adjacent_snv_pos": ("indel_adj", "pos"),
    "anchor_het_snv": ("anchor", "pos"),
    "anchor_hom_snv": ("anchor", "pos"),
    "anchor_clean_neg": ("anchor", "neg"),
}
HARD_STRATA = tuple(s for s in STRATA_R4 if not s.startswith("anchor_"))
ANCHOR_STRATA = tuple(s for s in STRATA_R4 if s.startswith("anchor_"))
PAIR_NAMES = ("mosaic", "tr", "segdup", "indel_adj")
# Negative stratum -> the positive stratum it is matched to (one-to-one, by bin).
MATCHED_TO = {
    "mosaic_caller_fp_neg": "mosaic_snv_pos",
    "tr_stutter_neg": "tr_indel_pos",
    "segdup_psv_neg": "segdup_cmrg_pos",
}

# Frozen quotas (design 3 / 7). Pilot: per hard stratum 6 at the default knob + 4 at the
# alternative; anchors 4 (1 het, 1 hom, 2 clean: half positive, as in the main set's
# 5 positive / 3 negative split rounded toward balance). Main: 15 per hard stratum; anchors
# 2 het / 3 hom / 3 clean.
QUOTAS_R4: dict[str, dict[str, dict[str, int]]] = {
    "main": {
        "default": {
            **{s: 15 for s in HARD_STRATA},
            "anchor_het_snv": 2,
            "anchor_hom_snv": 3,
            "anchor_clean_neg": 3,
        },
    },
    "pilot": {
        "default": {
            **{s: 6 for s in HARD_STRATA},
            "anchor_het_snv": 1,
            "anchor_hom_snv": 1,
            "anchor_clean_neg": 2,
        },
        "alt": {s: 4 for s in HARD_STRATA},
    },
}

# Knob grid (design 7): default -> the piloted alternative. Stratum 3 has two knobs (run
# floor and difficulty cut); its alternative level moves both, so the pilot observes one
# alternative cell per stratum.
KNOBS: dict[str, dict[str, dict[str, Any]]] = {
    "mosaic_snv_pos": {"default": {"alt_band": "split"}, "alt": {"alt_band": "ge7"}},
    "mosaic_caller_fp_neg": {
        "default": {"vaf_lo": 0.05, "vaf_hi": 0.20},
        "alt": {"vaf_lo": 0.08, "vaf_hi": 0.20},
    },
    "tr_indel_pos": {
        "default": {"hp_min": 10, "stutter_min": 0.08},
        "alt": {"hp_min": 12, "stutter_min": 0.12},
    },
    "tr_stutter_neg": {"default": {"hp_min": 10}, "alt": {"hp_min": 12}},
    "segdup_cmrg_pos": {"default": {"difficult_share": 0.5}, "alt": {"difficult_share": 0.75}},
    "segdup_psv_neg": {"default": {"vaf_min": 0.10}, "alt": {"vaf_min": 0.07}},
    "indel_adjacent_misalign_neg": {"default": {"vaf_min": 0.10}, "alt": {"vaf_min": 0.15}},
    "indel_adjacent_snv_pos": {
        "default": {"alt_min": 5, "alt_max": None},
        "alt": {"alt_min": 4, "alt_max": 8},
    },
}
KNOB_LEVELS = ("default", "alt")

# Section 3.0 rules.
MIN_SPACING_R4 = 10_000
MAX_PER_CLUSTER = 2
MIN_DEPTH = 10
STR_MIN = 20
CONTAMINATION_VAF_FLOOR = 0.08  # stratum 2's floor if contamination > 1% (design 3.2)
GNOMAD_MAX_AF = 0.01
# Per-stratum recoverability / signal constants.
MOSAIC_MIN_ALT = 4
MOSAIC_MIN_DEPTH = 25
MOSAIC_MIN_VAF = 0.06
HET_MIN_ALT = 3
HET_MIN_VAF = 0.10
HOM_MIN_FRAC = 0.60
FP_MIN_ALT = 4
# Amendment A1 (pre-pilot): 4 reads / 10% left 2 of 299 scanned confident runs (stutter in
# benchmark-confident homopolymers averages 1-3%); 3 reads / 5% keeps about 5%.
TR_NEG_MIN_READS = 3
TR_NEG_MIN_FRAC = 0.05
TR_HET_ALLELE_MAX = 0.40
TR_HOM_ALLELE_MAX = 0.85
PSV_MIN_READS = 3
ADJ_MIN_READS = 4
CMRG_DIFFICULT_VAF = 0.30
CMRG_DIFFICULT_MQ = 0.30
TR_POS_HOM_SHARE = 1 / 3  # stratum 3: 10 het + 5 hom in the main set


def pair_of(stratum: str) -> tuple[str, str]:
    return PAIRS[stratum]


def quotas(split: str) -> dict[str, dict[str, int]]:
    """``{knob_level: {stratum: n}}`` for ``split`` ("pilot" or "main")."""
    if split not in QUOTAS_R4:
        raise ValueError(f"unknown split {split!r}; choose pilot or main")
    return {k: dict(v) for k, v in QUOTAS_R4[split].items()}


# ---------------------------------------------------------------------------
# Labels
# ---------------------------------------------------------------------------


def label_for(stratum: str, row: dict[str, Any]) -> dict[str, Any]:
    """The answer-schema label of a row (design 3.0 "Labels").

    Positives are ``yes`` / ``real`` with the truth allele; mosaic positives are
    ``zygosity: low_vaf`` (an allowed answer value). Negatives are ``no`` / ``artifact``;
    ``artifact_type`` is ``unknown`` (unscored) for every negative, anchors' clean negative
    included (design 3.0, CONTRACT.md), except TR stutter (``homopolymer_indel``).
    """
    _, polarity = PAIRS[stratum]
    if polarity == "neg":
        atype = "homopolymer_indel" if stratum == "tr_stutter_neg" else "unknown"
        return {**NEGATIVE_LABEL, "artifact_type": atype}
    base = dict(row.get("label") or {})
    lab = {
        "variant_present": "yes",
        "variant_type": base.get("variant_type", "SNV"),
        "alt_allele": base.get("alt_allele", row.get("vcf_alt", "")),
        "zygosity": base.get("zygosity", "het"),
        "artifact_verdict": "real",
        "artifact_type": "none",
        "phase": "na",
    }
    if stratum == "mosaic_snv_pos":
        lab["zygosity"] = "low_vaf"
    return lab


# ---------------------------------------------------------------------------
# Recoverability and predicates
# ---------------------------------------------------------------------------


def _f(row: dict[str, Any], key: str, default: float = 0.0) -> float:
    v = (row.get("features") or {}).get(key)
    return float(default if v is None else v)


def recoverable(row: dict[str, Any]) -> tuple[bool, str | None]:
    """Design 3.0 recoverability: data only, applied to every positive.

    het: >= 3 alt reads (>= 4 mosaic) and VAF >= 0.10 (mosaic 0.06); hom: alt fraction >=
    0.60; any: depth >= 10. For repeat indels the counts are spanning reads carrying the
    truth length.
    """
    depth, alt, vaf = _f(row, "depth"), _f(row, "alt_reads"), _f(row, "vaf")
    if depth < MIN_DEPTH:
        return False, "depth_lt_10"
    zyg = (row.get("label") or {}).get("zygosity")
    if row.get("stratum") == "mosaic_snv_pos" or zyg == "low_vaf":
        if alt < MOSAIC_MIN_ALT or vaf < MOSAIC_MIN_VAF:
            return False, "mosaic_unrecoverable"
        return True, None
    if zyg == "hom_alt":
        return (vaf >= HOM_MIN_FRAC, None if vaf >= HOM_MIN_FRAC else "hom_frac_lt_0.60")
    if alt < HET_MIN_ALT or vaf < HET_MIN_VAF:
        return False, "het_unrecoverable"
    return True, None


def alt_band(alt: float) -> str:
    return "4-6" if alt <= 6 else "ge7"


def _alt_bin(alt: float) -> str:
    return "4-6" if alt <= 6 else "7-9" if alt <= 9 else "ge10"


def _depth_bin(depth: float) -> str:
    return "lt50" if depth < 50 else "50-70" if depth <= 70 else "gt70"


def _mq_bin(frac: float) -> str:
    return "lt0.10" if frac < 0.10 else "0.10-0.30" if frac < 0.30 else "ge0.30"


def _run_bin(row: dict[str, Any]) -> str:
    ctx = row.get("context") or {}
    unit = ctx.get("unit") or ""
    n = int(_f(row, "run_length"))
    if len(unit) > 1:
        return f"str{len(unit)}:" + ("20-29" if n < 30 else "ge30")
    return "hp:" + ("10-11" if n < 12 else "12-15" if n < 16 else "ge16")


def cmrg_difficult(row: dict[str, Any]) -> bool:
    """Design 3.5 difficulty cut: het VAF < 0.30 or >= 30% of reads at MAPQ < 20."""
    zyg = (row.get("label") or {}).get("zygosity")
    het_low = zyg == "het" and _f(row, "vaf") < CMRG_DIFFICULT_VAF
    return het_low or mapq_lt20_frac(row) >= CMRG_DIFFICULT_MQ


def mapq_lt20_frac(row: dict[str, Any]) -> float:
    """Fraction of all reads at the site with MAPQ < 20 (alt and ref classes pooled)."""
    d = _f(row, "depth")
    if not d:
        return 0.0
    a = _f(row, "alt_reads")
    return (_f(row, "alt_mapq_lt20_frac") * a + _f(row, "ref_mapq_lt20_frac") * max(d - a, 0)) / d


def _tr_run_ok(row: dict[str, Any], hp_min: int) -> bool:
    from .truth import is_tr_run

    unit = (row.get("context") or {}).get("unit") or ""
    return bool(unit) and is_tr_run(int(_f(row, "run_length")), unit, hp_min, STR_MIN)


Predicate = Callable[[dict[str, Any], dict[str, Any]], "tuple[bool, str | None]"]


def _p_mosaic_pos(row: dict[str, Any], knob: dict[str, Any]) -> tuple[bool, str | None]:
    if _f(row, "depth") < MOSAIC_MIN_DEPTH:
        return False, "depth_lt_25"
    ok, why = recoverable(row)
    if not ok:
        return ok, why
    if knob.get("alt_band") == "ge7" and _f(row, "alt_reads") < 7:
        return False, "knob_alt_band"
    return True, None


def _p_mosaic_fp(row: dict[str, Any], knob: dict[str, Any]) -> tuple[bool, str | None]:
    lo = max(float(knob["vaf_lo"]), float(knob.get("vaf_floor_override") or 0.0))
    if _f(row, "alt_reads") < FP_MIN_ALT:
        return False, "alt_lt_4"
    if not lo <= _f(row, "vaf") <= float(knob["vaf_hi"]):
        return False, "knob_vaf_band"
    if _f(row, "depth") < MIN_DEPTH:
        return False, "depth_lt_10"
    return True, None


def _p_tr_pos(row: dict[str, Any], knob: dict[str, Any]) -> tuple[bool, str | None]:
    if not _tr_run_ok(row, int(knob["hp_min"])):
        return False, "knob_run_floor"
    ok, why = recoverable(row)
    if not ok:
        return ok, why
    zyg = (row.get("label") or {}).get("zygosity")
    allele_max = TR_HOM_ALLELE_MAX if zyg == "hom_alt" else TR_HET_ALLELE_MAX
    if _f(row, "stutter_frac") >= float(knob["stutter_min"]) or _f(row, "vaf") <= allele_max:
        return True, None
    return False, "not_difficult"


def _p_tr_neg(row: dict[str, Any], knob: dict[str, Any]) -> tuple[bool, str | None]:
    if not _tr_run_ok(row, int(knob["hp_min"])):
        return False, "knob_run_floor"
    if _f(row, "alt_reads") < TR_NEG_MIN_READS or _f(row, "vaf") < TR_NEG_MIN_FRAC:
        return False, "stutter_signal_weak"
    return True, None


def _p_cmrg_pos(row: dict[str, Any], knob: dict[str, Any]) -> tuple[bool, str | None]:
    return recoverable(row)


def _p_psv_neg(row: dict[str, Any], knob: dict[str, Any]) -> tuple[bool, str | None]:
    if _f(row, "depth") < MIN_DEPTH:
        return False, "depth_lt_10"
    if _f(row, "alt_reads") < PSV_MIN_READS or _f(row, "vaf") < float(knob["vaf_min"]):
        return False, "knob_vaf_floor"
    return True, None


def _p_adj_neg(row: dict[str, Any], knob: dict[str, Any]) -> tuple[bool, str | None]:
    if _f(row, "depth") < MIN_DEPTH:
        return False, "depth_lt_10"
    if _f(row, "alt_reads") < ADJ_MIN_READS or _f(row, "vaf") < float(knob["vaf_min"]):
        return False, "knob_vaf_floor"
    return True, None


def _p_adj_pos(row: dict[str, Any], knob: dict[str, Any]) -> tuple[bool, str | None]:
    ok, why = recoverable(row)
    if not ok:
        return ok, why
    alt = _f(row, "alt_reads")
    hi = knob.get("alt_max")
    if alt < int(knob["alt_min"]) or (hi is not None and alt > int(hi)):
        return False, "knob_alt_reads"
    return True, None


def _p_anchor_pos(row: dict[str, Any], knob: dict[str, Any]) -> tuple[bool, str | None]:
    return recoverable(row)


def _p_anchor_neg(row: dict[str, Any], knob: dict[str, Any]) -> tuple[bool, str | None]:
    if _f(row, "depth") < 20:
        return False, "depth_lt_20"
    if _f(row, "alt_reads") > 0:
        return False, "not_clean"
    return True, None


PREDICATES: dict[str, Predicate] = {
    "mosaic_snv_pos": _p_mosaic_pos,
    "mosaic_caller_fp_neg": _p_mosaic_fp,
    "tr_indel_pos": _p_tr_pos,
    "tr_stutter_neg": _p_tr_neg,
    "segdup_cmrg_pos": _p_cmrg_pos,
    "segdup_psv_neg": _p_psv_neg,
    "indel_adjacent_misalign_neg": _p_adj_neg,
    "indel_adjacent_snv_pos": _p_adj_pos,
    "anchor_het_snv": _p_anchor_pos,
    "anchor_hom_snv": _p_anchor_pos,
    "anchor_clean_neg": _p_anchor_neg,
}


def knob_params(
    stratum: str, level: str, overrides: dict[str, Any] | None = None
) -> dict[str, Any]:
    params = dict(KNOBS.get(stratum, {}).get(level, {}))
    params.update((overrides or {}).get(stratum, {}))
    return params


def eligible(
    row: dict[str, Any], level: str, overrides: dict[str, Any] | None = None
) -> tuple[bool, str | None]:
    """Does ``row`` pass its stratum's predicate at knob ``level``?"""
    stratum = row["stratum"]
    return PREDICATES[stratum](row, knob_params(stratum, level, overrides))


# ---------------------------------------------------------------------------
# Sub-quotas within a stratum
# ---------------------------------------------------------------------------


def subkey(stratum: str, row: dict[str, Any]) -> str:
    """The within-stratum cell a row fills (alt band, zygosity, difficulty)."""
    if stratum == "mosaic_snv_pos":
        return alt_band(_f(row, "alt_reads"))
    if stratum == "tr_indel_pos":
        return "hom" if (row.get("label") or {}).get("zygosity") == "hom_alt" else "het"
    if stratum == "segdup_cmrg_pos":
        return "difficult" if cmrg_difficult(row) else "other"
    return "any"


def subquotas(
    stratum: str, level: str, n: int, overrides: dict[str, Any] | None = None
) -> dict[str, int]:
    """Split ``n`` items of ``stratum`` across its cells at knob ``level``."""
    knob = knob_params(stratum, level, overrides)
    if stratum == "mosaic_snv_pos":
        if knob.get("alt_band") == "ge7":
            return {"ge7": n}
        return {"4-6": n // 2, "ge7": n - n // 2}
    if stratum == "tr_indel_pos":
        hom = round(n * TR_POS_HOM_SHARE)
        return {"het": n - hom, "hom": hom}
    if stratum == "segdup_cmrg_pos":
        k = math.ceil(n * float(knob["difficult_share"]) - 1e-9)
        return {"difficult": k, "other": n - k}
    return {"any": n}


# ---------------------------------------------------------------------------
# Spacing and matching
# ---------------------------------------------------------------------------


class Spacing:
    """Design 3.0: items >= 10 kb apart on a chromosome, <= 2 per cluster (global)."""

    def __init__(self, min_bp: int = MIN_SPACING_R4, per_cluster: int = MAX_PER_CLUSTER):
        self.min_bp = min_bp
        self.per_cluster = per_cluster
        self.taken: list[dict[str, Any]] = []
        self.per: Counter[str] = Counter()

    def why_not(self, row: dict[str, Any]) -> str | None:
        if self.per[row["cluster"]] >= self.per_cluster:
            return "cluster_full"
        for t in self.taken:
            if (
                t["chrom"] == row["chrom"]
                and abs(int(t["site0"]) - int(row["site0"])) < self.min_bp
            ):
                return "spacing_10kb"
        return None

    def ok(self, row: dict[str, Any]) -> bool:
        return self.why_not(row) is None

    def add(self, row: dict[str, Any]) -> None:
        self.taken.append(row)
        self.per[row["cluster"]] += 1


def match_key(neg_stratum: str, row: dict[str, Any], *, level: int = 0) -> tuple[Any, ...]:
    """Matching bin of a row for ``neg_stratum``'s pair; ``level`` > 0 relaxes the key
    (drops the last dimension first), the pre-declared order when an exact match is absent.

    - mosaic: alt-read bin (4-6, 7-9, >= 10), depth bin (< 50, 50-70, > 70), segdup /
      low-mappability membership;
    - segdup: MAPQ < 20 fraction bin (< 10%, 10-30%, >= 30%), depth bin;
    - tr: run-length bin (homopolymer 10-11 / 12-15 / >= 16, STR by unit).
    """
    if neg_stratum == "mosaic_caller_fp_neg":
        key: tuple[Any, ...] = (
            _alt_bin(_f(row, "alt_reads")),
            _depth_bin(_f(row, "depth")),
            bool((row.get("context") or {}).get("segdup_or_lowmap")),
        )
    elif neg_stratum == "segdup_psv_neg":
        key = (_mq_bin(mapq_lt20_frac(row)), _depth_bin(_f(row, "depth")))
    elif neg_stratum == "tr_stutter_neg":
        key = (_run_bin(row),)
    else:
        key = ()
    return key[: max(0, len(key) - level)]


def seeded_order(rows: Iterable[dict[str, Any]], seed: int, tag: str) -> list[dict[str, Any]]:
    """Rows in a seeded, reproducible order (independent of input order)."""
    out = sorted(rows, key=lambda r: (r["chrom"], int(r["site0"]), r.get("stratum", "")))
    random.Random(f"{seed}:{tag}").shuffle(out)  # noqa: S311
    return out


def match_negatives(
    positives: list[dict[str, Any]],
    candidates: list[dict[str, Any]],
    neg_stratum: str,
    *,
    seed: int,
    accept: Callable[[dict[str, Any]], str | None] = lambda r: None,
    on_take: Callable[[dict[str, Any]], Any] | None = None,
) -> tuple[list[dict[str, Any]], dict[str, int]]:
    """One negative per positive with the same :func:`match_key`, drawn in seeded order.

    ``accept(row)`` returns None to take a row or an exclusion reason (lazy label checks,
    spacing); every rejected row is counted. When no exact match remains, the key is
    relaxed one dimension at a time and the row's ``match`` records the level used.
    ``on_take(row)`` runs as soon as a row is chosen (so spacing state seen by ``accept``
    includes earlier matches) and its return value replaces the row in the output.
    Returns (chosen, exclusion counts).
    """
    pool = seeded_order(candidates, seed, f"match:{neg_stratum}")
    used: set[int] = set()
    rejected: set[int] = set()
    counts: Counter[str] = Counter()
    chosen: list[dict[str, Any]] = []
    depth = len(match_key(neg_stratum, positives[0])) if positives else 0
    for pos in positives:
        got = None
        for level in range(depth + 1):
            want = match_key(neg_stratum, pos, level=level)
            for i, cand in enumerate(pool):
                if i in used or i in rejected or match_key(neg_stratum, cand, level=level) != want:
                    continue
                why = accept(cand)
                if why is not None:
                    counts[why] += 1
                    if why != "spacing_10kb" and why != "cluster_full":
                        rejected.add(i)
                    continue
                got = {**cand, "match": {"level": level, "to": f"{pos['chrom']}:{pos['site0']}"}}
                used.add(i)
                break
            if got is not None:
                break
        if got is None:
            counts["unmatched"] += 1
            continue
        chosen.append(on_take(got) if on_take is not None else got)
    return chosen, dict(counts)


def offset_in_run(row: dict[str, Any]) -> float | None:
    """POI offset within its repeat run as a fraction of the run (TR strata lint)."""
    ctx = row.get("context") or {}
    run = ctx.get("run")
    if not run:
        return None
    s, e = int(run[0]), int(run[1])
    return (int(row["site0"]) - s) / max(e - s, 1)


def ks_two_sample(a: list[float], b: list[float]) -> tuple[float, float]:
    """Two-sample Kolmogorov-Smirnov statistic and asymptotic p-value (no scipy)."""
    if not a or not b:
        return 0.0, 1.0
    xs = sorted(set(a) | set(b))
    sa, sb = sorted(a), sorted(b)
    import bisect

    d = max(
        abs(bisect.bisect_right(sa, x) / len(sa) - bisect.bisect_right(sb, x) / len(sb)) for x in xs
    )
    n = len(a) * len(b) / (len(a) + len(b))
    lam = (math.sqrt(n) + 0.12 + 0.11 / math.sqrt(n)) * d
    if lam < 0.2:  # the series is 1 to machine precision here (and converges badly)
        return d, 1.0
    p = 2 * sum((-1) ** (k - 1) * math.exp(-2 * k * k * lam * lam) for k in range(1, 101))
    return d, min(1.0, max(0.0, p))


def manifest_digest(text: str) -> str:
    return hashlib.sha256(text.encode()).hexdigest()
