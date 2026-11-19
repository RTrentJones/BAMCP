"""The one candidate-filter predicate: min_vaf / min_depth / min_alt_count.

Every tool-layer filter over candidate variants (get_variants' selection, the viewer and
region-summary alt-read floor, scan_variants' per-chunk floor, the requested-site
``filter_check``) goes through :func:`failed_candidate_filters`, so they cannot disagree.

Semantics, decided once:

- A threshold passed as ``None`` is *not applied* (e.g. the viewer path only applies the
  alt-read floor; the parser already applied VAF/depth).
- A candidate's missing ``alt_count`` means "count unknown" (a VCF record without per-sample
  AD, say) — it is never dropped by ``min_alt_count``; only a known count below the floor
  fails. Missing ``vaf``/``depth`` read as 0 (no evidence), as before.
"""

from __future__ import annotations

from collections.abc import Mapping
from typing import Any


def failed_candidate_filters(
    v: Mapping[str, Any],
    *,
    min_vaf: float | None = None,
    min_depth: int | None = None,
    min_alt_count: int | None = None,
) -> list[str]:
    """The thresholds ``v`` fails, in the order ``min_vaf``, ``min_depth``, ``min_alt_count``."""
    failed: list[str] = []
    if min_vaf is not None and float(v.get("vaf") or 0.0) < min_vaf:
        failed.append("min_vaf")
    if min_depth is not None and int(v.get("depth") or 0) < min_depth:
        failed.append("min_depth")
    alt_count = v.get("alt_count")
    if min_alt_count is not None and alt_count is not None and int(alt_count) < min_alt_count:
        failed.append("min_alt_count")
    return failed


def passes_candidate_filters(
    v: Mapping[str, Any],
    *,
    min_vaf: float | None = None,
    min_depth: int | None = None,
    min_alt_count: int | None = None,
) -> bool:
    """Whether a candidate passes the tool-layer thresholds. VCF records always pass:
    the caller's VCF is authoritative, so BAMCP's floors never hide its calls."""
    if v.get("source") == "vcf":
        return True
    return not failed_candidate_filters(
        v, min_vaf=min_vaf, min_depth=min_depth, min_alt_count=min_alt_count
    )
