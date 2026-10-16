"""Case types of the visual-vs-text eval and their ground-truth labels.

Each case type is planted at a seeded random position (see :mod:`.fixtures`). The label
fields mirror the answer schema in :mod:`.answers`. Label conventions (pre-registered in
``docs/eval/visual-preregistration.md``):

- ``variant_present`` — the simulated sample truly carries a variant at the site(s).
- ``artifact_verdict`` — ``real`` when it does; ``artifact`` when it does not (any
  non-reference signal is noise or a technical artifact).
- ``artifact_type`` — the systematic artifact mechanism, ``none`` for real variants and
  for negatives without a systematic signal.
- ``alt_allele`` — SNV: the alt base; MNV: the alt bases; deletion: the deleted
  reference base(s); insertion: the inserted bases; several alleles/sites: comma-joined
  and sorted. Empty when ``variant_present`` is ``no``.
- ``phase`` — only meaningful for the two-site phasing items, ``na`` otherwise.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any

DEFAULT_DEPTH = 30


@dataclass(frozen=True)
class CaseSpec:
    name: str
    variant_type: str  # SNV|ins|del|MNV|none
    zygosity: str  # het|hom_alt|low_vaf|none
    artifact_type: str  # none or an artifact mechanism
    phase: str = "na"
    depth: int = DEFAULT_DEPTH
    lenient_verdicts: tuple[str, ...] = ()

    @property
    def real(self) -> bool:
        return self.variant_type != "none"

    def label(self, alt_allele: str) -> dict[str, Any]:
        return {
            "variant_present": "yes" if self.real else "no",
            "variant_type": self.variant_type,
            "alt_allele": alt_allele if self.real else "",
            "zygosity": self.zygosity,
            "artifact_verdict": "real" if self.real else "artifact",
            "artifact_type": self.artifact_type,
            "phase": self.phase,
        }


CASE_TYPES: tuple[CaseSpec, ...] = (
    CaseSpec("het_snv", "SNV", "het", "none"),
    CaseSpec("hom_alt_snv", "SNV", "hom_alt", "none"),
    CaseSpec("low_vaf_snv", "SNV", "low_vaf", "none", depth=40),
    CaseSpec("low_baseq_alt", "none", "none", "low_baseq"),
    CaseSpec("strand_bias", "none", "none", "strand_bias"),
    CaseSpec("read_end_bias", "none", "none", "read_end"),
    CaseSpec("low_mapq_paralog", "none", "none", "low_mapq"),
    CaseSpec("homopolymer_del", "none", "none", "homopolymer_indel"),
    CaseSpec("real_del", "del", "het", "none"),
    CaseSpec("real_ins", "ins", "het", "none"),
    CaseSpec("multi_allelic", "SNV", "het", "none"),
    CaseSpec("pcr_duplicate", "none", "none", "pcr_duplicate"),
    CaseSpec("low_depth_het", "SNV", "het", "none", depth=3, lenient_verdicts=("uncertain",)),
    CaseSpec("clean_negative", "none", "none", "none"),
    CaseSpec("scattered_errors", "none", "none", "none"),
    CaseSpec("mnv", "MNV", "het", "none"),
    CaseSpec("phased_cis", "SNV", "het", "none", phase="cis"),
    CaseSpec("phased_trans", "SNV", "het", "none", phase="trans"),
)

CASE_BY_NAME = {c.name: c for c in CASE_TYPES}
