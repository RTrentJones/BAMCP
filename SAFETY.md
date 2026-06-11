# Safety & Intended Use

BAMCP makes genomic alignment evidence legible to an LLM and to a human
through an interactive viewer. That puts it one step away from clinical
decision-making, so the boundaries below are deliberate, documented, and —
where possible — **enforced by tests**, not just stated.

## Intended use

BAMCP is a **research and visualization tool**. It is **not** a diagnostic
device, not a clinical variant caller, and not a substitute for a validated
pipeline or expert review.

- Variant detection is a lightweight pileup caller for *exploration*. It is not
  benchmarked or regulated as a clinical caller (compare DeepVariant, GATK).
- The curation summary surfaces *quality evidence and risk flags* to support a
  human curator. It does not make a pathogenicity or actionability call.
- Every variant, regardless of how clean it looks, requires confirmation by a
  qualified human using a validated process.

This boundary travels with the data: `get_variant_curation_summary` includes an
`intended_use` field in its structured output and an `INTENDED USE` section in
its text output, so a downstream LLM sees the disclaimer in-band rather than
relying on it being in the docs.

## Why this matters

An LLM reading tool output can sound authoritative. The failure mode we most
want to avoid is **confident presentation of weak evidence** — e.g. calling a
1×-depth site, a strand-biased PCR artifact, or a homopolymer slip as if it
were a real, high-confidence variant. The design pushes against this in two
places:

1. **The tool hedges by construction.** `analysis/evidence.py` computes
   artifact risk (strand bias, read-position bias, low MAPQ, homopolymer
   context) and `compute_confidence` downgrades confidence accordingly. Clean
   "no flags" output explicitly says manual curation is still required rather
   than claiming clinical suitability.

2. **An eval guards the hedge.** The ground-truth gate
   (`bamcp.eval.truthset`, run on every PR) includes a **safety invariant**:
   no site flagged as high artifact risk may simultaneously be reported as
   high confidence. A regression that made the tool overconfident on
   likely-spurious evidence fails CI. See [EVALS.md](EVALS.md).

## Safety invariants under test

| Invariant | Enforced by |
| --------- | ----------- |
| No high-artifact site is high-confidence | `truthset` `overconfident_sites` floor (PR gate) |
| Low-coverage variants are flagged, not silently passed | `tests/unit/analysis/test_curation_safety.py` |
| "No concerns" output never claims clinical suitability | `tests/unit/analysis/test_curation_safety.py` |
| The intended-use disclaimer is present in every curation payload | `tests/unit/analysis/test_curation_rubric.py`, `test_curation_safety.py` |

## Data handling

Alignment files can contain identifiable human genomic data. BAMCP:

- does not persist read data — regions are fetched, serialized to the caller,
  and dropped; only remote *index* files are cached (24h TTL, session-isolated);
- treats remote access as opt-in and SSRF-filtered (see [SECURITY.md](SECURITY.md));
- emits sanitized errors that avoid echoing file paths or config.

Operators handling protected health information remain responsible for running
BAMCP within an environment that meets their regulatory obligations (e.g.
HIPAA, GDPR). BAMCP provides controls, not compliance.

## Reporting a safety concern

Safety issues that are not security vulnerabilities (e.g. a tool over-claiming
confidence, a misleading recommendation) can be raised as a GitHub issue. If a
concern has both safety and security dimensions, follow the private disclosure
path in [SECURITY.md](SECURITY.md).
