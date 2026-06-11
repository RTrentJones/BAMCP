# BAMCP Evaluation

BAMCP ships with two complementary evaluation tracks:

1. **A deterministic ground-truth gate** (`bamcp.eval.truthset`) — no LLM, no
   network. It scores the BAMCP tools against a labeled truth set and runs on
   every PR. This is the regression gate that keeps the tools honest.
2. **An LLM-driven harness** (`bamcp.eval` runner) — drives a real model
   (Anthropic / OpenAI) through the tools and grades the conversation, with an
   optional vision path that feeds rendered screenshots back to the model.

This document covers the methodology. To run things, see the [Commands](#commands).

---

## 1. Ground-truth gate

### What it measures

The gate loads a versioned **truth-set manifest** and computes three families
of deterministic metrics:

| Metric | Question it answers | Floor (synthetic_v1) |
|--------|---------------------|----------------------|
| Variant detection P / R / F1 | Does the detector recover exactly the planted variant alleles in a labeled region? | recall = 1.0, precision ≥ 0.95 |
| False-positive control | Does it stay silent in reference-only regions? | 0 calls across negative regions |
| Artifact-type recall | Does the curation tool surface the *expected* artifact at each known-artifact site (strand bias, low MAPQ, homopolymer)? | recall = 1.0 |
| Clean-control discrimination | Is a true clean variant scored *below* the artifact-prone sites and not labeled high-risk? | true |

Floors live in code (`TruthsetReport.meets_floors`) and in the CI step, so
tightening or loosening them is a reviewed change.

### Why synthetic ground truth first

`tests/eval/datasets/synthetic_v1/` is built entirely from the comprehensive
test fixture (`tests/create_fixtures.py`). Every variant and every artifact is
*planted*, so the expected labels are known by construction — there is no
second tool to disagree with. That makes the gate:

- **deterministic** — same inputs, same metrics, every run;
- **fast & hermetic** — no downloads, no API keys, runs in the normal test job;
- **meaningful** — it measures biological correctness (did we call the right
  variant, flag the right artifact), not just "the tool returned 200".

The current synthetic_v1 numbers: **8/8 variant alleles recovered (P=R=F1=1.0),
0 false positives across 3 negative regions, 3/3 artifact types surfaced,
clean control correctly discriminated.**

### Real ground truth (GIAB)

`tests/eval/datasets/giab/` is the extension to the
[Genome in a Bottle](https://www.nist.gov/programs-projects/genome-bottle)
HG001 high-confidence call set. It uses the **same manifest schema**, so the
scorer runs unchanged; only the provenance and floors differ (recall ≥ 0.95,
precision ≥ 0.90 — BAMCP is a lightweight pileup caller, not DeepVariant). GIAB
data is large and network-gated, so it is a manual / nightly job, not a PR
gate. See `tests/eval/datasets/giab/README.md`.

### Adding a truth set

A manifest is a single YAML file:

```yaml
dataset: my_set
version: 1
detection_region: chr1:1000-2600
provenance:
  bam: path/to/reads.bam
  reference: path/to/ref.fa
variant_sites:
  - {chrom: chr1, pos: 1050, ref: A, alt: T, expected_artifact: strand_bias}
  - {chrom: chr1, pos: 1150, ref: A, alt: G, is_clean_control: true}
negative_regions:
  - {region: chr1:2100-2150}
```

Then: `python -m bamcp.eval.truthset --manifest path/to/manifest.yaml`.

---

## 2. LLM harness

The runner (`bamcp.eval.runner`) drives one model through the tools per case,
captures tool-call telemetry, and grades two ways:

- **Deterministic** when a case declares `tools_expected` or inline rubric
  thresholds (`score_at_least: depth_quality=0.5`).
- **LLM judge** for free-form `expected` criteria, falling back to a small
  judge prompt.

Cases live in `tests/eval/test_cases.yaml` (text) and
`tests/eval/vision_cases.yaml` (vision-required). The harness is
MARRVEL-MCP-compatible so the two can be compared head-to-head.

---

## Commands

```bash
# Deterministic gate (CI runs this; no API key needed)
make eval-smoke                       # CLI report + nonzero exit on regression
pytest tests/eval/test_truthset_smoke.py   # same gate as pytest assertions

# LLM harness (needs ANTHROPIC_API_KEY or OPENAI_API_KEY + [eval] extras)
make eval                             # run text cases
make eval-compare                     # add vanilla + rendering-mode comparison
make eval-dry                         # list cases without calling a model
```

---

## Roadmap

The ground-truth gate is the foundation. Planned build-out, in priority order:

1. **GIAB wiring** — implement `fetch_giab.py`'s slice + VCF parse so a real
   truth set runs nightly.
2. **Judge validation** — human-label a 20-case subset and report
   judge–human concordance; add position-swap bias checks.
3. **Model comparison** — run the harness across Claude and OpenAI models and
   publish an accuracy-by-model table.
4. **Tool-vs-no-tool ablation** — quantify BAMCP's value-add by comparing a
   model with the tools against the same model answering from raw text.
