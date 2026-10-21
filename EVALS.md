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
| Artifact-type recall | Does the curation tool surface the *expected* artifact at each known-artifact site (strand bias, read-position bias, low MAPQ, homopolymer, low depth)? | recall = 1.0 |
| Artifact likelihood | Does each labeled site land on the right side of its `expected_likelihood` (artifact sites: at least it; non-artifact sites: at most it)? | 0 mismatches |
| VAF | Is each detected allele's VAF within `vaf_tolerance` (0.02) of `expected_vaf`? | 0 mismatches |
| Tool failures | Did every tool call succeed? A failed or `{"error": ...}` call never counts as "0 calls" / "risk 0". | 0 failures |
| Clean-control discrimination | Is a true clean variant scored *below* the artifact-prone sites and not labeled high-risk? | true |
| **Safety: overconfidence guard** | Is any high-artifact site *also* reported as high confidence? | 0 sites ([SAFETY.md](SAFETY.md)) |
| Confidence positive control | Does the one site engineered to be clean *actually* reach high confidence? | 0 mismatches |

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

The current synthetic_v1 numbers: **9/9 variant alleles recovered (P=R=F1=1.0),
0 false positives across 3 negative regions, 5/5 artifact types surfaced,
every likelihood and VAF label within bounds,
clean control correctly discriminated, 0 overconfident sites, and the
high-confidence positive control reaches high confidence.**

The positive control matters: a guard that says "no high-artifact site is
high-confidence" is vacuous if the tool can *never* report high confidence. The
chr1:2800 site (120bp reads, centered variant, balanced strands, high MAPQ, no
artifact context) is engineered to be the one site that legitimately earns it —
so the guard is discriminating, not just trivially satisfied.

### Real ground truth (GIAB)

`tests/eval/datasets/giab/` runs the same scorer against **real NIST
Genome-in-a-Bottle truth data** (HG001, v4.2.1, GRCh38). `fetch_giab.py`
streams the truth VCF, downloads the GRCh38 reference, and simulates reads over
it honoring the real genotypes — a semi-synthetic benchmark over real biology.
It uses the same manifest schema, so `bamcp.eval.truthset` scores it unchanged.

Headline result on a 60 kb chr20 slice (66 real NIST SNVs): **recall 1.000 at
both operating points; precision goes from 0.034 at BAMCP's default
curation-sensitivity thresholds to 1.000 when thresholded for calling** (VAF ≥
0.20, depth ≥ 8, MAPQ ≥ 20). This is **one 60 kb slice of simulated reads over real
truth genotypes** — evidence that the detector recovers planted SNVs on this slice, not a
general claim of caller correctness across genomes, indels, SVs, or real sequencer error.
That precision/recall tradeoff — and an honest statement of what a simulated-read benchmark
can and cannot show — is in
[`tests/eval/datasets/giab/GIAB_RESULTS.md`](tests/eval/datasets/giab/GIAB_RESULTS.md).
GIAB data is large and network-gated, so this is a manual / nightly job, not a
PR gate.

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
  judge prompt. The judge is **fixed and separate from the model under test**
  (default `anthropic` / `claude-haiku-5-5`, override with `--judge-provider` /
  `--judge-model`; a judge equal to the subject is refused) and runs at
  **temperature 0**, so every configuration is graded by the same deterministic
  judge. Verdict parsing tolerates markdown (`**VERDICT: PASS**`); no
  recognizable verdict is a fail.
- **Tool-use ablation** (`bamcp.eval.compare`) grades both arms on answer
  correctness only — the tools-off arm cannot call tools, so grading it on
  `tools_expected` would fail it by construction.

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

# Model comparison + tool-use ablation (EVAL_RESULTS.md)
make eval-matrix                      # real: accuracy-by-model + tool lift (needs key)
make eval-matrix-mock                 # offline: prove the pipeline with the mock provider

# Publish results to Tracer (best-effort; needs TRACER_INGEST_TOKEN to --post)
python -m bamcp.eval.truthset --manifest tests/eval/datasets/synthetic_v1/manifest.yaml --json > syn.json
python -m bamcp.eval.tracer truthset syn.json --out payload.json   # inspect the payload
python -m bamcp.eval.tracer llm .eval-results .eval-results/comparison --post
```

---

## Trending in Tracer

Results are published to [Tracer](https://tracer.rtrentjones.dev/projects/bamcp), the project
signals dashboard, by `bamcp.eval.tracer`:

- **Ground truth** (`mode=truthset`, one run per manifest, `suite` = dataset). The cases are the
  gate's floors, evaluated at the thresholds the report records under `floors` (so a floor case can
  never disagree with the gate), plus each artifact-type check and every missing or spurious call.
  Pass rate is the fraction of floors met; precision, recall, F1, TP/FP/FN and artifact recall are
  run metrics that Tracer charts per run.
- **LLM harness** (`mode=llm-eval`): one run per harness directory, with each case's input and
  expected answer joined from `test_cases.yaml`, the judge's rationale, and per-category accuracy as
  metrics. A comparison directory publishes one run per spec (model × tools/vision arm, `suite =
  comparison`), so Tracer's compare view sets the arms side by side.

Only the ground truth gates BAMCP's health in Tracer. LLM harness and comparison runs are published
with `gating: false` ("trended, not gated"): Tracer charts their pass rate and trend and labels the
signal *tracked*, but a model flicker never turns the project, badge or attention list red. Case
names longer than Tracer's 255-character limit (a long spurious indel, say) are clipped with `…`;
the full call stays in the case's output or expected value.

CI publishes from `ci.yml` (every push to `main`, ground truth only) and `eval-nightly.yml` (ground
truth nightly; the LLM harness and comparison when they run). Publishing never changes a job's
result: without the `TRACER_INGEST_TOKEN` repository secret it skips, and a failed POST is a
warning annotation. Commit, branch, CI run URL and trigger come from the Actions environment.

**Enabling it.** Set the repository secret once, with the same value as the wrapper repo's
`TF_VAR_TRACER_INGEST_TOKEN` (the token Tracer is deployed with):

```bash
gh secret set TRACER_INGEST_TOKEN -R RTrentJones/BAMCP
```

**Confirming it worked.** Each publish step writes its outcome to the job summary: `✅ Tracer:
published truthset/synthetic_v1/deterministic (201)` on success, or a `⚠️ Tracer publish failed` /
`⚠️ Tracer publish skipped: TRACER_INGEST_TOKEN unset` line otherwise. After the next push to
`main` (or a manual `eval-nightly` dispatch), the runs appear at
[tracer.rtrentjones.dev/projects/bamcp](https://tracer.rtrentjones.dev/projects/bamcp).

---

## Roadmap

The ground-truth gate is the foundation. Planned build-out, in priority order:

1. **Real reads for GIAB** — the GIAB benchmark now runs on real truth with
   simulated reads (done); swap in a real HG001 alignment (sliced remotely via
   BAMCP's remote-BAM support) to also capture real error/mapping behavior.
2. **Judge validation** — human-label a 20-case subset and report
   judge–human concordance; add position-swap bias checks.
3. **Model comparison + tool-use ablation** — *infrastructure built and
   validated offline* (`bamcp.eval.compare`, `make eval-matrix`): runs each
   model with and without the tools, aggregates accuracy by config and
   category, and emits a Markdown/JSON table. Real numbers pending an API key —
   a real run is a drop-in provider swap. See
   [EVAL_RESULTS.md](EVAL_RESULTS.md). Follow-up: add judge-graded reasoning
   cases so the ablation measures answer-quality lift, not just tool presence.
