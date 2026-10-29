# BAMCP model comparison & tool-use ablation harness

> **No real-model results are committed here.** This file documents the comparison harness and
> its mock-provider validation. BAMCP's committed evaluation results are elsewhere:
> the CI ground-truth gate ([EVALS.md](EVALS.md), truth sets in `tests/eval/datasets/`), the
> [GIAB benchmark](tests/eval/datasets/giab/GIAB_RESULTS.md), and the
> [visual-vs-text study](docs/eval/visual-results.md).

Two questions a deployment actually asks, answered by re-running the eval
harness across configurations (`bamcp.eval.compare`):

1. **Which model** best drives the BAMCP tools?
2. **Are the tools worth it** — how much accuracy does giving the model BAMCP
   buy over the same model answering from its own knowledge?

The harness runs each model **with** the tools and **without** them (tools
withheld, same cases) and reports the lift.

---

## Status

- **Pipeline: implemented and runs end-to-end offline.** The comparison runner,
  aggregation, and the Markdown/JSON tables run with the `mock` provider (see the
  validation run below and `tests/unit/eval/test_compare.py`).
- **Real-model results: not run.** No real-model comparison results are committed to this
  repository. A real run needs an API key and no code changes:

  ```bash
  export ANTHROPIC_API_KEY=...
  make eval-matrix   # models: COMPARISON_MODELS in src/bamcp/eval/models.py
  ```

  `make eval-matrix` runs `scripts/run_comparison.py --provider anthropic --results-md
  EVAL_RESULTS.md`; pass `--models` to choose others. The OpenAI provider is gated off in the
  nightly workflow until the runner translates its message format.

---

## Results

> _No real-model results. The tables below are the **mock-provider validation**,
> which proves the pipeline produces a well-formed, differentiated report without a
> key; they say nothing about any real model._

### Pipeline validation (`--provider mock`)

The `mock` model is scripted to behave like a capable agent: with tools
advertised it parses the prompt and calls the right BAMCP tools; with tools
withheld it can only answer from text. Measured on 2026-10-29 with
`make eval-matrix-mock` on the 4 cases in `tests/eval/test_cases.yaml`:

| Configuration | Tools | Passed | Cases | Accuracy |
|---|---|---:|---:|---:|
| mock-model (tools) | on | 1 | 4 | **25%** |
| mock-model (no tools) | off | 0 | 4 | **0%** |

| Model | With tools | Without tools | Lift from tools |
|---|---:|---:|---:|
| `mock-model` | 25% | 0% | **+25 pts** |

Only the case with rubric score thresholds (Variant Curation) is graded
deterministically offline. The other three need the LLM judge, which the mock
run does not have, and calling the right tools never passes a case on its own
(`grader.py`), so they score as not passed. This confirms the wiring: accuracy
aggregates per config and per category, the ablation pairs tools-on/tools-off by
model, and the tables render. It is not a measurement of any model.

---|---|---:|---:|---:|
| mock-model (tools) | on | 4 | 4 | **100%** |
| mock-model (no tools) | off | 0 | 4 | **0%** |

| Model | With tools | Without tools | Lift from tools |
|---|---:|---:|---:|
| `mock-model` | 100% | 0% | **+100 pts** |

This confirms the wiring: accuracy aggregates per config and per category, the
ablation pairs tools-on/tools-off by model, and the tables render. Swap in a
real provider and the same tables fill with real accuracy.

---

## Reading the ablation honestly

The current `tests/eval/test_cases.yaml` cases are **tool-targeted** — they ask
the model to call a specific tool, and their answer checks (rubric thresholds,
`expected_claims`, or the judge) depend on that tool's output. The ablation
grades answers only (`grade_tools=False`), but the "without tools" arm still
scores ~0 almost by construction: no tools, no tool-derived answer. That validates the *plumbing* but is not yet a fair
measurement of reasoning lift.

A meaningful ablation needs **judge-graded reasoning cases** — questions with a
correct biological answer that a strong model *might* get right from its own
knowledge, so the tool's contribution is the delta in answer quality, not the
mere presence of a tool call. Adding ~10–15 such cases (graded by the LLM judge,
which needs a key) is the natural next step; the comparison infrastructure
already supports them with no changes.

## How it works

- `bamcp.eval.compare.run_comparison` runs each `RunSpec` (provider, model,
  tools on/off) through the existing runner + grader, then aggregates.
- The "no tools" arm uses a router that advertises zero tools, so the model
  genuinely cannot call BAMCP — the honest counterfactual.
- Output: `comparison.json` (machine-readable) and a Markdown table
  (`--results-md`). The nightly workflow (`.github/workflows/eval-nightly.yml`) runs the
  comparison only when the `ANTHROPIC_API_KEY` secret is configured (otherwise it reports the
  skip), uploads `.eval-results/` as an artifact, and publishes to Tracer when
  `TRACER_INGEST_TOKEN` is set.
