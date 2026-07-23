# BAMCP model comparison & tool-use ablation

Two questions a deployment actually asks, answered by re-running the eval
harness across configurations (`bamcp.eval.compare`):

1. **Which model** best drives the BAMCP tools?
2. **Are the tools worth it** — how much accuracy does giving the model BAMCP
   buy over the same model answering from its own knowledge?

The harness runs each model **with** the tools and **without** them (tools
withheld, same cases) and reports the lift.

---

## Status

- **Pipeline: done and validated.** The comparison runner, aggregation, and the
  Markdown/JSON tables are implemented and proven end-to-end with the `mock`
  provider (see the validation table below and `tests/unit/eval/test_compare.py`).
- **Real numbers: pending an API key.** A real run is a drop-in swap of the
  provider/model — no code changes:

  ```bash
  export ANTHROPIC_API_KEY=...        # or OPENAI_API_KEY for --provider openai
  make eval-matrix                     # defaults to claude-opus-4-8 + claude-sonnet-4-5
  # or, explicitly:
  python scripts/run_comparison.py --provider anthropic \
    --models claude-opus-4-8 claude-sonnet-4-5 \
    --test-cases tests/eval/test_cases.yaml \
    --results-md EVAL_RESULTS.md
  ```

  Paste the generated table into the [Results](#results) section below.

---

## Results

> _Placeholder — populated by the first real run. The tables below are the
> **mock-provider validation**, which proves the pipeline produces a
> well-formed, differentiated report without a key._

### Pipeline validation (`--provider mock`)

The `mock` model is scripted to behave like a capable agent: with tools
advertised it parses the prompt and calls the right BAMCP tools; with tools
withheld it can only answer from text. On the 4 synthetic cases:

| Configuration | Tools | Passed | Cases | Accuracy |
|---|---|---:|---:|---:|
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
the model to call a specific tool and are graded on whether it did. So the
"without tools" arm scores ~0 almost by construction: no tools, no
tool-verified answer. That validates the *plumbing* but is not yet a fair
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
  (`--results-md`). Nightly CI runs this when `ANTHROPIC_API_KEY` is present and
  uploads both as artifacts.
