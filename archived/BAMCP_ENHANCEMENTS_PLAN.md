# BAMCP Phase 1 Enhancements — Rendering, Rubric, Telemetry, and Eval Harness

## Context

Recent landscape research (MARRVEL-MCP, MCP-Atlas/AgentBench, DeepVariant, AlphaMissense literature) identified six high-value enhancements for BAMCP. After scope discussion we are committing to a Phase 1 that delivers the substrate the other features will be evaluated on:

1. **DeepVariant-style pileup rendering** (two sub-modes: strips + composite) — gives us the comparative-eval substrate (IGV vs DV-strips vs DV-composite).
2. **Structured rubric output for `get_variant_curation_summary`** — exposes the structured dict already computed internally so eval harnesses can score curation deterministically.
3. **Tool-use telemetry (OTel + JSONL fallback)** — captures every tool call in a format the eval harness can consume.
4. **On-command eval harness** — compatible with MARRVEL-MCP's eval schema so we can run head-to-head benchmarks. `make eval` runs it, results land in `~/.cache/bamcp/evaluations/<run_id>/`.
5. **Full coverage discipline** — every new module ships with unit + integration tests; the existing 80% `fail_under` is preserved and applied to the new code paths.

Phase 2 (deferred): SpliceAI/AlphaMissense annotation clients, ACMG criterion scaffolding, MARRVEL-MCP interop manifest endpoint. Phase 1 is the prerequisite for Phase 2 evaluation.

Intended outcome: a single shippable increment where (a) the viewer can be flipped between IGV / DV-strips / DV-composite, (b) the curation tool produces machine-scorable rubric output, (c) every tool call is traced, (d) `make eval` runs a YAML benchmark against the live server and prints a head-to-head comparison report, and (e) coverage stays at or above 80% across new code.

---

## Feature 1 — DeepVariant-Style Pileup Rendering

Two new viewer modes, both selectable from the existing display-mode dropdown:

- **`dv-strips`** — each read renders as six stacked horizontal strips (base, base-quality, MAPQ, strand, supports-variant, base-differs-from-ref). Readable for humans.
- **`dv-composite`** — each read renders as a single row whose pixel color encodes all 6 channels (RGB-packed). Faithful to DV's `(height, width, 6)` tensor for the comparative eval.

### Backend changes

**[src/bamcp/core/serialization.py](src/bamcp/core/serialization.py)** — `_serialize_read()` currently omits `qualities`. The `AlignedRead.qualities` list is already extracted in [src/bamcp/core/parsers.py:281](src/bamcp/core/parsers.py#L281).
- Add `qualities: list[int]` to the serialized read, gated on `compact=False`.
- Payload impact: ~1 byte/base; bounded by the existing 500bp compact threshold.

No changes to `parsers.py` — data is already extracted.

### Viewer changes (TypeScript)

**[src/bamcp/static/types.ts](src/bamcp/static/types.ts)** — extend the `ReadDisplayMode` union to `'squished' | 'compact' | 'expanded' | 'dv-strips' | 'dv-composite'`. Add `qualities?: number[]` to the read type.

**[src/bamcp/static/constants.ts](src/bamcp/static/constants.ts)** — extend `DISPLAY_MODE_CONFIGS`:
- `dv-strips`: `{ readHeight: 36, readGap: 4 }`
- `dv-composite`: `{ readHeight: 8, readGap: 1 }`

Add two new palettes: `variantSupport` (gray → ref / alt color), `refDiffers` (black → cyan).

**[src/bamcp/static/viewer.html](src/bamcp/static/viewer.html)** — add two `<option>` entries to `#display-mode-select`.

**[src/bamcp/static/renderer.ts](src/bamcp/static/renderer.ts)** — split `renderReads()` (currently monolithic, lines 471–699). Branch by mode and add:
- `renderDeepVariantStrips(ctx, read, x, y, scale, activeVariantPos)` — six per-base strips.
- `renderDeepVariantComposite(ctx, read, x, y, scale, activeVariantPos)` — packed RGB row with a single packing helper for clarity.

**[src/bamcp/static/state.ts](src/bamcp/static/state.ts)** — add `activeVariantPosition: number | null` to settings. Variant-panel clicks set this; renderer reads it for the "supports-variant" channel.

**[src/bamcp/static/mcp-app.ts](src/bamcp/static/mcp-app.ts)** — extend the existing display-mode handler (lines 403–410) to pass through the new modes (no special-case logic; dimensions come from `DISPLAY_MODE_CONFIGS`).

### Build

`cd src/bamcp/static && npm run build` (existing flow). **Never** hand-edit `dist/`.

---

## Feature 2 — Structured Rubric Output for `get_variant_curation_summary`

The structured dict is already assembled at [src/bamcp/analysis/curation.py:209-237](src/bamcp/analysis/curation.py#L209) before being passed through `format_curation_summary()` (line 240). Expose it as opt-in.

### Changes

**[src/bamcp/analysis/curation.py](src/bamcp/analysis/curation.py)** — modify `handle_get_variant_curation_summary`:
- Add `format: Literal["text", "rubric"] = "text"` kwarg.
- When `format == "rubric"`: skip `format_curation_summary()`; return structured dict as JSON in the tool content envelope.
- Add `rubric_version: "1.0"` so the eval harness can pin behavior.
- Add a `scores` sub-dict with normalized 0–1 floats for: `vaf_quality`, `depth_quality`, `strand_balance`, `mapq_quality`, `position_quality`, `artifact_risk_inverse`. All computable from existing histograms — no new evidence pipeline.

**[src/bamcp/core/tools.py](src/bamcp/core/tools.py)** — propagate the `format` arg through the handler.

**[src/bamcp/server.py](src/bamcp/server.py)** — add `format: str = "text"` to the FastMCP wrapper (lines 211–231). FastMCP infers schema from signature.

---

## Feature 3 — Tool-Use Telemetry (OTel + JSONL Fallback)

Wrap every handler in `core/tools.py` with a telemetry decorator. JSONL is always written when enabled; OTel exporter layers on when configured.

### New file

**[src/bamcp/middleware/telemetry.py](src/bamcp/middleware/telemetry.py)** — new module exporting:

```python
def telemetry_wrapper(tool_name: str):
    """Decorator for handler functions in core/tools.py.
    Captures: tool_name, args (sanitized), result_summary, duration_ms, ok, error.
    Always writes JSONL when enabled; emits an OTel span when OTel is configured."""
```

Event schema (matches `tool_invocations` table from [archived/BAMCP_EVAL_HARNESS.md:848-857](archived/BAMCP_EVAL_HARNESS.md#L848-L857)):

```json
{
  "ts": "<ISO8601>",
  "session_id": "<str>",
  "tool_name": "<str>",
  "args": { /* sanitized: BAM paths hashed */ },
  "result_summary": { /* counts/status only, ≤1 KB */ },
  "duration_ms": <float>,
  "ok": <bool>,
  "error": "<str|null>",
  "trace_id": "<otel-trace-id-or-null>",
  "span_id": "<otel-span-id-or-null>"
}
```

Sanitization helper (single function):
- `bam_path`/file paths → `sha256[:12]` + filename suffix.
- Region strings/params: raw.
- Tokens/auth headers: never included.
- `result_summary` capped at ~1 KB.

### Config

**[src/bamcp/config.py](src/bamcp/config.py)** — extend `BAMCPConfig` and `from_env`, following the existing `*_enabled` pattern (lines 52–54, 136–137):
- `telemetry_enabled: bool = False` ← `BAMCP_TELEMETRY_ENABLED`
- `telemetry_path: str = ""` ← `BAMCP_TELEMETRY_PATH` (empty disables JSONL output)
- `telemetry_otel_enabled: bool = False` ← `BAMCP_TELEMETRY_OTEL`

### Wiring

**[src/bamcp/core/tools.py](src/bamcp/core/tools.py)** — apply `@telemetry_wrapper("<tool_name>")` to each `handle_*` handler (~11 handlers, mechanical edit).

**[src/bamcp/__main__.py](src/bamcp/__main__.py)** — initialize OTel tracer at startup when `config.telemetry_otel_enabled`. Uses `opentelemetry-sdk` + `opentelemetry-exporter-otlp-proto-http`, in a `[telemetry]` extras group so the base install stays slim.

---

## Feature 4 — On-Command Eval Harness (MARRVEL-MCP-Compatible)

A scriptable benchmark runner mirroring MARRVEL-MCP's eval CLI so the two systems can be compared head-to-head on the same harness conventions.

### Schema (test cases)

**[tests/eval/test_cases.yaml](tests/eval/test_cases.yaml)** — YAML adopting MARRVEL-MCP's exact field names (`name`, `category`, `input`, `expected`) plus two BAMCP-specific optional fields:

```yaml
- case:
    name: "TP53 R248W pathogenicity"
    category: "Variant Interpretation"
    input: "Classify chr17:7577538 G>A in the provided BAM. Use get_variant_curation_summary and lookup_clinvar."
    expected: "Likely Pathogenic or Pathogenic; ClinVar entry exists with ≥2 stars."
    # BAMCP-specific (optional):
    tools_expected: ["get_variant_curation_summary", "lookup_clinvar"]
    bam_fixture: "tests/fixtures/tp53_r248w.bam"
    rendering_modes: ["expanded", "dv-strips", "dv-composite"]  # if present, the case is run once per mode for comparison
```

`tools_expected` is checked against telemetry events (Feature 3) — passes if every expected tool was called at least once. `rendering_modes` triggers the comparative experiment.

### Runner

**[scripts/run_eval.py](scripts/run_eval.py)** — new CLI mirroring MARRVEL's flags:

```
python scripts/run_eval.py [--cache] [--subset N-M] \
  [--provider anthropic|openai|openrouter] [--model MODEL_ID] \
  [--with-vanilla] [--with-rendering-comparison] \
  [--output-dir DIR]
```

- `--cache` — reuse prior successful results; only re-run failures (matches MARRVEL).
- `--subset` — run a specific case range.
- `--provider`/`--model` — provider abstraction (default `anthropic` / `claude-opus-4-7` since this is a Claude-aligned project; `openai`/`openrouter` supported through a thin adapter).
- `--with-vanilla` — also run each case against the model with no MCP tools, for delta measurement (matches MARRVEL's `--with-vanilla`).
- `--with-rendering-comparison` — runs each case with `rendering_modes` once per mode and reports per-mode accuracy.
- Output: `~/.cache/bamcp/evaluations/<run_id>/` with `results.yaml`, `results.json`, `run_config.yaml`, `test_cases.yaml`, `report.html`, and `telemetry.jsonl` (the captured tool-call trace from Feature 3).

### Grading

A two-tier grader:
1. **Deterministic checks first** — when `tools_expected` is present, validate against telemetry; when the rubric tool was called, validate the `scores` dict against thresholds in `expected`.
2. **LLM judge fallback** — for free-form `expected` strings, use the same model (by default) or an override via `--judge-model` to grade pass/fail with a fixed rubric prompt. Cache judge verdicts by `(case_name, response_hash)` to keep runs cheap.

This matches what MARRVEL-MCP appears to use (its README leaves grading unspecified beyond pass/fail status), and lets us be deterministic where possible.

### Provider abstraction

**[src/bamcp/eval/](src/bamcp/eval/)** — new package:
- `runner.py` — main loop (case → LLM → tool calls → grading → result).
- `providers.py` — `AnthropicProvider`, `OpenAIProvider`, `OpenRouterProvider` (thin wrappers around their SDKs; each implements `chat_with_tools(messages, tools) -> response`).
- `grader.py` — deterministic + LLM-judge graders.
- `report.py` — HTML/JSON/YAML report writer matching MARRVEL's output layout.
- `comparison.py` — head-to-head report when both `--with-vanilla` and tool runs are present, plus per-rendering-mode delta when `--with-rendering-comparison` is set.

The runner talks to the running BAMCP server over MCP (stdio by default for local; configurable to HTTP). It does **not** invoke handlers directly — it goes through the MCP protocol so we're testing the integration the LLM actually sees.

### Dependencies

Add to `pyproject.toml`:
```toml
[project.optional-dependencies]
eval = [
    "anthropic>=0.40",
    "openai>=1.50",
    "pyyaml>=6",
    "jinja2>=3.1",       # HTML report
]
telemetry = [
    "opentelemetry-sdk>=1.27",
    "opentelemetry-exporter-otlp-proto-http>=1.27",
]
```

Install: `pip install -e ".[eval]"` (mirrors MARRVEL's `uv sync --extra eval`).

### Makefile

**[Makefile](Makefile)** — add:
```make
eval: build-viewer
    python scripts/run_eval.py --output-dir .eval-results

eval-cached:
    python scripts/run_eval.py --cache --output-dir .eval-results

eval-compare:
    python scripts/run_eval.py --with-vanilla --with-rendering-comparison --output-dir .eval-results

eval-marrvel-compat:
    # Sanity check that our schema parses MARRVEL test cases too
    python scripts/run_eval.py --test-cases external/marrvel_test_cases.yaml --dry-run
```

---

## Feature 5 — Full Coverage Discipline

Existing setup: `[tool.coverage.report] fail_under = 80` in [pyproject.toml:89](pyproject.toml#L89), `make test` runs everything except e2e, `make coverage` produces HTML.

This phase preserves that floor and *applies it to all new code* (which means real tests, not coverage-shy stubs).

### Required tests per feature

| Feature | Test files (new or extended) |
|---|---|
| Serialization adds `qualities` | `tests/unit/core/test_serialization.py` (extend) — qualities present when `compact=False`, absent when `True` |
| DV rendering modes | `tests/e2e/test_viewer_e2e.py` (extend) — switch to each new mode, assert canvas renders, assert variant-click sets active position; use existing Playwright pattern |
| Rubric mode | `tests/unit/analysis/test_curation_rubric.py` (new) — `format="rubric"` returns JSON with `scores` and `rubric_version`; `format="text"` unchanged; `scores` values bounded in [0,1] |
| Telemetry decorator | `tests/unit/middleware/test_telemetry.py` (new) — captures args/duration/error; JSONL output shape; sanitizer hashes BAM paths; OTel span emitted when configured (use in-memory exporter); decorator is a no-op when disabled |
| Telemetry config | `tests/unit/core/test_config.py` (extend) — three new env vars parsed correctly |
| Eval runner | `tests/unit/eval/test_runner.py`, `test_grader.py`, `test_providers.py` (new) — provider abstraction mocked; deterministic grader on telemetry; LLM-judge with stubbed model; subset/cache flags honored |
| Eval CLI smoke | `tests/integration/test_eval_cli.py` (new) — runs `scripts/run_eval.py --subset 1-1 --provider mock` end-to-end against a stub BAM fixture, asserts `results.yaml` exists and matches schema |

### Coverage gates

- Keep `fail_under = 80` global.
- Add `[tool.coverage.run] branch = true` to catch unhit branches in the new branching code (rendering mode switch, format param, telemetry on/off).
- Add `make coverage-strict` target that runs `pytest --cov=bamcp --cov-fail-under=85` so reviewers can verify the new modules don't dilute coverage:
```make
coverage-strict:
    python -m pytest tests/ --ignore=tests/e2e --cov=bamcp --cov-fail-under=85 --cov-report=term-missing
```
- CI (`.github/workflows/ci.yml`) — extend the existing test job to also run `make coverage-strict` on PRs touching `src/bamcp/eval/` or `src/bamcp/middleware/telemetry.py`.

### Test fixtures

Extend `tests/create_fixtures.py` to generate at least one small BAM with a known TP53 variant (for eval smoke + DV rendering tests). Keep size small (<100 reads).

---

## Critical Files Summary

| File | Change |
|------|--------|
| [src/bamcp/core/serialization.py](src/bamcp/core/serialization.py) | Add `qualities` field |
| [src/bamcp/static/types.ts](src/bamcp/static/types.ts) | Extend `ReadDisplayMode`; add `qualities?` to read type |
| [src/bamcp/static/constants.ts](src/bamcp/static/constants.ts) | DV mode configs + 2 new palettes |
| [src/bamcp/static/viewer.html](src/bamcp/static/viewer.html) | 2 dropdown options |
| [src/bamcp/static/renderer.ts](src/bamcp/static/renderer.ts) | Split `renderReads`; add DV strip + composite functions |
| [src/bamcp/static/state.ts](src/bamcp/static/state.ts) | `activeVariantPosition` |
| [src/bamcp/static/mcp-app.ts](src/bamcp/static/mcp-app.ts) | Wire variant click → active position |
| [src/bamcp/analysis/curation.py](src/bamcp/analysis/curation.py) | `format` param + `scores` dict |
| [src/bamcp/core/tools.py](src/bamcp/core/tools.py) | `@telemetry_wrapper` on all handlers; pass `format` |
| [src/bamcp/server.py](src/bamcp/server.py) | Add `format` to curation tool |
| [src/bamcp/middleware/telemetry.py](src/bamcp/middleware/telemetry.py) | **NEW** decorator + JSONL + OTel hook |
| [src/bamcp/config.py](src/bamcp/config.py) | 3 telemetry env vars |
| [src/bamcp/__main__.py](src/bamcp/__main__.py) | OTel tracer init |
| [src/bamcp/eval/](src/bamcp/eval/) | **NEW** package — `runner.py`, `providers.py`, `grader.py`, `report.py`, `comparison.py` |
| [scripts/run_eval.py](scripts/run_eval.py) | **NEW** CLI entry |
| [tests/eval/test_cases.yaml](tests/eval/test_cases.yaml) | **NEW** benchmark cases (MARRVEL-compatible schema) |
| [tests/create_fixtures.py](tests/create_fixtures.py) | Add TP53 variant fixture |
| [pyproject.toml](pyproject.toml) | `[eval]` + `[telemetry]` extras; `branch = true` |
| [Makefile](Makefile) | `eval`, `eval-cached`, `eval-compare`, `coverage-strict` targets |
| [.github/workflows/ci.yml](.github/workflows/ci.yml) | Run `make coverage-strict` on touched paths |

## Reused Patterns / Utilities

- **Display-mode dropdown precedent** — [viewer.html:602-606](src/bamcp/static/viewer.html#L602-L606) + [DISPLAY_MODE_CONFIGS](src/bamcp/static/constants.ts#L46-L50). New modes drop in.
- **Structured-dict-before-text pattern** — already at [curation.py:209-237](src/bamcp/analysis/curation.py#L209-L237). Rubric mode is "skip the formatter".
- **Config flag pattern** — `BAMCP_<FEATURE>_ENABLED` at [config.py:52-54, 136-137](src/bamcp/config.py#L52). Telemetry follows.
- **Handler signature** — `async def handle_<x>(args, config) -> dict` at [tools.py:197-217](src/bamcp/core/tools.py#L197). Telemetry decorator wraps this exact signature.
- **Async-safe TTL cache** — [clients/ttl_cache.py](src/bamcp/clients/ttl_cache.py) for any judge-verdict caching in eval grader.
- **OAuth-secured HTTP transport** already in [middleware/auth.py](src/bamcp/middleware/auth.py) — eval runner can target either stdio or HTTPS transports.

## Verification

End-to-end smoke (the full Phase 1 acceptance sequence):

```bash
# 1. Backend unit tests + coverage floor (80%)
make test

# 2. Strict coverage on new code (85%)
make coverage-strict

# 3. Lint + types
make lint && make typecheck

# 4. Rebuild viewer
cd src/bamcp/static && npm run build && cd -

# 5. E2E (Playwright) — exercises both new DV modes
make test-e2e

# 6. Manual: rubric mode through MCP inspector
#    call get_variant_curation_summary with format="rubric" → expect JSON with scores{} + rubric_version

# 7. Telemetry smoke
BAMCP_TELEMETRY_ENABLED=true BAMCP_TELEMETRY_PATH=/tmp/bamcp.jsonl \
  python -m bamcp <<< '<one visualize_region jsonrpc call>'
tail /tmp/bamcp.jsonl    # one event per tool call

# 8. Eval harness — minimal smoke
pip install -e ".[eval]"
make eval                 # uses default Anthropic provider; needs ANTHROPIC_API_KEY

# 9. Eval harness — full comparison run (rendering experiment)
make eval-compare         # produces report.html showing IGV vs DV-strips vs DV-composite accuracy

# 10. MARRVEL schema compatibility sanity check
make eval-marrvel-compat  # parses MARRVEL's test_cases.yaml without errors

# 11. Optional: OTel
OTEL_EXPORTER_OTLP_ENDPOINT=http://localhost:4318 \
BAMCP_TELEMETRY_OTEL=true BAMCP_TELEMETRY_ENABLED=true \
  python -m bamcp        # spans appear in local OTel collector
```

## Out of Scope (Phase 2)

- SpliceAI / AlphaMissense annotation clients (follow [clients/clinvar.py](src/bamcp/clients/clinvar.py) pattern).
- ACMG criterion scaffolding (will read evidence dict + external lookups).
- MARRVEL-MCP interop / `/manifest.json` discovery endpoint.
- Auto-mode selection for DV in `visualize_region` (Phase 1 is manual toggle only).
- Continuous-improvement loop from [archived/BAMCP_EVAL_HARNESS.md:591-654](archived/BAMCP_EVAL_HARNESS.md#L591-L654) (scaffold corrections fed back from eval data). Phase 1 produces the data; the loop is Phase 2.
