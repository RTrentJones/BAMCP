# BAMCP Implementation Plan

Full implementation schedule for MCP Apps viewer, ClinVar/gnomAD integration, and the Evaluation Harness.

**Current state**: Core server operational (FastMCP, 5 tools, OAuth 2.0, Docker, CI/CD, 192 tests passing). Documentation complete (CLAUDE.md, README.md, Strategy + Eval Harness design docs).

---

## Phase 1: MCP Apps Viewer Enhancement (Weeks 1-3)

**Goal**: Restructure the existing viewer to fully use MCP Apps SDK primitives so the LLM can proactively explain what the user is seeing.

### 1.1 Viewer SDK Integration — `src/bamcp/static/viewer.html`

- Replace raw `postMessage` bridge with `@modelcontextprotocol/ext-apps` App SDK
- Add `app.callServerTool()` calls for navigation (pan/zoom fetches new data from server without LLM inference)
- Add `app.updateModelContext()` on viewport changes — tells LLM what user is viewing (region, coverage, variants detected)
- Add `app.sendMessage()` on variant click — triggers LLM to explain the variant
- Add `app.requestDisplayMode("fullscreen")` toggle

### 1.2 App-Only Tools — `src/bamcp/server.py`, `src/bamcp/tools.py`

Add `visibility: ["app"]` tools for high-frequency viewer data fetching:

| Tool | Purpose |
|------|---------|
| `fetch_reads` | Get alignment data for a viewport region (chunked) |
| `fetch_coverage` | Get per-base coverage for a region |
| `fetch_variants` | Get detected variants in a region |
| `fetch_reference` | Get reference sequence bases |

These are called by the viewer directly via `callServerTool()`, not by the LLM.

### 1.3 LLM-Facing Tool Updates — `src/bamcp/server.py`

- Add `visualize_region` tool (replaces `browse_region` as the primary MCP App tool) with `_meta.ui.resourceUri`
- Add `get_region_summary` tool (text-only summary for LLM reasoning, no UI)
- Keep existing tools (`browse_region`, `get_variants`, `get_coverage`, `list_contigs`, `jump_to`)

### 1.4 Tests

- Update `tests/test_server.py` for new tools
- Update `tests/test_tools.py` for new handlers
- Update E2E tests for App SDK integration

### Success Criteria

- Viewer renders in Claude Desktop / claude.ai
- `updateModelContext` triggers proactive LLM explanations when user navigates
- App-only tools work for pan/zoom without LLM inference
- `sendMessage` on variant click produces LLM explanation

---

## Phase 2: External Database Clients (Weeks 2-4, parallel with Phase 1)

**Goal**: LLM can query ClinVar and gnomAD for variant annotation.

### 2.1 ClinVar Client — `src/bamcp/clinvar.py` (new)

- NCBI E-utilities API client (`esearch` + `esummary` for ClinVar)
- Rate limiting: respect 10 req/sec with API key, 3/sec without
- Response parsing: extract clinical_significance, review_status, stars, conditions, submitter_count
- Caching: in-memory LRU cache for repeated lookups in same session
- Error handling: timeout, rate limit, no result

### 2.2 gnomAD Client — `src/bamcp/gnomad.py` (new)

- gnomAD GraphQL API client
- Query: variant frequency by population, homozygote count, filter status
- Response parsing: global_af, max_pop_af, population breakdown
- Caching: same pattern as ClinVar

### 2.3 Tool Registration — `src/bamcp/server.py`, `src/bamcp/tools.py`

- `lookup_clinvar(chrom, pos, ref, alt)` → ClinVar significance, review status, conditions
- `lookup_gnomad(chrom, pos, ref, alt)` → Population allele frequency data
- Tool descriptions include: "research-grade information, not for clinical use"

### 2.4 Dependencies — `pyproject.toml`

- Add `httpx>=0.27` (async HTTP client for API calls)

### 2.5 Tests

- `tests/test_clinvar.py` — Mock API responses, test parsing, test error handling
- `tests/test_gnomad.py` — Mock API responses, test parsing
- Integration tests with real API (marked `@pytest.mark.slow`, skippable)

### Success Criteria

- `lookup_clinvar` and `lookup_gnomad` return structured annotation data
- Tests pass with mocked API responses
- LLM can explain variant significance with ClinVar/gnomAD context

---

## Phase 3: Evaluation Harness Infrastructure (Weeks 4-8)

**Goal**: Build the instrumentation, storage, and evaluation tools that let BAMCP assess LLM genomic reasoning accuracy.

### 3.1 SQLite Eval Store — `src/bamcp/eval_store.py` (new)

- Schema: `eval_sessions`, `tool_invocations`, `classifications`, `failure_modes` tables
- `failure_summary` view for aggregation
- Indexes on `failure_modes(mode)`, `failure_modes(criterion)`, `classifications(gene)`, `classifications(session_id)`
- CRUD operations: `create_session()`, `log_invocation()`, `log_classification()`, `get_invocation()`, `get_classifications()`
- DB path configurable via `BAMCP_EVAL_DB` env var (default: `~/.bamcp/eval.db`)
- Auto-create schema on first use

### 3.2 Instrumentation Layer — `src/bamcp/instrumentation.py` (new)

- Middleware/decorator that wraps every tool call
- Logs: tool name, arguments, response, timestamp, duration_ms, session context
- Ground truth comparison when available
- Session tracking (session_id propagated through tool calls)

### 3.3 Evidence Assembly — `src/bamcp/evidence.py` (new)

- `assemble_evidence(chrom, pos, ref, alt)` — aggregates data from BAM, ClinVar, gnomAD
- Returns structured evidence package:
  - observation: depth, VAF, strand bias, mapping quality, base quality, near-read-end fraction
  - clinvar: significance, review status, stars, conditions
  - population_frequency: global_af, max_pop_af, population breakdown
  - gene_context: gene, transcript, consequence, domain, constraint scores
  - predictions: SIFT, PolyPhen2, CADD, REVEL

### 3.4 ACMG Scaffold — `src/bamcp/acmg.py` (new)

- `build_acmg_scaffold(evidence)` — dynamically selects applicable ACMG criteria based on variant type
- Criteria definitions: PVS1, PS1, PM1, PM2, PP3, BP4, BA1, BS1 (extensible)
- Quality assessment checks: depth, strand bias, mapping quality, base quality, read-end fraction
- `load_failure_corrections()` — loads learned corrections from JSON file
- Scaffold injection: adds `known_pitfall` fields to criteria when corrections exist from prior evaluations

### 3.5 Failure Mode Detection — `src/bamcp/failure_modes.py` (new)

- `FailureMode` enum with 13 categories:
  - Classification errors: FALSE_PATHOGENIC, FALSE_BENIGN, OVERCONFIDENT_VUS, UNDERCONFIDENT_KNOWN
  - Reasoning errors: CRITERIA_MISAPPLICATION, FREQUENCY_MISINTERPRETATION, PREDICTION_OVERRELIANCE, EVIDENCE_IGNORED, EVIDENCE_FABRICATED
  - Data quality errors: ARTIFACT_MISSED, ARTIFACT_OVERCALLED
  - Communication errors: MISLEADING_CONFIDENCE, CLINICAL_LANGUAGE
- `detect_failure_modes(predicted, ground_truth, criteria_applied, evidence)` — categorizes *why* classification was wrong
- `score_classification(predicted, ground_truth, criteria_applied)` — exact match, within-one-step, per-criterion scoring

### 3.6 Eval Tools — `src/bamcp/server.py`, `src/bamcp/tools.py`

| Tool | Description |
|------|-------------|
| `classify_variant` | Assemble evidence + ACMG scaffold for a variant. Logs invocation with ground truth (hidden from LLM). |
| `submit_classification` | Capture LLM's structured classification, score against ground truth, detect failure modes. |
| `run_evaluation` | Start batch evaluation session over benchmark tiers. Presents variants one at a time. |
| `get_eval_report` | Generate accuracy report: exact match, within-one-step, confusion matrix, failure mode breakdown, per-criterion accuracy. |
| `get_expert_context` | Retrieve curated expert context for known variants (retrieval-augmented classification). |

### 3.7 Ground Truth Dataset — `src/bamcp/data/benchmark.json` (new)

Curated from ClinVar expert-reviewed variants (>= 2 stars):

| Tier | Count | Difficulty | Purpose |
|------|-------|------------|---------|
| Tier 1 (clear) | 50+ | Easy | Well-studied variants with unambiguous classifications |
| Tier 2 (nuanced) | 100+ | Medium | Requires careful evidence weighing |
| Tier 3 (adversarial) | 30+ | Hard | Trap variants targeting known LLM failure modes |

Tier 3 traps:
- `gene_name_bias`: Benign TP53 polymorphisms (tests if LLM overcalls due to gene reputation)
- `prediction_overreliance`: VUS with strong SIFT/PolyPhen predictions but no clinical evidence
- `frequency_confusion`: Edge-case AF thresholds around PM2/BA1 boundaries
- `review_status_ignored`: Low-star ClinVar entries
- `conflicting_evidence`: Contradictory ClinVar + gnomAD data
- `artifact_variant`: Low-quality observations (3 alt reads at 8x depth in homopolymer)

Format: JSON with `chrom`, `pos`, `ref`, `alt`, `gene`, `classification`, `expected_criteria`, `difficulty`, `trap`

### 3.8 Configuration — `src/bamcp/config.py`

Add eval-related fields:
- `eval_enabled: bool` (default False)
- `eval_db_path: str` (default `~/.bamcp/eval.db`)
- `corrections_file: str | None` (path to failure_corrections.json)

### 3.9 Tests

- `tests/test_eval_store.py` — SQLite CRUD operations, schema creation
- `tests/test_evidence.py` — Evidence assembly with mock data
- `tests/test_acmg.py` — Scaffold building, correction injection
- `tests/test_failure_modes.py` — All 13 failure mode detection scenarios
- `tests/test_eval_tools.py` — classify_variant, submit_classification, run_evaluation, get_eval_report integration

### Success Criteria

- SQLite eval store persists all tool invocations and classifications
- `classify_variant` returns structured evidence + ACMG scaffold without leaking ground truth
- `submit_classification` scores against ground truth and detects failure modes
- `run_evaluation` processes full benchmark tier and produces accuracy report
- `get_eval_report` shows exact match accuracy, within-one-step accuracy, confusion matrix, failure mode breakdown
- Baseline evaluation completes against Tier 1 variants

---

## Phase 4: Evaluation Execution & Improvement Loop (Weeks 8-10)

**Goal**: Run evaluations, identify failure patterns, build corrections, measure improvement.

### 4.1 Baseline Evaluation

- Run `run_evaluation(tier="all", include_corrections=False)` against Claude
- Capture baseline accuracy: exact match %, within-one-step %, per-criterion accuracy
- Identify top failure modes from `get_eval_report`

### 4.2 Scaffold Corrections — `src/bamcp/data/failure_corrections.json` (new)

- Analyze baseline failures to build targeted corrections for the ACMG scaffold
- Each correction: trigger failure mode, occurrence count, correction text, date added
- Initial corrections for top 5 failure modes (e.g., PM2 frequency thresholds, PP3 overreliance warnings)

### 4.3 Corrected Evaluation

- Run `run_evaluation(tier="all", include_corrections=True)` — scaffold now includes pitfall warnings
- Measure delta: baseline accuracy vs corrected accuracy
- Per-criterion improvement tracking
- The accuracy delta IS the publishable finding

### 4.4 Cross-Model Evaluation

- Run same benchmark against GPT-4o, Gemini (via their MCP clients)
- Compare accuracy across models on identical evidence packages
- Identify model-specific failure patterns

### Success Criteria

- Baseline accuracy numbers documented
- Correction-based improvement demonstrated (accuracy delta > 0)
- Cross-model comparison data for at least 2 models
- `improvement_over_time` metric shows trend

---

## Timeline Summary

| Phase | Work | Weeks | Depends On |
|-------|------|-------|------------|
| **1** | MCP Apps viewer enhancement | 1-3 | — |
| **2** | ClinVar + gnomAD clients | 2-4 | — (parallel with Phase 1) |
| **3** | Evaluation harness infrastructure | 4-8 | Phase 2 (needs API clients for evidence assembly) |
| **4** | Run evaluations + improvement loop | 8-10 | Phase 3 |

**Phases 1 and 2 run in parallel.** Phase 3 depends on Phase 2. Phase 4 depends on Phase 3.

**Total: ~10 weeks** from current state to publishable evaluation results.

---

## New Files Summary

| File | Phase | Purpose |
|------|-------|---------|
| `src/bamcp/clinvar.py` | 2 | ClinVar NCBI E-utilities API client |
| `src/bamcp/gnomad.py` | 2 | gnomAD GraphQL API client |
| `src/bamcp/eval_store.py` | 3 | SQLite evaluation data store |
| `src/bamcp/instrumentation.py` | 3 | Tool call logging middleware |
| `src/bamcp/evidence.py` | 3 | Evidence assembly for variant classification |
| `src/bamcp/acmg.py` | 3 | ACMG scaffold builder with correction injection |
| `src/bamcp/failure_modes.py` | 3 | Failure mode enum + detection logic |
| `src/bamcp/data/benchmark.json` | 3 | Ground truth variant benchmark dataset |
| `src/bamcp/data/failure_corrections.json` | 4 | Learned corrections from evaluation |
| `tests/test_clinvar.py` | 2 | ClinVar client tests |
| `tests/test_gnomad.py` | 2 | gnomAD client tests |
| `tests/test_eval_store.py` | 3 | Eval store CRUD tests |
| `tests/test_failure_modes.py` | 3 | Failure mode detection tests |
| `tests/test_eval_tools.py` | 3 | Eval tool integration tests |

## Modified Files Summary

| File | Phases | Change |
|------|--------|--------|
| `src/bamcp/server.py` | 1, 2, 3 | Register new tools (app-only, ClinVar, gnomAD, eval) |
| `src/bamcp/tools.py` | 1, 2, 3 | Add handlers for all new tools |
| `src/bamcp/config.py` | 3 | Add eval config fields |
| `src/bamcp/static/viewer.html` | 1 | App SDK integration |
| `pyproject.toml` | 2 | Add `httpx` dependency |

---

## Overall Success Criteria

1. **MCP Apps**: Viewer uses App SDK, `updateModelContext` triggers LLM explanations, fullscreen works
2. **API Clients**: `lookup_clinvar` and `lookup_gnomad` return structured data, tests pass with mocked APIs
3. **Eval Harness**: Full pipeline works: `classify_variant` → LLM reasons → `submit_classification` → scored + failure modes detected → `get_eval_report` produces accuracy summary
4. **Ground Truth**: 150+ curated variants across 3 tiers, including adversarial trap variants
5. **Measurable Improvement**: Scaffold corrections improve accuracy by measurable amount on repeated evaluation
6. **All tests pass**: `make test && make lint && make typecheck` green throughout
