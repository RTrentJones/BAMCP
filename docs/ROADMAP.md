# BAMCP Roadmap (living)

A living strengths/weaknesses assessment + prioritized focus, refreshed after the engineering-review
remediation ([`engineering-review-response.md`](engineering-review-response.md)) and 17 rounds of
adversarial review. Update the checkboxes and the "state" line as work lands.

**State (`main`, refreshed 2026-07-25):** ~820 Python tests + viewer vitest / ~90% coverage ·
~9.6k lines Python · 3.6k lines viewer TS · deployed to prod, verify-gated, healthy. The external
review scored it **7.5/10** (security + LLM-evals both **4.5**, the two weakest axes). A later
hands-on **peer review** then caught a class the 17 review rounds missed — a silent genome-build
mismatch — now closed (see "Biological correctness" below).

## Biological correctness (peer-review round) — closed

A peer review using the tool for real analysis found `search_gene` returned **GRCh38** coordinates
against a **GRCh37** BAM — a silent ~0.5–2 Mb offset, so the analysis was of anonymous intergenic
sequence, and the resulting empty ClinVar/gnomAD hits read as false "benign" evidence. The unifying
defect was **silent output that misleads the caller**. Fixed in four shipped phases, each with
regression tests:

- [x] **Phase 0 — build correctness** (#50): `search_gene` auto-detects the build from the BAM
  (explicit override allowed), states `genome_build`/`build_source`, and warns loudly on mismatch;
  fixed the assembly-selector that matched GRCh37's `.25` as GRCh38.
- [x] **Phase 1 — self-describing filters** (#51): every variant response echoes the `applied_filters`
  it used, so a variant's *absence* is interpretable (filtered vs. truly absent), not silent.
- [x] **Phase 2 — error vs. absence** (#52): ClinVar/gnomAD lookups carry an explicit `status`; a
  failed lookup is `unavailable`, never a null "absent" block — so a network blip can't become false
  rarity (PM2) evidence in the ACMG scaffold.
- [x] **Phase 3 — ruler rendering** (#53): position-ruler tick labels stay distinct at every zoom
  (were collapsing to `1.1K` repeated) and no longer collide with the contig label.

**The durable guardrail is a design rule, not a single test: _self-describing outputs — every
response states the configuration that shaped it_** (build + source, applied filters, lookup status).
When output carries its own provenance, a wrong assumption surfaces as a visible field instead of a
silent offset. Remaining deferred: surface the detected build in `visualize_region`/`get_region_summary`
payloads; a GRCh37-fixture biological-correctness eval (search a gene → visualize → assert reads land
on it); optional `BAMCP_GENOME_BUILD` in prod infra.

## Scorecard delta (external review → now)

| Axis | Then | Now | Note |
|---|---|---|---|
| Security posture | 4.5 | ~8 | Every P0 closed **and shipped**, with adversarial tests, verified live |
| LLM/vision evals | 4.5 | ~6.5 | Grader scores tool-selection vs answer-correctness; ablation reporting; nightly parity |
| Frontend quality | 6.5 | ~7.5 | tsc + ESLint + vitest + CI job |
| Test rigor | 7.5 | ~8.5 | Hermetic suite, lockfile-pinned CI gate |
| Architecture | 7 | ~6.5 ▼ | Regressed: `tools.py` grew to ~1,380 lines; god modules unsplit |

## Strengths (now)
- **The security boundary is real, not claimed** — the headline weakness is now the headline strength.
- **The verify-and-correct loop is evidenced** — 17 documented review rounds + `ENGINEERING_RETROSPECTIVE.md`.
- **Evals measure task success, not tool invocation** — the review's lead critique, closed at the root.
- **Deterministic-eval + safety-invariant discipline**, now behind a hermetic lockfile-pinned CI gate.

## Weaknesses (now)
- **`tools.py` is a ~1,380-line god module** (it *grew* during hardening); `parsers.py`, `mcp-app.ts`,
  `renderer.ts` unsplit.
- **Contracts still dynamically typed** — `dict[str, Any]` payload → hand-maintained TS interfaces
  (only a `schema_version` stamp + compat test so far).
- **Evals still run in-process, not through the real MCP transport** (`MCPStdioRouter` a stub); no
  human-labeled benchmark / judge-swap concordance.
- **Scientific breadth is narrow** — one 60 kb GIAB SNV slice; no indel/SV/real-read breadth.
- **Packaging/release incomplete** — no PyPI publish (`uvx bamcp` can't work), no clean-room wheel CI.
- **Coordinate semantics** convention-based (no dedicated type) — latent bug source.

## Focus next (prioritized)

### Tier 1 — depth + eval credibility (highest signal-per-effort)
- [x] **1. Port the ACMG feature** (#49) — `classify_variant` tool + evidence fusion + ACMG scoring
  eval. Shipped, and hardened in Phase 2 (error-vs-absence) so a failed lookup can't tilt a
  classification. Directly answers "scientific breadth is too narrow" with a differentiated clinical
  feature.
- [ ] **2. Finish eval credibility** (Deferred 3b/3c): implement `MCPStdioRouter` (real MCP transport)
  + a small human-labeled benchmark + judge-swap agreement. Port `tool_specs.py` (real tool schemas)
  alongside ACMG — it's a shared prerequisite.

### Tier 2 — pay down the debt that's growing
- [ ] **3. Split `tools.py`** (and the two big `.ts` files) by use case, before it grows further.
- [ ] **4. Typed payload contract** (Pydantic → JSON Schema → TS codegen) — kills the dict-boundary bug class.

### Tier 3 — release credibility
- [ ] **5. Publish to PyPI + clean-room wheel CI** (so `uvx bamcp` works, bundled viewer verified) +
  CHANGELOG/tag cadence.
- [ ] **6. Port the perf work** (read-free coverage, vectorized calling, benchmark gate) — solid, lower
  urgency; skip its caching (superseded by the MISS-sentinel work).

## Bottom line
Shifted from *"strong senior work undercut by P0 gaps"* to *"P0s closed and shipped; now flagship-polish
+ domain-depth."* The next unit of effort is best spent on **ACMG depth + finishing eval credibility**
(Tier 1) — turning "the security story is now honest" into "the scientific story is now deep."
