# BAMCP Roadmap (living)

A living strengths/weaknesses assessment + prioritized focus, refreshed after the engineering-review
remediation ([`engineering-review-response.md`](engineering-review-response.md)) and 17 rounds of
adversarial review. Update the checkboxes and the "state" line as work lands.

**State (`main`, refreshed 2026-10-16):** ~1,470 Python tests (94% coverage) + 68 vitest + 154 e2e. Previously (2026-07-25): ~820 Python tests + viewer vitest / ~90% coverage ·
~9.6k lines Python · 3.6k lines viewer TS · deployed to prod, verify-gated, healthy. The external
review scored it **7.5/10** (security + LLM-evals both **4.5**, the two weakest axes). A later
hands-on **peer review** then caught a class the 17 review rounds missed — a silent genome-build
mismatch — now closed (see "Biological correctness" below).

## Newcomer usability + visual-eval round (2026-10) — landed

A multi-agent round (PM live newcomer audit → principal-engineer eval design → TPM work packages →
six parallel implementation packages → QA acceptance against real public BAMs → blocker fixes):

- [x] **Works on a bare BAM URL.** Prod had no reference, so every newcomer saw "SNV detection
  disabled". The reference is now auto-selected from the BAM's detected build and contig style
  (`reference_source: auto`), for every tool including curation and classify.
- [x] **Demo path.** `demo:HG002` / `demo:HG00096` aliases with curated loci, `bamcp://demos`,
  `bamcp://glossary`, `bamcp://guide` resources, `explore_demo` / `check_variant` /
  `explain_terms` prompts, README "Try this in 2 minutes".
- [x] **Signal over noise.** Defaults min_mapq 20, min_baseq 13, min_vaf 0.10, min_alt_count 3;
  inline per-variant `artifact_risk`, `confidence`, `concerns`, `interpretation`;
  `low_base_quality` concern; requested-site probe (`filter_check`) so a filtered artifact is
  still explained rather than "not found".
- [x] **Fits an LLM context.** `list_contigs` 157 KB → <2 KB; viewer text <1.5 KB (payload stays
  in structuredContent for the UI); contig-style mapping (chr13 ↔ 13); single-base regions.
- [x] **Cold-start + long scans.** Separate cold budget, locally cached remote `.fai`/`.gzi`,
  demo prewarm, deadline-bounded resumable `scan_variants` (`next_region`).
- [x] **Viewer.** CIGAR-aware drawing (bases were shifted after soft clips/indels; indels were
  never drawn), strand/MAPQ cues at base zoom, legend, no-reference banner, hidden-reads count.
- [x] **Eval credibility.** Nightly LLM skip is now loud; model IDs centralized; vision ablation
  wired (`--vision`); pre-registered visual-vs-text eval
  ([`eval/visual-preregistration.md`](eval/visual-preregistration.md)) with leak-free
  synthetic fixtures, leak-lint, deterministic grader, `make eval-visual-dry` CI gate.
- [x] **Visual-vs-text eval run** ([`eval/visual-results.md`](eval/visual-results.md)): the
  pre-registered contrast (image + text vs text) is not significant at the ceiling (+0.038,
  p = 0.375). Image-only failures traced to rendering; the legibility fix (neutral read bodies,
  colour only on mismatches, hollow low-BQ tiles) raised image-only accuracy 0.82 → 0.91
  (exploratory, same items).
- [x] **Round 3, pre-registered confirmation** (fresh synthetic + real GIAB HG002 tier, 756 runs):
  all three confirmatory contrasts null after Holm. Raw read text beats the picture alone on
  both tiers; the picture alone ≈ structured JSON; adding the image to raw text changes nothing.
  The model keeps getting text evidence; the viewer is for the person.
- [x] **Round 4, harder real items** (design by a researcher/ML agent panel; 2 pilots + 540
  pre-registered runs): hard GIAB classes stay at 0.94–0.96 for text; only long-repeat indel
  genotyping is hard (0.84). There, BAMCP text and viewer were both inconclusive (−0.02/−0.03)
  and cost the small model sensitivity; the viewer doubled correct zygosity calls for the mid
  and large models (secondary, p = 0.004).
- [x] **From round 4:** the homopolymer concern is read-based for indels: with ≥ 8 reads
  spanning the run it is raised only when the called length is rare among them (< 20%, or no
  commoner than other non-reference lengths); the run-length rule remains the fallback. The
  curation summary and variant interpretations report the allele fraction among spanning reads
  and take the zygosity hint from it (mpileup-style counts made a het insertion look ~30%).
  Exploratory re-measure on the round-4 items: mpileup + BAMCP text 0.839 → 0.911 (+0.072,
  CI +0.011 to +0.131), correct genotypes 13/90 → 42/90.
- [ ] **Confirm on fresh items:** a pre-registered repeat-indel round on unused main-partition
  sites, mpileup vs mpileup + BAMCP text (read-based rule).
- [ ] **Open:** viewer: mark the position of interest / mismatch columns so a hom-alt column
  cannot be read as reference (round-3 failure mode); GRCh37 demo cold latency (~27 s);
  dependabot python 3.14 / mcp<3 (lockfile + API work); real-data (GIAB) eval tier.

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
