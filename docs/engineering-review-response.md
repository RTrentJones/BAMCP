# Engineering Review Response — Phased Remediation Plan

Response to an external engineering review of BAMCP (reviewed commit `6e91960`, scored ~7.5/10
as a hiring artifact). Every major claim below was verified against source before planning.
This is a **living checklist** — check items as they ship. Ship each item through the
Greenlight loop (see [Greenlight framing](#greenlight-framing)), not as one giant PR.

**Status:** Phases 0–1 complete (production-ready tier); Phases 2–5 advanced to bounded, tested
depth with large refactors scoped as follow-ups (see [Deferred](#deferred-follow-ups)).
**Owner:** @RTrentJones · **Created:** 2026-07-21 · Branch: `feat/eng-review-remediation`
(+ wrapper `feat/bamcp-remote-allowlist`)

---

## TL;DR — what actually matters

1. The review's major claims are **accurate** (verified in code). The two highest-value
   moves are (a) making the security claims *true* and (b) making the LLM eval measure
   *task success*, not tool invocation.
2. The single **live** production exposure is remote-fetch SSRF: prod runs with
   `BAMCP_ALLOW_REMOTE_FILES=true` and **no host allowlist**, and the reference-FASTA URL
   is completely unvalidated. Fastest fix = **populate** `BAMCP_ALLOWED_REMOTE_HOSTS` in the
   wrapper Terraform (one line) — **not** disabling remote files, which is a core feature.
3. The honest auth fix is to **delete the OAuth theater** and adopt the service-token model
   — which is *already* what the Greenlight verify contract asserts.

---

## Verified findings (checked against `6e91960`)

| Claim | Status | Evidence |
|---|---|---|
| OAuth is a token-minting service, not access control | ✅ Confirmed | `middleware/auth.py:80` `authorize()` issues a code with no user auth/consent; dynamic registration open; `required_scopes` defaults `None` |
| Rate limiter trusts `X-Forwarded-For`; unbounded dict | ✅ Confirmed | `middleware/ratelimit.py:33-37`; `_requests` is a `defaultdict(list)` cleaned only when the same IP returns |
| Remote-file validation is one-shot DNS; no redirect re-validation | ✅ Confirmed | `core/validation.py:88-100`; docstring admits the TOCTOU/rebind gap |
| **Reference-FASTA URL is unvalidated** | ✅ Confirmed | No validator for the `reference` arg; flows straight to `pysam.AlignmentFile(reference_filename=…)` — the one genuinely-open SSRF path |
| Index (`.bai`/`.crai`) URL *is* validated | ✅ Already handled | `core/tools.py:242` calls `validate_remote_url(index_url, config)` — review was wrong here; drop that item |
| Caches/registry unbounded, lazy expiry | ✅ Confirmed | `region_cache` plain dict (`tools.py:74`); `_services_registry` keyed by `id(config)`, never evicts |
| `id(config)` registry can alias stale services | ✅ Confirmed | `id()` reuse after GC → correctness bug, not just a leak |
| Grader passes on tool invocation alone | ✅ Confirmed | `eval/grader.py:48-92` — `tools_expected` satisfied + no score markers → deterministic PASS; answer text never checked |
| Eval router bypasses MCP; permissive schemas; stdio router is a stub | ✅ Confirmed | `eval/router.py:71,155,8` |
| Nightly runs only the original truth set; wires only Anthropic | ✅ Confirmed | `eval-nightly.yml` — no indel/GIAB; only `ANTHROPIC_API_KEY` despite a `provider` input |
| Frontend has no typecheck/lint/unit scripts; `uv.lock` unused by CI | ✅ Confirmed | `static/package.json` only `build`/`watch`; `ci.yml` uses `pip install -e` (broad ranges) |

## Reframes / judgment calls (where we diverge from the review)

- **"First to leverage MCP Apps"** → keep it, but **scope it**: *"the first genomics/genetic MCP
  to leverage MCP Apps"* — the verifiable form. Not a softening; a precise, checkable claim.
- **Auth**: don't integrate an IdP. **Delete** dynamic registration + the auth-code flow and
  adopt the **service-token trust model** — less code, more honest, and it's what the
  Greenlight verify contract already asserts (anon→401, service-token→authorized).
- **Remote SSRF**: "validate every redirect hop" is **not achievable while pysam/htslib does
  the fetch** (htslib re-resolves DNS + follows redirects invisibly to BAMCP). The **host
  allowlist checked pre-handoff is the real control**; residual = intra-allowlist rebinding
  (acceptable against trusted public genomics hosts). Full closure would require BAMCP to
  proxy the fetch (Option B) — deferred.
- **Streamable HTTP**: prod already runs `BAMCP_TRANSPORT=streamable-http`, so the review's
  "lead with Streamable HTTP" collapses to a **docs-only** fix.

---

## Greenlight framing

BAMCP is a **git submodule** in the `RTrentJones.dev` wrapper; deploy is **Greenlight-owned**.
The plan must respect that. Grounded in the wrapper's `infra/bamcp.tf` and
`verify/bamcp.config.ts`:

**Real prod env is Terraform, not `docker-compose`** (the compose `prod` profile is legacy):

| Prod env (OCI Container Instance) | Value | Implication |
|---|---|---|
| `BAMCP_TRANSPORT` | `streamable-http` | ✅ already non-deprecated → Phase 5c is docs-only |
| `BAMCP_AUTH_ENABLED` | `true` | on, but it's the theater (any client registers → mints a token) |
| `BAMCP_ALLOW_REMOTE_FILES` | `true` | remote fetch is live (a needed feature — see use cases) |
| `BAMCP_ALLOWED_REMOTE_HOSTS` | **unset** | ⚠️ unbounded SSRF blast radius |
| `required_scopes` | unset → `None` | no scopes enforced |
| `BAMCP_VERIFY_TOKEN` | `var.bamcp_verify_token` | the M2M service token = the Greenlight verify mechanism |

**Deploy loop (BAMCP is direct-to-prod, no beta/promote — free A1 cap, `envs=["prod"]`):**
```
greenlight preview bamcp (local container, mcp verify, auth OFF)
  → PR: ci.yml = make test + make eval-smoke            ← the gate
  → merge main → greenlight-build ship-gate (make test + eval-smoke)
  → build+push ghcr :prod / :<sha> → repository_dispatch → wrapper
  → oci-deploy-verify: OCI restart (re-pull :prod) → verify prod → commit status
```

**Prod verify contract (`verify/bamcp.config.ts`) — hard invariants any change must preserve:**
- `api` check: unauthenticated request to `/` returns **401** ("server up + OAuth enforced").
- `mcp` check: `Bearer $BAMCP_VERIFY_TOKEN` grants a working `tools/list`, `exactTools: true`
  (drift guard — a tool added/renamed/removed in code but not mirrored here **fails** the gate).
- `/__version` asserts artifact identity (v0.8).

**`make eval-smoke` is a deploy gate** (runs on every PR *and* inside the ship-gate). The
deterministic truth set must stay **hermetic, fast, green**. The stricter *LLM* grader
(answer-correctness, judge-swap) belongs in the **nightly** tier (needs API keys), never in
the blocking gate.

**Two-repo split — every phase touches both:**

| Lives in the **BAMCP submodule** (app code) | Lives in the **wrapper** (`RTrentJones.dev`) |
|---|---|
| auth/validation/ratelimit fixes, cache bounds, contracts, frontend, packaging | `infra/bamcp.tf` env changes (allowlist, remote-files, `required_scopes`) |
| adversarial auth e2e (rides `make test` → ship-gate) | `verify/bamcp.config.ts` if the auth contract or tool list changes |
| keep `make test` / `make eval-smoke` green | submodule pointer bump; `oci-deploy-verify` composite if verify semantics change |

Shared verify/deploy machinery (mcp-lane semantics, `/__version` identity) belongs **upstream
in Greenlight**, never pasted into the wrapper.

---

## Remote BAM use cases (why we don't disable remote files)

Remote streaming is a core differentiator. Four real use cases, all HTTP(S) only
(`REMOTE_FILE_SCHEMES = ("http://","https://")`; "S3" = virtual-hosted HTTPS, not `s3://`):

| # | Use case | Why remote fetch is required |
|---|---|---|
| 1 | **Public datasets** (GIAB, 1000 Genomes) | Multi-GB/TB BAMs — pysam does **region-scoped random access** over HTTP range requests; fetch only the reads in the window |
| 2 | **Cloud-hosted user data** (S3/GCS/Azure HTTPS URLs) | User's aligned BAM in a bucket; same range-scoped streaming, no local copy of a 100 GB file |
| 3 | **Remote reference FASTA** (CRAM decode) | CRAM is reference-compressed; `list_contigs` deliberately suggests remote UCSC FASTA URLs |
| 4 | **`.bai`/`.crai` index download** | **Mandatory for 1–3.** Random access over HTTP needs the index; htslib fetches `<url>.bai`, `BAMIndexCache` caches it (24 h TTL) |

Disabling `BAMCP_ALLOW_REMOTE_FILES` kills use cases 1–3. **The fix is a curated allowlist**, not a kill switch.

---

## The plan

Sizes: **S** ≤ half-day · **M** ≈ 1–2 days · **L** ≈ 3–5 days. Repo tags: **[sub]** BAMCP
submodule · **[wrap]** `RTrentJones.dev` wrapper · **[up]** upstream Greenlight.

### Phase 0 — Stop the credibility bleed  · ~2–3 days
Cheapest wins: align words with reality, de-flake CI. Do first — changes how everything reads.

- [x] **(S) [sub]** Docs truth pass: fixed the broken README table (VCF prose had split it);
  scoped the MCP-Apps claim to *"first genomics MCP server"*; noted `uvx bamcp` needs the
  pending PyPI publish; bounded the 60 kb GIAB headline; concrete SECURITY.md reporting route.
- [x] **(M) [sub]** Test hermeticity: added a `network` marker (excluded by default via
  `addopts -m 'not network'`), marked the 2 live gnomAD tests, mocked `socket.getaddrinfo` in
  the 4 `example.com` validation tests, added `make test-network`. Suite is green offline
  (696 passed, 2 deselected); `eval-smoke` untouched.
- [x] **(S) [sub]** `uv.lock` regenerated; added a CI `lockfile` job (`uv lock --check` +
  frozen install + ship-gate from the pinned env). The 3.10–3.12 matrix stays as the compat check.
- **Exit:** ✅ PR suite fully hermetic + green offline; every doc claim defensible.

### Phase 1 — Make "production-ready" honest (the P0 security boundary)  · ~1–1.5 weeks

- [x] **(S) [wrap] STOPGAP:** set `BAMCP_ALLOWED_REMOTE_HOSTS` (+ `BAMCP_RATE_LIMIT_TRUSTED_PROXIES`)
  in `infra/bamcp.tf` to curated genomics hosts. **Done on wrapper branch `feat/bamcp-remote-allowlist`
  (not applied — needs sign-off on user cloud-bucket hosts).**
- [x] **(M) [sub]** **Reference-URL validation** — `validate_reference_path()` closes the one open
  SSRF path (allowlist-aware, shares the `validate_remote_url` core). Wired into all 7 handlers.
- [x] **(M–L) [sub]** **Unified data-source validation** — `validate_data_sources()` is the one
  choke point over BAM/CRAM/VCF/reference/index. Index + BAM/VCF were already validated; reference
  folded in. Option A (allowlist-first) chosen and documented; Option B (BAMCP-proxied fetch) deferred.
  Adversarial tests: rebinding, metadata endpoint, CRAM-reference, allowlist. *(Temp-index atomic
  write/sweep — small follow-up; see Deferred.)*
- [x] **(M) [sub] + [wrap]** **Auth: deleted the theater** — dynamic registration is a config flag
  (default OFF): `/register` withheld + `register_client` refuses. Scope required when auth is on.
  Verify invariants preserved (anon→401, `BAMCP_VERIFY_TOKEN`→authorized). Adversarial tests added.
  `verify/bamcp.config.ts` unchanged (contract shape preserved).
- [x] **(M) [sub]** Rate limiting + bounded state: trusted-proxy XFF, bounded-LRU limiter table;
  `region_cache` LRU-bounded; **registry fixed** (bounded LRU + identity guard + explicit removal —
  GC-finalizer can't work, `services.config` pins the config); distinct cache-miss sentinel so
  negative ClinVar/gnomAD results aren't re-fetched.
- [x] **(M) [sub]** Parsing isolation: all pysam parses confined to a dedicated bounded
  `ThreadPoolExecutor` (backpressure + no starvation of the loop). Killable worker-process
  cancellation noted as a deliberate follow-up.
- **Exit:** ✅ security claims are now *true*, encoded by adversarial tests. Commits: reference SSRF,
  auth service-token, rate-limiter, caches/registry, parse pool.

### Phase 2 — Typed contracts & coordinates (kills a bug class)  · ~1 week

- [ ] **(L) [sub]** Payload contracts: replace `dict[str, Any]` at the server↔viewer boundary
  with Pydantic / precise `TypedDict`; generate JSON Schema → TS types; **version the payload**;
  schema-compat test against the compiled viewer.
- [ ] **(M) [sub]** Viewer reference/context preservation (P1): persist full data-source context
  (incl. explicit `reference`) as typed viewer state; assert every pan/zoom/detail refetch
  carries it. Fixes silent CRAM/mismatch-evidence loss on navigation — falls out of the typed contract.
- [ ] **(M) [sub]** Coordinate type: a dedicated 0-based-half-open vs 1-based boundary API; fix
  `jump_to`, NCBI gene conversion, end-inclusivity, and stale planted positions in eval prompts.
- **Exit:** the backend/frontend contract is machine-checked, not hand-maintained.

### Phase 3 — Evaluation credibility (weakest scorecard line, 4.5/10)  · ~1.5 weeks

- [ ] **(M) [sub]** Grader redesign: score **tool-selection, argument-correctness, and
  answer-correctness separately**; deterministically extract/grade required numeric/factual
  claims instead of passing on tool presence. *(Lands in the nightly tier — never the ship-gate.)*
- [ ] **(M) [sub]** Real MCP path: finish `MCPStdioRouter`; run model evals against the actual
  FastMCP schemas + transport, not `InProcessRouter` with permissive schemas.
- [ ] **(M) [sub]** Rigor: tool/no-tool + viewer/no-viewer ablations; a small human-labeled
  benchmark; judge-swap + human-concordance (avoid same-family judge); publish per-case results
  + failure examples, not just aggregate pass rates.
- [ ] **(S) [sub] + [wrap]** Nightly parity: run the indel truth set + documented GIAB eval;
  wire both providers to match the `provider` input; emit a durable trend report.
- [ ] **(L) [sub]** *(optional, high-signal)* At least one real-read/broader GIAB benchmark;
  narrow scientific claims to match.
- **Exit:** evals measure task success, run through the real server, and are reproducible/published.

### Phase 4 — Frontend quality & module decomposition  · ~1 week

- [ ] **(M) [sub]** Frontend static/unit layer: `tsc --noEmit` (strict), ESLint, fast unit tests
  (coordinate conversion, data store, render math); a smaller browser smoke suite + focused
  visual-regression tests split from the 3,181-line monolith.
- [ ] **(M–L) [sub] + [wrap]** Split god modules: `tools.py` (~1,234) by use case; `parsers.py`
  (~994) into source adapters + pipeline stages; `mcp-app.ts` (~1,247) into
  navigation/transport/selection/detail; `renderer.ts` (~1,250) into independent tracks; E2E by
  behavior. **Mirror any tool rename/add/remove into `verify/bamcp.config.ts` `TOOLS`** (exactTools guard).
- **Exit:** frontend has backend-grade static discipline; no single file carries the system.

### Phase 5 — Packaging, release engineering, deployment docs  · ~3–4 days

- [ ] **(M) [sub]** Clean-room packaging: CI builds a wheel, installs it in a fresh env, verifies
  the **bundled viewer** loads; integrate the frontend build into Python packaging robustly.
  **Keep `/__version` reporting the built SHA** (ship-gate identity check).
- [ ] **(S) [sub]** Release story: tags + CHANGELOG across the 150 commits.
- [ ] **(S) [sub]** Transport docs: lead remote-deploy docs with **Streamable HTTP**
  (SSE deprecated); keep SSE only as compat. *(Runtime already streamable-http — docs-only.)*
- [ ] **(S) [sub]** Engineering retrospective — the highest-signal hiring artifact: what you
  designed, what agents produced, what verification caught, which trade-offs remain. Write last.
- **Exit:** `pip install` from a clean env yields a working, viewer-complete server; story matches code.

---

## Sequencing

```
Phase 0 (credibility)  ──►  Phase 1 (security P0s)  ──►  Phase 2 (contracts)
                                                              │
                                            ┌─────────────────┴───────────────┐
                                        Phase 3 (evals)                 Phase 4 (frontend)
                                            └─────────────────┬───────────────┘
                                                          Phase 5 (packaging + retro)
```
0→1→2 is a hard chain. 3 and 4 parallelize once 2 lands. 5 is last. Review's own split:
**Phases 0–1 = "before production-ready"**, **Phases 2–5 = "before flagship interview artifact."**

## Phase 2–5 status (what landed this branch)

- **Phase 2** — [x] 2b viewer reference/context preservation; [x] 2a payload `schema_version`
  stamp + server↔viewer compat test. *Deferred: full Pydantic/TypedDict + TS codegen (2a-full),
  pervasive coordinate newtype (2c).*
- **Phase 3** — [x] 3a grader scores tool-selection vs answer-correctness separately
  (`expected_claims`); [x] 3d nightly runs both truth sets + both providers. *Deferred: 3b real-MCP
  router, 3c human-labeled benchmark + judge-swap, 3e real-read GIAB.*
- **Phase 4** — [x] 4a `tsc --noEmit` + ESLint + vitest + CI job (caught + fixed dead code).
  *Deferred: 4b god-module splits.*
- **Phase 5** — [x] new env vars documented; [x] transport docs lead Streamable HTTP; [x] CHANGELOG;
  [x] engineering retrospective. *Deferred: 5a clean-room wheel+viewer CI.*

## Deferred follow-ups (scoped, not rushed)

Each is genuinely multi-day and higher-quality with dedicated focus; rushing them would reproduce
the "large volume of AI-assisted code" the review flagged. Listed so they're tracked, not dropped:

- **2a-full / 2c** — generated typed payload contract (Pydantic → JSON Schema → TS) and a coordinate
  newtype threaded through the codebase.
- **3b / 3c / 3e** — finish `MCPStdioRouter` and run evals through the real transport; a small
  human-labeled benchmark with judge-swap concordance; a broader real-read GIAB benchmark.
  Also **OpenAI multi-turn tool support**: the runner records follow-up turns in Anthropic
  message shape, so `OpenAIProvider` needs provider-specific translation before the OpenAI eval
  path is runnable (it is gated off in the nightly workflow until then).
- **4b** — split `tools.py` / `parsers.py` / `mcp-app.ts` / `renderer.ts` / the E2E monolith.
- **5a** — CI builds a wheel, installs it clean, and verifies the bundled viewer + `/__version`.
- **1c tail** — atomic temp-index writes; Option B (BAMCP-proxied fetch) for full redirect closure.
- **1e tail** — killable worker-process parsing for true cancellation.

## Immediate next actions

1. **[wrap] Set `BAMCP_ALLOWED_REMOTE_HOSTS`** — decide whether prod needs specific cloud-bucket
   hosts beyond the curated public list. Ship the Terraform stopgap.
2. **[sub] Phase 0** — doc truth-pass (incl. scoped MCP-Apps claim), test hermeticity split, `uv.lock`.
3. **[sub] Phase 1** — reference validator → `resource_policy()` unification → auth-theater
   removal (+ adversarial e2e) → rate-limit/cache bounds, each as its own branch through the loop.
