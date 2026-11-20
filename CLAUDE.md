# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

BAMCP — interactive BAM/CRAM variant visualization MCP server with an MCP Apps UI.

## Architecture

Entry point: `src/bamcp/__main__.py` — reads `BAMCPConfig.from_env()`, creates the FastMCP
server, adds HTTP middleware, runs with the selected transport.

```
src/bamcp/
  server.py        — FastMCP setup, tool/resource registration, auth wiring
  config.py        — BAMCPConfig dataclass with from_env()
  constants.py     — shared defaults, thresholds, bins
  tool_specs.py    — TOOL_SPECS: each tool's description + JSON schema (the one source of
                     descriptions; server.py registers from it, the eval router reads it)
  resources.py     — get_viewer_html(): serves the bundled static/dist/viewer.html
  demos.py         — demo datasets + `demo:<NAME>` file_path alias (resolve_demo_alias)
  glossary.py      — plain-language glossary (bamcp://glossary)
  prompts.py       — SERVER_INSTRUCTIONS, bamcp://demos|glossary|guide resources, prompts
  viewer_build.py  — viewer source-hash build stamp (run by make build-viewer)
  core/
    tools.py         — region/variant/contig/gene tool handlers, RegionRequest preamble, probe_site
    services.py      — per-server runtime: parse pool, services registry, region cache
    remote.py        — remote index/reference-index download, cold-start warming
    lookups.py       — lookup_clinvar / lookup_gnomad handlers + status vocabulary
    scan.py          — scan_variants: time-bounded chunked SNV scan
    filters.py       — failed/passes_candidate_filters: the one min_vaf/min_depth/min_alt_count predicate
    parsers.py       — fetch_region(): pysam parsing, reads, coverage, variant detection
    validation.py    — region/path/URL validation, SSRF prevention
    serialization.py — serialize_region_data(): RegionData → JSON + evidence
    cache.py         — BAMIndexCache: file cache for remote BAM index files (TTL)
    reference.py     — genome build registry, detect_genome_build()
    errors.py        — one error envelope: ToolInputError, sanitize_tool_exception
    prewarm.py       — HTTP-only startup prewarm of demo datasets (BAMCP_PREWARM)
  analysis/
    evidence.py      — variant evidence, artifact risk scoring, confidence
    curation.py      — curation summaries + recommendations
    acmg.py          — ACMG criteria scoring
    fusion.py        — multi-source evidence fusion
  clients/
    clinvar.py, gnomad.py, genes.py — async NCBI / gnomAD GraphQL / NCBI gene clients
    ncbi.py          — shared NCBI E-utilities rate limiting + error detection, redact_error
    ttl_cache.py     — BoundedTTLCache (shared LRU + TTL)
  middleware/
    auth.py          — BAMCPAuthProvider (in-memory OAuth 2.0)
    bodylimit.py     — BodySizeLimitMiddleware: 413 on oversized request bodies
    ratelimit.py     — sliding-window IP rate limiting (IPv6 keyed by /64)
    security.py      — security response headers
    telemetry.py     — configure_telemetry() (optional OTel)
  eval/              — evaluation harness: truthset (deterministic gate), runner, router,
                       grader, metrics, report, compare, providers, renderer, schema, cli,
                       acmg, models (model-ID registry), tracer (Tracer publisher),
                       render_cli, visual/ (pre-registered visual-vs-text eval)
  static/            — viewer TS/Vite app (see "Viewer UI")
```

## Commands

```bash
make verify           # every ci.yml gate in one command (VERIFY_E2E=1 adds playwright) — run before pushing
make install          # build-viewer + pip install -e ".[dev]"
make build-viewer     # npm install && npm run build in static/, then python -m bamcp.viewer_build (source-hash stamp)
make test             # pytest (unit + integration, ignores e2e)
make test-network     # live-network tier only (excluded from the default suite)
make test-e2e         # playwright install chromium && pytest tests/e2e/
make lint             # ruff format --check AND ruff check
make format           # ruff format && ruff check --fix
make typecheck        # mypy src
make coverage         # HTML coverage report
make coverage-strict  # coverage with --cov-fail-under=85 (the CI bar)
make eval-smoke       # deterministic ground-truth gate (required ship-gate)
make eval-visual-dry  # visual-vs-text eval hermetic dry run (fixtures, stimuli, leak-lint, grader self-tests) — in ci.yml eval-smoke job
make giab-benchmark   # real GIAB benchmark — network-gated, not in CI
make render-viewer    # render the viewer to a PNG (.render-dev/last.png); finds an installed Chromium if Playwright's default is missing (override: BAMCP_CHROMIUM_EXECUTABLE)
make docker-build / docker-test / clean
```

`./dev_server.py` launches the server under the MCP Inspector (http://localhost:6274).

`make verify` (`scripts/verify-local.sh`) is the local regression loop: it runs against *its own
checkout's* `src/` (safe from a git worktree, unlike a bare `python -m bamcp...`, which imports
whichever editable install is on the path), restores the tracked fixtures that
`create_fixtures.py` rewrites, and ends with `OK — all local gates green`.

## CI gates (get these wrong and the push fails)

`.github/workflows/ci.yml` runs `lockfile`, `lint`, `test` (3.10/3.11/3.12), `eval-smoke`,
`frontend`, `e2e`, `docker`.

- **`uv.lock` must stay in sync.** Any `pyproject.toml` edit requires re-running `uv lock`
  — the `lockfile` job runs `uv lock --check` and fails otherwise.
- **`make lint` is two commands.** `ruff format --check` runs *before* `ruff check`; running
  only `ruff check` passes locally and then fails CI on formatting.
- **Coverage has two bars.** `fail_under = 80` in pyproject, but CI runs `make coverage-strict`
  (**85%**) on 3.12. Target 85.
- **`make eval-smoke` is a required ship-gate**, in both `ci.yml` and `greenlight-build.yml`.
  It needs `python tests/create_fixtures.py` to have run and the `[eval]` extra installed.
- **Eval results publish to Tracer** (`bamcp.eval.tracer`, see EVALS.md "Trending in Tracer") from
  `ci.yml` on `main` and from `eval-nightly.yml`. It is best-effort and never fails a job; it needs the
  `TRACER_INGEST_TOKEN` repo secret. If you change `TruthsetReport.meets_floors`, keep
  `tracer.DEFAULT_FLOORS` and `_floor_cases` in step (`test_tracer.py` pins both).
- **The frontend has its own gate.** In `src/bamcp/static/`: `npm run typecheck`, `npm run lint`
  (eslint), `npm test` (vitest), `npm run build` — all four run in CI.
- **mypy targets Python 3.12** (not the 3.10 floor). New imports from the `[eval]` or
  `[telemetry]` optional extras need an entry in `[[tool.mypy.overrides]]` or `make typecheck`
  fails on the CI lint job, which installs `[dev]` only.
- **The viewer bundle is part of the ship-gate.** `test_resources` asserts the built
  `static/dist/viewer.html` exists, so `make build-viewer` must have run.

`pre-commit` is configured (ruff-format + ruff `--fix`): `pip install pre-commit && pre-commit install`.

## Registered tools

Thirteen tools, all registered in `server.py`:

| Tool | Notes |
|------|-------|
| `visualize_region` | Reads + UI. Auto-detects compact mode for large regions. |
| `jump_to` | Jump to a position with a configurable window (UI + data). `contig` is required — never default it. |
| `get_variants` | `variant_source` (`auto`/`vcf`/`bamcp`) selects the source: `vcf` makes a caller's VCF authoritative and attaches BAMCP read-level evidence per site; `auto` overlays a VCF on local candidates; `bamcp` is local-only. |
| `scan_variants` | Time-bounded SNV scan of a contig or `region`; a truncated scan returns `next_region` (pass it back as `region` to resume). Output sorted by position. An unknown contig names the file's own spelling (`parsers.contig_not_found_message`). |
| `get_coverage` | Depth-of-coverage statistics. |
| `list_contigs` | Contigs + genome build detection + suggested public reference URL. |
| `get_region_summary` | Text-only summary for LLM reasoning (no UI). |
| `get_variant_curation_summary` | Curation summary with artifact risk assessment. The *requested* site is probed (`core/tools.py:probe_site`: min_baseq=0, no VAF/alt floor) when it fails the calling filters, and `filter_check` reports `passes_filters` + `failed_filters`; `classify_variant`'s observation does the same (`status: below_filters`). |
| `classify_variant` | ACMG classification via evidence fusion (`analysis/acmg.py`, `fusion.py`). |
| `lookup_clinvar` / `lookup_gnomad` | Clinical significance / population allele frequency. |
| `search_gene` | Gene symbol → coordinates. Build-aware: auto-detects from the BAM. |
| `cleanup_cache` | Clear downloaded index files + region cache, and forget which remote inputs are warm. stdio only (no-op on shared HTTP servers); never closes the shared HTTP clients. |

Also registered: the `ui://bamcp/viewer` resource; `bamcp://demos`, `bamcp://glossary`,
`bamcp://guide` resources; prompts `explore_demo`, `check_variant`, `explain_terms`; and
`GET /__version` (+ `/mcp/__version`), which serves the image's `GREENLIGHT_SHA` for Greenlight's
verify gate. Every `file_path` accepts a `demo:<NAME>` alias (`demo:HG002`, `demo:HG00096`),
resolved in `_run`.

**Every variant response echoes its `applied_filters`** — absence of a variant must stay
distinguishable from "filtered out". Preserve that when touching variant output. That includes
what *couldn't* be detected (`snv_detection` without a reference, `indel_detection` above the
size cap, `reads_truncated`, scan truncation), not just the thresholds.

**Coordinates.** Region strings at the tool boundary are **1-based inclusive, samtools-style**
(`chr1:1051-1051`, or just `chr1:1051`, is one base); `parse_region` converts once to 0-based half-open internally and
`format_region` is its inverse. Output positions are 1-based. The viewer's `ui/init` payload is
the one exception (`coordinate_system: "0-based-half-open"`). External lookups (ClinVar, gnomAD,
gene) must be queried in the **BAM's** assembly — pass `build`, never assume the server default.

## Configuration

All settings come from `BAMCP_*` env vars. The authoritative list is the dataclass itself:
`@src/bamcp/config.py` — read it rather than relying on a copy here.

Transport is `BAMCP_TRANSPORT`: `stdio` (default), `sse`, or `streamable-http`. HTTP transports
use `BAMCP_HOST` (default `127.0.0.1`; the Docker image sets `0.0.0.0`) / `BAMCP_PORT` (default `8000`).

Auth is opt-in via `BAMCP_AUTH_ENABLED=true` — an in-memory OAuth 2.0 authorization server
implementing `OAuthAuthorizationServerProvider` (auth code flow, refresh with rotation,
revocation). It also accepts a stateless service token (`BAMCP_VERIFY_TOKEN`). Dynamic client
registration is a separate opt-in (`BAMCP_ALLOW_DYNAMIC_REGISTRATION=true`, default off); auth
with neither a token nor registration fails at startup.

Security-relevant settings: `BAMCP_ALLOW_REMOTE_FILES` (default off), `BAMCP_ALLOWED_REMOTE_HOSTS`,
`BAMCP_ALLOWED_DIRECTORIES`, `BAMCP_ALLOW_LOCAL_FILES` (default on; off rejects every
caller-supplied local path, the configured `BAMCP_REFERENCE` excepted), `BAMCP_MAX_REQUEST_BYTES`
(default 1 MiB), `BAMCP_TRUSTED_HOSTS`, `BAMCP_RATE_LIMIT*`.

## Genome build

`list_contigs` auto-detects GRCh37 vs GRCh38 from chr1 length and suggests a public FASTA URL
when no local reference is configured (see `core/reference.py`).

**Any reference URL must be openable by pysam for random access**, which means uncompressed or
**bgzip**-compressed (never plain gzip) *and* carrying a faidx index (`.fai`, plus `.gzi` for
bgzip) at the same URL. The UCSC `bigZips` files (`hg19.fa.gz`, `hg38.fa.gz`) satisfy neither and
fail with a bare `error when opening file` — `GENOME_BUILDS` used to advertise them, which sent
callers into an unfixable error. Re-verify (HTTP 206 on the FASTA, `.fai`, and `.gzi`) before
changing any URL there. The suggestion is also **contig-naming-aware**: a `chr1`-style FASTA
cannot serve a BAM whose contigs are `1`, so the registry keeps one URL per naming style and
returns `None` (plus a `reference_note`) rather than suggesting a mismatch.

**Workflow rule: call `list_contigs` first on a new BAM to pin the build.** A build mismatch is
silent and shifts coordinates by ~0.5–2 Mb, which reads as false "benign" evidence downstream.
`search_gene` auto-detects the build from the BAM and warns loudly on mismatch — keep it that way.

## Security

HTTP middleware (SSE / streamable-http only), added in `__main__.py`:
`TrustedHostMiddleware` (DNS rebinding, `BAMCP_TRUSTED_HOSTS`) → `SecurityHeadersMiddleware` →
`BodySizeLimitMiddleware` (413; `BAMCP_MAX_REQUEST_BYTES`, default 1 MiB, 64 KiB on
`/register`, `/token`, `/revoke`, `/authorize`; checks `Content-Length` and counts streamed bytes) →
`RateLimitMiddleware` (default 120 req/min/IP — every MCP JSON-RPC message is a request; an IPv6
client is keyed by its /64). Behind
the Cloudflare tunnel the client IP comes from `CF-Connecting-IP`, then right-to-left
`X-Forwarded-For`, trusted only from `BAMCP_RATE_LIMIT_TRUSTED_PROXIES` (loopback).

The prod server is an **open public demo**: dynamic client registration with auto-approve is
intentional. It is off by default in code; prod turns it on in the wrapper's infra
(`RTrentJones.dev/infra/bamcp.tf`, `BAMCP_ALLOW_DYNAMIC_REGISTRATION = "true"`). Every in-memory store (OAuth clients/codes/tokens, index cache by bytes, region cache
by bytes) is therefore bounded — keep any new store bounded too. Registration also caps each
client record (`middleware/auth.py:check_client_metadata_size`: 10 `redirect_uris`, 2048 chars
each, 256-char `client_name`, 4 KiB of metadata).

**One scheme classifier.** Every caller-supplied data source (BAM/CRAM, VCF, reference) goes
through `core/validation.py:_is_remote_source`. Only an exact `http://` / `https://` prefix is
remote; any other value that starts like a URI scheme (`http:/` with one slash, `data:`, `ftp:`,
`file:`, `s3:`, …) is rejected, as is a local path that embeds `://` or `##idx##` — htslib
opens those itself, outside the SSRF checks. With `BAMCP_ALLOW_LOCAL_FILES=false` every local
path is rejected before any filesystem access (no existence oracle); the operator's
`BAMCP_REFERENCE` is exempt. `validate_variant_file_path` requires the VCF/BCF to exist.

`validate_remote_url()` in `core/validation.py` (http/https only, ports 80/443 only) rejects any
`#` (htslib reads `x.bam##idx##<url>` as a second, unvalidated URL), userinfo, and backslashes;
then checks `BAMCP_ALLOWED_REMOTE_HOSTS`; then resolves the host and blocks every non-global
address (`not ip.is_global`: loopback, RFC1918, link-local/metadata, CGNAT, multicast, …).
Index downloads follow redirects manually, re-validating each hop. Redirects inside
htslib/libcurl and DNS rebinding can't be checked per-URL, so **the exact-host allow-list is the
real SSRF control** — prod sets one, and startup warns without it. Validation resolves DNS, so
call it off the event loop (`await asyncio.to_thread(validate_path, ...)`).

At startup `__main__.py` sets `REF_PATH` / `REF_CACHE` (when the operator has not) to a local
directory under the cache dir, so htslib's CRAM reference-by-M5 lookup never goes to the network;
an operator value naming a URL is kept and logged with a warning. A CRAM header's `@SQ UR:` tag
is the other route: htslib decodes a contig against the file it names when no reference is given,
or when the given reference lacks that contig. With `BAMCP_ALLOW_LOCAL_FILES=false`, every
alignment open (`core/parsers.py:_open_alignment`, under `guard_cram_decoding`, which `server._run`
sets per call) refuses a CRAM with no resolved reference, or one whose `UR:`-tagged contigs the
reference does not cover. **It decides on the opened file's content** (`AlignmentFile.is_cram`),
never the name — htslib decodes a CRAM named `.bam` as CRAM; `is_cram(path)` is only for the
early `_require_cram_reference` refusal and for choosing an index suffix. The guard is a
ContextVar; `_run_in_parse_pool` copies the context into the parse thread. htslib itself refuses
a remote `UR:` (pinned by `test_cram_reference_gate.py`).

Input validation: regions against `REGION_PATTERN`, paths ≤2048 chars, regions ≤100 chars,
only `.bam`/`.cram` for local files, gene symbols against `GENE_SYMBOL_PATTERN` (1–40 of letters,
digits, `.`, `_`, `-`; they go into an Entrez query), error messages sanitized (never leak config).
Unavailable lookups carry a fixed `reason_code` (`timeout`, `network_error`, `upstream_error`,
`unsupported_build`) — upstream error text is logged, never returned, and logged only through
`clients/ncbi.py:redact_error` (an httpx error's text embeds the request URL, NCBI `api_key` and
all); a filter on the `httpx` logger scrubs the same parameters from its per-request INFO line.

**One error envelope.** Every tool runs through `_run()` in `server.py`, which maps any handler
exception via `core/errors.py:sanitize_tool_exception` to an `isError=true` result with a
client-safe message. Handlers raise `ToolInputError` for caller-facing failures — don't return
`{"error": ...}` as a success. Lookup `not_found` is an answer, not an error; invalid input (bad allele, position, build) is an
error. Argument-schema failures (wrong/missing parameter names) never reach `_run` — FastMCP
validates first — so `_install_argument_error_mapping` maps them to one short line listing the
tool's expected parameters instead of a raw pydantic dump.
Tool descriptions live only in `tool_specs.TOOL_SPECS` (`server.py` registers
`description=TOOL_SPECS[name]["description"]`); its schemas must match the server's parameter
sets (parity test in `test_tool_specs.py` checks parameters and descriptions).

`mcp>=1.23.0` is a hard floor (CVE-2025-66416 DNS rebinding fix).

Prod compose hardening: `cap_drop: [ALL]`, `no-new-privileges`, `read_only`, `memory: 2g`, non-root,
port bound to `127.0.0.1`, `BAMCP_ALLOWED_REMOTE_HOSTS` set, `BAMCP_ALLOW_LOCAL_FILES=false`.

**Resource bounds per call** (`constants.py`): `VCF_MAX_RECORDS` (5,000; `vcf_truncated` +
`vcf_max_records` in `applied_filters`), `VCF_MAX_SAMPLES_DETAIL` (16; GT only for the first 16,
`samples_truncated` + `sample_count` on the variant; `get_region_summary` lines end
`… +N more samples`), `READ_MAX_MISMATCHES` (1,000 per read;
`mismatches_truncated`), `READ_PAYLOAD_MAX_BYTES` (64 MiB; sets `reads_truncated`),
`GET_VARIANTS_MAX_RESULTS` (1,000; `truncated` / `total_candidates` / `max_results` in
`applied_filters`, VCF records kept first, then most alt reads; the viewer tools' payload uses
the same cap and reports `variants_truncated` / `total_candidates` / `max_results`, which the
viewer shows in its status banner), `SCAN_VARIANTS_MAX_RESULTS`
(500; `scan.py:_TopCandidates` keeps only the highest-VAF ones *while* scanning — the same picks,
ties included, as a sort at the end — and counts the rest for `total_candidates`).

## Viewer payload

The viewer tools (`visualize_region`, `jump_to`) return two things: the **text content** is the
compact LLM-facing summary (headline, top variant interpretations, applied filters — bounded by
`VIEWER_SUMMARY_MAX_CHARS`), while **`structuredContent`** (built from the handler's
`_meta["ui/init"]`) carries the full UI payload — reads, coverage, per-variant evidence — that
the MCP Apps viewer reads. Don't move the payload out of `structuredContent` (the viewer depends
on it) and don't grow the text back into a read dump.

## Viewer UI

TypeScript/Vite app in `src/bamcp/static/`. Sources: `viewer.html`, `client.ts` (MCP Apps SDK
wrapper), `mcp-app.ts` (app logic), `renderer.ts` (canvas), `layout.ts`, `palette.ts`,
`ui-sync.ts`, `state.ts`, `data-store.ts`, `types.ts`, `constants.ts` (vitest: `*.test.ts` for
constants/layout/palette/state/ui-sync).

```bash
cd src/bamcp/static
npm install && npm run build   # → dist/viewer.html
npm run watch                  # vite build --watch (rebuilds dist/ on change)
npm run typecheck && npm run lint && npm test   # what CI runs
```

**Never hand-edit anything in `src/bamcp/static/dist/`** — no `cat`, no `echo`, no manual copies.
Edit the source and run `npm run build`; Vite inlines the SDK and bundles everything.

After any `static/` edit, close the visual loop: run the **`render-viewer` skill** (or
`make render-viewer`) and `Read` the resulting PNG to confirm the change did what you intended.

### MCP Apps SDK

`callServerTool()` returns its result via the **Promise**, not a callback:

```typescript
// WRONG — ontoolresult will NOT fire for callServerTool
await app.callServerTool({ name: 'tool', arguments: {} });

// CORRECT
const result = await app.callServerTool({ name: 'tool', arguments: {} });
if (result.structuredContent) { handleData(result.structuredContent); }
```

`ontoolresult` / `ontoolinput` fire only for tool calls initiated by the **host** (the LLM),
never for calls the app makes itself.

### Non-obvious viewer behavior

- Regions >`COMPACT_REGION_THRESHOLD_BP` (500 bp) auto-omit sequences (compact mode) to shrink the
  payload. The server sends the threshold as `compact_threshold_bp` in `ui/init`; the viewer uses
  it and falls back to `DEFAULT_COMPACT_THRESHOLD_BP` (500) for older payloads.
- At zoom scale ≥10, `checkAndRequestSequences()` auto-fetches missing sequences via
  `callServerTool('visualize_region')`, with a "Load Detail" button as a 3s fallback.

## Tests

- Markers: `unit`, `integration`, `network`, `e2e`. Async tests use `pytest-asyncio`
  (`asyncio_mode = "auto"`).
- **The default suite is hermetic.** `addopts` carries `-m 'not network'`, so live-external-call
  tests are excluded by default and never flake offline. Run them explicitly with
  `make test-network` (a command-line `-m` overrides the default).
- E2E uses Playwright's **sync** API — it conflicts with pytest-asyncio's loop, so it runs
  separately (`tests/e2e/`, excluded from `addopts`).
- Fixtures in `tests/fixtures/` (`small.bam`, `ref.fa`, `empty.bam`, `comprehensive.bam` +
  `comprehensive_ref.fa`), generated by `tests/create_fixtures.py`.
- HTTP mocking via `pytest-httpx` for the ClinVar/gnomAD clients.
- Eval truth sets live in `tests/eval/datasets/` (`synthetic_v1`, `indel_v1`, `acmg_v1`, `giab`),
  each with a `manifest.yaml`. Floors are enforced in code (`TruthsetReport.meets_floors`) *and*
  in the CI step, so changing them is a reviewed change. Methodology: `@EVALS.md`.

## Repo etiquette

- **`main`-only — there is no `develop` branch here.** Work on a `feat/*`/`fix/*` branch, open a
  PR against `main`.
- Commits follow conventional style with the PR number: `fix(viewer): distinct ruler tick
  labels at zoom (#53)`.
- **Merging to `main` deploys straight to prod** (see below). There is no beta gate in this repo;
  the `beta`/`preview` docker-compose profiles are local-only.

## Deployment

Deploy is owned by Greenlight, not this repo — there is **no `deploy.yml` here**. On every push
to `main`:

1. **`.github/workflows/greenlight-build.yml`** — ship-gate (`make test` + `make eval-smoke`),
   then builds the arm64 image natively and pushes to GHCR
   (`ghcr.io/rtrentjones/bamcp:prod` moving tag + `:<sha>` immutable), then fires
   `repository_dispatch(deploy-bamcp)` at the wrapper repo (`RTrentJones/RTrentJones.dev`)
   via `GREENLIGHT_DISPATCH_TOKEN`.
2. The wrapper's **`greenlight-deploy-bamcp.yml`** runs the `oci-deploy-verify` composite:
   restart the OCI Container Instance to re-pull `:prod`, verify prod is serving the new image,
   then post a `greenlight/deploy-bamcp` commit status back here. OCI creds + the verify token
   live on the wrapper; this repo only needs `GREENLIGHT_DISPATCH_TOKEN`.

The container runs on a private OCI subnet with no public IP. External access is via a
**Cloudflare Tunnel** — a `cloudflared` sidecar in the same Container Instance, outbound-only
(no inbound rules, no public IP, no load balancer). Setup: `./scripts/setup-cloudflared.sh`
(recreates the instance with the sidecar); needs `CLOUDFLARE_TUNNEL_TOKEN`.

`docker/entrypoint.sh` respects a cloud-provider `PORT` env var, falling back to `BAMCP_PORT`.

## Docker

`Dockerfile` (prod, multi-stage slim, multi-arch; CI builds arm64) · `Dockerfile.dev` (adds pytest/playwright/ruff/mypy) ·
`docker-compose.yml` profiles: `dev` (test + e2e + lint services), `beta`, `preview`, `prod` ·
`docker/healthcheck.py` — the image `HEALTHCHECK` GETs `/__version` on HTTP transports (import
check for stdio); `docker/entrypoint.sh` runs it with `--preflight` (import/config/server-creation
check only, the server isn't listening yet). The viewer build's Node version is pinned once, in
`src/bamcp/static/.nvmrc`; every workflow reads it.

## Further reading

- `@docs/ROADMAP.md` — living strengths/weaknesses + current priorities. Update it as work lands.
- `@EVALS.md` — evaluation methodology · `@SAFETY.md` — the overconfidence guard.
- `@docs/engineering-review-response.md`, `@docs/ENGINEERING_RETROSPECTIVE.md` — why things are
  the way they are.
- `archived/` — superseded planning docs. Historical only; don't treat as current design.

## Greenlight loop

This repo is a Greenlight tool consumed as a submodule by `RTrentJones.dev`. The generic
deploy→verify→promote model is in `.claude/skills/deploy-verify-promote/SKILL.md`, and
per-provider detail in `.claude/skills/provider-*`. BAMCP's cell of that matrix is
**oci, direct-to-prod on `main`** — the ship-gate in `greenlight-build.yml` is the safety net
that a promote step would otherwise provide.
