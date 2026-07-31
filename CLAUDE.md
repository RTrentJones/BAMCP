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
  tool_specs.py    — JSON-schema specs per tool (get_tool_spec)
  resources.py     — get_viewer_html(): serves the bundled static/dist/viewer.html
  core/
    tools.py         — tool handlers + cache helpers
    parsers.py       — fetch_region(): pysam parsing, reads, coverage, variant detection
    validation.py    — region/path/URL validation, SSRF prevention
    serialization.py — serialize_region_data(): RegionData → JSON + evidence
    cache.py         — BAMIndexCache: file cache for remote BAM index files (TTL)
    reference.py     — genome build registry, detect_genome_build()
  analysis/
    evidence.py      — variant evidence, artifact risk scoring, confidence
    curation.py      — curation summaries + recommendations
    acmg.py          — ACMG criteria scoring
    fusion.py        — multi-source evidence fusion
  clients/
    clinvar.py, gnomad.py, genes.py — async NCBI / gnomAD GraphQL / NCBI gene clients
    ttl_cache.py     — BoundedTTLCache (shared LRU + TTL)
  middleware/
    auth.py          — BAMCPAuthProvider (in-memory OAuth 2.0)
    ratelimit.py     — sliding-window IP rate limiting
    security.py      — security response headers
    telemetry.py     — configure_telemetry() (optional OTel)
  eval/              — evaluation harness: truthset (deterministic gate), runner, router,
                       grader, metrics, report, compare, providers, renderer, schema, cli
  static/            — viewer TS/Vite app (see "Viewer UI")
```

## Commands

```bash
make install          # build-viewer + pip install -e ".[dev]"
make build-viewer     # cd src/bamcp/static && npm install && npm run build
make test             # pytest (unit + integration, ignores e2e)
make test-network     # live-network tier only (excluded from the default suite)
make test-e2e         # playwright install chromium && pytest tests/e2e/
make lint             # ruff format --check AND ruff check
make format           # ruff format && ruff check --fix
make typecheck        # mypy src
make coverage         # HTML coverage report
make coverage-strict  # coverage with --cov-fail-under=85 (the CI bar)
make eval-smoke       # deterministic ground-truth gate (required ship-gate)
make giab-benchmark   # real GIAB benchmark — network-gated, not in CI
make render-viewer    # render the viewer to a PNG (.render-dev/last.png)
make docker-build / docker-test / clean
```

`./dev_server.py` launches the server under the MCP Inspector (http://localhost:6274).

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
| `jump_to` | Jump to a position with a configurable window (UI + data). |
| `get_variants` | `variant_source` (`auto`/`vcf`/`bamcp`) selects the source: `vcf` makes a caller's VCF authoritative and attaches BAMCP read-level evidence per site; `auto` overlays a VCF on local candidates; `bamcp` is local-only. |
| `scan_variants` | Scan a broader region for candidate variants. |
| `get_coverage` | Depth-of-coverage statistics. |
| `list_contigs` | Contigs + genome build detection + suggested public reference URL. |
| `get_region_summary` | Text-only summary for LLM reasoning (no UI). |
| `get_variant_curation_summary` | Curation summary with artifact risk assessment. |
| `classify_variant` | ACMG classification via evidence fusion (`analysis/acmg.py`, `fusion.py`). |
| `lookup_clinvar` / `lookup_gnomad` | Clinical significance / population allele frequency. |
| `search_gene` | Gene symbol → coordinates. Build-aware: auto-detects from the BAM. |
| `cleanup_cache` | Clean the session's BAM index cache files. |

**Every variant response echoes its `applied_filters`** — absence of a variant must stay
distinguishable from "filtered out". Preserve that when touching variant output.

## Configuration

All settings come from `BAMCP_*` env vars. The authoritative list is the dataclass itself:
`@src/bamcp/config.py` — read it rather than relying on a copy here.

Transport is `BAMCP_TRANSPORT`: `stdio` (default), `sse`, or `streamable-http`. HTTP transports
use `BAMCP_HOST` (default `0.0.0.0`) / `BAMCP_PORT` (default `8000`).

Auth is opt-in via `BAMCP_AUTH_ENABLED=true` — an in-memory OAuth 2.0 authorization server
implementing `OAuthAuthorizationServerProvider` (dynamic client registration, auth code flow,
refresh, revocation).

## Genome build

`list_contigs` auto-detects GRCh37 vs GRCh38 from chr1 length and suggests a public UCSC FASTA
URL when no local reference is configured (see `core/reference.py`).

**Workflow rule: call `list_contigs` first on a new BAM to pin the build.** A build mismatch is
silent and shifts coordinates by ~0.5–2 Mb, which reads as false "benign" evidence downstream.
`search_gene` auto-detects the build from the BAM and warns loudly on mismatch — keep it that way.

## Security

HTTP middleware (SSE / streamable-http only), added in `__main__.py`:
`TrustedHostMiddleware` (DNS rebinding, `BAMCP_TRUSTED_HOSTS`) → `SecurityHeadersMiddleware` →
`RateLimitMiddleware` (default 60 req/min/IP).

`validate_remote_url()` in `core/validation.py` resolves hostnames and blocks private/internal IPs:
`127.0.0.0/8`, `10.0.0.0/8`, `172.16.0.0/12`, `192.168.0.0/16`, `169.254.0.0/16` (cloud metadata),
and IPv6 private/loopback/link-local. Optional allowlist via `BAMCP_ALLOWED_REMOTE_HOSTS`.

Input validation: regions against `REGION_PATTERN`, paths ≤2048 chars, regions ≤100 chars,
only `.bam`/`.cram` for local files, error messages sanitized (never leak config).

`mcp>=1.23.0` is a hard floor (CVE-2025-66416 DNS rebinding fix).

Prod compose hardening: `cap_drop: [ALL]`, `no-new-privileges`, `read_only`, `memory: 2g`, non-root.

## Viewer UI

TypeScript/Vite app in `src/bamcp/static/`. Sources: `viewer.html`, `client.ts` (MCP Apps SDK
wrapper), `mcp-app.ts` (app logic), `renderer.ts` (canvas), `state.ts`, `data-store.ts`,
`types.ts`, `constants.ts` (+ `constants.test.ts` for vitest).

```bash
cd src/bamcp/static
npm install && npm run build   # → dist/viewer.html
npm run dev                    # Vite dev server
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

- Regions >500 bp auto-omit sequences (compact mode) to shrink the payload.
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

`Dockerfile` (prod, slim, ARM64) · `Dockerfile.dev` (adds pytest/playwright/ruff/mypy) ·
`docker-compose.yml` profiles: `dev` (test + lint services), `beta`, `preview`, `prod` ·
`docker/healthcheck.py` (imports bamcp, verifies server creation).

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
