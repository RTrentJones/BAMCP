# BAMCP

Interactive BAM/CRAM variant visualization MCP server with MCP Apps UI.

## Architecture

```
FastMCP server (server.py)
  ├── Tool handlers (tools.py) → pysam parsers (parsers.py)
  ├── ClinVar client (clinvar.py) → NCBI E-utilities API
  ├── gnomAD client (gnomad.py) → gnomAD GraphQL API
  ├── UI resource (resources.py → static/viewer.html)
  ├── Auth provider (auth.py) — opt-in OAuth 2.0
  ├── Cache manager (cache.py) — remote BAM index caching
  └── Config (config.py) — all settings from env vars
```

**Entry point**: `src/bamcp/__main__.py` — reads `BAMCPConfig.from_env()`, creates FastMCP server, runs with selected transport.

## Key Modules

| Module | Role |
|--------|------|
| `server.py` | FastMCP setup, tool/resource registration, auth wiring |
| `tools.py` | All tool handlers (visualize, variants, coverage, contigs, jump, summary, ClinVar, gnomAD, curation) |
| `parsers.py` | `fetch_region()` — pysam BAM/CRAM parsing, read extraction, coverage, variant detection |
| `clinvar.py` | `ClinVarClient` — async NCBI E-utilities client for variant clinical significance |
| `gnomad.py` | `GnomadClient` — async gnomAD GraphQL client for population allele frequencies |
| `resources.py` | `get_viewer_html()` — serves bundled HTML from `static/viewer.html` |
| `config.py` | `BAMCPConfig` dataclass with `from_env()` classmethod |
| `auth.py` | `BAMCPAuthProvider` (in-memory OAuth 2.0), `build_auth_settings()` |
| `cache.py` | `BAMIndexCache` — file-based cache for remote BAM index files with TTL |
| `reference.py` | Genome build registry (`GENOME_BUILDS`), `detect_genome_build()`, `get_public_reference_url()` |

## Registered Tools

| Tool | Description |
|------|-------------|
| `visualize_region` | View aligned reads with interactive MCP Apps visualization (returns UI + data). Auto-detects compact mode for large regions. |
| `get_variants` | Detect and return variants in a region |
| `get_coverage` | Calculate depth of coverage statistics |
| `list_contigs` | List chromosomes/contigs and detect genome build (GRCh37/GRCh38) with suggested public reference URL |
| `jump_to` | Jump to a specific genomic position with configurable window (returns UI + data) |
| `get_region_summary` | Text-only region summary for LLM reasoning (no UI) |
| `lookup_clinvar` | Look up variant in ClinVar for clinical significance and conditions |
| `lookup_gnomad` | Look up variant in gnomAD for population allele frequency data |
| `get_variant_curation_summary` | Detailed curation summary for a variant with artifact risk assessment |
| `search_gene` | Search for a gene by symbol and return genomic coordinates (uses NCBI) |
| `cleanup_cache` | Clean up session's BAM index cache files |

## Transport Modes

Controlled by `BAMCP_TRANSPORT` env var:
- `stdio` (default) — standard MCP stdio transport
- `sse` — Server-Sent Events over HTTP
- `streamable-http` — Streamable HTTP transport

HTTP transports use `BAMCP_HOST` (default `0.0.0.0`) and `BAMCP_PORT` (default `8000`).

## Authentication

Opt-in via `BAMCP_AUTH_ENABLED=true`. In-memory OAuth 2.0 authorization server implementing `OAuthAuthorizationServerProvider`. Supports dynamic client registration, authorization code flow, token refresh, and revocation.

Related env vars: `BAMCP_ISSUER_URL`, `BAMCP_RESOURCE_SERVER_URL`, `BAMCP_REQUIRED_SCOPES`, `BAMCP_TOKEN_EXPIRY`.

## Configuration

All settings via environment variables. See `BAMCPConfig` in `src/bamcp/config.py`.

Core: `BAMCP_REFERENCE`, `BAMCP_MAX_READS`, `BAMCP_DEFAULT_WINDOW`, `BAMCP_MIN_VAF`, `BAMCP_MIN_DEPTH`, `BAMCP_MIN_MAPQ`

Transport: `BAMCP_TRANSPORT`, `BAMCP_HOST`, `BAMCP_PORT`

Auth: `BAMCP_AUTH_ENABLED`, `BAMCP_ISSUER_URL`, `BAMCP_RESOURCE_SERVER_URL`, `BAMCP_REQUIRED_SCOPES`, `BAMCP_TOKEN_EXPIRY`

External databases: `BAMCP_NCBI_API_KEY`, `BAMCP_CLINVAR_ENABLED`, `BAMCP_GNOMAD_ENABLED`, `BAMCP_GNOMAD_DATASET`, `BAMCP_GENOME_BUILD`

Cache: `BAMCP_CACHE_DIR` (default: `~/.cache/bamcp`), `BAMCP_CACHE_TTL` (default: 86400 seconds / 24 hours)

## Genome Build Detection

The `list_contigs` tool auto-detects genome build (GRCh37 or GRCh38) by comparing chr1 length:

| Build | chr1 Length | Notes |
|-------|-------------|-------|
| GRCh38 | 248,956,422 | Current human reference |
| GRCh37 | 249,250,621 | Older build, still common |

When no local reference is configured, the tool suggests public UCSC FASTA URLs that pysam can use remotely:
- GRCh38: `https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz`
- GRCh37: `https://hgdownload.soe.ucsc.edu/goldenPath/hg19/bigZips/hg19.fa.gz`

**Workflow**: Call `list_contigs` first on new BAM files to detect the build, then use the suggested reference URL or ask the user for their preferred reference.

## Commands

```bash
make install      # pip install -e ".[dev]"
make test         # pytest (unit + integration, ignores e2e)
make test-e2e     # playwright install chromium && pytest tests/e2e/
make lint         # ruff format --check && ruff check
make format       # ruff format && ruff check --fix
make typecheck    # mypy src
make docker-build # docker compose --profile prod build && --profile dev build
make docker-test  # docker compose --profile dev run --rm test && lint
make coverage     # pytest with --cov-report=html
make clean        # remove build artifacts
```

## Test Conventions

- Markers: `@pytest.mark.unit`, `@pytest.mark.integration`, `@pytest.mark.e2e`
- Async tests use `pytest-asyncio` with `asyncio_mode = "auto"`
- Fixtures in `tests/fixtures/` (small.bam, ref.fa, empty.bam)
- Fixture generator: `tests/create_fixtures.py` (creates test BAM/FASTA with known reads)
- HTTP mocking: `pytest-httpx` for ClinVar/gnomAD API tests
- Coverage threshold: 80% (fail_under in pyproject.toml)
- E2E tests use Playwright (sync API, run separately from async tests)

## Docker

- `Dockerfile` — production image (slim, no test tooling)
- `Dockerfile.dev` — development image (includes pytest, playwright, ruff, mypy)
- `docker-compose.yml` — profiles: `dev` (test+lint services), `beta`, `prod`
- `docker/entrypoint.sh` — entrypoint script
- `docker/healthcheck.py` — health check (imports bamcp, verifies server creation)

## Project Structure

```
src/bamcp/
  __init__.py, __main__.py, server.py, tools.py, parsers.py,
  clinvar.py, gnomad.py, cache.py, reference.py,
  resources.py, config.py, auth.py, static/viewer.html
tests/
  conftest.py, create_fixtures.py, fixtures/,
  test_parsers.py, test_tools.py, test_server.py, test_config.py,
  test_auth.py, test_resources.py, test_integration.py, test_docker.py,
  test_clinvar.py, test_gnomad.py, test_cache.py, test_reference.py,
  e2e/conftest.py, e2e/test_viewer_e2e.py
docker/
  entrypoint.sh, healthcheck.py
.github/
  workflows/ci.yml, workflows/release.yml,
  pull_request_template.md, ISSUE_TEMPLATE/
```

## Viewer UI Build Process

The interactive viewer is a TypeScript/Vite application in `src/bamcp/static/`.

**Build steps:**
```bash
cd src/bamcp/static
npm install          # Install dependencies (first time only)
npm run build        # Compile TypeScript + bundle with Vite → dist/viewer.html
```

**Development:**
```bash
npm run dev          # Start Vite dev server with hot reload
```

**CRITICAL - What NOT to do:**
- NEVER manually edit files in `src/bamcp/static/dist/` — always edit source files and rebuild
- NEVER copy files into `dist/` manually — Vite bundles everything automatically
- NEVER use `cat` or `echo` to modify dist files — use `npm run build`

**Source files:**
- `viewer.html` — HTML template with UI structure
- `client.ts` — MCP Apps SDK client wrapper (BAMCPClient)
- `mcp-app.ts` — Main application logic (BAMCPViewer)
- `renderer.ts` — Canvas-based read/coverage/variant rendering
- `state.ts` — Viewport and settings state management
- `types.ts` — TypeScript type definitions
- `constants.ts` — Color schemes and constants

**How bundling works:**
- Vite inlines the MCP Apps SDK into the HTML
- `get_viewer_html()` in `resources.py` serves `dist/viewer.html` if it exists
- Falls back to source `viewer.html` (won't work in sandboxed iframe)

## MCP Apps SDK Integration

The viewer uses `@modelcontextprotocol/ext-apps` to communicate with the MCP host.

**Key SDK methods:**
| Method | Purpose | Return |
|--------|---------|--------|
| `callServerTool()` | Directly invoke MCP server tool | `Promise<CallToolResult>` — result via Promise |
| `sendMessage()` | Send message to LLM (unreliable for tool calls) | `Promise<{isError?}>` |
| `updateModelContext()` | Update model context with viewer state | `Promise<{}>` |

**Common mistake:**
```typescript
// WRONG: Assumes result comes via callback
await app.callServerTool({ name: 'tool', arguments: {} });
// ontoolresult will NOT fire for callServerTool

// CORRECT: Use the Promise return value
const result = await app.callServerTool({ name: 'tool', arguments: {} });
if (result.structuredContent) { handleData(result.structuredContent); }
```

**Notification callbacks (`ontoolresult`, `ontoolinput`):**
- These fire for tool calls initiated by the HOST (LLM calling tools)
- They do NOT fire for `callServerTool` calls initiated by the APP

## Key Features & Test Cases

### 1. Region Visualization (`visualize_region` tool)
- **Feature**: Render aligned reads with coverage track, reference sequence, and variant calls
- **Compact mode**: Regions >500bp auto-omit sequences to reduce payload size
- **Test**: `tests/test_tools.py::test_visualize_region*`, `tests/e2e/test_viewer_e2e.py`

### 2. Auto-fetch Sequences on Zoom
- **Feature**: When zoomed to scale ≥10 (base-level), auto-fetch sequences if missing
- **Trigger**: `checkAndRequestSequences()` in `mcp-app.ts`
- **Uses**: `callServerTool('visualize_region')` for direct tool invocation
- **Fallback**: "Load Detail" button after 3s timeout
- **Test**: Manual — zoom into large region, bases should appear automatically

### 3. Variant Detection & Evidence
- **Feature**: Call SNVs/indels with VAF, depth, strand bias, MAPQ histograms
- **Artifact risk**: Position-in-read bias, strand bias detection
- **Test**: `tests/test_parsers.py::test_variant*`, `tests/test_tools.py::test_get_variant*`

### 4. External Database Lookups
- **ClinVar**: Clinical significance via NCBI E-utilities
- **gnomAD**: Population allele frequency via GraphQL API
- **Test**: `tests/test_clinvar.py`, `tests/test_gnomad.py` (uses `pytest-httpx` mocking)

### 5. Genome Build Detection
- **Feature**: Auto-detect GRCh37/GRCh38 from chr1 length
- **Suggests**: Public UCSC reference FASTA URLs
- **Test**: `tests/test_reference.py`

### 6. Remote BAM Support
- **Feature**: Stream BAM/CRAM from HTTP/S3 URLs with cached index files
- **Cache**: Session-isolated, 24h TTL by default
- **Test**: `tests/test_cache.py`

## Archived Planning Documents

Located in `archived/` folder:

- `BAMCP_Strategy.md` — MCP Apps architecture, viewer design, ClinVar/gnomAD integration roadmap
- `BAMCP_EVAL_HARNESS.md` — LLM genomic reasoning evaluation framework
- `BAMCP_IMPLEMENTATION_PLAN.md` — Original implementation schedule
