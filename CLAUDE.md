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
| `visualize_region` | View aligned reads with interactive MCP Apps visualization (returns UI + data) |
| `get_variants` | Detect and return variants in a region |
| `get_coverage` | Calculate depth of coverage statistics |
| `list_contigs` | List chromosomes/contigs and detect genome build (GRCh37/GRCh38) with suggested public reference URL |
| `jump_to` | Jump to a specific genomic position with configurable window |
| `get_region_summary` | Text-only region summary for LLM reasoning (no UI) |
| `lookup_clinvar` | Look up variant in ClinVar for clinical significance and conditions |
| `lookup_gnomad` | Look up variant in gnomAD for population allele frequency data |
| `get_variant_curation_summary` | Detailed curation summary with artifact risk assessment and recommendations |
| `search_gene` | Search for a gene by symbol and return genomic coordinates |
| `cleanup_cache` | Clean up expired BAM index cache files |

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

## Vision Documents

- `BAMCP_Strategy.md` — MCP Apps architecture, viewer design, ClinVar/gnomAD integration roadmap
- `BAMCP_EVAL_HARNESS.md` — LLM genomic reasoning evaluation framework (classify_variant, ACMG scaffolds, failure mode detection, ground truth benchmarks)
- `BAMCP_IMPLEMENTATION_PLAN.md` — Full implementation schedule for Phases 1-4
