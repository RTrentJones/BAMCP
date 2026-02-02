# BAMCP

Interactive BAM/CRAM variant visualization MCP server with MCP Apps UI.

## Architecture

```
FastMCP server (server.py)
  ├── Tool handlers (tools.py) → pysam parsers (parsers.py)
  ├── UI resource (resources.py → static/viewer.html)
  ├── Auth provider (auth.py) — opt-in OAuth 2.0
  └── Config (config.py) — all settings from env vars
```

**Entry point**: `src/bamcp/__main__.py` — reads `BAMCPConfig.from_env()`, creates FastMCP server, runs with selected transport.

## Key Modules

| Module | Role |
|--------|------|
| `server.py` | FastMCP setup, tool/resource registration, auth wiring |
| `tools.py` | `handle_browse_region`, `handle_get_variants`, `handle_get_coverage`, `handle_list_contigs`, `handle_jump_to` |
| `parsers.py` | `fetch_region()` — pysam BAM/CRAM parsing, read extraction, coverage, variant detection |
| `resources.py` | `get_viewer_html()` — serves bundled HTML from `static/viewer.html` |
| `config.py` | `BAMCPConfig` dataclass with `from_env()` classmethod |
| `auth.py` | `BAMCPAuthProvider` (in-memory OAuth 2.0), `build_auth_settings()` |

## Registered Tools

| Tool | Description |
|------|-------------|
| `browse_region` | View aligned reads with interactive visualization (returns UI + data) |
| `get_variants` | Detect and return variants in a region |
| `get_coverage` | Calculate depth of coverage statistics |
| `list_contigs` | List chromosomes/contigs in BAM/CRAM header |
| `jump_to` | Jump to a specific genomic position with configurable window |

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
  resources.py, config.py, auth.py, static/viewer.html
tests/
  conftest.py, create_fixtures.py, fixtures/,
  test_parsers.py, test_tools.py, test_server.py, test_config.py,
  test_auth.py, test_resources.py, test_integration.py, test_docker.py,
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
