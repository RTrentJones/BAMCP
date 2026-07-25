# BAMCP

**Interactive BAM/CRAM variant visualization for AI assistants via the Model Context Protocol**

[![CI](https://github.com/RTrentJones/BAMCP/actions/workflows/ci.yml/badge.svg)](https://github.com/RTrentJones/BAMCP/actions/workflows/ci.yml)
[![Security](https://github.com/RTrentJones/BAMCP/actions/workflows/security.yml/badge.svg)](https://github.com/RTrentJones/BAMCP/actions/workflows/security.yml)
[![MCP](https://img.shields.io/badge/MCP-Apps%20Extension-blue)](https://modelcontextprotocol.io)
[![Python](https://img.shields.io/badge/python-3.10+-green)](https://python.org)
[![License](https://img.shields.io/badge/license-MIT-blue.svg)](LICENSE)

BAMCP brings IGV-style alignment visualization directly into your AI conversations. Browse BAM/CRAM files, inspect variant evidence, and navigate genomic regions—all through natural language interaction with Claude or other MCP-enabled assistants.

<p align="center">
  <img src="docs/assets/demo.png" alt="BAMCP Demo" width="800">
</p>

---

## Table of Contents

- [Why BAMCP?](#why-bamcp)
- [Features](#features)
- [Architecture](#architecture)
- [Installation](#installation)
- [Usage](#usage)
- [Tools Reference](#tools-reference)
- [Configuration](#configuration)
- [Implementation Details](#implementation-details)
- [Development](#development)
- [Roadmap](#roadmap)
- [Contributing](#contributing)
- [License](#license)

---

## Why BAMCP?

Existing genomics MCP servers (bio-mcp-samtools, AWS HealthOmics MCP) provide command-line operations but return text/JSON. **BAMCP is the first genomics MCP server to leverage the [MCP Apps Extension](https://github.com/modelcontextprotocol/ext-apps) to render interactive visualizations inline** — verifiable against the known genomics-MCP landscape above.

| Feature | CLI-based MCPs | BAMCP |
|---------|---------------|-------|
| Read BAM/CRAM files | ✅ | ✅ |
| Coverage statistics | ✅ | ✅ |
| **Interactive visualization** | ❌ | ✅ |
| **Variant browsing UI** | ❌ | ✅ |
| **Pan/zoom navigation** | ❌ | ✅ |
| **Inline in chat** | ❌ | ✅ |

---

## Features

- **Alignment Viewer** — Visualize reads with color-coded mismatches, insertions, deletions, and soft clips
- **Coverage Track** — Real-time depth of coverage across the viewing window
- **Variant Highlighting** — Automatic detection and highlighting of positions with non-reference alleles
- **VCF/BCF Overlay** — Use caller-produced variants as a primary/overlay source, including sample genotypes and structural-variant metadata
- **Region Navigation** — Jump to coordinates, genes, or specific variants
- **Format Support** — BAM, CRAM (with reference), and indexed remote files (HTTP/S3)
- **Canvas Rendering** — High-performance visualization of thousands of reads via HTML5 Canvas
- **Multiple Transports** — stdio, Streamable HTTP (preferred for network), and legacy SSE
- **OAuth 2.0 Authentication** — Optional token-based auth for network deployments
- **Docker Ready** — Production and development Docker images with compose profiles

---

## Architecture

```
┌─────────────────────────────────────────────────────────────────────────────┐
│                        MCP Client (Claude, Cursor, etc.)                    │
└─────────────────────────────────┬───────────────────────────────────────────┘
                                  │
                                  │ MCP Protocol (JSON-RPC over stdio/SSE/HTTP)
                                  ▼
┌─────────────────────────────────────────────────────────────────────────────┐
│                       BAMCP MCP Server (Python + FastMCP)                    │
│                                                                             │
│  ┌─────────────────┐  ┌─────────────────┐  ┌─────────────────────────────┐  │
│  │                 │  │                 │  │                             │  │
│  │  pysam Layer    │  │  Tool Handlers  │  │  UI Resource Provider       │  │
│  │                 │  │                 │  │                             │  │
│  │  • AlignmentFile│  │  • visualize    │  │  • ui://bamcp/viewer        │  │
│  │  • fetch()      │  │  • get_variants │  │  • text/html+mcp            │  │
│  │  • pileup()     │  │  • lookup_*     │  │  • Sandboxed iframe content │  │
│  │  • CRAM support │  │  • search_gene  │  │                             │  │
│  │                 │  │  • +6 more      │  │                             │  │
│  └─────────────────┘  └─────────────────┘  └─────────────────────────────┘  │
│                                                                             │
│  ┌─────────────────┐  ┌─────────────────────────────────────────────────┐  │
│  │  Auth (opt-in)  │  │  Config (env vars)                              │  │
│  │  • OAuth 2.0 AS │  │  • Transport, host, port                        │  │
│  │  • Token mgmt   │  │  • Reference path, read limits                  │  │
│  │  • Client reg   │  │  • Auth settings                                │  │
│  └─────────────────┘  └─────────────────────────────────────────────────┘  │
└─────────────────────────────────────────────────────────────────────────────┘
```

### Data Flow

```
User: "Show me reads at chr17:7577000-7577500 in tumor.bam"
                              │
                              ▼
                    ┌─────────────────┐
                    │  MCP Client     │
                    │  (Claude)       │
                    └────────┬────────┘
                             │ tools/call: visualize_region
                             ▼
                    ┌─────────────────┐
                    │  BAMCP Server   │
                    │                 │
                    │  1. Parse region│
                    │  2. pysam.fetch │
                    │  3. Serialize   │
                    │  4. Return UI   │
                    └────────┬────────┘
                             │ CallToolResult with ui/resourceUri
                             ▼
                    ┌─────────────────┐
                    │  MCP Client     │
                    │  renders iframe │
                    └────────┬────────┘
                             │ postMessage: init with read data
                             ▼
                    ┌─────────────────┐
                    │  Canvas UI      │
                    │                 │
                    │  1. Parse reads │
                    │  2. Pack rows   │
                    │  3. Render      │
                    └─────────────────┘
```

---

## Installation

### Prerequisites

- Python 3.10+
- [pysam](https://pysam.readthedocs.io/) (requires htslib)

### From PyPI (coming soon)

```bash
pip install bamcp
```

### From Source

```bash
git clone https://github.com/RTrentJones/BAMCP.git
cd BAMCP
python -m venv venv
source venv/bin/activate  # or `venv\Scripts\activate` on Windows
pip install -e ".[dev]"
```

### Docker

```bash
# Build production image
docker compose --profile prod build

# Run production server (Streamable HTTP on port 8000). The prod profile enables auth, so it
# needs a service token — the server refuses to start without one:
BAMCP_VERIFY_TOKEN="$(openssl rand -hex 32)" docker compose --profile prod up

# Prefer to run without auth? Use the beta profile (auth disabled) instead:
docker compose --profile beta up

# Run tests in Docker
docker compose --profile dev run --rm test
```

### MCP Client Configuration

#### Claude Desktop (stdio)

Add to `~/.config/claude/claude_desktop_config.json`:

```json
{
  "mcpServers": {
    "bamcp": {
      "command": "python",
      "args": ["-m", "bamcp"],
      "env": {
        "BAMCP_REFERENCE": "/path/to/hg38.fa"
      }
    }
  }
}
```

#### Claude Desktop (Streamable HTTP — network, preferred)

Run the server with `BAMCP_TRANSPORT=streamable-http` and point the client at `/mcp`:

```json
{
  "mcpServers": {
    "bamcp": {
      "url": "http://localhost:8000/mcp"
    }
  }
}
```

> Legacy SSE (`BAMCP_TRANSPORT=sse`, URL `http://localhost:8000/sse`) is still supported for
> compatibility, but the MCP spec deprecates the HTTP+SSE transport in favor of Streamable HTTP.
> Prod runs Streamable HTTP.

#### Cursor / VS Code

> **Note:** `uvx bamcp` requires the PyPI release, which is not yet published. Until
> then, run from a local checkout — `"command": "uv"`, `"args": ["run", "bamcp"]` with
> `"cwd"` set to your clone — or use the stdio `python -m bamcp` form above.

```json
{
  "mcpServers": {
    "bamcp": {
      "command": "uvx",
      "args": ["bamcp"]
    }
  }
}
```

---

## Usage

### Browse a Region

> "Show me the reads at chr17:7577000-7577500 in /data/tumor.bam"

The assistant will call `visualize_region` and render an interactive alignment viewer inline.

### Jump to a Position

> "Jump to position 7674220 on chr17 in my BAM"

Centers the viewer on a specific position with a configurable window size.

### Inspect a Variant

> "What's the evidence for the variant at chr7:140453136 in my sample?"

BAMCP will center on the position and highlight reads supporting reference vs. alternate alleles.

### Get Coverage Stats

> "What's the average coverage across TP53 in my exome?"

Returns coverage statistics without visualization.

### List Available Chromosomes

> "What contigs are in this BAM file?"

---

## Tools Reference

| Tool | Description | Required Args | Optional Args |
|------|-------------|---------------|---------------|
| `visualize_region` | View aligned reads with interactive UI | `file_path`, `region` | `reference`, `vcf_path`, `variant_source` |
| `get_variants` | Detect and return candidate variants in a region | `file_path`, `region` | `reference`, `vcf_path`, `variant_source`, `min_vaf`, `min_depth` |
| `get_coverage` | Calculate depth statistics | `file_path`, `region` | `reference` |
| `list_contigs` | List chromosomes and detect genome build | `file_path` | `reference` |
| `jump_to` | Jump to a specific genomic position | `file_path`, `position` | `contig`, `window`, `reference`, `vcf_path`, `variant_source` |
| `get_region_summary` | Text summary for LLM reasoning (no UI) | `file_path`, `region` | `reference`, `vcf_path`, `variant_source` |
| `lookup_clinvar` | Look up variant in ClinVar | `chrom`, `pos`, `ref`, `alt` | — |
| `lookup_gnomad` | Look up variant in gnomAD | `chrom`, `pos`, `ref`, `alt` | — |
| `get_variant_curation_summary` | Detailed curation with artifact risk | `file_path`, `chrom`, `pos`, `ref`, `alt` | `window`, `reference` |
| `search_gene` | Search gene by symbol (NCBI); returns build-consistent coordinates | `symbol` | `file_path`, `build` |
| `cleanup_cache` | Clean up session's index cache | — | — |

> **Genome build & `search_gene`.** Gene coordinates are build-specific — a GRCh38 locus in a
> GRCh37 BAM is silently ~0.5–2 Mb off the real gene. Pass the BAM's **`file_path`** to
> `search_gene` and it auto-detects the build from the header and returns matching coordinates
> (or set **`build`** to override). The response states the coordinates' `genome_build` and warns
> on any mismatch. Without a `file_path`, coordinates use the server default (`BAMCP_GENOME_BUILD`,
> default GRCh38) — confirm your BAM's build with `list_contigs` first.

#### VCF as the primary variant source

BAMCP is at its best as an **evidence viewer** for a trusted caller's variants rather than
an ad-hoc caller of its own. Pass `vcf_path` (VCF/BCF) plus `variant_source` to control the
relationship:

- **`variant_source="vcf"`** — the VCF is the authoritative variant set. BAMCP does *not* add
  its own candidates; it reads the BAM and attaches read-level evidence (strand balance, mean
  quality, artifact risk, confidence) to each VCF site so you can review support. Requires `vcf_path`.
- **`variant_source="auto"`** (default) — BAMCP's read-level candidates with the VCF overlaid and
  de-duplicated.
- **`variant_source="bamcp"`** — BAMCP's read-level candidates only; any `vcf_path` is ignored.

Enrich further with `lookup_clinvar` / `lookup_gnomad` on the returned 1-based positions.

### Region Format

Regions can be specified as:
- `chr1:1000-2000` — Standard format
- `chr1:1,000-2,000` — Commas allowed
- `1:1000-2000` — Without "chr" prefix (depends on BAM header)

---

## Configuration

### Core Settings

| Environment Variable | Description | Default |
|---------------------|-------------|---------|
| `BAMCP_REFERENCE` | Path to reference FASTA (required for CRAM) | None |
| `BAMCP_MAX_READS` | Maximum reads to fetch per region | `10000` |
| `BAMCP_DEFAULT_WINDOW` | Default viewing window size (bp) | `500` |
| `BAMCP_MIN_VAF` | Minimum allele fraction for candidate variants | `0.02` |
| `BAMCP_MIN_DEPTH` | Minimum depth for candidate variants | `2` |
| `BAMCP_MIN_MAPQ` | Minimum mapping quality filter | `0` |
| `BAMCP_MIN_BASEQ` | Minimum base quality for coverage and candidate variant counts | `0` |
| `BAMCP_MAX_REMOTE_INDEX_BYTES` | Maximum remote BAM/CRAM index download size | `67108864` |

### Transport Settings

| Environment Variable | Description | Default |
|---------------------|-------------|---------|
| `BAMCP_TRANSPORT` | Transport mode: `stdio`, `sse`, `streamable-http` | `stdio` |
| `BAMCP_HOST` | HTTP bind address | `0.0.0.0` |
| `BAMCP_PORT` | HTTP port | `8000` |

### Authentication Settings

| Environment Variable | Description | Default |
|---------------------|-------------|---------|
| `BAMCP_AUTH_ENABLED` | Enable auth (`true`/`false`) | `false` |
| `BAMCP_ISSUER_URL` | OAuth issuer URL | `http://localhost:8000` |
| `BAMCP_RESOURCE_SERVER_URL` | OAuth resource server URL | `http://localhost:8000` |
| `BAMCP_REQUIRED_SCOPES` | Comma-separated required scopes (defaults to `bamcp:read` when auth is on) | None |
| `BAMCP_TOKEN_EXPIRY` | Access token lifetime (seconds) | `3600` |
| `BAMCP_ALLOW_DYNAMIC_REGISTRATION` | Allow OAuth dynamic client registration | `false` |
| `BAMCP_VERIFY_TOKEN` | M2M service token (a single bearer accepted statelessly) | None |

**Trust model.** BAMCP is a single-instance server, not a multi-tenant identity provider. Access
is a **service token** (`BAMCP_VERIFY_TOKEN`) — a stateless bearer that survives restarts.
Dynamic client registration is **off by default**: an open `/register` endpoint would let any
client self-register and mint a token, which is not access control. Enable it
(`BAMCP_ALLOW_DYNAMIC_REGISTRATION=true`) only for interactive OAuth clients in a trusted/dev
context. With auth enabled, an unauthenticated request is rejected (401) and only the service
token (or an out-of-band-registered client) is authorized.

### Security Settings

| Environment Variable | Description | Default |
|---------------------|-------------|---------|
| `BAMCP_ALLOW_REMOTE_FILES` | Allow HTTP(S) BAM/CRAM/VCF/reference URLs | `false` |
| `BAMCP_ALLOWED_REMOTE_HOSTS` | Comma-separated host allow-list for remote fetch (the real SSRF control) | None |
| `BAMCP_ALLOWED_DIRECTORIES` | Comma-separated local dirs that file paths must be under | None |
| `BAMCP_RATE_LIMIT` | Requests per minute per client IP | `60` |
| `BAMCP_RATE_LIMIT_TRUSTED_PROXIES` | IPs/CIDRs whose `X-Forwarded-For` is trusted | loopback |
| `BAMCP_TRUSTED_HOSTS` | Allowed `Host` headers (DNS-rebinding protection) | None |

Remote fetch is SSRF-guarded: private/internal/metadata IPs are always blocked, and — because the
actual fetch is performed by pysam/htslib (which follows redirects outside BAMCP) — the
**host allow-list is the primary control**. Set `BAMCP_ALLOWED_REMOTE_HOSTS` in any
internet-reachable deployment. `X-Forwarded-For` is honored only from a trusted proxy so a client
cannot spoof its IP to bypass rate limits.

---

## Implementation Details

### Why Canvas over React?

React's virtual DOM reconciliation struggles with thousands of elements. Genomic regions typically contain 1,000-10,000+ reads:

| Approach | Max Reads | DOM Nodes | Frame Rate |
|----------|-----------|-----------|------------|
| React/SVG | ~500 | O(n) | Degrades |
| Canvas 2D | ~10,000 | O(1) | 60fps |
| WebGL | ~100,000+ | O(1) | 60fps |

BAMCP uses Canvas 2D for the rendering layer with vanilla JS for controls.

### Read Packing Algorithm

Reads are packed into rows using a greedy algorithm:

```
For each read (sorted by start position):
    Find first row where read.start > row.end
    If found: add to that row
    Else: create new row
```

This minimizes vertical space while preventing overlaps.

### Variant Detection

BAMCP performs simple pileup-based variant detection:

1. Count bases at each position from aligned reads
2. Calculate variant allele frequency (VAF) = alt_count / total_depth
3. Report positions where VAF >= threshold and depth >= minimum

This is intentionally basic—for production variant calling, use dedicated tools (GATK, DeepVariant, etc.).

### CRAM Support

CRAM files require a reference FASTA for decoding. Set `BAMCP_REFERENCE` or pass `reference` to each tool call.

---

## Development

### Project Structure

```
BAMCP/
├── src/
│   └── bamcp/
│       ├── __init__.py
│       ├── __main__.py       # Entry point (transport selection)
│       ├── server.py         # FastMCP server setup
│       ├── tools.py          # Tool handlers
│       ├── parsers.py        # pysam wrappers
│       ├── resources.py      # UI resource provider
│       ├── config.py         # Environment-based configuration
│       ├── auth.py           # OAuth 2.0 authorization server
│       └── static/
│           └── viewer.html   # Bundled alignment viewer UI
├── tests/
│   ├── conftest.py
│   ├── create_fixtures.py    # Generate test BAM/FASTA files
│   ├── test_parsers.py
│   ├── test_tools.py
│   ├── test_server.py
│   ├── test_config.py
│   ├── test_auth.py
│   ├── test_resources.py
│   ├── test_integration.py
│   ├── test_docker.py
│   ├── fixtures/
│   └── e2e/
│       ├── conftest.py
│       └── test_viewer_e2e.py
├── docker/
│   ├── entrypoint.sh
│   └── healthcheck.py
├── .github/
│   └── workflows/
│       ├── ci.yml
│       ├── greenlight-build.yml   # ship-gate + build/push image + trigger deploy
│       ├── eval-nightly.yml
│       └── security.yml
├── Dockerfile              # Production (multi-stage)
├── Dockerfile.dev          # Development (with test tooling)
├── docker-compose.yml      # dev/beta/prod profiles
├── pyproject.toml
├── Makefile
├── README.md
├── CONTRIBUTING.md
├── SECURITY.md
└── LICENSE
```

### Running Locally

```bash
# stdio mode (default)
python -m bamcp

# SSE mode (network)
BAMCP_TRANSPORT=sse python -m bamcp

# Streamable HTTP mode
BAMCP_TRANSPORT=streamable-http python -m bamcp

# With MCP inspector
npx @modelcontextprotocol/inspector python -m bamcp
```

### Running Tests

```bash
# Unit + integration tests
make test

# With HTML coverage report
make coverage

# E2E tests (requires Playwright)
make test-e2e

# Tests in Docker
make docker-test
```

### Code Style

```bash
# Check
make lint

# Auto-format
make format

# Type check
make typecheck
```

Tools: [ruff](https://github.com/astral-sh/ruff) (linting + formatting), [mypy](https://mypy.readthedocs.io/) (type checking).

### Docker Development

```bash
# Build all images
make docker-build

# Run tests in container
docker compose --profile dev run --rm test

# Run linting in container
docker compose --profile dev run --rm lint

# Start beta server
docker compose --profile beta up

# Start production server (with auth)
docker compose --profile prod up
```

---

## Roadmap

- [x] **v0.1** — Core viewer with BAM support, Canvas renderer, variant detection
- [x] **v0.2** — CRAM support, coverage track, MCP Apps UI resource
- [x] **v0.3** — FastMCP migration, SSE/HTTP transports, OAuth 2.0
- [x] **v0.4** — Docker infrastructure (prod/dev/beta), CI/CD
- [x] **v0.5** — MCP Apps SDK integration (`updateModelContext`, display modes, auto-fetch)
- [x] **v0.6** — ClinVar + gnomAD API clients, gene search, variant curation tools
- [x] **v0.7** — Evaluation harness: deterministic ground-truth gate (variant P/R/F1, artifact-type recall) wired into CI, plus LLM + vision harness ([EVALS.md](EVALS.md)). GIAB real-data set scaffolded.
- [ ] **v0.8** — Gene annotation track (RefSeq)
- [ ] **v0.9** — VCF overlay, multi-sample comparison
- [ ] **v1.0** — Stable release

---

## Related Projects

- [bio-mcp-samtools](https://github.com/bio-mcp/bio-mcp-samtools) — CLI-based samtools MCP wrapper
- [igv.js](https://github.com/igvteam/igv.js) — JavaScript genome visualization library
- [pysam](https://github.com/pysam-developers/pysam) — Python interface to htslib
- [MCP Apps Extension](https://github.com/modelcontextprotocol/ext-apps) — UI extension for MCP

---

## Contributing

Contributions are welcome! Please see [CONTRIBUTING.md](CONTRIBUTING.md) for guidelines.

1. Fork the repository
2. Create a feature branch (`git checkout -b feature/amazing-feature`)
3. Commit your changes (`git commit -m 'Add amazing feature'`)
4. Push to the branch (`git push origin feature/amazing-feature`)
5. Open a Pull Request

---

## Security & Safety

BAMCP reads identifiable genomic data and feeds an LLM, so both are first-class:

- **[SECURITY.md](SECURITY.md)** — threat model, in-tree controls (SSRF
  filtering, opt-in remote files, OAuth, rate limiting, container hardening),
  and automated scanning (`pip-audit` + CodeQL in CI).
- **[SAFETY.md](SAFETY.md)** — intended-use boundary (research/visualization,
  **not** a diagnostic device) and the safety invariants enforced by the eval
  gate, including a guard against overconfident calls on artifact-prone sites.

---

## License

MIT License — see [LICENSE](LICENSE) for details.

---

## Acknowledgments

- The [Model Context Protocol](https://modelcontextprotocol.io) team at Anthropic
- [MCP-UI](https://github.com/idosal/mcp-ui) for pioneering interactive MCP interfaces
- The [pysam](https://github.com/pysam-developers/pysam) and [htslib](https://github.com/samtools/htslib) maintainers
- [IGV](https://igv.org/) for inspiration on genomics visualization

---

<p align="center">
  <b>BAMCP</b> — BAM files + Model Context Protocol<br>
  <sub>Built for the computational biology community</sub>
</p>
