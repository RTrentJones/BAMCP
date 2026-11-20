# Security Policy

## Supported Versions

| Version | Supported |
| ------- | --------- |
| 0.1.x   | Yes       |

## Reporting a Vulnerability

If you discover a security vulnerability, please report it responsibly:

1. **Do not** open a public GitHub issue.
2. **Preferred:** open a private report via GitHub Security Advisories —
   <https://github.com/RTrentJones/BAMCP/security/advisories/new>. This keeps the
   report private and lets us collaborate on a fix and CVE in-thread.
3. **Alternative:** email <rtrentjones@gmail.com> with a description of the
   vulnerability and reproduction steps.
4. Allow reasonable time for a fix before public disclosure.

We acknowledge receipt within 48 hours and aim to release a patch within 7
days for critical issues.

---

## Threat Model

BAMCP reads alignment files (BAM/CRAM) — potentially from remote URLs — and
exposes tools to an LLM over MCP, optionally over the network (SSE /
streamable-HTTP transports). The primary risks are therefore:

- **Server-Side Request Forgery (SSRF)** via attacker-influenced remote file
  URLs reaching internal services or cloud metadata endpoints.
- **Path traversal / arbitrary file read** via local file-path arguments.
- **DNS rebinding** against the HTTP transports.
- **Resource exhaustion** from large regions / high request volume.
- **Sensitive-data exposure** — alignment files can contain identifiable human
  genomic data (see [SAFETY.md](SAFETY.md) for the data-handling stance).

## Controls

These are implemented in-tree; the module is named for each.

| Control | Where | Notes |
| ------- | ----- | ----- |
| Scheme classification | `core/validation.py::_is_remote_source` | Every caller-supplied BAM/CRAM, VCF and reference passes through it. Only an exact `http://` / `https://` prefix is remote; any other value that starts like a URI scheme (`http:/` with one slash, `data:`, `ftp:`, `file:`, `s3:`, …) is rejected, as is a local path embedding `://` or `##idx##`, so no other htslib URL handler is reachable. A local VCF/BCF must exist. |
| SSRF prevention | `core/validation.py::validate_remote_url` | http/https on ports 80/443 only. Rejects `#` (htslib `##idx##` smuggling), userinfo, and backslash/control characters. Resolves the host and blocks every non-global address (loopback, RFC1918, link-local/metadata, CGNAT, multicast, reserved; IPv4-mapped IPv6 judged by the embedded IPv4). Redirects inside htslib can't be checked per-URL, so the exact-host allow-list `BAMCP_ALLOWED_REMOTE_HOSTS` is the real SSRF control; startup warns when remote files are enabled and it is unset. |
| Path validation | `core/validation.py::validate_path` | When `BAMCP_ALLOWED_DIRECTORIES` is set, resolved paths must fall inside it (unset means no directory restriction). `BAMCP_ALLOW_LOCAL_FILES=false` rejects every caller-supplied local path before any filesystem access (the configured `BAMCP_REFERENCE` is exempt). Only `.bam`/`.cram` extensions; paths ≤2048 chars (regions ≤100 chars, via `validate_region`). |
| Input validation | `core/validation.py` | Region/variant strings checked against strict patterns before reaching pysam; `search_gene` symbols must match `GENE_SYMBOL_PATTERN` (1–40 letters, digits, `.`, `_`, `-`) before they enter an NCBI Entrez query. |
| CRAM reference lookup | `__main__.py::_configure_htslib_reference_lookup` | Unless the operator sets them, `REF_PATH` / `REF_CACHE` point at a local directory under the cache dir, so htslib never looks a CRAM reference up over the network by M5 checksum. An operator value naming a URL is kept and logged with a warning. |
| CRAM header references | `core/parsers.py::_open_alignment` (`guard_cram_decoding`, `_check_cram_decode`) | htslib decodes a contig against the file in its `@SQ UR:` tag when no reference is given or the reference lacks that contig. With `BAMCP_ALLOW_LOCAL_FILES=false`, a CRAM needs a resolved reference covering every `UR:`-tagged contig. Checked on every alignment open by the detected format, not the file name (a CRAM named `.bam` is still gated). htslib refuses remote `UR:` values itself; a test pins that. |
| DNS-rebinding protection | `__main__.py` (`TrustedHostMiddleware`) | Configured via `BAMCP_TRUSTED_HOSTS`. |
| Security headers | `middleware/security.py` | X-Content-Type-Options, X-Frame-Options, CSP, etc. |
| Request body limit | `middleware/bodylimit.py` | 413 above `BAMCP_MAX_REQUEST_BYTES` (default 1 MiB; 64 KiB on `/register`, `/token`, `/revoke`, `/authorize`), on the declared `Content-Length` and on streamed bytes; a malformed or conflicting `Content-Length` is a 400. Runs after the Host check and before the rate limiter. |
| Rate limiting | `middleware/ratelimit.py` | Sliding-window per-IP, default 120 req/min (`BAMCP_RATE_LIMIT`; every MCP JSON-RPC message counts). Client IP comes from `CF-Connecting-IP` / `X-Forwarded-For` only when the peer is in `BAMCP_RATE_LIMIT_TRUSTED_PROXIES` (default loopback); an IPv6 client is keyed by its /64; the tracking table is a bounded LRU. |
| AuthN/Z (opt-in) | `middleware/auth.py` | OAuth 2.0 authorization server; enable with `BAMCP_AUTH_ENABLED`. Accepts a stateless service token (`BAMCP_VERIFY_TOKEN`); dynamic client registration is off unless `BAMCP_ALLOW_DYNAMIC_REGISTRATION=true` (the hosted demo turns it on, with size-bounded in-memory stores and per-client caps: 10 `redirect_uris` of ≤2048 chars, `client_name` ≤256 chars, ≤16 KiB of metadata). Scopes are required when auth is on (default `bamcp:read`). |
| Container hardening | `docker-compose.yml` (prod) | `cap_drop: [ALL]`, `no-new-privileges`, `read_only`, memory cap, non-root user; port bound to `127.0.0.1`; `BAMCP_ALLOWED_REMOTE_HOSTS` set and `BAMCP_ALLOW_LOCAL_FILES=false`. |
| Error sanitization | `server.py::_run` → `core/errors.py::sanitize_tool_exception` | Every tool exception becomes an `isError` result with a client-safe message; config/path internals are never echoed. An `unavailable` ClinVar/gnomAD answer carries a fixed `reason_code`; upstream error text is logged, not returned. |

### Remote files are opt-in

Remote BAM/CRAM access is **disabled by default**. It must be explicitly
enabled with `BAMCP_ALLOW_REMOTE_FILES=true`, and even then every URL passes
through `validate_remote_url`. Prefer also setting `BAMCP_ALLOWED_REMOTE_HOSTS`
to an allowlist in any networked deployment. A remote-only public server
should also set `BAMCP_ALLOW_LOCAL_FILES=false`.

### Residual risks

- Redirects followed inside htslib/libcurl and DNS rebinding between BAMCP's
  check and htslib's own connection cannot be validated per URL. The
  exact-host allow-list (`BAMCP_ALLOWED_REMOTE_HOSTS`) is the real SSRF
  control; startup warns when remote files are enabled without it.
- An operator-set `REF_PATH` / `REF_CACHE` that names a URL is honoured (with a
  warning) and is fetched outside the allow-list.

## Automated Scanning (CI)

| Check | Workflow | Catches |
| ----- | -------- | ------- |
| `ruff` security lints (`S` ruleset) | `ci.yml::lint` | Insecure code patterns (subprocess, eval, hardcoded creds, etc.). |
| `pip-audit --strict` | `security.yml::pip-audit` | Known CVEs in resolved dependencies, on PRs and weekly. |
| CodeQL (`security-extended`) | `security.yml::codeql` | Data-flow vulnerabilities (injection, SSRF, path traversal) in `python` and `javascript-typescript` (the viewer). |
| Dependabot | `dependabot.yml` | Outdated/vulnerable pip, GitHub Actions, and Docker base-image (digest-pinned) versions. |

To accept a specific advisory that has no available fix, add
`--ignore-vuln <GHSA-id>` to the `pip-audit` step with a short justification and
keep the list reviewed.

## Dependency Pinning

`mcp>=1.23.0` is required (CVE-2025-66416, DNS-rebinding fix). Security-relevant
floors like this are documented at the point of use; do not lower them without
a corresponding note here.

## Deployment Notes

The production deployment runs on a private subnet with **no public IP**;
external access is via an outbound-only Cloudflare Tunnel (`cloudflared`
sidecar), so there are no inbound firewall rules or load balancers to harden.
TLS terminates at the Cloudflare edge.
