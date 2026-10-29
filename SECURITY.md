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
| SSRF prevention | `core/validation.py::validate_remote_url` | http/https on ports 80/443 only. Rejects `#` (htslib `##idx##` smuggling), userinfo, and backslash/control characters. Resolves the host and blocks every non-global address (loopback, RFC1918, link-local/metadata, CGNAT, multicast, reserved; IPv4-mapped IPv6 judged by the embedded IPv4). Redirects inside htslib can't be checked per-URL, so the exact-host allow-list `BAMCP_ALLOWED_REMOTE_HOSTS` is the real SSRF control; startup warns when it is unset. |
| Path validation | `core/validation.py::validate_path` | When `BAMCP_ALLOWED_DIRECTORIES` is set, resolved paths must fall inside it (unset means no directory restriction). Only `.bam`/`.cram` extensions; paths ≤2048 chars, regions ≤100 chars. |
| Input validation | `core/validation.py` | Region/variant strings checked against strict patterns before reaching pysam. |
| DNS-rebinding protection | `__main__.py` (`TrustedHostMiddleware`) | Configured via `BAMCP_TRUSTED_HOSTS`. |
| Security headers | `middleware/security.py` | X-Content-Type-Options, X-Frame-Options, CSP, etc. |
| Rate limiting | `middleware/ratelimit.py` | Sliding-window per-IP, default 120 req/min (`BAMCP_RATE_LIMIT`; every MCP JSON-RPC message counts). Client IP comes from `CF-Connecting-IP` / `X-Forwarded-For` only when the peer is in `BAMCP_RATE_LIMIT_TRUSTED_PROXIES` (default loopback); the tracking table is a bounded LRU. |
| AuthN/Z (opt-in) | `middleware/auth.py` | OAuth 2.0 authorization server; enable with `BAMCP_AUTH_ENABLED`. Accepts a stateless service token (`BAMCP_VERIFY_TOKEN`); dynamic client registration is off unless `BAMCP_ALLOW_DYNAMIC_REGISTRATION=true` (the hosted demo turns it on, with size-bounded in-memory stores). Scopes are required when auth is on (default `bamcp:read`). |
| Container hardening | `docker-compose.yml` (prod) | `cap_drop: [ALL]`, `no-new-privileges`, `read_only`, memory cap, non-root user. |
| Error sanitization | `server.py::_run` → `core/errors.py::sanitize_tool_exception` | Every tool exception becomes an `isError` result with a client-safe message; config/path internals are never echoed. |

### Remote files are opt-in

Remote BAM/CRAM access is **disabled by default**. It must be explicitly
enabled with `BAMCP_ALLOW_REMOTE_FILES=true`, and even then every URL passes
through `validate_remote_url`. Prefer also setting `BAMCP_ALLOWED_REMOTE_HOSTS`
to an allowlist in any networked deployment.

## Automated Scanning (CI)

| Check | Workflow | Catches |
| ----- | -------- | ------- |
| `ruff` security lints (`S` ruleset) | `ci.yml::lint` | Insecure code patterns (subprocess, eval, hardcoded creds, etc.). |
| `pip-audit --strict` | `security.yml::pip-audit` | Known CVEs in resolved dependencies, on PRs and weekly. |
| CodeQL (`security-extended`) | `security.yml::codeql` | Data-flow vulnerabilities (injection, SSRF, path traversal). |
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
