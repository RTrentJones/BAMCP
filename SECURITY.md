# Security Policy

## Supported Versions

| Version | Supported |
| ------- | --------- |
| 0.1.x   | Yes       |

## Reporting a Vulnerability

If you discover a security vulnerability, please report it responsibly:

1. **Do not** open a public GitHub issue.
2. Email the maintainers with a description of the vulnerability.
3. Allow reasonable time for a fix before public disclosure.

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
| SSRF prevention | `core/validation.py::validate_remote_url` | Resolves hostnames and blocks RFC1918, loopback, link-local, and `169.254.0.0/16` (cloud metadata); IPv6 equivalents included. Optional host allowlist via `BAMCP_ALLOWED_REMOTE_HOSTS`. |
| Path validation | `core/validation.py::validate_path` | Restricts to `BAMCP_ALLOWED_DIRECTORIES`; only `.bam`/`.cram` extensions; length caps. |
| Input validation | `core/validation.py` | Region/variant strings checked against strict patterns before reaching pysam. |
| DNS-rebinding protection | `__main__.py` (`TrustedHostMiddleware`) | Configured via `BAMCP_TRUSTED_HOSTS`. |
| Security headers | `middleware/security.py` | X-Content-Type-Options, X-Frame-Options, CSP, etc. |
| Rate limiting | `middleware/ratelimit.py` | Sliding-window per-IP, default 60 req/min (`BAMCP_RATE_LIMIT`). |
| AuthN/Z (opt-in) | `middleware/auth.py` | OAuth 2.0 authorization server; enable with `BAMCP_AUTH_ENABLED`. |
| Container hardening | `docker-compose.yml` (prod) | `cap_drop: [ALL]`, `no-new-privileges`, `read_only`, memory cap, non-root user. |
| Error sanitization | tool handlers | Error messages avoid leaking config/path internals. |

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
| Dependabot | `dependabot.yml` | Outdated/vulnerable pip + GitHub Actions versions. |

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
