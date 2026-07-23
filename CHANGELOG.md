# Changelog

All notable changes to BAMCP are documented here. The format is based on
[Keep a Changelog](https://keepachangelog.com/en/1.1.0/), and this project aims to follow
[Semantic Versioning](https://semver.org/spec/v2.0.0.html) once it reaches a first tagged release.

## [Unreleased]

Response to an external engineering review — a phased hardening of the security boundary,
evaluation credibility, and frontend quality. See `docs/engineering-review-response.md` for the
plan and `docs/ENGINEERING_RETROSPECTIVE.md` for the write-up.

### Security

- **Reference FASTA URL is now SSRF-validated** — the `reference` argument was handed to pysam
  unvalidated, so a remote reference could target a private/internal host or the cloud metadata
  endpoint. All data sources (BAM/CRAM/VCF/reference/index) now pass through one policy.
- **Retired the OAuth theater** — dynamic client registration is off by default (an open
  `/register` that mints tokens is not access control); access is a stateless service token.
  Scopes are required when auth is on. Anonymous → 401 invariant preserved.
- **Proxy-aware, bounded rate limiter** — `X-Forwarded-For` is honored only behind a trusted
  proxy (no IP spoofing) and the tracking table is a bounded LRU (no memory exhaustion).

### Fixed

- Viewer lost the reference genome on pan/zoom/detail refetches (mismatch evidence vanished,
  CRAM failed to decode after navigation) — the full data-source context is now preserved.
- Negative ClinVar/gnomAD results (variant-not-found) were re-fetched every lookup — a distinct
  cache-miss sentinel now caches absences.
- Region cache and the services registry were unbounded — both are bounded now; the registry
  also guards against `id()` reuse and removes entries on close.
- pysam parses ran on the event loop's default executor, so a hung parse could starve unrelated
  async work — parsing is confined to a dedicated bounded pool.

### Added

- `schema_version` stamp on the viewer payload + a server↔viewer compatibility test.
- Eval grader scores **tool selection and answer correctness separately** (a case can no longer
  pass on the tool call alone); nightly runs both truth sets and both LLM providers.
- Frontend static/test layer: `tsc --noEmit` (strict), ESLint, and vitest, wired into CI.

### Changed

- Docs lead remote deployment with **Streamable HTTP** (SSE is deprecated by the MCP spec).
- Test suite is hermetic by default (live-network tests behind a `network` marker); CI installs
  from the frozen `uv.lock`.
