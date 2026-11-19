# Changelog

All notable changes to BAMCP are documented here. The format is based on
[Keep a Changelog](https://keepachangelog.com/en/1.1.0/), and this project aims to follow
[Semantic Versioning](https://semver.org/spec/v2.0.0.html) once it reaches a first tagged release.

## [Unreleased]

Response to an external engineering review — a phased hardening of the security boundary,
evaluation credibility, and frontend quality. See `docs/engineering-review-response.md` for the
plan and `docs/ENGINEERING_RETROSPECTIVE.md` for the write-up. Entries that cite a PR number
landed after that review.

### Security

- **Reference FASTA URL is now SSRF-validated** — the `reference` argument was handed to pysam
  unvalidated, so a remote reference could target a private/internal host or the cloud metadata
  endpoint. All data sources (BAM/CRAM/VCF/reference/index) now pass through one policy.
- **Retired the OAuth theater** — dynamic client registration is off by default (an open
  `/register` that mints tokens is not access control); access is a stateless service token.
  Scopes are required when auth is on. Anonymous → 401 invariant preserved. (The public prod
  demo turns registration back on in the wrapper's infra, `BAMCP_ALLOW_DYNAMIC_REGISTRATION`.)
- **Proxy-aware, bounded rate limiter** — `X-Forwarded-For` is honored only behind a trusted
  proxy (no IP spoofing) and the tracking table is a bounded LRU (no memory exhaustion).
- **Deep-review security fixes (#63)** — remote URLs containing `#` are rejected (htslib reads
  `x.bam##idx##<url>` as a second, unvalidated URL); index-download redirects are re-validated
  at every hop; every non-global address is blocked; DNS resolution runs off the event loop;
  the OAuth stores and caches are bounded; the rate limiter reads `CF-Connecting-IP` behind the
  tunnel (default 120 requests/min per IP); viewer strings are HTML-escaped (XSS).

### Fixed

- An indel in a homopolymer run was flagged from the run's reference length alone. With at
  least 8 reads spanning the whole run, the homopolymer risk now comes from the called
  length's share among those reads, which also give the allele fraction and a genotype hint;
  with fewer, the reference-length rule still applies (#71).
- `list_contigs` suggested UCSC `bigZips` FASTAs that pysam cannot open; suggestions are now
  verified bgzip-or-uncompressed URLs with an index, matched to the BAM's contig naming, or
  `None` with a `reference_note`. Parse timeouts carry the region and mode instead of an empty
  message, and an empty `BAMCP_REFERENCE=` reads as unset (#59).
- VCF header and AF parsing crashes; external lookups now run in the BAM's build; gnomAD
  not-found handling; `max_pop_af` gated so it cannot produce a false BA1; read-position bias
  corrected (#63).
- The Docker image pinned `BAMCP_MIN_VAF`/`_DEPTH`/`_MAPQ`, so prod ran `min_mapq=0` against a
  code default of 20; the image no longer overrides the calling thresholds (#68).
- Eval comparison runs now persist each spec's results, and Tracer floor verdicts are computed
  from raw counts rather than rounded report values.
- GIAB benchmark manifest wrote its region strings 0-based after region strings became 1-based
  (#63), so two negative regions started on a truth SNV and `make giab-benchmark` failed with 2
  false positives; `fetch_giab.py` now writes 1-based strings, and both published operating
  points reproduce.
- Viewer lost the reference genome on pan/zoom/detail refetches (mismatch evidence vanished,
  CRAM failed to decode after navigation) — the full data-source context is now preserved.
- Negative ClinVar/gnomAD results (variant-not-found) were re-fetched every lookup — a distinct
  cache-miss sentinel now caches absences.
- Region cache and the services registry were unbounded — both are bounded now; the registry
  also guards against `id()` reuse and removes entries on close.
- pysam parses ran on the event loop's default executor, so a hung parse could starve unrelated
  async work — parsing is confined to a dedicated bounded pool.

### Added

- **`/__version` deploy identity** — serves `{"sha": …}` (the commit the image was built from,
  baked in via the `GREENLIGHT_SHA` build-arg) at `/__version` and `/mcp/__version`, so the
  Greenlight deploy gate can prove the running container is the build CI just shipped rather than
  the still-running previous image. Unauthenticated; a build without the arg reports `null`.
- `schema_version` stamp on the viewer payload + a server↔viewer compatibility test.
- Eval grader scores **tool selection and answer correctness separately** (a case can no longer
  pass on the tool call alone); nightly runs both truth sets and, when `ANTHROPIC_API_KEY` is
  set, the Anthropic provider (the OpenAI path is gated until the runner translates
  provider-specific messages).
- Frontend static/test layer: `tsc --noEmit` (strict), ESLint, and vitest, wired into CI.
- Newcomer features (#67): an automatically selected reference (the registry FASTA matching the
  BAM's build and contig naming), `demo:<NAME>` aliases with the `bamcp://demos`, `glossary` and
  `guide` resources and three prompts, per-variant interpretation, a probe of the requested site
  (`filter_check`), a time-bounded `scan_variants` that returns `next_region` to resume, and a
  CIGAR-aware viewer.
- Eval results publish to Tracer from CI on `main` and nightly (best-effort; needs
  `TRACER_INGEST_TOKEN`).
- Evaluation: a real-data GIAB HG002 tier and pre-registered visual-vs-text rounds 3–5
  (#69, #70, #74); results in `docs/eval/visual-results.md`.
- `make verify` runs every `ci.yml` gate in one command (#63).

### Changed

- Region strings at the tool boundary are 1-based inclusive (samtools-style); every tool error
  goes through one sanitized `isError` envelope (#63).
- Viewer rendering: neutral read bodies, colour only on mismatches, strand shown by arrow and
  underline, hollow tiles for low-base-quality mismatches (#68).
- Docs lead remote deployment with **Streamable HTTP** (SSE is deprecated by the MCP spec).
- Test suite is hermetic by default (live-network tests behind a `network` marker); CI's
  `lockfile` job checks `uv.lock` is fresh and runs the ship-gate from the frozen lock, while the
  3.10–3.12 `test` matrix installs broad ranges as the compatibility check.

### Biological correctness (peer-review round)

A hands-on peer review caught a class the automated review rounds missed: silent output that
misleads the caller. The unifying fix is **self-describing outputs — every response states the
configuration that shaped it**. See `docs/ENGINEERING_RETROSPECTIVE.md`.

- **Genome-build correctness** — `search_gene` returned coordinates for a static server build with
  no tie to the BAM (a GRCh38 locus in a GRCh37 BAM is silently ~0.5–2 Mb off the real gene). It
  now auto-detects the build from the BAM (`file_path`), allows an explicit `build` override, states
  `genome_build`/`build_source`, and warns loudly on any mismatch. Also fixed the assembly selector
  that matched GRCh37's `GCF_000001405.25` as GRCh38.
- **Self-describing filters** — `get_variants`/`visualize_region`/`jump_to`/`get_region_summary`/
  `scan_variants` now echo the `applied_filters` (VAF/depth/MAPQ/baseQ, and `max_reads` where reads
  are shown), so a variant's *absence* is interpretable (filtered vs. truly absent), not silent.
- **Error vs. absence** — `lookup_clinvar`/`lookup_gnomad` carry an explicit `status`
  (`found`/`not_found`/`unavailable`/…); a transient gnomAD GraphQL error now propagates instead of
  collapsing to `None`. In the ACMG scaffold a failed lookup is marked `unavailable` rather than a
  null "absent" block, so a network blip can't become false rarity (PM2) evidence.
- **Ruler rendering** — position-ruler tick labels stay distinct at every zoom (were collapsing to
  `1.1K` repeated) and no longer collide with the contig label; extracted a tested
  `formatRulerLabel`.
