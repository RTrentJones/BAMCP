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
- **One URI-scheme classifier** — a `vcf_path` such as `http:/127.0.0.1:<port>/x.vcf.gz` (one
  slash, so no `://`) was treated as a local path and reached htslib, which fetched it over HTTP
  without the SSRF checks. Every data source now goes through `_is_remote_source`: only an exact
  `http://` / `https://` prefix is remote, and any other scheme-looking value (`data:`, `ftp:`,
  `file:`, `s3:`, `http:/`, …) or a local path embedding `://` / `##idx##` is rejected. A local
  VCF/BCF must exist.
- **`BAMCP_ALLOW_LOCAL_FILES`** (default `true`) — `false` rejects every caller-supplied local
  path before any filesystem access; the operator's `BAMCP_REFERENCE` keeps working. The compose
  `prod` profile sets it to `false`, sets `BAMCP_ALLOWED_REMOTE_HOSTS`, and binds the port to
  `127.0.0.1`.
- **Request body limit** — `BodySizeLimitMiddleware` answers 413 above `BAMCP_MAX_REQUEST_BYTES`
  (default 1 MiB; 64 KiB on `/register`, `/token`, `/revoke`, `/authorize`), on the declared
  `Content-Length` and on streamed bytes; a malformed or conflicting `Content-Length` is a 400. It
  sits after `SecurityHeadersMiddleware` and before the rate limiter.
- **Client-registration caps** — at most 10 `redirect_uris` of at most 2048 characters, a
  `client_name` of at most 256 characters, and 4 KiB of metadata per client (a real registration
  is well under 1 KiB).
- **IPv6 rate-limit buckets** — an IPv6 client is keyed by its /64 (IPv4-mapped addresses by the
  embedded IPv4), so rotating source addresses inside one /64 does not dodge the limit.
- **CRAM reference lookup stays local** — at startup `REF_PATH` / `REF_CACHE` default to a
  directory under the cache dir, so htslib never fetches a reference by M5 checksum from the
  network; an operator-set URL value is kept and logged with a warning.
- **CRAM header references gated** — htslib decodes a contig against the path in its `@SQ UR:`
  tag when no reference is given or the given reference lacks that contig, so a caller-supplied
  CRAM could return a server file's bytes as read sequence. With `BAMCP_ALLOW_LOCAL_FILES=false`
  a CRAM needs a resolved reference that covers every `UR:`-tagged contig in its header. The
  gate is applied on the opened file's detected format, not its name: a CRAM renamed `.bam`
  (htslib detects CRAM from the bytes) used to skip it and decode against its `UR:` file.
- **NCBI API key kept out of logs** — `str()` of an httpx error embeds the request URL, so a failed
  ClinVar/gnomAD/gene lookup logged `api_key=…`, as did httpx's own per-request INFO line. Every
  log site now goes through `redact_error` (type, status, URL with credential parameters
  replaced), and a filter scrubs the `httpx` logger.
- **`scan_variants` memory bounded while scanning** — every candidate of a scan (up to 250 Mbp,
  at caller floors as low as `min_vaf=0`, `min_alt_count=1`) was held until the end and only then
  capped at 500. The highest-VAF candidates are now kept as the scan goes (same selection, ties
  included); `total_candidates` / `results_capped` still count everything seen.
- **Viewer variant list capped** — `visualize_region` / `jump_to` put every passing candidate of
  a window (up to 1 Mbp) into the viewer payload. They now keep `GET_VARIANTS_MAX_RESULTS`,
  chosen as `get_variants` chooses, record `variants_truncated` / `total_candidates` /
  `max_results` in `applied_filters`, and the viewer says how many are shown.
- **`search_gene` symbol validation** — symbols must be 1–40 letters, digits, `.`, `_` or `-`
  before they are put into an NCBI Entrez query.
- **`BAMCP_HOST` defaults to `127.0.0.1`** in code; the Docker image still sets `0.0.0.0`.
- **Lookup answers no longer echo upstream error text** — an `unavailable` ClinVar/gnomAD answer
  carries a fixed `reason_code` (`timeout`, `network_error`, `upstream_error`,
  `unsupported_build`); the upstream detail is logged server-side.
- CI: `packages: write` is granted only to the image `build` job; CodeQL also scans
  `javascript-typescript` (the viewer).

### Fixed

- **Mismatch-capped reads no longer weaken SNV evidence silently** — a read past
  `READ_MAX_MISMATCHES` stops listing mismatches, so it stopped supporting later candidates
  while the pileup alt count stayed the same. Responses that report variants now set
  `applied_filters.mismatches_truncated` (+ `mismatches_truncated_reads`), and candidates such a
  read covers past its kept mismatches are marked `evidence_incomplete`.
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
- Region cache and the services registry were unbounded — both are bounded now; each registry
  entry holds its own config (so its `id()` cannot be reused while registered) and is removed on
  close.
- Per-call resource bounds: at most `VCF_MAX_RECORDS` (5,000) VCF records per call
  (`vcf_truncated` + `vcf_max_records` in `applied_filters`); records with more than
  `VCF_MAX_SAMPLES_DETAIL` (16) samples keep the GT of the first 16 (`samples_truncated`,
  `sample_count`; the viewer shows the full count and `get_region_summary` appends
  `… +N more samples`); at most `READ_MAX_MISMATCHES` (1,000)
  mismatches per read (`mismatches_truncated`); read extraction stops past an estimated
  `READ_PAYLOAD_MAX_BYTES` (64 MiB, `reads_truncated`); `get_variants` returns at most
  `GET_VARIANTS_MAX_RESULTS` (1,000) variants (`truncated`, `total_candidates`, `max_results`),
  and `get_region_summary` lists at most as many.
- ClinVar/gnomAD timeouts were reported as network errors on Python 3.11+ (`asyncio.TimeoutError`
  is an `OSError` there); the timeout branch is now reached, and every lookup outcome carries
  `genome_build` / `build_source`.
- Index URLs for signed BAM/CRAM/reference URLs put `.bai`/`.crai`/`.fai`/`.gzi` inside the query
  string; the suffix now goes on the URL path and the query is kept.
- `scan_variants` re-raised pysam's `ValueError` text unsanitized; an unknown contig is now a
  `ToolInputError` naming the file's own spelling, and other errors go through the sanitizer.
- `get_variant_curation_summary` returned invalid input, timeouts, unknown contigs and unreadable
  files as ordinary results; invalid input now raises `ToolInputError` (checked before reference
  resolution, as `classify_variant` now also does) and other failures take the `isError` envelope.
- 10 of 13 tool descriptions in `tool_specs.py` had drifted from the registered text; the eval
  router sent placeholder descriptions to the eval model.
- `cleanup_cache` now also forgets which remote inputs are warm, so the next open of a remote
  input gets the cold-start budget again.
- A candidate with an unknown `alt_count` was dropped by the viewer/region-summary floor and the
  site probe but kept by `get_variants`; all of them now keep it.
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
- `compact_threshold_bp` in the viewer's `ui/init` payload (`schema_version` unchanged); the
  viewer's sequence auto-fetch uses it, with a 500 bp fallback for older servers.

### Changed

- Region strings at the tool boundary are 1-based inclusive (samtools-style); every tool error
  goes through one sanitized `isError` envelope (#63).
- Viewer rendering: neutral read bodies, colour only on mismatches, strand shown by arrow and
  underline, hollow tiles for low-base-quality mismatches (#68).
- Docs lead remote deployment with **Streamable HTTP** (SSE is deprecated by the MCP spec).
- Test suite is hermetic by default (live-network tests behind a `network` marker); CI's
  `lockfile` job checks `uv.lock` is fresh and runs the ship-gate from the frozen lock, while the
  3.10–3.12 `test` matrix installs broad ranges as the compatibility check.
- `tool_specs.TOOL_SPECS` is the one source of tool descriptions: `server.py` registers from it,
  the eval router sends its descriptions and schemas, and the parity test checks descriptions as
  well as parameters.
- `core/tools.py` split into `core/services.py`, `remote.py`, `lookups.py`, `scan.py` and
  `filters.py` (one candidate-filter predicate); the region tools share one `RegionRequest`
  preamble; one genome-build alias table in `constants.py`.
- Removed dead code: `validate_data_sources`, `get_tool_spec`, telemetry `is_enabled` /
  `get_session_id`, the test-only `parsers.scan_contig` / `scan_variants_chunked`, and unused
  visual-eval constants.
- Docker: the runtime `HEALTHCHECK` GETs `/__version` on HTTP transports; the entrypoint's
  pre-flight runs the import check (`--preflight`).
- One Node version for the viewer build, in `src/bamcp/static/.nvmrc`, read by every workflow.

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
