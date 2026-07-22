# Engineering Retrospective

BAMCP explores how an AI assistant, a typed scientific tool server, and an interactive evidence
viewer can operate as one verifiable system. This note is the honest account the code alone can't
give: what I designed, what agents produced, what verification caught, and which trade-offs remain.
It was prompted by an external engineering review; the response plan is in
[`engineering-review-response.md`](engineering-review-response.md).

## What I set out to build

Not "a variant caller behind MCP." The thesis is that a domain agent tool is one system —
**model + structured genomic evidence + a human verification surface** — and that all three should
be *evaluated together*. The differentiators follow from that: an interactive alignment viewer that
returns as an MCP App (not text/JSON), progressive retrieval instead of dumping large payloads, a
clean separation between evidence display and clinical interpretation, and deterministic biological
truth sets rather than leaning only on a model judge.

## What I designed vs. what agents produced

Much of the implementation was written with heavy agent assistance; that is visible in the commit
history and I don't hide it. What I supplied was the judgment: the architecture (thin FastMCP
handlers over a `core`/`analysis`/`clients`/`middleware` split), the evaluation design (planted
ground truth, negative controls, a positive control so safety checks can't pass vacuously, an
artifact-region recall channel), the safety invariants (in-band disclaimers verified by tests), and
— importantly — deciding **which claims the evidence does not support**.

The most useful thing I did was not accept agent output at face value. Live verification against
real external schemas (ClinVar, SPDI, VCF quality fields, gnomAD GraphQL) repeatedly contradicted
plausible-looking code, and each contradiction drove a real change.

## What verification caught

A few examples where a check falsified something that looked correct:

- **An unvalidated SSRF path.** The `reference` argument flowed straight to
  `pysam.AlignmentFile(reference_filename=…)` with no validation, even though the BAM, VCF, and
  index paths were all validated. Writing the unified data-source policy surfaced the one input
  that had slipped through — a remote reference could have reached the cloud metadata endpoint.

- **A memory-leak "fix" that couldn't work.** My first attempt to bound the services registry used
  a `weakref` finalizer to drop an entry when its config was garbage-collected. A test proved it
  never fired: `BAMCPServices.config` holds a *strong* reference to the config, so the config is
  pinned alive as long as the entry exists. The finalizer was theater; the real fix is a bounded
  LRU plus explicit removal on close. The test is the reason I know that.

- **Dead code the type checker had never seen.** The viewer TypeScript had only ever been
  transpiled, never type-checked. Adding `tsc --noEmit` immediately flagged an unused private
  method and an unused destructured field in `renderer.ts` — small, but exactly the kind of rot a
  build-that-doesn't-typecheck lets accumulate.

- **A grader that measured the wrong thing.** The eval grader passed a case whenever the expected
  tool was invoked, regardless of whether the answer was correct. It scored *tool selection* and
  called it *task success*. The fix separates the two dimensions and requires the answer to carry
  the required factual claims.

## Trade-offs I made deliberately

- **Service-token auth over an identity provider.** For a single-instance server behind a
  Cloudflare Tunnel, integrating a real IdP would be more code defending a threat model that
  doesn't exist. I removed the dynamic-registration theater and adopted a stateless service token
  instead. Interactive OAuth remains available behind a flag for dev.

- **Allow-list over redirect-revalidation for SSRF.** Because pysam/htslib performs the actual
  fetch (and follows redirects invisibly to the application), "validate every redirect hop" isn't
  achievable without re-implementing htslib's HTTP layer. The host allow-list checked before
  hand-off is the honest control; the residual (intra-allow-list DNS rebinding) is documented
  rather than papered over.

- **Backpressure over true cancellation for parsing.** `asyncio.wait_for` abandons a timed-out
  parse but can't cancel the underlying C call. Confining parses to a dedicated bounded pool stops
  them starving the event loop; a genuinely killable parse would need a worker *process*, which is
  a deliberate follow-up, not a claim I'm making today.

## What the evidence does *not* support (yet)

- The biological benchmark is one 60 kb slice of simulated reads over real GIAB genotypes. It shows
  the detector recovers planted SNVs on that slice — **not** general caller correctness across
  genomes, indels, SVs, or real sequencer error.
- The LLM eval now grades answer correctness, but a broad, human-labeled benchmark with judge-swap
  concordance is still ahead. Running model evals through the real MCP transport (rather than the
  in-process router) is scoped but not done.
- The payload contract is versioned and compat-tested, but not yet a generated, fully-typed
  schema. That migration is planned, not finished.

## The honest summary

The strongest signal here isn't the volume of implementation — it's the evaluation and
verification discipline: deterministic truth sets, negative and positive controls, adversarial auth
tests, and a habit of writing the check that proves the agent-generated code wrong. The remaining
gaps are named above rather than hidden, which is the point.
