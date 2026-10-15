# Visual-vs-text eval — pre-registration

**Status:** registered before any model answer was generated (2026-10-15).
Any change to hypotheses, endpoints, scoring or exclusions after the first model run must be
recorded in the *Amendments* section with a date and a reason, and reported alongside results.

Code: `src/bamcp/eval/visual/` (`python -m bamcp.eval.visual`). Hermetic gate: `make eval-visual-dry`.

## Question

Does showing a model the BAMCP pileup **render** improve its analysis of a candidate site over
giving it the same evidence as **text**? We report the result honestly either way: "no
improvement" and "worse" are publishable outcomes.

## Hypotheses

- **H1 (additive value of the picture).** `V+Traw` > `Traw` on the primary endpoint.
- **H2 (picture vs. information-equivalent text).** `V` vs `Tjson` — two-sided; no direction is
  predicted. `Tjson` carries the same information as the render (it *is* the render payload).
- **H3 (secondary).** `Tfull+V` > `Tfull`; `V` vs `Traw`; `Tjson` vs `Traw` (two-sided).
- **H4 (secondary, exploratory).** Visual benefit concentrates in case types whose evidence is
  spatial (strand bias, read-end bias, low-MAPQ paralog, PCR-duplicate stack, phasing).

## Stimulus set

Deterministic generator `create_visual_eval_fixtures(seed)` (`tests/create_visual_eval_fixtures.py`
→ `tests/fixtures/vt_eval.bam` + `vt_eval_ref.fa`, generated, not committed). Default seed `1`.

- Non-periodic seeded reference, 3 contigs (`chr1..chr3`, "GRCh38-like" names), GC ≈ 41%,
  natural runs capped at 3 so the only homopolymers are planted.
- Single-end 100–150 bp reads, ~30x (40x for low-VAF, 3x for low-depth), jittered starts,
  both orientations, base qualities decaying along the read with Gaussian noise and occasional
  dips, ~0.3% random substitution errors (lower BQ on errors), hashed 16-hex read names.
- **18 case types × 3 seeded instances = 54 items**, each a planted site at a seeded random
  slot (slot order shuffled, so position does not encode case type):

| case type | truth | artifact_type | planting |
|---|---|---|---|
| het_snv | real SNV, het | none | hap-1 SNV, VAF 0.42–0.58 |
| hom_alt_snv | real SNV, hom_alt | none | both haplotypes |
| low_vaf_snv | real SNV, low_vaf | none | VAF 0.10–0.14 at 40x (≥3 alt reads) |
| low_baseq_alt | no variant | low_baseq | 25–35% of reads alt at BQ 4–12 |
| strand_bias | no variant | strand_bias | 35–50% of forward reads alt, reverse none |
| read_end_bias | no variant | read_end | 8–10 extra reads, alt ≤4 bp from the 3′ aligned end, ⅔ with 3′ soft clips |
| low_mapq_paralog | no variant | low_mapq | 10–12 extra MAPQ 0–3 reads with alt + 3 linked paralog mismatches |
| homopolymer_del | no variant | homopolymer_indel | poly-A 9, 12–18% of spanning reads with a left-aligned 1 bp deletion |
| real_del | real del, het | none | 1 bp deletion, non-repeat context |
| real_ins | real ins, het | none | 2–3 bp non-homopolymer insertion |
| multi_allelic | real SNV (two alts), het | none | alt1 on hap 1, alt2 on hap 0 |
| pcr_duplicate | no variant | pcr_duplicate | 7–9 identical (start, end, orientation, sequence) alt reads, unflagged |
| low_depth_het | real SNV, het | none | exactly 3 reads at the site (2 alt) |
| clean_negative | no variant | none | nothing planted |
| scattered_errors | no variant | none | 3% error rate in ±30 bp |
| mnv | real MNV, het | none | two adjacent substitutions on hap 1 |
| phased_cis | two real SNVs, het, cis | none | both on hap 1, 25–45 bp apart |
| phased_trans | two real SNVs, het, trans | none | one per haplotype |

Not included (generator complexity): soft-clip breakpoint pileup, het large-deletion depth dip,
misalignment cluster; also no paired-end reads. `misalignment` and `low_depth` stay in the answer
enum but are never a true label.

`tests/eval/visual/test_planted_signals.py` checks every planted pattern from the BAM itself.

## Conditions

Each (item, condition) is an opaque directory `stim/<sha8>/` containing `caption.md` (uniform
template: region, position(s) of interest 1-based, "simulated diploid sample on a synthetic
GRCh38-like reference"), `prompt.md` (fixed task + schema; identical for identical file sets)
and the condition's evidence:

| condition | evidence files |
|---|---|
| `Traw` | `traw.md` — pileup numbers at the site: ±15 bp reference context and per-position allele counts, per-allele counts by orientation, distinct alignments, MAPQ / BQ / 3′-distance histograms, 3′ soft-clip counts, and for two-site items the read-level allele combinations. **No** verdicts, risk, confidence, recommendations or flags. |
| `Tjson` | `tjson.json` — the viewer payload for the ~200 bp window passed through `renderer.blind_payload()` (drawable keys only, `variants: []`, hashed read names). Information-equivalent to `V`. |
| `Tfull` | `tfull.md` — shipped `get_region_summary` text for the window + `get_variant_curation_summary` for the strongest candidate within 2 bp of each site (or a note that none was reported). |
| `V` | `view_wide.png` (~200 bp) + `view_zoom.png` (~60 bp), expanded mode, rendered with `ViewerRenderer.capture(..., blind=True)`: pinned display settings, clipped to `#viewer`, `[data-blind-hide]` overlays hidden, viewport height auto-fitted to depth. `prompt.md` adds a viewer key (the legend is outside the clip). |
| `V+Traw` | `V` + `Traw` |
| `Tfull+V` | `Tfull` + `V` |

`Tfull` is the only condition allowed to contain tool verdict wording — that wording is what it
measures. Note for interpretation: the render is the viewer *as built at stimulus time*; the
viewer key in `stimuli.VIEWER_KEY` must be kept in sync with it.

Labels live in `labels/labels.jsonl`, outside `stim/`; the run manifest `runs.jsonl` maps
`run_id → (item, condition, model, repeat, files, answer_path)` and is never shown to a model.

**Known structural cue:** two-site (phasing) items list two positions of interest in every
condition; this reveals "a phasing item" equally across conditions, so contrasts are
unaffected, and we additionally report the primary endpoint excluding phasing items.

## Blinding

- Opaque stimulus directory names and item ids; labels stored elsewhere.
- **Leak-lint** (`python -m bamcp.eval.visual lint`, part of every build): fails on label
  vocabulary in evidence files (all but `prompt.md`/`tfull.md`), BAM read names, item ids or
  case-type names in any stimulus (PNG text chunks included), non-opaque or nested dirs,
  unexpected files, a caption off the template, or a `prompt.md` differing from the uniform text.
- Model-under-test runs are Claude subagents told to Read only the files listed for their run
  and to write only the answer JSON; transcripts are audited for any other file access.

## Answer schema

One JSON object: `variant_present` yes|no · `variant_type` SNV|ins|del|MNV|none · `alt_allele`
string · `zygosity` het|hom_alt|low_vaf|none · `artifact_verdict` real|artifact|uncertain ·
`artifact_type` none|strand_bias|read_end|low_mapq|low_baseq|homopolymer_indel|pcr_duplicate|misalignment|low_depth ·
`phase` cis|trans|na · `confidence` 0–1 (probability the variant_present + artifact_verdict call
is right) · `rationale` ≤ 40 words.

Label conventions: `variant_present` = the simulated sample truly carries a variant;
`artifact_verdict` = `real` when it does, `artifact` when it does not (non-reference signal is
noise or artifact); `alt_allele` = alt base(s) / deleted ref base(s) / inserted bases,
comma-joined and sorted for several alleles or sites, empty when absent.

## Endpoints and metrics

- **Primary endpoint:** joint correctness of (`variant_present`, `artifact_verdict`) per run.
  `uncertain` is never correct on the primary endpoint.
- Secondary: per-field accuracy (`phase` only on phasing items; `alt_allele` normalized:
  uppercase, `+`/`-` stripped, sorted comma list, none-like → empty); "lenient joint" accepting
  `uncertain` for `low_depth_het`; artifact precision/recall (positive = `artifact`);
  macro-F1 over true `artifact_type` classes + confusion matrix; Brier score and 10-bin ECE of
  `confidence` against the primary endpoint (valid answers only); per-case-type accuracy.

## Analysis plan

- Unit of resampling = **item** (cluster); repeats are averaged within item. Contrasts are
  paired within model: delta = mean over items of (acc_A − acc_B), **95% percentile CI from
  2000 paired cluster bootstrap resamples** (seed 0). Exact two-sided **McNemar** on runs
  matched by (item, model, repeat).
- Primary contrasts: `V+Traw` vs `Traw`, `V` vs `Tjson`, reported per model and pooled
  across models. Decision rule: H1 supported if the pooled CI lies above 0; H2 reported as
  the estimate + CI (two-sided). No multiplicity correction on the two primary contrasts
  (reported as such); secondary contrasts and per-case-type deltas are descriptive.
- Charts: accuracy dot plot (condition × model, 95% CI), forest plot of paired deltas by case
  type, reliability diagrams, artifact-type confusion matrices, stimulus gallery.

## Pilot

26 items (all 18 case types' instance 0, then 8 instance-1 items, stratified by
`--n-items 26`) × {Traw, Tjson, V, V+Traw} × 3 models (haiku, sonnet, opus) × 1 repeat.
The pilot estimates variance and validates the pipeline; it is not powered for H1/H2.

## Exclusion rules (decided now)

- An **invalid** answer (not JSON after stripping one optional code fence, missing field,
  enum violation, confidence outside [0, 1], rationale > 60 words) or a **missing** answer file
  is scored **wrong on every field** — never excluded.
- A run is **excluded** only when the transcript audit shows the subagent read a file outside
  its listed paths (or the labels), or an infrastructure failure meant the model never saw
  the stimulus (e.g. image failed to load). Excluded run ids are passed via `--exclude` and
  the count is reported.
- No item is dropped after seeing results.

## Self-tests (hermetic, CI)

`make eval-visual-dry`: fixtures → all conditions (PNG step skipped with a message when no
Chromium) → leak-lint → oracle answers score 100% on every field, uniformly random answers
score ≈ chance (1/6 on the primary endpoint, within max(0.08, 4 SE)), one corrupted answer is
counted invalid and wrong.

## Amendments

_None._
