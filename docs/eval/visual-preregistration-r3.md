# Visual-vs-text eval, round 3 — pre-registration

**Status:** registered before any round-3 model answer was generated (2026-10-16); the
commit that adds this file precedes the first run. Changes
after the first model run go in *Amendments* with a date and a reason, and are reported
alongside the results.

Round 1 (pilot, pre-registered in [`visual-preregistration.md`](visual-preregistration.md))
and round 2 (exploratory re-run of the image conditions after the viewer legibility fix) are
reported in [`visual-results.md`](visual-results.md). Round 2 was measured on the same items
whose failures motivated the fix, so its gain (`V` 0.821 → 0.910) is a lead, not a result.
This round tests it on items nobody has seen, and adds real sequencing data, where the text
conditions are not expected to sit at the ceiling.

Everything not stated here is unchanged from round 1: answer schema, prompt and file
descriptions, `VIEWER_KEY`, the primary endpoint, scoring and normalization, invalid-answer
rule, the subagent harness (one fresh blind subagent per run, only `Read` of its listed files
and the structured answer allowed), label vaulting, and transcript audit.

Viewer: the renderer as of BAMCP `ae5111b` (neutral read bodies, colour only on mismatches,
hollow low-quality mismatches). Models: the same three tiers as before (haiku, sonnet, opus),
one repeat.

## Tiers

### S — fresh synthetic items

`python -m bamcp.eval.visual build --seed 2 --instances 2` (the round-1 generator, a seed never
used before): **36 items**, all 18 case types × 2 instances. Conditions **`Traw`, `Tjson`,
`V`**. 36 × 3 × 3 = **324 runs**.

### R — real HG002 reads (GIAB v4.2.1 truth)

Reads: the public GIAB HG002 Illumina 60x GRCh38 BAM (the `demo:HG002` dataset). Truth: GIAB
HG002 v4.2.1 benchmark VCF and the `noinconsistent` confident BED. Context: the GIAB v3.3
`AllTandemRepeatsandHomopolymers_slop5` stratification. Code: `src/bamcp/eval/visual/real.py`.

Site selection (`python -m bamcp.eval.visual select-real --seed 3`, chr20 only): seeded random
2 kb windows inside confident intervals ≥ 2 kb are scanned until each stratum's candidate pool
holds 3x its quota (or 600 windows), then each quota is drawn with the same seed, every chosen
site ≥ 1 kb from every other.

| stratum | n | definition |
|---|---|---|
| `real_snv_het` | 6 | truth SNV 0/1, outside repeats |
| `real_snv_hom` | 6 | truth SNV 1/1, outside repeats |
| `real_snv_in_repeat` | 4 | truth SNV (0/1 or 1/1) inside the repeat/homopolymer stratification |
| `real_indel_het` | 6 | truth indel 1–10 bp, 0/1, any context |
| `real_indel_hom` | 4 | truth indel 1–10 bp, 1/1, any context |
| `real_noisy_negative` | 8 | no truth variant within 10 bp; ≥ 3 reads and ≥ 5% carry the same non-reference allele; half inside repeats |
| `real_clean_negative` | 2 | no truth variant within 10 bp; no non-reference read; depth ≥ 20; outside repeats |

Indels are taken in any context because isolated indels outside repeats are too rare on chr20
to fill a stratum (a first 100 kb scan found none); most selected indels will sit in
homopolymers or tandem repeats, as real indels do. This was decided before any model run.

Truth variants are used only when biallelic, fully inside a confident interval, and with no
other truth variant within 30 bp, so the position of interest has one answer. Labels follow
round 1's conventions (deletion: first deleted base, `alt_allele` = deleted bases; insertion:
anchor base, `alt_allele` = inserted bases). Negatives are `variant_present: no`,
`artifact_verdict: artifact`; a noisy negative's `artifact_type` is `unknown` and **not
scored** (the truth set says the sample is reference, not why reads disagree).

Blinding: each site's reference (3 kb) and reads are copied into a local BAM on re-indexed
contigs `chr1..chr3` (one 3 kb slot per item, seeded slot order), so no coordinate can be
looked up. Kept reads are primary, non-duplicate, non-QC-fail alignments overlapping the site
± 450 bp; they are written single-end with hashed names and no tags, so every condition counts
the same reads. The caption says "sequenced human sample (Illumina short reads); coordinates
are re-indexed" and the prompt says "a sequenced sample" in place of "a simulated sample".

Conditions **`Traw`, `Tjson`, `V`, `V+Traw`**. 36 × 4 × 3 = **432 runs**.

Frozen site manifest: [`data/r3_real_sites.jsonl`](data/r3_real_sites.jsonl) (36 sites, pools filled after 111 windows), sha256 **`a3e85c3dcbe91511c33577212d6bac4897fc638a2b1dc2e2ce961ab645f3232b`** (recorded in each
run set's `build.json`).

## Run sets as built

Both run sets were built and leak-linted (0 problems) before this document was committed; the
viewer source hash (`6ccf323b…`) is recorded in each `build.json`. Labels are sealed in an archive outside the
run directories for the whole run and restored only for grading.

## Hypotheses

Primary endpoint as in round 1: joint correctness of `variant_present` and `artifact_verdict`
(`uncertain` is never correct).

Confirmatory contrasts, Holm-corrected as one family at α = 0.05 (two-sided):

- **C1 (tier S): `V` vs `Tjson`.** Predicted direction `V` > `Tjson` (round 2: 0.910 vs 0.859).
- **C2 (tier R): `V+Traw` vs `Traw`.** No direction predicted.
- **C3 (tier R): `V` vs `Tjson`.** No direction predicted.

Secondary, descriptive: tier S `V` vs `Traw`; tier R `V` vs `Traw` and `Tjson` vs `Traw`;
every contrast per model; per-stratum accuracy; among real variants, `zygosity` and
`alt_allele` accuracy; Brier score and mean confidence on right vs wrong answers.

## Analysis

Unchanged from round 1 (`python -m bamcp.eval.visual grade`): unit of resampling = item;
delta = mean over items of (acc_A − acc_B), each item's accuracy averaged over the three
models; 95% percentile CI from 2000 paired cluster-bootstrap resamples (seed 0); exact
two-sided McNemar on runs matched by (item, model). Holm is applied to the
three McNemar p-values; a contrast is "supported" only if its Holm-adjusted p < 0.05 and its
CI excludes 0. Tiers are analysed separately and never pooled.

## Exclusions

As in round 1: invalid or missing answers are wrong, never excluded. A run is excluded only
for an audit violation (a read outside its listed files) or an infrastructure failure that
kept the model from seeing the stimulus; such runs are re-run once with the identical prompt
and the count is reported. Runs that fail on a usage limit before producing an answer are
re-run with the identical prompt. No item is dropped after seeing results.

## Amendments

- **2026-10-16, after the runs (code review of the selection code; no effect on the data).**
  Two defects in `real.py` were fixed for future seeds: (1) truth proximity was checked against
  VCF anchor positions only, so a long truth deletion anchored more than 10 bp upstream could
  cover a "negative"; it now uses each record's full reference span. (2) The noisy-negative
  draw could silently fill a short in-repeat pool from outside repeats; it now fails instead.
  The frozen manifest was checked against the truth VCF: none of its 10 negatives lies within
  10 bp of any truth span, and its noisy negatives are 4 inside / 4 outside repeats, so the
  registered item set and every reported result are unchanged.
- **2026-10-16, after the runs (second review round; no effect on the data).** (3) An SNV's
  repeat flag was computed over its own base plus the next one, so an SNV directly before a
  repeat could be counted as inside it; it now uses the variant's own span (insertions also
  the base after the anchor). All four `real_snv_in_repeat` sites overlap a repeat on their
  own base and no other selected SNV is adjacent to one. (4) Scanning could stop when the
  combined noisy-negative pool was full but one repeat sub-pool was short; the stop rule now
  checks both. (5) The self-test oracle answered the unscored `unknown` type, which is not a
  valid answer value; it now answers `none`. The oracle self-test passes on the real run set.
- **2026-10-16, after the runs (later review rounds; no effect on the data).** Further
  selection-code fixes for future seeds: results.csv leaves an unscored `artifact_type` cell
  empty; a zero quota takes nothing; the scan stops only when each pool holds enough mutually
  spaced candidates; truth isolation excludes the current record by index, not by start
  position; only scanned windows count toward `max_windows`; `build-real` takes
  `--png/--no-png`; and a source read is copied into one item only. Checks on the frozen set:
  no selected variant shares a POS with another truth record; its closest pair of sites is
  1,502 bp apart, so no read can reach two items. **Reproduction:** re-running
  `select-real --seed 3` with the registered code (commit 4cb4048) reproduces
  `data/r3_real_sites.jsonl` byte for byte (sha256 a3e85c3d…, pools full after 111 windows).
  The current code selects a different set for seed 3 (pools fill later under the
  spacing-aware stop rule), as expected after these fixes; the registered set stays the one
  defined by the registered code.
