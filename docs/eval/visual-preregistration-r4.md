# Visual-vs-text eval, round 4 — pre-registration (focused: repeat-run indels)

**Status:** registered 2026-10-17 before any main-set model answer; the commit adding this file
precedes the first main run. Later changes go in *Amendments* with a date and a reason.

## Why this round is narrower than the design

The full round-4 design ([`visual-design-r4.md`](visual-design-r4.md)) failed its baseline gate
twice ([`visual-pilot-r4.md`](visual-pilot-r4.md)): on real HG002 sites from the classes where
short-read callers fail most, current models scored 0.96 with an expert text summary and 0.94
with raw `samtools mpileup`. The one class with a baseline near 0.70 under both was indel
genotyping in long repeat runs at the alternative difficulty setting (0.67 with the summary,
0.58 with mpileup; stutter negatives at that setting 0.83 with mpileup). This round tests
whether BAMCP helps there, against the raw-input baseline a user without BAMCP would have.

## Items

60 items from the **main** partition (no cluster shared with any pilot item), drawn by
`select-r4 --split main --seed 501 --pairs tr --knob tr_indel_pos=alt --knob tr_stutter_neg=alt
--quota tr_indel_pos=30 --quota tr_stutter_neg=30` over the scanned candidates of chr19–22, at
harness commit `1b20d44`. Frozen manifest [`data/r4_tr_main_sites.jsonl`](data/r4_tr_main_sites.jsonl),
sha256 **`aeb16948e02f8925908064a27295c530c9f2b0ea8184f7465eded67b32140185`**.

- **`tr_indel_pos` (30):** v5.0q PASS indels of 1–4 bp, whole units of a homopolymer ≥ 12 bp or a
  di/tri-nucleotide repeat ≥ 20 bp, with ±1-unit stutter ≥ 12% or allele fraction ≤ 0.40 (het)
  / ≤ 0.85 (hom) on spanning reads; HiFi modal lengths at the truth genotype. 18 insertions and
  12 deletions; 20 het, 10 hom. Label `yes` / `real`.
- **`tr_stutter_neg` (30):** runs of the same kinds inside the v5.0q BED with no truth record in
  the run ± 30 bp, a ±1-unit allele on ≥ 3 spanning reads and ≥ 5%, reference as HiFi's modal
  length with the candidate ≤ 25% of HiFi reads, and the allele recurring in the parents
  (amendments A1, A3, A4 of the design, made before any pilot run). Matched one-to-one to the
  positives on run-length bin. Label `no` / `artifact`, `artifact_type: homopolymer_indel`.

Known composition difference: all 30 negatives are one-unit deletions (stutter is mostly
contraction), while 18 of 30 positives are insertions, so the position of interest is the run's
first base for every negative but the base before the run for the insertions (POI-offset KS
lint p < 0.001 after 5 redraws). This can move absolute accuracy but not the paired contrasts,
which compare arms on the same items; accuracy is also reported separately for insertions and
deletions.

## Conditions, models, runs

| code | files | role |
|---|---|---|
| `Tmp` | `mpileup.txt` | baseline: unfiltered `samtools mpileup -B -Q 0 -q 0 -s --output-BP`, site ± 30 bp |
| `Tmp+Tcur` | `mpileup.txt`, `tcur.md` | + BAMCP's text tools as shipped (`get_region_summary`, `get_variant_curation_summary` at the top non-reference allele) |
| `V+Tmp` | `mpileup.txt`, `view_wide.png`, `view_zoom.png` | + the BAMCP viewer as shipped (viewer source hash `6ccf323b…`, unchanged since round 3) |

BAMCP product code is as on `main` (`4ed4cab`); this branch changes only `src/bamcp/eval/`.
None of the design's product changes (J1–J6, V1–V6) were made. Caption, prompt, answer schema
and harness as in the pilots (caption sample line and repeat-run anchoring sentence of the r4
design). Models: the three tiers (small, mid, large), one run each.
60 × 3 × 3 = **540 runs**, one fresh blind subagent per run; labels sealed during the runs;
transcript audit as before.

## Hypotheses

Primary endpoint: strict joint correctness (`variant_present` and `artifact_verdict`;
`uncertain` is wrong), aggregated as **balanced accuracy** (mean of the positive- and
negative-stratum accuracies), equal-weighted over the three models
(`python -m bamcp.eval.visual grade-r4`).

Confirmatory, Holm-corrected as one family at α = 0.05, two-sided:

- **C1: `Tmp+Tcur` vs `Tmp`.** Predicted direction: **no gain or harm**. As shipped,
  `get_variant_curation_summary` raises a HIGH "homopolymer run" concern for every indel in a
  run ≥ 4 bp, true or not, which should push answers toward "artifact": specificity up,
  sensitivity down.
- **C2: `V+Tmp` vs `Tmp`.** No direction predicted.

Test: paired sign-flip permutation test on item-level differences (model-averaged, weighted by
1 / (2 · n_class)), flips by cluster, 10,000 flips; 95% CI by test inversion; supported = Holm p
< 0.05 and CI excludes 0. Equivalence: TOST at ± 0.07; a null is read as "no meaningful gain"
only if TOST succeeds. Expected baseline from the pilots: about 0.71. With 60 items the SE of
a paired difference is about 0.04, so this round can detect gains or losses of about 0.10 and
above, and is inconclusive for smaller ones.

Secondary, descriptive: sensitivity and specificity per arm; per model; insertions vs
deletions; het vs hom; genotype (zygosity) accuracy among positives; utility (+1 / 0 / −1);
Brier and confidence on right vs wrong.

## Exclusions

As in round 3: invalid or missing answers are wrong. A run is re-run once with the identical
prompt only for an audit violation, or an infrastructure failure (including a usage limit)
that kept the model from producing an answer. No item is dropped after seeing results.

## Amendments

- **2026-10-18, after the runs (grading code; no effect on the data).** `grade-r4` had its
  contrast family hard-wired to the design's `Traw-r4` arms, so it reported no contrasts for
  this run set; it now uses the family of the run set's baseline arm (`Tmp`). The test,
  weights and Holm family are as registered above. 240 runs first stopped on a usage limit
  before answering and were re-run with the identical prompt (Exclusions); one run produced two
  identical answers and is counted once.
