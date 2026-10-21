# Visual-vs-text eval, round 5 — pre-registration (confirming the spanning-read rule)

**Status:** registered 2026-10-21. The commit adding this file precedes the item draw; the
manifest hash is added in a second commit after the draw and before any model answer. Later
changes go in *Amendments* with a date and a reason.

## Why this round

The round-4 follow-up ([`visual-results.md`](visual-results.md#round-4-follow-up--read-based-homopolymer-risk-exploratory))
changed BAMCP's curation evidence for indels in homopolymer runs (BAMCP #71): the homopolymer
concern and the allele fraction now come from reads that span the whole run. Re-measured on the
60 round-4 items, BAMCP's text tools raised balanced accuracy over raw `samtools mpileup` from
0.839 to 0.911 (+0.072, 95% CI +0.011 to +0.131) and correct genotypes on true indels from
13/90 to 42/90. That result is exploratory: the items are the ones whose round-4 answers
motivated the change. This round tests it on items no earlier round has seen.

## Items

- **Pools.** Round 4's scanned repeat-run candidates on chr19–22 plus new scans of chr15–18,
  built by the same code (`select-r4 --pairs tr`, harness at this commit), main partition only.
- **Exclusion of earlier rounds.** `--exclude-manifest docs/eval/data/r4_tr_main_sites.jsonl`:
  every candidate within 10 kb of a round-4 site, or in one of their segmental-duplication
  clusters, is dropped before the draw. Pilot-partition clusters are excluded by the partition.
- **Strata and rules.** As registered for round 4 at the alternative difficulty setting
  ([`visual-preregistration-r4.md`](visual-preregistration-r4.md#items)): `tr_indel_pos`
  (v5.0q PASS indels of 1–4 bp in a homopolymer ≥ 12 bp or a di/tri-nucleotide repeat ≥ 20 bp,
  stutter ≥ 12% or het allele fraction ≤ 0.40, HiFi modal lengths at the truth genotype) and
  `tr_stutter_neg` (matched one-to-one on run-length bin; HiFi and trio checks of amendments
  A1, A3, A4). Spacing ≥ 10 kb, ≤ 2 items per cluster.
- **Size rule (fixed now).** Draw seed **701**, quota 45 per stratum. If the matched negatives
  fall short, the draw is re-run with the same seed and quota n for both strata, where n is the
  number of matched negatives the first draw found, provided n ≥ 30 (the strata stay equal and
  the het/hom cells keep their shares). Below 30 the round does not run and the pools are
  extended first (recorded as an amendment).
- The frozen manifest is committed as `docs/eval/data/r5_tr_confirm_sites.jsonl` with its
  sha256 in *Manifest* below before the first run.

## Conditions, models, runs

| code | files | role |
|---|---|---|
| `Tmp` | `mpileup.txt` | baseline: unfiltered `samtools mpileup -B -Q 0 -q 0 -s --output-BP`, site ± 30 bp |
| `Tmp+Tcur` | `mpileup.txt`, `tcur.md` | + BAMCP's text tools (`get_region_summary`, `get_variant_curation_summary` at the top non-reference allele) |

BAMCP product code is `main` at `24538b8` (BAMCP #71); nothing in `src/bamcp/` outside
`eval/` changes before the runs. Caption, prompt, answer schema, run harness and the three
model tiers (small, mid, large) are as in round 4, one run per item × condition × model, one
fresh blind subagent per run with Read access to its listed files only; labels sealed during
the runs; transcript audit as before. The viewer arm is not re-run.

## Hypotheses

Both are confirmatory, Holm-corrected as one family of two at α = 0.05, two-sided; the
registered direction is `Tmp+Tcur` better.

- **H1 (primary): balanced accuracy.** Strict joint correctness (`variant_present` and
  `artifact_verdict`; `uncertain` is wrong) as balanced accuracy over the two strata,
  equal-weighted over the three models. Paired cluster sign-flip permutation test (10,000
  flips), 95% CI by test inversion, TOST at ± 0.07. Supported = Holm p < 0.05, CI excludes 0,
  Δ > 0, and the guardrail holds (neither sensitivity nor specificity falls by more than 0.05).
- **H2: genotype on true indels.** Among `tr_indel_pos` items, a run is genotype-correct when
  it is jointly correct and its zygosity, variant type and allele match the truth
  (`grade.genotype_correct`). Exact McNemar on item × model pairs. Supported = Holm p < 0.05
  and `Tmp+Tcur` wins more discordant pairs.

Grading: `python -m bamcp.eval.visual grade-r4 --run-dir <set> --genotype-confirmatory`.

**Power** (resampling the follow-up's 44 per-cluster differences, two-sided 0.05): for H1,
about 0.82 with 45 + 45 items and 0.59 with 30 + 30 if the follow-up effect holds in full, and
about 0.44 and 0.26 at three quarters of it. H2 is well powered at either size (the follow-up
had 33 vs 4 discordant pairs of 90). A null H1 with a supported H2 would mean the change
improves genotypes without a detectable gain on the real-vs-artifact call.

Secondary, descriptive: sensitivity and specificity per arm; per model; insertions vs
deletions; het vs hom; utility; Brier and confidence on right vs wrong.

## Exclusions

As in round 4: invalid or missing answers are wrong. A run is re-run once with the identical
prompt only for an audit violation, or an infrastructure failure (including a usage limit)
that kept the model from producing an answer. No item is dropped after seeing results.

## Manifest

(added after the draw, before the first run)

## Amendments

(none)
