# Visual-vs-text eval, round 4 — pilot (text baseline only)

**Status:** pilot batch 1 of the round-4 design ([`visual-design-r4.md`](visual-design-r4.md)), run
2026-10-17. Only the baseline condition `B` (`Traw-r4`) was run, as the design requires: no
treatment arm, so no product change has been tuned against round-4 model output.

## What was run

- Harness commit `16f7540` (frozen before the first run). Site manifest
  [`data/r4_pilot_b1_sites.jsonl`](data/r4_pilot_b1_sites.jsonl) (sha256 `064b067ab270269a…`): 52 pilot-partition
  items from four of the five pairs (mosaic 20, tandem-repeat 20, segdup 8, anchors 4). The
  indel-adjacent pair and the rest of the segdup pair were still being scanned.
- 52 items × 3 models (small / mid / large) × 1 = **156 blind runs**, one fresh subagent per run
  with only `Read` of its three listed files. Audit: 156/156 answered, 0 violations, 0 runs that
  skipped a file. 155 runs first failed on a usage limit before producing any answer and were
  re-run with the identical prompt (design section 6).

## Result

| | primary (balanced) | sensitivity | specificity |
|---|---|---|---|
| all models | **0.961** | 0.922 | 1.000 |
| small | 0.933 | 0.867 | 1.000 |
| mid | 0.983 | 0.967 | 1.000 |
| large | 0.967 | 0.933 | 1.000 |

Per stratum (joint accuracy over the three models):

| stratum | knob default | knob alt |
|---|---|---|
| `mosaic_snv_pos` (4–12 clean alt reads of 50–80) | 0.94 (18 runs) | 0.92 (12) |
| `mosaic_caller_fp_neg` | 1.00 (18) | 1.00 (12) |
| `tr_indel_pos` | 0.94 (18) | **0.67** (12) |
| `tr_stutter_neg` | 1.00 (18) | 1.00 (12) |
| `segdup_cmrg_pos` | 1.00 (9) | 1.00 (3) |
| `segdup_psv_neg` | 1.00 (9) | 1.00 (3) |
| anchors | 1.00 (12) | — |

Mean confidence was 0.78 on right answers and 0.68 on wrong ones.

## Go/no-go

**G1 fails: NO-GO.** The predicted main primary `B` is about 0.96, against the 0.62–0.78 window.
Every alternative knob level was piloted and none moves a stratum below 0.92 except
`tr_indel_pos` at the alternative setting, so the one permitted knob move cannot reach the window.
The indel-adjacent pair cannot change this: its pool's negatives carry the same low-base-quality
signature (mean alt base quality about 12) that made every other negative stratum easy.

The design predicted 0.69 (prediction interval ±0.14). It was wrong because the strengthened
baseline already contains what decides these sites: the per-allele strand, base-quality, read-end,
mate and linkage panel. Given it, all three models separate true low-fraction variants from
caller false positives, and stutter from real repeat indels, almost perfectly. The one
residual difficulty is genotyping repeat indels at low allele fraction in long homopolymers
(`tr_indel_pos` alt, 0.67; 4 of 5 errors called a true insertion stutter).

## Amendments made before the pilot (selection only)

- **A1.** Stutter-negative signal floor 3 reads / 5% (design: 4 / 10%). Benchmark-confident
  homopolymers carry 1–3% stutter; 4 / 10% left 2 of 299 scanned runs.
- **A2 → A3.** HiFi check for stutter negatives: HiFi's modal length is reference and the
  candidate length is on ≤ 25% of HiFi spanning reads. HiFi's own ±1 error in these runs is
  about 20%, so the SNV negatives' cap and the design's "≤ half of Illumina's" rule (A2) rejected
  every clean run.
- **A4.** Trio recurrence also pools the applicable parents (≥ 2 reads and ≥ 2% together); a
  3–5% artifact shows on 1–2 reads of a 40–60x parent.
- Scan-driven strata were extended from chr20 to chr19–22.
- A selection review fixed five defects before the pilot (HiFi read-and-fraction cap,
  isolation against every truth set, `anchor_clean_neg` artifact type unscored, insertion spans
  in the variable-copy exclusion, anchoring reach); see the commit log.

## What this means

With an expert-style text summary of the site, real hard-class sites in HG002 do not produce a
0.70 baseline for current models; the round-3 ceiling was not an artifact of easy items. The
headroom seen in rounds 1–3 lives in the *input format* instead: raw read payloads (`Tjson`,
0.82–0.93) and the picture alone (`V`, 0.87–0.89) scored below the summary panel (0.98–0.99).

## Pilot 2: raw-input baseline (`samtools mpileup`)

After the NO-GO the baseline was changed to what a user without BAMCP would paste: unfiltered
`samtools mpileup -B -Q 0 -q 0 -s --output-BP` of the position of interest ± 30 bp, with a fixed
column key (condition `Tmp`, harness `179b604`). Same 52 items, same three models, 156 blind
runs; audit 156/156 answered, 0 violations.

| | primary (balanced) | sensitivity | specificity |
|---|---|---|---|
| all models | **0.944** | 0.911 | 0.978 |
| small | 0.900 | 0.833 | 0.967 |
| mid | 0.967 | 0.967 | 0.967 |
| large | 0.967 | 0.933 | 1.000 |

| stratum | knob default | knob alt |
|---|---|---|
| `mosaic_snv_pos` | 1.00 | 0.92 |
| `mosaic_caller_fp_neg` | 1.00 | 1.00 |
| `tr_indel_pos` | 0.89 | **0.58** |
| `tr_stutter_neg` | 1.00 | 0.83 |
| `segdup_cmrg_pos` / `segdup_psv_neg` | 1.00 / 1.00 | 1.00 / 1.00 |
| anchors | 0.83–1.00 | — |

**Still NO-GO (0.94).** Removing the expert summary cost 0.017 overall: the models decode raw
mpileup strings (strand case, Phred characters, indel tokens, read positions) almost as well as
the pre-digested panel. The one consistently hard class, under both baselines, is
**genotyping repeat indels in long homopolymers at the alternative knob** (run ≥ 12 bp, stutter
≥ 12% or allele fraction ≤ 0.40 het): 0.67 with the panel, 0.58 with mpileup, mostly true
insertions called stutter by the small model.

**Conclusion of the round-4 pilots.** On single-site review of real HG002 reads, from the classes
where short-read callers fail most (low-fraction mosaic, caller false positives, segmental
duplications, homopolymer indels), current models are at 0.94–0.96 with either a text summary or
raw mpileup. A 0.70 text baseline exists only for long-homopolymer indel genotyping.
