# Visual-vs-text eval — results (2026-10)

Method and analysis plan: [`visual-preregistration.md`](visual-preregistration.md) (committed before
any model run). Harness: `python -m bamcp.eval.visual build|grade|lint|dry`.

**Setup.** 26 synthetic items over 18 case types (real SNV/indel/MNV/multi-allelic/phased variants,
seven artifact classes, clean negatives). Three Claude model tiers (small/mid/large), one blind run
per item × condition × model. Each run was a fresh subagent allowed only to read its own stimulus
files and return the fixed answer schema; labels were sealed until all runs finished. Transcript
audit: 0 violations across 468 runs. Metric: joint accuracy (both "is a variant present" and "real or artifact" answered correctly), pooled
over models, 95% CIs by item-cluster bootstrap.

## Pilot (pre-registered, original viewer) — 312 runs

| Condition | Joint accuracy | 95% CI | Brier |
|-----------|---------------:|--------|------:|
| `Traw` — raw read text | 0.949 | 0.885–1.000 | 0.053 |
| `Tjson` — BAMCP structured evidence | 0.859 | 0.744–0.949 | 0.116 |
| `V` — viewer image only | 0.821 | 0.705–0.923 | 0.147 |
| `V+Traw` — image + raw text | 0.987 | 0.962–1.000 | 0.051 |

**Primary contrast (V+Traw − Traw): +0.038, CI −0.026 to +0.103, p = 0.375. Not significant.**
Text alone is near the ceiling on these items.

The image-only arm failed hom-alt SNVs (0.17) and low-base-quality artifacts (0.33). Root cause
was rendering: every base was coloured and base colours collided with strand colours, so a column
of identical mismatches did not stand out; low-BQ mismatches were drawn like high-BQ ones.

## Iteration 2 (exploratory, fixed viewer) — 156 runs

Fix: neutral read bodies, colour only on mismatches, strand by arrow/underline, low-BQ mismatches
drawn hollow (`static/renderer.ts`, `VIEWER_KEY`). Same items and run ids, re-rendered.

| Condition | Pilot | Fixed viewer | Paired change |
|-----------|------:|-------------:|---------------|
| `V` | 0.821 | **0.910** | +0.090 (CI 0.000 to +0.192), 9 up / 2 down, McNemar p = 0.065 |
| `V+Traw` | 0.987 | 0.962 | −0.026 (CI −0.077 to +0.026), 1 up / 3 down, p = 0.63 |

- Image-only by model tier: small 0.692 → 0.769, mid 0.885 → 0.962, large 0.885 → **1.000**.
- Case types: hom-alt SNV 0.17 → 0.83, low-BQ artifact 0.33 → 0.67, phased-trans and scattered
  errors 0.67 → 1.00; homopolymer deletion 0.83 → 0.67.
- Image-only Brier 0.147 → 0.094, but it stays overconfident when wrong (mean confidence 0.66 on
  wrong answers vs 0.78 on right ones). Raw text, with or without the image, calibrates best.

## Reading these results

- Legibility of the rendering mattered more than the presence of an image.
- With the fix, a picture alone beats BAMCP's structured JSON (0.910 vs 0.859) and approaches raw
  text (0.949); the largest model is perfect from the picture alone.
- Adding the image to raw text gives no measurable gain at this ceiling. BAMCP therefore keeps
  returning text evidence alongside the viewer.

**Limits.** Iteration 2 was designed after seeing which pilot items failed and was measured on
those items, so its gain is an upper estimate. All data is synthetic. n is small (26 items, one
run per cell), several case types have a single item, and text arms are near the ceiling.

**Next.** A pre-registered test of the fixed viewer on fresh items, including GIAB
difficult-region sites from real data, with repeats per item.

## Round 3 (pre-registered confirmation) — 756 runs

Pre-registration: [`visual-preregistration-r3.md`](visual-preregistration-r3.md), committed
(4cb4048) before any run. Fixed viewer (`ae5111b`). Two tiers, analysed separately:

- **S — fresh synthetic:** 36 never-seen items (seed 2, all 18 case types × 2), `Traw`, `Tjson`,
  `V`; 324 runs.
- **R — real reads:** 36 sites from the GIAB HG002 60x BAM on chr20, truth from GIAB v4.2.1
  (26 truth variants, half of them indels or SNVs in repeats; 8 noisy negatives where ≥ 3 reads
  disagree with a reference call; 2 clean negatives), blinded onto re-indexed contigs; `Traw`,
  `Tjson`, `V`, `V+Traw`; 432 runs. Site manifest [`data/r3_real_sites.jsonl`](data/r3_real_sites.jsonl).

Audit: 756/756 answered, 1 invalid answer (scored wrong), one rule deviation (a run used Grep on
its own listed `tjson.json`; not a read outside its files, so kept per the exclusion rule;
excluding it changes nothing reported below).

| Tier | `Traw` | `Tjson` | `V` | `V+Traw` |
|------|------:|------:|----:|------:|
| S (fresh synthetic) | 0.981 | 0.824 | 0.889 | — |
| R (real HG002) | 0.991 | 0.926 | 0.870 | 0.981 |

**Confirmatory contrasts (Holm over three):**

| | delta | 95% CI | McNemar p | Holm p | |
|---|---:|---|---:|---:|---|
| C1 S: `V` vs `Tjson` | +0.065 | −0.019 to +0.148 | 0.167 | 0.50 | not supported |
| C2 R: `V+Traw` vs `Traw` | −0.009 | −0.028 to 0.000 | 1.0 | 1.0 | not supported |
| C3 R: `V` vs `Tjson` | −0.056 | −0.139 to +0.019 | 0.210 | 0.50 | not supported |

**Secondary (descriptive):** raw read text beats the picture alone on both tiers (S −0.093,
CI −0.167 to −0.028, p = 0.006; R −0.120, CI −0.204 to −0.056, p = 0.001), and beats structured
JSON (S −0.157; R −0.065). Image-only accuracy by model tier: S 0.806 / 0.861 / 1.000,
R 0.694 / 0.944 / 0.972 (small / mid / large). Among real variants, `V` got zygosity right 85%
and the allele 80% of the time, against 96% and 99% for `Traw`. Brier: `Traw` 0.042 / 0.034,
`V+Traw` 0.033 (R), `V` 0.115 / 0.107. Image-only answers stayed overconfident when wrong
(mean confidence 0.74 on S and 0.69 on R for wrong answers, against 0.78–0.81 for right ones).

Where the picture alone failed on real data: homopolymer insertions read as polymerase stutter
(one item missed by all three models), and three homozygous SNVs read as "every read matches
the reference". The blinded reference and render are correct there (a solid column of alt tiles
in every read), so this is a reading error, plausibly because the reference track is fully
coloured and a column coloured in every read looks like more reference.

## Conclusion across rounds

- On these items, **a picture of the pileup does not beat the reads as text**. Raw read text is
  the strongest single input on synthetic and on real data, and adding the image to it neither
  helps nor hurts measurably (two pre-registered tests, both null at the ceiling).
- **The picture alone is a usable but weaker channel**: 0.87–0.89 joint accuracy, close to
  BAMCP's structured JSON and well behind raw text, and overconfident when wrong. The round-2
  gain from the legibility fix is real in direction (fresh items: `V` 0.889 vs `Tjson` 0.824)
  but did not reach significance, so it is not confirmed.
- **Model size matters most for images**: the largest model was near-perfect from pictures
  (1.000 synthetic, 0.972 real); the smallest was not (0.81, 0.69).
- For BAMCP this supports the current design: the viewer is for the person; the model should
  keep getting the read evidence as text.
