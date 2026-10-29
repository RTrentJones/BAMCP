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

**Primary contrasts (both pre-registered):**
- **H1 (V+Traw − Traw): +0.038, CI −0.026 to +0.103, p = 0.375. Not significant.**
- **H2 (V − Tjson): −0.038, CI −0.167 to +0.077, McNemar 6 vs 9 discordant, p = 0.61. Not
  significant.**

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
  (26 truth variants: 10 indels, 4 SNVs in repeats and 12 SNVs outside repeats; 8 noisy negatives where ≥ 3 reads
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
CI −0.167 to −0.028, p = 0.006; R −0.120, CI −0.204 to −0.046, p = 0.001), and beats structured
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

## Round 4 — harder items (pre-registered, focused on repeat-run indels)

Goal: find real cases where the text baseline sits near 0.70, so a gain has room to show.
Design: [`visual-design-r4.md`](visual-design-r4.md) (built by a panel of researcher and ML
agents, two adversarial reviews). Pilots: [`visual-pilot-r4.md`](visual-pilot-r4.md).
Pre-registration: [`visual-preregistration-r4.md`](visual-preregistration-r4.md).

**Pilots (text baseline only, 52 items, 156 runs each).** On real HG002 sites from the classes
where short-read callers fail most — low-fraction mosaic SNVs against caller false positives,
segmental-duplication variants against paralog look-alikes, repeat-run indels against stutter —
the models scored **0.961** with an expert text summary and **0.944** with raw
`samtools mpileup`. The 0.70 regime does not exist there. The one class near it under both
inputs was indel genotyping in long homopolymers / STRs (0.58–0.67).

**Main run (60 fresh items, 540 runs, 3 models).** 30 true indels (1–4 bp, runs ≥ 12 bp, stutter
≥ 12% or allele fraction ≤ 0.40 het) and 30 matched stutter negatives; baseline = unfiltered
mpileup of the site ± 30 bp. Audit: 540/540 answered, 0 violations (240 runs first stopped on a
usage limit before answering and were re-run with the identical prompt; one run was answered
twice, identically).

| arm | primary (balanced) | sensitivity | specificity | small / mid / large |
|---|---:|---:|---:|---|
| `Tmp` (mpileup) | **0.839** | 0.711 | 0.967 | 0.783 / 0.900 / 0.833 |
| `Tmp+Tcur` (+ BAMCP text tools) | 0.817 | 0.633 | 1.000 | 0.700 / 0.900 / 0.850 |
| `V+Tmp` (+ BAMCP viewer) | 0.811 | 0.633 | 0.989 | 0.667 / 0.917 / 0.850 |

**Confirmatory contrasts (Holm over two, cluster sign-flip test):**

| | delta | 95% CI | p | Holm p | TOST ±0.07 | |
|---|---:|---|---:|---:|---|---|
| C1 `Tmp+Tcur` vs `Tmp` | −0.022 | −0.086 to +0.040 | 0.60 | 0.82 | p = 0.069 | inconclusive |
| C2 `V+Tmp` vs `Tmp` | −0.028 | −0.083 to +0.025 | 0.41 | 0.82 | p = 0.064 | inconclusive |

Both null, and neither is shown equivalent (TOST just misses). Both additions cost sensitivity
(0.711 → 0.633) and both fail the pre-registered guardrail; the loss is concentrated in the
small model (true-indel accuracy 0.60 → 0.40 with the BAMCP text, 0.37 with the viewer), while
the mid model's true-indel accuracy falls a little (0.87 → 0.80 and 0.83) and the large model's
rises (0.67 → 0.70 in both). C1's predicted direction (no gain or harm,
because the shipped curation summary raises a HIGH homopolymer concern on every indel in a run
≥ 4 bp) is consistent with the result. The baseline landed at 0.84, not the 0.71 the pilot
cells suggested: the pilot's 12-run cells overstated the difficulty.

**Secondary: the viewer improved genotype calls.** Counting a true indel as fully right only
when it is called real with the right zygosity, variant type and allele (`genotype_correct`), the viewer arm got **27 / 90** against
**13 / 90** for mpileup alone (exact McNemar on item × model pairs, 18 vs 4 discordant,
p = 0.004; the BAMCP text arm 13 / 90). The gain is in the mid and large models (6 → 12,
4 → 13). Mechanism: in a long run, mpileup's per-column counts make a heterozygous insertion
look like a ~30% minority, and models called it `low_vaf` (30 of the 41 `Tmp` runs that gave a
true het indel a zygosity); with the picture they called it `het` (23 of 32). Homozygous indels
were mostly called `het` in every arm (22 of 30 runs in each; `hom_alt` only 2, 0 and 4 times),
because reads that stop inside the run dilute the fraction. This is a pre-specified secondary
endpoint, not a confirmatory one, and it is the first picture-specific gain in four rounds.

### Round 4 follow-up — read-based homopolymer risk (exploratory)

Round 4 pointed at two product changes, now made in `analysis/evidence.py` and
`analysis/curation.py`: an indel's homopolymer concern is decided from the reads that span the
whole run (raised only when the called length is rare among them), and the curation summary
reports the allele fraction among those reads, with the zygosity hint taken from it. The
thresholds (≥ 8 spanning reads; concern below 20%, HIGH below 10%, or when no commoner than
other non-reference lengths) were set from the mechanism before any model run below.

Re-measured on the same 60 round-4 items: the `Tmp+Tcur` arm only, 180 blind runs on the final
code (audit 180/180, 0 violations), with each item's `mpileup.txt` byte-identical to round 4 and
the round-4 `Tmp` answers reused, so the only change between arms is `tcur.md`. Review fixes made
after the first pass changed the text on 56 of the 60 items (mostly the spanning-read counts, once
spanning reads had to pass the same base-quality floor as the VAF); every run whose text changed
was redone on the final text, and the rest keep answers given on byte-identical text. The first
pass, before those fixes, scored 0.933.

| arm | primary | sensitivity | specificity | true indels: called real *and* right genotype (`genotype_correct`) |
|---|---:|---:|---:|---:|
| `Tmp` (mpileup) | 0.839 | 0.711 | 0.967 | 13 / 90 |
| `Tmp+Tcur`, round 4 (length-only rule) | 0.817 | 0.633 | 1.000 | 13 / 90 |
| `V+Tmp`, round 4 (viewer) | 0.811 | 0.633 | 0.989 | 27 / 90 |
| **`Tmp+Tcur`, read-based rule** | **0.911** | **0.833** | **0.989** | **42 / 90** |

Against mpileup alone: +0.072 (95% CI +0.011 to +0.131, sign-flip p = 0.036), with the
guardrail holding (sensitivity +0.122, specificity +0.022); all three models improve (small
0.783 → 0.867, mid 0.900 → 0.917, large 0.833 → 0.950). Against the round-4 text arm the new
text was right on 21 runs the old one missed and wrong on 4. Genotype calls on true indels:
42/90 against 13/90 (McNemar 33 vs 4 discordant, p < 0.001).

**Caveats.** This is exploratory, not pre-registered: the items are the ones whose round-4
results motivated the change, the comparison arms come from an earlier session, and it is one
run per item and model. It is the first BAMCP output in four rounds to beat the text baseline;
a pre-registered run on fresh main-partition items (the repeat pools hold several hundred
unused positives) is the way to confirm it.

## Round 5 — confirmation on fresh items (pre-registered)

Pre-registration: [`visual-preregistration-r5.md`](visual-preregistration-r5.md), committed
(36513ba) before the item draw; manifest frozen (7866bbb) before any run. 90 sites no earlier
round has seen: 45 true repeat-run indels and 45 matched stutter negatives from chr15–22, none
within 10 kb of a round-4 site and none in a round-4 segmental-duplication cluster. Same rules,
prompt, caption and models as round 4; BAMCP at `main` (24538b8). Arms: `Tmp` (mpileup) and
`Tmp+Tcur` (+ BAMCP text tools with the spanning-read rule). 540 blind runs; audit 540/540
answered, 0 violations, no re-runs.

| arm | primary (balanced) | sensitivity | specificity | small / mid / large | true indels: real *and* right genotype |
|---|---:|---:|---:|---|---:|
| `Tmp` (mpileup) | 0.826 | 0.659 | 0.993 | 0.756 / 0.867 / 0.856 | 30 / 135 |
| `Tmp+Tcur` (spanning-read rule) | **0.885** | **0.770** | **1.000** | 0.822 / 0.922 / 0.911 | **54 / 135** |

**Confirmatory (Holm over two):**

| | effect | 95% CI | p | Holm p | |
|---|---:|---|---:|---:|---|
| H1 balanced accuracy | +0.059 | +0.009 to +0.109 | 0.027 | 0.027 | **supported** (guardrail holds) |
| H2 genotype on true indels | 54 vs 30 of 135 | discordant 33 vs 9 | 0.0003 | 0.0005 | **supported** |

Both registered hypotheses are supported on fresh items. The primary effect is smaller than
the follow-up's +0.072, as expected for a re-test, and lies inside the follow-up's interval;
TOST at ±0.07 does not apply to a supported effect. Each model improves by a similar amount
(+0.067, +0.056, +0.056), and none of the per-model differences is significant on its own.
Sensitivity rises by 0.111 and specificity does not fall (0.993 → 1.000).

**Where the gain is.** It is in heterozygous true indels: joint accuracy 0.52 → 0.69, with
correct genotypes 20 → 38 of 90. Among runs that called a true het indel real, models without
BAMCP's text gave it `low_vaf` 27 times and `het` 20 times; with it, 24 and 38. Homozygous indels stay at 0.93 joint
accuracy in both arms, and most correct calls on them still say `het` (29 and 24 of 42),
because reads that stop inside a long run dilute every fraction the models see. Insertions
gain more than deletions in genotype (18 → 39 of 93 against 12 → 15 of 42); joint accuracy
rises by a similar amount in both (0.70 → 0.80 and 0.57 → 0.71).

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
- **Harder real items do not open a gap (round 4).** On hard classes of real HG002 sites the
  text baseline stays at 0.94–0.96; only long-repeat indel genotyping is hard (0.84 on fresh
  items). There, neither BAMCP's text tools nor the viewer raised joint accuracy (both −0.02 to
  −0.03, inconclusive), both cost the small model sensitivity, and the viewer raised
  correct genotype calls for the mid and large models from 10 to 25 of 60 (secondary, p = 0.004). After the
  product fix those results pointed to (read-based homopolymer risk, spanning-read allele
  fraction), BAMCP's text tools raised accuracy on the same items from 0.839 to 0.911 and
  correct genotypes from 13/90 to 42/90 (exploratory).
- **Confirmed on fresh items (round 5).** On 90 sites no earlier round had seen, BAMCP's text
  tools beat raw mpileup on both pre-registered endpoints: balanced accuracy +0.059 (0.826 →
  0.885, Holm p = 0.027) and correct genotypes on true indels 30 → 54 of 135 (Holm p < 0.001),
  with specificity not lower (0.993 → 1.000). This is the first confirmed gain from a BAMCP output in five
  rounds, and it comes from the text tools, not the viewer.
- For BAMCP this supports the current design: the viewer is for the person; the model should
  keep getting the read evidence as text. Two product changes followed from round 4 and are
  now made (#71): the curation summary's homopolymer concern depends on the reads that span the
  run, not only its length (design item J2), and the text tools report allele fraction among
  those spanning reads, which is the information the picture was supplying for zygosity.

## Corrections (2026-10-29)

An audit of this document against the graded run sets and the committed manifests changed the
following; no result or conclusion changed. The pilot's second pre-registered primary contrast
(H2, `V` vs `Tjson`) was missing and is now reported; the round-3 real-tier `V` vs `Traw` upper
CI bound was −0.056 and is −0.046; round 4's mid-model true-indel accuracy fell slightly rather
than staying flat; the round-4 zygosity counts and the round-5 het counts now state the runs
they are counted over; "called `het` in every arm" now gives the counts (22 of 30 per arm);
"roughly doubled" is replaced by the counts (10 → 25 of 60); the insertion-vs-deletion sentence
no longer claims a larger joint-accuracy gain for insertions; the genotype endpoint is named as
`genotype_correct` (zygosity, type and allele), which gives the same counts as zygosity alone in
every arm; the round-3 truth composition is stated exactly; and the round-4 product changes are
described as made.

