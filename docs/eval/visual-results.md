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
