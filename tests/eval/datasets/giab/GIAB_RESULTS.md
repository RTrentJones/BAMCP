# GIAB benchmark — results

Running BAMCP's variant detection against **real NIST Genome-in-a-Bottle truth
data** for HG001 / NA12878 (v4.2.1, GRCh38). This is the real-biology
counterpart to the deterministic `synthetic_v1` gate: `synthetic_v1` proves the
harness and guards regressions on every PR; this shows the caller's behavior on
authentic variant coordinates, alleles, and zygosity.

Reproduce point B (or just `make giab-benchmark`, which pins every generation
parameter and scores point B only):

```bash
python tests/eval/datasets/giab/fetch_giab.py --region chr20:1000000-1060000 \
  --depth 30 --read-len 100 --error-rate 0.001 --seed 1
python -m bamcp.eval.truthset --manifest tests/eval/datasets/giab/manifest.yaml \
  --min-vaf 0.2 --min-depth 8 --min-mapq 20 --max-reads 50000 \
  --min-recall 0.90 --min-precision 0.90
```

Point A scores the same reads at the curation-sensitivity thresholds (the manifest
pins no operating point, so every threshold is given explicitly; `--min-alt-count 1`
is required because the current default is 3):

```bash
python -m bamcp.eval.truthset --manifest tests/eval/datasets/giab/manifest.yaml \
  --min-vaf 0.02 --min-depth 2 --min-mapq 0 --min-alt-count 1 \
  --min-recall 0.90 --min-precision 0.0
```

`fetch_giab.py`'s `--region` is a 0-based half-open window; the region strings it
writes into the manifest are 1-based inclusive, as BAMCP reads them.

## Coordinate-basis correction (2026-07-26)

**Every number below was unreproducible between #29 and this correction.** The
numbers themselves were correct when first published — at that time both the
manifest and `get_variants` reported 0-based positions, so they agreed. #29
("report 1-based positions from the variant tools") moved the tool output to
1-based and migrated `synthetic_v1`'s manifest, but not this dataset's. Because
the scorer matches `(chrom, pos, ref, alt)` **exactly, with no tolerance**, the
two sides then disagreed by one at every single site.

The failure was silent in the worst way — it looked like a caller regression
rather than a coordinate-space mismatch:

```
variant detection  P=0.000 R=0.000 F1=0.000  (tp=0 fp=66 fn=66)
  missing:  chr20:1002294:G>A ...
  spurious: chr20:1002295:G>A ...
```

Same 66 variants on both lists, one apart. The manifest was regenerated
1-based (`TruthSNV.pos1`), the results reproduced exactly as originally
published, and `test_manifest_positions_are_one_based` fails against the old
writer so this cannot recur.

## Region-string correction (2026-10-29)

The same mismatch then recurred for the region strings. #63 made BAMCP's region
*strings* 1-based inclusive (`parse_region` subtracts one from the start), but
`fetch_giab.py` still wrote `detection_region` and `negative_regions` 0-based, so
each window started one base early. Two negative regions then began on a
homozygous truth SNV (chr20:1008797 A>T, chr20:1011613 C>T), and `make
giab-benchmark` reported `false positives 2 (negative regions)` and failed, with
precision and recall unchanged. The generator now writes 1-based strings
(`chr20:1000001-1060000`; negatives `chr20:1000001-1002000`,
`chr20:1008798-1010797`, `chr20:1011614-1013613`), the regenerated manifest differs
only in those four strings, and `test_manifest_region_strings_are_one_based` pins
them. Re-measured on 2026-10-29, both rows below reproduce exactly.

## Dataset

| | |
|---|---|
| Sample | HG001 / NA12878 |
| Truth set | NIST GIAB v4.2.1, GRCh38, high-confidence BED |
| Reference | UCSC hg38 chr20 |
| Region | chr20:1,000,000–1,060,000 (60 kb) |
| Truth SNVs (in high-conf BED) | 66 (59 hom, 7 het) |
| Indel truth records in BED | 7 — skipped, SNV-only (reported by the generator) |
| Reads | simulated, 30× coverage, 100 bp, base-error rate 0.001 (~Q30), seed 1 |

**This is semi-synthetic over real biology.** The reference sequence and the
variant catalog (positions, alleles, genotypes) are real GIAB data; the reads
are simulated honoring each site's genotype (het ≈ 50 % alt, hom ≈ 100 % alt)
with a uniform Illumina-like base-error rate. See [Limitations](#limitations).

## Results

The headline finding is a **precision/recall tradeoff driven by the caller's
operating point**, measured on the same reads:

| Operating point | min VAF | min depth | min MAPQ | Recall | Precision | F1 | FP (neg. regions) |
|---|---|---|---|---|---|---|---|
| **A — curation sensitivity (pre-2026-10 defaults)** | 0.02 | 2 | 0 | **1.000** | 0.034 | 0.066 | 171 |
| **B — variant calling** | 0.20 | 8 | 20 | **1.000** | **1.000** | **1.000** | **0** |

- **Recall is 1.000 at both points**: all 66 real NIST SNVs are recovered,
  including the 7 heterozygous sites at ~50 % VAF.
- **Precision depends entirely on the threshold.** Point A (`min_vaf=0.02`,
  `min_depth=2`, `min_mapq=0`, `min_alt_count=1`, BAMCP's defaults before
  2026-10) is a curation-sensitivity point — surface everything, let a human
  review — so at Q30 noise it emits many low-VAF noise calls (precision 0.034).
  Thresholded for *calling* (VAF ≥ 0.20, depth ≥ 8, MAPQ ≥ 20), the same tool on
  the same reads reaches **precision 1.000 with no loss of recall**.
- The current server defaults (`min_vaf=0.10`, `min_depth=2`, `min_mapq=20`,
  `min_baseq=13`, `min_alt_count=3`) are neither point and were not measured here.

The practical takeaway: on this slice the caller recovers every truth SNV, and
precision is set by the operating point. The scorer's `--min-vaf` /
`--min-depth` / `--min-mapq` / `--min-alt-count` flags make the operating point
explicit and measurable.

## Limitations

Be precise about what these numbers do and do not say:

- **They measure** recovery of real variant structures on real sequence context
  (recall) and rejection of uniform sequencing noise at a given threshold
  (precision).
- **They do not measure** robustness to *real* read data: real base-quality and
  error profiles, mismapping in segmental duplications, reference bias, or
  alignment artifacts. Those require a real HG001 alignment (a ~300 GB file),
  which is out of scope for a portable, reproducible benchmark.
- **SNVs only.** Indel read-simulation is not implemented, so the **7** indel
  truth records inside this BED slice are skipped rather than mis-simulated.
  The generator now prints that count — previously the skip was real but the
  count was not, so the truth set read as complete when it is SNV-only.
- Because reads are simulated cleanly, precision at point B is an *upper bound*
  — real reads would introduce failure modes this benchmark cannot see.

## Next step toward real reads

`fetch_giab.py` already emits a manifest the scorer consumes unchanged. Swapping
the simulated BAM for a real HG001 alignment (via BAMCP's remote-BAM support,
slicing just this region through the `.bai`) would turn this into a
fully-real-data benchmark. That is the natural follow-up; the harness is ready
for it.
