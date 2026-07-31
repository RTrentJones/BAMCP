# GIAB benchmark — results

Running BAMCP's variant detection against **real NIST Genome-in-a-Bottle truth
data** for HG001 / NA12878 (v4.2.1, GRCh38). This is the real-biology
counterpart to the deterministic `synthetic_v1` gate: `synthetic_v1` proves the
harness and guards regressions on every PR; this shows the caller's behavior on
authentic variant coordinates, alleles, and zygosity.

Reproduce (or just `make giab-benchmark`, which pins every generation
parameter):

```bash
python tests/eval/datasets/giab/fetch_giab.py --region chr20:1000000-1060000 \
  --depth 30 --read-len 100 --error-rate 0.001 --seed 1
python -m bamcp.eval.truthset --manifest tests/eval/datasets/giab/manifest.yaml \
  --min-vaf 0.2 --min-depth 8 --min-mapq 20 --max-reads 50000 \
  --min-recall 0.90 --min-precision 0.90
```

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

Same 66 variants on both lists, one apart. The manifest is now regenerated
1-based (`TruthSNV.pos1`), the results below reproduce exactly as originally
published, and `test_manifest_positions_are_one_based` fails against the old
writer so this cannot recur.

Note the asymmetry that made this easy to get wrong, which still holds: BAMCP's
region *strings* are 0-based half-open (`parse_region` hands `start`/`end`
straight to pysam), while reported variant *positions* are 1-based. So
`detection_region` and `negative_regions` in this manifest are correct as
0-based; only `variant_sites[].pos` needed the shift.

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
| **A — default (curation sensitivity)** | 0.02 | 2 | 0 | **1.000** | 0.034 | 0.066 | 171 |
| **B — variant calling** | 0.20 | 8 | 20 | **1.000** | **1.000** | **1.000** | **0** |

- **Recall is 1.000 at both points**: all 66 real NIST SNVs are recovered,
  including the 7 heterozygous sites at ~50 % VAF.
- **Precision depends entirely on the threshold.** BAMCP's *default* thresholds
  (`min_vaf=0.02`, `min_depth=2`, `min_mapq=0`) are deliberately tuned for
  curation sensitivity — surface everything, let a human review — so at Q30
  noise they emit many low-VAF noise calls (precision 0.034). Thresholded for
  *calling* (VAF ≥ 0.20, depth ≥ 8, MAPQ ≥ 20), the same tool on the same reads
  reaches **precision 1.000 with no loss of recall**.

The practical takeaway: BAMCP's caller is correct, but its defaults are a
max-sensitivity operating point. Anything doing automated calling (rather than
human-in-the-loop curation) should raise the thresholds — the `--min-vaf` /
`--min-depth` / `--min-mapq` knobs on the scorer exist to make the operating
point explicit and measurable.

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
