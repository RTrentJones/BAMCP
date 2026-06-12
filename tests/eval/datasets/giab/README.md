# GIAB truth set (real-data extension)

The `synthetic_v1` dataset proves the scorer works and gates regressions
deterministically in CI. This directory is the path to *real* ground truth:
the [Genome in a Bottle](https://www.nist.gov/programs-projects/genome-bottle)
consortium's high-confidence small-variant calls for the NA12878 / HG001
reference sample.

Real GIAB data is large and network-gated, so it is **not** fetched in CI.
The synthetic gate runs on every PR; a GIAB run is a manual / nightly job.

## What you need

| File | Source | Purpose |
|------|--------|---------|
| HG001 truth VCF | GIAB v4.2.1 release | gold-standard variant calls |
| HG001 high-confidence BED | GIAB v4.2.1 release | regions where the truth set is callable |
| An aligned BAM/CRAM for HG001 | e.g. a downsampled NIST/Illumina run | the reads BAMCP actually inspects |
| Matching reference FASTA | GRCh38 (or GRCh37) | required for variant calling |

GIAB release index:
`https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/NA12878_HG001/`

## How it plugs in

`fetch_giab.py` downloads a release, verifies checksums, and slices a single
small region (default: a few hundred kb on chr20) so the scorer stays fast.
It then emits a `manifest.yaml` in the *same schema* as
`../synthetic_v1/manifest.yaml`:

- `variant_sites` ← every truth VCF record inside the slice that also falls in
  the high-confidence BED. The `expected_artifact` / `is_clean_control` fields
  are left null (GIAB labels truth, not artifacts), so the GIAB run reports
  variant-detection precision/recall/F1 while the synthetic run additionally
  covers artifact-type recall.
- `negative_regions` ← high-confidence stretches in the slice with no truth
  variant, used for false-positive control.

Once the manifest exists, the existing scorer runs unchanged:

```bash
python -m bamcp.eval.truthset \
  --manifest tests/eval/datasets/giab/manifest.yaml \
  --min-recall 0.95 --min-precision 0.90
```

The lower floors (recall ≥ 0.95, precision ≥ 0.90 rather than 1.0/0.95)
reflect that BAMCP is a lightweight pileup caller, not a production variant
caller — the point is to track its real-data behavior over time and catch
regressions, not to beat DeepVariant.

## Status

`fetch_giab.py` is a scaffold: it documents the exact URLs, the checksum
verification step, and the manifest-emit logic. Wiring the slice + VCF parse is
the first task of the GIAB milestone (see `EVALS.md`).
