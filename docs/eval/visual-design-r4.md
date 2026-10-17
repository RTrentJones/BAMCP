# Visual-vs-text eval, round 4: final design (hard real strata, baseline about 0.70)

**Status:** FINAL DESIGN (revision of `design_draft.md` after two adversarial reviews). This becomes the pre-registration (`docs/eval/visual-preregistration-r4.md`) only after the pilot passes its go/no-go gate (section 7) and the freeze list (section 12) is committed.

Design-phase evidence: the power simulation is committed as [`r4/power_r4b.py`](r4/power_r4b.py) with its output [`r4/power_r4b_out.txt`](r4/power_r4b_out.txt) (balanced-accuracy primary, stratum-specific effects, 2,000 simulations). The site-count probes cited below were run against the GIAB v5.0q smvar, CMRG v1.00 and Mosaic SNV v1.1 releases and the v3.x stratifications; their scripts are rebuilt as `truth.py`/`strata_r4.py` in the implementation.

**What changed from the draft (summary).** Orthogonal-technology (HiFi) and trio label checks on every hard stratum; mosaic concordance rule fixed; B gains mate/pairing fields, 5′ distance and a linkage table, so it no longer withholds read-level information; viewer gains soft clips and pair information (V5, V6); primary endpoint becomes the mean over 4 pos/neg pairs of balanced accuracy, with a sensitivity/specificity guardrail; quotas frozen before the pilot; larger B-only pilot with two knob levels; no treatment runs in the pilot; recoverability gate G8; cluster-aware sampling and testing; power recomputed for the predicted (mixed-sign) effects; negatives sourced from caller false positives where possible.

---

## 1. Why round 3 hit the ceiling, and what round 4 changes

1. **`Traw` is not a raw read listing.** `stimuli.traw_text` hands the model a site panel resembling bam-readcount output: per-allele counts split by strand, MAPQ, base quality and 3′-distance histograms, soft-clip counts, and ±15 bp context. That is an expert's artifact checklist.
2. **`Tjson` is the blinded read payload, not BAMCP's evidence.** `blind_payload` empties `variants` and strips `variant_evidence`. BAMCP's evidence layer (`get_variant_curation_summary`, `get_region_summary`) has never been tested.
3. **v4.2.1 confident regions are easy by construction.** In 46 reference-confident homopolymers of 12 bp or more, worst stutter was below 8%. Segdup truth SNVs inside v4.2.1 had MAPQ 60 and VAF 0.30-0.58. Of 210 v5.0q-only segdup reference positions on chr20, none had ≥15% non-reference reads.

Round 4 takes its difficulty from **real sites where site-local evidence is genuinely ambiguous**: truth sets covering what v4.2.1 excludes (Mosaic SNV v1.1, CMRG v1.00, v5.0q), real low-VAF biology, and negatives that mimic signal (caller false positives, PSVs, misalignment beside an indel).

**Round 4 does not weaken the text; it strengthens it.** No truncation, downsampling or depth budget. Because difficulty must not come from information withheld from text, the baseline `B` = **`Traw-r4`** adds what a careful command-line user gets cheaply (section 4). The deferred hard cases (low depth, soft-clip-only long indels, window triage) are listed at the end of section 13.

## 2. Research question and claim framing

**RQ.** On real HG002 short-read candidate sites from technically challenging classes, where the best realistic plain-text site summary lands near 0.70 on the primary endpoint, does adding (C1) BAMCP's structured curation evidence or (C2) the BAMCP viewer image raise LLM joint accuracy on variant presence plus the real-vs-artifact verdict?

Low-complexity plus segdup sites make up about 77% of technically challenging P/LP findings in Lincoln et al. 2021.

**Claim scope.**
- "On candidate sites with a pileup signal in these hard classes (caller-flagged where noted in section 3), adding X to an expert site summary improves balanced joint correctness by Δ [95% CI]." Not a general claim that images help variant review.
- Rounds 1-3 remain the record for easy sites.
- **External yardstick:** Barnell et al. 2019, manual IGV review 77.4% → 94.1% with a structured SOP (+16.7 points).
- **Prevalence:** class balance is set **per pair** to model review of flagged candidates, not genome-wide prevalence. There is no "50:50 overall" claim; sensitivity and specificity are reported separately.
- **Comparability:** new caption sentence and a strengthened `B`; absolute accuracies are not comparable to rounds 1-3.
- **Provenance of the features:** J1, J3 and J5 were motivated partly by round-3 cells of n = 6 (e.g. strand_bias Tjson 3/6, homopolymer_del V 2/6) that may be noise. Stratum-level predictions are hypotheses, not confirmations.

## 3. Strata (main set: 128 items = 4 pairs × (15 pos + 15 neg) + 8 anchors)

### 3.0 Rules common to every stratum

- **Reads.** HG002 60x GRCh38 Illumina (`HG002.GRCh38.60x.1.bam`, NHGRI Novoalign, remote, local `.bai`). Primary, non-duplicate, non-QC-fail. Blinded and re-indexed as in round 3 (`write_real_fixtures`: 3 kb slots, hashed names), **but paired-end information is now kept** (V6): FLAG proper-pair / mate-unmapped / mate-reverse bits, the `MQ` (mate MAPQ) tag, and mate position re-indexed if the mate lies in the slot, otherwise a sentinel "mate elsewhere" contig with no real coordinate. `SA` reduced to a presence flag. Before the pilot, confirm from the BAM header/README that the library derives from NIST RM 8391 (the mosaic VAFs differ in Coriell DNA).
- **Orthogonal label check (new, frozen before the pilot).** Every item in strata 1-8 must pass a check against **HG002 PacBio HiFi** (`PacBio_CCS_15kb_20kb_chemistry2/GRCh38/HG002.SequelII.merged_15kb_20kb.pbmm2.GRCh38.haplotag.10x.bam` on ftp-trace; probed: opens in 6.5 s, about 55-63x depth, 2-3 s per site; all 4 mosaic sites probed showed the alt: 8/55, 12/60, 4/51, 9/63):
  - SNV positives: HiFi alt reads at MAPQ ≥ 20 ≥ 3 and fraction ≥ 0.03 (mosaic), ≥ 0.25 (het), ≥ 0.80 (hom).
  - TR-indel positives: the HiFi spanning-read length histogram has its modal length(s) at the truth genotype.
  - Negatives: HiFi carries the candidate allele in ≤ 1 read and ≤ 2%. For TR stutter negatives (stratum 4), HiFi's ±1-unit fraction must be ≤ half of the Illumina fraction.
  - **Trio recurrence** (HG003/HG004 `NIST_HiSeq_*_Homogeneity` 60x Novoalign BAMs, same pipeline), at sites where the parents' v4.2.1 truth is reference or the parental BED does not cover the site: required for strata 4 and 6 (the same non-reference allele at ≥ 5% in at least one parent marks a locus-specific technical artifact, not HG002 mosaicism); recorded for all negatives.
  - Exclusion counts per stratum and rule are reported.
- **Label-error audit (residual).** Before freezing the main set, two reviewers adjudicate 3 random main items per hard stratum (24 items) blind, from HiFi plus trio pileups. If either finds ≥ 2 confirmed errors in a stratum, the whole stratum is reviewed the same way and confirmed-error items are replaced by the next item in the seeded draw order. This happens before any main run and never uses model output. The confirmed-error rate is reported.
- **Truth-set coverage is per variant class.** A v4.2.1 or v5.0q reference call is *not* evidence about variants below 30% VAF, so the v4.2.1 "covers and agrees" check is not applied to mosaic sites. Where two truth sets both assert a **germline** genotype at a site they must agree after `bcftools norm` left-normalisation; discordant sites are dropped and counted.
- **v5.0q FILTER handling.** Positives use PASS records inside the v5.0q BED only. Isolation and negative clearance count records of **any** FILTER, inside or outside the BED (chr20 has 59,606 HET1/HET2/GAP1/GAP2 records and 25,831 PASS records outside the BED).
- **Isolation.** No other truth record (any set, any FILTER) within 30 bp of the POI; for TR strata, none anywhere in the repeat run ±30 bp, and the run holds exactly one truth record. Indel-adjacent strata use their own rule (3.7-3.8).
- **Spacing and clusters.** Cluster = connected component of the GRCh38 `genomicSuperDups` pairing (a segment and all its partners); outside segdups, a 1 Mb bin. Items must be ≥ 10 kb apart, with ≤ 2 items per cluster. The pilot/main partition is a seeded hash of the **cluster ID** (25% pilot), so paralog families such as SMN1/SMN2, PMS2/PMS2CL, GBA/GBAP1 and NCF1/NCF1B/NCF1C never straddle pilot and main.
- **Recoverability and evidence sufficiency (data only, frozen).** A truth het needs ≥ 3 alt reads (≥ 4 for mosaic) and VAF ≥ 0.10 among the reads shown (mosaic: ≥ 0.06). A truth hom needs alt fraction ≥ 0.60. Depth ≥ 10. Failing candidates are never sampled and are counted.
- **Labels.** Positive = `yes` / `real`; negative = `no` / `artifact`, `artifact_type` unscored except stratum 4.
- **Never select on outcomes.** No item is selected, kept or dropped on the basis of any model's answer (selecting failed items produced a spurious +0.06 to +0.14 with no true effect in the ML simulation).
- **Signature report.** For every pair, the distribution of strand, base-quality, MAPQ, read-end, soft-clip and linked-mismatch features of selected negatives vs positives is computed on the pool and published before the main runs.

| # | stratum | n main | n pilot (default / alt knob) | truth source + orthogonal check | expected `B` |
|---|---|---|---|---|---|
| 1 | `mosaic_snv_pos` | 15 | 6 / 4 | Mosaic SNV v1.1; HiFi alt | 0.62 |
| 2 | `mosaic_caller_fp_neg` | 15 | 6 / 4 | inside mosaic BED, not in mosaic VCF, caller-flagged; HiFi clean | 0.68 |
| 3 | `tr_indel_pos` (10 het, 5 hom) | 15 | 6 / 4 | v5.0q PASS (+v4.2.1 where it covers); HiFi length histogram | 0.72 |
| 4 | `tr_stutter_neg` | 15 | 6 / 4 | v5.0q (+v4.2.1) reference; HiFi + trio recurrence | 0.66 |
| 5 | `segdup_cmrg_pos` | 15 | 6 / 4 | CMRG v1.00 (+v5.0q); HiFi alt | 0.68 |
| 6 | `segdup_psv_neg` | 15 | 6 / 4 | CMRG BED ∩ v5.0q BED reference; HiFi clean + trio recurrence; not in variable-copy regions | 0.62 |
| 7 | `indel_adjacent_misalign_neg` | 15 | 6 / 4 | v4.2.1 ∩ v5.0q; HiFi clean | 0.70 |
| 8 | `indel_adjacent_snv_pos` | 15 | 6 / 4 | v4.2.1 ∩ v5.0q; HiFi alt | 0.82 |
| 9 | `anchor_controls` (2 het SNV, 3 hom SNV, 3 clean neg) | 8 | 4 | v4.2.1 ∩ v5.0q | 0.98 |
| | **total** | **128** | **84** | | |

**Expected primary B** (mean over pairs of balanced accuracy, anchors excluded): mosaic (0.62+0.68)/2 = 0.65; TR 0.69; segdup 0.65; indel-adjacent 0.76 → **0.69**. By model, from round-3 offsets: about 0.58 haiku, 0.69 sonnet, 0.80 opus. These expectations are guesses; the honest 95% prediction interval for the main-set primary B before the pilot is about **±0.14** (pilot SE about 0.05-0.06 with correlated model runs, main item sampling about 0.04, regression to the mean after knob moves). After the pilot it is about ±0.09.

### 3.1 `mosaic_snv_pos` (n = 15)
- **Truth.** `HG002_GRCh38_MosaicSNVv1.1.vcf.gz` + BED: 85 SNVs, VAF about 5-30% at 300x (median 0.12). Label `yes`/`real`, `zygosity: low_vaf`.
- **Selection.** ≥ 4 truth-alt reads, depth ≥ 25 at 60x. **Drop any site with a v5.0q record (any FILTER) within 30 bp**: 6 sites, including 2 that v5.0q calls 0/1 (chr1:79144396 — HiFi shows 8/55 = 15%, so the v5.0q germline call is itself suspect — and chr3:21601448). The v4.2.1 check is not applied (all 85 sites lie in the v4.2.1 BED with no record; v4.2.1 cannot see < 30% VAF). Alt-count bands: half 4-6, half ≥ 7. HiFi check.
- **Pool.** 85 − 6 (v5.0q within 30 bp) − about 7 (≤ 3 alt at 60x) ≈ 72, minus HiFi failures; about 10 go to pilot by cluster hash, leaving about 55-60 for a main quota of 15 (G5 requires ≥ 30).
- **Why hard.** 4-12 clean alt reads out of 50-80 with no strand, base-quality or MAPQ flag; the germline reflex ("not half the reads, so artifact") is wrong.

### 3.2 `mosaic_caller_fp_neg` (n = 15)
- **Truth.** Inside the mosaic v1.1 BED (asserts VAF < 2% at any unlisted SNV) and not in the mosaic VCF; no v4.2.1, v5.0q or mosaic record within 30 bp; HiFi clean.
- **Source (caller-flagged, matching the claim scope).** Candidates rejected from the benchmark: `var_cat` = putative/candidate rows of `SupplementaryFiles/HG002_GRCh38_mosaic_database.tsv` (which also carries per-platform alt depths: PacBio, Element, BGI, Duplex), `HG002_GRCh38_Strelka2-Ill300X.vcf.gz`, DRAGEN `HG002_300x.mosaic-only-hard-filtered.vcf.gz` and DeepSomatic Illumina-300x Q0 calls (all in the release's `external_validator_vcfs/`), each minus the mosaic VCF. Fallback if fewer than 2× quota survive: a pileup-signal scan inside the mosaic BED (≥ 4 reads, VAF 0.05-0.20), with the claim scope for that share restated as "pileup-signal".
- **Matching.** One-to-one to positives on alt-read bin (4-6, 7-9, ≥ 10), depth bin (< 50, 50-70, > 70) **and segdup/low-mappability membership** (the same about 4% share as the positives, 3 of 85), so neither VAF nor MAPQ separates the classes for free.
- **Population filter.** Exclude negatives at gnomAD v4 AF ≥ 1% (looked up before the draw, never shown): clean low-level signal at a common SNP is indistinguishable from contamination without population data, which the blind setting removes.
- **Contamination check.** Before the pilot, estimate HG002 60x cross-sample contamination from the reference fraction at 2,000 v4.2.1 hom-alt sites of common SNPs; report it, and if above 1% move the VAF floor of 3.2 to 0.08.
- **No filter on artifact signatures** (recorded, reported, not selected on).
- **Why hard.** Counts match the positives; remaining cues (strand locking, low base quality, read ends, linked mismatches, motif) are partly in `B` and partly read-level.

### 3.3 `tr_indel_pos` (n = 15: 10 het, 5 hom)
- **Truth.** v5.0q PASS inside its BED; v4.2.1 must agree where it covers.
- **Selection.** Indel of 1-4 bp in a homopolymer ≥ 10 bp or di/tri STR ≥ 20 bp; biallelic (1/2 excluded); isolation over the whole run ±30 bp. Difficulty as a data property on spanning reads: ±1-unit non-truth stutter ≥ 8%, or truth-allele fraction ≤ 0.40 (het) / ≤ 0.85 (hom). HiFi length check: this is what protects against v5.0q genotype errors, which concentrate exactly where the reads disagree with the label.
- **Pool.** On chr20, v5.0q-only PASS in-BED records in homopolymers ≥ 12 bp: 4,481, of which 785 are 1/2 → 3,696 biallelic (2,440 het, 1,256 hom). The post-isolation, post-difficulty, post-HiFi pool is reported before the pilot.
- **Why hard.** Stutter spreads reads over several lengths; round 3's "homopolymer insertion read as stutter".

### 3.4 `tr_stutter_neg` (n = 15)
- **Truth.** Inside the v5.0q BED (and v4.2.1 where it covers), no truth record within the run ±30 bp; HiFi ±1 fraction ≤ half the Illumina fraction; the same ±1 allele at ≥ 5% in at least one parent. Label `no`/`artifact`, `artifact_type: homopolymer_indel` (secondary).
- **Selection.** Homopolymer ≥ 10 bp or STR ≥ 20 bp with a ±1-unit allele on ≥ 4 spanning reads and ≥ 10% of them; run length distribution matched to 3.3. **POI placement:** the left-normalised anchor of the ±1 allele (first deleted base or insertion anchor), the same rule as 3.3; lint checks that the POI-offset-within-run distributions of 3.3 and 3.4 do not differ (KS p > 0.2), else redraw by seed (data-only).

### 3.5 `segdup_cmrg_pos` (n = 15)
- **Truth.** CMRG v1.00 (curated hifiasm assembly). SNV or indel ≤ 10 bp inside `GRCh38_segdups` or `lowmappabilityall`, **outside the v4.2.1 BED** (2,401 records, 1,631 SNVs); v5.0q agrees where it covers; HiFi alt.
- **Selection.** Depth ≥ 10; recoverability; no other CMRG record within 30 bp; half meet a difficulty cut (het VAF < 0.30 or ≥ 30% of reads at MAPQ < 20).
- **Why hard.** Paralog reads lower MAPQ and skew VAF; "low MAPQ means artifact" is wrong here.

### 3.6 `segdup_psv_neg` (n = 15)
- **Truth.** Inside the CMRG BED, segdup/low-mappability **and** the v5.0q BED; no CMRG or v5.0q record within 50 bp; ≥ 50 bp inside both BEDs; **outside the v4.2.1 BED** (same rule as 3.5); HiFi clean; trio recurrence; excluded if in an HG002 variable-copy region (GIAB HG002 CNV/SV v0.6 and the CMRG SV BED). The two assembly truths share a method, so HiFi + trio is what makes the negative label credible.
- **Matching to 3.5** on MAPQ < 20 fraction bin (< 10%, 10-30%, ≥ 30%) and depth bin, so MAPQ is not a free cue in either direction.
- **Source.** Primary: short-read caller false positives inside the CMRG BED, i.e. DRAGEN `HG002_300x.hard-filtered.vcf.gz` minus CMRG truth (what a user would actually review). Supplement: a pileup scan of the territory, which is **CMRG ∩ v5.0q ∩ segdup/low-mappability = 1.56 Mb** (about 800 windows of 2 kb; remote fetch 0.4 s each). Selection signal: one non-reference base at ≥ 3 reads and ≥ 10%.
- **Fallbacks, pre-registered in order, balance-preserving:** (1) VAF floor 7%; (2) extend to v5.0q ∩ segdup/low-mappability outside CMRG with the same HiFi/trio checks; (3) reduce strata 5 **and** 6 equally, to no fewer than 10 each. Under the per-pair balanced primary the pair weight does not change.
- **Why hard.** A PSV looks like a clean het; the tell is a band of reads sharing several mismatches, or discordant/low-MQ mates, which `B` now partly shows (linkage table, mate fields).

### 3.7 `indel_adjacent_misalign_neg` (n = 15)
- **Truth.** A truth indel ≥ 3 bp in both v4.2.1 and v5.0q. Distance is measured from the indel's **full reference span** (as `near_truth` in `real.py`): gap > 2 bp and ≤ 12 bp; the POI lies outside the repeat run containing the indel; no other truth record within 3 bp; inside both BEDs; HiFi clean. An automated check fails any item where the prompt's anchoring sentence could map the POI onto a truth record.
- **Selection.** At the POI, one non-reference base on ≥ 4 reads and ≥ 10%.
- **Why hard.** Good base quality, MAPQ 60, both strands; the cue is that the same reads end or soft-clip next to the indel, which V5 (soft clips) and `B`'s linkage table and 5′/3′ distances now expose.

### 3.8 `indel_adjacent_snv_pos` (n = 15)
- **Truth.** Truth SNV (0/1 or 1/1) 2-12 bp from the full span of a truth indel ≥ 3 bp, same span/run rule as 3.7, in v4.2.1 and v5.0q, agreeing; HiFi alt. Pool: `complexindel10bp`, `snpswithin10bp` (chr20 has 635 SNVs in `complexindel10bp`). ≥ 5 alt reads.
- **Purpose.** The mirror of 3.7, so "mismatch near indel = artifact" is not a free heuristic.

### 3.9 `anchor_controls` (n = 8, excluded from the primary)
2 clean het SNVs, 3 clean hom SNVs ("hom SNV read as reference" was a round-3 V failure), 3 clean negatives. Harm check, descriptive.

## 4. Conditions

| code | files | role |
|---|---|---|
| **B** = `Traw-r4` | `traw.md` | Baseline. The round-3 generator plus three additions any command-line user gets from `samtools view`/bam-readcount: (a) per-allele mate fields (proper-pair fraction, mate-unmapped count, "mate elsewhere" count, mate-MAPQ histogram); (b) distance to the 5′ end as well as the 3′ end; (c) an always-printed **linkage table**: co-occurrence of the POI allele with every non-reference column (mismatch, indel, soft-clip boundary) within ±15 bp. Generator hash frozen before the pilot. |
| **B+J** = `Traw-r4+Tcur` | `traw.md`, `tcur.md` | **C1.** BAMCP's text tools at shipped defaults: `get_region_summary` over the 200 bp window and `get_variant_curation_summary` always run at the POI's top non-reference allele via `probe_site` (J5), with risk score and reasons, `filter_check`, `applied_filters`, J3, J4, J6. |
| **B+V** = `V+Traw-r4` | `traw.md`, `view_wide.png`, `view_zoom.png` | **C2.** Viewer at the frozen round-4 commit (V1, V2, V5, V6). |

**Dropped:** `V` alone and `Tjson` alone (answered in rounds 1-3); the exploratory B+J+V arm (budget moved to items). A samtools-mpileup "unaided user" arm stays future work; it is not a second baseline.

**Why B is strengthened rather than kept identical to round 3.** Strata 5-7 depend on mates and read-level linkage. Withholding them from text would make text lose on information no real user lacks, the artificial handicap the user asked to avoid. Any C2 gain on 6-7 then measures perception and integration, not access. Comparability with round 3 is already broken by the caption.

**Common to all arms (frozen before the pilot):**
- **Caption sample line:** `sequenced human sample (Illumina short reads, ~60x whole-genome); the sample may carry germline variants or variants present in only a fraction of cells (mosaic); coordinates are re-indexed`. Its label risk on strata 4, 6, 7 (a low-VAF "real" answer could be biologically right) is closed by the HiFi and trio checks in 3.0.
- **Indel anchoring sentence:** "An insertion or deletion may be anchored anywhere within the repeat run that contains the position of interest, or within 2 bp of it." Strata 7-8 are constructed so this rule never maps the POI onto a truth record.
- Answer schema, invalid-answer rule and harness unchanged; `VIEWER_KEY` extended only to describe V1, V2, V5, V6.

**Models and repeats.** haiku, sonnet, opus (exact IDs in `build.json`, identical between pilot and main, else the calibration is void); one run per (item, condition, model); one fresh blind subagent per run. Two-model alternatives were simulated ([`r4/power_r4b_out.txt`](r4/power_r4b_out.txt): 176 items × 2 models gives power 0.88 vs 0.90 at uniform +0.10 and 0.12 vs 0.13 under the C1 prediction) and give no gain while raising the baseline, so all three models stay.

## 5. Endpoints

- **Primary endpoint.** Per run, strict joint correctness (`variant_present` AND `artifact_verdict`; `uncertain` = wrong). Aggregated as the **mean over the 4 pos/neg pairs of balanced accuracy** (per pair: mean of positive-stratum and negative-stratum accuracy), equal-weighted over the 3 models, **anchors excluded**. This estimand does not depend on the quotas and cannot be moved by a calling-threshold shift that trades sensitivity for specificity.
- **Guardrail (co-primary).** A contrast is "supported" only if, in addition to the test, neither mean sensitivity (over strata 1, 3, 5, 8) nor mean specificity (over 2, 4, 6, 7) falls by more than 0.05 (point estimate).
- **Key secondary: utility** (+1 correct, 0 `uncertain`, −1 wrong), with the same aggregation and test. This tells an abstention-driven null or harm (an arm that makes models hedge honestly) apart from a real loss. Sensitivity analysis at −2.
- `zygosity` het vs `low_vaf` does not affect the primary endpoint.
- **Descriptive secondaries:**
  1. sensitivity and specificity per arm;
  2. per-pair balanced accuracy;
  3. per-stratum and per-model deltas;
  4. genotype (joint + zygosity + allele, left-normalised with repeat-shift equivalence; `artifact_type` on stratum 4);
  5. calibration: Brier, 5-bin ECE, AUROC of confidence, confident-wrong rate, selective accuracy, AURC, abstention rate;
  6. accuracy by truth source;
  7. pooled item-weighted joint accuracy (the rounds 1-3 metric) for continuity;
  8. feature sufficiency: the G8 classifier (section 7) applied to main items;
  9. the "information headroom" gap: G8 AUROC on J/V features minus AUROC on `B`-visible features (reported, not gated; see Rejected critiques).

## 6. Exclusions and deviations

- **Invalid answers** count as wrong; the rate is reported by arm.
- **Infrastructure deviations** (did not `Read` every listed file, or read an unlisted one): rerun once with a fresh subagent, capped at 5% of runs; sensitivity analysis excludes the cell.
- **Adjudication before any main run, two raters.** Rater 1: the engineer, blind to labels, from `traw.md` plus PNGs, with main items shuffled together with 20 decoys (pilot-partition items of known label). Rater 2 (mandatory): a second human who has not seen the strata definitions, same protocol; if none is available, an automated proxy (G8 pooled posterior between 0.4 and 0.6 = "undecidable"). κ between raters is reported. Flagged items stay in the primary and are excluded in a sensitivity analysis. Flags are frozen before any model run.
- **Label errors found after the runs:** amendment; the primary keeps frozen labels, a sensitivity analysis uses corrected ones.
- Exclusion counts (recoverability, concordance, HiFi, trio, gnomAD, CNV, spacing) are reported per stratum.

## 7. Difficulty calibration: pilot and frozen rule

**Principle.** Freeze the selection rule and the **quotas** before the pilot; the pilot may move only within-stratum knobs; main items are drawn fresh afterwards with a new seed; nothing is ever selected on outcomes. With the per-pair balanced primary, the estimand weights are fixed by design (4 equal pairs), so pilot-driven changes affect precision and difficulty, not the target.

**Before any pilot run (firewall).**
- Freeze and hash the **eval-harness commit**: `Traw-r4` generator, fixture writer (V6 fields), strata predicates, knob grid, HiFi/trio rules, `calibrate.py`.
- `calibrate.py` ships with a **simulated operating characteristic**: under a prior on stratum B (each stratum's expected value ± 0.12), how often each branch (GO, knob move, re-pilot, NO-GO) fires and the resulting distribution of main-set primary B. Committed with its output.
- Leak-lint: 0 problems. Label audit of every pilot item against truth, normalisation, displayed reads and the HiFi/trio pileups.
- G8 computed (below).

**Pilot set.** 84 items from the pilot partition: per hard stratum 6 at the default knob setting and 4 at the pre-registered alternative setting (both levels are observed, so a knob move never extrapolates), plus 4 anchors.

**Pilot runs: B only.** 84 items × 3 models = **252 runs**, plus a **caption check** of the 18 round-3 real items under the new caption, B only × 3 models = **54 runs** (measures the caption's prior shift against round 3's 0.99). **No treatment arm is run in the pilot**: there is nothing to tune J or V against, and no discordance gate (main n is fixed in advance, see section 9).

**Knob grid (default → alternative piloted).**

| stratum | knob | default → alt |
|---|---|---|
| 1 | alt-count band | 4-6/≥7 halves → ≥7 only |
| 2 | VAF band | 0.05-0.20 → 0.08-0.20 |
| 3, 4 | homopolymer run floor | 10 → 12 bp |
| 3 | difficulty cut | stutter ≥ 8% → ≥ 12% |
| 5 | share of "difficult" items | ½ → ¾ |
| 6 | VAF floor | 10% → 7% (also pool fallback 1) |
| 7 | VAF floor | 10% → 15% |
| 8 | alt reads | ≥ 5 → 4-8 band |

**Estimation.** Beta-binomial empirical-Bayes shrinkage of each stratum × knob-level accuracy toward the pooled mean, with the within-item correlation of the 3 model runs modelled (effective n per 6-item cell about 8-9), then the predicted main primary B for the frozen quotas.

**Recoverability gate G8 (model-free, no forking path).** For each pair, on the **whole labelled candidate pool** of the pilot partition (all candidates passing the stratum predicates, not only items: hundreds per pair at zero model cost; for stratum 1, all about 10 pilot-partition mosaic sites plus about 60 main-partition sites, which are labelled sites, not drawn items), fit a 5-fold cross-validated logistic regression on read-level features (J3/J4 fields, linked mismatches, read-end stacking, soft-clip co-occurrence, mate fields). Requirements:
- pair AUROC ≥ 0.80, else the pair is redesigned (alternative knob) or, if still failing, replaced by its pre-registered reserve pair (`low_mapq_real_vs_mismap` from CMRG, built with the same rules);
- within stratum 1, any alt-count band with band-level AUROC < 0.70 is dropped from the sampling frame (this targets the "4 clean alt reads out of 60" undecidable case).

**Go/no-go.**

| criterion | GO | if failed |
|---|---|---|
| G1 predicted main primary B (3-model mean) | **0.62-0.78** | apply the knob rule once using piloted levels; re-predict; if still outside, one **B-only re-pilot** of 20 items × 3 models (60 runs) at the new settings; if still outside, **NO-GO** |
| G2 stratum knob moves | only for a stratum whose shrunken 90% interval lies wholly outside 0.45-0.88, or (to fix G1) the strata farthest from target in the needed direction | moves only to the piloted alternative level |
| G3 opus predicted primary B | ≤ 0.90 | treated under G1/G2 |
| G4 label audit | 0 confirmed label errors; ≤ 1 representation problem per stratum | ≥ 2 problems: review the full stratum with HiFi/trio; replace confirmed errors by seed order |
| G5 pools | every main stratum pool ≥ 2× quota after all exclusions | stratum 6 fallbacks (3.6), stratum 2 fallback (3.2); otherwise reduce the pair equally, never below 10 per stratum |
| G6 harness | ≥ 95% valid answers; 0 leak-lint problems; every `tcur.md` (generated, not run) has a curation block at the POI and parses; every PNG renders | fix the harness and re-run affected **B** pilot cells; G6 inspects file and parse properties only, never answer content |
| G8 recoverability | as above | as above |

**After GO.** Freeze the BAMCP product commit (J1-J6, V1-V2, V5-V6; thresholds set from round-3 items and pilot-partition candidate-pool features only, never from any model output of round 4). Draw main items (new seed), label-error audit, leak-lint, adjudication, signature report; hash the manifest; commit the pre-registration **before any main run**. Accept whatever baseline is realised; never redraw.

## 8. Confirmatory contrasts

Holm over two contrasts, two-sided α = 0.05, estimand = the primary endpoint (section 5):
- **C1.** B+J vs B.
- **C2.** B+V vs B.

**Test.** A cluster-level paired sign-flip permutation test (10,000 flips) on each item's model-averaged difference weighted by its balanced-accuracy weight (1 / (2 · n_class · 4)), with flips applied to whole clusters (section 3.0). The 95% CI is obtained by inverting the same test. Supported = Holm-adjusted p < 0.05, CI excludes 0, and the guardrail holds. The design effect (cluster vs item level) is reported. Sensitivity analyses: exact McNemar per model; a logistic GLMM (item and cluster random intercepts, model fixed).

**Equivalence.** Pre-registered TOST with margin ±0.07: a null may be read as "no meaningful gain" only if TOST succeeds. At n = 120 pair items TOST succeeds with probability 0.40 when Δ = 0 ([`r4/power_r4b_out.txt`](r4/power_r4b_out.txt)); ±0.05 would be essentially unreachable, so ±0.07 is the smallest margin with useful operating characteristics. When TOST fails the result is reported as inconclusive.

**Pre-specified directions (secondary; one summary statistic).** The only directional test is a sign test of how many of the 16 stratum × contrast cells marked + or − below are matched in sign (0/+ and 0/− cells are not counted; 0 cells are not counted). Stratum-level claims are forbidden.

| stratum | C1 (B+J) | C2 (B+V) |
|---|---|---|
| 1 mosaic pos | **0/−** (region summary at default `min_vaf = 0.10`, `min_alt_count = 3` shows 0 candidates for 4-6 alt reads; J6 framing mitigates) | +/0 (V2 fraction bar) |
| 2 mosaic caller-FP neg | + (J3 mismatch burden, J4 local noise) | + (scattered vs strand-locked) |
| 3 TR indel pos | 0 (J2 ships, so the context-only flag is gone) | 0 |
| 4 TR stutter neg | + | 0 |
| 5 segdup pos | − (low-MAPQ risk) unless J3 offsets | 0/+ (V6 mate shading) |
| 6 PSV neg | +/0 | + (co-occurring mismatch band) |
| 7 misalign neg | + (region summary lists the sibling indel) | + (V5 clips stacked at the indel) |
| 8 indel-adjacent pos | 0/− | 0/− |

**Exploratory:** model × condition interaction ("helps small models").

## 9. Power ([`r4/power_r4b.py`](r4/power_r4b.py), 2,000 simulations, 1,000 flips each)

**Model.** Logistic: stratum intercepts set to the section-3 expected B; model offsets −0.9 / 0 / +0.9; item SD 1.2; cluster SD 0.6 in clusters of 2; item × condition SD 1.0; nested arms (ρ = 0.5, shared run randomness q = 0.3); balanced-accuracy weights; one-sided count at the Holm rank-1 α of 0.025.

| scenario (true Δ on primary) | 120 pair items × 3 models | 136 × 3 | 160 × 3 | 176 × 2 |
|---|---|---|---|---|
| null (TOST ±0.07 success) | 0.01 (0.40) | 0.01 (0.50) | 0.01 (0.62) | 0.01 (0.52) |
| uniform +0.10 | **0.90** | 0.94 | 0.97 | 0.88 |
| concentrated +0.15 in 2 pairs (Δ = +0.075) | 0.60 | 0.66 | 0.77 | 0.55 |
| C1 as predicted in section 8 (Δ = +0.037) | **0.13** | 0.15 | 0.18 | 0.12 |
| C2 as predicted in section 8 (Δ = +0.043) | **0.20** | 0.23 | 0.26 | 0.20 |

The SE of the primary Δ is about 0.028 at 120 pair items.

**What this means, stated in the pre-registration.**
- The design detects a broad +0.10 gain (power 0.90) and has moderate power for a +0.15 gain concentrated in half the pairs.
- **The effects the design itself predicts are mixed in sign and pool to about +0.04, where power is 0.13-0.20 at any feasible n.** Even 160 pair items would raise this only to 0.18-0.26, and the pools (about 55-60 main-partition mosaic sites; PSV pool unknown) cap n. A null at the predicted effect sizes is therefore the most likely single outcome. It will be reported as inconclusive unless TOST succeeds, with per-pair balanced accuracies and the G8 sufficiency results showing *where* information existed but was not used.
- n = 120 pair items is chosen over 136-160: the power gain is small (0.90 → 0.94-0.97 at +0.10) against tighter pools and +144-432 runs.

## 10. BAMCP changes needed before the freeze

All through normal BAMCP PRs to `main` (each passing `make verify`: lint, mypy, 85% coverage, `make eval-smoke`, frontend gate, `make render-viewer` checked by eye for viewer changes). Merging deploys to prod, which is acceptable because these are product improvements. **Product thresholds are set from round-3 items and pilot-partition candidate-pool features only, and the product commit is frozen before the main draw; no round-4 model output is ever seen before that freeze, because the pilot runs B only.**

**Required (C1 tests the harness rather than the product without these):**
- **J5 (eval harness, `tcur_text`).** The curation summary is always produced at the POI's top non-reference allele with `probe_site` semantics (`min_baseq = 0`, no VAF or alt floor), reporting `filter_check`. Today `tfull_text` prints "not called" on 7 of 8 round-3 noisy negatives.
- **J1 (consistency bug).** One strand-bias rule (`>=` the threshold) in `is_balanced`, risks and recommendations, measured relative to the ref reads' strand split (Fisher/SOR-style), in `analysis/evidence.py` and `analysis/curation.py`. Round 3 printed "IMBALANCED" next to "Risk Score 0.00" at a 1F/9R site.
- **J2 (evidence-based homopolymer risk) — hard gate.** Replace the context-only risk (`ARTIFACT_HOMOPOLYMER_LENGTH_THRESHOLD = 4` fires on any run ≥ 4 bp; it flagged 6 of 10 true round-3 indels) with a spanning-read rule: raise the risk only when the alt fraction is < 0.2 or within about 2× the ±1 stutter fraction; run length stays as information. If J2 is not merged, round 4 does not start.
- **J3 (per-allele comparison).** For alt vs ref reads: mismatch burden excluding the site (per 100 aligned bp, computed against the reference), soft-clip fraction, median MAPQ, **median mate MAPQ and proper-pair fraction**, near-end fraction vs expectation, mean base quality. In `analysis/curation.py`, using the read loop of `core/tools.probe_site`.
- **J4 (local context), computed over ±300 bp** (not the 2 kb window: fixtures keep only reads within ±450 bp of the site, so depth falls to zero about 600 bp out and a 2 kb median is a fixture artifact). Depth relative to the ±300 bp median, other non-reference candidates within 100 bp with their strand and base-quality profile, repeat unit and run length.
- **J6 (neutral low-VAF framing).** When the probed site fails only the VAF or alt-count floor, `filter_check` states the observed VAF and that it lies "below the default calling floor (0.10); fractions of 0.02-0.10 are consistent with either low-level mosaicism or noise", instead of a bare failure. `get_region_summary` runs at shipped defaults (a real user would), so the "0 candidates" line remains, which is why stratum 1 is still predicted 0/−.

**Required (viewer, `src/bamcp/static/`, `npm run build`, never hand-editing `dist/`):**
- **V1.** Neutral POI marker (ruler caret + outlined column) via a new blind-safe setting separate from `activeVariantPosition`.
- **V2.** IGV-style allele-fraction coverage bar (columns ≥ 5% non-reference split by base colour; insertion counts marked).
- **V5.** Soft clips drawn in blind mode (stratum 7's predicted C2 mechanism depends on them).
- **V6.** Paired-end display: reads shaded by mate MAPQ, discordant / mate-elsewhere reads coloured; needs the fixture writer to keep pair flags and `MQ` (section 3.0).

**Eval harness (not product):** `Traw-r4` additions (section 4), fixture writer V6 fields, HiFi/trio label checks.

**Deferred (limitations):** V3 muted reference tiles; V4 insertion-length labels.

**Leak checks for new text.** `tcur.md` and `traw.md` must not contain truth-set or stratification names (GIAB, CMRG, v4.2.1, v5.0q, benchmark, segdup, mosaic-benchmark, HiFi, PacBio, trio); "mosaic" only in `caption.md` and `prompt.md`. External lookups off (no ClinVar or gnomAD). The "mate elsewhere" sentinel contig must carry no real coordinate.

## 11. Implementation plan (`src/bamcp/eval/visual/`)

| module | change |
|---|---|
| `truth.py` (new) | `TruthSet` (VCF + BED, FILTER-aware), loaders for v4.2.1, v5.0q, CMRG, Mosaic v1.1 and the mosaic release's candidate/validator files; `covers`, `records_near`, class-aware `concordant`, `stratum_tags`; `orthogonal_check(site, label)` over the HiFi BAM; `trio_recurrence(site)` over the HG003/HG004 60x BAMs; `cluster_id(site)` from `genomicSuperDups` components. |
| `strata_r4.py` (new) | Frozen `QUOTAS_R4`; one predicate per stratum; `KNOBS` (default + piloted alternative); `partition(site)` by cluster hash; spacing (≥ 10 kb, ≤ 2 per cluster); `label_for`; recoverability; matching (`match_negatives` on alt, depth, segdup, MAPQ bins). |
| `real.py` | `_allele_features` (strand, base quality, MAPQ < 20 fraction, mate fields, spanning-read indel-length histogram); generators `iter_mosaic_sites`, `iter_mosaic_caller_fps`, `iter_cmrg_windows`, `iter_dragen_fps`, `iter_indel_adjacent` (full-span distance, run rule), `iter_tr_runs`; JSONL scan cache; `write_real_fixtures` keeps pair flags/`MQ` with re-indexed or sentinel mates; `SAMPLE_LINE_R4`. |
| `stimuli.py` | `traw_r4_text()` (round-3 panel + mate fields + 5′ distance + linkage table); `tcur_text()` (J5); conditions `"Traw-r4"`, `"Traw-r4+Tcur"`, `"V+Traw-r4"`; `FILE_DESCRIPTIONS`; `prompt_text(..., repeat_anchor=True)`; `VIEWER_KEY` for V1, V2, V5, V6. |
| `leaklint.py` | Vocabulary checks for `tcur.md` and the new `traw.md` fields; anchoring check for strata 7-8; POI-offset lint for 3/4. |
| `grade.py` | r4 `PRIMARY_CONTRASTS`; per-pair balanced accuracy; guardrail; utility; `allele_equivalent()`; accuracy by truth source. |
| `stats.py` | `cluster_signflip_test(weighted_diffs, clusters)`, inverted CI, `holm`, `tost`, `aurc`, 5-bin ECE. |
| `sufficiency.py` (new) | G8: pool feature extraction and cross-validated logistic AUROC per pair and band; posterior for the adjudication proxy. |
| `calibrate.py` (new) | Reads only the pilot B results; shrinkage with within-item correlation; G1-G8; frozen knob rule; writes `calibration.json`; `--simulate` prints the operating characteristic. |
| `cli.py` | `select-r4 --split pilot|main --seed N`, `sufficiency`, `calibrate`, `build-real --conditions …`; `build.json` records model IDs, generator hashes (`Traw-r4`, `tcur`), viewer hash, harness and product commits, manifest sha. |
| `selftest.py` | Oracle = 1.0, random ≈ chance, always-real / always-artifact → per-pair balanced accuracy 0.5 exactly. |
| tests | `truth.py` concordance/equivalence/FILTER handling, predicates on tiny synthetic BAMs, cluster spacing, `cluster_signflip_test`, `tost`, `calibrate.py` on canned inputs, J1-J4/J6 in curation, V6 fixture fields. |

**Schedule (one engineer coordinating subagents):**
- **Days 1-2.** `truth.py`, `strata_r4.py`, generators, `Traw-r4`, fixture V6 fields. Background: mosaic caller-FP extraction, DRAGEN FP extraction in CMRG, the 1.56 Mb PSV scan, TR runs, matched pools, HiFi and trio checks (about 3 s per site; several thousand candidates ≈ 3-6 h), contamination estimate. In parallel open the product PRs (J1-J6, V1, V2, V5, V6).
- **Day 3.** Freeze the harness commit; G8; pilot draw, label audit, leak-lint; the 306 B pilot runs; `calibrate`; GO/NO-GO (re-pilot if needed).
- **Day 4.** Merge product PRs with CI green; freeze the product commit; main draw; label-error audit; adjudication; signature report; commit the pre-registration.
- **Days 5-6.** 1,152 main runs.
- **Day 7.** Grade, analyse, sensitivity analyses; append the round-4 section to `docs/eval/visual-results.md`.

## 12. Freeze list

Before the pilot: caption, prompt and anchoring sentence; `FILE_DESCRIPTIONS`; answer schema; `Traw-r4` generator hash; fixture writer; strata predicates, quotas, knob grid; HiFi/trio/gnomAD/CNV/spacing rules; seeds; `calibrate.py` and its simulated operating characteristic; `sufficiency.py`; `grade.py`, `stats.py`, analysis scripts; model IDs.
Before the main runs (with the pre-registration): BAMCP product commit, `viewer_source_hash`, `tcur` generator hash, J thresholds in `constants.py`, `VIEWER_KEY`, `calibration.json`, G8 output, manifest sha256, label-audit results, adjudication flags, signature report.
Any later change is a dated amendment.

## 13. Compute budget

| phase | runs |
|---|---|
| Pilot: B × 3 models × 84 items | 252 |
| Caption check: B × 3 models × 18 round-3 real items | 54 |
| Conditional B-only re-pilot (G1) | ≤ 60 |
| Pilot reruns (reserve) | ≤ 20 |
| **Pilot total** | **≤ 386** |
| Main: {B, B+J, B+V} × 3 models × 128 items | 1,152 |
| Main reruns (5% cap) | ≤ 58 |
| **Main total** | **≤ 1,210** |
| **Overall** | **≤ 1,596** (draft: ≤ 1,258) |

Background CPU (no model runs): about 10-16 h of remote BAM work (Illumina scans plus HiFi and trio checks).

**Deferred to round 5, with reasons:**
- **Natural low depth** (6-14x, GC-extreme): mainly a calibration question that needs leniency for `uncertain`.
- **Long indels of 15-50 bp shown only as soft clips:** representation mismatch risk; BAMCP's indel size cap interferes.
- **2-3 kb window triage** ("what is in this region?"): a different task with per-window scoring; the honest test of `scan_variants` and `get_region_summary` against a text dump.
- **Agentic arm** (the model calls BAMCP tools): different estimand; image tool results in the restricted harness unsolved.
- **Multi-allelic 1/2 TR genotypes:** mainly a genotype-endpoint issue.

## 14. Main risks and mitigations

1. **Truth error concentrated in the hard regions selected for difficulty.** Mitigations: HiFi check on every item, trio recurrence on strata 4 and 6, CNV exclusion, class-aware concordance, the residual label-error audit, accuracy by truth source.
2. **Irreducible items.** Mitigations: G8 (pool AUROC ≥ 0.80 per pair, band drop < 0.70), two-rater adjudication with κ, the utility secondary.
3. **Baseline misses 0.70.** The pre-pilot prediction interval is ±0.14; two knob levels piloted per stratum; G1 window 0.62-0.78; one B-only re-pilot; NO-GO otherwise.
4. **Mixed-sign feature effects.** Predicted and pre-registered; power at the predicted effects is low (section 9), stated up front; TOST distinguishes "no meaningful gain" from "inconclusive".
5. **BAMCP features harm** (low-MAPQ flag on stratum 5, region-summary filter on stratum 1): predicted, therefore product findings.
6. **Tuning the product to the test.** Product thresholds come only from round-3 items and pilot-partition pool features; no round-4 model output exists before the product freeze.
7. **Caption prior shift.** Measured by the 54-run caption check plus anchors and the sensitivity/specificity split.
8. **Small pools** (about 72 eligible mosaic sites; PSV pool unknown until the 1.56 Mb scan plus DRAGEN FPs). Pre-registered balance-preserving fallbacks; the minimum per stratum is 10.
9. **Non-independence** (paralogs, shared reads). Cluster partition, spacing, cluster-level test, design effect reported.
10. **Model drift.** IDs pinned; a mismatch voids the calibration.

## 15. Rejected or modified critiques

- **C1-#2, second half (gate on J/V-feature AUROC exceeding `B`-visible AUROC by ≥ 0.05): rejected as a gate, kept as a report.** The headroom we want to measure includes model *uptake*: a baseline at 0.70 whose `B` text already contains separable information still leaves room for a pre-digested J summary or an image to help integration. Gating on information asymmetry would also conflict with strengthening `B` (section 4), which deliberately removes the asymmetry.
- **C2-#6 (AUROC ≥ 0.85 per stratum): modified.** AUROC is undefined within a single-class stratum, so it is applied per pair at 0.80 (critic 1's level), with the band-level < 0.70 drop inside stratum 1 targeting the undecidable low-alt case critic 2 identified.
- **C1-#1 ("≥ 2 confirmed errors in a stratum → drop the pair"): modified.** Errors trigger a full-stratum HiFi/trio review and replacement by seed order. This uses labels only, before any model run, so it is not outcome selection, and it avoids losing a pair to a fixable labelling problem.
- **C1-#7, "confirmatory on sonnet and opus only, haiku exploratory": rejected.** The simulation shows no power gain (0.88 vs 0.90 at +0.10; 0.12 vs 0.13 at the C1 prediction), and dropping haiku raises the baseline away from 0.70 and loses the small-model question.
- **C1-#7, TOST margin ±0.05: modified to ±0.07.** At ±0.05 equivalence is essentially unreachable at feasible n; ±0.07 succeeds in 40% of null worlds at n = 120 and is stated as such.
- **C1-#3, abstention: partially rejected.** The primary stays strict (`uncertain` = wrong) so that "joint correctness" keeps its rounds 1-3 meaning and the 0.70 calibration target is well defined. Critic 1's remedy is adopted in its stated form: utility is the key secondary with the same test, so abstention effects are identified.
- **C2-#10 (separate `Traw+link` control arm on strata 6-7): rejected in favour of a stronger fix.** The linkage table, 5′ distance and mate fields are folded into `B` itself for every stratum, which removes the text-only handicap everywhere instead of measuring it in a side arm.
- **C2-#2, alternative "germline caption for those strata": not needed.** The HiFi and trio checks make the negative labels robust under the mosaic caption, and a per-stratum caption would leak the stratum.
- **C1-#9, VerifyBamID: replaced** by a cheaper reference-fraction contamination estimate at hom-alt common-SNP sites, which serves the same purpose.
- **C2-#15, sampling positives from one caller's output: modified.** Negatives in strata 2 and 6 are now caller false positives; positives are not restricted to caller-found sites (that would select easier positives). The claim scope is restated as "candidate sites with a pileup signal (caller-flagged where noted)".
- **C1-#5, option of dropping strata 5-7 to exploratory: not taken.** The preferred option (V5 and V6, with pair information in every arm including `B` and `Tcur`) is implemented.
- **C1-#19 (wrong section reference): fixed.** Deferred cases are now referenced as "end of section 13".

All other critique points (C1 #4, 6, 8, 10-18; C2 #1, 3-5, 7-9, 11-14, 16-22) are adopted as described in the sections above.
