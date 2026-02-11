# BAMCP Evaluation Harness: Architecture for LLM Genomic Reasoning Assessment

## The Core Idea

Every time an LLM interprets a variant through BAMCP, two things happen simultaneously:
1. The user gets an explanation (the product)
2. The system records whether the LLM got it right (the research)

The MCP server sits between the LLM and the genomic data. It controls what
evidence the LLM sees, captures the LLM's reasoning, and can compare that
reasoning against ground truth. This position is uniquely powerful for
evaluation — you're not scraping chat logs after the fact, you're instrumenting
the reasoning pipeline at the tool level.

---

## Architecture

```
┌────────────────────────────────────────────────────────────┐
│                    MCP Host (Claude, etc.)                  │
│                                                            │
│  ┌──────────────────────────────────────────────────────┐  │
│  │                       LLM                             │  │
│  │                                                       │  │
│  │  1. Receives variant evidence from tools              │  │
│  │  2. Reasons about pathogenicity                       │  │
│  │  3. Returns classification + reasoning                │  │
│  │  4. Gets corrective context if reasoning was flawed   │  │
│  └───────────┬───────────────────────────▲───────────────┘  │
│              │tool calls                 │tool results       │
│              │                           │                   │
│  ┌───────────▼───────────────────────────┴───────────────┐  │
│  │              MCP App Viewer (iframe)                    │  │
│  │  - Alignment visualization                             │  │
│  │  - Variant table with eval indicators                  │  │
│  │  - updateModelContext() on every interaction            │  │
│  │  - sendMessage() to trigger classification             │  │
│  └───────────┬───────────────────────────▲───────────────┘  │
│              │callServerTool             │results            │
└──────────────┼───────────────────────────┼──────────────────┘
               │                           │
    ┌──────────▼───────────────────────────┴──────────────┐
    │                  BAMCP MCP Server                     │
    │                                                      │
    │  ┌─────────────────────────────────────────────┐     │
    │  │            Instrumentation Layer              │     │
    │  │                                               │     │
    │  │  Every tool call logged:                      │     │
    │  │  - tool name, arguments, timestamp            │     │
    │  │  - full response returned                     │     │
    │  │  - session context (what region, what BAM)    │     │
    │  │  - ground truth comparison (when available)   │     │
    │  └──────────────────┬──────────────────────────┘     │
    │                     │                                 │
    │  ┌──────────────────▼──────────────────────────┐     │
    │  │           Tool Layer                          │     │
    │  │                                               │     │
    │  │  User-facing:                                 │     │
    │  │  - visualize_region (MCP App)                 │     │
    │  │  - classify_variant (structured reasoning)    │     │
    │  │  - lookup_clinvar / lookup_gnomad             │     │
    │  │                                               │     │
    │  │  App-only (viewer data):                      │     │
    │  │  - fetch_reads, fetch_coverage, etc.          │     │
    │  │                                               │     │
    │  │  Eval-specific:                               │     │
    │  │  - run_evaluation (batch mode)                │     │
    │  │  - get_eval_report (results summary)          │     │
    │  └──────────────────┬──────────────────────────┘     │
    │                     │                                 │
    │  ┌──────────────────▼──────────────────────────┐     │
    │  │        Evidence & Scaffold Layer              │     │
    │  │                                               │     │
    │  │  - Ground truth DB (ClinVar gold standard)    │     │
    │  │  - ACMG criteria templates                    │     │
    │  │  - Known failure mode catalog                 │     │
    │  │  - Corrective context generators              │     │
    │  └──────────────────┬──────────────────────────┘     │
    │                     │                                 │
    │  ┌──────────────────▼──────────────────────────┐     │
    │  │           Eval Store (SQLite)                 │     │
    │  │                                               │     │
    │  │  Tables:                                      │     │
    │  │  - eval_sessions                              │     │
    │  │  - tool_invocations                           │     │
    │  │  - classifications (LLM output vs ground truth│     │
    │  │  - failure_modes (categorized errors)         │     │
    │  │  - reasoning_chains (step-by-step logic)      │     │
    │  └─────────────────────────────────────────────┘     │
    │                                                      │
    │  ┌─────────────────────────────────────────────┐     │
    │  │           Data Layer (existing)               │     │
    │  │  - pysam BAM/CRAM parsing                     │     │
    │  │  - ClinVar API client                         │     │
    │  │  - gnomAD API client                          │     │
    │  │  - Reference genome                           │     │
    │  └─────────────────────────────────────────────┘     │
    └──────────────────────────────────────────────────────┘
```

---

## The Key Tool: `classify_variant`

This is the tool that does double duty as product AND evaluation. When the LLM
encounters a variant (either because a user clicked it, or because an eval suite
is running), it calls `classify_variant`. The tool:

1. Assembles the evidence package for the variant
2. Provides an ACMG reasoning scaffold
3. Tells the LLM to classify and show its work
4. Captures the structured output
5. Compares against ground truth (if available)
6. Returns the result to the user AND logs the evaluation data

### Tool Definition

```python
@server.tool(
    name="classify_variant",
    description="""Classify a genomic variant using ACMG/AMP guidelines.

    Provide your classification and reasoning in the structured format below.
    Apply each ACMG criterion explicitly. State the evidence for and against.
    Assign a final classification: Pathogenic, Likely Pathogenic, VUS,
    Likely Benign, or Benign.

    IMPORTANT: This is a research-grade analysis, not a clinical report.
    Always state your confidence level and what additional evidence would
    change your assessment."""
)
async def classify_variant(
    chrom: str,
    pos: int,
    ref: str,
    alt: str,
    context: dict | None = None  # optional: gene, transcript, BAM stats
) -> dict:
    """
    Returns structured evidence AND captures LLM classification for eval.
    """

    # 1. Gather all available evidence
    evidence = await assemble_evidence(chrom, pos, ref, alt)

    # 2. Check ground truth (don't return to LLM — only for scoring)
    ground_truth = await lookup_ground_truth(chrom, pos, ref, alt)

    # 3. Build the reasoning scaffold
    scaffold = build_acmg_scaffold(evidence)

    # 4. Log the invocation (pre-classification)
    invocation_id = await eval_store.log_invocation(
        tool="classify_variant",
        variant=f"{chrom}:{pos}{ref}>{alt}",
        evidence_provided=evidence,
        ground_truth=ground_truth,  # stored but not returned
        timestamp=datetime.utcnow(),
        session_id=current_session_id(),
    )

    # 5. Return evidence + scaffold to LLM
    # The scaffold shapes HOW the LLM reasons without giving away the answer
    return {
        "content": [
            {
                "type": "text",
                "text": format_evidence_report(evidence, scaffold)
            }
        ],
        "structuredContent": {
            "invocation_id": invocation_id,  # track this classification
            "evidence": evidence,
            "scaffold": scaffold,
            "has_ground_truth": ground_truth is not None,
        }
    }
```

### The Evidence Package

What `assemble_evidence` returns — this is what the LLM gets to reason with:

```python
async def assemble_evidence(chrom, pos, ref, alt) -> dict:
    return {
        # From the BAM (observed data)
        "observation": {
            "depth": 45,
            "alt_reads": 21,
            "ref_reads": 24,
            "vaf": 0.467,
            "strand_bias": {"forward": 10, "reverse": 11},  # alt reads
            "avg_mapping_quality": 58.2,
            "avg_base_quality": 34.7,
            "near_read_end_fraction": 0.05,  # artifact indicator
        },

        # From ClinVar
        "clinvar": {
            "variation_id": "17661",
            "clinical_significance": "Pathogenic",
            "review_status": "criteria_provided,_multiple_submitters",
            "stars": 3,
            "conditions": ["Li-Fraumeni syndrome", "Hereditary cancer"],
            "last_evaluated": "2024-08-15",
            "submitter_count": 12,
        },

        # From gnomAD
        "population_frequency": {
            "global_af": 0.00001,
            "max_pop_af": 0.00003,  # highest in any population
            "max_pop": "European (non-Finnish)",
            "homozygote_count": 0,
            "filter": "PASS",
        },

        # From gene annotation
        "gene_context": {
            "gene": "TP53",
            "transcript": "NM_000546.6",
            "consequence": "missense_variant",
            "protein_change": "p.Arg248Trp",
            "exon": "7/11",
            "domain": "DNA-binding domain",
            "constraint": {  # from gnomAD
                "loeuf": 0.19,  # very constrained
                "missense_z": 3.72,
            },
        },

        # Computational predictions
        "predictions": {
            "sift": "Deleterious (0.0)",
            "polyphen2": "Probably_damaging (1.0)",
            "cadd_phred": 35.0,
            "revel": 0.95,
        },
    }
```

### The ACMG Scaffold

This is the reasoning structure the tool provides. It doesn't give the answer —
it gives the framework for arriving at one:

```python
def build_acmg_scaffold(evidence: dict) -> dict:
    """
    Returns applicable ACMG criteria with the evidence needed to evaluate each.
    The LLM must apply each criterion and state whether it's met.
    """
    scaffold = {
        "instructions": (
            "Evaluate each applicable ACMG/AMP criterion below. "
            "For each, state: (1) the criterion code, (2) whether it is met, "
            "(3) the specific evidence supporting your decision, (4) your "
            "confidence (high/medium/low). Then provide a final classification."
        ),
        "criteria_to_evaluate": [],
        "classification_options": [
            "Pathogenic", "Likely Pathogenic", "Uncertain Significance (VUS)",
            "Likely Benign", "Benign"
        ],
        "response_format": {
            "criteria_applied": [
                {"code": "PVS1/PS1/PM2/etc", "met": True, "evidence": "...", "confidence": "high"}
            ],
            "final_classification": "...",
            "confidence": "high/medium/low",
            "key_uncertainties": ["..."],
            "additional_evidence_needed": ["..."],
        }
    }

    # Dynamically include relevant criteria based on variant type
    if evidence["gene_context"]["consequence"] == "missense_variant":
        scaffold["criteria_to_evaluate"].extend([
            {
                "code": "PS1",
                "description": "Same amino acid change as an established pathogenic variant",
                "evaluate_with": "clinvar data, protein_change field",
            },
            {
                "code": "PM1",
                "description": "Located in a mutational hot spot or well-established functional domain",
                "evaluate_with": "gene_context.domain",
            },
            {
                "code": "PM2",
                "description": "Absent from controls (or extremely low frequency)",
                "evaluate_with": "population_frequency — note: PM2_supporting if AF < 0.0001",
            },
            {
                "code": "PP3",
                "description": "Multiple computational evidence supports deleterious effect",
                "evaluate_with": "predictions (SIFT, PolyPhen2, CADD, REVEL)",
            },
            {
                "code": "BP4",
                "description": "Computational evidence suggests no impact (counter-evidence)",
                "evaluate_with": "predictions — apply if tools suggest benign",
            },
            {
                "code": "BA1",
                "description": "Allele frequency > 5% in any population (standalone benign)",
                "evaluate_with": "population_frequency.max_pop_af",
            },
        ])

    # Add observation quality criteria (custom — not ACMG but critical)
    scaffold["quality_assessment"] = {
        "instructions": (
            "Before classifying, assess the quality of the observed data. "
            "Flag any concerns about the reliability of the variant call itself."
        ),
        "checks": [
            {"check": "depth", "threshold": "≥20x for reliable het call", "value": evidence["observation"]["depth"]},
            {"check": "strand_bias", "threshold": "both strands represented", "value": evidence["observation"]["strand_bias"]},
            {"check": "mapping_quality", "threshold": "≥40 average", "value": evidence["observation"]["avg_mapping_quality"]},
            {"check": "base_quality", "threshold": "≥30 average", "value": evidence["observation"]["avg_base_quality"]},
            {"check": "read_end_fraction", "threshold": "<0.1 (artifact if high)", "value": evidence["observation"]["near_read_end_fraction"]},
        ]
    }

    return scaffold
```

### Why This Scaffold Matters for Evaluation

The scaffold does three things simultaneously:

1. **For the user**: It makes the LLM's reasoning transparent and structured,
   so even a layman can follow the logic chain.

2. **For evaluation**: It forces the LLM to show its work in a parseable format.
   You can extract each criterion application and score it independently.

3. **For improvement**: When you find the LLM consistently misapplies PM2
   (population frequency), you can add more specific guidance to that criterion's
   scaffold entry. The scaffold is your primary lever for improving LLM accuracy
   without fine-tuning.

---

## Capturing the LLM's Classification

Here's the critical piece: how do you get the LLM's structured classification
back into the server for scoring? Two approaches:

### Approach 1: Second Tool Call (Recommended)

After the LLM reasons through the evidence, it calls a `submit_classification`
tool with its structured result:

```python
@server.tool(
    name="submit_classification",
    description="""Submit your variant classification after reasoning through
    the evidence. You MUST call this after classify_variant to record your
    assessment. Provide your classification in the structured format."""
)
async def submit_classification(
    invocation_id: str,
    classification: str,           # Pathogenic / Likely Pathogenic / VUS / etc.
    confidence: str,               # high / medium / low
    criteria_applied: list[dict],  # [{code, met, evidence, confidence}, ...]
    key_uncertainties: list[str],
    additional_evidence_needed: list[str],
    reasoning_summary: str,
) -> dict:

    # 1. Retrieve the ground truth for this invocation
    invocation = await eval_store.get_invocation(invocation_id)
    ground_truth = invocation.ground_truth

    # 2. Score the classification
    score = score_classification(
        predicted=classification,
        ground_truth=ground_truth,
        criteria_applied=criteria_applied,
    )

    # 3. Detect failure modes
    failure_modes = detect_failure_modes(
        predicted=classification,
        ground_truth=ground_truth,
        criteria_applied=criteria_applied,
        evidence=invocation.evidence_provided,
    )

    # 4. Log everything
    await eval_store.log_classification(
        invocation_id=invocation_id,
        classification=classification,
        confidence=confidence,
        criteria_applied=criteria_applied,
        reasoning_summary=reasoning_summary,
        score=score,
        failure_modes=failure_modes,
    )

    # 5. Return acknowledgment (and corrective feedback if wrong)
    result = {
        "recorded": True,
        "invocation_id": invocation_id,
    }

    # In eval mode: return ground truth comparison
    # In user mode: return disclaimer only
    if is_eval_mode():
        result["ground_truth"] = ground_truth
        result["score"] = score
        result["failure_modes"] = failure_modes
        result["feedback"] = generate_corrective_feedback(
            classification, ground_truth, failure_modes
        )
    else:
        result["disclaimer"] = (
            "Classification recorded. This is a research-grade analysis. "
            "Clinical decisions require validated testing by a certified laboratory."
        )

    return {"content": [{"type": "text", "text": json.dumps(result, indent=2)}]}
```

### Approach 2: Model Context Capture (Supplementary)

The MCP App viewer can also capture classification data via `updateModelContext`
round-trips. When the user interacts with a variant in the viewer and the LLM
responds, the viewer can parse the response and push structured data back:

```typescript
// In the viewer UI — after LLM provides classification in chat
app.ontoolresult = (result) => {
  if (result.structuredContent?.invocation_id) {
    // Track that this variant has been classified
    currentInvocationId = result.structuredContent.invocation_id;
  }
};

// When user clicks "Accept classification" or moves to next variant
acceptButton.addEventListener('click', async () => {
  await app.updateModelContext({
    content: [{
      type: "text",
      text: `User reviewed classification for invocation ${currentInvocationId}. ` +
            `User feedback: accepted / rejected / modified. ` +
            `Now call submit_classification with your final assessment.`
    }]
  });
});
```

---

## Failure Mode Detection

This is where the research value lives. `detect_failure_modes` categorizes
*why* the LLM got something wrong, not just *that* it was wrong.

```python
class FailureMode(Enum):
    # Classification errors
    FALSE_PATHOGENIC = "false_pathogenic"          # Called pathogenic, actually benign/VUS
    FALSE_BENIGN = "false_benign"                  # Called benign, actually pathogenic
    OVERCONFIDENT_VUS = "overconfident_vus"         # Classified definitively, should be VUS
    UNDERCONFIDENT_KNOWN = "underconfident_known"   # Called VUS, actually has clear classification

    # Reasoning errors
    CRITERIA_MISAPPLICATION = "criteria_misapplication"  # Applied ACMG criterion incorrectly
    FREQUENCY_MISINTERPRETATION = "frequency_misinterpretation"  # Misread population data
    PREDICTION_OVERRELIANCE = "prediction_overreliance"  # Too much weight on computational
    EVIDENCE_IGNORED = "evidence_ignored"                # Had evidence, didn't use it
    EVIDENCE_FABRICATED = "evidence_fabricated"           # Cited evidence not in the data

    # Data quality errors
    ARTIFACT_MISSED = "artifact_missed"              # Didn't flag low-quality data
    ARTIFACT_OVERCALLED = "artifact_overcalled"       # Flagged good data as artifact

    # Communication errors
    MISLEADING_CONFIDENCE = "misleading_confidence"   # Stated high confidence incorrectly
    CLINICAL_LANGUAGE = "clinical_language"            # Used clinical diagnostic language


def detect_failure_modes(predicted, ground_truth, criteria_applied, evidence) -> list[dict]:
    """
    Analyzes the LLM's classification and reasoning to identify specific failure modes.
    Returns a list of detected failures with explanations.
    """
    failures = []

    if not ground_truth:
        return failures  # Can't evaluate without ground truth

    gt_class = ground_truth["classification"]
    gt_criteria = ground_truth.get("expected_criteria", {})

    # --- Classification-level failures ---

    severity_order = ["Benign", "Likely Benign", "Uncertain Significance",
                      "Likely Pathogenic", "Pathogenic"]

    pred_idx = severity_order.index(predicted) if predicted in severity_order else -1
    gt_idx = severity_order.index(gt_class) if gt_class in severity_order else -1

    if pred_idx > gt_idx + 1:  # Predicted more severe than truth by >1 step
        failures.append({
            "mode": FailureMode.FALSE_PATHOGENIC,
            "severity": "high",
            "detail": f"Classified as {predicted}, ground truth is {gt_class}",
            "direction": "overcall",
        })
    elif pred_idx < gt_idx - 1:  # Predicted less severe
        failures.append({
            "mode": FailureMode.FALSE_BENIGN,
            "severity": "high",
            "detail": f"Classified as {predicted}, ground truth is {gt_class}",
            "direction": "undercall",
        })

    # --- Criteria-level failures ---

    applied_codes = {c["code"]: c for c in criteria_applied}

    for code, expected in gt_criteria.items():
        if code not in applied_codes:
            failures.append({
                "mode": FailureMode.EVIDENCE_IGNORED,
                "severity": "medium",
                "detail": f"Did not evaluate criterion {code}, expected: {expected['met']}",
                "criterion": code,
            })
        elif applied_codes[code]["met"] != expected["met"]:
            failures.append({
                "mode": FailureMode.CRITERIA_MISAPPLICATION,
                "severity": "medium",
                "detail": (
                    f"Criterion {code}: LLM said {'met' if applied_codes[code]['met'] else 'not met'}, "
                    f"expected {'met' if expected['met'] else 'not met'}"
                ),
                "criterion": code,
                "llm_reasoning": applied_codes[code].get("evidence", ""),
                "expected_reasoning": expected.get("evidence", ""),
            })

    # --- Hallucination detection ---

    for criterion in criteria_applied:
        evidence_text = criterion.get("evidence", "").lower()
        # Check if LLM references data not in the evidence package
        if "clinvar" in evidence_text and not evidence.get("clinvar"):
            failures.append({
                "mode": FailureMode.EVIDENCE_FABRICATED,
                "severity": "critical",
                "detail": f"Referenced ClinVar data in {criterion['code']} but none was provided",
                "criterion": criterion["code"],
            })

    # --- Population frequency specific checks ---

    for criterion in criteria_applied:
        if criterion["code"] in ("PM2", "BA1", "BS1"):
            pop_data = evidence.get("population_frequency", {})
            if pop_data:
                af = pop_data.get("global_af", 0)
                # PM2: should be absent/extremely rare
                if criterion["code"] == "PM2" and criterion["met"] and af > 0.01:
                    failures.append({
                        "mode": FailureMode.FREQUENCY_MISINTERPRETATION,
                        "severity": "high",
                        "detail": f"Applied PM2 (absent from controls) but AF is {af}",
                        "criterion": "PM2",
                    })
                # BA1: should be > 5%
                if criterion["code"] == "BA1" and criterion["met"] and af < 0.05:
                    failures.append({
                        "mode": FailureMode.FREQUENCY_MISINTERPRETATION,
                        "severity": "high",
                        "detail": f"Applied BA1 (AF>5%) but AF is {af}",
                        "criterion": "BA1",
                    })

    return failures
```

---

## The Improvement Loop: Corrective Context

Here's how evaluation findings feed back into better responses WITHOUT
fine-tuning. The mechanism is the ACMG scaffold itself.

### Step 1: Aggregate Failure Patterns

```python
async def get_failure_summary() -> dict:
    """Analyze accumulated evaluation data for systematic patterns."""
    return await eval_store.query("""
        SELECT
            failure_mode,
            criterion,
            COUNT(*) as count,
            AVG(CASE WHEN severity = 'critical' THEN 3
                      WHEN severity = 'high' THEN 2
                      ELSE 1 END) as avg_severity
        FROM failure_modes
        GROUP BY failure_mode, criterion
        ORDER BY count DESC
    """)
```

### Step 2: Generate Targeted Corrections

When the data shows the LLM consistently misapplies a specific criterion,
add targeted guidance to the scaffold:

```python
# failure_corrections.json — built from evaluation data
{
    "PM2": {
        "trigger": "frequency_misinterpretation",
        "occurrences": 23,
        "correction": "IMPORTANT: PM2 (absent from controls) requires the variant to be absent or at extremely low frequency (< 0.0001) in gnomAD. Do NOT apply PM2 if the global allele frequency is above 0.01%. PM2_supporting (moderate→supporting) is appropriate for AF between 0.0001-0.001%. Check max_pop_af, not just global_af.",
        "added_after_eval_session": "2026-02-15"
    },
    "PP3": {
        "trigger": "prediction_overreliance",
        "occurrences": 17,
        "correction": "CAUTION: PP3 (computational evidence) is supporting-level evidence only. Multiple in-silico predictions (SIFT, PolyPhen2, REVEL, CADD) suggesting deleterious effect can satisfy PP3, but this alone is never sufficient for Pathogenic or Likely Pathogenic. Do not let strong computational predictions override weak clinical evidence.",
        "added_after_eval_session": "2026-02-20"
    }
}
```

### Step 3: Inject Corrections Into the Scaffold

```python
def build_acmg_scaffold(evidence: dict) -> dict:
    scaffold = build_base_scaffold(evidence)

    # Load corrections learned from evaluation data
    corrections = load_failure_corrections()

    for criterion in scaffold["criteria_to_evaluate"]:
        code = criterion["code"]
        if code in corrections:
            criterion["known_pitfall"] = corrections[code]["correction"]
            criterion["pitfall_frequency"] = corrections[code]["occurrences"]

    return scaffold
```

Now when the LLM sees the PM2 criterion in the scaffold, it also sees:
"IMPORTANT: PM2 requires the variant to be absent or at extremely low
frequency..." This is prompt engineering informed by evaluation data. The
scaffold evolves as you collect more failure data.

### Step 4: Measure Improvement

```python
@server.tool(name="get_eval_report")
async def get_eval_report(
    after: str | None = None,  # ISO date — only evals after this date
    before: str | None = None,
    group_by: str = "failure_mode",  # or "criterion", "gene", "variant_type"
) -> dict:
    """Generate an evaluation report showing LLM accuracy and failure patterns."""

    results = await eval_store.get_classifications(after=after, before=before)

    report = {
        "total_evaluations": len(results),
        "accuracy": {
            "exact_match": sum(1 for r in results if r.predicted == r.ground_truth) / len(results),
            "within_one_step": sum(1 for r in results if abs(r.pred_idx - r.gt_idx) <= 1) / len(results),
        },
        "confusion_matrix": build_confusion_matrix(results),
        "failure_modes": aggregate_failures(results, group_by),
        "criteria_accuracy": per_criterion_accuracy(results),
        "improvement_over_time": accuracy_by_date(results),  # the key metric
    }

    return {"content": [{"type": "text", "text": format_eval_report(report)}]}
```

The `improvement_over_time` metric is what you publish. It shows: as the
scaffold incorporates more corrections, LLM accuracy improves from X% to Y%
on the same benchmark set. That's a research finding.

---

## The Ground Truth Dataset

The evaluation harness is only as good as the ground truth. Build this
carefully:

### Source: ClinVar Expert-Reviewed Variants

```python
def build_ground_truth_db():
    """
    Curate a benchmark set from ClinVar with these requirements:
    - Review status >= 2 stars (multiple submitters, no conflict)
    - Has detailed evidence summary
    - Covers the spectrum: pathogenic, likely pathogenic, VUS, likely benign, benign
    - Covers multiple gene categories: tumor suppressors, oncogenes,
      metabolic enzymes, structural proteins
    - Includes common LLM failure bait:
      - Pathogenic variants in non-constrained genes
      - Benign variants in cancer genes (TP53 benign polymorphisms exist)
      - VUS with strong computational predictions (tests overreliance)
      - Common variants misclassified as rare
    """
    benchmark_set = {
        # Easy — clear pathogenic, well-studied
        "tier1_clear": [
            {"chrom": "17", "pos": 7674220, "ref": "C", "alt": "T",
             "gene": "TP53", "classification": "Pathogenic",
             "expected_criteria": {"PS1": {"met": True}, "PM1": {"met": True}, "PM2": {"met": True}},
             "difficulty": "easy", "trap": None},
            # ... 50+ clear cases
        ],

        # Medium — requires careful evidence weighing
        "tier2_nuanced": [
            {"chrom": "7", "pos": 55249071, "ref": "C", "alt": "T",
             "gene": "EGFR", "classification": "Likely Pathogenic",
             "expected_criteria": {"PM1": {"met": True}, "PM2": {"met": True}, "PP3": {"met": True}},
             "difficulty": "medium", "trap": "computational_overreliance"},
            # ... 100+ nuanced cases
        ],

        # Hard — designed to catch specific LLM failure modes
        "tier3_adversarial": [
            # Benign variant in a scary gene
            {"chrom": "17", "pos": 7676154, "ref": "G", "alt": "C",
             "gene": "TP53", "classification": "Benign",
             "expected_criteria": {"BA1": {"met": True}},
             "difficulty": "hard", "trap": "gene_name_bias",
             "note": "Common polymorphism in TP53 — tests if LLM overcalls due to gene reputation"},

            # VUS with strong computational predictions
            {"chrom": "2", "pos": 215595164, "ref": "A", "alt": "G",
             "gene": "BARD1", "classification": "Uncertain Significance",
             "expected_criteria": {"PP3": {"met": True}, "PM2": {"met": True}},
             "difficulty": "hard", "trap": "prediction_overreliance",
             "note": "SIFT/PolyPhen say deleterious but no clinical evidence — should remain VUS"},
        ],
    }
    return benchmark_set
```

### Tier 3 "Adversarial" Variants Are the Most Valuable

These are specifically designed to test known LLM weaknesses:

| Trap | What It Tests | Example |
|------|--------------|---------|
| `gene_name_bias` | Does the LLM overcall variants in "scary" genes? | Benign TP53 polymorphism |
| `prediction_overreliance` | Does computational evidence override lack of clinical data? | VUS with deleterious SIFT/PolyPhen but no clinical evidence |
| `frequency_confusion` | Can the LLM correctly interpret population AF thresholds? | Variant at 0.1% AF (below BA1, above PM2_supporting) |
| `review_status_ignored` | Does the LLM notice ClinVar review quality? | 1-star "pathogenic" from single submitter |
| `conflicting_evidence` | How does the LLM handle contradictory data? | ClinVar pathogenic + gnomAD AF of 1% |
| `artifact_variant` | Does the LLM flag low-quality observations? | 3 alt reads at 8x depth in homopolymer |

---

## Batch Evaluation Mode

For systematic benchmarking, run the full suite:

```python
@server.tool(name="run_evaluation")
async def run_evaluation(
    tier: str = "all",           # tier1_clear, tier2_nuanced, tier3_adversarial, all
    sample_size: int | None = None,  # random sample from tier, or all
    include_corrections: bool = True,  # use scaffold corrections or baseline
) -> dict:
    """
    Run the evaluation benchmark. This will present variants one at a time.
    For each variant, classify it using the evidence provided.

    After all variants are classified, a summary report will be generated.
    """

    variants = load_benchmark_set(tier, sample_size)

    session_id = await eval_store.create_session(
        tier=tier,
        variant_count=len(variants),
        corrections_enabled=include_corrections,
    )

    return {
        "content": [{
            "type": "text",
            "text": (
                f"Evaluation session {session_id} started.\n"
                f"Tier: {tier}\n"
                f"Variants to classify: {len(variants)}\n"
                f"Corrections enabled: {include_corrections}\n\n"
                f"I will now present the first variant. For each, call "
                f"classify_variant with the coordinates, reason through "
                f"the evidence, then call submit_classification with your "
                f"assessment.\n\n"
                f"First variant: {variants[0]['gene']} "
                f"{variants[0]['chrom']}:{variants[0]['pos']} "
                f"{variants[0]['ref']}>{variants[0]['alt']}"
            )
        }],
        "structuredContent": {
            "session_id": session_id,
            "variant_queue": [
                {"chrom": v["chrom"], "pos": v["pos"], "ref": v["ref"], "alt": v["alt"]}
                for v in variants
            ],
        }
    }
```

This creates a conversation where the LLM classifies variants in sequence.
Each classification is scored. At the end, `get_eval_report` produces the
summary.

Run it twice: once with `include_corrections=False` (baseline), once with
`include_corrections=True` (improved). The delta IS the paper.

---

## Eval Store Schema

```sql
CREATE TABLE eval_sessions (
    id TEXT PRIMARY KEY,
    created_at TIMESTAMP DEFAULT CURRENT_TIMESTAMP,
    tier TEXT,
    variant_count INTEGER,
    corrections_enabled BOOLEAN,
    model_id TEXT,            -- "claude-opus-4-5-20251101" etc
    completed_at TIMESTAMP,
    accuracy_exact REAL,
    accuracy_within_one REAL
);

CREATE TABLE tool_invocations (
    id TEXT PRIMARY KEY,
    session_id TEXT REFERENCES eval_sessions(id),
    tool_name TEXT NOT NULL,
    arguments JSON,
    response JSON,
    ground_truth JSON,         -- NULL if no ground truth available
    timestamp TIMESTAMP DEFAULT CURRENT_TIMESTAMP,
    duration_ms INTEGER
);

CREATE TABLE classifications (
    id TEXT PRIMARY KEY,
    invocation_id TEXT REFERENCES tool_invocations(id),
    session_id TEXT REFERENCES eval_sessions(id),
    variant TEXT NOT NULL,      -- "chr17:7674220C>T"
    gene TEXT,
    predicted TEXT NOT NULL,    -- LLM's classification
    ground_truth TEXT,          -- from benchmark
    confidence TEXT,            -- LLM's stated confidence
    criteria_applied JSON,     -- full structured criteria
    reasoning_summary TEXT,
    score JSON,                -- exact match, within-one, per-criterion
    timestamp TIMESTAMP DEFAULT CURRENT_TIMESTAMP
);

CREATE TABLE failure_modes (
    id INTEGER PRIMARY KEY AUTOINCREMENT,
    classification_id TEXT REFERENCES classifications(id),
    session_id TEXT REFERENCES eval_sessions(id),
    mode TEXT NOT NULL,          -- FailureMode enum value
    severity TEXT NOT NULL,      -- critical, high, medium, low
    criterion TEXT,              -- which ACMG criterion, if applicable
    detail TEXT,
    llm_reasoning TEXT,
    expected_reasoning TEXT
);

-- Indexes for common queries
CREATE INDEX idx_failures_mode ON failure_modes(mode);
CREATE INDEX idx_failures_criterion ON failure_modes(criterion);
CREATE INDEX idx_classifications_gene ON classifications(gene);
CREATE INDEX idx_classifications_session ON classifications(session_id);

-- View: failure mode summary
CREATE VIEW failure_summary AS
SELECT
    mode,
    criterion,
    COUNT(*) as occurrences,
    COUNT(DISTINCT classification_id) as affected_classifications,
    GROUP_CONCAT(DISTINCT gene) as affected_genes,
    ROUND(AVG(CASE severity
        WHEN 'critical' THEN 4 WHEN 'high' THEN 3
        WHEN 'medium' THEN 2 ELSE 1 END), 2) as avg_severity
FROM failure_modes
GROUP BY mode, criterion
ORDER BY occurrences DESC;
```

---

## What You Can Publish

After running evaluations across models and scaffold iterations:

### Findings Template

1. **Baseline LLM accuracy on ACMG variant classification**: X% exact match,
   Y% within one classification step, across N variants in M genes.

2. **Systematic failure modes identified**:
   - Gene name bias: LLMs overcall pathogenicity in well-known cancer genes by Z%
   - Computational prediction overreliance: X% of VUS misclassified when SIFT/PolyPhen
     both predict deleterious
   - Population frequency misinterpretation: X% error rate on PM2/BA1 threshold
     application
   - [etc]

3. **Structured prompting intervention**: By providing ACMG reasoning scaffolds
   with targeted corrections for known failure modes, accuracy improved from
   A% to B% (p < 0.05) on the adversarial tier.

4. **Cross-model comparison**: Claude Opus vs Sonnet vs GPT-4o vs Gemini on
   the same benchmark. (MCP works with multiple hosts.)

5. **Recommendation**: Specific prompt engineering strategies for improving
   LLM genomic reasoning, backed by quantitative evaluation.

### Artifact: The Benchmark + Evaluation Harness

Release the benchmark set, the evaluation harness (the MCP server), and the
scaffold correction catalog as open-source. Other researchers can:
- Run the benchmark against new models
- Add variants to the benchmark
- Contribute corrections to the scaffold
- Extend to new domains (pharmacogenomics, somatic variants, CNVs)

---

## "Can I expose another tool for fine-tuned responses?"

Not fine-tuned in the ML sense (you can't fine-tune Claude from within MCP),
but you can build something that functions similarly: a **guided reasoning
tool** that provides calibrated, evidence-informed context that shapes how the
LLM responds.

The scaffold IS this tool. It's effectively a dynamic system prompt per variant
that incorporates lessons from evaluation. As you accumulate evaluation data,
the scaffold gets better, and the LLM's accuracy improves — without any model
training. You're doing prompt engineering at scale, informed by structured
evaluation.

If you wanted to go further, you could build a `get_expert_context` tool that
retrieves curated expert reasoning for known variants:

```python
@server.tool(name="get_expert_context")
async def get_expert_context(
    chrom: str, pos: int, ref: str, alt: str
) -> dict:
    """
    Returns curated expert-level context for a known variant.
    Includes: correct classification, evidence summary, common
    misinterpretations, and key distinguishing factors.

    Use this to verify your reasoning or when confidence is low.
    """
    # Lookup in curated knowledge base
    expert = await knowledge_base.lookup(chrom, pos, ref, alt)

    if expert:
        return {
            "content": [{
                "type": "text",
                "text": (
                    f"Expert context for {expert['gene']} {expert['hgvs']}:\n\n"
                    f"Established classification: {expert['classification']} "
                    f"({expert['review_stars']}-star review)\n\n"
                    f"Key evidence: {expert['key_evidence']}\n\n"
                    f"Common misinterpretation: {expert['common_error']}\n\n"
                    f"This context is from curated expert sources. Use it to "
                    f"calibrate your assessment, not as a substitute for "
                    f"evaluating the evidence yourself."
                )
            }]
        }
    else:
        return {
            "content": [{
                "type": "text",
                "text": "No curated expert context available for this variant. "
                        "Rely on the evidence provided by classify_variant."
            }]
        }
```

This creates a retrieval-augmented classification pipeline: the LLM gets raw
evidence from `classify_variant`, can optionally check `get_expert_context`
for known variants, and submits its classification via `submit_classification`.
The whole chain is logged and scored.

---

## Implementation Priority

1. **SQLite eval store + instrumentation layer** (the foundation — 2 days)
2. **classify_variant + submit_classification tools** (the core loop — 3 days)
3. **Ground truth dataset: Tier 1 clear cases** (50 variants — 2 days of curation)
4. **ACMG scaffold builder** (criterion selection logic — 2 days)
5. **Failure mode detection** (scoring logic — 2 days)
6. **get_eval_report tool** (aggregation and reporting — 1 day)
7. **Run baseline evaluation** (first real data — 1 day)
8. **Tier 3 adversarial variants** (the interesting test cases — 3 days)
9. **Scaffold corrections from baseline failures** (the improvement — 2 days)
10. **Run corrected evaluation, measure delta** (the result — 1 day)
11. **Cross-model evaluation** (run against GPT-4o, Gemini — 1 day per model)
12. **Write up findings** (the publishable artifact)

Total: ~3-4 weeks of focused work to have a publishable evaluation with
quantitative results.