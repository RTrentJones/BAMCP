"""Grading logic for eval cases.

Two-tier:

1. Deterministic — **tool selection and answer correctness are scored separately**:
   - ``tools_expected`` (tool selection): every named tool must appear in telemetry.
   - ``expected_claims`` (answer correctness): required factual/numeric strings must appear
     in the assistant's answer.
   - ``score_at_least: <key>=<float>`` in ``expected`` (answer correctness): checked against
     the rubric payload the curation tool emits.
   A case fails if *any* declared check fails, so satisfying ``tools_expected`` no longer
   passes a case whose answer is wrong. Sub-scores (``tool_score``/``answer_score``) are
   surfaced on the verdict.

2. LLM judge: free-form ``expected`` strings are compared against the assistant response via a
   small judge prompt. The judge model defaults to the same provider/model running the eval;
   users can override it.
"""

from __future__ import annotations

import json
import re
from collections.abc import Iterable

from .providers import LLMProvider
from .schema import EvalCase, GraderVerdict

_JUDGE_PROMPT = (
    "You are grading an AI assistant's answer to a genomics question.\n"
    "Question:\n{question}\n\n"
    "Expected criteria:\n{expected}\n\n"
    "Assistant's response:\n{response}\n\n"
    "Reply with a single line in this exact format:\n"
    "VERDICT: PASS — <one-sentence rationale>\n"
    "or\n"
    "VERDICT: FAIL — <one-sentence rationale>\n"
)

_SCORE_THRESHOLD_PATTERN = re.compile(r"score_at_least:\s*(\w+)\s*=\s*([0-9]*\.?[0-9]+)")


def _missing_claims(expected_claims: Iterable[str], response_text: str) -> list[str]:
    """Return the required claims not present in the response as a bounded token.

    Matches on alphanumeric boundaries (not raw substring), so ``chr1`` does not spuriously
    match inside ``chr10`` and ``43.5%`` does not match inside ``143.5%``. This checks literal
    token *presence* only — it does not model negation (e.g. "not pathogenic" still contains
    "pathogenic"); use the LLM judge or distinct claims for negation-sensitive cases.
    """
    text = response_text or ""
    missing: list[str] = []
    for claim in expected_claims:
        c = claim.strip()
        if not c:
            continue
        pattern = r"(?<![A-Za-z0-9])" + re.escape(c) + r"(?![A-Za-z0-9])"
        if not re.search(pattern, text, re.IGNORECASE):
            missing.append(claim)
    return missing


def grade_case(
    case: EvalCase,
    response_text: str,
    tool_calls: Iterable[str],
    rubric_payloads: Iterable[dict] = (),
) -> GraderVerdict | None:
    """Run deterministic checks. Return None if no deterministic rule applies.

    Tool selection and answer correctness are graded **separately**: satisfying
    ``tools_expected`` no longer passes a case on its own when ``expected_claims`` (required
    factual/numeric content) are declared — those must appear in the answer too. Sub-scores are
    surfaced on the verdict so a report can show which dimension failed.
    """
    tool_set = set(tool_calls)
    deterministic_signals = []
    failures = []
    tool_score: float | None = None
    answer_score: float | None = None

    # --- Tool selection -----------------------------------------------------
    if case.tools_expected:
        missing = [t for t in case.tools_expected if t not in tool_set]
        tool_score = 1.0 if not missing else 0.0
        if missing:
            failures.append(f"missing expected tool calls: {missing}")
        else:
            deterministic_signals.append("tools_expected satisfied")

    # --- Answer correctness: required literal/numeric claims ----------------
    if case.expected_claims:
        missing_claims = _missing_claims(case.expected_claims, response_text)
        answer_score = 1.0 if not missing_claims else 0.0
        if missing_claims:
            failures.append(f"answer missing required claims: {missing_claims}")
        else:
            deterministic_signals.append(f"all {len(case.expected_claims)} expected claims present")

    # --- Answer correctness: rubric score thresholds -----------------------
    score_thresholds = _SCORE_THRESHOLD_PATTERN.findall(case.expected)
    if score_thresholds:
        merged: dict[str, float] = {}
        for payload in rubric_payloads:
            scores = payload.get("scores") if isinstance(payload, dict) else None
            if isinstance(scores, dict):
                for k, v in scores.items():
                    if isinstance(v, (int, float)):
                        merged[k] = max(merged.get(k, 0.0), float(v))
        threshold_ok = True
        for key, threshold_str in score_thresholds:
            threshold = float(threshold_str)
            actual = merged.get(key)
            if actual is None:
                failures.append(f"score '{key}' not present in any rubric payload")
                threshold_ok = False
            elif actual < threshold:
                failures.append(f"score '{key}'={actual:.3f} below threshold {threshold:.3f}")
                threshold_ok = False
        answer_score = min(
            answer_score if answer_score is not None else 1.0, 1.0 if threshold_ok else 0.0
        )
        if threshold_ok:
            deterministic_signals.append(
                "rubric thresholds satisfied: "
                + ", ".join(f"{k}>={v}" for k, v in score_thresholds)
            )

    # No deterministic rule applied → defer to the LLM judge.
    if not deterministic_signals and not failures:
        return None

    passed = not failures
    rationale = "; ".join(failures) if failures else "; ".join(deterministic_signals)
    return GraderVerdict(
        passed=passed,
        rationale=rationale,
        deterministic=True,
        score=1.0 if passed else 0.0,
        tool_score=tool_score,
        answer_score=answer_score,
    )


async def judge_with_llm(
    case: EvalCase,
    response_text: str,
    judge: LLMProvider,
) -> GraderVerdict:
    """Ask the LLM judge for a pass/fail verdict on the assistant response."""
    prompt = _JUDGE_PROMPT.format(
        question=case.input, expected=case.expected, response=response_text or "(empty)"
    )
    resp = await judge.chat_with_tools(
        messages=[{"role": "user", "content": prompt}],
        tools=[],
    )
    text = (resp.text or "").strip()
    verdict_line = text.splitlines()[0] if text else ""
    passed = verdict_line.upper().startswith("VERDICT: PASS")
    return GraderVerdict(
        passed=passed,
        rationale=verdict_line or "(no judge output)",
        deterministic=False,
        score=1.0 if passed else 0.0,
    )


def collect_rubric_payloads(tool_results: Iterable[tuple[str, str]]) -> list[dict]:
    """Extract JSON rubric payloads from a sequence of (tool_name, text) pairs.

    Only ``get_variant_curation_summary`` rubric responses are recognized;
    other tool outputs are ignored.
    """
    out: list[dict] = []
    for name, text in tool_results:
        if name != "get_variant_curation_summary":
            continue
        try:
            parsed = json.loads(text)
        except (ValueError, TypeError):
            continue
        if isinstance(parsed, dict) and parsed.get("rubric_version"):
            out.append(parsed)
    return out
