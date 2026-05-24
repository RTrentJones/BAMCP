"""Grading logic for eval cases.

Two-tier:

1. Deterministic: when ``tools_expected`` is set, require those tools to
   appear in the captured telemetry. When ``expected`` contains structured
   markers (e.g. ``"score_at_least: vaf_quality=0.5"``), check them against
   the recorded tool output. Pass deterministically.

2. LLM judge: free-form ``expected`` strings are compared against the
   assistant response via a small judge prompt. The judge model defaults to
   the same provider/model running the eval; users can override it.
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


def grade_case(
    case: EvalCase,
    response_text: str,
    tool_calls: Iterable[str],
    rubric_payloads: Iterable[dict] = (),
) -> GraderVerdict | None:
    """Run deterministic checks. Return None if no deterministic rule applies."""
    tool_set = set(tool_calls)
    deterministic_signals = []

    if case.tools_expected:
        missing = [t for t in case.tools_expected if t not in tool_set]
        if missing:
            return GraderVerdict(
                passed=False,
                rationale=f"Missing expected tool calls: {missing}",
                deterministic=True,
            )
        deterministic_signals.append("tools_expected satisfied")

    score_thresholds = _SCORE_THRESHOLD_PATTERN.findall(case.expected)
    if score_thresholds:
        merged: dict[str, float] = {}
        for payload in rubric_payloads:
            scores = payload.get("scores") if isinstance(payload, dict) else None
            if isinstance(scores, dict):
                for k, v in scores.items():
                    if isinstance(v, (int, float)):
                        merged[k] = max(merged.get(k, 0.0), float(v))
        for key, threshold_str in score_thresholds:
            threshold = float(threshold_str)
            actual = merged.get(key)
            if actual is None:
                return GraderVerdict(
                    passed=False,
                    rationale=f"Score '{key}' not present in any rubric payload",
                    deterministic=True,
                )
            if actual < threshold:
                return GraderVerdict(
                    passed=False,
                    rationale=f"Score '{key}'={actual:.3f} below threshold {threshold:.3f}",
                    deterministic=True,
                )
        deterministic_signals.append(
            "rubric thresholds satisfied: " + ", ".join(f"{k}>={v}" for k, v in score_thresholds)
        )

    if deterministic_signals:
        return GraderVerdict(
            passed=True,
            rationale="; ".join(deterministic_signals),
            deterministic=True,
            score=1.0,
        )
    return None


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
