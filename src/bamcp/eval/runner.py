"""Main eval loop: case → LLM → tool calls → grading → result.

The runner drives one EvalCase to completion:

1. Configure telemetry to write to a per-run JSONL file.
2. Hand the LLM provider the BAMCP tool descriptors plus the case input.
3. Process the model's tool-use loop (up to ``MAX_TURNS`` rounds).
4. Read the captured telemetry to enumerate which tools were called.
5. Grade deterministically; fall back to the LLM judge when no deterministic
   rule applies.
6. Return an :class:`EvalResult` for the case.
"""

from __future__ import annotations

import json
import time
import uuid
from collections.abc import Iterable
from pathlib import Path
from typing import Any

from ..middleware import telemetry
from .grader import collect_rubric_payloads, grade_case, judge_with_llm
from .providers import LLMProvider, ToolCall, ToolDescriptor
from .router import ToolRouter, tool_descriptors
from .schema import EvalCase, EvalResult, GraderVerdict

MAX_TURNS = 8


def _read_jsonl(path: Path) -> list[dict[str, Any]]:
    if not path.exists():
        return []
    out: list[dict[str, Any]] = []
    for line in path.read_text(encoding="utf-8").splitlines():
        if not line.strip():
            continue
        try:
            out.append(json.loads(line))
        except json.JSONDecodeError:
            continue
    return out


def _format_user_message(case: EvalCase, rendering_mode: str | None) -> str:
    parts = [case.input]
    if case.bam_fixture:
        parts.append(f"\n(BAM fixture: {case.bam_fixture})")
    if rendering_mode:
        parts.append(f"\n(Rendering mode under test: {rendering_mode})")
    return "".join(parts)


async def _run_tool_loop(
    provider: LLMProvider,
    router: ToolRouter,
    tools: list[ToolDescriptor],
    initial_user: str,
) -> tuple[str, list[tuple[str, str]]]:
    """Drive a single conversation through up to MAX_TURNS rounds.

    Returns ``(final_text, tool_results)`` where ``tool_results`` is the list
    of ``(tool_name, response_text)`` tuples in call order.
    """
    messages: list[dict[str, Any]] = [{"role": "user", "content": initial_user}]
    tool_results: list[tuple[str, str]] = []
    final_text = ""

    for _ in range(MAX_TURNS):
        resp = await provider.chat_with_tools(messages=messages, tools=tools)
        final_text = resp.text
        if not resp.tool_calls:
            break
        # Append the assistant turn carrying the tool_use blocks.
        messages.append(
            {
                "role": "assistant",
                "content": _assistant_content_blocks(resp.text, resp.tool_calls),
            }
        )
        # Execute each tool and append a user turn with the results.
        user_content: list[dict[str, Any]] = []
        for call in resp.tool_calls:
            outcome = await router.call(call.name, call.arguments)
            tool_results.append((call.name, outcome.text))
            user_content.append(
                {
                    "type": "tool_result",
                    "tool_use_id": call.id,
                    "content": outcome.text or (outcome.error or ""),
                    "is_error": not outcome.ok,
                }
            )
        messages.append({"role": "user", "content": user_content})

    return final_text, tool_results


def _assistant_content_blocks(text: str, tool_calls: list[ToolCall]) -> list[dict[str, Any]]:
    """Build the Anthropic-style content-block list for an assistant turn."""
    blocks: list[dict[str, Any]] = []
    if text:
        blocks.append({"type": "text", "text": text})
    for call in tool_calls:
        blocks.append(
            {"type": "tool_use", "id": call.id, "name": call.name, "input": call.arguments}
        )
    return blocks


async def run_case(
    case: EvalCase,
    provider: LLMProvider,
    router: ToolRouter,
    judge: LLMProvider | None,
    output_dir: Path,
    rendering_mode: str | None = None,
) -> EvalResult:
    """Run one case to completion and grade it."""
    case_slug = case.name.replace(" ", "_").replace("/", "_").replace(":", "_").replace(">", "_")
    suffix = f"__{rendering_mode}" if rendering_mode else ""
    telemetry_path = output_dir / "telemetry" / f"{case_slug}{suffix}.jsonl"
    telemetry_path.parent.mkdir(parents=True, exist_ok=True)
    if telemetry_path.exists():
        telemetry_path.unlink()

    session_id = f"eval-{case_slug}-{uuid.uuid4().hex[:8]}"
    telemetry.configure_telemetry(
        enabled=True, jsonl_path=str(telemetry_path), session_id=session_id
    )

    start = time.monotonic()
    error: str | None = None
    final_text = ""
    tool_results: list[tuple[str, str]] = []
    try:
        tools = tool_descriptors(router)
        user_message = _format_user_message(case, rendering_mode)
        final_text, tool_results = await _run_tool_loop(provider, router, tools, user_message)
    except Exception as e:  # noqa: BLE001 — capture and surface to caller
        error = f"{type(e).__name__}: {e}"

    duration_ms = (time.monotonic() - start) * 1000
    telemetry.reset_telemetry()

    events = _read_jsonl(telemetry_path)
    captured_tools: list[str] = [
        str(e["tool_name"]) for e in events if isinstance(e.get("tool_name"), str)
    ]
    rubric_payloads = collect_rubric_payloads(tool_results)

    verdict = grade_case(case, final_text, captured_tools, rubric_payloads)
    if verdict is None:
        if judge is None:
            verdict = GraderVerdict(
                passed=False,
                rationale="No deterministic rule applied and no judge provider configured",
                deterministic=True,
            )
        else:
            verdict = await judge_with_llm(case, final_text, judge)

    return EvalResult(
        case_name=case.name,
        category=case.category,
        rendering_mode=rendering_mode,
        response_text=final_text,
        tool_calls=captured_tools,
        duration_ms=round(duration_ms, 3),
        verdict=verdict,
        error=error,
        telemetry_path=str(telemetry_path),
    )


async def run_cases(
    cases: Iterable[EvalCase],
    provider: LLMProvider,
    router: ToolRouter,
    judge: LLMProvider | None,
    output_dir: Path,
    with_rendering_comparison: bool = False,
) -> list[EvalResult]:
    """Run a sequence of cases. Returns one result per (case, rendering_mode)."""
    results: list[EvalResult] = []
    for case in cases:
        modes: list[str | None] = (
            list(case.rendering_modes)
            if (with_rendering_comparison and case.rendering_modes)
            else [None]
        )
        for mode in modes:
            results.append(await run_case(case, provider, router, judge, output_dir, mode))
    return results
