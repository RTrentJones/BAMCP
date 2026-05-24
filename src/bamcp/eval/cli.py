"""CLI entry point for the BAMCP eval harness.

Mirrors MARRVEL-MCP's CLI flags so the two harnesses are directly
comparable. Designed to be invoked from ``scripts/run_eval.py`` or via
``python -m bamcp.eval.cli``.
"""

from __future__ import annotations

import argparse
import asyncio
import json
import sys
from pathlib import Path

from .providers import LLMProvider, get_provider
from .report import write_reports, write_summary_line
from .router import InProcessRouter, ToolRouter
from .runner import run_cases
from .schema import RunConfig, apply_subset, load_cases


def parse_args(argv: list[str] | None = None) -> RunConfig:
    p = argparse.ArgumentParser(
        prog="bamcp-eval",
        description=(
            "Run the BAMCP eval harness. Compatible with MARRVEL-MCP's "
            "test_cases.yaml schema for head-to-head comparison."
        ),
    )
    p.add_argument(
        "--test-cases",
        default="tests/eval/test_cases.yaml",
        help="Path to a YAML file of cases (default: tests/eval/test_cases.yaml).",
    )
    p.add_argument(
        "--output-dir",
        default=".eval-results",
        help="Where to write results (default: .eval-results).",
    )
    p.add_argument(
        "--provider",
        default="anthropic",
        choices=["anthropic", "openai", "mock"],
        help="LLM provider (default: anthropic).",
    )
    p.add_argument("--model", default="claude-opus-4-7", help="Provider-specific model id.")
    p.add_argument(
        "--judge-model",
        default=None,
        help="Override model for the LLM judge (defaults to the same model).",
    )
    p.add_argument(
        "--cache",
        action="store_true",
        help="Reuse passing results from a previous run; only re-test failures.",
    )
    p.add_argument(
        "--subset",
        default=None,
        help="1-indexed inclusive range, e.g. '1-5' to run only the first 5 cases.",
    )
    p.add_argument(
        "--with-vanilla",
        action="store_true",
        help="Also run each case without any MCP tools, for delta measurement.",
    )
    p.add_argument(
        "--with-rendering-comparison",
        action="store_true",
        help=(
            "For cases that declare 'rendering_modes', run once per mode and "
            "report per-mode accuracy."
        ),
    )
    p.add_argument(
        "--router",
        default="in-process",
        choices=["in-process", "mcp-stdio"],
        help="How to dispatch tool calls (default: in-process).",
    )
    p.add_argument(
        "--dry-run",
        action="store_true",
        help="Parse cases and exit without running. Validates schema.",
    )
    args = p.parse_args(argv)

    subset: tuple[int, int] | None = None
    if args.subset:
        try:
            a, b = (int(x) for x in args.subset.split("-", 1))
            subset = (a, b)
        except ValueError:
            p.error(f"--subset must look like 'N-M', got {args.subset!r}")

    return RunConfig(
        test_cases_path=Path(args.test_cases),
        output_dir=Path(args.output_dir),
        provider=args.provider,
        model=args.model,
        judge_model=args.judge_model,
        cache=args.cache,
        subset=subset,
        with_vanilla=args.with_vanilla,
        with_rendering_comparison=args.with_rendering_comparison,
        router=args.router,
        dry_run=args.dry_run,
    )


def _filter_cached(results_path: Path, cases):  # noqa: ANN001, ANN201
    """Drop cases that have a passing result in ``results.json``."""
    if not results_path.exists():
        return cases
    try:
        prior = json.loads(results_path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError):
        return cases
    passed = {r.get("case_name") for r in prior if r.get("passed")}
    return [c for c in cases if c.name not in passed]


def _build_router(config: RunConfig) -> ToolRouter:
    if config.router == "in-process":
        return InProcessRouter()
    if config.router == "mcp-stdio":
        raise NotImplementedError("mcp-stdio router is Phase 2. Use --router in-process for now.")
    raise ValueError(f"Unknown router: {config.router}")


def _build_judge(config: RunConfig) -> LLMProvider | None:
    if config.provider == "mock":
        return None  # No judge in fully-offline runs.
    model = config.judge_model or config.model
    return get_provider(config.provider, model)


async def _run(config: RunConfig) -> int:
    if not config.test_cases_path.exists():
        print(f"Test cases file not found: {config.test_cases_path}", file=sys.stderr)
        return 2

    cases = load_cases(config.test_cases_path)
    cases = apply_subset(cases, config.subset)

    if config.dry_run:
        print(f"Parsed {len(cases)} cases from {config.test_cases_path}")
        return 0

    if config.cache:
        cases = _filter_cached(config.output_dir / "results.json", cases)
        if not cases:
            print("All cases cached. Nothing to do.")
            return 0

    provider = get_provider(config.provider, config.model)
    judge = _build_judge(config)
    router = _build_router(config)

    results = await run_cases(
        cases=cases,
        provider=provider,
        router=router,
        judge=judge,
        output_dir=config.output_dir,
        with_rendering_comparison=config.with_rendering_comparison,
    )

    write_reports(results, config, config.output_dir)
    print(write_summary_line(results))
    return 0 if all(r.verdict.passed for r in results) else 1


def main(argv: list[str] | None = None) -> int:
    config = parse_args(argv)
    return asyncio.run(_run(config))


if __name__ == "__main__":  # pragma: no cover
    sys.exit(main())
