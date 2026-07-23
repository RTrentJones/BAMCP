"""Model comparison + tool-use ablation for the BAMCP eval harness.

The single-run harness (`bamcp.eval.runner`) answers "did this model pass these
cases". This module answers the two questions a deployment actually cares about:

1. **Which model?** Run the same cases across several (provider, model) configs
   and rank them by accuracy — overall and per category.
2. **Are the tools worth it?** Run each model *with* the BAMCP tools and *without*
   them (same cases, tools withheld), and report the accuracy lift the tools buy.

It reuses the existing runner and grader unchanged; it only orchestrates
multiple runs and aggregates them into a comparison table (Markdown + JSON).

Real runs need ``ANTHROPIC_API_KEY`` / ``OPENAI_API_KEY`` and the ``[eval]``
extras. The ``mock`` provider makes the whole pipeline runnable offline (see
``tests/unit/eval/test_compare.py`` and the ``--provider mock`` self-check), so
the table generation is proven without a key — a real run is a drop-in swap of
the spec's provider/model.
"""

from __future__ import annotations

import json
import re
from collections.abc import Iterable
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

from .providers import (
    LLMProvider,
    MockProvider,
    ProviderResponse,
    ToolCall,
    ToolDescriptor,
    get_provider,
)
from .router import InProcessRouter, RouterResult, ToolRouter
from .runner import run_cases
from .schema import EvalCase, EvalResult, apply_subset, load_cases

# Local (offline, no-network) tools the mock is willing to "call". Excludes
# ClinVar / gnomAD / gene search, which need the network.
_LOCAL_TOOLS = frozenset(
    {
        "get_coverage",
        "get_variants",
        "list_contigs",
        "jump_to",
        "visualize_region",
        "get_region_summary",
        "scan_variants",
        "get_variant_curation_summary",
    }
)


@dataclass(frozen=True)
class RunSpec:
    """One configuration to evaluate: a provider/model, with tools on or off."""

    label: str
    provider: str
    model: str
    tools: bool = True

    @property
    def slug(self) -> str:
        return re.sub(r"[^a-z0-9]+", "-", self.label.lower()).strip("-")


@dataclass
class SpecResult:
    """Aggregated outcome of running one spec over the case set."""

    spec: RunSpec
    passed: int
    total: int
    by_category: dict[str, tuple[int, int]] = field(default_factory=dict)  # cat -> (passed,total)

    @property
    def accuracy(self) -> float:
        return self.passed / self.total if self.total else 0.0

    def as_dict(self) -> dict[str, Any]:
        return {
            "label": self.spec.label,
            "provider": self.spec.provider,
            "model": self.spec.model,
            "tools": self.spec.tools,
            "passed": self.passed,
            "total": self.total,
            "accuracy": round(self.accuracy, 4),
            "by_category": {
                cat: {"passed": p, "total": t, "accuracy": round(p / t, 4) if t else 0.0}
                for cat, (p, t) in sorted(self.by_category.items())
            },
        }


@dataclass
class ComparisonReport:
    """The comparison across all specs."""

    specs: list[SpecResult]

    @property
    def categories(self) -> list[str]:
        cats: set[str] = set()
        for s in self.specs:
            cats.update(s.by_category)
        return sorted(cats)

    def ablation_pairs(self) -> list[tuple[str, SpecResult, SpecResult]]:
        """Match tools-on / tools-off specs by (provider, model) for the lift table."""
        by_key: dict[tuple[str, str], dict[bool, SpecResult]] = {}
        for s in self.specs:
            by_key.setdefault((s.spec.provider, s.spec.model), {})[s.spec.tools] = s
        pairs = []
        for (_provider, model), variants in by_key.items():
            if True in variants and False in variants:
                pairs.append((model, variants[True], variants[False]))
        return pairs

    def as_dict(self) -> dict[str, Any]:
        return {
            "specs": [s.as_dict() for s in self.specs],
            "ablation": [
                {
                    "model": model,
                    "with_tools": round(on.accuracy, 4),
                    "without_tools": round(off.accuracy, 4),
                    "lift": round(on.accuracy - off.accuracy, 4),
                }
                for model, on, off in self.ablation_pairs()
            ],
        }


# ── Aggregation ─────────────────────────────────────────────────────────────


def aggregate(spec: RunSpec, results: Iterable[EvalResult]) -> SpecResult:
    """Fold a spec's per-case results into pass counts overall and by category."""
    passed = total = 0
    by_cat: dict[str, list[int]] = {}
    for r in results:
        if r.skipped:
            continue
        total += 1
        ok = 1 if r.verdict.passed else 0
        passed += ok
        bucket = by_cat.setdefault(r.category, [0, 0])
        bucket[0] += ok
        bucket[1] += 1
    return SpecResult(
        spec=spec,
        passed=passed,
        total=total,
        by_category={c: (p, t) for c, (p, t) in by_cat.items()},
    )


# ── Markdown rendering ──────────────────────────────────────────────────────


def _pct(x: float) -> str:
    return f"{x * 100:.0f}%"


def format_markdown(report: ComparisonReport, *, title: str = "BAMCP model comparison") -> str:
    lines = [f"# {title}", ""]

    lines += ["## Accuracy by configuration", ""]
    lines += ["| Configuration | Provider | Model | Tools | Passed | Cases | Accuracy |"]
    lines += ["|---|---|---|---|---:|---:|---:|"]
    for s in sorted(report.specs, key=lambda s: s.accuracy, reverse=True):
        lines.append(
            f"| {s.spec.label} | {s.spec.provider} | `{s.spec.model}` | "
            f"{'on' if s.spec.tools else 'off'} | {s.passed} | {s.total} | "
            f"**{_pct(s.accuracy)}** |"
        )
    lines.append("")

    pairs = report.ablation_pairs()
    if pairs:
        lines += ["## Tool-use ablation (does the tool help?)", ""]
        lines += ["| Model | With tools | Without tools | Lift from tools |"]
        lines += ["|---|---:|---:|---:|"]
        for model, on, off in pairs:
            lift = on.accuracy - off.accuracy
            lines.append(
                f"| `{model}` | {_pct(on.accuracy)} | {_pct(off.accuracy)} | "
                f"**{lift * 100:+.0f} pts** |"
            )
        lines.append("")

    cats = report.categories
    if cats:
        lines += ["## Accuracy by category", ""]
        header = "| Category | " + " | ".join(s.spec.label for s in report.specs) + " |"
        sep = "|---|" + "|".join("---:" for _ in report.specs) + "|"
        lines += [header, sep]
        for cat in cats:
            cells = []
            for s in report.specs:
                p, t = s.by_category.get(cat, (0, 0))
                cells.append(_pct(p / t) if t else "—")
            lines.append(f"| {cat} | " + " | ".join(cells) + " |")
        lines.append("")

    return "\n".join(lines)


# ── Providers / routers per spec ────────────────────────────────────────────


class _NoToolRouter:
    """Router that advertises no tools — the 'without tools' ablation arm."""

    def list_tools(self) -> list[str]:
        return []

    async def call(self, name: str, arguments: dict[str, Any]) -> RouterResult:
        return RouterResult(ok=False, text="", error="tools disabled for this run")


_REGION_RE = re.compile(r"(chr[\w]+:\d+-\d+)")
_BAM_RE = re.compile(r"(\S+\.(?:bam|cram))")
_REF_RE = re.compile(r"reference\s+(\S+\.fa)")
_VARIANT_RE = re.compile(r"(chr[\w]+):(\d+)[:\s]+([ACGT]+)\s*>\s*([ACGT]+)", re.IGNORECASE)


def _extract_args(text: str) -> dict[str, Any]:
    """Best-effort parse of tool arguments from a case prompt (what a model does)."""
    args: dict[str, Any] = {}
    if m := _BAM_RE.search(text):
        args["file_path"] = m.group(1)
    if m := _REGION_RE.search(text):
        args["region"] = m.group(1)
    if m := _REF_RE.search(text):
        args["reference"] = m.group(1)
    elif args.get("file_path", "").endswith("comprehensive.bam"):
        args["reference"] = "tests/fixtures/comprehensive_ref.fa"
    if m := _VARIANT_RE.search(text):
        args["chrom"] = m.group(1)
        args["pos"] = int(m.group(2))
        args["ref"] = m.group(3).upper()
        args["alt"] = m.group(4).upper()
    if "format=" in text or "rubric" in text:
        args["format"] = "rubric"
    return args


def realistic_mock() -> MockProvider:
    """A mock that behaves like a capable model: parse the prompt, call the
    advertised local tools with extracted args, then summarize.

    With tools advertised it satisfies ``tools_expected``; with none advertised
    (the ablation arm) it can only answer from text — exactly the contrast the
    ablation measures.
    """

    def script(messages: list[dict[str, Any]], tools: list[ToolDescriptor]) -> ProviderResponse:
        # Second turn onward (tool results already appended) -> final answer.
        already_called = any(
            m["role"] == "user" and isinstance(m.get("content"), list) for m in messages[1:]
        )
        prompt = messages[0]["content"] if messages else ""
        prompt = prompt if isinstance(prompt, str) else ""
        if already_called or not tools:
            answer = "Based on the available information, here is the analysis." + (
                "" if tools else " (No tools were available for this run.)"
            )
            return ProviderResponse(text=answer, tool_calls=[])
        args = _extract_args(prompt)
        calls = [
            ToolCall(id=f"call_{i}", name=t.name, arguments=dict(args))
            for i, t in enumerate(tools)
            if t.name in _LOCAL_TOOLS
        ]
        return ProviderResponse(text="Let me use the BAMCP tools.", tool_calls=calls)

    return MockProvider(script)


def build_provider(spec: RunSpec) -> LLMProvider:
    if spec.provider == "mock":
        return realistic_mock()
    return get_provider(spec.provider, spec.model)


def build_judge(spec: RunSpec) -> LLMProvider | None:
    if spec.provider == "mock":
        return None
    return get_provider(spec.provider, spec.model)


def build_router(spec: RunSpec) -> ToolRouter:
    return InProcessRouter() if spec.tools else _NoToolRouter()


# ── Orchestration ───────────────────────────────────────────────────────────


async def run_comparison(
    specs: list[RunSpec],
    cases: list[EvalCase],
    output_dir: Path,
) -> ComparisonReport:
    """Run every spec over the shared case set and aggregate the outcomes."""
    spec_results: list[SpecResult] = []
    for spec in specs:
        results = await run_cases(
            cases=cases,
            provider=build_provider(spec),
            router=build_router(spec),
            judge=build_judge(spec),
            output_dir=output_dir / spec.slug,
            no_vision=True,
        )
        spec_results.append(aggregate(spec, results))
    return ComparisonReport(spec_results)


def _default_specs(provider: str, models: list[str]) -> list[RunSpec]:
    """Build tools-on + tools-off specs for each model (the ablation grid)."""
    specs: list[RunSpec] = []
    for model in models:
        short = model if provider == "mock" else model
        specs.append(RunSpec(label=f"{short} (tools)", provider=provider, model=model, tools=True))
        specs.append(
            RunSpec(label=f"{short} (no tools)", provider=provider, model=model, tools=False)
        )
    return specs


def main(argv: list[str] | None = None) -> int:
    """CLI: run a model comparison + tool ablation and write the results table."""
    import argparse
    import asyncio

    p = argparse.ArgumentParser(
        prog="bamcp-eval-compare",
        description="Compare models and measure the tool-use lift for BAMCP.",
    )
    p.add_argument("--provider", default="anthropic", choices=["anthropic", "openai", "mock"])
    p.add_argument(
        "--models",
        nargs="+",
        default=["claude-opus-4-8", "claude-sonnet-4-5"],
        help="Models to compare (each run with and without tools).",
    )
    p.add_argument("--test-cases", default="tests/eval/test_cases.yaml")
    p.add_argument("--output-dir", default=".eval-results/comparison")
    p.add_argument("--subset", help="1-indexed inclusive range, e.g. 1-5")
    p.add_argument("--results-md", default=None, help="Write the Markdown table here.")
    args = p.parse_args(argv)

    subset = None
    if args.subset:
        lo, hi = args.subset.split("-")
        subset = (int(lo), int(hi))
    cases = apply_subset(load_cases(args.test_cases), subset)

    specs = _default_specs(args.provider, args.models)
    output_dir = Path(args.output_dir)
    report = asyncio.run(run_comparison(specs, cases, output_dir))

    output_dir.mkdir(parents=True, exist_ok=True)
    (output_dir / "comparison.json").write_text(
        json.dumps(report.as_dict(), indent=2), encoding="utf-8"
    )
    md = format_markdown(report)
    md_path = Path(args.results_md) if args.results_md else output_dir / "comparison.md"
    md_path.write_text(md, encoding="utf-8")
    print(md)
    print(f"\nWrote {md_path} and {output_dir / 'comparison.json'}")
    return 0


if __name__ == "__main__":  # pragma: no cover
    import sys

    sys.exit(main())
