"""Schema definitions for the BAMCP eval harness.

The YAML format mirrors MARRVEL-MCP's ``mcp_llm_test/test_cases.yaml`` exactly
(``name``/``category``/``input``/``expected``) and adds two optional
BAMCP-specific fields:

- ``tools_expected``: list of tool names that must appear in the captured
  telemetry for the case to pass deterministic grading.
- ``bam_fixture``: path to a BAM/CRAM fixture the eval runner injects into
  the prompt.
- ``rendering_modes``: list of viewer display modes; when set, the case is
  run once per mode for head-to-head comparison.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path
from typing import Any


@dataclass
class EvalCase:
    """A single benchmark case."""

    name: str
    input: str
    expected: str
    category: str = "uncategorized"
    tools_expected: list[str] = field(default_factory=list)
    bam_fixture: str | None = None
    rendering_modes: list[str] = field(default_factory=list)
    # When True, the case can only be answered with vision. Non-vision runs
    # skip it (and the report marks it as "skipped: vision unavailable").
    vision_required: bool = False

    @classmethod
    def from_dict(cls, d: dict[str, Any]) -> EvalCase:
        """Build an EvalCase from a parsed YAML entry."""
        # MARRVEL nests cases under a top-level ``case:`` key.
        if "case" in d and isinstance(d["case"], dict):
            d = d["case"]
        return cls(
            name=str(d["name"]),
            input=str(d["input"]),
            expected=str(d["expected"]),
            category=str(d.get("category", "uncategorized")),
            tools_expected=list(d.get("tools_expected", []) or []),
            bam_fixture=d.get("bam_fixture"),
            rendering_modes=list(d.get("rendering_modes", []) or []),
            vision_required=bool(d.get("vision_required", False)),
        )


@dataclass
class GraderVerdict:
    """The grader's verdict for a single case."""

    passed: bool
    rationale: str
    deterministic: bool  # True if a deterministic check decided; False = LLM-judged
    score: float = 0.0  # 0.0 to 1.0 for partial credit (deterministic = 0 or 1)


@dataclass
class EvalResult:
    """The full result of running one case end-to-end."""

    case_name: str
    category: str
    rendering_mode: str | None
    response_text: str
    tool_calls: list[str]
    duration_ms: float
    verdict: GraderVerdict
    error: str | None = None
    telemetry_path: str | None = None  # path to the JSONL trace for this case
    # Paths to screenshots captured during this case (when vision was active).
    screenshots: list[str] = field(default_factory=list)
    # True when at least one image was passed to the LLM for this case.
    vision_used: bool = False
    # True when the case was skipped because vision was unavailable but
    # ``vision_required`` was set on the case.
    skipped: bool = False
    skip_reason: str | None = None


@dataclass
class RunConfig:
    """Configuration for an eval run (parsed from CLI flags)."""

    test_cases_path: Path
    output_dir: Path
    provider: str = "anthropic"
    model: str = "claude-opus-4-7"
    judge_model: str | None = None
    cache: bool = False
    subset: tuple[int, int] | None = None  # 1-indexed inclusive range
    with_vanilla: bool = False
    with_rendering_comparison: bool = False
    router: str = "in-process"  # "in-process" or "mcp-stdio"
    dry_run: bool = False
    # Force text-only runs even when cases declare rendering_modes / provider
    # supports vision. Useful for cost control.
    no_vision: bool = False
    vision_screenshot_dir: Path | None = None


def load_cases(path: str | Path) -> list[EvalCase]:
    """Parse a YAML file of benchmark cases.

    Accepts both flat list (``- name: ...``) and MARRVEL-style nested
    (``- case: { name: ... }``) layouts.
    """
    try:
        import yaml  # type: ignore[import-untyped]
    except ImportError as e:  # pragma: no cover — install hint
        raise RuntimeError(
            "pyyaml is required for the eval harness. Install with 'pip install \".[eval]\"'."
        ) from e
    text = Path(path).read_text(encoding="utf-8")
    parsed = yaml.safe_load(text) or []
    if not isinstance(parsed, list):
        raise ValueError(f"Top-level YAML must be a list of cases; got {type(parsed).__name__}")
    return [EvalCase.from_dict(item) for item in parsed]


def apply_subset(cases: list[EvalCase], subset: tuple[int, int] | None) -> list[EvalCase]:
    """Filter cases by a 1-indexed inclusive subset (e.g. (1, 5))."""
    if subset is None:
        return cases
    start, end = subset
    return cases[max(0, start - 1) : end]
