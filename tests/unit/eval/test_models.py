"""Model IDs live only in bamcp.eval.models — guard against drift back into literals."""

from __future__ import annotations

import re
from pathlib import Path

import pytest

from bamcp.eval import models

REPO = Path(__file__).resolve().parents[3]

# A concrete model ID: a Claude family name followed by a version, or an OpenAI gpt-N id.
MODEL_ID_RE = re.compile(r"\bclaude-(?:opus|sonnet|haiku|fable|mythos)-\d[\w.-]*|\bgpt-\d[\w.-]*")

# Code and CI config that must take model IDs from models.py. Docs (EVAL_RESULTS.md records
# historical runs; archived/ is history) and tests (which use IDs as fixtures) are exempt.
_SCANNED = ("src", "scripts", "Makefile", ".github", "pyproject.toml", "docker-compose.yml")
_ALLOWED = {"src/bamcp/eval/models.py"}


_SKIP_DIRS = {"node_modules", "dist", "__pycache__", ".mypy_cache", ".ruff_cache"}
_TEXT_SUFFIXES = {"", ".py", ".ts", ".js", ".html", ".yml", ".yaml", ".toml", ".sh", ".json"}


def _scanned_files() -> list[Path]:
    files: list[Path] = []
    for entry in _SCANNED:
        root = REPO / entry
        candidates = [root] if root.is_file() else root.rglob("*") if root.is_dir() else []
        for f in candidates:
            if not f.is_file() or _SKIP_DIRS.intersection(f.relative_to(REPO).parts):
                continue
            if f.suffix in _TEXT_SUFFIXES and "egg-info" not in str(f):
                files.append(f)
    return files


@pytest.mark.unit
def test_no_model_id_literal_outside_models_py():
    offenders = []
    files = _scanned_files()
    assert any(f.name == "compare.py" for f in files), "scan found no eval sources"
    for path in files:
        rel = path.relative_to(REPO).as_posix()
        if rel in _ALLOWED:
            continue
        try:
            text = path.read_text(encoding="utf-8")
        except (UnicodeDecodeError, OSError):
            continue  # binary / unreadable
        for lineno, line in enumerate(text.splitlines(), 1):
            if m := MODEL_ID_RE.search(line):
                offenders.append(f"{rel}:{lineno}: {m.group(0)}")
    assert not offenders, "model-ID literals outside eval/models.py:\n" + "\n".join(offenders)


@pytest.mark.unit
def test_scan_pattern_catches_ids():
    for literal in (models.CLAUDE_OPUS, models.CLAUDE_FABLE, models.OPENAI_DEFAULT):
        assert MODEL_ID_RE.search(f'default="{literal}"')
    assert not MODEL_ID_RE.search("claude-x / claude-y")


@pytest.mark.unit
def test_roles_are_consistent():
    assert models.DEFAULT_SUBJECT_MODEL in models.COMPARISON_MODELS
    # The judge never grades a comparison subject.
    assert models.JUDGE_MODEL not in models.COMPARISON_MODELS
    assert len(set(models.COMPARISON_MODELS)) == len(models.COMPARISON_MODELS)


@pytest.mark.unit
def test_defaults_wired_through_models():
    from bamcp.eval import providers, schema

    assert providers.DEFAULT_ANTHROPIC_MODEL == models.DEFAULT_SUBJECT_MODEL
    assert providers.DEFAULT_JUDGE_MODEL == models.JUDGE_MODEL
    assert providers.DEFAULT_OPENAI_MODEL == models.OPENAI_DEFAULT
    assert schema.RunConfig(test_cases_path=Path("x"), output_dir=Path("y")).model == (
        models.DEFAULT_SUBJECT_MODEL
    )


@pytest.mark.unit
@pytest.mark.parametrize(
    ("model", "accepts"),
    [
        (models.CLAUDE_OPUS, False),
        (models.CLAUDE_SONNET, False),
        (models.CLAUDE_HAIKU, False),
        (models.CLAUDE_FABLE, False),
        ("claude-opus-5", False),
        ("claude-opus-4-8", False),
        ("claude-opus-4-7", False),
        ("claude-opus-4-6", True),
        ("claude-sonnet-4-6", True),
        ("claude-haiku-4-5", True),
        ("claude-x", True),
        ("gpt-4o-mini", True),
    ],
)
def test_accepts_sampling_params(model, accepts):
    assert models.accepts_sampling_params(model) is accepts
