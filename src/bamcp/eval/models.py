"""The single source of truth for every model ID the eval harness uses.

Model IDs used to drift across providers.py, schema.py, acmg.py, compare.py, the
Makefile and the nightly workflow. Everything now imports from here (the Makefile
and nightly workflow rely on ``compare.py``'s default, :data:`COMPARISON_MODELS`), and
``tests/unit/eval/test_models.py`` fails if a model-ID literal appears in code or CI
config anywhere else.

Current Claude family (Opus 5.5 / Sonnet 5.5 / Haiku 5.5 / Fable 5.1). IDs are exact —
never append a date suffix.
"""

from __future__ import annotations

import re

# ── Claude ──────────────────────────────────────────────────────────────────
CLAUDE_OPUS = "claude-opus-5-5"
CLAUDE_SONNET = "claude-sonnet-5-5"
CLAUDE_HAIKU = "claude-haiku-5-5"
CLAUDE_FABLE = "claude-fable-5-1"

# ── OpenAI (provider kept for parity; gated in the nightly run) ─────────────
OPENAI_DEFAULT = "gpt-4o-mini"

# ── Roles ───────────────────────────────────────────────────────────────────
# Model under test when none is given (single runs, ACMG eval, RunConfig default).
DEFAULT_SUBJECT_MODEL = CLAUDE_OPUS
# The model comparison grid (`make eval-matrix`, nightly comparison step).
COMPARISON_MODELS: tuple[str, ...] = (CLAUDE_OPUS, CLAUDE_SONNET)
# The fixed, separate LLM judge. Deliberately NOT one of the comparison subjects, so a
# model never grades its own answers (providers.resolve_judge refuses that anyway).
JUDGE_PROVIDER = "anthropic"
JUDGE_MODEL = CLAUDE_HAIKU

# Models that reject explicit sampling parameters (temperature/top_p/top_k): sending one
# returns HTTP 400. The 5.x family and Opus 4.7/4.8 are in this set; determinism there
# comes from the fixed model + prompt, not temperature 0.
_NO_SAMPLING_RE = re.compile(r"^claude-(?:(?:opus|sonnet|haiku|fable|mythos)-5|opus-4-[78])(?:-|$)")


def accepts_sampling_params(model: str) -> bool:
    """Whether ``model`` accepts an explicit ``temperature`` (False → omit it)."""
    return not _NO_SAMPLING_RE.match(model)
