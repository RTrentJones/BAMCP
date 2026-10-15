"""LLM provider abstraction for the eval harness.

Each provider exposes a single async method:

    chat_with_tools(messages, tools) -> ProviderResponse

Implementations:

- ``MockProvider``: scripts a deterministic conversation from a callback.
  Used in tests and ``--provider mock``.
- ``AnthropicProvider``: real Anthropic API with tool-use turn loop.
- ``OpenAIProvider``: real OpenAI API with tool-call turn loop.

Both real providers require the ``[eval]`` extras to be installed and an API
key in the environment (``ANTHROPIC_API_KEY``, ``OPENAI_API_KEY``).
"""

from __future__ import annotations

from collections.abc import Awaitable, Callable
from dataclasses import dataclass, field
from typing import Any, Protocol

from . import models


@dataclass
class ToolDescriptor:
    """A tool exposed to the LLM."""

    name: str
    description: str
    input_schema: dict[str, Any]


@dataclass
class ToolCall:
    """A tool-call request emitted by the LLM."""

    id: str
    name: str
    arguments: dict[str, Any]


@dataclass
class ProviderResponse:
    """One assistant turn from the provider."""

    text: str
    tool_calls: list[ToolCall] = field(default_factory=list)


class LLMProvider(Protocol):
    """Common provider interface."""

    # Whether this provider can accept image content blocks (vision).
    vision_capable: bool
    # Short tag identifying the provider's content-block shape:
    # ``"anthropic"``, ``"openai"``, or ``"mock"``. Used by
    # :func:`build_tool_result_content`.
    kind: str

    async def chat_with_tools(
        self,
        messages: list[dict[str, Any]],
        tools: list[ToolDescriptor],
    ) -> ProviderResponse: ...


# -----------------------------------------------------------------------------
# Mock provider — used by tests and `--provider mock`
# -----------------------------------------------------------------------------


class MockProvider:
    """Deterministic scripted provider.

    Construct with a ``script`` callable that, given (messages, tools), returns
    the next ProviderResponse. The callback receives the full conversation so
    far so it can implement multi-turn behavior. Set ``vision_capable=True`` to
    exercise the image-in-tool-result code path from tests.
    """

    kind: str = "mock"
    vision_capable: bool = False

    def __init__(
        self,
        script: Callable[[list[dict[str, Any]], list[ToolDescriptor]], ProviderResponse]
        | Callable[[list[dict[str, Any]], list[ToolDescriptor]], Awaitable[ProviderResponse]],
        *,
        vision_capable: bool = False,
    ) -> None:
        self._script = script
        self.vision_capable = vision_capable

    async def chat_with_tools(
        self,
        messages: list[dict[str, Any]],
        tools: list[ToolDescriptor],
    ) -> ProviderResponse:
        result = self._script(messages, tools)
        if hasattr(result, "__await__"):
            return await result  # type: ignore[no-any-return]
        return result  # type: ignore[return-value]


def make_one_shot_mock(
    text: str,
    tool_calls: list[tuple[str, dict[str, Any]]] | None = None,
) -> MockProvider:
    """Convenience: a mock that emits ``text`` plus optional tool calls once,
    then a plain text turn on subsequent calls.
    """
    emitted = {"done": False}

    def _script(messages, tools):  # noqa: ANN001
        if emitted["done"]:
            return ProviderResponse(text=text, tool_calls=[])
        emitted["done"] = True
        calls = [
            ToolCall(id=f"call_{i}", name=name, arguments=args)
            for i, (name, args) in enumerate(tool_calls or [])
        ]
        return ProviderResponse(text=text, tool_calls=calls)

    return MockProvider(_script)


# -----------------------------------------------------------------------------
# Anthropic provider
# -----------------------------------------------------------------------------


# Per-provider default model ids (single source of truth for construction + run records).
# The IDs themselves live in models.py.
DEFAULT_ANTHROPIC_MODEL = models.DEFAULT_SUBJECT_MODEL
DEFAULT_OPENAI_MODEL = models.OPENAI_DEFAULT


# The LLM judge is FIXED and SEPARATE from the model under test: a model grading its own
# answers inflates its score, and a judge that changes with the subject makes runs
# incomparable. Default: models.JUDGE_MODEL (no comparison subject uses it), at
# temperature 0 where the model accepts sampling parameters. Override with
# --judge-provider / --judge-model; a judge identical to the model under test is refused
# (see resolve_judge).
DEFAULT_JUDGE_PROVIDER = models.JUDGE_PROVIDER
DEFAULT_JUDGE_MODEL = models.JUDGE_MODEL
JUDGE_TEMPERATURE = 0.0


def resolve_judge(
    subject_provider: str,
    subject_model: str | None,
    judge_provider: str | None = None,
    judge_model: str | None = None,
) -> tuple[str, str] | None:
    """Return the (provider, model) of the judge, or None for fully offline mock runs.

    Raises ValueError when the judge would be the model under test (self-grading).
    """
    if subject_provider == "mock" and judge_provider in (None, "mock"):
        return None
    jp = judge_provider or DEFAULT_JUDGE_PROVIDER
    jm = judge_model or (DEFAULT_JUDGE_MODEL if jp == DEFAULT_JUDGE_PROVIDER else default_model(jp))
    subject = (subject_provider, subject_model or default_model(subject_provider))
    if (jp, jm) == subject:
        raise ValueError(
            f"LLM judge {jp}/{jm} is the model under test; pick a separate judge "
            "with --judge-provider/--judge-model"
        )
    return jp, jm


def build_judge_provider(provider: str, model: str) -> LLMProvider:
    """Construct the judge: the given provider/model at temperature 0 (deterministic).

    For models that reject sampling parameters the provider drops the temperature.
    """
    return get_provider(provider, model, temperature=JUDGE_TEMPERATURE)


def default_model(provider: str) -> str:
    """The model id a provider uses when none is given — recorded so runs stay reproducible."""
    return {
        "anthropic": DEFAULT_ANTHROPIC_MODEL,
        "openai": DEFAULT_OPENAI_MODEL,
        "mock": "mock",
    }.get(provider, "")


class AnthropicProvider:
    """Anthropic Claude provider using the tool-use API."""

    kind: str = "anthropic"
    # Every model in models.py accepts image content blocks.
    vision_capable: bool = True

    def __init__(
        self,
        model: str = DEFAULT_ANTHROPIC_MODEL,
        # Current models think by default and thinking spends output tokens; 4096 truncated.
        max_tokens: int = 16000,
        temperature: float | None = None,
    ) -> None:
        try:
            import anthropic  # noqa: F401
        except ImportError as e:  # pragma: no cover — install hint
            raise RuntimeError(
                "anthropic SDK not installed. Install with 'pip install \".[eval]\"'."
            ) from e
        self.model = model
        self.max_tokens = max_tokens
        self.temperature = temperature

    async def chat_with_tools(
        self,
        messages: list[dict[str, Any]],
        tools: list[ToolDescriptor],
    ) -> ProviderResponse:
        import anthropic

        client = anthropic.AsyncAnthropic()
        tool_specs = [
            {"name": t.name, "description": t.description, "input_schema": t.input_schema}
            for t in tools
        ]
        extra: dict[str, Any] = {}
        # The 5.x family rejects explicit sampling params with a 400 — omit, don't send.
        if self.temperature is not None and models.accepts_sampling_params(self.model):
            extra["temperature"] = self.temperature
        resp = await client.messages.create(
            model=self.model,
            max_tokens=self.max_tokens,
            messages=messages,  # type: ignore[arg-type]
            tools=tool_specs,  # type: ignore[arg-type]
            **extra,
        )
        text_parts: list[str] = []
        tool_calls: list[ToolCall] = []
        for block in resp.content:
            if block.type == "text":
                text_parts.append(block.text)
            elif block.type == "tool_use":
                tool_calls.append(
                    ToolCall(id=block.id, name=block.name, arguments=dict(block.input))
                )
        return ProviderResponse(text="".join(text_parts), tool_calls=tool_calls)


# -----------------------------------------------------------------------------
# OpenAI provider
# -----------------------------------------------------------------------------


class OpenAIProvider:
    """OpenAI Responses API tool-calling provider."""

    kind: str = "openai"
    # The OpenAI default (models.OPENAI_DEFAULT) accepts image content blocks.
    vision_capable: bool = True

    def __init__(
        self,
        model: str = DEFAULT_OPENAI_MODEL,
        max_tokens: int = 4096,
        temperature: float | None = None,
    ) -> None:
        try:
            import openai  # noqa: F401
        except ImportError as e:  # pragma: no cover — install hint
            raise RuntimeError(
                "openai SDK not installed. Install with 'pip install \".[eval]\"'."
            ) from e
        self.model = model
        self.max_tokens = max_tokens
        self.temperature = temperature

    async def chat_with_tools(
        self,
        messages: list[dict[str, Any]],
        tools: list[ToolDescriptor],
    ) -> ProviderResponse:
        import openai

        client = openai.AsyncOpenAI()
        tool_specs = [
            {
                "type": "function",
                "function": {
                    "name": t.name,
                    "description": t.description,
                    "parameters": t.input_schema,
                },
            }
            for t in tools
        ]
        extra: dict[str, Any] = {}
        if self.temperature is not None:
            extra["temperature"] = self.temperature
        if tool_specs:
            extra["tools"] = tool_specs
        resp = await client.chat.completions.create(
            model=self.model,
            max_tokens=self.max_tokens,
            messages=messages,  # type: ignore[arg-type]
            **extra,
        )
        choice = resp.choices[0].message
        tool_calls: list[ToolCall] = []
        if choice.tool_calls:
            import json as _json

            for tc in choice.tool_calls:
                tool_calls.append(
                    ToolCall(
                        id=tc.id,
                        name=tc.function.name,  # type: ignore[union-attr]
                        arguments=_json.loads(tc.function.arguments or "{}"),  # type: ignore[union-attr]
                    )
                )
        return ProviderResponse(text=choice.content or "", tool_calls=tool_calls)


def get_provider(
    name: str, model: str | None = None, *, temperature: float | None = None
) -> LLMProvider:
    """Factory: resolve a provider name + model to an LLMProvider instance.

    ``model`` may be ``None`` — each provider then uses its own default, so a caller that
    selects ``openai`` without a model does not inherit an Anthropic model id (and vice versa).
    ``mock`` is supported for offline runs (echoes the case input back as text).
    """
    if name == "mock":
        return make_one_shot_mock(text="[mock response]")
    if name == "anthropic":
        return AnthropicProvider(model=model or DEFAULT_ANTHROPIC_MODEL, temperature=temperature)
    if name == "openai":
        return OpenAIProvider(model=model or DEFAULT_OPENAI_MODEL, temperature=temperature)
    raise ValueError(f"Unknown provider {name!r}. Supported: anthropic, openai, mock.")


# -----------------------------------------------------------------------------
# Tool-result content construction (with optional images)
# -----------------------------------------------------------------------------


# How a UI tool's render reaches the model under test:
#   off         — text only (no render)
#   text+image  — the tool's text result plus the viewer PNG
#   image-only  — just the viewer PNG; the tool's text is withheld (the visual arm)
VISION_OFF = "off"
VISION_TEXT_IMAGE = "text+image"
VISION_IMAGE_ONLY = "image-only"
VISION_MODES = (VISION_OFF, VISION_TEXT_IMAGE, VISION_IMAGE_ONLY)

# Sent in place of the withheld text when an image-only render could not be produced, so
# the visual arm never silently falls back to the text it is meant to exclude.
IMAGE_UNAVAILABLE_TEXT = "[viewer render unavailable; tool text withheld in image-only mode]"
# What a non-vision (mock) provider sees for an image-only result.
IMAGE_ONLY_PLACEHOLDER = "[viewer image]"


def build_tool_result_content(
    text: str,
    image_bytes: bytes | None,
    provider_kind: str,
    *,
    image_only: bool = False,
) -> Any:
    """Build the right content-block shape for a tool_result, per provider.

    - Anthropic accepts a list of ``image`` and ``text`` blocks in tool_result.
    - OpenAI's chat-completions tool_result is a plain string; image content
      must be attached on a subsequent user turn. The eval runner handles that
      split — this helper produces an OpenAI-shaped image+text bundle that the
      runner can serialize appropriately.
    - Mock providers just receive plain text (images are ignored).

    When ``image_bytes`` is None or the provider doesn't support vision, the
    return value is the original ``text`` string.

    ``image_only=True`` (the ``image-only`` vision mode) withholds ``text``: the result is
    the image block alone, or :data:`IMAGE_UNAVAILABLE_TEXT` when there is no image —
    never the tool text the visual arm is meant to exclude.
    """
    if image_bytes is None:
        return IMAGE_UNAVAILABLE_TEXT if image_only else text

    import base64

    encoded = base64.b64encode(image_bytes).decode("ascii")

    image_block: dict[str, Any]
    if provider_kind == "anthropic":
        image_block = {
            "type": "image",
            "source": {
                "type": "base64",
                "media_type": "image/png",
                "data": encoded,
            },
        }
    elif provider_kind == "openai":
        image_block = {
            "type": "image_url",
            "image_url": {"url": f"data:image/png;base64,{encoded}"},
        }
    else:
        # Mock / unknown: no image support — plain text (or a placeholder in image-only).
        return IMAGE_ONLY_PLACEHOLDER if image_only else text
    if image_only:
        return [image_block]
    return [image_block, {"type": "text", "text": text}]
