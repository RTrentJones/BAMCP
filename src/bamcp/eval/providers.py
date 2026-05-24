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


class AnthropicProvider:
    """Anthropic Claude provider using the tool-use API."""

    kind: str = "anthropic"
    # All current claude-opus-4-* and claude-sonnet-4-* models accept vision.
    vision_capable: bool = True

    def __init__(self, model: str = "claude-opus-4-7", max_tokens: int = 4096) -> None:
        try:
            import anthropic  # noqa: F401
        except ImportError as e:  # pragma: no cover — install hint
            raise RuntimeError(
                "anthropic SDK not installed. Install with 'pip install \".[eval]\"'."
            ) from e
        self.model = model
        self.max_tokens = max_tokens

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
        resp = await client.messages.create(
            model=self.model,
            max_tokens=self.max_tokens,
            messages=messages,  # type: ignore[arg-type]
            tools=tool_specs,  # type: ignore[arg-type]
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
    # gpt-4o family accepts image content blocks.
    vision_capable: bool = True

    def __init__(self, model: str = "gpt-4o-mini", max_tokens: int = 4096) -> None:
        try:
            import openai  # noqa: F401
        except ImportError as e:  # pragma: no cover — install hint
            raise RuntimeError(
                "openai SDK not installed. Install with 'pip install \".[eval]\"'."
            ) from e
        self.model = model
        self.max_tokens = max_tokens

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
        resp = await client.chat.completions.create(
            model=self.model,
            max_tokens=self.max_tokens,
            messages=messages,  # type: ignore[arg-type]
            tools=tool_specs,  # type: ignore[arg-type]
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


def get_provider(name: str, model: str) -> LLMProvider:
    """Factory: resolve a provider name + model to an LLMProvider instance.

    ``mock`` is supported for offline runs (echoes the case input back as text).
    """
    if name == "mock":
        return make_one_shot_mock(text="[mock response]")
    if name == "anthropic":
        return AnthropicProvider(model=model)
    if name == "openai":
        return OpenAIProvider(model=model)
    raise ValueError(f"Unknown provider {name!r}. Supported: anthropic, openai, mock.")


# -----------------------------------------------------------------------------
# Tool-result content construction (with optional images)
# -----------------------------------------------------------------------------


def build_tool_result_content(
    text: str,
    image_bytes: bytes | None,
    provider_kind: str,
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
    """
    if image_bytes is None:
        return text

    import base64

    encoded = base64.b64encode(image_bytes).decode("ascii")

    if provider_kind == "anthropic":
        return [
            {
                "type": "image",
                "source": {
                    "type": "base64",
                    "media_type": "image/png",
                    "data": encoded,
                },
            },
            {"type": "text", "text": text},
        ]
    if provider_kind == "openai":
        return [
            {
                "type": "image_url",
                "image_url": {"url": f"data:image/png;base64,{encoded}"},
            },
            {"type": "text", "text": text},
        ]
    # Mock / unknown: fall back to plain text.
    return text
