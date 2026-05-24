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
    far so it can implement multi-turn behavior.
    """

    def __init__(
        self,
        script: Callable[[list[dict[str, Any]], list[ToolDescriptor]], ProviderResponse]
        | Callable[[list[dict[str, Any]], list[ToolDescriptor]], Awaitable[ProviderResponse]],
    ) -> None:
        self._script = script

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
