"""Unit tests for provider abstractions (mock-only — no API calls)."""

from __future__ import annotations

import pytest

from bamcp.eval.providers import MockProvider, ProviderResponse, ToolCall, get_provider


@pytest.mark.unit
async def test_mock_provider_returns_scripted_response():
    def script(messages, tools):
        return ProviderResponse(text="hi", tool_calls=[ToolCall(id="1", name="x", arguments={})])

    provider = MockProvider(script)
    resp = await provider.chat_with_tools(messages=[], tools=[])
    assert resp.text == "hi"
    assert resp.tool_calls[0].name == "x"


@pytest.mark.unit
async def test_mock_provider_accepts_async_script():
    async def script(messages, tools):
        return ProviderResponse(text="async hi")

    provider = MockProvider(script)
    resp = await provider.chat_with_tools(messages=[], tools=[])
    assert resp.text == "async hi"


@pytest.mark.unit
def test_get_provider_mock_works():
    provider = get_provider("mock", "any-model")
    assert provider is not None


@pytest.mark.unit
def test_get_provider_rejects_unknown():
    with pytest.raises(ValueError, match="Unknown provider"):
        get_provider("not-a-provider", "x")


@pytest.mark.unit
def test_get_provider_anthropic_requires_sdk(monkeypatch):
    """Importing anthropic must succeed for AnthropicProvider construction."""
    # The SDK was pulled in by `[eval]` extras; just smoke-test construction.
    from bamcp.eval.providers import AnthropicProvider

    p = AnthropicProvider(model="claude-opus-4-7")
    assert p.model == "claude-opus-4-7"


@pytest.mark.unit
def test_get_provider_openai_construction():
    from bamcp.eval.providers import OpenAIProvider

    p = OpenAIProvider(model="gpt-4o-mini")
    assert p.model == "gpt-4o-mini"
