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


@pytest.mark.unit
def test_get_provider_none_model_uses_provider_default():
    """A None model must not leak an Anthropic id into the OpenAI provider (and vice versa)."""
    assert get_provider("openai", None).model.startswith("gpt")
    assert get_provider("anthropic", None).model.startswith("claude")
    # An explicit model is still honored.
    assert get_provider("openai", "gpt-4o").model == "gpt-4o"


# ── Fixed, separate, deterministic judge ────────────────────────────────────


@pytest.mark.unit
def test_resolve_judge_defaults_to_fixed_separate_model():
    from bamcp.eval.providers import DEFAULT_JUDGE_MODEL, DEFAULT_JUDGE_PROVIDER, resolve_judge

    for subject in ("claude-opus-4-7", "claude-opus-4-8", "claude-sonnet-4-5"):
        assert resolve_judge("anthropic", subject) == (DEFAULT_JUDGE_PROVIDER, DEFAULT_JUDGE_MODEL)
    assert resolve_judge("mock", "mock") is None


@pytest.mark.unit
def test_resolve_judge_refuses_self_grading():
    from bamcp.eval.providers import DEFAULT_JUDGE_MODEL, resolve_judge

    with pytest.raises(ValueError, match="model under test"):
        resolve_judge("anthropic", DEFAULT_JUDGE_MODEL)
    with pytest.raises(ValueError, match="model under test"):
        resolve_judge("anthropic", "claude-x", "anthropic", "claude-x")
    # An explicit override that differs is fine, and openai falls back to its default.
    assert resolve_judge("anthropic", "claude-x", "anthropic", "claude-y") == (
        "anthropic",
        "claude-y",
    )
    assert resolve_judge("anthropic", "claude-x", "openai")[0] == "openai"


@pytest.mark.unit
def test_judge_runs_at_temperature_zero():
    pytest.importorskip("anthropic")
    from bamcp.eval.providers import build_judge_provider

    judge = build_judge_provider("anthropic", "claude-y")
    assert judge.temperature == 0.0
    assert judge.model == "claude-y"


@pytest.mark.unit
async def test_anthropic_provider_forwards_temperature(monkeypatch):
    anthropic = pytest.importorskip("anthropic")
    from bamcp.eval.providers import AnthropicProvider

    seen = {}

    class _Resp:
        content = []

    class _Messages:
        async def create(self, **kwargs):
            seen.update(kwargs)
            return _Resp()

    class _Client:
        messages = _Messages()

    monkeypatch.setattr(anthropic, "AsyncAnthropic", lambda: _Client())
    await AnthropicProvider("m", temperature=0.0).chat_with_tools([], [])
    assert seen["temperature"] == 0.0
    seen.clear()
    await AnthropicProvider("m").chat_with_tools([], [])
    assert "temperature" not in seen
