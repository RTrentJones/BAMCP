"""Unit tests for the visual-eval integration points.

Covers:
- Router surfacing ``_meta.ui/init`` into RouterResult.ui_payload.
- Provider vision_capable / kind attributes.
- build_tool_result_content shape per provider.
- Runner capturing screenshots and embedding image content blocks.
- Skip path for vision_required cases when vision is unavailable.
"""

from __future__ import annotations

import json
from pathlib import Path

import pytest

from bamcp.config import BAMCPConfig
from bamcp.eval import runner as runner_mod
from bamcp.eval.providers import (
    AnthropicProvider,
    MockProvider,
    OpenAIProvider,
    ProviderResponse,
    ToolCall,
    build_tool_result_content,
)
from bamcp.eval.router import InProcessRouter
from bamcp.eval.runner import run_case
from bamcp.eval.schema import EvalCase


@pytest.fixture
def cfg(ref_fasta_path):
    return BAMCPConfig(reference=ref_fasta_path)


# ---------------------------------------------------------------------------
# Router surfacing _meta.ui/init
# ---------------------------------------------------------------------------


class TestRouterUiPayload:
    @pytest.mark.unit
    async def test_visualize_region_exposes_ui_payload(self, small_bam_path, cfg):
        router = InProcessRouter(config=cfg)
        result = await router.call(
            "visualize_region", {"file_path": small_bam_path, "region": "chr1:90-200"}
        )
        assert result.ok is True
        assert result.ui_payload is not None
        assert result.ui_payload.get("contig") == "chr1"
        assert "reads" in result.ui_payload

    @pytest.mark.unit
    async def test_non_ui_tool_has_no_payload(self, small_bam_path, cfg):
        router = InProcessRouter(config=cfg)
        result = await router.call(
            "get_coverage", {"file_path": small_bam_path, "region": "chr1:90-200"}
        )
        assert result.ok is True
        assert result.ui_payload is None


# ---------------------------------------------------------------------------
# Provider vision_capable / kind
# ---------------------------------------------------------------------------


class TestProviderVisionAttrs:
    @pytest.mark.unit
    def test_anthropic_is_vision_capable(self):
        p = AnthropicProvider()
        assert p.vision_capable is True
        assert p.kind == "anthropic"

    @pytest.mark.unit
    def test_openai_is_vision_capable(self):
        p = OpenAIProvider()
        assert p.vision_capable is True
        assert p.kind == "openai"

    @pytest.mark.unit
    def test_mock_defaults_to_text_only(self):
        p = MockProvider(lambda m, t: ProviderResponse(text=""))
        assert p.vision_capable is False
        assert p.kind == "mock"

    @pytest.mark.unit
    def test_mock_can_opt_into_vision(self):
        p = MockProvider(lambda m, t: ProviderResponse(text=""), vision_capable=True)
        assert p.vision_capable is True


# ---------------------------------------------------------------------------
# build_tool_result_content
# ---------------------------------------------------------------------------


class TestBuildToolResultContent:
    @pytest.mark.unit
    def test_no_image_returns_plain_text(self):
        out = build_tool_result_content("hello", None, "anthropic")
        assert out == "hello"

    @pytest.mark.unit
    def test_anthropic_shape(self):
        png = b"\x89PNG\r\n\x1a\nbytes"
        out = build_tool_result_content("summary", png, "anthropic")
        assert isinstance(out, list)
        assert out[0]["type"] == "image"
        assert out[0]["source"]["type"] == "base64"
        assert out[0]["source"]["media_type"] == "image/png"
        assert out[1] == {"type": "text", "text": "summary"}

    @pytest.mark.unit
    def test_openai_shape(self):
        png = b"\x89PNG\r\n\x1a\nbytes"
        out = build_tool_result_content("summary", png, "openai")
        assert isinstance(out, list)
        assert out[0]["type"] == "image_url"
        assert out[0]["image_url"]["url"].startswith("data:image/png;base64,")
        assert out[1] == {"type": "text", "text": "summary"}

    @pytest.mark.unit
    def test_unknown_kind_falls_back_to_text(self):
        png = b"\x89PNG\r\n\x1a\nbytes"
        out = build_tool_result_content("summary", png, "some-unknown-kind")
        assert out == "summary"


# ---------------------------------------------------------------------------
# Runner: vision-required skip path
# ---------------------------------------------------------------------------


class TestVisionRequiredSkip:
    @pytest.mark.unit
    async def test_skipped_when_provider_lacks_vision(self, tmp_path: Path, cfg):
        case = EvalCase(
            name="vision-only",
            input="show me",
            expected="ok",
            rendering_modes=["expanded"],
            vision_required=True,
        )
        provider = MockProvider(lambda m, t: ProviderResponse(text="should not run"))
        result = await run_case(
            case,
            provider,
            InProcessRouter(config=cfg),
            judge=None,
            output_dir=tmp_path,
            rendering_mode="expanded",
        )
        assert result.skipped is True
        assert "vision-capable" in (result.skip_reason or "")
        assert result.verdict.passed is False
        assert result.tool_calls == []

    @pytest.mark.unit
    async def test_skipped_when_no_vision_flag_set(self, tmp_path: Path, cfg):
        case = EvalCase(
            name="vision-only",
            input="show me",
            expected="ok",
            rendering_modes=["expanded"],
            vision_required=True,
        )
        provider = MockProvider(lambda m, t: ProviderResponse(text="ok"), vision_capable=True)
        result = await run_case(
            case,
            provider,
            InProcessRouter(config=cfg),
            judge=None,
            output_dir=tmp_path,
            rendering_mode="expanded",
            no_vision=True,
        )
        assert result.skipped is True
        assert "no-vision" in (result.skip_reason or "").lower() or "disabled" in (
            result.skip_reason or ""
        )

    @pytest.mark.unit
    async def test_skipped_when_no_rendering_mode(self, tmp_path: Path, cfg):
        case = EvalCase(
            name="vision-only",
            input="show me",
            expected="ok",
            vision_required=True,
        )
        provider = MockProvider(lambda m, t: ProviderResponse(text="ok"), vision_capable=True)
        # No rendering_mode passed in → vision can't activate.
        result = await run_case(
            case,
            provider,
            InProcessRouter(config=cfg),
            judge=None,
            output_dir=tmp_path,
            rendering_mode=None,
        )
        assert result.skipped is True


# ---------------------------------------------------------------------------
# Runner: actual image-in-tool-result flow via FakeRenderer
# ---------------------------------------------------------------------------


class TestRunnerVisionFlow:
    @pytest.mark.unit
    async def test_vision_capture_attached_to_tool_result(
        self, tmp_path: Path, small_bam_path, cfg, fake_renderer, monkeypatch
    ):
        """When the LLM calls visualize_region, the runner should capture a
        screenshot, pass it to the LLM, and record vision_used=True."""

        # Stub _make_renderer to return our FakeRenderer.
        monkeypatch.setattr(runner_mod, "_make_renderer", lambda: fake_renderer)

        # Script the provider: first turn calls visualize_region, second turn
        # produces a plain-text response. Record the messages the second turn
        # receives so we can assert an image block reached the LLM.
        seen_messages: list = []

        def script(messages, tools):
            seen_messages.append(messages)
            if len(seen_messages) == 1:
                return ProviderResponse(
                    text="visualizing",
                    tool_calls=[
                        ToolCall(
                            id="c1",
                            name="visualize_region",
                            arguments={
                                "file_path": small_bam_path,
                                "region": "chr1:90-200",
                            },
                        )
                    ],
                )
            return ProviderResponse(text="done looking")

        provider = MockProvider(script, vision_capable=True)
        # MockProvider defaults to kind="mock" — switch to anthropic so the
        # image is actually built into the tool_result content list.
        provider.kind = "anthropic"

        case = EvalCase(
            name="vision smoke",
            input="render the region",
            expected="any",
            rendering_modes=["expanded"],
            vision_required=False,
        )
        result = await run_case(
            case,
            provider,
            InProcessRouter(config=cfg),
            judge=None,
            output_dir=tmp_path,
            rendering_mode="expanded",
        )

        # Vision was used, one screenshot persisted, fake renderer entered/exited
        assert result.vision_used is True
        assert len(result.screenshots) == 1
        screenshot = Path(result.screenshots[0])
        assert screenshot.exists()
        assert fake_renderer.entered and fake_renderer.exited

        # The second LLM turn saw a tool_result with image content.
        second_turn_messages = seen_messages[1]
        last_user = second_turn_messages[-1]
        assert last_user["role"] == "user"
        tool_result = last_user["content"][0]
        assert tool_result["type"] == "tool_result"
        # Image content block present.
        image_block = tool_result["content"][0]
        assert image_block["type"] == "image"
        assert image_block["source"]["type"] == "base64"

    @pytest.mark.unit
    async def test_max_captures_per_case_respected(
        self, tmp_path: Path, small_bam_path, cfg, fake_renderer, monkeypatch
    ):
        """The runner caps captures per case at MAX_CAPTURES_PER_CASE."""

        monkeypatch.setattr(runner_mod, "_make_renderer", lambda: fake_renderer)
        # Tighten the cap so we only need a few turns to hit it.
        monkeypatch.setattr(runner_mod, "MAX_CAPTURES_PER_CASE", 2)
        # And widen MAX_TURNS so the loop can keep producing tool calls.
        monkeypatch.setattr(runner_mod, "MAX_TURNS", 10)

        turn = {"n": 0}

        def script(messages, tools):
            turn["n"] += 1
            if turn["n"] <= 4:
                return ProviderResponse(
                    text=f"turn {turn['n']}",
                    tool_calls=[
                        ToolCall(
                            id=f"c{turn['n']}",
                            name="visualize_region",
                            arguments={
                                "file_path": small_bam_path,
                                "region": "chr1:90-200",
                            },
                        )
                    ],
                )
            return ProviderResponse(text="done")

        provider = MockProvider(script, vision_capable=True)
        provider.kind = "anthropic"

        case = EvalCase(
            name="cap test",
            input="x",
            expected="x",
            rendering_modes=["expanded"],
        )
        result = await run_case(
            case,
            provider,
            InProcessRouter(config=cfg),
            judge=None,
            output_dir=tmp_path,
            rendering_mode="expanded",
        )
        # Cap respected even though the LLM kept asking for more visuals.
        assert len(result.screenshots) == 2

    @pytest.mark.unit
    async def test_capture_failure_degrades_gracefully(
        self, tmp_path: Path, small_bam_path, cfg, fake_renderer, monkeypatch
    ):
        """A renderer crash shouldn't abort the case — the tool_result just
        falls back to text-only."""

        monkeypatch.setattr(runner_mod, "_make_renderer", lambda: fake_renderer)
        fake_renderer.make_fail_next()

        def script(messages, tools):
            if not hasattr(script, "called"):
                script.called = True
                return ProviderResponse(
                    text="trying",
                    tool_calls=[
                        ToolCall(
                            id="c1",
                            name="visualize_region",
                            arguments={
                                "file_path": small_bam_path,
                                "region": "chr1:90-200",
                            },
                        )
                    ],
                )
            return ProviderResponse(text="recovered")

        provider = MockProvider(script, vision_capable=True)
        provider.kind = "anthropic"

        result = await run_case(
            EvalCase(
                name="degrade",
                input="x",
                expected="x",
                rendering_modes=["expanded"],
            ),
            provider,
            InProcessRouter(config=cfg),
            judge=None,
            output_dir=tmp_path,
            rendering_mode="expanded",
        )
        # No screenshots, but the run still completed.
        assert result.screenshots == []
        assert result.vision_used is False
        assert result.error is None
        assert "visualize_region" in result.tool_calls


# ---------------------------------------------------------------------------
# Reports: vision_used count + thumbnails
# ---------------------------------------------------------------------------


class TestReportsWithScreenshots:
    @pytest.mark.unit
    def test_html_includes_thumbnail_when_screenshots_present(self, tmp_path: Path):
        from bamcp.eval.report import write_reports
        from bamcp.eval.schema import EvalResult, GraderVerdict, RunConfig

        shot = tmp_path / "screenshots" / "x__expanded__step0.png"
        shot.parent.mkdir(parents=True)
        shot.write_bytes(b"\x89PNG fake")

        result = EvalResult(
            case_name="x",
            category="vision",
            rendering_mode="expanded",
            response_text="ok",
            tool_calls=["visualize_region"],
            duration_ms=10.0,
            verdict=GraderVerdict(passed=True, rationale="ok", deterministic=True, score=1.0),
            screenshots=[str(shot)],
            vision_used=True,
        )
        cfg = RunConfig(
            test_cases_path=tmp_path / "cases.yaml",
            output_dir=tmp_path,
            provider="mock",
            model="m",
        )
        paths = write_reports([result], cfg, tmp_path)
        html = paths["html"].read_text()
        assert "<img class='thumb'" in html
        assert "screenshots/x__expanded__step0.png" in html
        # Vision summary line.
        assert "vision: 1 used" in html
        # The companion side-by-side report exists too.
        assert "vision_html" in paths
        vision_html = paths["vision_html"].read_text()
        assert "<img" in vision_html

    @pytest.mark.unit
    def test_html_marks_skipped_rows(self, tmp_path: Path):
        from bamcp.eval.report import write_reports
        from bamcp.eval.schema import EvalResult, GraderVerdict, RunConfig

        result = EvalResult(
            case_name="vision-only",
            category="vision",
            rendering_mode="expanded",
            response_text="",
            tool_calls=[],
            duration_ms=0.0,
            verdict=GraderVerdict(
                passed=False,
                rationale="Skipped: provider not vision-capable",
                deterministic=True,
            ),
            skipped=True,
            skip_reason="provider not vision-capable",
        )
        cfg = RunConfig(
            test_cases_path=tmp_path / "cases.yaml",
            output_dir=tmp_path,
            provider="mock",
            model="m",
        )
        paths = write_reports([result], cfg, tmp_path)
        html = paths["html"].read_text()
        assert "tr class='skipped'" in html
        assert "SKIP" in html


# ---------------------------------------------------------------------------
# Schema: vision_required parsing
# ---------------------------------------------------------------------------


class TestSchemaVisionRequired:
    @pytest.mark.unit
    def test_vision_required_parsed_from_yaml(self, tmp_path: Path):
        from bamcp.eval.schema import load_cases

        p = tmp_path / "cases.yaml"
        p.write_text(
            "- case:\n    name: vis\n    input: q\n    expected: e\n    vision_required: true\n"
        )
        cases = load_cases(p)
        assert cases[0].vision_required is True

    @pytest.mark.unit
    def test_vision_required_defaults_to_false(self, tmp_path: Path):
        from bamcp.eval.schema import load_cases

        p = tmp_path / "cases.yaml"
        p.write_text("- case:\n    name: x\n    input: q\n    expected: e\n")
        cases = load_cases(p)
        assert cases[0].vision_required is False

    @pytest.mark.unit
    def test_seed_vision_cases_file_parses(self):
        """The shipped tests/eval/vision_cases.yaml must parse cleanly."""
        from bamcp.eval.schema import load_cases

        path = Path(__file__).parent.parent.parent / "eval" / "vision_cases.yaml"
        cases = load_cases(path)
        assert len(cases) >= 5
        # Every shipped vision case must declare vision_required.
        for c in cases:
            assert c.vision_required is True
            assert c.rendering_modes, f"{c.name}: should declare rendering_modes"


# ---------------------------------------------------------------------------
# CLI: --no-vision and --vision-screenshot-dir flags
# ---------------------------------------------------------------------------


class TestCliVisionFlags:
    @pytest.mark.unit
    def test_no_vision_defaults_to_false(self):
        from bamcp.eval.cli import parse_args

        cfg = parse_args(["--provider", "mock"])
        assert cfg.no_vision is False
        assert cfg.vision_screenshot_dir is None

    @pytest.mark.unit
    def test_no_vision_flag_sets_config(self):
        from bamcp.eval.cli import parse_args

        cfg = parse_args(["--provider", "mock", "--no-vision"])
        assert cfg.no_vision is True

    @pytest.mark.unit
    def test_screenshot_dir_flag_sets_path(self, tmp_path: Path):
        from bamcp.eval.cli import parse_args

        cfg = parse_args(["--provider", "mock", "--vision-screenshot-dir", str(tmp_path)])
        assert cfg.vision_screenshot_dir == tmp_path


# ---------------------------------------------------------------------------
# Sanity check: results JSON round-trips the new fields
# ---------------------------------------------------------------------------


@pytest.mark.unit
def test_results_json_includes_new_fields(tmp_path: Path):
    from bamcp.eval.report import write_reports
    from bamcp.eval.schema import EvalResult, GraderVerdict, RunConfig

    result = EvalResult(
        case_name="x",
        category="vision",
        rendering_mode="expanded",
        response_text="ok",
        tool_calls=["visualize_region"],
        duration_ms=10.0,
        verdict=GraderVerdict(passed=True, rationale="ok", deterministic=True, score=1.0),
        screenshots=[str(tmp_path / "x.png")],
        vision_used=True,
    )
    cfg = RunConfig(
        test_cases_path=tmp_path / "cases.yaml",
        output_dir=tmp_path,
        provider="mock",
        model="m",
    )
    paths = write_reports([result], cfg, tmp_path)
    parsed = json.loads(paths["json"].read_text())
    assert parsed[0]["vision_used"] is True
    assert parsed[0]["screenshots"] == [str(tmp_path / "x.png")]
