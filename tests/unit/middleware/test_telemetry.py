"""Unit tests for tool-use telemetry middleware."""

from __future__ import annotations

import json
from pathlib import Path

import pytest

from bamcp.middleware import telemetry


@pytest.fixture(autouse=True)
def reset_telemetry_state():
    """Reset module-global telemetry state around every test."""
    telemetry.reset_telemetry()
    yield
    telemetry.reset_telemetry()


@pytest.fixture
def jsonl_path(tmp_path: Path) -> Path:
    return tmp_path / "events.jsonl"


def _read_events(p: Path) -> list[dict]:
    if not p.exists():
        return []
    return [json.loads(line) for line in p.read_text().splitlines() if line]


class TestSanitization:
    @pytest.mark.unit
    def test_file_path_is_hashed(self):
        out = telemetry.sanitize_args({"file_path": "/data/sample.bam", "region": "chr1:1-100"})
        assert out["region"] == "chr1:1-100"
        assert out["file_path"] != "/data/sample.bam"
        assert out["file_path"].endswith(":sample.bam")
        prefix = out["file_path"].split(":", 1)[0]
        assert len(prefix) == 12  # sha256[:12]

    @pytest.mark.unit
    def test_reference_path_is_hashed(self):
        out = telemetry.sanitize_args({"reference": "/refs/hg38.fa"})
        assert out["reference"].endswith(":hg38.fa")

    @pytest.mark.unit
    def test_auth_keys_are_dropped(self):
        out = telemetry.sanitize_args(
            {"token": "secret", "api_key": "key", "auth": "x", "region": "chr1:1-100"}
        )
        assert "token" not in out
        assert "api_key" not in out
        assert "auth" not in out
        assert out == {"region": "chr1:1-100"}

    @pytest.mark.unit
    def test_url_path_is_hashed(self):
        out = telemetry.sanitize_args({"file_path": "https://example.com/data/x.bam"})
        assert out["file_path"].endswith(":x.bam")

    @pytest.mark.unit
    def test_non_dict_returns_empty(self):
        assert telemetry.sanitize_args("not a dict") == {}  # type: ignore[arg-type]
        assert telemetry.sanitize_args(None) == {}  # type: ignore[arg-type]


class TestResultSummary:
    @pytest.mark.unit
    def test_extracts_count_field(self):
        result = {"content": [{"type": "text", "text": json.dumps({"count": 5, "variants": []})}]}
        s = telemetry.summarize_result(result)
        assert s["count"] == 5
        assert s["variants_len"] == 0

    @pytest.mark.unit
    def test_marks_ui_meta(self):
        result = {
            "content": [{"type": "text", "text": "summary"}],
            "_meta": {"ui/resourceUri": "ui://x"},
        }
        s = telemetry.summarize_result(result)
        assert s["has_ui"] is True

    @pytest.mark.unit
    def test_non_json_text_preview(self):
        result = {"content": [{"type": "text", "text": "Region chr1: 5 reads"}]}
        s = telemetry.summarize_result(result)
        assert "text_preview" in s
        assert "Region chr1" in s["text_preview"]

    @pytest.mark.unit
    def test_truncates_large_summaries(self):
        # Force a result that summarizes to a huge text_preview.
        big = "x" * 5000
        result = {"content": [{"type": "text", "text": big}]}
        s = telemetry.summarize_result(result)
        # text_preview is capped at 200 chars so size will be small; the cap
        # path is exercised separately when many fields are present.
        assert len(s["text_preview"]) == 200

    @pytest.mark.unit
    def test_truncation_triggers_on_oversized_payload(self, monkeypatch):
        # Pretend the summarizer produced a huge dict by raising the threshold.
        monkeypatch.setattr(telemetry, "_RESULT_SUMMARY_MAX_BYTES", 10)
        result = {"content": [{"type": "text", "text": json.dumps({"count": 1234567890})}]}
        s = telemetry.summarize_result(result)
        assert s.get("truncated") is True


class TestDecoratorDisabledByDefault:
    @pytest.mark.unit
    async def test_passes_through_when_disabled(self, jsonl_path):
        # Telemetry not configured — decorator must be a no-op.
        @telemetry.telemetry_wrapper("test_tool")
        async def handler(args, config):
            return {"content": [{"type": "text", "text": "ok"}]}

        result = await handler({"file_path": "/data/x.bam"}, None)
        assert result["content"][0]["text"] == "ok"
        # No file written.
        assert not jsonl_path.exists()


class TestDecoratorEnabled:
    @pytest.mark.unit
    async def test_writes_jsonl_event(self, jsonl_path):
        telemetry.configure_telemetry(enabled=True, jsonl_path=str(jsonl_path))

        @telemetry.telemetry_wrapper("get_variants")
        async def handler(args, config):
            return {
                "content": [
                    {"type": "text", "text": json.dumps({"count": 3, "variants": [1, 2, 3]})}
                ]
            }

        await handler({"file_path": "/data/x.bam", "region": "chr1:1-100"}, None)
        events = _read_events(jsonl_path)
        assert len(events) == 1
        evt = events[0]
        assert evt["tool_name"] == "get_variants"
        assert evt["ok"] is True
        assert evt["error"] is None
        assert evt["session_id"] == telemetry.get_session_id()
        assert evt["args"]["region"] == "chr1:1-100"
        assert evt["args"]["file_path"].endswith(":x.bam")
        assert evt["result_summary"]["count"] == 3
        assert evt["result_summary"]["variants_len"] == 3
        assert evt["duration_ms"] >= 0
        assert evt["trace_id"] is None
        assert evt["span_id"] is None

    @pytest.mark.unit
    async def test_captures_exception(self, jsonl_path):
        telemetry.configure_telemetry(enabled=True, jsonl_path=str(jsonl_path))

        @telemetry.telemetry_wrapper("get_coverage")
        async def handler(args, config):
            raise ValueError("boom")

        with pytest.raises(ValueError, match="boom"):
            await handler({"region": "chr1:1-100"}, None)

        events = _read_events(jsonl_path)
        assert len(events) == 1
        assert events[0]["ok"] is False
        assert events[0]["error"] == "ValueError: boom"

    @pytest.mark.unit
    async def test_multiple_calls_append(self, jsonl_path):
        telemetry.configure_telemetry(enabled=True, jsonl_path=str(jsonl_path))

        @telemetry.telemetry_wrapper("list_contigs")
        async def handler(args, config):
            return {"content": [{"type": "text", "text": "{}"}]}

        await handler({}, None)
        await handler({}, None)
        await handler({}, None)
        assert len(_read_events(jsonl_path)) == 3

    @pytest.mark.unit
    async def test_no_jsonl_path_skips_file_write(self, tmp_path):
        telemetry.configure_telemetry(enabled=True, jsonl_path="")

        @telemetry.telemetry_wrapper("list_contigs")
        async def handler(args, config):
            return {"content": [{"type": "text", "text": "{}"}]}

        # Should not raise even with no path configured.
        await handler({}, None)
        # Nothing in tmp_path either.
        assert not any(tmp_path.iterdir())

    @pytest.mark.unit
    async def test_session_id_is_stable_across_calls(self, jsonl_path):
        telemetry.configure_telemetry(
            enabled=True, jsonl_path=str(jsonl_path), session_id="fixed-session"
        )

        @telemetry.telemetry_wrapper("list_contigs")
        async def handler(args, config):
            return {"content": [{"type": "text", "text": "{}"}]}

        await handler({}, None)
        await handler({}, None)
        events = _read_events(jsonl_path)
        assert all(e["session_id"] == "fixed-session" for e in events)


class TestOtelIntegration:
    @pytest.mark.unit
    async def test_otel_span_emitted_when_configured(self, jsonl_path, monkeypatch):
        """When an OTel tracer is configured, span_id and trace_id are captured."""
        try:
            from opentelemetry.sdk.trace import TracerProvider
            from opentelemetry.sdk.trace.export import SimpleSpanProcessor
            from opentelemetry.sdk.trace.export.in_memory_span_exporter import (
                InMemorySpanExporter,
            )
        except ImportError:
            pytest.skip("opentelemetry SDK not installed")

        # Build an in-memory tracer that captures spans.
        exporter = InMemorySpanExporter()
        provider = TracerProvider()
        provider.add_span_processor(SimpleSpanProcessor(exporter))
        # Set the module-global tracer directly (bypass OTLP exporter init).
        telemetry.configure_telemetry(enabled=True, jsonl_path=str(jsonl_path))
        monkeypatch.setattr(telemetry, "_otel_tracer", provider.get_tracer("bamcp-test"))

        @telemetry.telemetry_wrapper("visualize_region")
        async def handler(args, config):
            return {"content": [{"type": "text", "text": "ok"}]}

        await handler({"region": "chr1:1-100"}, None)

        events = _read_events(jsonl_path)
        assert len(events) == 1
        assert events[0]["trace_id"] is not None
        assert events[0]["span_id"] is not None
        # Span captured by the exporter:
        spans = exporter.get_finished_spans()
        assert len(spans) == 1
        assert spans[0].name == "bamcp.visualize_region"

    @pytest.mark.unit
    async def test_otel_records_error_attributes(self, jsonl_path, monkeypatch):
        try:
            from opentelemetry.sdk.trace import TracerProvider
            from opentelemetry.sdk.trace.export import SimpleSpanProcessor
            from opentelemetry.sdk.trace.export.in_memory_span_exporter import (
                InMemorySpanExporter,
            )
        except ImportError:
            pytest.skip("opentelemetry SDK not installed")

        exporter = InMemorySpanExporter()
        provider = TracerProvider()
        provider.add_span_processor(SimpleSpanProcessor(exporter))
        telemetry.configure_telemetry(enabled=True, jsonl_path=str(jsonl_path))
        monkeypatch.setattr(telemetry, "_otel_tracer", provider.get_tracer("bamcp-test"))

        @telemetry.telemetry_wrapper("lookup_clinvar")
        async def handler(args, config):
            raise RuntimeError("nope")

        with pytest.raises(RuntimeError):
            await handler({}, None)

        spans = exporter.get_finished_spans()
        assert len(spans) == 1
        attrs = dict(spans[0].attributes or {})
        assert attrs.get("error") is True
        assert "nope" in str(attrs.get("error.message", ""))


class TestConfigure:
    @pytest.mark.unit
    def test_disabled_is_default(self):
        assert telemetry.is_enabled() is False
        assert telemetry.get_session_id() == ""

    @pytest.mark.unit
    def test_enable_sets_session(self):
        telemetry.configure_telemetry(enabled=True, session_id="abc")
        assert telemetry.is_enabled() is True
        assert telemetry.get_session_id() == "abc"

    @pytest.mark.unit
    def test_creates_parent_dir(self, tmp_path):
        target = tmp_path / "nested" / "deep" / "events.jsonl"
        telemetry.configure_telemetry(enabled=True, jsonl_path=str(target))
        assert target.parent.is_dir()
