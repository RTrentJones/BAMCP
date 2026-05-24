"""Unit tests for the render-viewer CLI helpers (no real Playwright)."""

from __future__ import annotations

import time
from pathlib import Path

import pytest

from bamcp.eval.render_cli import (
    _load_ui_init_json,
    _parse_active_variant,
    _parse_args,
    _resolve_region,
)


class TestArgParsing:
    @pytest.mark.unit
    def test_bam_and_json_are_mutually_exclusive(self):
        with pytest.raises(SystemExit):
            _parse_args(["--bam", "x.bam", "--ui-init-json", "y.json"])

    @pytest.mark.unit
    def test_requires_source(self):
        with pytest.raises(SystemExit):
            _parse_args(["--mode", "expanded"])

    @pytest.mark.unit
    def test_mode_validated(self):
        with pytest.raises(SystemExit):
            _parse_args(["--bam", "x.bam", "--mode", "not-a-real-mode"])

    @pytest.mark.unit
    def test_defaults(self):
        args = _parse_args(["--bam", "x.bam"])
        assert args.mode == "expanded"
        assert args.width == 1024
        assert args.height == 800
        assert args.out == ".render-dev/last.png"
        assert args.print_path_only is False


class TestActiveVariantParsing:
    @pytest.mark.unit
    def test_parses_valid_spec(self):
        pos, alt = _parse_active_variant("chr1:1050:A>T")
        assert pos == 1050
        assert alt == "T"

    @pytest.mark.unit
    def test_uppercases_alt(self):
        _pos, alt = _parse_active_variant("chr1:50:a>t")
        assert alt == "T"

    @pytest.mark.unit
    def test_none_returns_nones(self):
        assert _parse_active_variant(None) == (None, None)
        assert _parse_active_variant("") == (None, None)

    @pytest.mark.unit
    def test_rejects_malformed(self):
        with pytest.raises(SystemExit, match="active-variant"):
            _parse_active_variant("chr1-50-A-T")


class TestUiInitJsonLoader:
    @pytest.mark.unit
    def test_loads_dict(self, tmp_path: Path):
        p = tmp_path / "x.json"
        p.write_text('{"contig": "chr1", "start": 0, "end": 100}')
        payload = _load_ui_init_json(str(p))
        assert payload["contig"] == "chr1"

    @pytest.mark.unit
    def test_rejects_non_object(self, tmp_path: Path):
        p = tmp_path / "x.json"
        p.write_text("[1, 2, 3]")
        with pytest.raises(ValueError, match="JSON object"):
            _load_ui_init_json(str(p))


class TestResolveRegion:
    @pytest.mark.unit
    def test_explicit_region_wins(self, small_bam_path):
        assert _resolve_region(small_bam_path, "chr1:200-500") == "chr1:200-500"

    @pytest.mark.unit
    def test_default_uses_first_contig(self, small_bam_path):
        region = _resolve_region(small_bam_path, None)
        assert ":" in region
        # Should start at 1 and end somewhere reasonable
        contig, span = region.split(":")
        start, end = span.split("-")
        assert int(start) == 1
        assert 1 <= int(end) <= 1000


class TestRebuildGuard:
    @pytest.mark.unit
    def test_returns_none_when_dist_fresh(self, monkeypatch, tmp_path: Path):
        from bamcp.eval import render_cli

        static = tmp_path / "static"
        static.mkdir()
        (static / "dist").mkdir()
        dist_html = static / "dist" / "viewer.html"
        ts_file = static / "renderer.ts"
        ts_file.write_text("// source")
        # Touch dist after sources so dist is newer.
        time.sleep(0.01)
        dist_html.write_text("<html>")

        monkeypatch.setattr(render_cli, "_STATIC_DIR", static)
        assert render_cli._check_rebuild_needed() is None

    @pytest.mark.unit
    def test_warns_when_sources_newer(self, monkeypatch, tmp_path: Path):
        from bamcp.eval import render_cli

        static = tmp_path / "static"
        static.mkdir()
        (static / "dist").mkdir()
        dist_html = static / "dist" / "viewer.html"
        dist_html.write_text("<html>")
        time.sleep(0.01)
        ts_file = static / "renderer.ts"
        ts_file.write_text("// updated")

        monkeypatch.setattr(render_cli, "_STATIC_DIR", static)
        warning = render_cli._check_rebuild_needed()
        assert warning is not None
        assert "renderer.ts" in warning

    @pytest.mark.unit
    def test_warns_when_dist_missing(self, monkeypatch, tmp_path: Path):
        from bamcp.eval import render_cli

        static = tmp_path / "static"
        static.mkdir()
        monkeypatch.setattr(render_cli, "_STATIC_DIR", static)
        warning = render_cli._check_rebuild_needed()
        assert warning is not None
        assert "does not exist" in warning


class TestPrintPathOnly:
    @pytest.mark.unit
    def test_main_prints_only_path_when_flag_set(self, monkeypatch, capsys, tmp_path: Path):
        """When --print-path-only is set, stdout is exactly one line — the output path."""
        from bamcp.eval import render_cli

        # Stub out the async renderer so we don't launch Chromium.
        async def fake_capture(self, *a, **kw):
            return b"\x89PNG\r\n\x1a\n" + b"fake"

        async def fake_aenter(self):
            return self

        async def fake_aexit(self, *exc):
            return None

        monkeypatch.setattr(render_cli.ViewerRenderer, "capture", fake_capture)
        monkeypatch.setattr(render_cli.ViewerRenderer, "__aenter__", fake_aenter)
        monkeypatch.setattr(render_cli.ViewerRenderer, "__aexit__", fake_aexit)

        # And the BAM build path; pretend it returns a tiny payload.
        async def fake_build(*a, **kw):
            return {"contig": "chr1", "start": 1, "end": 10}

        monkeypatch.setattr(render_cli, "_build_payload_from_bam", fake_build)
        monkeypatch.setattr(render_cli, "_resolve_region", lambda bam, region: "chr1:1-10")

        out = tmp_path / "out.png"
        rc = render_cli.main(
            [
                "--bam",
                "fake.bam",
                "--mode",
                "expanded",
                "--out",
                str(out),
                "--print-path-only",
            ]
        )
        captured = capsys.readouterr()
        assert rc == 0
        assert captured.out.strip() == str(out)
        # No rebuild warning on stderr either (print-path-only suppresses it).
        assert "WARNING" not in captured.err

    @pytest.mark.unit
    def test_main_prints_size_summary_by_default(self, monkeypatch, capsys, tmp_path: Path):
        from bamcp.eval import render_cli

        async def fake_capture(self, *a, **kw):
            return b"\x89PNG\r\n\x1a\n" + b"x" * 1000

        async def fake_aenter(self):
            return self

        async def fake_aexit(self, *exc):
            return None

        monkeypatch.setattr(render_cli.ViewerRenderer, "capture", fake_capture)
        monkeypatch.setattr(render_cli.ViewerRenderer, "__aenter__", fake_aenter)
        monkeypatch.setattr(render_cli.ViewerRenderer, "__aexit__", fake_aexit)

        async def fake_build(*a, **kw):
            return {"contig": "chr1", "start": 1, "end": 10}

        monkeypatch.setattr(render_cli, "_build_payload_from_bam", fake_build)
        monkeypatch.setattr(render_cli, "_resolve_region", lambda bam, region: "chr1:1-10")
        # Suppress the rebuild guard so we don't pollute stderr in this test.
        monkeypatch.setattr(render_cli, "_check_rebuild_needed", lambda: None)

        out = tmp_path / "out.png"
        rc = render_cli.main(["--bam", "fake.bam", "--mode", "dv-strips", "--out", str(out)])
        captured = capsys.readouterr()
        assert rc == 0
        assert "Rendered" in captured.out
        assert "dv-strips" in captured.out
