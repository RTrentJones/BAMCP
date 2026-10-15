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
        assert args.height is None
        assert args.blind is False
        assert args.clip is None
        assert args.out == ".render-dev/last.png"
        assert args.print_path_only is False

    @pytest.mark.unit
    def test_height_resolution(self):
        from bamcp.eval.render_cli import _resolve_height

        assert _resolve_height(_parse_args(["--bam", "x.bam"])) == 800
        assert _resolve_height(_parse_args(["--bam", "x.bam", "--blind"])) is None
        assert _resolve_height(_parse_args(["--bam", "x.bam", "--height", "auto"])) is None
        assert _resolve_height(_parse_args(["--bam", "x.bam", "--blind", "--height", "600"])) == 600

    @pytest.mark.unit
    @pytest.mark.parametrize("bad", ["tall", "0", "-5"])
    def test_height_rejects_garbage(self, bad):
        with pytest.raises(SystemExit):
            _parse_args(["--bam", "x.bam", "--height", bad])

    @pytest.mark.unit
    def test_clip_validated(self):
        with pytest.raises(SystemExit):
            _parse_args(["--bam", "x.bam", "--clip", "#table"])


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


class TestHashStaleness:
    def _static(self, tmp_path: Path) -> Path:
        static = tmp_path / "static"
        (static / "dist").mkdir(parents=True)
        (static / "renderer.ts").write_text("// v1")
        (static / "viewer.html").write_text("<html>src</html>")
        (static / "dist" / "viewer.html").write_text("<html>built</html>")
        return static

    @pytest.mark.unit
    def test_fresh_stamp_ignores_newer_mtimes(self, monkeypatch, tmp_path: Path):
        """A checkout/branch switch bumps mtimes without changing content: no warning."""
        import os

        from bamcp.eval import render_cli

        static = self._static(tmp_path)
        render_cli.record_build_hash(static)
        future = time.time() + 3600
        os.utime(static / "renderer.ts", (future, future))
        monkeypatch.setattr(render_cli, "_STATIC_DIR", static)
        assert render_cli._check_rebuild_needed() is None

    @pytest.mark.unit
    def test_content_change_warns_even_with_old_mtime(self, monkeypatch, tmp_path: Path):
        import os

        from bamcp.eval import render_cli

        static = self._static(tmp_path)
        render_cli.record_build_hash(static)
        (static / "renderer.ts").write_text("// v2")
        os.utime(static / "renderer.ts", (0, 0))
        monkeypatch.setattr(render_cli, "_STATIC_DIR", static)
        warning = render_cli._check_rebuild_needed()
        assert warning is not None and "content hash" in warning

    @pytest.mark.unit
    def test_test_files_do_not_affect_hash(self, tmp_path: Path):
        from bamcp.eval import render_cli

        static = self._static(tmp_path)
        before = render_cli.viewer_source_hash(static)
        (static / "constants.test.ts").write_text("// test only")
        assert render_cli.viewer_source_hash(static) == before
        (static / "state.ts").write_text("// new source")
        assert render_cli.viewer_source_hash(static) != before

    @pytest.mark.unit
    def test_stamp_for_older_dist_is_ignored(self, monkeypatch, tmp_path: Path):
        """Rebuilt without re-stamping: the stamp no longer describes dist → mtime fallback."""
        import os

        from bamcp.eval import render_cli

        static = self._static(tmp_path)
        render_cli.record_build_hash(static)
        (static / "renderer.ts").write_text("// v2")
        (static / "dist" / "viewer.html").write_text("<html>rebuilt</html>")
        os.utime(static / "renderer.ts", (0, 0))
        monkeypatch.setattr(render_cli, "_STATIC_DIR", static)
        assert render_cli._check_rebuild_needed() is None

    @pytest.mark.unit
    def test_corrupt_stamp_falls_back(self, monkeypatch, tmp_path: Path):
        from bamcp.eval import render_cli

        static = self._static(tmp_path)
        (static / "dist" / render_cli.BUILD_STAMP_NAME).write_text("not json")
        monkeypatch.setattr(render_cli, "_STATIC_DIR", static)
        assert render_cli._read_stamp(static, static / "dist" / "viewer.html") is None

    @pytest.mark.unit
    def test_record_requires_dist(self, tmp_path: Path):
        from bamcp.eval import render_cli

        with pytest.raises(SystemExit, match="does not exist"):
            render_cli.record_build_hash(tmp_path)

    @pytest.mark.unit
    def test_main_record_build_hash(self, monkeypatch, capsys, tmp_path: Path):
        from bamcp.eval import render_cli

        static = self._static(tmp_path)
        monkeypatch.setattr(render_cli, "_STATIC_DIR", static)
        assert render_cli.main(["--record-build-hash"]) == 0
        assert render_cli.BUILD_STAMP_NAME in capsys.readouterr().out
        assert render_cli._check_rebuild_needed() is None


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


class TestBlindCli:
    def _stub(self, monkeypatch, calls: list):
        from bamcp.eval import render_cli

        async def fake_capture(self, payload, mode, **kw):
            calls.append({"mode": mode, **kw})
            return b"\x89PNG"

        async def fake_aenter(self):
            return self

        async def fake_aexit(self, *exc):
            return None

        async def fake_build(*a, **kw):
            return {"contig": "chr1", "start": 1, "end": 10}

        monkeypatch.setattr(render_cli.ViewerRenderer, "capture", fake_capture)
        monkeypatch.setattr(render_cli.ViewerRenderer, "__aenter__", fake_aenter)
        monkeypatch.setattr(render_cli.ViewerRenderer, "__aexit__", fake_aexit)
        monkeypatch.setattr(render_cli, "_build_payload_from_bam", fake_build)
        monkeypatch.setattr(render_cli, "_resolve_region", lambda bam, region: "chr1:1-10")
        monkeypatch.setattr(render_cli, "_check_rebuild_needed", lambda: None)
        return render_cli

    @pytest.mark.unit
    def test_blind_flag_reaches_capture(self, monkeypatch, tmp_path: Path):
        calls: list = []
        render_cli = self._stub(monkeypatch, calls)
        out = tmp_path / "b.png"
        assert render_cli.main(["--bam", "f.bam", "--blind", "--out", str(out)]) == 0
        assert calls[0]["blind"] is True
        assert calls[0]["viewport"] == (1024, None)
        assert calls[0]["clip"] is None

    @pytest.mark.unit
    def test_default_is_unblinded_fixed_viewport(self, monkeypatch, tmp_path: Path):
        calls: list = []
        render_cli = self._stub(monkeypatch, calls)
        render_cli.main(["--bam", "f.bam", "--out", str(tmp_path / "x.png")])
        assert calls[0]["blind"] is False
        assert calls[0]["viewport"] == (1024, 800)

    @pytest.mark.unit
    def test_blind_rejects_active_variant(self, monkeypatch, tmp_path: Path):
        render_cli = self._stub(monkeypatch, [])
        with pytest.raises(SystemExit, match="blind"):
            render_cli.main(
                ["--bam", "f.bam", "--blind", "--active-variant", "chr1:5:A>T", "--out", "x.png"]
            )
