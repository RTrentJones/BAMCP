"""Remote index download hardening: per-hop redirect validation and temp-file cleanup."""

from __future__ import annotations

from pathlib import Path
from unittest.mock import patch

import httpx
import pytest

from bamcp.config import BAMCPConfig
from bamcp.core import tools as tools_module

_HOST_IPS = {
    "data.example.com": "93.184.216.34",
    "cdn.example.com": "93.184.216.35",
    "internal.example.com": "10.0.0.5",
}


def _fake_getaddrinfo(host, port, proto=0):
    return [(2, 1, 6, "", (_HOST_IPS.get(host, "93.184.216.36"), port))]


@pytest.fixture(autouse=True)
def _hermetic_dns():
    with patch("bamcp.core.validation.socket.getaddrinfo", side_effect=_fake_getaddrinfo):
        yield


def _client(handler) -> httpx.AsyncClient:
    return httpx.AsyncClient(transport=httpx.MockTransport(handler), follow_redirects=False)


def _tmp_files(directory: Path) -> list[Path]:
    return list(directory.glob(".index-*.tmp"))


@pytest.mark.unit
class TestIndexRedirects:
    async def test_follows_validated_relative_redirect(self, tmp_path):
        seen: list[str] = []

        def handler(request: httpx.Request) -> httpx.Response:
            seen.append(str(request.url))
            if request.url.path == "/x.bam.bai":
                return httpx.Response(302, headers={"location": "/moved/x.bam.bai"})
            return httpx.Response(200, content=b"BAI\x01")

        out = tmp_path / "x.bai"
        async with _client(handler) as client:
            result = await tools_module._download_index_streaming(
                client, "https://data.example.com/x.bam.bai", str(out), BAMCPConfig()
            )

        assert result == str(out)
        assert out.read_bytes() == b"BAI\x01"
        assert seen == [
            "https://data.example.com/x.bam.bai",
            "https://data.example.com/moved/x.bam.bai",
        ]

    async def test_redirect_to_internal_host_blocked(self, tmp_path):
        fetched: list[str] = []

        def handler(request: httpx.Request) -> httpx.Response:
            fetched.append(request.url.host)
            return httpx.Response(
                301, headers={"location": "http://internal.example.com/latest/meta-data"}
            )

        async with _client(handler) as client:
            with pytest.raises(ValueError, match="private/internal"):
                await tools_module._download_index_streaming(
                    client,
                    "https://data.example.com/x.bam.bai",
                    str(tmp_path / "x.bai"),
                    BAMCPConfig(),
                )
        # The internal host was never contacted.
        assert fetched == ["data.example.com"]

    async def test_redirect_off_allow_list_blocked(self, tmp_path):
        def handler(request: httpx.Request) -> httpx.Response:
            return httpx.Response(307, headers={"location": "https://cdn.example.com/x.bai"})

        config = BAMCPConfig(allowed_remote_hosts=["data.example.com"])
        async with _client(handler) as client:
            with pytest.raises(ValueError, match="not in the allowed remote hosts"):
                await tools_module._download_index_streaming(
                    client, "https://data.example.com/x.bam.bai", str(tmp_path / "x.bai"), config
                )

    async def test_redirect_to_non_http_scheme_blocked(self, tmp_path):
        def handler(request: httpx.Request) -> httpx.Response:
            return httpx.Response(302, headers={"location": "file:///etc/passwd"})

        async with _client(handler) as client:
            with pytest.raises(ValueError, match="Scheme not supported"):
                await tools_module._download_index_streaming(
                    client,
                    "https://data.example.com/x.bam.bai",
                    str(tmp_path / "x.bai"),
                    BAMCPConfig(),
                )

    async def test_redirect_loop_bounded(self, tmp_path):
        calls = 0

        def handler(request: httpx.Request) -> httpx.Response:
            nonlocal calls
            calls += 1
            return httpx.Response(302, headers={"location": f"/hop{calls}.bai"})

        async with _client(handler) as client:
            with pytest.raises(ValueError, match="redirects"):
                await tools_module._download_index_streaming(
                    client,
                    "https://data.example.com/x.bam.bai",
                    str(tmp_path / "x.bai"),
                    BAMCPConfig(),
                )
        assert calls == tools_module._INDEX_MAX_REDIRECTS + 1

    async def test_redirect_without_location_is_a_miss(self, tmp_path):
        async with _client(lambda request: httpx.Response(302)) as client:
            result = await tools_module._download_index_streaming(
                client,
                "https://data.example.com/x.bam.bai",
                str(tmp_path / "x.bai"),
                BAMCPConfig(),
            )
        assert result is None

    async def test_non_200_is_a_miss(self, tmp_path):
        async with _client(lambda request: httpx.Response(500)) as client:
            result = await tools_module._download_index_streaming(
                client,
                "https://data.example.com/x.bam.bai",
                str(tmp_path / "x.bai"),
                BAMCPConfig(),
            )
        assert result is None


@pytest.mark.unit
class TestTempFileCleanup:
    async def test_oversized_stream_leaves_no_tmp(self, tmp_path):
        async def body():
            yield b"abc"
            yield b"def"

        def handler(request: httpx.Request) -> httpx.Response:
            return httpx.Response(200, content=body())  # no Content-Length

        async with _client(handler) as client:
            with pytest.raises(ValueError, match="exceeds configured maximum"):
                await tools_module._download_index_streaming(
                    client,
                    "https://data.example.com/x.bam.bai",
                    str(tmp_path / "x.bai"),
                    BAMCPConfig(max_remote_index_bytes=4),
                )
        assert _tmp_files(tmp_path) == []
        assert not (tmp_path / "x.bai").exists()

    async def test_network_error_mid_stream_leaves_no_tmp(self, tmp_path):
        async def body():
            yield b"abc"
            raise httpx.ReadError("connection reset")

        def handler(request: httpx.Request) -> httpx.Response:
            return httpx.Response(200, content=body())

        async with _client(handler) as client:
            with pytest.raises(httpx.ReadError):
                await tools_module._download_index_streaming(
                    client,
                    "https://data.example.com/x.bam.bai",
                    str(tmp_path / "x.bai"),
                    BAMCPConfig(),
                )
        assert _tmp_files(tmp_path) == []
        assert not (tmp_path / "x.bai").exists()


@pytest.mark.unit
class TestEnsureCachedIndex:
    async def test_download_enforces_budget_and_prunes_lock(self, tmp_path, monkeypatch):
        from bamcp.core.tools import _ensure_cached_index, get_services

        def handler(request: httpx.Request) -> httpx.Response:
            return httpx.Response(200, content=b"BAI\x01")

        real_client = httpx.AsyncClient

        def client_factory(**kwargs):
            assert kwargs.get("follow_redirects") is False
            return real_client(transport=httpx.MockTransport(handler), **kwargs)

        monkeypatch.setattr(tools_module.httpx, "AsyncClient", client_factory)
        config = BAMCPConfig(cache_dir=str(tmp_path), allow_remote_files=True)
        services = get_services(config)
        budget_calls: list[str | None] = []
        monkeypatch.setattr(
            services.cache, "enforce_budget", lambda keep=None: budget_calls.append(keep) or 0
        )

        result = await _ensure_cached_index("https://data.example.com/S1.CRAM?sig=1", config)

        assert result is not None and result.endswith(".crai")
        assert Path(result).read_bytes() == b"BAI\x01"
        assert budget_calls == [result]
        assert services.index_download_locks == {}

    async def test_failed_download_prunes_lock(self, tmp_path, monkeypatch):
        from bamcp.core.tools import _ensure_cached_index, get_services

        real_client = httpx.AsyncClient

        def client_factory(**kwargs):
            return real_client(
                transport=httpx.MockTransport(lambda r: httpx.Response(404)), **kwargs
            )

        monkeypatch.setattr(tools_module.httpx, "AsyncClient", client_factory)
        config = BAMCPConfig(cache_dir=str(tmp_path), allow_remote_files=True)

        assert await _ensure_cached_index("https://data.example.com/x.bam", config) is None
        assert get_services(config).index_download_locks == {}
        assert _tmp_files(get_services(config).cache.cache_dir) == []
