"""No NCBI API key (or any credential parameter) reaches a log line via an httpx error."""

from __future__ import annotations

import logging
import re

import httpx
import pytest

from bamcp.clients.ncbi import redact_error
from bamcp.config import BAMCPConfig
from bamcp.core.errors import sanitize_tool_exception

pytestmark = pytest.mark.unit

SECRET = "SECRETKEY0123456789"


def _status_error(status: int = 400) -> httpx.HTTPStatusError:
    request = httpx.Request(
        "GET",
        f"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=clinvar&api_key={SECRET}",
    )
    response = httpx.Response(status, request=request)
    return httpx.HTTPStatusError(
        f"Client error for url {request.url}", request=request, response=response
    )


class TestRedactError:
    def test_status_error_keeps_type_status_and_url_without_key(self):
        text = redact_error(_status_error(403))
        assert SECRET not in text
        assert "HTTPStatusError" in text and "HTTP 403" in text
        assert "esearch.fcgi" in text and "db=clinvar" in text
        assert "api_key=REDACTED" in text

    def test_request_error_url_redacted(self):
        request = httpx.Request("GET", f"https://gnomad.example/api?token={SECRET}")
        text = redact_error(httpx.ConnectError(f"failed {request.url}", request=request))
        assert SECRET not in text and "ConnectError" in text

    @pytest.mark.parametrize("param", ["api_key", "key", "token", "access_token", "API_KEY"])
    def test_other_exception_text_scrubbed(self, param):
        text = redact_error(ValueError(f"bad url https://h/x?{param}={SECRET}&db=gene"))
        assert SECRET not in text and "db=gene" in text

    def test_plain_message_unchanged(self):
        assert redact_error(ValueError("no key here")) == "ValueError: no key here"


class TestLogSites:
    async def test_clinvar_lookup_failure_log_has_no_key(self, httpx_mock, caplog):
        from bamcp.core.lookups import handle_lookup_clinvar

        httpx_mock.add_response(url=re.compile(r".*/esearch\.fcgi.*"), status_code=400)
        config = BAMCPConfig(ncbi_api_key=SECRET)
        with caplog.at_level(logging.DEBUG):
            result = await handle_lookup_clinvar(
                {"chrom": "chr17", "pos": 7674220, "ref": "G", "alt": "A", "build": "GRCh38"},
                config,
            )
        assert "unavailable" in result["content"][0]["text"]
        assert "HTTP 400" in caplog.text
        assert SECRET not in caplog.text
        # The request really carried the key: the redaction is what kept it out of the log.
        assert SECRET in str(httpx_mock.get_requests()[0].url)

    async def test_fusion_lookup_failure_log_has_no_key(self, caplog, monkeypatch):
        from bamcp.analysis import fusion
        from bamcp.core import tools

        class Boom:
            async def lookup(self, *a, **k):
                raise _status_error()

        monkeypatch.setattr(tools, "get_clinvar_client", lambda config: Boom())
        with caplog.at_level(logging.DEBUG):
            await fusion._lookup_clinvar_safe(BAMCPConfig(), "chr1", 1, "A", "G")
        assert "HTTP 400" in caplog.text and SECRET not in caplog.text

    def test_unexpected_httpx_error_in_a_tool_logs_no_key(self, caplog):
        with caplog.at_level(logging.DEBUG):
            err = sanitize_tool_exception(_status_error(), "lookup_clinvar", BAMCPConfig())
        assert SECRET not in caplog.text and SECRET not in str(err)
        assert all(r.exc_info is None for r in caplog.records)
