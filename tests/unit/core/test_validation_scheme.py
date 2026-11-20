"""Scheme routing and local-file policy for the data-source validators.

Regression for an SSRF: ``http:/127.0.0.1:40409/b.vcf.gz`` (single slash) has no ``://`` so it
used to be treated as a *local* path, but htslib parses the scheme and fetches it over HTTP.
"""

import socket
import threading
from unittest.mock import patch

import pytest

from bamcp.config import BAMCPConfig
from bamcp.core.validation import (
    LOCAL_FILES_DISABLED_MESSAGE,
    validate_path,
    validate_reference_path,
    validate_variant_file_path,
)

_PUBLIC = [(2, 1, 6, "", ("93.184.216.34", 443))]

_VALIDATORS = [
    (validate_path, "x.bam"),
    (validate_variant_file_path, "x.vcf.gz"),
    (validate_reference_path, "x.fa"),
]

_SCHEME_PAYLOADS = [
    "http:/127.0.0.1:40409/{name}",
    "https:/127.0.0.1/{name}",
    "http:127.0.0.1/{name}",
    "HTTP://example.com/{name}",
    "data:,{name}",
    "ftp://example.com/{name}",
    "ftp:/example.com/{name}",
    "file:///etc/{name}",
    "file:/etc/{name}",
    "s3://bucket/{name}",
    "gs://bucket/{name}",
    "s3+https://bucket/{name}",
    "1x:/host/{name}",
]


@pytest.mark.unit
@pytest.mark.parametrize(("validator", "name"), _VALIDATORS)
@pytest.mark.parametrize("template", _SCHEME_PAYLOADS)
@pytest.mark.parametrize("remote", [False, True])
def test_non_http_schemes_rejected(validator, name, template, remote):
    config = BAMCPConfig(allow_remote_files=remote)
    with (
        patch("bamcp.core.validation.socket.getaddrinfo", return_value=_PUBLIC) as dns,
        pytest.raises(ValueError, match="Scheme not supported"),
    ):
        validator(template.format(name=name), config)
    dns.assert_not_called()


@pytest.mark.unit
@pytest.mark.parametrize(("validator", "name"), _VALIDATORS)
@pytest.mark.parametrize(
    "template", ["/data/http://evil/{name}", "local.bam##idx##/tmp/{name}", "a/b://c/{name}"]
)
def test_embedded_urls_in_local_paths_rejected(validator, name, template):
    with pytest.raises(ValueError, match="embedded URLs"):
        validator(template.format(name=name), BAMCPConfig(allow_remote_files=True))


@pytest.mark.unit
@pytest.mark.parametrize(("validator", "name"), _VALIDATORS)
@pytest.mark.parametrize("scheme", ["http", "https"])
def test_http_urls_go_through_validate_remote_url(validator, name, scheme):
    url = f"{scheme}://example.com/{name}"
    config = BAMCPConfig(allow_remote_files=True)
    with patch("bamcp.core.validation.validate_remote_url") as remote:
        validator(url, config)
    remote.assert_called_once_with(url, config)


@pytest.mark.unit
@pytest.mark.parametrize(("validator", "name"), _VALIDATORS)
def test_http_url_still_rejected_when_remote_disabled(validator, name):
    with pytest.raises(ValueError, match="Remote files are disabled"):
        validator(f"https://example.com/{name}", BAMCPConfig())


@pytest.mark.unit
def test_real_local_files_still_accepted(tmp_path):
    vcf = tmp_path / "calls.vcf.gz"
    bam = tmp_path / "reads.bam"
    ref = tmp_path / "ref.fa"
    for f in (vcf, bam, ref):
        f.write_bytes(b"x")
    config = BAMCPConfig(allowed_directories=[str(tmp_path)])
    validate_variant_file_path(str(vcf), config)
    validate_path(str(bam), config)
    validate_reference_path(str(ref), config)
    # Relative paths with ':' later in the string are not schemes.
    validate_variant_file_path(str(tmp_path / "." / "calls.vcf.gz"), config)


@pytest.mark.unit
def test_missing_local_vcf_is_not_found(tmp_path):
    with pytest.raises(FileNotFoundError):
        validate_variant_file_path(str(tmp_path / "nope.vcf.gz"), BAMCPConfig())


@pytest.mark.unit
def test_vcf_allow_list_checked_before_existence(tmp_path):
    allowed = tmp_path / "allowed"
    allowed.mkdir()
    with pytest.raises(ValueError, match="not in allowed directories"):
        validate_variant_file_path(
            str(tmp_path / "nope.vcf.gz"), BAMCPConfig(allowed_directories=[str(allowed)])
        )


class TestLocalFilesDisabled:
    @pytest.mark.unit
    @pytest.mark.parametrize(("validator", "name"), _VALIDATORS)
    def test_local_paths_rejected_before_filesystem_access(self, validator, name, tmp_path):
        existing = tmp_path / name
        existing.write_bytes(b"x")
        config = BAMCPConfig(allow_local_files=False)
        with patch("bamcp.core.validation.Path") as path_cls:
            for candidate in (str(existing), str(tmp_path / ("missing-" + name))):
                with pytest.raises(ValueError, match=LOCAL_FILES_DISABLED_MESSAGE):
                    validator(candidate, config)
        path_cls.assert_not_called()  # no existence oracle

    @pytest.mark.unit
    @pytest.mark.parametrize(("validator", "name"), _VALIDATORS)
    def test_remote_urls_still_work(self, validator, name):
        config = BAMCPConfig(allow_local_files=False, allow_remote_files=True)
        with patch("bamcp.core.validation.socket.getaddrinfo", return_value=_PUBLIC):
            validator(f"https://example.com/{name}", config)

    @pytest.mark.unit
    def test_operator_reference_still_accepted(self, tmp_path):
        ref = tmp_path / "reference.fa"
        ref.write_bytes(b">chr1\nACGT\n")
        config = BAMCPConfig(allow_local_files=False, reference=str(ref))
        validate_reference_path(str(ref), config)
        # ...but a different caller-chosen local FASTA is not.
        other = tmp_path / "other.fa"
        other.write_bytes(b">chr1\nACGT\n")
        with pytest.raises(ValueError, match=LOCAL_FILES_DISABLED_MESSAGE):
            validate_reference_path(str(other), config)

    @pytest.mark.unit
    @pytest.mark.parametrize("allow_local", [True, False])
    def test_configured_reference_url_passed_explicitly(self, allow_local):
        """A configured *URL* reference named by the caller goes through remote validation."""
        url = "https://example.com/hg38.fa.gz"
        config = BAMCPConfig(allow_local_files=allow_local, allow_remote_files=True, reference=url)
        with patch("bamcp.core.validation.socket.getaddrinfo", return_value=_PUBLIC):
            validate_reference_path(url, config)  # not "configured reference file not found"
        # ...and is still SSRF-checked like any other remote reference.
        private = [(2, 1, 6, "", ("10.0.0.5", 443))]
        with (
            patch("bamcp.core.validation.socket.getaddrinfo", return_value=private),
            pytest.raises(ValueError),
        ):
            validate_reference_path(url, config)

    @pytest.mark.unit
    async def test_config_reference_resolution_not_rejected(self, tmp_path):
        from bamcp.core.tools import _resolve_reference_info

        ref = tmp_path / "reference.fa"
        ref.write_bytes(b">chr1\nACGT\n")
        config = BAMCPConfig(allow_local_files=False, reference=str(ref))
        resolved = await _resolve_reference_info({"file_path": "https://h/x.bam"}, config)
        assert resolved.effective == str(ref)
        resolved = await _resolve_reference_info({"reference": str(ref)}, config)
        assert resolved.effective == str(ref)

    @pytest.mark.unit
    def test_demo_aliases_resolve_to_https(self):
        from bamcp.demos import DEMO_DATASETS, resolve_demo_alias

        for dataset in DEMO_DATASETS.values():
            resolved = resolve_demo_alias({"file_path": dataset.alias})["file_path"]
            assert resolved.startswith("https://")
            if dataset.reference_url:
                assert dataset.reference_url.startswith("https://")

    @pytest.mark.unit
    def test_env_flag_parsed_strictly(self, monkeypatch):
        monkeypatch.setenv("BAMCP_ALLOW_LOCAL_FILES", "false")
        assert BAMCPConfig.from_env().allow_local_files is False
        monkeypatch.setenv("BAMCP_ALLOW_LOCAL_FILES", "nope")
        with pytest.raises(ValueError, match="BAMCP_ALLOW_LOCAL_FILES"):
            BAMCPConfig.from_env()
        monkeypatch.delenv("BAMCP_ALLOW_LOCAL_FILES")
        assert BAMCPConfig.from_env().allow_local_files is True


@pytest.mark.unit
async def test_get_variants_single_slash_vcf_never_fetched():
    """The reproduced exploit: no request may reach a loopback listener."""
    from bamcp.core.tools import handle_get_variants

    hits: list[bytes] = []
    srv = socket.socket()
    srv.bind(("127.0.0.1", 0))
    srv.listen(4)
    srv.settimeout(0.5)
    port = srv.getsockname()[1]
    stop = threading.Event()

    def _accept() -> None:
        while not stop.is_set():
            try:
                conn, _ = srv.accept()
            except OSError:
                continue
            hits.append(conn.recv(200))
            conn.close()

    t = threading.Thread(target=_accept, daemon=True)
    t.start()
    try:
        with pytest.raises(ValueError, match="Scheme not supported"):
            await handle_get_variants(
                {
                    "file_path": "tests/fixtures/small.bam",
                    "region": "chr1:100-200",
                    "vcf_path": f"http:/127.0.0.1:{port}/b.vcf.gz",
                },
                BAMCPConfig(allow_remote_files=True),
            )
    finally:
        stop.set()
        t.join(2)
        srv.close()
    assert hits == []
