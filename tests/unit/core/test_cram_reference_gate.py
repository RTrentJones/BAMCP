"""CRAM decode on a remote-only server must not fall back to the header's ``@SQ UR:`` tag."""

from __future__ import annotations

import hashlib
import re
from pathlib import Path

import pysam
import pytest

from bamcp.config import BAMCPConfig
from bamcp.core.errors import ToolInputError
from bamcp.core.tools import (
    CRAM_REFERENCE_REQUIRED_MESSAGE,
    _check_cram_header_references,
    _fetch_region_with_timeout,
    _require_cram_reference,
)

pytestmark = pytest.mark.unit

FIXTURES = Path(__file__).resolve().parents[2] / "fixtures"


def _cram_with_header_ur(tmp_path: Path, ur: str | None = None) -> tuple[Path, Path]:
    """Write a CRAM whose only @SQ carries ``UR:`` (default ``file://<local fasta>``).

    Returns ``(cram, fasta)``; the FASTA is the one the reads were encoded against.
    """
    ref = tmp_path / "header_ref.fa"
    ref.write_text((FIXTURES / "ref.fa").read_text())
    pysam.faidx(str(ref))
    with pysam.FastaFile(str(ref)) as fa:
        seq = fa.fetch("chr1").upper()
        length = fa.get_reference_length("chr1")
    header = {
        "HD": {"VN": "1.6", "SO": "coordinate"},
        "SQ": [
            {
                "SN": "chr1",
                "LN": length,
                "M5": hashlib.md5(seq.encode(), usedforsecurity=False).hexdigest(),
                "UR": f"file://{ref}",
            }
        ],
    }
    cram = tmp_path / "header_ur.cram"
    with (
        pysam.AlignmentFile(str(FIXTURES / "small.bam")) as src,
        pysam.AlignmentFile(str(cram), "wc", header=header, reference_filename=str(ref)) as out,
    ):
        for read in src.fetch("chr1"):
            d = read.to_dict()
            if d["next_ref_name"] not in ("=", "chr1"):
                d["next_ref_name"] = "*"
            out.write(pysam.AlignedSegment.from_dict(d, out.header))
    if ur is not None:
        # htslib rewrites UR: to the encoding reference on write, so set the attacker's
        # value afterwards, as `samtools reheader` would.
        with pysam.AlignmentFile(str(cram)) as f:
            text = str(f.header)
        header_sam = tmp_path / "header.sam"
        header_sam.write_text(re.sub(r"UR:[^\t\n]+", f"UR:{ur}", text))
        pysam.samtools.reheader("-i", str(header_sam), str(cram))
    pysam.index(str(cram))
    return cram, ref


class TestRequireCramReference:
    def test_remote_only_cram_without_reference_refused(self):
        config = BAMCPConfig(allow_local_files=False)
        with pytest.raises(ToolInputError, match="CRAM decoding needs a reference"):
            _require_cram_reference("https://example.org/x.cram", None, config)

    def test_resolved_reference_passes(self):
        config = BAMCPConfig(allow_local_files=False)
        _require_cram_reference("https://example.org/x.cram", "https://example.org/r.fa", config)

    def test_bam_without_reference_passes(self):
        config = BAMCPConfig(allow_local_files=False)
        _require_cram_reference("https://example.org/x.bam", None, config)

    def test_local_files_allowed_is_unchanged(self):
        _require_cram_reference("https://example.org/x.cram", None, BAMCPConfig())


class TestFetchGate:
    async def test_header_ur_reference_is_used_when_local_files_allowed(self, tmp_path):
        # Documents the htslib behavior the gate exists for: with no reference, the CRAM
        # decodes against the file its header's UR: tag names.
        cram, _ = _cram_with_header_ur(tmp_path)
        data = await _fetch_region_with_timeout(
            str(cram), "chr1:1-200", None, BAMCPConfig(cache_dir=str(tmp_path / "c"))
        )
        assert data.reads

    async def test_remote_only_refuses_before_opening(self, tmp_path):
        cram, _ = _cram_with_header_ur(tmp_path)
        config = BAMCPConfig(allow_local_files=False, cache_dir=str(tmp_path / "c"))
        with pytest.raises(ToolInputError) as exc:
            await _fetch_region_with_timeout(str(cram), "chr1:1-200", None, config)
        assert str(exc.value) == CRAM_REFERENCE_REQUIRED_MESSAGE


class TestHeaderReferenceFallback:
    """A resolved reference lacking a header contig still lets htslib follow that @SQ UR:."""

    @pytest.fixture(autouse=True)
    def _clear_checks(self):
        from bamcp.core import tools

        tools._CRAM_HEADER_REF_CHECKS.clear()
        yield
        tools._CRAM_HEADER_REF_CHECKS.clear()

    def _other_ref(self, tmp_path: Path) -> Path:
        other = tmp_path / "other.fa"
        other.write_text(">other\nACGTACGTACGT\n")
        pysam.faidx(str(other))
        return other

    async def test_htslib_falls_back_to_header_ur(self, tmp_path):
        cram, _ = _cram_with_header_ur(tmp_path)
        other = self._other_ref(tmp_path)
        with pysam.AlignmentFile(str(cram), "rc", reference_filename=str(other)) as f:
            seqs = [r.query_sequence for r in f.fetch(until_eof=True)]
        assert seqs and all(seqs)

    async def test_remote_only_refuses_uncovered_ur_contig(self, tmp_path):
        cram, _ = _cram_with_header_ur(tmp_path)
        other = self._other_ref(tmp_path)
        config = BAMCPConfig(allow_local_files=False, cache_dir=str(tmp_path / "c"))
        with pytest.raises(ToolInputError, match=r"\(chr1\) at an external reference \(@SQ UR:\)"):
            await _fetch_region_with_timeout(str(cram), "chr1:1-200", str(other), config)

    async def test_remote_only_covering_reference_decodes(self, tmp_path):
        cram, ref = _cram_with_header_ur(tmp_path)
        config = BAMCPConfig(allow_local_files=False, cache_dir=str(tmp_path / "c"))
        data = await _fetch_region_with_timeout(str(cram), "chr1:1-200", str(ref), config)
        assert data.reads

    async def test_local_files_allowed_permits_local_ur(self, tmp_path):
        # With local files allowed the caller could name that FASTA directly, so a local
        # UR: is no escalation and decoding proceeds.
        cram, _ = _cram_with_header_ur(tmp_path)
        other = self._other_ref(tmp_path)
        await _check_cram_header_references(str(cram), str(other), BAMCPConfig())


class TestNetworkHeaderReference:
    """htslib itself refuses a URL in @SQ UR: ("UR tags pointing to remote files are not
    supported"), so BAMCP's header gate only needs to handle local paths. Pin that: if a
    pysam/htslib upgrade starts following remote UR: tags, this fails before any release."""

    @pytest.fixture
    def listener(self):
        import http.server
        import threading

        hits: list[str] = []

        class Handler(http.server.BaseHTTPRequestHandler):
            def do_GET(self):  # noqa: N802 - stdlib hook name
                hits.append(self.path)
                self.send_response(404)
                self.end_headers()

            def log_message(self, *args):
                pass

        server = http.server.HTTPServer(("127.0.0.1", 0), Handler)
        threading.Thread(target=server.serve_forever, daemon=True).start()
        yield server.server_port, hits
        server.shutdown()

    def test_htslib_never_fetches_a_remote_ur(self, tmp_path, listener, monkeypatch):
        port, hits = listener
        for var in ("HTTP_PROXY", "HTTPS_PROXY", "http_proxy", "https_proxy"):
            monkeypatch.delenv(var, raising=False)  # reach the loopback listener directly
        monkeypatch.setenv("REF_PATH", str(tmp_path / "refs" / "%s"))  # no M5 hit either
        cram, _ = _cram_with_header_ur(tmp_path, ur=f"http://127.0.0.1:{port}/internal.fa")
        with pytest.raises(OSError), pysam.AlignmentFile(str(cram), "rc") as f:
            for _read in f.fetch(until_eof=True):
                pass
        assert hits == []
