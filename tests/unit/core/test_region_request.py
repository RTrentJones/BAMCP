"""The shared region-tool preamble (tools._region_request)."""

from __future__ import annotations

import pytest

from bamcp.config import BAMCPConfig
from bamcp.core import tools

pytestmark = pytest.mark.unit


async def test_defaults_come_from_config(small_bam_path, ref_fasta_path):
    cfg = BAMCPConfig(reference=ref_fasta_path, min_vaf=0.2, min_depth=4, min_alt_count=2)
    req = await tools._region_request({"file_path": small_bam_path, "region": "chr1:1-50"}, cfg)
    assert (req.min_vaf, req.min_depth, req.min_alt_count) == (0.2, 4, 2)
    assert req.reference == ref_fasta_path
    assert req.resolved.source == tools.REFERENCE_SOURCE_CONFIG
    assert (req.variant_source, req.vcf_path, req.vcf_primary) == ("auto", None, False)


async def test_overrides_and_variant_source(small_bam_path):
    args = {
        "file_path": small_bam_path,
        "region": "chr1:1-50",
        "min_vaf": 0.5,
        "min_depth": 9,
        "min_alt_count": 7,
        "variant_source": "bamcp",
        "vcf_path": "/ignored.vcf.gz",
    }
    req = await tools._region_request(args, BAMCPConfig())
    assert (req.min_vaf, req.min_depth, req.min_alt_count) == (0.5, 9, 7)
    assert (req.vcf_path, req.vcf_primary) == (None, False)
    with pytest.raises(ValueError, match="requires a vcf_path"):
        await tools._region_request(
            {**args, "variant_source": "vcf", "vcf_path": None}, BAMCPConfig()
        )


async def test_region_override_is_validated(small_bam_path):
    args = {"file_path": small_bam_path}
    req = await tools._region_request(args, BAMCPConfig(), region="chr1:10-20")
    assert req.region == "chr1:10-20"
    with pytest.raises(ValueError, match="Invalid region format"):
        await tools._region_request(args, BAMCPConfig(), region="chr1;bad")
