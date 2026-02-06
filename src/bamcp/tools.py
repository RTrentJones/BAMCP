"""MCP tool handlers for BAMCP."""

from __future__ import annotations

import json
import logging
from dataclasses import asdict
from typing import Any

import pysam

from .clinvar import ClinVarClient
from .config import BAMCPConfig
from .gnomad import GnomadClient
from .parsers import RegionData, fetch_region

logger = logging.getLogger(__name__)


async def handle_browse_region(args: dict[str, Any], config: BAMCPConfig) -> dict:
    """
    Handle browse_region tool call.

    Returns serialized region data with UI metadata.
    """
    file_path = args["file_path"]
    region = args["region"]
    reference = args.get("reference", config.reference)

    data = fetch_region(
        file_path,
        region,
        reference,
        max_reads=config.max_reads,
        min_mapq=config.min_mapq,
    )

    payload = _serialize_region_data(data)

    return {
        "content": [{"type": "text", "text": json.dumps(payload)}],
        "_meta": {
            "ui/resourceUri": "ui://bamcp/viewer",
            "ui/init": payload,
        },
    }


async def handle_get_variants(args: dict[str, Any], config: BAMCPConfig) -> dict:
    """Return variants without UI."""
    file_path = args["file_path"]
    region = args["region"]
    reference = args.get("reference", config.reference)
    min_vaf = args.get("min_vaf", config.min_vaf)
    min_depth = args.get("min_depth", config.min_depth)

    data = fetch_region(
        file_path,
        region,
        reference,
        max_reads=config.max_reads,
        min_mapq=config.min_mapq,
    )

    variants = [v for v in data.variants if v["vaf"] >= min_vaf and v["depth"] >= min_depth]

    return {
        "content": [
            {"type": "text", "text": json.dumps({"variants": variants, "count": len(variants)})}
        ]
    }


async def handle_get_coverage(args: dict[str, Any], config: BAMCPConfig) -> dict:
    """Return coverage statistics."""
    file_path = args["file_path"]
    region = args["region"]
    reference = args.get("reference", config.reference)

    data = fetch_region(
        file_path,
        region,
        reference,
        max_reads=config.max_reads,
        min_mapq=config.min_mapq,
    )

    coverage = data.coverage
    stats = {
        "region": f"{data.contig}:{data.start}-{data.end}",
        "mean": round(sum(coverage) / len(coverage), 2) if coverage else 0,
        "min": min(coverage) if coverage else 0,
        "max": max(coverage) if coverage else 0,
        "median": sorted(coverage)[len(coverage) // 2] if coverage else 0,
        "bases_covered": sum(1 for c in coverage if c > 0),
        "total_bases": len(coverage),
    }

    return {"content": [{"type": "text", "text": json.dumps(stats)}]}


async def handle_list_contigs(args: dict[str, Any], config: BAMCPConfig) -> dict:
    """List contigs in a BAM/CRAM file."""
    file_path = args["file_path"]
    reference = args.get("reference", config.reference)

    mode: str = "rc" if file_path.endswith(".cram") else "rb"
    samfile = pysam.AlignmentFile(file_path, mode, reference_filename=reference)  # type: ignore[arg-type]

    contigs = [
        {"name": name, "length": length}
        for name, length in zip(samfile.references, samfile.lengths, strict=True)
    ]

    samfile.close()

    return {"content": [{"type": "text", "text": json.dumps({"contigs": contigs})}]}


async def handle_jump_to(args: dict[str, Any], config: BAMCPConfig) -> dict:
    """
    Handle jump_to tool call.

    Centers the viewer on a specific genomic position with a configurable window.
    """
    file_path = args["file_path"]
    position = args["position"]
    contig = args.get("contig", "chr1")
    window = args.get("window", config.default_window)
    reference = args.get("reference", config.reference)

    start = max(0, position - window // 2)
    end = position + window // 2
    region = f"{contig}:{start}-{end}"

    data = fetch_region(
        file_path,
        region,
        reference,
        max_reads=config.max_reads,
        min_mapq=config.min_mapq,
    )

    payload = _serialize_region_data(data)

    return {
        "content": [{"type": "text", "text": json.dumps(payload)}],
        "_meta": {
            "ui/resourceUri": "ui://bamcp/viewer",
            "ui/init": payload,
        },
    }


async def handle_visualize_region(args: dict[str, Any], config: BAMCPConfig) -> dict:
    """
    Handle visualize_region tool call.

    Primary MCP Apps tool — returns serialized region data with UI metadata.
    Equivalent to browse_region but named to clarify its App-centric purpose.
    """
    file_path = args["file_path"]
    region = args["region"]
    reference = args.get("reference", config.reference)

    data = fetch_region(
        file_path,
        region,
        reference,
        max_reads=config.max_reads,
        min_mapq=config.min_mapq,
    )

    payload = _serialize_region_data(data)

    return {
        "content": [{"type": "text", "text": json.dumps(payload)}],
        "_meta": {
            "ui/resourceUri": "ui://bamcp/viewer",
            "ui/init": payload,
        },
    }


async def handle_get_region_summary(args: dict[str, Any], config: BAMCPConfig) -> dict:
    """
    Handle get_region_summary tool call.

    Text-only summary for LLM reasoning — no UI metadata.
    """
    file_path = args["file_path"]
    region = args["region"]
    reference = args.get("reference", config.reference)

    data = fetch_region(
        file_path,
        region,
        reference,
        max_reads=config.max_reads,
        min_mapq=config.min_mapq,
    )

    coverage = data.coverage
    mean_cov = round(sum(coverage) / len(coverage), 2) if coverage else 0
    max_cov = max(coverage) if coverage else 0

    summary_lines = [
        f"Region: {data.contig}:{data.start}-{data.end}",
        f"Reads: {len(data.reads)}",
        f"Coverage: mean={mean_cov}x, max={max_cov}x",
        f"Variants detected: {len(data.variants)}",
    ]

    for v in data.variants:
        summary_lines.append(
            f"  {v['contig']}:{v['position']} {v['ref']}>{v['alt']} "
            f"VAF={v['vaf']:.1%} depth={v['depth']}"
        )

    return {"content": [{"type": "text", "text": "\n".join(summary_lines)}]}


_CLINVAR_DISCLAIMER = (
    "Note: This is research-grade information from ClinVar and is not intended "
    "for clinical diagnostic use."
)

_GNOMAD_DISCLAIMER = (
    "Note: This is research-grade population frequency data from gnomAD and is "
    "not intended for clinical diagnostic use."
)


async def handle_lookup_clinvar(args: dict[str, Any], config: BAMCPConfig) -> dict:
    """
    Look up a variant in ClinVar via NCBI E-utilities.

    Returns clinical significance, review status, and associated conditions.
    """
    chrom = args["chrom"]
    pos = args["pos"]
    ref = args["ref"]
    alt = args["alt"]

    client = ClinVarClient(api_key=config.ncbi_api_key)

    try:
        result = await client.lookup(chrom, pos, ref, alt)
    except Exception:
        logger.exception("ClinVar lookup failed for %s:%d %s>%s", chrom, pos, ref, alt)
        return {
            "content": [
                {
                    "type": "text",
                    "text": json.dumps(
                        {"error": "ClinVar lookup failed", "disclaimer": _CLINVAR_DISCLAIMER}
                    ),
                }
            ]
        }

    if result is None:
        return {
            "content": [
                {
                    "type": "text",
                    "text": json.dumps(
                        {
                            "found": False,
                            "message": f"No ClinVar entry found for {chrom}:{pos} {ref}>{alt}",
                            "disclaimer": _CLINVAR_DISCLAIMER,
                        }
                    ),
                }
            ]
        }

    payload = asdict(result)
    payload["disclaimer"] = _CLINVAR_DISCLAIMER

    return {"content": [{"type": "text", "text": json.dumps(payload)}]}


async def handle_lookup_gnomad(args: dict[str, Any], config: BAMCPConfig) -> dict:
    """
    Look up a variant in gnomAD for population allele frequency data.

    Returns global and per-population allele frequencies.
    """
    chrom = args["chrom"]
    pos = args["pos"]
    ref = args["ref"]
    alt = args["alt"]

    client = GnomadClient(dataset=config.gnomad_dataset)

    try:
        result = await client.lookup(chrom, pos, ref, alt)
    except Exception:
        logger.exception("gnomAD lookup failed for %s:%d %s>%s", chrom, pos, ref, alt)
        return {
            "content": [
                {
                    "type": "text",
                    "text": json.dumps(
                        {"error": "gnomAD lookup failed", "disclaimer": _GNOMAD_DISCLAIMER}
                    ),
                }
            ]
        }

    if result is None:
        return {
            "content": [
                {
                    "type": "text",
                    "text": json.dumps(
                        {
                            "found": False,
                            "message": f"No gnomAD entry found for {chrom}:{pos} {ref}>{alt}",
                            "disclaimer": _GNOMAD_DISCLAIMER,
                        }
                    ),
                }
            ]
        }

    payload = asdict(result)
    payload["disclaimer"] = _GNOMAD_DISCLAIMER

    return {"content": [{"type": "text", "text": json.dumps(payload)}]}


def _serialize_region_data(data: RegionData) -> dict:
    """Serialize RegionData to a JSON-compatible dict."""
    return {
        "contig": data.contig,
        "start": data.start,
        "end": data.end,
        "reads": [
            {
                "name": r.name,
                "sequence": r.sequence,
                "qualities": r.qualities,
                "cigar": r.cigar,
                "position": r.position,
                "end_position": r.end_position,
                "mapping_quality": r.mapping_quality,
                "is_reverse": r.is_reverse,
                "mismatches": r.mismatches,
            }
            for r in data.reads
        ],
        "coverage": data.coverage,
        "variants": data.variants,
        "reference_sequence": data.reference_sequence,
    }
