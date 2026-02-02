"""MCP tool handlers for BAMCP."""

import json
from typing import Any

import pysam

from .config import BAMCPConfig
from .parsers import RegionData, fetch_region


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
