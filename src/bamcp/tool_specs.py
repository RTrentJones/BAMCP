"""Canonical tool descriptions and JSON schemas.

Single source of truth for the LLM-facing description and parameter schema of
each BAMCP tool. The FastMCP server (``server.py``) derives its schemas from
function signatures; this registry mirrors those so the eval harness exposes the
*same* descriptions and parameter schemas a real MCP client would see — instead
of the placeholder ``{additionalProperties: true}`` it used before. Keeping the
two in sync is a documentation/fidelity concern; a test asserts the registry
covers every routed tool.

Each entry is ``{"description": str, "input_schema": <JSON Schema object>}``.
"""

from __future__ import annotations

from typing import Any

_STR = {"type": "string"}
_INT = {"type": "integer"}
_NUM = {"type": "number"}


def _obj(props: dict[str, Any], required: list[str]) -> dict[str, Any]:
    return {
        "type": "object",
        "properties": props,
        "required": required,
        "additionalProperties": False,
    }


TOOL_SPECS: dict[str, dict[str, Any]] = {
    "get_variants": {
        "description": "Detect and return variants in a genomic region. Positions are "
        "1-based (VCF/dbSNP convention). Each variant carries strand counts, mean base "
        "quality, artifact-risk assessment, and a confidence level.",
        "input_schema": _obj(
            {
                "file_path": _STR,
                "region": {**_STR, "description": "e.g. chr1:1000-2000"},
                "reference": _STR,
                "min_vaf": _NUM,
                "min_depth": _INT,
            },
            ["file_path", "region"],
        ),
    },
    "get_coverage": {
        "description": "Calculate depth-of-coverage statistics (mean, min, max, median, "
        "bases covered) for a genomic region.",
        "input_schema": _obj(
            {"file_path": _STR, "region": _STR, "reference": _STR},
            ["file_path", "region"],
        ),
    },
    "list_contigs": {
        "description": "List chromosomes/contigs in a BAM/CRAM file and detect genome build "
        "(GRCh37 vs GRCh38). Use this first on new BAM files. Returns detected build, "
        "confidence, and a suggested public reference URL if available.",
        "input_schema": _obj({"file_path": _STR, "reference": _STR}, ["file_path"]),
    },
    "jump_to": {
        "description": "Jump to a specific genomic position and view the surrounding reads in "
        "the interactive viewer.",
        "input_schema": _obj(
            {
                "file_path": _STR,
                "position": _INT,
                "contig": _STR,
                "window": _INT,
                "reference": _STR,
            },
            ["file_path", "position"],
        ),
    },
    "visualize_region": {
        "description": "Visualize aligned reads in a genomic region with the interactive MCP "
        "Apps viewer. Auto-detects compact mode for large regions.",
        "input_schema": _obj(
            {"file_path": _STR, "region": _STR, "reference": _STR},
            ["file_path", "region"],
        ),
    },
    "get_region_summary": {
        "description": "Get a text summary of a genomic region (read count, coverage, and each "
        "variant with confidence, artifact risk, and strand balance) for LLM reasoning. No UI.",
        "input_schema": _obj(
            {"file_path": _STR, "region": _STR, "reference": _STR},
            ["file_path", "region"],
        ),
    },
    "lookup_clinvar": {
        "description": "Look up a variant in ClinVar for clinical significance, review status "
        "(star rating), and associated conditions. Research-grade, not for clinical use. "
        "Positions are 1-based.",
        "input_schema": _obj(
            {"chrom": _STR, "pos": _INT, "ref": _STR, "alt": _STR},
            ["chrom", "pos", "ref", "alt"],
        ),
    },
    "lookup_gnomad": {
        "description": "Look up a variant in gnomAD for population allele-frequency data "
        "(global and per-population AF). Research-grade, not for clinical use. Positions are "
        "1-based.",
        "input_schema": _obj(
            {"chrom": _STR, "pos": _INT, "ref": _STR, "alt": _STR},
            ["chrom", "pos", "ref", "alt"],
        ),
    },
    "get_variant_curation_summary": {
        "description": "Get a detailed curation summary for a specific variant with artifact-risk "
        "assessment, quality metrics, and curator recommendations. Set format='rubric' for a "
        "machine-scorable JSON with 0-1 quality scores. Positions are 1-based.",
        "input_schema": _obj(
            {
                "file_path": _STR,
                "chrom": _STR,
                "pos": _INT,
                "ref": _STR,
                "alt": _STR,
                "window": _INT,
                "reference": _STR,
                "format": {**_STR, "enum": ["text", "rubric"]},
                "include_clinical": {"type": "boolean"},
            },
            ["file_path", "chrom", "pos", "ref", "alt"],
        ),
    },
    "scan_variants": {
        "description": "Scan an entire contig for variants using fast coverage-based detection. "
        "Requires a reference genome. Returns up to 500 variants ranked by VAF. Use list_contigs "
        "first to detect the genome build and get a reference URL.",
        "input_schema": _obj(
            {
                "file_path": _STR,
                "contig": _STR,
                "reference": _STR,
                "min_vaf": _NUM,
                "min_depth": _INT,
            },
            ["file_path"],
        ),
    },
    "search_gene": {
        "description": "Search for a gene by symbol (e.g. BRCA1, TP53) and return its genomic "
        "coordinates. Use this to navigate to a gene's location by name.",
        "input_schema": _obj({"symbol": _STR}, ["symbol"]),
    },
    "cleanup_cache": {
        "description": "Clean up this session's cached BAM index files.",
        "input_schema": _obj({}, []),
    },
    "classify_variant": {
        "description": "Assemble a multi-source evidence package (BAM observation + ClinVar + "
        "gnomAD) for a variant and return an ACMG reasoning scaffold. Apply each applicable "
        "criterion, then give a classification (Pathogenic / Likely Pathogenic / VUS / Likely "
        "Benign / Benign) with confidence. Research-grade, not a clinical determination. "
        "Positions are 1-based.",
        "input_schema": _obj(
            {
                "file_path": _STR,
                "chrom": _STR,
                "pos": _INT,
                "ref": _STR,
                "alt": _STR,
                "gene": _STR,
                "window": _INT,
                "reference": _STR,
            },
            ["file_path", "chrom", "pos", "ref", "alt"],
        ),
    },
}


def get_tool_spec(name: str) -> dict[str, Any] | None:
    """Return the ``{description, input_schema}`` spec for a tool, or None."""
    return TOOL_SPECS.get(name)
