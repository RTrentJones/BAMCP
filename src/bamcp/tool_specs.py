"""Canonical tool descriptions and JSON schemas.

Single source of truth for the LLM-facing description and parameter schema of
each BAMCP tool. The FastMCP server (``server.py``) derives its schemas from
function signatures; this registry mirrors those so the eval harness exposes the
*same* descriptions and parameter schemas a real MCP client would see — instead
of the placeholder ``{additionalProperties: true}`` it used before. Because each
schema sets ``additionalProperties: False``, a parameter missing here is one the
eval model can never pass — so ``tests/unit/test_tool_specs.py`` asserts, for every
tool, that the property and required sets match the schema FastMCP derives from
``server.py``. Add or rename a parameter in both places.

Each entry is ``{"description": str, "input_schema": <JSON Schema object>}``.
"""

from __future__ import annotations

from typing import Any

_STR = {"type": "string"}
_INT = {"type": "integer"}
_NUM = {"type": "number"}
_VCF_PATH = {**_STR, "description": "Optional VCF/BCF of the caller's variant calls"}
_VARIANT_SOURCE = {
    **_STR,
    "enum": ["auto", "vcf", "bamcp"],
    "description": "auto: BAMCP candidates with the VCF overlaid; vcf: the VCF is "
    "authoritative (requires vcf_path); bamcp: BAMCP candidates only",
}


# The region-string contract every region-taking tool shares (parse_region converts it).
REGION_FORMAT_NOTE = (
    "region is 1-based and inclusive, samtools-style ('chr1:1000-2000' covers bases "
    "1000..2000; 'chr1:1051-1051' or just 'chr1:1051' is the single base 1051)."
)
_REGION = {**_STR, "description": "1-based inclusive, samtools-style, e.g. chr1:1000-2000"}

# Reference resolution shared by the variant tools (see core/tools._resolve_reference_info).
AUTO_REFERENCE_NOTE = (
    "If reference is omitted, a public FASTA matching the BAM's detected build and contig "
    "naming is used automatically when available (applied_filters.reference_source says "
    "explicit/config/auto/none)."
)
# The detection operating point (constants.DEFAULT_MIN_*), loosenable per call.
DETECTION_DEFAULTS_NOTE = (
    "Defaults: min_vaf 0.10, min_alt_count 3 alt reads, MAPQ>=20, base quality>=13; all "
    "echoed in applied_filters. Lower min_vaf/min_alt_count to see weaker candidates."
)
_MIN_ALT_COUNT = {
    **_INT,
    "description": "Minimum alt-supporting reads for a BAMCP candidate (default 3)",
}


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
        "1-based (VCF/dbSNP convention); "
        + REGION_FORMAT_NOTE
        + " Each variant carries strand counts, mean base "
        "quality, artifact-risk assessment, a confidence level, concerns, and a plain-language "
        "interpretation with a zygosity hint. " + DETECTION_DEFAULTS_NOTE,
        "input_schema": _obj(
            {
                "file_path": _STR,
                "region": _REGION,
                "reference": _STR,
                "vcf_path": _VCF_PATH,
                "variant_source": _VARIANT_SOURCE,
                "min_vaf": _NUM,
                "min_depth": _INT,
                "min_alt_count": _MIN_ALT_COUNT,
            },
            ["file_path", "region"],
        ),
    },
    "get_coverage": {
        "description": "Calculate depth-of-coverage statistics (mean, min, max, median, "
        "bases covered) for a genomic region. " + REGION_FORMAT_NOTE,
        "input_schema": _obj(
            {"file_path": _STR, "region": _REGION, "reference": _STR},
            ["file_path", "region"],
        ),
    },
    "list_contigs": {
        "description": "Detect a BAM/CRAM's genome build (GRCh37 vs GRCh38) and contig naming "
        "style, and list its primary chromosomes. Use this first on new BAM files. The summary "
        "(build, style, suggested public reference URL, contig counts) comes first; alt/decoy/"
        "HLA/unplaced contigs are counted, not listed, unless include_all=true.",
        "input_schema": _obj(
            {
                "file_path": _STR,
                "reference": _STR,
                "include_all": {
                    "type": "boolean",
                    "description": "List every contig, not just the primary chromosomes",
                },
            },
            ["file_path"],
        ),
    },
    "jump_to": {
        "description": "Jump to a specific genomic position and view the surrounding reads in "
        "the interactive viewer. contig and position are both required; position is 1-based "
        "(VCF/samtools convention). Contig names follow the BAM (list_contigs: 'chr17' vs '17').",
        "input_schema": _obj(
            {
                "file_path": _STR,
                "position": _INT,
                "contig": _STR,
                "window": _INT,
                "reference": _STR,
                "vcf_path": _VCF_PATH,
                "variant_source": _VARIANT_SOURCE,
            },
            ["file_path", "contig", "position"],
        ),
    },
    "visualize_region": {
        "description": "Visualize aligned reads in a genomic region with the interactive MCP "
        "Apps viewer. Auto-detects compact mode for large regions. " + REGION_FORMAT_NOTE,
        "input_schema": _obj(
            {
                "file_path": _STR,
                "region": _REGION,
                "reference": _STR,
                "vcf_path": _VCF_PATH,
                "variant_source": _VARIANT_SOURCE,
            },
            ["file_path", "region"],
        ),
    },
    "get_region_summary": {
        "description": "Get a text summary of a genomic region (read count, coverage, and each "
        "variant with a plain-language interpretation: zygosity hint, artifact risk, confidence) "
        "for LLM reasoning. No UI. " + REGION_FORMAT_NOTE + " " + AUTO_REFERENCE_NOTE,
        "input_schema": _obj(
            {
                "file_path": _STR,
                "region": _REGION,
                "reference": _STR,
                "vcf_path": _VCF_PATH,
                "variant_source": _VARIANT_SOURCE,
                "min_alt_count": _MIN_ALT_COUNT,
            },
            ["file_path", "region"],
        ),
    },
    "lookup_clinvar": {
        "description": "Look up a variant in ClinVar for clinical significance, review status "
        "(star rating), and associated conditions. Research-grade, not for clinical use. "
        "Positions are 1-based, in 'build' (GRCh37/hg19 or GRCh38/hg38; default: server "
        "build) — pass the BAM's build from list_contigs, or a GRCh37 variant is searched "
        "at the wrong locus.",
        "input_schema": _obj(
            {"chrom": _STR, "pos": _INT, "ref": _STR, "alt": _STR, "build": _STR},
            ["chrom", "pos", "ref", "alt"],
        ),
    },
    "lookup_gnomad": {
        "description": "Look up a variant in gnomAD for population allele-frequency data "
        "(global and per-population AF). Research-grade, not for clinical use. Positions are "
        "1-based, in 'build' (GRCh37/hg19 → gnomAD v2.1, GRCh38/hg38 → gnomAD v4; default: "
        "server build) — pass the BAM's build from list_contigs.",
        "input_schema": _obj(
            {"chrom": _STR, "pos": _INT, "ref": _STR, "alt": _STR, "build": _STR},
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
        "Requires a reference genome (" + AUTO_REFERENCE_NOTE + ") Returns up to 500 variants "
        "ranked by VAF. If contig is omitted, the BAM's first contig is scanned. "
        + DETECTION_DEFAULTS_NOTE,
        "input_schema": _obj(
            {
                "file_path": _STR,
                "contig": _STR,
                "reference": _STR,
                "min_vaf": _NUM,
                "min_depth": _INT,
                "min_alt_count": _MIN_ALT_COUNT,
            },
            ["file_path"],
        ),
    },
    "search_gene": {
        "description": "Search for a gene by symbol (e.g. BRCA1, TP53) and return its genomic "
        "coordinates (1-based). Use this to navigate to a gene's location by name. Pass "
        "file_path to return coordinates in that BAM's build (auto-detected), or build "
        "(GRCh37/hg19 or GRCh38/hg38) to choose it explicitly. The returned region is "
        "1-based inclusive, ready to pass to visualize_region.",
        "input_schema": _obj(
            {
                "symbol": _STR,
                "file_path": _STR,
                "build": {**_STR, "description": "GRCh37/hg19 or GRCh38/hg38"},
            },
            ["symbol"],
        ),
    },
    "cleanup_cache": {
        "description": "Clear the server's downloaded BAM index files and region cache "
        "(single-client stdio servers only; shared HTTP servers expire their cache by TTL).",
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
