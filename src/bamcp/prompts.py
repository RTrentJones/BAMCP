"""Newcomer surface: server instructions, demo/glossary/guide resources, and MCP prompts.

Registered by :func:`register_demo_surface` (called once from ``create_server``). Nothing here
adds a *tool* — the tool list is pinned (the deploy verify has an exact-match drift guard) — so
discovery goes through MCP resources (``bamcp://demos``, ``bamcp://glossary``,
``bamcp://guide``) and prompts (``explore_demo``, ``check_variant``, ``explain_terms``).
"""

from __future__ import annotations

import json

from mcp.server.fastmcp import FastMCP

from .demos import DEMO_DATASETS, available_aliases, demos_as_dict, get_demo
from .glossary import GLOSSARY, glossary_markdown, lookup_terms

SERVER_INSTRUCTIONS = (
    "BAMCP inspects aligned sequencing reads (BAM/CRAM) and candidate variants, with an "
    "interactive pileup viewer. Research and education only — not for clinical use.\n"
    "- No file of your own? Use a demo alias as file_path: "
    + ", ".join(available_aliases())
    + " (curated loci in the bamcp://demos resource; walkthrough prompt: explore_demo).\n"
    "- Call list_contigs first on any new file to pin the genome build (GRCh37 vs GRCh38) and "
    "contig naming ('chr17' vs '17'); pass that build to lookup_clinvar / lookup_gnomad / "
    "search_gene.\n"
    "- Regions and positions are 1-based inclusive (samtools-style: 'chr17:43099914-43099914' "
    "is one base). Keep visualize_region/get_variants windows small (≤200 bp); for a whole "
    "gene use get_coverage, then zoom in.\n"
    "- Every variant response echoes applied_filters: a missing variant may have been filtered "
    "out, not absent. Plain-language term definitions: bamcp://glossary."
)

GUIDE_MARKDOWN = (
    """# BAMCP quick guide

BAMCP lets an assistant look at the actual sequencing reads behind a variant, the way a
genome analyst would in IGV. It is a **research and education tool, not a diagnostic test**.

## 1. Pick a file
No BAM of your own? Use a demo alias as `file_path`:

"""
    + "\n".join(
        f"- `{d.alias}` — {d.description} Build {d.build}, contigs `{d.contig_style}`."
        for d in DEMO_DATASETS.values()
    )
    + """

Curated regions are in `bamcp://demos`.

## 2. Pin the build
Call `list_contigs` first. It tells you GRCh37 vs GRCh38 and whether contigs are named
`chr17` or `17`. Coordinates from the other build are silently ~0.5-2 Mb off.

## 3. Look
- `visualize_region` — interactive pileup for a small window (≤200 bp is easiest to read).
- `get_variants` — candidate variants with their `applied_filters` (also a small window).
- `get_coverage` — depth over a whole gene; use it (or `search_gene`) before zooming in.
- `get_variant_curation_summary` — artifact risk and confidence for one variant.
- `lookup_clinvar` / `lookup_gnomad` — what others have reported (pass the BAM's build).

## 4. Interpret carefully
Read `bamcp://glossary`. Rules of thumb: a real heterozygous variant is ~50% of many reads on
both strands with good quality; a lopsided, low-quality, one-strand cluster is usually an
artifact; a variant common in gnomAD is almost never the cause of a rare disease; and nothing
here replaces confirmation by a qualified professional.

Coordinates are 1-based inclusive throughout.
"""
)


def register_demo_surface(mcp: FastMCP) -> None:
    """Register the newcomer resources and prompts on ``mcp`` (no tools)."""

    # -- Resources -----------------------------------------------------------

    @mcp.resource(
        "bamcp://demos",
        name="BAMCP demo datasets",
        description="Public demo BAMs (use their demo: alias as file_path) with curated loci.",
        mime_type="application/json",
    )
    def demos_resource() -> str:
        return json.dumps(demos_as_dict(), indent=2)

    @mcp.resource(
        "bamcp://glossary",
        name="BAMCP glossary",
        description="Plain-language meaning of VAF, depth, MAPQ, artifact risk, ClinVar stars...",
        mime_type="text/markdown",
    )
    def glossary_resource() -> str:
        return glossary_markdown()

    @mcp.resource(
        "bamcp://guide",
        name="BAMCP quick guide",
        description="How to use BAMCP in four steps, for newcomers.",
        mime_type="text/markdown",
    )
    def guide_resource() -> str:
        return GUIDE_MARKDOWN

    # -- Prompts -------------------------------------------------------------

    @mcp.prompt(
        name="explore_demo",
        description="Guided tour of a public demo BAM: a real variant vs. an artifact.",
    )
    def explore_demo(dataset: str = "demo:HG002") -> str:
        demo = get_demo(dataset) or DEMO_DATASETS["HG002"]
        small = "\n".join(
            f"   - {loc.region} — {loc.label}" for loc in demo.loci if not loc.gene_scale
        )
        genes = "\n".join(f"   - {loc.region} — {loc.label}" for loc in demo.loci if loc.gene_scale)
        ref = (
            f" If SNV detection is reported disabled, pass reference='{demo.reference_url}'."
            if demo.reference_url
            else ""
        )
        steps = [
            f"Give me a beginner-friendly tour of the BAMCP demo dataset {demo.alias} "
            f"({demo.description})\n",
            f"1. Call list_contigs with file_path='{demo.alias}' and tell me the genome build "
            "and contig naming in one sentence.",
            "2. For each small curated locus below (each ≤ ~200 bp), call get_variants and "
            f"visualize_region on exactly that window with file_path='{demo.alias}'.{ref}\n"
            + small,
            "3. For each variant you find, call get_variant_curation_summary, then "
            f"lookup_clinvar and lookup_gnomad using build='{demo.build}', and explain in plain "
            "language whether the read evidence looks real or like an artifact (depth, VAF, "
            "strand balance, base quality), and what the population and ClinVar data add."
            + (
                " For the artifact locus, call get_variant_curation_summary on the named site "
                "even if get_variants did not report it — it evaluates that site and says "
                "which filters it fails."
                if any("artifact" in loc.label for loc in demo.loci)
                else ""
            ),
        ]
        if genes:
            steps.append(
                "4. Gene-scale loci below are tens of kb: do NOT call get_variants, "
                "visualize_region or get_region_summary on the whole gene (slow, and the "
                "pileup is unreadable). Call get_coverage on the gene region (fast) for depth, "
                "or search_gene with file_path to get its coordinates, then pick a ≤200 bp "
                "window inside it for any closer look.\n" + genes
            )
        steps.append(
            f"{len(steps)}. Finish with what a curious person can and cannot conclude (see "
            "bamcp://glossary). This is research/education only, not a medical result."
        )
        return "\n".join(steps)

    @mcp.prompt(
        name="check_variant",
        description="Check whether the reads support a variant at a position or region.",
    )
    def check_variant(file_path: str = "demo:HG002", region: str = "chr17:43099914") -> str:
        target = region.strip()
        # Split on the LAST colon, like parse_region: contigs such as HLA-A*01:01:01:01
        # contain colons themselves.
        if ":" in target and "-" not in target.rsplit(":", 1)[1]:
            contig, pos = target.rsplit(":", 1)
            target = f"{contig}:{pos}-{pos}"
        return (
            f"Check the read evidence for a variant at {target} in file_path='{file_path}'.\n\n"
            "1. Call list_contigs first to pin the genome build and contig naming; adjust the "
            "contig name if needed ('chr17' vs '17').\n"
            f"2. Call get_variants for a ~100 bp window around {target} (1-based inclusive) "
            "and visualize_region for the same window. Report applied_filters, and say so if "
            "SNV detection was disabled.\n"
            "3. Call get_variant_curation_summary for the variant (it evaluates the named "
            "site even when get_variants filtered it out, and says which filters it fails), "
            "then lookup_clinvar and lookup_gnomad with the BAM's build.\n"
            "4. Summarize in plain language: is it likely real or an artifact (and why), "
            "likely heterozygous or homozygous, how common it is, and what ClinVar says with "
            "its star rating. State what cannot be concluded. Research use only."
        )

    @mcp.prompt(
        name="explain_terms",
        description="Explain BAMCP output terms (VAF, MAPQ, artifact risk...) in plain language.",
    )
    def explain_terms(terms: str = "") -> str:
        names = [t for t in terms.replace(";", ",").split(",") if t.strip()]
        entries = lookup_terms(names) if names else GLOSSARY
        if not entries:
            entries = GLOSSARY
        body = "\n".join(
            f"- {t.term}: {t.meaning} Can conclude: {t.conclude} "
            f"Do not conclude: {t.do_not_conclude}"
            for t in entries
        )
        return (
            "Explain these BAMCP terms to a curious non-specialist, with a short everyday "
            "analogy for each, and keep the 'can / cannot conclude' distinction. "
            "BAMCP is research-only, not a diagnostic test.\n\n" + body
        )
