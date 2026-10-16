"""Plain-language glossary of the terms BAMCP output uses.

Served as the ``bamcp://glossary`` MCP resource and used by the ``explain_terms`` prompt. Each
entry says what a term means, and — following SAFETY.md — what a curious person *should* and
*should not* conclude from it. BAMCP is a research/visualization tool, not a diagnostic device.
"""

from __future__ import annotations

from dataclasses import dataclass


@dataclass(frozen=True)
class GlossaryTerm:
    term: str
    meaning: str
    conclude: str
    do_not_conclude: str


GLOSSARY: list[GlossaryTerm] = [
    GlossaryTerm(
        "VAF (variant allele fraction)",
        "The share of reads at a position that show the alternative base: alt_count / depth. "
        "0.5 is what one copy out of two (heterozygous) usually looks like; ~1.0 means both "
        "copies (homozygous).",
        "Around 0.4-0.6 with good depth fits a heterozygous variant; near 1.0 fits homozygous.",
        "A VAF alone is not a genotype. At low depth it is very noisy (2 of 4 reads = 0.5), "
        "and low VAFs (<0.1-0.2) are often sequencing errors, not mosaicism.",
    ),
    GlossaryTerm(
        "Depth (coverage)",
        "How many reads overlap a position. More reads means more evidence.",
        "~30x or more gives a reasonably stable picture of a germline variant.",
        "A variant seen in 1-3 reads is a hint, not a finding. Low depth cannot rule a variant "
        "in or out.",
    ),
    GlossaryTerm(
        "MAPQ (mapping quality)",
        "The aligner's confidence that a read belongs where it was placed (Phred scale: 20 "
        "≈ 1% chance it is misplaced, 60 = very confident). 0 means it fits equally well "
        "elsewhere.",
        "Evidence from high-MAPQ reads is about this locus.",
        "Mismatches carried mostly by low-MAPQ reads often come from a similar sequence "
        "elsewhere in the genome (paralog/pseudogene), not from this gene.",
    ),
    GlossaryTerm(
        "Base quality",
        "The sequencer's confidence in each individual base (Phred: 20 ≈ 1% error, 30 ≈ "
        "0.1%). Grey/faded bases in the viewer are low quality.",
        "An alternative base backed by high-quality bases is more believable.",
        "A cluster of mismatches with low base quality (mean below ~13-20) is the typical "
        "look of a sequencing error, not a variant.",
    ),
    GlossaryTerm(
        "Strand bias",
        "Reads come from both DNA strands (forward and reverse). A real variant usually shows "
        "up on both; an artifact often appears on only one.",
        "Alt reads split across both strands support a real variant.",
        "Alt reads on one strand only (e.g. 11 forward / 0 reverse) are a red flag, though "
        "not proof of an artifact — check depth and base quality too.",
    ),
    GlossaryTerm(
        "Read-position bias",
        "Whether the alternative base sits mostly at the ends of reads, where sequencing and "
        "alignment errors concentrate.",
        "Alt bases spread across read positions support a real variant.",
        "Alt bases piled at read ends suggest an alignment or chemistry artifact.",
    ),
    GlossaryTerm(
        "Homopolymer",
        "A run of the same base (e.g. AAAAAAAA). Sequencers often miscount run length.",
        "Indels outside repetitive context are easier to trust.",
        "A one-base insertion/deletion inside a long homopolymer is frequently an error, even "
        "at moderate VAF.",
    ),
    GlossaryTerm(
        "Artifact risk",
        "BAMCP's combined score (low/medium/high) of the warning signs above: strand bias, "
        "read-position bias, low MAPQ, low base quality, homopolymer context, low depth.",
        "High artifact risk means: treat the call as probably not real until proven otherwise.",
        "Low artifact risk does not mean the variant is real or important — it means none of "
        "the checked warning signs fired.",
    ),
    GlossaryTerm(
        "Confidence",
        "BAMCP's overall confidence that the read evidence supports a real variant, "
        "downgraded by artifact risk and low depth. By design it is never 'high' when artifact "
        "risk is high.",
        "High confidence means the reads look like a clean variant.",
        "Confidence is about the read evidence, not about whether the variant causes disease. "
        "It is not a clinical call.",
    ),
    GlossaryTerm(
        "gnomAD allele frequency (AF)",
        "How common the variant is in the gnomAD population database (0.2 = 20% of "
        "chromosomes carry it). popmax / faf95 give the highest frequency across populations.",
        "A variant carried by more than a few percent of people is very unlikely to cause a "
        "rare, severe disease (that is what ACMG BA1/BS1 encode).",
        "Absence from gnomAD does not make a variant harmful; many harmless variants are rare.",
    ),
    GlossaryTerm(
        "gnomAD reliability (frequency_reliable)",
        "Whether gnomAD's own QC passed at that site. false means gnomAD flagged it (e.g. "
        "discrepant exome/genome frequencies) or too few people were sequenced.",
        "With frequency_reliable=true the AF is a usable population estimate.",
        "With frequency_reliable=false, do not lean on the exact AF number.",
    ),
    GlossaryTerm(
        "ClinVar stars (review status)",
        "ClinVar collects labs' interpretations (Benign ... Pathogenic). Stars (0-4) say how "
        "well reviewed: 0 = no assertion criteria, 1 = single submitter, 2 = multiple "
        "agreeing, 3 = expert panel, 4 = practice guideline.",
        "A 3-4 star classification is a strong, reviewed opinion.",
        "A 0-1 star entry may be one lab's old opinion; 'not found' only means nobody "
        "submitted it, not that it is benign.",
    ),
    GlossaryTerm(
        "applied_filters",
        "The thresholds a variant response used (min_vaf, min_depth, min_mapq, min_baseq, "
        "min_alt_count, excluded read types) and anything that could not be checked "
        "(snv_detection disabled without a reference, truncated reads).",
        "A variant listed passed these filters.",
        "A missing variant may simply have been filtered out or undetectable — absence is not "
        "proof the site is reference. Loosen filters to look again, or ask "
        "get_variant_curation_summary about the exact site: it evaluates it anyway and lists "
        "the filters it fails (filter_check).",
    ),
    GlossaryTerm(
        "1-based coordinates",
        "Positions count from 1, like VCF, samtools and genome browsers. 'chr17:100-200' "
        "includes both base 100 and base 200; 'chr17:150-150' (or just 'chr17:150') is one base.",
        "Positions from BAMCP can be pasted into ClinVar, gnomAD or a genome browser as-is "
        "(same genome build).",
        "Coordinates from a different genome build (GRCh37 vs GRCh38) do not line up — "
        "they are ~0.5-2 Mb off. Call list_contigs first to pin the build.",
    ),
    GlossaryTerm(
        "reference_source",
        "Where the reference genome used for SNV detection came from: supplied by you "
        "('explicit'), configured by the server operator ('config'), or chosen automatically "
        "to match the BAM's build and contig naming ('auto'). 'none' means there was no "
        "reference, so SNV/indel detection is disabled.",
        "With a matching reference, mismatches against it are real differences in the reads.",
        "Without one, a lack of SNVs means nothing; with the wrong build, everything is off.",
    ),
]

_SAFETY = (
    "BAMCP is a research and visualization tool, not a diagnostic device. Variant calls are "
    "exploratory, and every finding needs confirmation by a qualified person using a "
    "validated pipeline. Common benign polymorphisms (like those in the demos) are normal — "
    "everyone carries thousands."
)


def glossary_markdown() -> str:
    """Render the glossary as markdown (the ``bamcp://glossary`` resource body)."""
    lines = [
        "# BAMCP glossary",
        "",
        f"> {_SAFETY}",
        "",
    ]
    for t in GLOSSARY:
        lines += [
            f"## {t.term}",
            "",
            t.meaning,
            "",
            f"- **You can conclude:** {t.conclude}",
            f"- **Do not conclude:** {t.do_not_conclude}",
            "",
        ]
    return "\n".join(lines).rstrip() + "\n"


def lookup_terms(names: list[str]) -> list[GlossaryTerm]:
    """Glossary entries whose term contains any of ``names`` (case-insensitive)."""
    wanted = [n.strip().lower() for n in names if n.strip()]
    return [t for t in GLOSSARY if any(w in t.term.lower() for w in wanted)]
