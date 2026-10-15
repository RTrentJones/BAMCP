"""Curated public demo datasets and the ``demo:<NAME>`` file_path alias.

A newcomer should be able to try BAMCP without knowing a BAM URL: every tool that takes
``file_path`` also accepts ``demo:HG002`` / ``demo:HG00096``. :func:`resolve_demo_alias` swaps
the alias for the real public URL *before* the handler runs, so the normal path/URL validation
(SSRF checks, the host allow-list) still applies to the resolved URL — the alias is a
convenience, never a bypass.

Every locus below was confirmed against the live server (get_variants with the public
reference, lookup_clinvar / lookup_gnomad in the BAM's build). Regions are 1-based inclusive,
samtools-style, like every region at the tool boundary.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any

from .core.errors import ToolInputError
from .core.reference import CONTIG_STYLE_CHR, CONTIG_STYLE_NOCHR, get_public_reference_url

DEMO_PREFIX = "demo:"


@dataclass(frozen=True)
class DemoLocus:
    """A curated region inside a demo dataset (1-based inclusive ``region``)."""

    region: str
    label: str
    note: str


@dataclass(frozen=True)
class DemoDataset:
    """A public, index-backed alignment file reachable through a ``demo:<NAME>`` alias."""

    alias: str
    file_path: str
    build: str
    contig_style: str
    reference_url: str | None
    description: str = ""
    loci: list[DemoLocus] = field(default_factory=list)


_HG002_URL = (
    "https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/data/AshkenazimTrio/"
    "HG002_NA24385_son/NIST_HiSeq_HG002_Homogeneity-10953946/"
    "NHGRI_Illumina300X_AJtrio_novoalign_bams/HG002.GRCh38.60x.1.bam"
)
_HG00096_URL = (
    "https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/phase3/data/HG00096/alignment/"
    "HG00096.mapped.ILLUMINA.bwa.GBR.low_coverage.20120522.bam"
)

DEMO_DATASETS: dict[str, DemoDataset] = {
    "HG002": DemoDataset(
        alias="demo:HG002",
        file_path=_HG002_URL,
        build="GRCh38",
        contig_style=CONTIG_STYLE_CHR,
        reference_url=get_public_reference_url("GRCh38", CONTIG_STYLE_CHR),
        description=(
            "Genome in a Bottle HG002 (NA24385, Ashkenazi son), Illumina ~60x, GRCh38, "
            "chr-style contigs. A well-characterised public reference sample."
        ),
        loci=[
            DemoLocus(
                region="chr17:43099864-43099964",
                label="BRCA1 real heterozygous SNV chr17:43099914 G>A",
                note=(
                    "Clean het SNV: ~49% of ~57 reads carry A. ClinVar Benign (3 stars, "
                    "expert panel; BRCA1 c.442-34C>T); gnomAD v4 AF ~0.21 — a common "
                    "polymorphism, not a disease finding."
                ),
            ),
            DemoLocus(
                region="chr17:43096101-43096200",
                label="BRCA1 artifact cluster around chr17:43096151 T>G",
                note=(
                    "A cluster of low-VAF mismatches (43096151 T>G ~22%, plus 43096165, "
                    "43096197 ...) carried by low-base-quality bases on one strand: the "
                    "classic look of a sequencing artifact. Compare it with the het SNV."
                ),
            ),
            DemoLocus(
                region="chr17:43109038-43109138",
                label="BRCA1 common polymorphism chr17:43109088 G>C",
                note=(
                    "Clean het SNV (~50% of ~68 reads). ClinVar Benign (3 stars; BRCA1 "
                    "c.135-2555C>G); gnomAD v4 AF ~0.023 overall, ~0.07 in Ashkenazi "
                    "Jewish — note gnomAD flags this site (frequency_reliable=false)."
                ),
            ),
            DemoLocus(
                region="chr17:43044295-43125483",
                label="BRCA1 gene (GRCh38)",
                note=(
                    "Gene-scale region: use get_coverage / get_region_summary here, and "
                    "zoom to a ≤200 bp window before visualize_region."
                ),
            ),
        ],
    ),
    "HG00096": DemoDataset(
        alias="demo:HG00096",
        file_path=_HG00096_URL,
        build="GRCh37",
        contig_style=CONTIG_STYLE_NOCHR,
        reference_url=get_public_reference_url("GRCh37", CONTIG_STYLE_NOCHR),
        description=(
            "1000 Genomes phase 3 HG00096 (GBR), low coverage (~4-8x), GRCh37, contigs "
            "named '1', '13', ... (no 'chr'). Shows build and naming differences."
        ),
        loci=[
            DemoLocus(
                region="13:32911838-32911938",
                label="BRCA2 common homozygous SNV 13:32911888 A>G",
                note=(
                    "All ~7 reads carry G (homozygous alt, but at low depth). ClinVar "
                    "Benign (3 stars; BRCA2 c.3396A>G p.Lys1132=); gnomAD v2.1 AF ~0.29."
                ),
            ),
            DemoLocus(
                region="13:32889611-32973805",
                label="BRCA2 gene (GRCh37)",
                note=(
                    "Gene-scale region on a nochr GRCh37 BAM: contig is '13', not 'chr13'. "
                    "Low coverage means many single-read candidates — mostly noise."
                ),
            ),
        ],
    ),
}


def available_aliases() -> list[str]:
    """Every accepted alias, e.g. ``['demo:HG002', 'demo:HG00096']``."""
    return [d.alias for d in DEMO_DATASETS.values()]


def get_demo(alias: str) -> DemoDataset | None:
    """Look up a dataset by ``demo:NAME`` (or bare ``NAME``), case-insensitively."""
    name = alias[len(DEMO_PREFIX) :] if alias.lower().startswith(DEMO_PREFIX) else alias
    for key, dataset in DEMO_DATASETS.items():
        if key.lower() == name.strip().lower():
            return dataset
    return None


def resolve_demo_alias(args: dict[str, Any]) -> dict[str, Any]:
    """Return ``args`` with a ``demo:<NAME>`` ``file_path`` replaced by the dataset's URL.

    Anything else passes through untouched (same dict object), so this is safe to apply to
    every tool call. An unknown ``demo:`` alias raises :class:`ToolInputError` naming the
    available demos.
    """
    file_path = args.get("file_path")
    if not isinstance(file_path, str) or not file_path.lower().startswith(DEMO_PREFIX):
        return args
    dataset = get_demo(file_path)
    if dataset is None:
        raise ToolInputError(
            f"Unknown demo dataset '{file_path[:64]}'. Available demos: "
            + ", ".join(available_aliases())
            + " (see the bamcp://demos resource)."
        )
    return {**args, "file_path": dataset.file_path}


def demos_as_dict() -> dict[str, Any]:
    """JSON-ready listing for the ``bamcp://demos`` resource."""
    return {
        "coordinate_system": "1-based-inclusive",
        "usage": (
            "Pass the alias as file_path to any tool (e.g. file_path='demo:HG002'). Call "
            "list_contigs first; pass reference_url as reference to enable SNV detection "
            "if the server does not pick one automatically."
        ),
        "datasets": [
            {
                "alias": d.alias,
                "file_path": d.file_path,
                "build": d.build,
                "contig_style": d.contig_style,
                "reference_url": d.reference_url,
                "description": d.description,
                "loci": [
                    {"region": loc.region, "label": loc.label, "note": loc.note} for loc in d.loci
                ],
            }
            for d in DEMO_DATASETS.values()
        ],
    }


__all__ = [
    "CONTIG_STYLE_CHR",
    "CONTIG_STYLE_NOCHR",
    "DEMO_DATASETS",
    "DEMO_PREFIX",
    "DemoDataset",
    "DemoLocus",
    "available_aliases",
    "demos_as_dict",
    "get_demo",
    "resolve_demo_alias",
]
