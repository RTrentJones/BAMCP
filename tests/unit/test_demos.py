"""Demo datasets and the ``demo:<NAME>`` file_path alias."""

from __future__ import annotations

import re

import pytest

from bamcp.core.errors import ToolInputError
from bamcp.core.reference import get_public_reference_url
from bamcp.core.validation import REGION_PATTERN
from bamcp.demos import (
    DEMO_DATASETS,
    available_aliases,
    demos_as_dict,
    get_demo,
    resolve_demo_alias,
)
from bamcp.glossary import GLOSSARY, glossary_markdown, lookup_terms

_RANGE = re.compile(r"^(?P<contig>[^:]+):(?P<start>\d+)-(?P<end>\d+)$")


@pytest.mark.unit
class TestDemoDatasets:
    def test_expected_demos(self):
        assert set(DEMO_DATASETS) == {"HG002", "HG00096"}
        assert available_aliases() == ["demo:HG002", "demo:HG00096"]

    @pytest.mark.parametrize("name", list(DEMO_DATASETS))
    def test_reference_url_matches_registry(self, name):
        d = DEMO_DATASETS[name]
        assert d.reference_url is not None
        assert d.reference_url == get_public_reference_url(d.build, d.contig_style)

    @pytest.mark.parametrize("name", list(DEMO_DATASETS))
    def test_loci_are_valid_1_based_regions_in_bam_naming(self, name):
        d = DEMO_DATASETS[name]
        assert d.alias == f"demo:{name}"
        assert d.file_path.startswith("https://") and d.file_path.endswith(".bam")
        assert d.loci
        small = 0
        for locus in d.loci:
            assert REGION_PATTERN.match(locus.region), locus.region
            m = _RANGE.match(locus.region)
            assert m, locus.region
            start, end = int(m["start"]), int(m["end"])
            assert 1 <= start <= end
            assert m["contig"].startswith("chr") == (d.contig_style == "chr")
            small += end - start + 1 <= 200
            assert locus.label and locus.note
        assert small >= 1  # at least one viewer-sized locus

    def test_hg002_curated_sites(self):
        regions = " ".join(loc.label for loc in DEMO_DATASETS["HG002"].loci)
        assert "43099914 G>A" in regions
        assert "43096151 T>G" in regions

    def test_demos_as_dict_shape(self):
        data = demos_as_dict()
        assert data["coordinate_system"] == "1-based-inclusive"
        aliases = [d["alias"] for d in data["datasets"]]
        assert aliases == available_aliases()
        assert all(d["loci"] for d in data["datasets"])


@pytest.mark.unit
class TestResolveDemoAlias:
    def test_resolves_alias(self):
        out = resolve_demo_alias({"file_path": "demo:HG002", "region": "chr17:1-10"})
        assert out == {"file_path": DEMO_DATASETS["HG002"].file_path, "region": "chr17:1-10"}

    def test_case_insensitive(self):
        out = resolve_demo_alias({"file_path": "DEMO:hg00096"})
        assert out["file_path"] == DEMO_DATASETS["HG00096"].file_path

    def test_does_not_mutate_input(self):
        args = {"file_path": "demo:HG002"}
        resolve_demo_alias(args)
        assert args == {"file_path": "demo:HG002"}

    @pytest.mark.parametrize(
        "args",
        [{}, {"file_path": "/data/x.bam"}, {"file_path": "https://h/x.bam"}, {"chrom": "chr1"}],
    )
    def test_passthrough(self, args):
        assert resolve_demo_alias(args) is args

    def test_unknown_demo_is_clear_error(self):
        with pytest.raises(ToolInputError) as exc:
            resolve_demo_alias({"file_path": "demo:NOPE"})
        msg = str(exc.value)
        assert "demo:NOPE" in msg
        assert "demo:HG002" in msg and "demo:HG00096" in msg

    def test_get_demo_bare_name(self):
        assert get_demo("HG002") is DEMO_DATASETS["HG002"]
        assert get_demo("demo:missing") is None


@pytest.mark.unit
class TestGlossary:
    REQUIRED = [
        "VAF",
        "Depth",
        "MAPQ",
        "Base quality",
        "Strand bias",
        "Read-position bias",
        "Homopolymer",
        "Artifact risk",
        "Confidence",
        "gnomAD allele frequency",
        "frequency_reliable",
        "ClinVar stars",
        "applied_filters",
        "1-based",
        "reference_source",
    ]

    @pytest.mark.parametrize("term", REQUIRED)
    def test_term_present(self, term):
        assert lookup_terms([term]), term

    def test_every_term_has_conclusions(self):
        for t in GLOSSARY:
            assert t.meaning and t.conclude and t.do_not_conclude

    def test_markdown_carries_safety_framing(self):
        md = glossary_markdown()
        assert "not a diagnostic" in md
        assert md.count("**Do not conclude:**") == len(GLOSSARY)

    def test_lookup_ignores_blank(self):
        assert lookup_terms(["", "  "]) == []


@pytest.mark.unit
def test_gene_scale_loci_are_flagged_by_span():
    from bamcp.demos import DemoLocus

    loci = [loc for d in DEMO_DATASETS.values() for loc in d.loci]
    genes = {loc.region for loc in loci if loc.gene_scale}
    assert genes == {"chr17:43044295-43125483", "13:32889214-32974405"}
    assert DemoLocus("chr1:1-200", "x", "y").span == 200
