"""load_vcf_variants against REAL pysam-written VCFs (no VariantFile stubs).

Stubbed records hid pysam's ``ValueError: Invalid header`` on undeclared INFO keys,
which every bcftools/GATK/DeepVariant small-variant VCF triggers (none declare
SVTYPE/END/SVLEN).
"""

from __future__ import annotations

import pysam
import pytest

from bamcp.core.parsers import load_vcf_variants


def _header(
    info: tuple[str, ...] = (),
    fmt: tuple[str, ...] = (),
    samples: tuple[str, ...] = (),
    contigs: tuple[tuple[str, int], ...] = (("chr1", 1000), ("chr2", 5000), ("chr3", 5000)),
) -> pysam.VariantHeader:
    h = pysam.VariantHeader()
    for name, length in contigs:
        h.add_line(f"##contig=<ID={name},length={length}>")
    defs = {
        "DP": '##INFO=<ID=DP,Number=1,Type=Integer,Description="Depth">',
        "AF": '##INFO=<ID=AF,Number=A,Type=Float,Description="AF">',
        "SVTYPE": '##INFO=<ID=SVTYPE,Number=1,Type=String,Description="SV type">',
        "END": '##INFO=<ID=END,Number=1,Type=Integer,Description="End">',
        "SVLEN": '##INFO=<ID=SVLEN,Number=A,Type=Integer,Description="SV length">',
    }
    for key in info:
        h.add_line(defs[key])
    fdefs = {
        "GT": '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">',
        "AD": '##FORMAT=<ID=AD,Number=R,Type=Integer,Description="Allele depths">',
        "DP": '##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Sample depth">',
    }
    for key in fmt:
        h.add_line(fdefs[key])
    for s in samples:
        h.add_sample(s)
    return h


def _write(tmp_path, header, records, name="calls.vcf.gz", index=True) -> str:
    """Write records (dicts) to a VCF; bgzip + tabix it unless ``index`` is False."""
    path = str(tmp_path / name)
    mode = "wz" if path.endswith(".gz") else "w"
    with pysam.VariantFile(path, mode, header=header) as out:
        for spec in records:
            rec = out.new_record(
                contig=spec["contig"], start=spec["pos"] - 1, alleles=spec["alleles"]
            )
            info = dict(spec.get("info", {}))
            end = info.pop("END", None)
            for k, v in info.items():
                rec.info[k] = v
            if end is not None:
                # pysam reserves INFO/END (written via record.stop); set it last since
                # htslib re-derives the span from SVLEN when that is set afterwards.
                rec.stop = end
            for sample, fields in spec.get("samples", {}).items():
                for k, v in fields.items():
                    rec.samples[sample][k] = v
            out.write(rec)
    if index and path.endswith(".gz"):
        pysam.tabix_index(path, preset="vcf", force=True)
    return path


class TestUndeclaredInfo:
    @pytest.mark.unit
    def test_small_variant_vcf_without_sv_header_lines_loads(self, tmp_path):
        """Typical bcftools output: only DP/AF declared — SVTYPE/END/SVLEN lookups must
        not raise ``Invalid header``."""
        vcf = _write(
            tmp_path,
            _header(info=("DP", "AF")),
            [{"contig": "chr1", "pos": 134, "alleles": ("T", "G"), "info": {"DP": 20, "AF": 0.25}}],
        )
        [v] = load_vcf_variants(vcf, "chr1:101-200")
        assert (v["position"], v["ref"], v["alt"]) == (133, "T", "G")
        assert v["variant_kind"] == "snv"
        assert (v["vaf"], v["depth"], v["alt_count"]) == (0.25, 20, 5)
        assert (v["sv_type"], v["sv_end"], v["sv_len"]) == (None, None, None)

    @pytest.mark.unit
    def test_no_info_declared_at_all(self, tmp_path):
        vcf = _write(tmp_path, _header(), [{"contig": "chr1", "pos": 134, "alleles": ("T", "G")}])
        [v] = load_vcf_variants(vcf, "chr1:101-200")
        assert (v["vaf"], v["depth"], v["alt_count"]) == (0.0, 0, 0)

    @pytest.mark.unit
    def test_missing_af_entry_is_skipped_not_crashing(self, tmp_path):
        """AF=.,0.3 arrives from pysam as (None, 0.3) — float(None) used to TypeError."""
        vcf = _write(
            tmp_path,
            _header(info=("DP", "AF")),
            [
                {
                    "contig": "chr1",
                    "pos": 134,
                    "alleles": ("T", "G", "C"),
                    "info": {"DP": 10, "AF": (None, 0.3)},
                }
            ],
        )
        g, c = load_vcf_variants(vcf, "chr1:101-200")
        assert (g["alt"], g["vaf"], g["alt_count"]) == ("G", 0.0, 0)
        assert (c["alt"], c["vaf"], c["alt_count"]) == ("C", 0.3, 3)


class TestAlleleClassification:
    @pytest.mark.unit
    def test_star_allele_is_not_an_snv(self, tmp_path):
        """'*' (spanning deletion) is carried by the upstream deletion record, not here."""
        vcf = _write(
            tmp_path,
            _header(info=("AF",)),
            [
                {
                    "contig": "chr1",
                    "pos": 134,
                    "alleles": ("T", "G", "*"),
                    "info": {"AF": (0.4, 0.2)},
                }
            ],
        )
        variants = load_vcf_variants(vcf, "chr1:101-200")
        assert [v["alt"] for v in variants] == ["G"]

    @pytest.mark.unit
    def test_mnp_is_labelled_mnp_not_indel(self, tmp_path):
        vcf = _write(tmp_path, _header(), [{"contig": "chr1", "pos": 134, "alleles": ("AC", "GT")}])
        [v] = load_vcf_variants(vcf, "chr1:101-200")
        assert v["variant_kind"] == "mnp"
        assert (v["ref"], v["alt"]) == ("AC", "GT")

    @pytest.mark.unit
    def test_multiallelic_alleles_trimmed_independently(self, tmp_path):
        """REF=ACG ALT=ATG,A: the first ALT is really C>T one base later (an SNV that
        must dedup against BAMCP's own candidate); the second is a deletion."""
        vcf = _write(
            tmp_path, _header(), [{"contig": "chr1", "pos": 134, "alleles": ("ACG", "ATG", "A")}]
        )
        snv, dele = load_vcf_variants(vcf, "chr1:101-200")
        assert (snv["position"], snv["ref"], snv["alt"], snv["variant_kind"]) == (
            134,
            "C",
            "T",
            "snv",
        )
        assert (dele["position"], dele["ref"], dele["alt"], dele["variant_kind"]) == (
            133,
            "ACG",
            "A",
            "indel",
        )

    @pytest.mark.unit
    def test_padded_insertion_trimmed_to_anchor(self, tmp_path):
        vcf = _write(
            tmp_path, _header(), [{"contig": "chr1", "pos": 134, "alleles": ("ATT", "ATTT", "AT")}]
        )
        ins, dele = load_vcf_variants(vcf, "chr1:101-200")
        assert (ins["position"], ins["ref"], ins["alt"]) == (133, "A", "AT")
        assert (dele["position"], dele["ref"], dele["alt"]) == (133, "AT", "A")


class TestSampleSupport:
    @pytest.mark.unit
    def test_single_sample_format_ad_dp_preferred_over_info(self, tmp_path):
        vcf = _write(
            tmp_path,
            _header(info=("DP", "AF"), fmt=("GT", "AD", "DP"), samples=("S1",)),
            [
                {
                    "contig": "chr1",
                    "pos": 134,
                    "alleles": ("T", "G"),
                    # Cohort-level INFO disagrees with the sample's own read counts.
                    "info": {"DP": 500, "AF": 0.01},
                    "samples": {"S1": {"GT": (0, 1), "AD": (6, 4), "DP": 10}},
                }
            ],
        )
        [v] = load_vcf_variants(vcf, "chr1:101-200")
        assert (v["depth"], v["alt_count"], v["vaf"]) == (10, 4, 0.4)
        assert v["samples"]["S1"]["AD"] == [6, 4]

    @pytest.mark.unit
    def test_ad_without_format_dp_sums_ad(self, tmp_path):
        vcf = _write(
            tmp_path,
            _header(fmt=("AD",), samples=("S1",)),
            [
                {
                    "contig": "chr1",
                    "pos": 134,
                    "alleles": ("T", "G", "C"),
                    "samples": {"S1": {"AD": (5, 3, 2)}},
                }
            ],
        )
        g, c = load_vcf_variants(vcf, "chr1:101-200")
        assert (g["depth"], g["alt_count"], g["vaf"]) == (10, 3, 0.3)
        assert (c["depth"], c["alt_count"], c["vaf"]) == (10, 2, 0.2)

    @pytest.mark.unit
    def test_multisample_falls_back_to_info(self, tmp_path):
        """With several samples, which one matches the BAM is unknowable → INFO."""
        vcf = _write(
            tmp_path,
            _header(info=("DP", "AF"), fmt=("GT", "AD"), samples=("tumor", "normal")),
            [
                {
                    "contig": "chr1",
                    "pos": 134,
                    "alleles": ("T", "G"),
                    "info": {"DP": 40, "AF": 0.25},
                    "samples": {
                        "tumor": {"GT": (0, 1), "AD": (10, 10)},
                        "normal": {"GT": (0, 0), "AD": (20, 0)},
                    },
                }
            ],
        )
        [v] = load_vcf_variants(vcf, "chr1:101-200")
        assert (v["depth"], v["alt_count"], v["vaf"]) == (40, 10, 0.25)
        assert v["sample_names"] == ["tumor", "normal"]
        assert v["samples"]["tumor"]["GT"] == [0, 1]


class TestStructural:
    @pytest.mark.unit
    def test_multisample_structural_variant_metadata(self, tmp_path):
        vcf = _write(
            tmp_path,
            _header(
                info=("DP", "AF", "SVTYPE", "END", "SVLEN"),
                fmt=("GT", "AD"),
                samples=("tumor", "normal"),
            ),
            [
                {
                    "contig": "chr2",
                    "pos": 501,
                    "alleles": ("N", "<DEL>"),
                    "info": {"DP": 42, "AF": 0.5, "SVTYPE": "DEL", "END": 900, "SVLEN": -399},
                    "samples": {
                        "tumor": {"GT": (0, 1), "AD": (20, 22)},
                        "normal": {"GT": (0, 0), "AD": (40, 0)},
                    },
                }
            ],
        )
        [v] = load_vcf_variants(vcf, "chr2:401-1000")
        assert v["variant_kind"] == "sv"
        assert (v["ref"], v["alt"]) == ("N", "<DEL>")  # symbolic alleles are never trimmed
        assert (v["sv_type"], v["sv_end"], v["sv_len"]) == ("DEL", 900, -399)
        assert v["samples"]["tumor"]["AD"] == [20, 22]

    @pytest.mark.unit
    def test_per_alt_svlen_selected_by_alt_index(self, tmp_path):
        vcf = _write(
            tmp_path,
            _header(info=("SVTYPE", "END", "SVLEN")),
            [
                {
                    "contig": "chr3",
                    "pos": 1001,
                    "alleles": ("N", "<DEL>", "<DUP>"),
                    "info": {"SVTYPE": "DEL", "END": 2000, "SVLEN": (-500, 900)},
                }
            ],
        )
        variants = load_vcf_variants(vcf, "chr3:901-3000")
        assert [v["alt"] for v in variants] == ["<DEL>", "<DUP>"]
        assert [v["sv_len"] for v in variants] == [-500, 900]

    @pytest.mark.unit
    def test_symbolic_alt_without_svtype_is_still_sv(self, tmp_path):
        vcf = _write(
            tmp_path, _header(), [{"contig": "chr2", "pos": 501, "alleles": ("N", "<INV>")}]
        )
        [v] = load_vcf_variants(vcf, "chr2:401-1000")
        assert v["variant_kind"] == "sv"


class TestRegionQuery:
    @pytest.mark.unit
    def test_one_based_region_bounds_are_inclusive(self, tmp_path):
        vcf = _write(
            tmp_path,
            _header(),
            [
                {"contig": "chr1", "pos": 100, "alleles": ("A", "C")},
                {"contig": "chr1", "pos": 101, "alleles": ("A", "C")},
                {"contig": "chr1", "pos": 200, "alleles": ("A", "C")},
                {"contig": "chr1", "pos": 201, "alleles": ("A", "C")},
            ],
        )
        positions = [v["position"] + 1 for v in load_vcf_variants(vcf, "chr1:101-200")]
        assert positions == [101, 200]

    @pytest.mark.unit
    def test_unindexed_plain_vcf_is_a_clear_error(self, tmp_path):
        vcf = _write(
            tmp_path,
            _header(),
            [{"contig": "chr1", "pos": 134, "alleles": ("T", "G")}],
            name="calls.vcf",
            index=False,
        )
        with pytest.raises(ValueError, match="bgzip-compressed and indexed"):
            load_vcf_variants(vcf, "chr1:101-200")
