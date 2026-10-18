"""A tiny synthetic "GIAB" for the round-4 selection tests (hermetic: no network).

One 60 kb ``chr20`` with paired reads (100 bp, insert 300) at about 40x and a handful of
planted loci, each standing in for one round-4 stratum. ``build_world(tmp)`` writes the
reference, a "60x Illumina" BAM, a "HiFi" BAM, two parent BAMs, truth VCFs + BEDs, a caller
VCF, stratification BEDs and a ``genomicSuperDups`` table, and returns their paths plus
the planted coordinates.
"""

from __future__ import annotations

import gzip
import random
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

import pysam

CHROM = "chr20"
LEN = 60_000
READ = 100
INSERT = 300
STEP = 5  # one pair every 5 bp per strand pattern -> ~40x

# Planted loci (0-based). Each comment names the stratum it feeds.
HP_POS = 5000  # 12 x A run [5000, 5012): truth 1 bp deletion, het      (tr_indel_pos)
HP_NEG = 15000  # 12 x T run [15000, 15012): no truth, -1 stutter     (tr_stutter_neg)
MOSAIC = 25000  # truth mosaic SNV, ~12% alt                          (mosaic_snv_pos)
FP = 30000  # caller-flagged, ~12% alt, not in mosaic VCF            (mosaic_caller_fp_neg)
FP_NEAR_TRUTH = 31000  # caller-flagged but a v5.0q record 10 bp away (excluded)
ADJ_DEL = 35000  # truth 3 bp deletion (anchor 35000, deletes 35001-35003)
ADJ_SNV = 35009  # truth het SNV, gap 6 from the deletion span       (indel_adjacent_snv_pos)
ADJ_NOISE = 34994  # misalignment noise, gap 6 upstream               (indel_adjacent_misalign_neg)
CMRG_SNV = 21000  # CMRG het SNV in segdup, outside v4.2.1 BED       (segdup_cmrg_pos)
PSV = 21500  # paralog-like noise, ~15% alt, in DRAGEN VCF            (segdup_psv_neg)
ANCHOR_HET = 10000  # clean het SNV in both sets                      (anchor_het_snv)
ANCHOR_HOM = 12000  # clean hom SNV in both sets                      (anchor_hom_snv)
SEGDUP = (20000, 23000)
SEGDUP_PARTNER = (45000, 48000)
V421_BED = [(1000, 20000), (23000, 59000)]  # excludes the segdup
V5Q_BED = [(1000, 59000)]
MOSAIC_BED = [(1000, 59000)]
CMRG_BED = [(19000, 24000)]


def _ref() -> str:
    rng = random.Random(11)  # noqa: S311
    seq: list[str] = []
    while len(seq) < LEN:
        b = rng.choice("ACGT")
        if len(seq) >= 2 and seq[-1] == seq[-2] == b:
            continue
        seq.append(b)
    s = "".join(seq)
    # Homopolymers: flank with non-matching bases so runs are exactly 12.
    s = s[: HP_POS - 1] + "C" + "A" * 12 + "G" + s[HP_POS + 13 :]
    s = s[: HP_NEG - 1] + "G" + "T" * 12 + "C" + s[HP_NEG + 13 :]
    return s


def other(b: str) -> str:
    return {"A": "C", "C": "G", "G": "T", "T": "A"}[b]


@dataclass
class Edit:
    """One haplotype change applied to every read drawn from that haplotype."""

    kind: str  # "snv" | "del" | "ins"
    pos: int
    data: Any  # base / deleted length / inserted bases


@dataclass
class Locus:
    """Reads covering ``pos`` get ``edits`` when ``pick(i)`` says so (i = read index)."""

    pos: int
    edits: list[Edit]
    pick: Any
    low_bq: bool = False
    mapq: int | None = None


def _build_read(ref: str, start: int, edits: list[Edit]) -> tuple[str, list[tuple[int, int]]]:
    """Sequence + CIGAR of a read starting at ``start`` spanning READ reference bases."""
    seq: list[str] = []
    cig: list[tuple[int, int]] = []

    def push(op: int, n: int) -> None:
        if cig and cig[-1][0] == op:
            cig[-1] = (op, cig[-1][1] + n)
        else:
            cig.append((op, n))

    by_pos = {e.pos: e for e in edits}
    p = start
    end = start + READ
    while p < end:
        e = by_pos.get(p)
        if e is not None and e.kind == "del" and p + e.data <= end:
            push(2, e.data)
            p += e.data
            continue
        if e is not None and e.kind == "ins":
            seq.extend(e.data)
            push(1, len(e.data))
        base = ref[p]
        if e is not None and e.kind == "snv":
            base = e.data
        seq.append(base)
        push(0, 1)
        p += 1
    return "".join(seq), cig


@dataclass
class World:
    d: Path
    ref: str
    paths: dict[str, Path] = field(default_factory=dict)


def _write_bam(
    path: Path, ref: str, loci: list[Locus], *, seed: int, pair: bool = True, depth_step: int = STEP
) -> None:
    rng = random.Random(seed)  # noqa: S311
    header = {"HD": {"VN": "1.6", "SO": "coordinate"}, "SQ": [{"SN": CHROM, "LN": LEN}]}
    rows: list[tuple[int, str, int, str, list[tuple[int, int]], int, int, int, list[int]]] = []
    for n0, start in enumerate(range(500, LEN - 1000, depth_step)):
        n = n0 + 1
        name = f"p{n0}"
        s1 = start
        s2 = start + INSERT - READ
        mates = [(s1, False, True), (s2, True, False)] if pair else [(s1, n % 2 == 1, True)]
        mapqs = []
        built = []
        for s, rev, r1 in mates:
            edits: list[Edit] = []
            mq = 60
            low = False
            for loc in loci:
                if s <= loc.pos < s + READ and loc.pick(n):
                    edits.extend(loc.edits)
                    low = low or loc.low_bq
                if s <= loc.pos < s + READ and loc.mapq is not None:
                    mq = loc.mapq
            seq, cig = _build_read(ref, s, edits)
            quals = [35] * len(seq)
            if low:
                quals = [12] * len(seq)
            built.append((s, rev, r1, seq, cig, quals))
            mapqs.append(mq)
        for k, (s, rev, r1, seq, cig, quals) in enumerate(built):
            flag = 0
            if pair:
                flag = 0x1 | 0x2 | (0x40 if r1 else 0x80)
                if rev:
                    flag |= 0x10
                else:
                    flag |= 0x20
                mate = built[1 - k]
                rows.append((s, name, flag, seq, cig, mapqs[k], mate[0], 0, quals))
            else:
                rows.append(
                    (s, name, flag | (0x10 if rev else 0), seq, cig, mapqs[k], -1, 0, quals)
                )
    _ = rng
    rows.sort(key=lambda t: (t[0], t[1]))
    with pysam.AlignmentFile(str(path), "wb", header=header) as out:
        for s, name, flag, seq, cig, mq, mstart, _t, quals in rows:
            a = pysam.AlignedSegment(out.header)
            a.query_name = name
            a.query_sequence = seq
            a.flag = flag
            a.reference_id = 0
            a.reference_start = s
            a.mapping_quality = mq
            a.cigartuples = cig
            a.query_qualities = pysam.qualitystring_to_array("".join(chr(q + 33) for q in quals))
            if flag & 0x1:
                a.next_reference_id = 0
                a.next_reference_start = mstart
                a.template_length = INSERT if not flag & 0x10 else -INSERT
            out.write(a)
    pysam.index(str(path))


def _vcf(path: Path, records: list[tuple[int, str, str, str, str]], *, index: bool = True) -> Path:
    """records: (POS 1-based, REF, ALT, FILTER, GT)."""
    body = "".join(
        f"{CHROM}\t{p}\t.\t{r}\t{a}\t50\t{f}\t.\tGT\t{g}\n" for p, r, a, f, g in sorted(records)
    )
    txt = (
        "##fileformat=VCFv4.2\n"
        f"##contig=<ID={CHROM},length={LEN}>\n"
        '##FILTER=<ID=HET1,Description="x">\n'
        '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n'
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tHG002\n" + body
    )
    if not index:
        path.write_text(txt)
        return path
    plain = path.with_suffix("")
    plain.write_text(txt)
    pysam.tabix_compress(str(plain), str(path), force=True)
    pysam.tabix_index(str(path), preset="vcf", force=True)
    return path


def _bed(path: Path, ivs: list[tuple[int, int]], chrom: str = CHROM) -> Path:
    opener: Any = gzip.open if path.suffix == ".gz" else open
    with opener(path, "wt") as fh:
        for s, e in ivs:
            fh.write(f"{chrom}\t{s}\t{e}\n")
    return path


def build_world(d: Path) -> World:
    ref = _ref()
    w = World(d, ref)
    fa = d / "ref.fa"
    fa.write_text(f">{CHROM}\n" + "\n".join(ref[i : i + 80] for i in range(0, LEN, 80)) + "\n")
    pysam.faidx(str(fa))
    w.paths["ref"] = fa

    snv = lambda pos: Edit("snv", pos, other(ref[pos]))  # noqa: E731
    illumina = [
        # TR pos: het 1 bp deletion of the A run (left-aligned at its first base), plus
        # 2 bp stutter on some reads.
        Locus(HP_POS, [Edit("del", HP_POS, 1)], lambda i: i % 2 == 0),
        Locus(HP_POS, [Edit("del", HP_POS, 2)], lambda i: i % 7 == 1),
        # TR neg: -1 stutter on ~1/6 of reads.
        Locus(HP_NEG, [Edit("del", HP_NEG, 1)], lambda i: i % 6 == 1),
        Locus(MOSAIC, [snv(MOSAIC)], lambda i: i % 8 == 3),
        Locus(FP, [snv(FP)], lambda i: i % 8 == 5, low_bq=True),
        Locus(FP_NEAR_TRUTH, [snv(FP_NEAR_TRUTH)], lambda i: i % 8 == 5),
        Locus(ADJ_DEL + 1, [Edit("del", ADJ_DEL + 1, 3)], lambda i: i % 2 == 0),
        Locus(ADJ_SNV, [snv(ADJ_SNV)], lambda i: i % 2 == 1),
        Locus(ADJ_NOISE, [snv(ADJ_NOISE)], lambda i: i % 5 == 2),
        Locus(CMRG_SNV, [snv(CMRG_SNV)], lambda i: i % 2 == 0, mapq=15),
        Locus(PSV, [snv(PSV)], lambda i: i % 6 == 0, mapq=30),
        Locus(ANCHOR_HET, [snv(ANCHOR_HET)], lambda i: i % 2 == 0),
        Locus(ANCHOR_HOM, [snv(ANCHOR_HOM)], lambda i: True),
    ]
    _write_bam(d / "ill.bam", ref, illumina, seed=1)
    w.paths["bam"] = d / "ill.bam"
    hifi = [
        Locus(HP_POS, [Edit("del", HP_POS, 1)], lambda i: i % 2 == 0),
        Locus(MOSAIC, [snv(MOSAIC)], lambda i: i % 8 == 3),
        Locus(ADJ_DEL + 1, [Edit("del", ADJ_DEL + 1, 3)], lambda i: i % 2 == 0),
        Locus(ADJ_SNV, [snv(ADJ_SNV)], lambda i: i % 2 == 1),
        Locus(CMRG_SNV, [snv(CMRG_SNV)], lambda i: i % 2 == 0),
        Locus(ANCHOR_HET, [snv(ANCHOR_HET)], lambda i: i % 2 == 0),
        Locus(ANCHOR_HOM, [snv(ANCHOR_HOM)], lambda i: True),
    ]
    _write_bam(d / "hifi.bam", ref, hifi, seed=2, pair=False)
    w.paths["hifi"] = d / "hifi.bam"
    # Parents: both carry the TR stutter and the PSV signal (locus-specific artifacts).
    parent = [
        Locus(HP_NEG, [Edit("del", HP_NEG, 1)], lambda i: i % 6 == 1),
        Locus(PSV, [snv(PSV)], lambda i: i % 6 == 0),
    ]
    for nm, seed in (("HG003", 3), ("HG004", 4)):
        _write_bam(d / f"{nm}.bam", ref, parent if nm == "HG003" else [], seed=seed)
        w.paths[nm] = d / f"{nm}.bam"

    def rec(pos0: int, alt: str | None = None, gt: str = "0/1", filt: str = "PASS"):
        return (pos0 + 1, ref[pos0], alt or other(ref[pos0]), filt, gt)

    hp_del = (HP_POS, ref[HP_POS - 1 : HP_POS + 1], ref[HP_POS - 1], "PASS", "0/1")
    adj_del = (ADJ_DEL + 1, ref[ADJ_DEL : ADJ_DEL + 4], ref[ADJ_DEL], "PASS", "0/1")
    germ = [
        hp_del,
        adj_del,
        rec(ADJ_SNV),
        rec(ANCHOR_HET),
        rec(ANCHOR_HOM, gt="1/1"),
    ]
    _vcf(d / "v421.vcf.gz", germ)
    _bed(d / "v421.bed", V421_BED)
    v5q = [*germ, rec(FP_NEAR_TRUTH + 10), rec(40000, filt="HET1"), rec(CMRG_SNV)]
    _vcf(d / "v5q.vcf.gz", v5q)
    _bed(d / "v5q.bed", V5Q_BED)
    _vcf(d / "cmrg.vcf", [rec(CMRG_SNV)], index=False)
    _bed(d / "cmrg.bed", CMRG_BED)
    _vcf(d / "mosaic.vcf.gz", [rec(MOSAIC)])
    _bed(d / "mosaic.bed", MOSAIC_BED)
    _vcf(d / "caller.vcf.gz", [rec(FP), rec(FP_NEAR_TRUTH), rec(MOSAIC)])
    _vcf(d / "dragen.vcf.gz", [rec(PSV), rec(CMRG_SNV)])
    _vcf(d / "sv.vcf.gz", [])
    for k, v in {
        "v421": "v421.vcf.gz",
        "v421_bed": "v421.bed",
        "v5q": "v5q.vcf.gz",
        "v5q_bed": "v5q.bed",
        "cmrg": "cmrg.vcf",
        "cmrg_bed": "cmrg.bed",
        "mosaic": "mosaic.vcf.gz",
        "mosaic_bed": "mosaic.bed",
        "caller": "caller.vcf.gz",
        "dragen": "dragen.vcf.gz",
        "sv": "sv.vcf.gz",
    }.items():
        w.paths[k] = d / v
    strat = d / "strat"
    strat.mkdir()
    w.paths["strat"] = strat
    w.paths["segdups"] = _bed(strat / "segdups.bed.gz", [SEGDUP, SEGDUP_PARTNER])
    w.paths["lowmappability"] = _bed(strat / "lowmap.bed.gz", [])
    w.paths["homopolymer_ge12"] = _bed(
        strat / "hp12.bed.gz", [(HP_POS - 5, HP_POS + 17), (HP_NEG - 5, HP_NEG + 17)]
    )
    w.paths["homopolymer_7to11"] = _bed(strat / "hp7.bed.gz", [])
    w.paths["diTR"] = _bed(strat / "ditr.bed.gz", [])
    sd = d / "genomicSuperDups.txt.gz"
    with gzip.open(sd, "wt") as fh:
        fh.write(
            f"585\t{CHROM}\t{SEGDUP[0]}\t{SEGDUP[1]}\tx\t0\t+\t{CHROM}\t"
            f"{SEGDUP_PARTNER[0]}\t{SEGDUP_PARTNER[1]}\n"
        )
    w.paths["superdups"] = sd
    return w


def resources(w: World, *, hifi: bool = True, parents: bool = True) -> Any:
    """An :class:`R4Resources` over the world (no network: every path is local)."""
    from bamcp.eval.visual import real

    p = w.paths
    return real.open_r4_resources(
        chroms=[CHROM],
        truth_chroms=[CHROM],
        bam_path=str(p["bam"]),
        ref_path=str(p["ref"]),
        v421=(str(p["v421"]), str(p["v421_bed"])),
        v5q=(str(p["v5q"]), str(p["v5q_bed"])),
        cmrg=(str(p["cmrg"]), str(p["cmrg_bed"])),
        mosaic=(str(p["mosaic"]), str(p["mosaic_bed"])),
        strat={
            k: str(p[k])
            for k in ("segdups", "lowmappability", "homopolymer_ge12", "homopolymer_7to11", "diTR")
        }
        | {"complexindel10bp": str(p["diTR"]), "snpswithin10bp": str(p["diTR"])},
        superdups=str(p["superdups"]),
        sv_vcfs=[str(p["sv"])],
        hifi=(str(p["hifi"]), None) if hifi else None,
        parents={n: (str(p[n]), None) for n in ("HG003", "HG004")} if parents else None,
        log=lambda *_: None,
    )
