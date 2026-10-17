"""Real-data tier of the visual-vs-text eval: GIAB HG002 sites on re-indexed coordinates.

Items are real Illumina reads at sites whose truth comes from the GIAB v4.2.1 benchmark
(``docs/eval/visual-preregistration-r3.md``). Each site's reads and reference are cut out of
the public data and written into a small local BAM on re-indexed contigs (``chr1``..``chr3``,
one 3 kb slot per item, the same layout as the synthetic tier), so a model cannot look up a
known coordinate, and every condition reads the same local files.

Strata (``case_type``):

- ``real_snv_het`` / ``real_snv_hom`` — truth SNVs outside tandem repeats and homopolymers.
- ``real_snv_in_repeat`` — truth SNVs inside the GIAB tandem-repeat/homopolymer
  stratification.
- ``real_indel_het`` / ``real_indel_hom`` — truth indels (1-10 bp) in any context (real
  indels sit almost entirely in repeats).
- ``real_noisy_negative`` — no truth variant within 10 bp, but >= 3 primary non-duplicate
  reads (>= 5%) carry the same non-reference allele: what a naive caller would flag.
- ``real_clean_negative`` — no truth variant within 10 bp and no non-reference read.

Truth variants are kept only when no other truth variant lies within 30 bp, so the position
of interest has one unambiguous answer. Negatives' ``artifact_type`` is ``unknown`` (the
truth set says the sample is reference there, not why reads disagree) and is not scored.
"""

from __future__ import annotations

import bisect
import gzip
import hashlib
import json
import random
from collections import Counter
from collections.abc import Callable, Iterable
from dataclasses import asdict, dataclass, field
from pathlib import Path
from typing import Any

from .fixtures import Item, _derive_rng, opaque_id
from .stimuli import SAMPLE_LINE_R4

GIAB_BASE = (
    "https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/AshkenazimTrio/"
    "HG002_NA24385_son/NISTv4.2.1/GRCh38/"
)
TRUTH_VCF = GIAB_BASE + "HG002_GRCh38_1_22_v4.2.1_benchmark.vcf.gz"
CONFIDENT_BED = GIAB_BASE + "HG002_GRCh38_1_22_v4.2.1_benchmark_noinconsistent.bed"
REPEAT_BED = (
    "https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/genome-stratifications/"
    "v3.3/GRCh38@all/LowComplexity/GRCh38_AllTandemRepeatsandHomopolymers_slop5.bed.gz"
)
HG002_BAM = (
    "https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/data/AshkenazimTrio/"
    "HG002_NA24385_son/NIST_HiSeq_HG002_Homogeneity-10953946/"
    "NHGRI_Illumina300X_AJtrio_novoalign_bams/HG002.GRCh38.60x.1.bam"
)
GRCH38_FASTA = (
    "https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/GRCh38_reference_genome/"
    "GRCh38_full_analysis_set_plus_decoy_hla.fa"
)

QUOTAS: dict[str, int] = {
    "real_snv_het": 6,
    "real_snv_hom": 6,
    "real_snv_in_repeat": 4,
    "real_indel_het": 6,
    "real_indel_hom": 4,
    "real_noisy_negative": 8,
    "real_clean_negative": 2,
}
REAL_CASE_TYPES = tuple(QUOTAS)
REAL_BAM_NAME = "vt_real.bam"
REAL_REF_NAME = "vt_real_ref.fa"
WINDOW_BP = 2000
ISOLATION_BP = 30
NEG_CLEAR_BP = 10
MAX_INDEL = 10
NOISY_MIN_READS = 3
NOISY_MIN_FRAC = 0.05
CLEAN_MIN_DEPTH = 20
N_CONTIGS = 3
SLOT_BP = 3000
CONTIG_PAD = 1000
READ_HALF_SPAN = 450  # keep reads overlapping site +/- this many bp
MIN_SPACING = 1000  # between selected sites; reads near two sites go to the first only
SAMPLE_LINE = "sequenced human sample (Illumina short reads); coordinates are re-indexed"


@dataclass
class RealSite:
    """One selected site. ``site0`` is the 0-based position of interest on ``chrom``."""

    case_type: str
    chrom: str
    site0: int
    label: dict[str, Any]
    in_repeat: bool
    vcf_pos: int | None = None
    vcf_ref: str = ""
    vcf_alt: str = ""
    alt_reads: int | None = None
    depth: int | None = None
    extra: dict[str, Any] = field(default_factory=dict)


# ---------------------------------------------------------------------------
# Intervals
# ---------------------------------------------------------------------------


class Intervals:
    """Sorted, merged half-open intervals on one contig with O(log n) lookups."""

    def __init__(self, ivs: Iterable[tuple[int, int]]):
        merged: list[list[int]] = []
        for s, e in sorted(ivs):
            if merged and s <= merged[-1][1]:
                merged[-1][1] = max(merged[-1][1], e)
            else:
                merged.append([s, e])
        self.starts = [s for s, _ in merged]
        self.ends = [e for _, e in merged]

    def __len__(self) -> int:
        return len(self.starts)

    def contains(self, pos0: int) -> bool:
        i = bisect.bisect_right(self.starts, pos0) - 1
        return i >= 0 and pos0 < self.ends[i]

    def overlaps(self, s: int, e: int) -> bool:
        i = bisect.bisect_right(self.starts, e - 1) - 1
        return i >= 0 and self.ends[i] > s

    def covers(self, s: int, e: int) -> bool:
        i = bisect.bisect_right(self.starts, s) - 1
        return i >= 0 and self.ends[i] >= e

    def pairs(self) -> list[tuple[int, int]]:
        return list(zip(self.starts, self.ends, strict=True))


def load_bed(path: str | Path, chrom: str) -> Intervals:
    p = Path(path)
    opener = gzip.open if p.suffix == ".gz" else open
    ivs: list[tuple[int, int]] = []
    with opener(p, "rt") as fh:
        for line in fh:
            if not line.strip() or line.startswith(("#", "track", "browser")):
                continue
            f = line.split("\t")
            if f[0] == chrom:
                ivs.append((int(f[1]), int(f[2])))
    return Intervals(ivs)


# ---------------------------------------------------------------------------
# Truth records
# ---------------------------------------------------------------------------


def classify_record(ref: str, alt: str, gt: tuple[int | None, ...]) -> dict[str, Any] | None:
    """Label + 0-based offset of the position of interest from a biallelic truth record.

    Conventions match the synthetic tier: a deletion's position of interest is its first
    deleted base and its ``alt_allele`` the deleted bases; an insertion's is the anchor base
    and its ``alt_allele`` the inserted bases. Returns ``None`` for anything else
    (multi-allelic genotypes, MNVs/complex alleles, indels over ``MAX_INDEL``).
    """
    calls = sorted(a for a in gt if a is not None)
    if len(gt) != 2 or len(calls) != 2 or calls not in ([0, 1], [1, 1]):
        return None
    zyg = "hom_alt" if gt[0] == gt[1] == 1 else "het"
    if len(ref) == len(alt) == 1:
        vtype, allele, offset = "SNV", alt, 0
    elif len(ref) > len(alt) == 1 and ref[0] == alt[0] and len(ref) - 1 <= MAX_INDEL:
        vtype, allele, offset = "del", ref[1:], 1
    elif len(alt) > len(ref) == 1 and ref[0] == alt[0] and len(alt) - 1 <= MAX_INDEL:
        vtype, allele, offset = "ins", alt[1:], 0
    else:
        return None
    return {
        "offset": offset,
        "label": {
            "variant_present": "yes",
            "variant_type": vtype,
            "alt_allele": allele,
            "zygosity": zyg,
            "artifact_verdict": "real",
            "artifact_type": "none",
            "phase": "na",
        },
    }


NEGATIVE_LABEL = {
    "variant_present": "no",
    "variant_type": "none",
    "alt_allele": "",
    "zygosity": "none",
    "artifact_verdict": "artifact",
    "phase": "na",
}


def _usable(read: Any) -> bool:
    return not (
        read.is_unmapped
        or read.is_secondary
        or read.is_supplementary
        or read.is_duplicate
        or read.is_qcfail
    )


# ---------------------------------------------------------------------------
# Scanning
# ---------------------------------------------------------------------------


def _allele_counts(bam: Any, chrom: str, a: int, b: int) -> dict[int, Counter[str]]:
    """Per-position allele counts (base, ``-`` deletion, ``X+SEQ`` insertion) over [a, b)."""
    counts: dict[int, Counter[str]] = {}
    for col in bam.pileup(
        chrom,
        a,
        b,
        truncate=True,
        min_base_quality=0,
        stepper="nofilter",
        ignore_orphans=False,
        ignore_overlaps=False,
        max_depth=100_000,
    ):
        c: Counter[str] = Counter()
        for pr in col.pileups:
            if not _usable(pr.alignment):
                continue
            if pr.is_del or pr.is_refskip:
                c["-"] += 1
                continue
            q = pr.alignment.query_sequence
            base = q[pr.query_position] if q else "N"
            if pr.indel > 0:  # a following deletion is counted at the deleted positions
                base += "+" + q[pr.query_position + 1 : pr.query_position + 1 + pr.indel]
            c[base.upper()] += 1
        counts[col.reference_pos] = c
    return counts


def scan_window(
    bam: Any,
    fasta: Any,
    vcf: Any,
    chrom: str,
    a: int,
    b: int,
    confident: Intervals,
    repeats: Intervals,
) -> dict[str, list[RealSite]]:
    """Candidate sites in [a, b) for every stratum (truth from ``vcf``, signal from ``bam``)."""
    pools: dict[str, list[RealSite]] = {k: [] for k in QUOTAS}
    truth = list(vcf.fetch(chrom, max(0, a - 100), b + 100))
    # Whole reference spans, not anchors: a long deletion anchored far upstream still covers
    # (and makes non-reference reads at) every base it deletes.
    spans = [(r.pos - 1, r.pos - 1 + max(len(r.ref), 1)) for r in truth]

    def near_truth(s0: int, e0: int, bp: int, exclude: int | None = None) -> bool:
        """Does any truth span (other than record index ``exclude``) come within ``bp`` of
        [s0, e0)? By index, so a second record at the same POS still counts."""
        return any(
            start < e0 + bp and end > s0 - bp
            for i, (start, end) in enumerate(spans)
            if i != exclude
        )

    for idx, rec in enumerate(truth):
        p0 = rec.pos - 1
        if not (a <= p0 < b) or len(rec.alts or ()) != 1:
            continue
        span_end = p0 + max(len(rec.ref), 1)
        if not confident.covers(p0, span_end) or near_truth(
            p0, span_end, ISOLATION_BP, exclude=idx
        ):
            continue
        gt = next(iter(rec.samples.values()))["GT"]
        cls = classify_record(rec.ref, rec.alts[0], tuple(gt))
        if cls is None:
            continue
        label = cls["label"]
        # The variant's own reference span; an insertion also touches the base after its
        # anchor (the inserted bases sit between the two).
        is_ins = label["variant_type"] == "ins"
        in_rep = repeats.overlaps(p0, span_end + 1 if is_ins else span_end)
        if label["variant_type"] != "SNV":  # real indels sit mostly in repeats: any context
            stratum = "real_indel_het" if label["zygosity"] == "het" else "real_indel_hom"
        elif in_rep:
            stratum = "real_snv_in_repeat"
        else:
            stratum = "real_snv_het" if label["zygosity"] == "het" else "real_snv_hom"
        pools[stratum].append(
            RealSite(
                stratum,
                chrom,
                p0 + cls["offset"],
                dict(label),
                in_rep,
                vcf_pos=rec.pos,
                vcf_ref=rec.ref,
                vcf_alt=rec.alts[0],
            )
        )

    ref_seq = fasta.fetch(chrom, a, b).upper()
    counts = _allele_counts(bam, chrom, a, b)
    for p0, c in counts.items():
        if not confident.covers(p0 - NEG_CLEAR_BP, p0 + NEG_CLEAR_BP + 1):
            continue
        if near_truth(p0, p0 + 1, NEG_CLEAR_BP):
            continue
        ref_base = ref_seq[p0 - a]
        depth = sum(c.values())
        nonref = {k: n for k, n in c.items() if k != ref_base}
        top = max(nonref.items(), key=lambda kv: (kv[1], kv[0]), default=("", 0))
        in_rep = repeats.contains(p0)
        if top[1] >= NOISY_MIN_READS and depth and top[1] / depth >= NOISY_MIN_FRAC:
            pools["real_noisy_negative"].append(
                RealSite(
                    "real_noisy_negative",
                    chrom,
                    p0,
                    {**NEGATIVE_LABEL, "artifact_type": "unknown"},
                    in_rep,
                    alt_reads=top[1],
                    depth=depth,
                    extra={"top_allele": top[0]},
                )
            )
        elif not nonref and depth >= CLEAN_MIN_DEPTH and not in_rep:
            pools["real_clean_negative"].append(
                RealSite(
                    "real_clean_negative",
                    chrom,
                    p0,
                    {**NEGATIVE_LABEL, "artifact_type": "none"},
                    in_rep,
                    alt_reads=0,
                    depth=depth,
                )
            )
    return pools


def _pick(
    pool: list[RealSite], k: int, rng: random.Random, taken: list[RealSite]
) -> list[RealSite]:
    """``k`` sites from ``pool``, each >= ``MIN_SPACING`` bp from every site already taken
    (any stratum), so no two items share reads."""
    rng.shuffle(pool)
    chosen: list[RealSite] = []
    for s in pool:
        if len(chosen) >= k:
            break
        if all(s.chrom != c.chrom or abs(s.site0 - c.site0) >= MIN_SPACING for c in taken):
            chosen.append(s)
            taken.append(s)
    return chosen


def _spaced(pool: list[RealSite]) -> int:
    """How many of ``pool`` can be taken at once, ``MIN_SPACING`` apart (greedy on position):
    adjacent candidates from one locus count once."""
    n, last = 0, None
    for s in sorted(pool, key=lambda s: (s.chrom, s.site0)):
        if last is None or s.chrom != last.chrom or s.site0 - last.site0 >= MIN_SPACING:
            n, last = n + 1, s
    return n


def _pools_full(pools: dict[str, list[RealSite]], quotas: dict[str, int], factor: int) -> bool:
    """Every stratum's pool holds ``factor`` x its quota of mutually spaced candidates; noisy
    negatives need that inside and outside repeats separately (drawn half from each). The
    ``factor`` margin absorbs cross-stratum spacing conflicts; a draw that still falls short
    raises rather than returning a short manifest."""
    for k, q in quotas.items():
        pool = pools.get(k, [])
        if k == "real_noisy_negative":
            inside = [s for s in pool if s.in_repeat]
            outside = [s for s in pool if not s.in_repeat]
            if _spaced(inside) < factor * (q // 2) or _spaced(outside) < factor * (q - q // 2):
                return False
        elif _spaced(pool) < factor * q:
            return False
    return True


def select_sites(
    seed: int,
    *,
    chrom: str = "chr20",
    bam_path: str = HG002_BAM,
    ref_path: str = GRCH38_FASTA,
    vcf_path: str = TRUTH_VCF,
    confident_bed: str | Path,
    repeat_bed: str | Path,
    quotas: dict[str, int] | None = None,
    max_windows: int = 600,
    pool_factor: int = 3,
    log: Any = print,
) -> list[RealSite]:
    """Scan seeded random windows of ``chrom`` until every stratum's pool holds
    ``pool_factor`` x its quota (or ``max_windows`` is reached), then draw each quota with
    the same seed. Noisy negatives are drawn half inside repeats, half outside."""
    import pysam

    quotas = dict(quotas or QUOTAS)
    rng = _derive_rng(seed, "real-sites")
    confident = load_bed(confident_bed, chrom)
    repeats = load_bed(repeat_bed, chrom)
    big = [(s, e) for s, e in confident.pairs() if e - s >= WINDOW_BP]
    if not big:
        raise ValueError(f"no confident interval >= {WINDOW_BP} bp on {chrom}")
    pools: dict[str, list[RealSite]] = {k: [] for k in quotas}
    used: list[int] = []
    with (
        pysam.AlignmentFile(bam_path) as bam,
        pysam.FastaFile(ref_path) as fasta,
        pysam.VariantFile(vcf_path) as vcf,
    ):
        # Only scanned windows count toward max_windows; a start that overlaps an earlier
        # window is redrawn (attempts are capped so a crowded contig cannot loop forever).
        scanned = attempts = 0
        full = False
        while scanned < max_windows and attempts < 20 * max_windows:
            attempts += 1
            s, e = rng.choice(big)
            a = rng.randint(s, e - WINDOW_BP)
            if any(abs(a - u) < WINDOW_BP for u in used):
                continue
            used.append(a)
            scanned += 1
            for k, v in scan_window(
                bam, fasta, vcf, chrom, a, a + WINDOW_BP, confident, repeats
            ).items():
                pools.setdefault(k, []).extend(v)
            if scanned % 50 == 0:
                log(f"{scanned} windows: { {k: len(v) for k, v in pools.items()} }")
            if _pools_full(pools, quotas, pool_factor):
                log(f"pools full after {scanned} windows ({attempts} draws)")
                full = True
                break
        if not full:
            sizes = {k: len(v) for k, v in pools.items()}
            log(f"scanned {scanned} windows in {attempts} draws; pool sizes: {sizes}")

    chosen: list[RealSite] = []
    pick_rng = _derive_rng(seed, "real-pick")
    for k, q in quotas.items():
        pool = sorted(pools.get(k, []), key=lambda s: (s.chrom, s.site0))
        if k == "real_noisy_negative":
            rep = [s for s in pool if s.in_repeat]
            non = [s for s in pool if not s.in_repeat]
            in_rep = _pick(rep, q // 2, pick_rng, chosen)
            out_rep = _pick(non, q - q // 2, pick_rng, chosen)
            if len(in_rep) < q // 2 or len(out_rep) < q - q // 2:
                raise ValueError(
                    f"stratum {k}: need {q // 2} inside repeats and {q - q // 2} outside, "
                    f"found {len(in_rep)} and {len(out_rep)}"
                )
            got = in_rep + out_rep
        else:
            got = _pick(pool, q, pick_rng, chosen)
        if len(got) < q:
            raise ValueError(f"stratum {k}: only {len(got)} of {q} sites available")
    return sorted(chosen, key=lambda s: (s.case_type, s.chrom, s.site0))


def sites_to_jsonl(sites: list[RealSite], path: str | Path) -> str:
    """Write the frozen site manifest; returns its sha256 (recorded in the pre-registration)."""
    text = "".join(json.dumps(asdict(s), sort_keys=True) + "\n" for s in sites)
    Path(path).write_text(text)
    return hashlib.sha256(text.encode()).hexdigest()


def sites_from_jsonl(path: str | Path) -> list[RealSite]:
    return [RealSite(**json.loads(line)) for line in Path(path).read_text().splitlines() if line]


# ---------------------------------------------------------------------------
# Blinded fixtures
# ---------------------------------------------------------------------------


def _layout(sites: list[RealSite], seed: int) -> list[tuple[str, int]]:
    """A seeded (contig, slot-start) for every site, shuffled so order carries no stratum."""
    per = -(-len(sites) // N_CONTIGS)
    slots = [
        (f"chr{c + 1}", CONTIG_PAD + k * SLOT_BP) for c in range(N_CONTIGS) for k in range(per)
    ]
    _derive_rng(seed, "real-layout").shuffle(slots)
    return slots[: len(sites)]


MATE_SENTINEL = "chrU"  # "mate elsewhere": carries no real coordinate
MATE_SENTINEL_LEN = 1000
# Paired, proper pair, mate unmapped, reverse, mate reverse, first / second in pair.
PAIR_FLAG_MASK = 0x1 | 0x2 | 0x8 | 0x10 | 0x20 | 0x40 | 0x80


@dataclass
class _OutRead:
    tid: int
    start: int
    name: str
    flag: int
    seq: str
    quals: Any
    cigar: Any
    mapq: int
    mate_tid: int | None = None
    mate_start: int = -1
    tlen: int = 0
    mate_mapq: int | None = None


def _pair_fields(
    rec: _OutRead, r: Any, chrom: str, src0: int, shift: int, tid: int, sentinel: int
) -> None:
    """Copy a source read's pairing into ``rec`` (``write_real_fixtures(keep_pairs=True)``)."""
    rec.flag = r.flag & PAIR_FLAG_MASK
    if not r.is_paired:
        return
    if r.mate_is_unmapped:  # SAM convention: an unmapped mate sits at the read's position
        rec.mate_tid, rec.mate_start = tid, rec.start
    elif r.next_reference_name == chrom and src0 <= r.next_reference_start < src0 + SLOT_BP:
        rec.mate_tid, rec.mate_start = tid, r.next_reference_start + shift
        rec.tlen = r.template_length
    else:
        rec.mate_tid, rec.mate_start = sentinel, 0


def write_real_fixtures(
    sites: list[RealSite],
    out_dir: str | Path,
    seed: int,
    *,
    bam_path: str = HG002_BAM,
    ref_path: str = GRCH38_FASTA,
    keep_pairs: bool = False,
    bam_index: str | None = None,
) -> tuple[dict[str, Path], list[Item], list[dict[str, Any]]]:
    """Write the blinded BAM + reference and return (paths, items, mapping).

    Kept reads: primary, non-duplicate, non-QC-fail alignments overlapping the site +/-
    ``READ_HALF_SPAN`` whose whole alignment fits the slot. By default they are written
    single-end (flag 0/16), with hashed names and no tags, so every condition counts the
    same reads. ``keep_pairs`` (round 4, viewer V6) also keeps the pairing FLAG bits
    (``PAIR_FLAG_MASK``), the ``MQ`` (mate MAPQ) tag and the mate position: re-indexed when
    the mate starts inside the same slot, else on the sentinel contig ``MATE_SENTINEL``
    (all N, no real coordinate) at position 1 with template length 0. ``MQ`` is the source
    tag, else the mate's MAPQ when the mate is among the fetched reads (Novoalign writes no
    ``MQ``), else absent. ``mapping`` (original
    -> blinded coordinates) is label-bearing and stays out of stim/.
    """
    import pysam

    out = Path(out_dir)
    out.mkdir(parents=True, exist_ok=True)
    layout = _layout(sites, seed)
    per = -(-len(sites) // N_CONTIGS)
    contig_len = 2 * CONTIG_PAD + per * SLOT_BP
    names = [f"chr{c + 1}" for c in range(N_CONTIGS)]
    refs = {n: ["N"] * contig_len for n in names}
    lengths = {n: contig_len for n in names}
    if keep_pairs:
        names.append(MATE_SENTINEL)
        refs[MATE_SENTINEL] = ["N"] * MATE_SENTINEL_LEN
        lengths[MATE_SENTINEL] = MATE_SENTINEL_LEN
    records: list[_OutRead] = []
    items: list[Item] = []
    mapping: list[dict[str, Any]] = []
    # A read is copied into one item only: two nearby sites must never share evidence.
    assigned: set[tuple[str | None, bool, int, bool]] = set()
    with (
        pysam.AlignmentFile(bam_path, index_filename=bam_index) as bam,
        pysam.FastaFile(ref_path) as fasta,
    ):
        for idx, (site, (contig, slot0)) in enumerate(zip(sites, layout, strict=True)):
            src0 = site.site0 - SLOT_BP // 2
            shift = slot0 - src0
            seq = fasta.fetch(site.chrom, src0, src0 + SLOT_BP).upper()
            refs[contig][slot0 : slot0 + len(seq)] = list(seq)
            item_id = "i" + opaque_id(seed, "real", site.chrom, site.site0)
            n = shared = 0
            fetched = list(
                bam.fetch(site.chrom, site.site0 - READ_HALF_SPAN, site.site0 + READ_HALF_SPAN)
            )
            mates = _mate_mapqs(fetched) if keep_pairs else {}
            for r in fetched:
                if not _usable(r) or r.reference_end is None:
                    continue
                if r.reference_start < src0 or r.reference_end > src0 + SLOT_BP:
                    continue
                key = (r.query_name, bool(r.is_read1), r.reference_start, bool(r.is_reverse))
                if key in assigned:
                    shared += 1
                    continue
                assigned.add(key)
                name = hashlib.sha1(f"{seed}:{item_id}:{n}".encode()).hexdigest()[:16]  # noqa: S324
                n += 1
                tid = names.index(contig)
                rec = _OutRead(
                    tid,
                    r.reference_start + shift,
                    name,
                    16 if r.is_reverse else 0,
                    r.query_sequence or "",
                    r.query_qualities,
                    r.cigartuples,
                    r.mapping_quality,
                )
                if keep_pairs:
                    _pair_fields(rec, r, site.chrom, src0, shift, tid, names.index(MATE_SENTINEL))
                    rec.mate_mapq = _mate_mapq(r, mates)
                records.append(rec)
            new_site = site.site0 + shift
            item = Item(
                item_id=item_id,
                case_type=site.case_type,
                instance=idx,
                contig=contig,
                sites=[new_site],
                depth=n,
                label=dict(site.label),
            )
            items.append(item)
            mapping.append(
                {
                    "item": item_id,
                    "case_type": site.case_type,
                    "source": f"{site.chrom}:{site.site0 + 1}",
                    "blinded": f"{contig}:{new_site + 1}",
                    "vcf": [site.vcf_pos, site.vcf_ref, site.vcf_alt],
                    "reads": n,
                    "shared_reads_dropped": shared,
                }
            )

    ref_out = out / REAL_REF_NAME
    with ref_out.open("w") as f:
        for nm in names:
            s = "".join(refs[nm])
            f.write(f">{nm}\n")
            for i in range(0, len(s), 80):
                f.write(s[i : i + 80] + "\n")
    fai = ref_out.with_name(REAL_REF_NAME + ".fai")
    if fai.exists():
        fai.unlink()
    pysam.faidx(str(ref_out))

    bam_out = out / REAL_BAM_NAME
    header = {
        "HD": {"VN": "1.6", "SO": "coordinate"},
        "SQ": [{"SN": nm, "LN": lengths[nm]} for nm in names],
    }
    records.sort(key=lambda t: (t.tid, t.start, t.name))
    with pysam.AlignmentFile(str(bam_out), "wb", header=header) as w:
        for rec in records:
            a = pysam.AlignedSegment(w.header)
            a.query_name = rec.name
            a.query_sequence = rec.seq
            a.flag = rec.flag
            a.reference_id = rec.tid
            a.reference_start = rec.start
            a.mapping_quality = rec.mapq
            a.cigartuples = rec.cigar
            a.query_qualities = rec.quals
            if rec.mate_tid is not None:
                a.next_reference_id = rec.mate_tid
                a.next_reference_start = rec.mate_start
                a.template_length = rec.tlen
            if rec.mate_mapq is not None:
                a.set_tag("MQ", rec.mate_mapq, "i")
            w.write(a)
    bai = bam_out.with_name(REAL_BAM_NAME + ".bai")
    if bai.exists():
        bai.unlink()
    pysam.index(str(bam_out))
    return {"bam": bam_out, "reference": ref_out}, items, mapping


def real_label_record(item: Item) -> dict[str, Any]:
    lenient = [item.label["artifact_verdict"]]
    return {
        "item": item.item_id,
        "case_type": item.case_type,
        "instance": item.instance,
        "contig": item.contig,
        "positions": item.positions_1based(),
        **item.label,
        "lenient_verdicts": lenient,
    }


def build_real_run_set(
    out: str | Path,
    sites_path: str | Path,
    *,
    seed: int,
    conditions: list[str],
    models: list[str],
    repeats: int = 1,
    png: bool | None = None,
    fixture_dir: str | Path | None = None,
    bam_path: str = HG002_BAM,
    ref_path: str = GRCH38_FASTA,
    log: Any = print,
) -> Any:
    """Blinded fixtures from a frozen site manifest, then stimuli + labels + runs."""
    from .stimuli import build_items_run_set

    sites = sites_from_jsonl(sites_path)
    fx = Path(fixture_dir) if fixture_dir else Path(out) / "fixtures"
    paths, items, mapping = write_real_fixtures(
        sites, fx, seed, bam_path=bam_path, ref_path=ref_path
    )
    res = build_items_run_set(
        out,
        bam=str(paths["bam"]),
        ref=str(paths["reference"]),
        items=items,
        label_records=[real_label_record(it) for it in items],
        seed=seed,
        conditions=conditions,
        models=models,
        repeats=repeats,
        png=png,
        log=log,
        meta={
            "tier": "real",
            "sites": str(sites_path),
            "sites_sha256": hashlib.sha256(Path(sites_path).read_bytes()).hexdigest(),
        },
        sample_line=SAMPLE_LINE,
        sample_kind="sequenced",
    )
    with (Path(out) / "labels" / "mapping.jsonl").open("w") as f:
        for m in mapping:
            f.write(json.dumps(m, sort_keys=True) + "\n")
    return res


# ===========================================================================
# Round 4 (docs/eval/visual-design-r4.md): hard real strata
# ===========================================================================
#
# Candidates are manifest rows (dicts; ``CONTRACT.md``). One generator per stratum family
# yields rows whose ``features`` are computed from the 60x Illumina reads for the candidate
# allele (:func:`_allele_features`, cached in a JSONL scan cache). :func:`collect_pools_r4`
# runs the generators and labels each row with its cluster and partition;
# :func:`select_r4` applies the frozen predicates and draws a split; label checks against
# HiFi / the parents / gnomAD run lazily, only on rows a draw reaches, and are cached too.

R4_CONDITIONS = ("Traw-r4", "Traw-r4+Tcur", "V+Traw-r4")
NEAR_END_BP = 10
LINK_BP = 15
LOCAL_BP = 100
LOCAL_MIN_READS = 3
LOCAL_MIN_FRAC = 0.10
ISOLATION_R4 = 30
PSV_CLEAR_BP = 50
ADJ_MIN_GAP = 3  # "gap > 2 bp"
ADJ_MAX_GAP = 12
ADJ_CLEAR_BP = 3
MIN_TR_INDEL_BP = 3  # indel-adjacent strata: truth indel >= 3 bp
TR_MAX_INDEL = 4
CMRG_MAX_INDEL = 10
ANCHOR_CLEAR_BP = 10
REMOTE_RETRIES = 4
FETCH_BLOCK = 16_384
FETCH_BLOCKS_MAX = 512  # bounded: about 8 MB of reference
TR_NEG_MIN_READS = 4  # = strata_r4.TR_NEG_MIN_READS / TR_NEG_MIN_FRAC (knob-independent)
TR_NEG_MIN_FRAC = 0.10


# ---------------------------------------------------------------------------
# Scan cache
# ---------------------------------------------------------------------------


class ScanCache:
    """Append-only JSONL ``{"key": ..., "value": ...}`` cache of remote-BAM work.

    Keys name the BAM, site and allele, so a cache is reusable across draws and splits;
    ``path=None`` keeps it in memory only.
    """

    def __init__(self, path: str | Path | None = None):
        self.path = Path(path) if path else None
        self.data: dict[str, Any] = {}
        if self.path and self.path.exists():
            for line in self.path.read_text().splitlines():
                if line.strip():
                    rec = json.loads(line)
                    self.data[rec["key"]] = rec["value"]

    def get(self, key: str) -> Any:
        return self.data.get(key)

    def __contains__(self, key: str) -> bool:
        return key in self.data

    def put(self, key: str, value: Any) -> Any:
        self.data[key] = value
        if self.path:
            self.path.parent.mkdir(parents=True, exist_ok=True)
            with self.path.open("a") as f:
                f.write(json.dumps({"key": key, "value": value}, sort_keys=True) + "\n")
        return value

    def memo(self, key: str, fn: Any) -> Any:
        if key in self.data:
            return self.data[key]
        return self.put(key, fn())


# ---------------------------------------------------------------------------
# Read-level features
# ---------------------------------------------------------------------------


def _median(xs: list[float]) -> float | None:
    if not xs:
        return None
    s = sorted(xs)
    m = len(s) // 2
    return float(s[m]) if len(s) % 2 else (s[m - 1] + s[m]) / 2


def _mean(xs: list[float]) -> float | None:
    return sum(xs) / len(xs) if xs else None


def _frac(xs: list[bool]) -> float | None:
    return sum(xs) / len(xs) if xs else None


MATE_SEARCH_BP = 800


def _mate_mapqs(reads: Iterable[Any]) -> dict[tuple[str, bool], int]:
    """``(name, is_read1) -> MAPQ`` of primary alignments, to look a read's mate up by."""
    return {
        (r.query_name, bool(r.is_read1)): r.mapping_quality
        for r in reads
        if not (r.is_secondary or r.is_supplementary)
    }


def _mate_mapq(read: Any, mates: dict[tuple[str, bool], int]) -> int | None:
    if read.has_tag("MQ"):
        return int(read.get_tag("MQ"))
    if not read.is_paired or read.mate_is_unmapped:
        return None
    return mates.get((read.query_name, not read.is_read1))


@dataclass
class _ReadObs:
    reverse: bool
    bq: int | None
    mapq: int
    near_end: bool
    softclip: bool
    mismatches: float  # per 100 aligned bp, excluding the site
    proper: bool
    mate_mapq: int | None
    nonref_cols: set[int]


def _observe(
    read: Any, ref: str, ref0: int, skip: tuple[int, int]
) -> tuple[dict[int, str], set[int], float, int]:
    """Per-read reference comparison over its aligned bases.

    Returns ({ref_pos: base} for aligned bases, positions where the read differs from the
    reference (mismatch, deletion, insertion point, soft-clip boundary), mismatches per 100
    aligned bp outside ``skip``, aligned bp)."""
    bases: dict[int, str] = {}
    nonref: set[int] = set()
    mism = aligned = 0
    q = read.query_sequence or ""
    rpos = read.reference_start
    qpos = 0
    cig = read.cigartuples or []
    for k, (op, n) in enumerate(cig):
        if op in (0, 7, 8):
            for j in range(n):
                p = rpos + j
                b = q[qpos + j].upper() if qpos + j < len(q) else "N"
                bases[p] = b
                rb = ref[p - ref0] if 0 <= p - ref0 < len(ref) else "N"
                if b != rb and rb != "N" and b != "N":
                    nonref.add(p)
                    if not skip[0] <= p < skip[1]:
                        mism += 1
                if not skip[0] <= p < skip[1]:
                    aligned += 1
            rpos += n
            qpos += n
        elif op == 1:
            nonref.add(rpos)
            qpos += n
        elif op in (2, 3):
            for j in range(n):
                bases[rpos + j] = "-"
                nonref.add(rpos + j)
            rpos += n
        elif op == 4:
            nonref.add(rpos if k > 0 else rpos - 1)  # clip boundary column
            qpos += n
    return bases, nonref, (100.0 * mism / aligned if aligned else 0.0), aligned


def _allele_features(
    bam: Any,
    fasta: Any,
    chrom: str,
    pos0: int,
    *,
    allele: str | None = None,
    run: tuple[int, int] | None = None,
    delta: int | None = None,
    unit: str = "",
    truth_deltas: tuple[int, ...] = (0,),
) -> dict[str, Any]:
    """The contract's ``features`` dict for one candidate allele at ``pos0``.

    Two modes: a base allele (``allele``: SNV candidates; alt = reads with that base at
    ``pos0``) or a length allele (``run`` + ``delta``: reads spanning the run whose net indel
    length inside it equals ``delta``; ref = net length 0). Per class: strand counts, mean
    base quality at the site, MAPQ < 20 fraction, near-read-end fraction (< 10 bp from an
    aligned end), soft-clip fraction, mismatch burden (per 100 aligned bp, excluding the
    site, against the reference), proper-pair fraction, median mate MAPQ (``MQ``). Then the
    linkage maximum (the largest share of alt reads that also differ from the reference at
    one other column within +/-15 bp), the number of other non-reference columns within
    100 bp (>= 3 reads and >= 10%), the reference run length at ``pos0``, and for length
    alleles the +/-1-unit stutter fraction around the truth lengths ``truth_deltas``.
    """
    from .truth import read_run_delta, repeat_run

    lo, hi = pos0 - LOCAL_BP - 300, pos0 + LOCAL_BP + 300
    lo = max(0, lo)
    ref = fasta.fetch(chrom, lo, hi).upper()
    skip = run if run is not None else (pos0, pos0 + 1)
    ref_base = ref[pos0 - lo] if 0 <= pos0 - lo < len(ref) else "N"
    cols: Counter[int] = Counter()
    cover: Counter[int] = Counter()
    obs: dict[str, list[_ReadObs]] = {"alt": [], "ref": [], "other": []}
    hist: Counter[int] = Counter()
    # Mate MAPQ: the MQ tag when the aligner wrote one (Novoalign does not), else the mate
    # itself when it aligns within MATE_SEARCH_BP (most inserts are 300-700 bp).
    reads = list(bam.fetch(chrom, max(0, pos0 - MATE_SEARCH_BP), pos0 + MATE_SEARCH_BP))
    # The candidate's own columns (the site, or the whole run) are not "linked" mismatches.
    own_cols = set(range(skip[0] - 1, skip[1] + 1)) if run is not None else {pos0}
    mates = _mate_mapqs(reads)
    for r in reads:
        if not _usable(r) or r.reference_end is None:
            continue
        if r.reference_end <= pos0 - LOCAL_BP or r.reference_start > pos0 + LOCAL_BP:
            continue
        bases, nonref, mism, _ = _observe(r, ref, lo, skip)
        for p in bases:
            if abs(p - pos0) <= LOCAL_BP:
                cover[p] += 1
        for p in nonref:
            if abs(p - pos0) <= LOCAL_BP:
                cols[p] += 1
        if run is not None:
            d = read_run_delta(r, run[0], run[1])
            if d is None:
                continue
            hist[d] += 1
            cls = "alt" if d == delta else "ref" if d == 0 else "other"
            site_pos = run[0] - 1  # the anchor base: always aligned, unlike a deleted one
        else:
            if not r.reference_start <= pos0 < r.reference_end:
                continue
            b = bases.get(pos0)
            if b is None:
                continue
            cls = "alt" if b == (allele or "").upper() else "ref" if b == ref_base else "other"
            site_pos = pos0
        qpos = None
        for qp, rp in r.get_aligned_pairs(matches_only=True):
            if rp == site_pos:
                qpos = qp
                break
        quals = r.query_qualities
        bq = int(quals[qpos]) if qpos is not None and quals is not None else None
        if qpos is not None:
            near = min(qpos - r.query_alignment_start, r.query_alignment_end - 1 - qpos)
        else:
            near = min(site_pos - r.reference_start, r.reference_end - 1 - site_pos)
        obs[cls].append(
            _ReadObs(
                reverse=bool(r.is_reverse),
                bq=bq,
                mapq=r.mapping_quality,
                near_end=near < NEAR_END_BP,
                softclip=any(op == 4 for op, _ in r.cigartuples or ()),
                mismatches=mism,
                proper=bool(r.is_proper_pair),
                mate_mapq=_mate_mapq(r, mates),
                nonref_cols={
                    p for p in nonref if abs(p - site_pos) <= LINK_BP and p not in own_cols
                },
            )
        )
    alt, refo = obs["alt"], obs["ref"]
    depth = len(alt) + len(refo) + len(obs["other"])
    feats: dict[str, Any] = {"depth": depth, "alt_reads": len(alt)}
    if depth:
        feats["vaf"] = round(len(alt) / depth, 4)
    for name, group in (("alt", alt), ("ref", refo)):
        feats[f"{name}_fwd"] = sum(not o.reverse for o in group)
        feats[f"{name}_rev"] = sum(o.reverse for o in group)
        vals: dict[str, float | None] = {
            "mean_bq": _mean([o.bq for o in group if o.bq is not None]),
            "mapq_lt20_frac": _frac([o.mapq < 20 for o in group]),
            "near_end_frac": _frac([o.near_end for o in group]),
            "softclip_frac": _frac([o.softclip for o in group]),
            "mismatch_burden": _mean([o.mismatches for o in group]),
            "proper_pair_frac": _frac([o.proper for o in group]),
            "mate_mapq_median": _median([o.mate_mapq for o in group if o.mate_mapq is not None]),
        }
        for k, v in vals.items():
            if v is not None:
                feats[f"{name}_{k}"] = round(v, 4)
    if alt:
        link: Counter[int] = Counter()
        for o in alt:
            link.update(o.nonref_cols)
        feats["linked_mismatch_max_frac"] = round(max(link.values(), default=0) / len(alt), 4)
    skip_cols = set(range(skip[0] - 1, skip[1] + 1))
    feats["local_nonref_sites_100bp"] = sum(
        1
        for p, n in cols.items()
        if p not in skip_cols
        and n >= LOCAL_MIN_READS
        and cover.get(p, 0)
        and n / cover[p] >= LOCAL_MIN_FRAC
    )
    i = pos0 - lo
    if 0 <= i < len(ref):
        s, e, _u = repeat_run(ref, i)
        feats["run_length"] = e - s
    if run is not None:
        feats["run_length"] = run[1] - run[0]
        u = len(unit) or 1
        n = sum(hist.values())
        truth = set(truth_deltas)
        stut = {t + k for t in truth for k in (-u, u)} - truth
        feats["stutter_frac"] = round(sum(hist[d] for d in stut) / n, 4) if n else 0.0
    return feats


# ---------------------------------------------------------------------------
# Resources
# ---------------------------------------------------------------------------


@dataclass
class R4Resources:
    """Open handles and loaded truth for one selection run (see :func:`open_r4_resources`)."""

    bam: Any
    fasta: Any
    v421: Any
    v5q: Any
    cmrg: Any = None
    mosaic: Any = None
    segdups: Any = None
    lowmap: Any = None
    repeats: Any = None  # homopolymer + diTR stratification union (anchors' "not in repeat")
    variable_copy: Any = None  # SV spans: stratum 6 excludes sites in or near them
    tr_beds: dict[str, Any] = field(default_factory=dict)  # run sources for stratum 4
    clusters: Any = None
    cache: ScanCache = field(default_factory=ScanCache)
    hifi: Any = None
    parents: dict[str, Any] = field(default_factory=dict)
    parent_truth: dict[str, Any] | None = None
    bam_id: str = HG002_BAM
    log: Any = print
    # name ("bam", "fasta", "hifi", "parent:<sample>") -> re-open function, for retries
    openers: dict[str, Any] = field(default_factory=dict)
    _blocks: dict[tuple[str, int], str] = field(default_factory=dict, repr=False)

    def retry(self, fn: Any, tries: int = REMOTE_RETRIES) -> Any:
        """``fn()``, re-opening every remote handle and retrying when a remote read fails
        (:func:`truth.transient`: htslib surfaces a dropped connection as ``OSError:
        truncated file`` or a pileup ``ValueError: error during iteration``)."""
        import time

        for attempt in range(tries):
            try:
                return fn()
            except (OSError, ValueError) as exc:
                from .truth import transient

                if attempt == tries - 1 or not transient(exc):
                    raise
                self.log(f"remote read failed ({exc}); re-opening and retrying")
                time.sleep(2 * (attempt + 1))
                self.reopen()
        raise AssertionError("unreachable")  # pragma: no cover

    def reopen(self) -> None:
        for name, opener in self.openers.items():
            if name.startswith("parent:"):
                self.parents[name.split(":", 1)[1]] = opener()
            else:
                setattr(self, name, opener())

    def fetch(self, chrom: str, s: int, e: int) -> str:
        """Upper-case reference [s, e), served from cached 16 kb blocks: normalisation and
        run detection make many tiny reads, each a round trip on a remote FASTA."""
        s = max(0, s)
        if e <= s:
            return ""
        out = []
        for b in range(s // FETCH_BLOCK, (e - 1) // FETCH_BLOCK + 1):
            key = (chrom, b)
            blk = self._blocks.get(key)
            if blk is None:
                a0 = b * FETCH_BLOCK
                blk = str(
                    self.retry(lambda a0=a0: self.fasta.fetch(chrom, a0, a0 + FETCH_BLOCK))
                ).upper()
                if len(self._blocks) >= FETCH_BLOCKS_MAX:
                    self._blocks.pop(next(iter(self._blocks)))
                self._blocks[key] = blk
            out.append(blk)
        joined = "".join(out)
        off = s - (s // FETCH_BLOCK) * FETCH_BLOCK
        return joined[off : off + (e - s)]

    def features(self, chrom: str, pos0: int, **kw: Any) -> dict[str, Any]:
        run = kw.get("run")
        key = "|".join(
            str(x)
            for x in (
                "feat-v1",
                self.bam_id,
                chrom,
                pos0,
                kw.get("allele"),
                tuple(run) if run else None,
                kw.get("delta"),
                tuple(kw.get("truth_deltas", (0,))),
            )
        )
        return dict(
            self.cache.memo(
                key,
                lambda: self.retry(
                    lambda: _allele_features(self.bam, self.fasta, chrom, pos0, **kw)
                ),
            )
        )

    def in_segdup_or_lowmap(self, chrom: str, pos0: int) -> bool:
        return bool(
            (self.segdups is not None and self.segdups.contains(chrom, pos0))
            or (self.lowmap is not None and self.lowmap.contains(chrom, pos0))
        )


def _sv_spans(path: str | Path, chroms: Iterable[str] | None, pad: int = 50) -> Any:
    """Reference spans (+/- ``pad``) of every record in an SV VCF, as a BedSet."""
    from .truth import BedSet, TruthSet, cached

    # One whole-file pass over a local copy: cheaper than genome-wide remote region reads.
    ts = TruthSet.load("sv", cached(path), None, chroms=chroms)
    by: dict[str, Intervals] = {}
    for c, recs in ts.by_chrom.items():
        # GIAB SV VCFs are sequence-resolved: REF spans a deletion's deleted bases.
        by[c] = Intervals([(max(0, r.start0 - pad), r.end0 + pad) for r in recs])
    return BedSet(by)


def open_r4_resources(
    *,
    chroms: Iterable[str],
    truth_chroms: Iterable[str] | None = None,
    bam_path: str = HG002_BAM,
    bam_index: str | None = None,
    ref_path: str = GRCH38_FASTA,
    v421: tuple[str, str] | None = None,
    v5q: tuple[str, str] | None = None,
    cmrg: tuple[str, str] | None = None,
    mosaic: tuple[str, str] | None = None,
    strat: dict[str, str] | None = None,
    superdups: str | None = None,
    superdup_min_identity: float = 0.0,
    sv_vcfs: Iterable[str] | None = None,
    hifi: tuple[str, str | None] | None = None,
    parents: dict[str, tuple[str, str | None]] | None = None,
    parent_truth: dict[str, tuple[str, str]] | None = None,
    cache_path: str | Path | None = None,
    log: Any = print,
) -> R4Resources:
    """Open the BAMs/FASTA and load every truth set for ``chroms`` (scan-driven strata)
    plus ``truth_chroms`` (truth-list strata: mosaic, CMRG; default all autosomes).

    Every path defaults to the public URL in :mod:`.truth`; pass local copies to avoid
    re-downloads. ``hifi=None`` / ``parents=None`` skip those checks' handles (a draw then
    records the check as not run); pass ``(url, local_bai_or_None)`` to enable them.
    """
    import pysam

    from . import truth as T

    scan = list(chroms)
    tchroms = list(truth_chroms) if truth_chroms is not None else list(T.AUTOSOMES)
    allc = sorted(set(scan) | set(tchroms), key=T._chrom_order)
    strat = {**T.STRAT_URLS, **(strat or {})}
    log(f"loading truth sets for {len(allc)} chromosome(s) ...")
    res = R4Resources(
        bam=T.open_bam(bam_path, bam_index),
        fasta=pysam.FastaFile(ref_path),
        # v4.2.1 / v5.0q: scan chromosomes held in memory, others by region query;
        # CMRG and Mosaic are small and loaded whole.
        v421=T.TruthSet.load("v4.2.1", *(v421 or (T.V421_VCF, T.V421_BED)), chroms=scan),
        v5q=T.TruthSet.load("v5.0q", *(v5q or (T.V5Q_VCF, T.V5Q_BED)), chroms=scan),
        cmrg=T.TruthSet.load("cmrg_v1.00", *(cmrg or (T.CMRG_VCF, T.CMRG_BED))),
        mosaic=T.TruthSet.load("mosaic_v1.1", *(mosaic or (T.MOSAIC_VCF, T.MOSAIC_BED))),
        segdups=T.BedSet.load(strat["segdups"], allc),
        lowmap=T.BedSet.load(strat["lowmappability"], allc),
        cache=ScanCache(cache_path),
        bam_id=bam_path,
        log=log,
    )
    rep: dict[str, Intervals] = {}
    for key in ("homopolymer_ge12", "homopolymer_7to11", "diTR"):
        for c, iv in T.load_bed_multi(strat[key], allc).items():
            rep[c] = Intervals(rep.get(c, Intervals([])).pairs() + iv.pairs())
    res.repeats = T.BedSet(rep)
    res.tr_beds = {
        k: T.BedSet.load(strat[k], scan) for k in ("homopolymer_ge12", "homopolymer_7to11", "diTR")
    }
    sv = list(sv_vcfs) if sv_vcfs is not None else [T.CMRG_SV_VCF, T.V5Q_SV_VCF]
    by: dict[str, Intervals] = {}
    for p in sv:
        for c, iv in _sv_spans(p, allc).by_chrom.items():
            by[c] = Intervals(by.get(c, Intervals([])).pairs() + iv.pairs())
    res.variable_copy = T.BedSet(by)
    res.clusters = T.cluster_index(superdups, min_identity=superdup_min_identity)
    res.openers["bam"] = lambda: T.open_bam(bam_path, bam_index)
    res.openers["fasta"] = lambda: pysam.FastaFile(ref_path)
    if hifi is not None:
        h_url, h_bai = hifi
        res.hifi = T.open_bam(h_url, h_bai)
        res.openers["hifi"] = lambda: T.open_bam(h_url, h_bai)
    for name, (url, bai) in (parents or {}).items():
        res.parents[name] = T.open_bam(url, bai)
        res.openers[f"parent:{name}"] = lambda url=url, bai=bai: T.open_bam(url, bai)
    if parent_truth:
        res.parent_truth = {
            n: T.TruthSet.load(n, v, b, chroms=()) for n, (v, b) in parent_truth.items()
        }
    return res


# ---------------------------------------------------------------------------
# Row construction
# ---------------------------------------------------------------------------


def _new_row(
    stratum: str,
    chrom: str,
    site0: int,
    *,
    truth_source: str,
    label: dict[str, Any] | None = None,
    vcf: tuple[int | None, str, str] = (None, "", ""),
    context: dict[str, Any] | None = None,
) -> dict[str, Any]:
    from .strata_r4 import PAIRS

    pair, polarity = PAIRS[stratum]
    return {
        "chrom": chrom,
        "site0": int(site0),
        "stratum": stratum,
        "pair": pair,
        "polarity": polarity,
        "cluster": None,
        "partition": None,
        "knob": None,
        "truth_source": truth_source,
        "vcf_pos": vcf[0],
        "vcf_ref": vcf[1],
        "vcf_alt": vcf[2],
        "label": dict(label or {}),
        "features": {},
        "checks": {},
        "context": dict(context or {}),
        "depth": 0,
        "alt_reads": 0,
        "exclude": None,
    }


def _finish(row: dict[str, Any], feats: dict[str, Any]) -> dict[str, Any]:
    row["features"] = feats
    row["depth"] = int(feats.get("depth", 0))
    row["alt_reads"] = int(feats.get("alt_reads", 0))
    return row


def _excluded(row: dict[str, Any], reason: str) -> dict[str, Any]:
    row["exclude"] = reason
    return row


def _germline_label(rec: Any) -> tuple[dict[str, Any], int] | None:
    cls = classify_record(rec.ref, rec.alts[0], tuple(rec.gt)) if rec.biallelic else None
    if cls is None:
        return None
    return cls["label"], cls["offset"]


def _span_of(rec: Any, fetch: Any) -> tuple[int, int, str]:
    """(run start, run end, unit) for an indel record: the reference repeat run holding the
    first affected base, widened to the indel's own deleted span."""
    from .truth import repeat_run

    p0 = rec.start0
    w0 = max(0, p0 - 200)
    seq = fetch(rec.chrom, w0, p0 + 300)
    s, e, unit = repeat_run(seq, p0 + 1 - w0)
    s, e = s + w0, e + w0
    dl = len(rec.ref) - 1
    return min(s, p0 + 1), max(e, p0 + 1 + max(dl, 0)), unit


def _truth_sets(res: R4Resources) -> list[Any]:
    """Every loaded truth set: isolation and clearance count records of any of them
    (design 3.0, "no other truth record (any set, any FILTER)")."""
    return [ts for ts in (res.v421, res.v5q, res.cmrg, res.mosaic) if ts is not None]


def _same_variant(res: R4Resources, a: Any, b: Any) -> bool:
    """Do two truth records describe the same allele after left-normalisation (two sets'
    representations of one variant are not "other" records)?"""
    from .truth import normalized_key

    if not a.alts or not b.alts:
        return False
    want = {normalized_key(res.fetch, a, x) for x in a.alts}
    return any(normalized_key(res.fetch, b, x) in want for x in b.alts)


def _truth_near(
    res: R4Resources,
    c: str,
    s0: int,
    e0: int,
    bp: int,
    *,
    sets: Iterable[Any] | None = None,
    own: Any = None,
    keep: Callable[[Any], bool] | None = None,
) -> Any:
    """The first truth set (of ``sets``, default every loaded one) holding a record of any
    FILTER within ``bp`` of [s0, e0) that is not ``own`` (or the same variant as it) and
    passes ``keep``; None when the window is clear."""
    for ts in sets if sets is not None else _truth_sets(res):
        for o in ts.records_near(c, s0, e0, bp):
            if own is not None and (o is own or o == own or _same_variant(res, own, o)):
                continue
            if keep is not None and not keep(o):
                continue
            return ts
    return None


def _allele_in_unit(allele: str, unit: str) -> bool:
    if not unit or len(allele) % len(unit):
        return False
    return allele.upper() in (unit * (len(allele) // len(unit) + 2)).upper()


# ---------------------------------------------------------------------------
# Generators
# ---------------------------------------------------------------------------


def iter_mosaic_sites(res: R4Resources, chroms: Iterable[str]) -> Iterable[dict[str, Any]]:
    """Stratum 1: Mosaic SNV v1.1 records (genome-wide by default). A site with any v5.0q
    record (any FILTER) within 30 bp is dropped; v4.2.1 is not consulted (it cannot see
    < 30% VAF)."""
    for c in chroms:
        for rec in res.mosaic.records(c):
            if rec.kind != "SNV":
                continue
            alt = rec.alts[0].upper()
            label = {
                "variant_present": "yes",
                "variant_type": "SNV",
                "alt_allele": alt,
                "zygosity": "low_vaf",
                "artifact_verdict": "real",
                "artifact_type": "none",
                "phase": "na",
            }
            row = _new_row(
                "mosaic_snv_pos",
                c,
                rec.start0,
                truth_source="mosaic_v1.1",
                label=label,
                vcf=(rec.pos, rec.ref, alt),
                context={
                    "candidate_allele": alt,
                    "segdup_or_lowmap": res.in_segdup_or_lowmap(c, rec.start0),
                },
            )
            if res.v5q.records_near(c, rec.start0, rec.end0, ISOLATION_R4):
                yield _excluded(row, "v5q_within_30bp")
                continue
            yield _finish(row, res.features(c, rec.start0, allele=alt))


def _caller_records(sources: dict[str, str], chroms: Iterable[str]) -> Iterable[tuple[str, Any]]:
    from .truth import TruthSet

    for name, path in sources.items():
        ts = TruthSet.load(name, path, None, chroms=list(chroms))
        for c in ts.chroms():
            for rec in ts.records(c):
                yield name, rec


def _database_records(path: str | Path, chroms: Iterable[str]) -> Iterable[tuple[str, Any]]:
    """``putative`` / ``candidate`` SNV rows of the mosaic release's annotated database."""
    from .truth import TruthRecord, cached

    want = set(chroms)
    with cached(path).open() as fh:
        header = fh.readline().rstrip("\n").split("\t")
        col = {h: i for i, h in enumerate(header)}
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if f[col["CHROM"]] not in want:
                continue
            if f[col["var_cat"]] not in ("putative", "candidate"):
                continue
            if f[col["var_type"]].lower() != "snv":
                continue
            yield (
                "mosaic_database",
                TruthRecord(
                    f[col["CHROM"]],
                    int(f[col["POS"]]),
                    f[col["REF"]],
                    (f[col["ALT"]],),
                    (),
                    "PASS",
                    "mosaic_database",
                ),
            )


def iter_mosaic_caller_fps(
    res: R4Resources,
    chroms: Iterable[str],
    *,
    sources: dict[str, str] | None = None,
    database: str | Path | None = None,
    seed: int = 0,
    limit: int | None = None,
) -> Iterable[dict[str, Any]]:
    """Stratum 2: caller-flagged SNVs the mosaic benchmark rejected.

    Sources (design 3.2): the release's DRAGEN mosaic-only, DeepSomatic Illumina-300x Q0
    and Strelka2 Ill300X calls (``truth.MOSAIC_CALLER_VCFS``) and, when ``database`` is
    given, the putative/candidate rows of ``HG002_GRCh38_mosaic_database.tsv`` -- each
    minus the mosaic VCF. A candidate must lie inside the mosaic BED (which asserts < 2%
    VAF at any unlisted SNV) with no v4.2.1, v5.0q or mosaic record within 30 bp.
    Candidates are read from the BAM in seeded order, at most ``limit`` of them (chr20
    alone carries about 12,000 caller records, mostly Strelka2's).
    """
    from . import truth as T

    chroms = list(chroms)
    seen: dict[tuple[str, int, str], dict[str, Any]] = {}
    streams: list[Iterable[tuple[str, Any]]] = [
        _caller_records(sources if sources is not None else T.MOSAIC_CALLER_VCFS, chroms)
    ]
    if database is not None:
        streams.append(_database_records(database, chroms))
    order: list[tuple[str, int, str]] = []
    for stream in streams:
        for name, rec in stream:
            if rec.kind != "SNV":
                continue
            key = (rec.chrom, rec.start0, rec.alts[0].upper())
            if key in seen:
                seen[key]["context"]["fp_sources"].append(name)
                continue
            row = _new_row(
                "mosaic_caller_fp_neg",
                rec.chrom,
                rec.start0,
                truth_source="mosaic_v1.1",
                label=None,
                vcf=(rec.pos, rec.ref, key[2]),
                context={
                    "candidate_allele": key[2],
                    "fp_sources": [name],
                    "segdup_or_lowmap": res.in_segdup_or_lowmap(rec.chrom, rec.start0),
                },
            )
            seen[key] = row
            order.append(key)
    order.sort()
    _derive_rng(seed, "mosaic-fp").shuffle(order)
    scanned = 0
    for key in order:
        if limit is not None and scanned >= limit:
            return
        row = seen[key]
        c, p0, alt = key
        if not res.mosaic.covers(c, p0):
            yield _excluded(row, "outside_mosaic_bed")
            continue
        order = [res.mosaic, *(t for t in _truth_sets(res) if t is not res.mosaic)]
        near = _truth_near(res, c, p0, p0 + 1, ISOLATION_R4, sets=order)
        if near is not None:
            yield _excluded(row, f"{near.name}_within_30bp")
        else:
            scanned += 1
            yield _finish(row, res.features(c, p0, allele=alt))


def iter_tr_sites(
    res: R4Resources,
    chroms: Iterable[str],
    *,
    polarity: str,
    seed: int,
    limit: int | None = None,
    hp_min: int = 10,
) -> Iterable[dict[str, Any]]:
    """Strata 3 (``polarity="pos"``) and 4 (``"neg"``): repeat-run indels / stutter.

    pos: v5.0q PASS in-BED biallelic indels of 1-4 bp whose allele is whole units of a
    homopolymer >= ``hp_min`` bp or di/tri STR >= 20 bp; the run +/- 30 bp holds no other
    truth record (v5.0q any FILTER, v4.2.1); v4.2.1 must agree where it covers.
    neg: stratification runs inside the v5.0q BED with no truth record within the run +/-
    30 bp; the candidate is the commoner +/-1-unit length on spanning reads. The POI is the
    left-normalised anchor in both: a deletion's first deleted base (the run's first
    base), an insertion's anchor (the base before the run). ``limit`` caps how many loci
    are read from the BAM, in seeded order.
    """
    from . import truth as T

    rng = _derive_rng(seed, "tr", polarity)
    scanned = 0
    if polarity == "pos":
        recs = [
            r
            for c in chroms
            for r in res.v5q.positives(c)
            if r.kind in ("ins", "del") and r.zygosity in ("het", "hom_alt")
        ]
        recs.sort(key=lambda r: (r.chrom, r.pos))
        rng.shuffle(recs)
        for rec in recs:
            if limit is not None and scanned >= limit:
                return
            got = _germline_label(rec)
            if got is None:
                continue
            label, offset = got
            allele = label["alt_allele"]
            if not 1 <= len(allele) <= TR_MAX_INDEL:
                continue
            s, e, unit = _span_of(rec, res.fetch)
            if not _allele_in_unit(allele, unit) or not T.is_tr_run(e - s, unit, hp_min):
                continue
            delta = len(allele) if label["variant_type"] == "ins" else -len(allele)
            v4_cov = res.v421.covers(rec.chrom, rec.start0, rec.end0)
            row = _new_row(
                "tr_indel_pos",
                rec.chrom,
                rec.start0 + offset,
                truth_source="v5.0q+v4.2.1" if v4_cov else "v5.0q",
                label=label,
                vcf=(rec.pos, rec.ref, rec.alts[0]),
                context={"run": [s, e], "unit": unit, "delta": delta},
            )
            # The run +/- 30 bp holds no truth record of any set other than this variant
            # (another set's representation of the same allele is not "other").
            if _truth_near(res, rec.chrom, s, e, ISOLATION_R4, own=rec) is not None:
                yield _excluded(row, "truth_in_run_30bp")
                continue
            if T.concordant(rec, res.v421, res.fetch) is False:
                yield _excluded(row, "discordant_v4.2.1")
                continue
            scanned += 1
            truth = (0, delta) if label["zygosity"] == "het" else (delta,)
            yield _finish(
                row,
                res.features(
                    rec.chrom,
                    row["site0"],
                    run=(s, e),
                    delta=delta,
                    unit=unit,
                    truth_deltas=truth,
                ),
            )
        return

    beds = res.tr_beds
    loci: list[tuple[str, int, int]] = []
    for c in chroms:
        for key in ("homopolymer_ge12", "homopolymer_7to11", "diTR"):
            if key in beds:
                loci.extend((c, s, e) for s, e in beds[key].get(c).pairs())
    loci = sorted(set(loci))
    rng.shuffle(loci)
    for c, a, b in loci:
        if limit is not None and scanned >= limit:
            return
        seq = res.fetch(c, a, b)
        best = None
        i = 0
        while i < len(seq):
            s, e, unit = T.repeat_run(seq, i)
            if T.is_tr_run(e - s, unit, hp_min) and (best is None or e - s > best[1] - best[0]):
                best = (s + a, e + a, unit)
            i = max(e, i + 1)
        if best is None:
            continue
        s, e, unit = best
        if not res.v5q.covers(c, s, e):
            continue
        u = len(unit)
        row = _new_row(
            "tr_stutter_neg",
            c,
            s,
            truth_source="v5.0q+v4.2.1" if res.v421.covers(c, s, e) else "v5.0q",
            label=None,
            context={"run": [s, e], "unit": unit},
        )
        if _truth_near(res, c, s, e, ISOLATION_R4) is not None:
            yield _excluded(row, "truth_in_run_30bp")
            continue
        if T.reference_concordant(res.v421, c, s, e, ISOLATION_R4) is False:
            yield _excluded(row, "discordant_v4.2.1")
            continue
        scanned += 1
        hist = res.cache.memo(
            f"hist|{res.bam_id}|{c}|{s}|{e}",
            lambda c=c, s=s, e=e: res.retry(
                lambda: {str(k): v for k, v in T.spanning_hist(res.bam, c, s, e).items()}
            ),
        )
        h = {int(k): v for k, v in hist.items()}
        delta = -u if h.get(-u, 0) >= h.get(u, 0) else u
        n_span = sum(h.values())
        if h.get(delta, 0) < TR_NEG_MIN_READS or h.get(delta, 0) < TR_NEG_MIN_FRAC * n_span:
            # Below every piloted level of the stratum-4 signal rule: skip the full read pass.
            yield _excluded(row, "stutter_signal_weak")
            continue
        row["site0"] = s if delta < 0 else s - 1
        row["context"].update(
            delta=delta,
            candidate_allele=("del:" if delta < 0 else "ins:") + unit,
        )
        yield _finish(
            row,
            res.features(c, row["site0"], run=(s, e), delta=delta, unit=unit, truth_deltas=(0,)),
        )


def iter_cmrg_sites(res: R4Resources, chroms: Iterable[str]) -> Iterable[dict[str, Any]]:
    """Stratum 5: CMRG v1.00 PASS in-BED SNVs and indels <= 10 bp inside segdups or
    low-mappability, outside the v4.2.1 BED, with no other CMRG record within 30 bp; v5.0q
    must agree where it covers (germline-vs-germline, left-normalised)."""
    from . import truth as T

    for c in chroms:
        for rec in res.cmrg.positives(c):
            if not res.in_segdup_or_lowmap(c, rec.start0):
                continue
            got = _germline_label(rec)
            if got is None:
                continue
            label, offset = got
            if res.v421.bed is not None and res.v421.bed.overlaps(c, rec.start0, rec.end0):
                continue
            src = "cmrg_v1.00+v5.0q" if res.v5q.covers(c, rec.start0, rec.end0) else "cmrg_v1.00"
            ctx: dict[str, Any] = {"segdup_or_lowmap": True}
            row = _new_row(
                "segdup_cmrg_pos",
                c,
                rec.start0 + offset,
                truth_source=src,
                label=label,
                vcf=(rec.pos, rec.ref, rec.alts[0]),
                context=ctx,
            )
            if _truth_near(res, c, rec.start0, rec.end0, ISOLATION_R4, sets=[res.cmrg], own=rec):
                yield _excluded(row, "cmrg_within_30bp")
                continue
            others = [ts for ts in _truth_sets(res) if ts is not res.cmrg]
            if _truth_near(res, c, rec.start0, rec.end0, ISOLATION_R4, sets=others, own=rec):
                yield _excluded(row, "truth_within_30bp")
                continue
            if T.concordant(rec, res.v5q, res.fetch) is False:
                yield _excluded(row, "discordant_v5.0q")
                continue
            if label["variant_type"] == "SNV":
                feats = res.features(c, rec.start0, allele=label["alt_allele"])
            else:
                s, e, unit = _span_of(rec, res.fetch)
                n = len(label["alt_allele"])
                delta = n if label["variant_type"] == "ins" else -n
                ctx.update(run=[s, e], unit=unit, delta=delta)
                truth = (0, delta) if label["zygosity"] == "het" else (delta,)
                feats = res.features(
                    c, row["site0"], run=(s, e), delta=delta, unit=unit, truth_deltas=truth
                )
            yield _finish(row, feats)


def _intersect(a: list[tuple[int, int]], b: list[tuple[int, int]]) -> list[tuple[int, int]]:
    out: list[tuple[int, int]] = []
    i = j = 0
    while i < len(a) and j < len(b):
        s, e = max(a[i][0], b[j][0]), min(a[i][1], b[j][1])
        if s < e:
            out.append((s, e))
        if a[i][1] < b[j][1]:
            i += 1
        else:
            j += 1
    return out


def _subtract(a: list[tuple[int, int]], b: list[tuple[int, int]]) -> list[tuple[int, int]]:
    out: list[tuple[int, int]] = []
    j = 0
    for s, e in a:
        cur = s
        while j < len(b) and b[j][1] <= cur:
            j += 1
        k = j
        while k < len(b) and b[k][0] < e:
            if b[k][0] > cur:
                out.append((cur, b[k][0]))
            cur = max(cur, b[k][1])
            k += 1
        if cur < e:
            out.append((cur, e))
    return out


def psv_territory(res: R4Resources, chrom: str) -> list[tuple[int, int]]:
    """CMRG BED ∩ v5.0q BED ∩ (segdups ∪ low-mappability), minus the v4.2.1 BED."""

    def pairs(bs: Any) -> list[tuple[int, int]]:
        return bs.get(chrom).pairs() if bs is not None else []

    mapp = Intervals(pairs(res.segdups) + pairs(res.lowmap)).pairs()
    t = _intersect(_intersect(pairs(res.cmrg.bed), pairs(res.v5q.bed)), mapp)
    return _subtract(t, pairs(res.v421.bed))


def _psv_exclusion(res: R4Resources, c: str, p0: int) -> str | None:
    if not (
        res.cmrg.covers(c, p0 - PSV_CLEAR_BP, p0 + PSV_CLEAR_BP + 1)
        and res.v5q.covers(c, p0 - PSV_CLEAR_BP, p0 + PSV_CLEAR_BP + 1)
    ):
        return "not_50bp_inside_beds"
    if res.v421.bed is not None and res.v421.bed.contains(c, p0):
        return "inside_v4.2.1_bed"
    if not res.in_segdup_or_lowmap(c, p0):
        return "not_segdup_or_lowmap"
    order = [res.cmrg, res.v5q, *(t for t in _truth_sets(res) if t not in (res.cmrg, res.v5q))]
    near = _truth_near(res, c, p0, p0 + 1, PSV_CLEAR_BP, sets=order)
    if near is not None:
        return f"{near.name}_within_50bp"
    if res.variable_copy is not None and res.variable_copy.contains(c, p0):
        return "variable_copy"
    return None


def iter_psv_candidates(
    res: R4Resources,
    chroms: Iterable[str],
    *,
    seed: int,
    dragen_vcf: str | None = None,
    scan_windows: int | None = 800,
    window_bp: int = 2000,
    min_frac: float = 0.07,
) -> Iterable[dict[str, Any]]:
    """Stratum 6: paralogous-sequence-variant negatives in the PSV territory
    (:func:`psv_territory`).

    Primary source: DRAGEN 300x germline calls (``truth.DRAGEN_300X_VCF``) minus CMRG truth
    (what a user would actually review). Supplement: a pileup scan of up to
    ``scan_windows`` seeded 2 kb windows of the territory (one non-reference base at >= 3
    reads and >= ``min_frac``, the lowest piloted floor). Both require >= 50 bp inside the
    CMRG and v5.0q BEDs, no CMRG or v5.0q record within 50 bp, and no variable-copy span.
    """
    from . import truth as T

    chroms = list(chroms)
    seen: set[tuple[str, int]] = set()

    def make(
        c: str, p0: int, alt: str, src: str, vcf: tuple[int | None, str, str]
    ) -> dict[str, Any]:
        row = _new_row(
            "segdup_psv_neg",
            c,
            p0,
            truth_source="cmrg_v1.00+v5.0q",
            vcf=vcf,
            context={"candidate_allele": alt, "psv_source": src, "segdup_or_lowmap": True},
        )
        why = _psv_exclusion(res, c, p0)
        if why:
            return _excluded(row, why)
        return _finish(row, res.features(c, p0, allele=alt))

    terr = {c: psv_territory(res, c) for c in chroms}
    dragen = dragen_vcf or T.DRAGEN_300X_VCF
    holder = [T.open_vcf(dragen)]
    contigs = set(holder[0].header.contigs)
    for c in chroms:
        if c not in contigs:
            continue
        for s, e in terr[c]:
            for r in T.fetch_vcf(holder, dragen, c, s, e):
                if len(r.alts or ()) != 1 or len(r.ref) != 1 or len(r.alts[0]) != 1:
                    continue
                p0 = r.pos - 1
                if (c, p0) in seen or res.cmrg.find(c, r.pos, r.ref, r.alts[0]):
                    continue
                seen.add((c, p0))
                yield make(c, p0, r.alts[0].upper(), "dragen_fp", (r.pos, r.ref, r.alts[0]))
    if not scan_windows:
        return
    windows = [(c, a) for c in chroms for s, e in terr[c] for a in range(s, e, window_bp)]
    _derive_rng(seed, "psv-scan").shuffle(windows)
    for n, (c, a) in enumerate(windows[:scan_windows]):
        ends = [e for s, e in terr[c] if s <= a < e]
        b = min(a + window_bp, ends[0] if ends else a + window_bp)
        ref = res.fetch(c, a, b)
        counts = res.cache.memo(
            f"pile|{res.bam_id}|{c}|{a}|{b}",
            lambda c=c, a=a, b=b: res.retry(
                lambda: {str(p): dict(v) for p, v in _allele_counts(res.bam, c, a, b).items()}
            ),
        )
        for ps, cnt in counts.items():
            p0 = int(ps)
            if (c, p0) in seen or not a <= p0 < b:
                continue
            rb = ref[p0 - a]
            depth = sum(cnt.values())
            nonref = {k: v for k, v in cnt.items() if k != rb and len(k) == 1 and k in "ACGT"}
            if not nonref or depth < 10:
                continue
            top, k = max(nonref.items(), key=lambda kv: (kv[1], kv[0]))
            if k < 3 or k / depth < min_frac:
                continue
            seen.add((c, p0))
            yield make(c, p0, top, "pileup_scan", (None, "", ""))
        if (n + 1) % 100 == 0:
            res.log(f"psv scan: {n + 1} windows")


def _gap(p: int, s: int, e: int) -> int:
    """Bases from ``p`` to the span [s, e): 1 = adjacent; 0 inside."""
    if p < s:
        return s - p
    if p >= e:
        return p - e + 1
    return 0


def anchoring_safe(fetch: Any, chrom: str, poi0: int, spans: Iterable[tuple[int, int]]) -> bool:
    """Can the prompt's anchoring sentence ("an indel may be anchored anywhere within the
    repeat run that contains the position of interest, or within 2 bp of it") map the POI
    onto a truth record? Safe when no span touches the POI's run (any unit 1-3, >= 2
    copies) widened by 2 bp."""
    from .truth import repeat_run

    w0 = max(0, poi0 - 200)
    seq = fetch(chrom, w0, poi0 + 200)
    s, e, _ = repeat_run(seq, poi0 - w0)
    lo, hi = min(s + w0, poi0) - 2, max(e + w0, poi0 + 1) + 2
    return not any(a < hi and b > lo for a, b in spans)


def _adjacent_problem(
    res: R4Resources,
    c: str,
    p: int,
    indel: tuple[int, int],
    indel_run: tuple[int, int],
    spans: list[tuple[int, int]],
    *,
    own_snv: bool = False,
) -> str | None:
    """Why POI ``p`` cannot sit next to the truth indel spanning ``indel`` (None = it can).

    Records overlapping the indel's own span (either set's representation of it) are not
    "other" truth records; with ``own_snv`` (stratum 8) neither is the SNV at the POI.
    """
    s, e = indel
    if not ADJ_MIN_GAP <= _gap(p, s, e) <= ADJ_MAX_GAP:
        return "gap"
    if indel_run[0] <= p < indel_run[1]:
        return "in_indel_run"
    if not (res.v421.covers(c, p) and res.v5q.covers(c, p)):
        return "outside_beds"

    def other(o: Any) -> bool:
        return not o.overlaps(s, e) and not (own_snv and o.start0 == p and o.end0 == p + 1)

    if _truth_near(res, c, p, p + 1, ADJ_CLEAR_BP, keep=other) is not None:
        return "truth_within_3bp"
    sp = [x for x in spans if not (own_snv and x[0] == p and x[1] == p + 1)]
    if not anchoring_safe(res.fetch, c, p, sp):
        return "anchoring_unsafe"
    return None


def iter_indel_adjacent(
    res: R4Resources,
    chroms: Iterable[str],
    *,
    polarity: str,
    seed: int,
    limit: int | None = None,
    min_frac: float = 0.10,
) -> Iterable[dict[str, Any]]:
    """Strata 7 (``polarity="neg"``, misalignment next to a truth indel) and 8 (``"pos"``,
    a truth SNV next to one).

    Anchor indels: >= 3 bp, PASS in both v4.2.1 and v5.0q BEDs, concordant. Distances are
    from the indel's full reference span: 2 < gap <= 12 bp; the POI lies outside the
    indel's repeat run and passes :func:`anchoring_safe`; no other truth record (any set,
    any FILTER) within 3 bp; inside both BEDs. pos: a truth SNV there, in both sets and
    agreeing. neg: one non-reference base on >= 4 reads and >= ``min_frac``.
    """
    from . import truth as T

    rng = _derive_rng(seed, "adj", polarity)
    indels = [
        r
        for c in chroms
        for r in res.v421.positives(c)
        if r.kind in ("ins", "del")
        and abs(len(r.ref) - len(r.alts[0])) >= MIN_TR_INDEL_BP
        and r.zygosity in ("het", "hom_alt")
    ]
    indels.sort(key=lambda r: (r.chrom, r.pos))
    rng.shuffle(indels)
    stratum = "indel_adjacent_snv_pos" if polarity == "pos" else "indel_adjacent_misalign_neg"
    scanned = 0
    seen: set[tuple[str, int]] = set()
    for ind in indels:
        if limit is not None and scanned >= limit:
            return
        c = ind.chrom
        if (
            not res.v5q.covers(c, ind.start0, ind.end0)
            or T.concordant(ind, res.v5q, res.fetch) is not True
        ):
            continue
        s, e = ind.start0, ind.end0
        rs, re_, _unit = _span_of(ind, res.fetch)
        spans = [
            (o.start0, o.end0) for ts in (res.v421, res.v5q) for o in ts.records_near(c, s, e, 40)
        ]
        ctx = {"indel": [s, e], "indel_run": [rs, re_]}

        scanned += 1
        if polarity == "pos":
            for snv in res.v421.records_near(c, s, e, ADJ_MAX_GAP):
                if snv.kind != "SNV" or snv.zygosity not in ("het", "hom_alt") or not snv.is_pass:
                    continue
                p = snv.start0
                if (c, p) in seen:
                    continue
                seen.add((c, p))
                got = _germline_label(snv)
                if got is None:
                    continue
                row = _new_row(
                    stratum,
                    c,
                    p,
                    truth_source="v4.2.1+v5.0q",
                    label=got[0],
                    vcf=(snv.pos, snv.ref, snv.alts[0]),
                    context={**ctx, "gap": _gap(p, s, e)},
                )
                why = _adjacent_problem(res, c, p, (s, e), (rs, re_), spans, own_snv=True)
                if why is None and T.concordant(snv, res.v5q, res.fetch) is not True:
                    why = "discordant_v5.0q"
                if why:
                    yield _excluded(row, why)
                    continue
                yield _finish(row, res.features(c, p, allele=snv.alts[0].upper()))
            continue
        a, b = s - ADJ_MAX_GAP - 1, e + ADJ_MAX_GAP + 1
        ref = res.fetch(c, a, b)
        counts = res.cache.memo(
            f"pile|{res.bam_id}|{c}|{a}|{b}",
            lambda c=c, a=a, b=b: res.retry(
                lambda: {str(p): dict(v) for p, v in _allele_counts(res.bam, c, a, b).items()}
            ),
        )
        for ps, cnt in counts.items():
            p = int(ps)
            if (c, p) in seen or not ADJ_MIN_GAP <= _gap(p, s, e) <= ADJ_MAX_GAP:
                continue
            rb = ref[p - a]
            depth = sum(cnt.values())
            nonref = {k: v for k, v in cnt.items() if k != rb and len(k) == 1 and k in "ACGT"}
            if not nonref or not depth:
                continue
            top, k = max(nonref.items(), key=lambda kv: (kv[1], kv[0]))
            if k < 4 or k / depth < min_frac:
                continue
            seen.add((c, p))
            row = _new_row(
                stratum,
                c,
                p,
                truth_source="v4.2.1+v5.0q",
                context={**ctx, "gap": _gap(p, s, e), "candidate_allele": top},
            )
            why = _adjacent_problem(res, c, p, (s, e), (rs, re_), spans)
            if why:
                yield _excluded(row, why)
                continue
            yield _finish(row, res.features(c, p, allele=top))


def iter_anchor_sites(
    res: R4Resources, chroms: Iterable[str], *, seed: int, limit: int = 200
) -> Iterable[dict[str, Any]]:
    """Stratum 9 anchors: clean het / hom SNVs (v4.2.1 and v5.0q agree, isolated 30 bp,
    outside repeats and segdups) and clean negatives (both BEDs, no record within 10 bp,
    no non-reference read, depth >= 20), up to ``limit`` of each kind in seeded order."""
    from . import truth as T

    rng = _derive_rng(seed, "anchor")
    snvs = [
        r
        for c in chroms
        for r in res.v421.positives(c)
        if r.kind == "SNV" and r.zygosity in ("het", "hom_alt")
    ]
    snvs.sort(key=lambda r: (r.chrom, r.pos))
    rng.shuffle(snvs)
    n_snv = {"het": 0, "hom_alt": 0}
    for rec in snvs:
        z = rec.zygosity
        if n_snv[z] >= limit:
            if min(n_snv.values()) >= limit:
                break
            continue
        c, p = rec.chrom, rec.start0
        if res.repeats.contains(c, p) or res.in_segdup_or_lowmap(c, p):
            continue
        if _truth_near(res, c, p, p + 1, ISOLATION_R4, own=rec) is not None:
            continue
        if T.concordant(rec, res.v5q, res.fetch) is not True:
            continue
        got = _germline_label(rec)
        if got is None:
            continue
        n_snv[z] += 1
        stratum = "anchor_het_snv" if z == "het" else "anchor_hom_snv"
        row = _new_row(
            stratum,
            c,
            p,
            truth_source="v4.2.1+v5.0q",
            label=got[0],
            vcf=(rec.pos, rec.ref, rec.alts[0]),
        )
        yield _finish(row, res.features(c, p, allele=rec.alts[0].upper()))
    made = tries = 0
    both = {c: _intersect(res.v421.bed.get(c).pairs(), res.v5q.bed.get(c).pairs()) for c in chroms}
    flat = [(c, s, e) for c, ivs in both.items() for s, e in ivs if e - s > 2 * ANCHOR_CLEAR_BP]
    while made < limit and flat and tries < 20 * limit:
        tries += 1
        c, s, e = rng.choice(flat)
        p = rng.randint(s + ANCHOR_CLEAR_BP, e - ANCHOR_CLEAR_BP - 1)
        if res.repeats.contains(c, p) or res.in_segdup_or_lowmap(c, p):
            continue
        if _truth_near(res, c, p, p + 1, ANCHOR_CLEAR_BP) is not None:
            continue
        rb = res.fetch(c, p, p + 1)
        cnt = res.retry(lambda c=c, p=p: T.base_counts(res.bam, c, p))
        depth = sum(cnt.values())
        alt = depth - cnt.get(rb, 0)
        if depth < 20 or alt:
            continue
        made += 1
        row = _new_row("anchor_clean_neg", c, p, truth_source="v4.2.1+v5.0q")
        yield _finish(row, {"depth": depth, "alt_reads": 0, "vaf": 0.0})


# ---------------------------------------------------------------------------
# Pools, lazy checks, the draw
# ---------------------------------------------------------------------------

REQUIRED_TRIO = frozenset({"tr_stutter_neg", "segdup_psv_neg"})


def collect_pools_r4(
    res: R4Resources,
    *,
    chroms: Iterable[str],
    truth_chroms: Iterable[str],
    seed: int,
    pairs: Iterable[str] = ("mosaic", "tr", "segdup", "indel_adj", "anchor"),
    limits: dict[str, int | None] | None = None,
    fp_sources: dict[str, str] | None = None,
    mosaic_database: str | Path | None = None,
    dragen_vcf: str | None = None,
) -> list[dict[str, Any]]:
    """Run the generators for ``pairs`` and return every candidate row (excluded ones keep
    their ``exclude`` reason), each with its ``cluster`` and ``partition``.

    The pool seed is fixed per pool (not per draw), so the pilot and main draws see the same
    candidates. Mosaic positives, CMRG positives and PSV negatives are truth-list driven and
    use ``truth_chroms``; everything else scans ``chroms``.
    """
    from .truth import partition

    chroms, truth_chroms = list(chroms), list(truth_chroms)
    lim = {
        "tr": 1500,
        "indel_adj": 1500,
        "anchor": 60,
        "psv_windows": 800,
        "mosaic_fp": 3000,
        **(limits or {}),
    }
    gens: list[tuple[str, Iterable[dict[str, Any]]]] = []
    want = set(pairs)
    if "mosaic" in want:
        gens.append(("mosaic_snv_pos", iter_mosaic_sites(res, truth_chroms)))
        gens.append(
            (
                "mosaic_caller_fp_neg",
                iter_mosaic_caller_fps(
                    res,
                    chroms,
                    sources=fp_sources,
                    database=mosaic_database,
                    seed=seed,
                    limit=lim["mosaic_fp"],
                ),
            )
        )
    if "tr" in want:
        gens.append(
            ("tr_indel_pos", iter_tr_sites(res, chroms, polarity="pos", seed=seed, limit=lim["tr"]))
        )
        gens.append(
            (
                "tr_stutter_neg",
                iter_tr_sites(res, chroms, polarity="neg", seed=seed, limit=lim["tr"]),
            )
        )
    if "segdup" in want:
        gens.append(("segdup_cmrg_pos", iter_cmrg_sites(res, truth_chroms)))
        gens.append(
            (
                "segdup_psv_neg",
                iter_psv_candidates(
                    res,
                    truth_chroms,
                    seed=seed,
                    dragen_vcf=dragen_vcf,
                    scan_windows=lim["psv_windows"],
                ),
            )
        )
    if "indel_adj" in want:
        for pol in ("pos", "neg"):
            gens.append(
                (
                    f"indel_adjacent_{pol}",
                    iter_indel_adjacent(
                        res, chroms, polarity=pol, seed=seed, limit=lim["indel_adj"]
                    ),
                )
            )
    if "anchor" in want:
        gens.append(
            ("anchor", iter_anchor_sites(res, chroms, seed=seed, limit=lim["anchor"] or 60))
        )
    rows: list[dict[str, Any]] = []
    import time

    for name, gen in gens:
        t0 = time.time()
        n0 = len(rows)
        for row in gen:
            row["cluster"] = res.clusters.cluster_id(row["chrom"], row["site0"])
            row["partition"] = partition(row["cluster"])
            rows.append(row)
            if (len(rows) - n0) % 200 == 0:
                res.log(f"  {name}: {len(rows) - n0} candidates ...")
        kept = sum(1 for r in rows[n0:] if not r["exclude"])
        res.log(
            f"{name}: {len(rows) - n0} candidates, {kept} not excluded ({time.time() - t0:.0f}s)"
        )
    return rows


def _gnomad_af(chrom: str, pos: int, ref: str, alt: str) -> float | None:
    import asyncio

    from ...clients.gnomad import GnomadClient

    async def go() -> float | None:
        client = GnomadClient()
        try:
            r = await client.lookup(chrom, pos, ref, alt, build="GRCh38")
            return None if r is None else r.global_af
        finally:
            await client.close()

    return asyncio.run(go())


class LabelChecks:
    """Lazy, cached label checks for a row a draw reaches (design 3.0, 3.2, 3.6).

    ``__call__(row)`` returns None when the row may be taken, else the exclusion reason.
    Results are written into ``row["checks"]``. A check whose data handle is absent
    (``hifi`` / ``parents`` not opened) is recorded as ``"not_run"`` and does not exclude,
    so a dry run without them still draws; the manifest then says so per row.
    """

    def __init__(self, res: R4Resources, *, gnomad: bool = True, gnomad_fn: Any = None):
        self.res = res
        self.gnomad = gnomad
        self.gnomad_fn = gnomad_fn or _gnomad_af
        self.counts: Counter[str] = Counter()

    def _key(self, name: str, row: dict[str, Any]) -> str:
        ctx = row.get("context") or {}
        return "|".join(
            str(x)
            for x in (
                name,
                row["stratum"],
                row["chrom"],
                row["site0"],
                ctx.get("candidate_allele") or (row.get("label") or {}).get("alt_allele"),
                ctx.get("delta"),
            )
        )

    def __call__(self, row: dict[str, Any]) -> str | None:
        from .truth import orthogonal_check, trio_recurrence

        res, checks = self.res, row.setdefault("checks", {})
        stratum = row["stratum"]
        if stratum.startswith("anchor_"):
            return None
        if res.hifi is not None:
            hifi = res.cache.memo(
                self._key("hifi", row), lambda: res.retry(lambda: orthogonal_check(row, res.hifi))
            )
            checks["hifi"] = hifi
            if not hifi["pass"]:
                return self._fail(hifi.get("reason") or "hifi")
        else:
            checks["hifi"] = "not_run"
        if row["polarity"] == "neg":
            if res.parents:
                trio = res.cache.memo(
                    self._key("trio", row),
                    lambda: res.retry(lambda: trio_recurrence(row, res.parents, res.parent_truth)),
                )
                checks["trio"] = trio
                if stratum in REQUIRED_TRIO and not trio["recurrent"]:
                    return self._fail("trio_not_recurrent")
            else:
                checks["trio"] = "not_run"
        if stratum == "mosaic_caller_fp_neg" and self.gnomad and row.get("vcf_pos"):
            try:
                af = res.cache.memo(
                    self._key("gnomad", row),
                    lambda: self.gnomad_fn(
                        row["chrom"], row["vcf_pos"], row["vcf_ref"], row["vcf_alt"]
                    ),
                )
            except Exception as exc:  # noqa: BLE001 - an outage must not read as "rare"
                checks["gnomad_af"] = f"unavailable: {type(exc).__name__}"
                return self._fail("gnomad_unavailable")
            checks["gnomad_af"] = af
            from .strata_r4 import GNOMAD_MAX_AF

            if af is not None and af >= GNOMAD_MAX_AF:
                return self._fail("gnomad_af_ge_1pct")
        return None

    def _fail(self, reason: str) -> str:
        self.counts[reason] += 1
        return reason


def pools_report_r4(
    rows: list[dict[str, Any]], *, overrides: dict[str, Any] | None = None
) -> dict[str, Any]:
    """Per stratum: candidates, exclusion counts by rule, and eligible counts by partition
    and knob level, plus the G5 check (main pool >= 2x the main quota) and mean features
    of eligible rows (the pool side of the signature report)."""
    from .strata_r4 import KNOB_LEVELS, QUOTAS_R4, STRATA_R4, eligible

    out: dict[str, Any] = {}
    main_q = QUOTAS_R4["main"]["default"]
    feat_keys = (
        "depth",
        "vaf",
        "alt_reads",
        "alt_mean_bq",
        "alt_mapq_lt20_frac",
        "alt_near_end_frac",
        "alt_softclip_frac",
        "alt_mismatch_burden",
        "alt_proper_pair_frac",
        "linked_mismatch_max_frac",
        "local_nonref_sites_100bp",
        "stutter_frac",
    )
    for s in STRATA_R4:
        rs = [r for r in rows if r["stratum"] == s]
        excl: Counter[str] = Counter(r["exclude"] for r in rs if r["exclude"])
        elig: dict[str, dict[str, int]] = {p: {} for p in ("pilot", "main")}
        fails: Counter[str] = Counter()
        good: list[dict[str, Any]] = []
        for r in rs:
            if r["exclude"]:
                continue
            for lvl in KNOB_LEVELS if not s.startswith("anchor_") else ("default",):
                ok, why = eligible(r, lvl, overrides)
                if ok:
                    elig[r["partition"]][lvl] = elig[r["partition"]].get(lvl, 0) + 1
                    if lvl == "default":
                        good.append(r)
                elif lvl == "default":
                    fails[why or "predicate"] += 1
        means = {}
        for k in feat_keys:
            vals = [float(r["features"][k]) for r in good if r["features"].get(k) is not None]
            if vals:
                means[k] = round(sum(vals) / len(vals), 4)
        main_pool = elig["main"].get("default", 0)
        out[s] = {
            "candidates": len(rs),
            "excluded": dict(excl),
            "predicate_fail_default": dict(fails),
            "eligible": elig,
            "main_quota": main_q[s],
            "g5_main_pool_ge_2x": main_pool >= 2 * main_q[s],
            "clusters_main": len({r["cluster"] for r in good if r["partition"] == "main"}),
            "feature_means_default": means,
        }
    return out


def select_r4(
    rows: list[dict[str, Any]],
    *,
    split: str,
    seed: int,
    checks: Any = None,
    knob_overrides: dict[str, str] | None = None,
    param_overrides: dict[str, Any] | None = None,
    pairs: Iterable[str] | None = None,
    ks_min_p: float = 0.2,
    max_redraws: int = 5,
    log: Any = print,
) -> tuple[list[dict[str, Any]], dict[str, Any]]:
    """Draw ``split`` ("pilot"/"main") from candidate rows with ``seed``.

    Per pair: positives fill their within-stratum cells (:func:`strata_r4.subquotas`) in
    seeded order, then matched negatives are drawn one-to-one (:func:`strata_r4.
    match_negatives`); unmatched strata fill directly. Every take passes the global spacing
    rule (>= 10 kb, <= 2 per cluster) and ``checks`` (lazy label checks; None = none).
    Pilot draws both knob levels (6 default + 4 alt); main uses ``knob_overrides``
    (``{stratum: "alt"}``) on top of the default. The TR pair is redrawn by seed while the
    POI-offset-in-run distributions of strata 3 and 4 differ (KS p <= ``ks_min_p``).
    Returns (manifest rows, report); shortfalls are reported, not raised.
    """
    import functools

    from .strata_r4 import (
        ANCHOR_STRATA,
        MATCHED_TO,
        PAIRS,
        Spacing,
        eligible,
        ks_two_sample,
        label_for,
        match_negatives,
        offset_in_run,
        quotas,
        seeded_order,
        subkey,
        subquotas,
    )

    q = quotas(split)
    knob_overrides = dict(knob_overrides or {})
    check = checks or (lambda r: None)
    spacing = Spacing()
    chosen: list[dict[str, Any]] = []
    report: dict[str, Any] = {
        "split": split,
        "seed": seed,
        "short": {},
        "exclusions": {},
        "matching": {},
    }
    want_pairs = (
        list(pairs) if pairs is not None else ["mosaic", "tr", "segdup", "indel_adj", "anchor"]
    )
    pool = [r for r in rows if not r["exclude"] and r["partition"] == split]
    taken_ids: set[tuple[str, str, int]] = set()

    def rid(r: dict[str, Any]) -> tuple[str, str, int]:
        return (r["stratum"], r["chrom"], int(r["site0"]))

    def level_of(stratum: str, level: str) -> str:
        return knob_overrides.get(stratum, level) if split == "main" else level

    def accept_fn(excl: Counter[str]) -> Any:
        def accept(r: dict[str, Any]) -> str | None:
            if rid(r) in taken_ids:
                return "taken"
            why = spacing.why_not(r)
            if why:
                excl[why] += 1
                return why
            why = check(r)
            if why:
                excl[why] += 1
            return why

        return accept

    def finalize(r: dict[str, Any], level: str) -> dict[str, Any]:
        out = {**r, "knob": level, "label": label_for(r["stratum"], r)}
        out.pop("exclude", None)
        spacing.add(out)
        taken_ids.add(rid(out))
        chosen.append(out)
        return out

    def fill(stratum: str, level: str, n: int, tag: str) -> list[dict[str, Any]]:
        excl = report["exclusions"].setdefault(stratum, Counter())
        accept = accept_fn(excl)
        cells = subquotas(stratum, level, n, param_overrides)
        cand = [
            r for r in pool if r["stratum"] == stratum and eligible(r, level, param_overrides)[0]
        ]
        got: list[dict[str, Any]] = []
        for cell, k in cells.items():
            have = 0
            for r in seeded_order(cand, seed, f"{tag}:{stratum}:{level}:{cell}"):
                if have >= k:
                    break
                if subkey(stratum, r) != cell and cell != "any":
                    continue
                if accept(r) is None:
                    got.append(finalize(r, level))
                    have += 1
            if have < k:
                report["short"][f"{stratum}|{level}|{cell}"] = {"need": k, "got": have}
        return got

    def draw_pair(pair: str, tag: str) -> None:
        strata = [s for s, (p, _) in PAIRS.items() if p == pair]
        pos_s = [s for s in strata if PAIRS[s][1] == "pos"]
        neg_s = [s for s in strata if PAIRS[s][1] == "neg"]
        for level, per in q.items():
            for s in pos_s:
                if s in per:
                    pos_rows = fill(s, level_of(s, level), per[s], tag)
                    for ns in neg_s:
                        if MATCHED_TO.get(ns) != s or ns not in per:
                            continue
                        lvl = level_of(ns, level)
                        cand = [
                            r
                            for r in pool
                            if r["stratum"] == ns and eligible(r, lvl, param_overrides)[0]
                        ]
                        excl = report["exclusions"].setdefault(ns, Counter())
                        need = per[ns]
                        # One negative per drawn positive: a short positive stratum shortens
                        # its negative equally (balance-preserving, design 3.6 fallback 3).
                        target = pos_rows[:need]
                        got, counts = match_negatives(
                            target,
                            cand,
                            ns,
                            seed=_seed_int(seed, tag, level),
                            accept=accept_fn(excl),
                            on_take=functools.partial(finalize, level=lvl),
                        )
                        report["matching"][f"{ns}|{level}"] = {
                            **counts,
                            "levels": Counter(g["match"]["level"] for g in got),
                        }
                        if len(got) < need:
                            report["short"][f"{ns}|{lvl}|matched"] = {"need": need, "got": len(got)}
            for s in neg_s:
                if s in per and s not in MATCHED_TO:
                    fill(s, level_of(s, level), per[s], tag)

    for pair in want_pairs:
        if pair == "anchor":
            for s in ANCHOR_STRATA:
                if s in q["default"]:
                    fill(s, "default", q["default"][s], "anchor")
            continue
        if pair != "tr":
            draw_pair(pair, "draw")
            continue
        snap = (list(spacing.taken), Counter(spacing.per), list(chosen), set(taken_ids))
        for attempt in range(max_redraws + 1):
            draw_pair(pair, f"draw-tr{attempt}" if attempt else "draw")
            got = [r for r in chosen if r["pair"] == "tr"]
            a = [
                x
                for x in (offset_in_run(r) for r in got if r["polarity"] == "pos")
                if x is not None
            ]
            b = [
                x
                for x in (offset_in_run(r) for r in got if r["polarity"] == "neg")
                if x is not None
            ]
            d, p = ks_two_sample(a, b)
            report["tr_offset_ks"] = {"attempt": attempt, "D": round(d, 4), "p": round(p, 4)}
            if p > ks_min_p or attempt == max_redraws:
                break
            spacing.taken, spacing.per = list(snap[0]), Counter(snap[1])
            chosen[:] = list(snap[2])
            taken_ids.clear()
            taken_ids.update(snap[3])
            for k in [k for k in report["short"] if k.startswith("tr_")]:
                del report["short"][k]
    report["exclusions"] = {k: dict(v) for k, v in report["exclusions"].items()}
    for v in report["matching"].values():
        v["levels"] = {str(k): n for k, n in v["levels"].items()}
    report["n"] = Counter(r["stratum"] for r in chosen)
    if report["short"]:
        log(f"WARNING: short strata: {report['short']}")
    chosen.sort(key=lambda r: (r["stratum"], r["chrom"], int(r["site0"])))
    return chosen, report


def _seed_int(seed: int, *parts: object) -> int:
    return int(hashlib.sha256(":".join(map(str, (seed, *parts))).encode()).hexdigest()[:12], 16)


def manifest_to_jsonl(rows: list[dict[str, Any]], path: str | Path) -> str:
    """Write the r4 site manifest (one row per line, sorted keys); returns its sha256."""
    text = "".join(json.dumps(r, sort_keys=True, default=_jsonable) + "\n" for r in rows)
    Path(path).parent.mkdir(parents=True, exist_ok=True)
    Path(path).write_text(text)
    return hashlib.sha256(text.encode()).hexdigest()


def manifest_from_jsonl(path: str | Path) -> list[dict[str, Any]]:
    return [json.loads(line) for line in Path(path).read_text().splitlines() if line.strip()]


def _jsonable(o: Any) -> Any:
    if isinstance(o, (set, frozenset, tuple)):
        return list(o)
    if isinstance(o, Counter):
        return dict(o)
    raise TypeError(f"not JSON serialisable: {type(o).__name__}")


# ---------------------------------------------------------------------------
# Build
# ---------------------------------------------------------------------------

R4_LABEL_KEYS = ("stratum", "pair", "polarity", "cluster", "truth_source", "knob")


def row_to_site(row: dict[str, Any]) -> RealSite:
    return RealSite(
        case_type=row["stratum"],
        chrom=row["chrom"],
        site0=int(row["site0"]),
        label=dict(row["label"]),
        in_repeat=bool((row.get("context") or {}).get("run")),
        vcf_pos=row.get("vcf_pos"),
        vcf_ref=row.get("vcf_ref") or "",
        vcf_alt=row.get("vcf_alt") or "",
        alt_reads=row.get("alt_reads"),
        depth=row.get("depth"),
        extra={k: row.get(k) for k in R4_LABEL_KEYS},
    )


def r4_label_record(item: Item, row: dict[str, Any]) -> dict[str, Any]:
    """Round-3 label record plus the contract's round-4 keys (``stratum`` = case_type)."""
    return {**real_label_record(item), **{k: row.get(k) for k in R4_LABEL_KEYS}}


def build_r4_run_set(
    out: str | Path,
    manifest_path: str | Path,
    *,
    seed: int,
    conditions: list[str] | tuple[str, ...] = R4_CONDITIONS,
    models: list[str],
    repeats: int = 1,
    png: bool | None = None,
    fixture_dir: str | Path | None = None,
    bam_path: str = HG002_BAM,
    bam_index: str | None = None,
    ref_path: str = GRCH38_FASTA,
    log: Any = print,
) -> Any:
    """Blinded fixtures (pairs kept) from a frozen r4 manifest, then stimuli + labels + runs
    via ``stimuli.build_items_run_set(..., r4=True)``."""
    from . import stimuli

    rows = manifest_from_jsonl(manifest_path)
    sites = [row_to_site(r) for r in rows]
    r4_kw: dict[str, Any] = {"r4": True}  # the stimuli agent's r4 switch (CONTRACT.md)
    fx = Path(fixture_dir) if fixture_dir else Path(out) / "fixtures"
    if bam_index is None and bam_path.startswith(("http://", "https://")):
        from .truth import cached

        bam_index = str(cached(bam_path + ".bai"))
    paths, items, mapping = write_real_fixtures(
        sites, fx, seed, bam_path=bam_path, ref_path=ref_path, keep_pairs=True, bam_index=bam_index
    )
    res = stimuli.build_items_run_set(
        out,
        bam=str(paths["bam"]),
        ref=str(paths["reference"]),
        items=items,
        label_records=[r4_label_record(it, r) for it, r in zip(items, rows, strict=True)],
        seed=seed,
        conditions=list(conditions),
        models=models,
        repeats=repeats,
        png=png,
        log=log,
        meta={
            "tier": "real-r4",
            "sites": str(manifest_path),
            "sites_sha256": hashlib.sha256(Path(manifest_path).read_bytes()).hexdigest(),
        },
        sample_line=SAMPLE_LINE_R4,
        sample_kind="sequenced",
        **r4_kw,
    )
    with (Path(out) / "labels" / "mapping.jsonl").open("w") as f:
        for m, r in zip(mapping, rows, strict=True):
            f.write(
                json.dumps({**m, **{k: r.get(k) for k in R4_LABEL_KEYS}}, sort_keys=True) + "\n"
            )
    return res
