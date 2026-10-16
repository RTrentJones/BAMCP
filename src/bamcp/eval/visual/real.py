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
from collections.abc import Iterable
from dataclasses import asdict, dataclass, field
from pathlib import Path
from typing import Any

from .fixtures import Item, _derive_rng, opaque_id

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
MIN_SPACING = 1000  # between selected sites; > 2 * READ_HALF_SPAN
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
    truth_pos = sorted(r.pos - 1 for r in truth)

    def near_truth(p0: int, bp: int, exclude: int | None = None) -> bool:
        i = bisect.bisect_left(truth_pos, p0 - bp)
        while i < len(truth_pos) and truth_pos[i] <= p0 + bp:
            if truth_pos[i] != exclude:
                return True
            i += 1
        return False

    for rec in truth:
        p0 = rec.pos - 1
        if not (a <= p0 < b) or len(rec.alts or ()) != 1:
            continue
        span_end = p0 + max(len(rec.ref), 1)
        if not confident.covers(p0, span_end) or near_truth(p0, ISOLATION_BP, exclude=p0):
            continue
        gt = next(iter(rec.samples.values()))["GT"]
        cls = classify_record(rec.ref, rec.alts[0], tuple(gt))
        if cls is None:
            continue
        label = cls["label"]
        in_rep = repeats.overlaps(p0, span_end + 1)
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
        if near_truth(p0, NEG_CLEAR_BP):
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
        if all(s.chrom != c.chrom or abs(s.site0 - c.site0) >= MIN_SPACING for c in taken):
            chosen.append(s)
            taken.append(s)
        if len(chosen) == k:
            break
    return chosen


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
        for n in range(max_windows):
            s, e = rng.choice(big)
            a = rng.randint(s, e - WINDOW_BP)
            if any(abs(a - u) < WINDOW_BP for u in used):
                continue
            used.append(a)
            for k, v in scan_window(
                bam, fasta, vcf, chrom, a, a + WINDOW_BP, confident, repeats
            ).items():
                pools.setdefault(k, []).extend(v)
            if (n + 1) % 50 == 0:
                log(f"{n + 1} windows: { {k: len(v) for k, v in pools.items()} }")
            if all(len(pools.get(k, [])) >= pool_factor * q for k, q in quotas.items()):
                log(f"pools full after {n + 1} windows")
                break
        else:
            sizes = {k: len(v) for k, v in pools.items()}
            log(f"scanned {max_windows} windows; pool sizes: {sizes}")

    chosen: list[RealSite] = []
    pick_rng = _derive_rng(seed, "real-pick")
    for k, q in quotas.items():
        pool = sorted(pools.get(k, []), key=lambda s: (s.chrom, s.site0))
        if k == "real_noisy_negative":
            rep = [s for s in pool if s.in_repeat]
            non = [s for s in pool if not s.in_repeat]
            half = _pick(rep, q // 2, pick_rng, chosen)
            got = half + _pick(non, q - len(half), pick_rng, chosen)
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


def write_real_fixtures(
    sites: list[RealSite],
    out_dir: str | Path,
    seed: int,
    *,
    bam_path: str = HG002_BAM,
    ref_path: str = GRCH38_FASTA,
) -> tuple[dict[str, Path], list[Item], list[dict[str, Any]]]:
    """Write the blinded BAM + reference and return (paths, items, mapping).

    Kept reads: primary, non-duplicate, non-QC-fail alignments overlapping the site +/-
    ``READ_HALF_SPAN`` whose whole alignment fits the slot. They are written single-end
    (flag 0/16), with hashed names and no tags, so every condition counts the same reads.
    ``mapping`` (original -> blinded coordinates) is label-bearing and stays out of stim/.
    """
    import pysam

    out = Path(out_dir)
    out.mkdir(parents=True, exist_ok=True)
    layout = _layout(sites, seed)
    per = -(-len(sites) // N_CONTIGS)
    contig_len = 2 * CONTIG_PAD + per * SLOT_BP
    names = [f"chr{c + 1}" for c in range(N_CONTIGS)]
    refs = {n: ["N"] * contig_len for n in names}
    records: list[tuple[int, int, str, int, str, Any, Any, int]] = []
    items: list[Item] = []
    mapping: list[dict[str, Any]] = []
    with pysam.AlignmentFile(bam_path) as bam, pysam.FastaFile(ref_path) as fasta:
        for idx, (site, (contig, slot0)) in enumerate(zip(sites, layout, strict=True)):
            src0 = site.site0 - SLOT_BP // 2
            shift = slot0 - src0
            seq = fasta.fetch(site.chrom, src0, src0 + SLOT_BP).upper()
            refs[contig][slot0 : slot0 + len(seq)] = list(seq)
            item_id = "i" + opaque_id(seed, "real", site.chrom, site.site0)
            n = 0
            for r in bam.fetch(
                site.chrom, site.site0 - READ_HALF_SPAN, site.site0 + READ_HALF_SPAN
            ):
                if not _usable(r) or r.reference_end is None:
                    continue
                if r.reference_start < src0 or r.reference_end > src0 + SLOT_BP:
                    continue
                name = hashlib.sha1(f"{seed}:{item_id}:{n}".encode()).hexdigest()[:16]  # noqa: S324
                n += 1
                records.append(
                    (
                        names.index(contig),
                        r.reference_start + shift,
                        name,
                        16 if r.is_reverse else 0,
                        r.query_sequence or "",
                        r.query_qualities,
                        r.cigartuples,
                        r.mapping_quality,
                    )
                )
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
        "SQ": [{"SN": nm, "LN": contig_len} for nm in names],
    }
    records.sort(key=lambda t: (t[0], t[1], t[2]))
    with pysam.AlignmentFile(str(bam_out), "wb", header=header) as w:
        for tid, start, name, flag, seq, quals, cigar, mapq in records:
            a = pysam.AlignedSegment(w.header)
            a.query_name = name
            a.query_sequence = seq
            a.flag = flag
            a.reference_id = tid
            a.reference_start = start
            a.mapping_quality = mapq
            a.cigartuples = cigar
            a.query_qualities = quals
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
