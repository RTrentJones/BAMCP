"""Round-4 truth handling: truth sets, concordance, orthogonal (HiFi) and trio label checks.

Round 4 (``docs/eval/visual-design-r4.md`` section 3) takes its sites from several GIAB
releases at once, each authoritative for a different variant class:

- **v4.2.1** and **v5.0q** small-variant benchmarks (germline; blind below about 30% VAF);
- **CMRG v1.00** (curated assembly truth in challenging medically relevant genes);
- **Mosaic SNV v1.1** (85 SNVs at 5-30% VAF, plus a BED in which any unlisted SNV is
  asserted to be below 2% VAF).

:class:`TruthSet` loads one VCF + BED pair per chromosome. It is FILTER-aware: positives use
PASS records inside the BED only, while isolation and negative clearance count records of
any FILTER, inside or outside the BED (v5.0q carries tens of thousands of non-PASS records).
Distances are always measured from a record's **full reference span**, as ``near_truth``
in :mod:`.real`: a long deletion anchored upstream still covers every base it deletes.

Concordance is class-aware (:func:`concordant`): two truth sets are compared only where
both assert a germline genotype (the site lies in both BEDs), after left-normalisation, so a
v4.2.1 reference call is never read as evidence about a mosaic site.

Label checks against independent data, frozen before the pilot:

- :func:`orthogonal_check` — HG002 PacBio HiFi at the site (alt support for positives,
  absence for negatives, spanning-read length histograms for repeat indels);
- :func:`trio_recurrence` — the same non-reference allele in the parents' 60x Illumina
  BAMs (same NHGRI Novoalign pipeline) marks a locus-specific technical artifact.

:func:`cluster_index` builds paralog clusters from the UCSC ``genomicSuperDups`` table
(connected components of segments and their partners; a 1 Mb bin outside segdups), and
:func:`partition` assigns pilot/main by a seeded hash of the cluster ID, so paralog
families never straddle the two splits. Downloads go to a user cache directory
(:func:`cache_dir`), never into the repository.
"""

from __future__ import annotations

import bisect
import gzip
import hashlib
import os
import shutil
import urllib.error
import urllib.request
from collections import Counter
from collections.abc import Callable, Iterable, Iterator
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

from .real import GIAB_BASE, HG002_BAM, TRUTH_VCF, Intervals, _usable

# ---------------------------------------------------------------------------
# Public data
# ---------------------------------------------------------------------------

FTP = "https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/"
HG002_RELEASE = FTP + "release/AshkenazimTrio/HG002_NA24385_son/"
V421_VCF = TRUTH_VCF
V421_BED = GIAB_BASE + "HG002_GRCh38_1_22_v4.2.1_benchmark_noinconsistent.bed"
V5Q_VCF = HG002_RELEASE + "v5.0q/HG002_GRCh38_v5.0q_smvar.vcf.gz"
V5Q_BED = HG002_RELEASE + "v5.0q/HG002_GRCh38_v5.0q_smvar.benchmark.bed"
V5Q_SV_VCF = HG002_RELEASE + "v5.0q/HG002_GRCh38_v5.0q_stvar.vcf.gz"
CMRG_BASE = HG002_RELEASE + "CMRG_v1.00/GRCh38/"
CMRG_VCF = CMRG_BASE + "SmallVariant/HG002_GRCh38_CMRG_smallvar_v1.00.vcf.gz"
CMRG_BED = CMRG_BASE + "SmallVariant/HG002_GRCh38_CMRG_smallvar_v1.00.bed"
CMRG_SV_VCF = CMRG_BASE + "StructuralVariant/HG002_GRCh38_CMRG_SV_v1.00.vcf.gz"
MOSAIC_BASE = HG002_RELEASE + "mosaic_v1.10/GRCh38/SNV/"
MOSAIC_VCF = MOSAIC_BASE + "HG002_GRCh38_MosaicSNVv1.1.vcf.gz"
MOSAIC_BED = MOSAIC_BASE + "HG002_GRCh38_MosaicSNVv1.1.bed"
_MOSAIC_SUPP = MOSAIC_BASE + "SupplementaryFiles/"
_VALIDATORS = _MOSAIC_SUPP + "external_validator_vcfs/"
# Caller output the mosaic benchmark rejected or never listed: the sources of stratum 2's
# caller-flagged negatives (each minus the mosaic VCF). All but Strelka2 are tabix-indexed
# on the server; the Strelka2 VCF and the annotated database are whole-file downloads.
MOSAIC_CALLER_VCFS = {
    "dragen_mosaic": _VALIDATORS + "dragen/HG002_300x.mosaic-only-hard-filtered.vcf.gz",
    "deepsomatic_ill300x_q0": _VALIDATORS
    + "deepsomatic_google/HG002.illumina.300x.deepsomatic.mosaic.Q0.vcf.gz",
    "strelka2_ill300x": _MOSAIC_SUPP + "HG002_GRCh38_Strelka2-Ill300X.vcf.gz",
}
MOSAIC_DATABASE_TSV = _MOSAIC_SUPP + "HG002_GRCh38_mosaic_database.tsv"
# DRAGEN germline calls on the 300x data: inside the CMRG BED, minus CMRG truth, these are
# the short-read false positives a user would actually review (stratum 6's primary source).
DRAGEN_300X_VCF = _VALIDATORS + "dragen/HG002_300x.hard-filtered.vcf.gz"

_DATA = FTP + "data/AshkenazimTrio/"
HIFI_BAM = (
    _DATA + "HG002_NA24385_son/PacBio_CCS_15kb_20kb_chemistry2/GRCh38/"
    "HG002.SequelII.merged_15kb_20kb.pbmm2.GRCh38.haplotag.10x.bam"
)
PARENT_BAMS = {
    "HG003": _DATA + "HG003_NA24149_father/NIST_HiSeq_HG003_Homogeneity-12389378/"
    "NHGRI_Illumina300X_AJtrio_novoalign_bams/HG003.GRCh38.60x.1.bam",
    "HG004": _DATA + "HG004_NA24143_mother/NIST_HiSeq_HG004_Homogeneity-14572558/"
    "NHGRI_Illumina300X_AJtrio_novoalign_bams/HG004.GRCh38.60x.1.bam",
}
_PARENT_TRUTH = (
    FTP + "release/AshkenazimTrio/{d}/NISTv4.2.1/GRCh38/{s}_GRCh38_1_22_v4.2.1_benchmark"
)
PARENT_TRUTH = {
    s: (
        _PARENT_TRUTH.format(d=d, s=s) + ".vcf.gz",
        _PARENT_TRUTH.format(d=d, s=s) + "_noinconsistent.bed",
    )
    for s, d in (("HG003", "HG003_NA24149_father"), ("HG004", "HG004_NA24143_mother"))
}
SUPERDUPS_URL = "https://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/genomicSuperDups.txt.gz"
_STRAT = FTP + "release/genome-stratifications/v3.3/GRCh38@all/"
STRAT_URLS = {
    "segdups": _STRAT + "SegmentalDuplications/GRCh38_segdups.bed.gz",
    "lowmappability": _STRAT + "Mappability/GRCh38_lowmappabilityall.bed.gz",
    "homopolymer_ge12": _STRAT + "LowComplexity/GRCh38_SimpleRepeat_homopolymer_ge12_slop5.bed.gz",
    "homopolymer_7to11": _STRAT
    + "LowComplexity/GRCh38_SimpleRepeat_homopolymer_7to11_slop5.bed.gz",
    "diTR": _STRAT + "LowComplexity/GRCh38_SimpleRepeat_diTR_10to49_slop5.bed.gz",
    "complexindel10bp": _STRAT
    + "GenomeSpecific/GRCh38_HG002_v4.2.1_complexindel10bp_slop50.bed.gz",
    "snpswithin10bp": _STRAT + "GenomeSpecific/GRCh38_HG002_v4.2.1_snpswithin10bp_slop50.bed.gz",
}
AUTOSOMES = tuple(f"chr{i}" for i in range(1, 23))

# Orthogonal-check thresholds (design 3.0), frozen before the pilot.
HIFI_MIN_MAPQ = 20
HIFI_MIN_ALT = 3
HIFI_ALT_FRAC = {"low_vaf": 0.03, "het": 0.25, "hom_alt": 0.80}
HIFI_NEG_MAX_READS = 1
HIFI_NEG_MAX_FRAC = 0.02
HIFI_STUTTER_RATIO = 0.5  # HiFi +/-1-unit fraction must be <= this x the Illumina fraction
HIFI_MIN_DEPTH = 10
HIFI_HOM_REF_RATIO = 0.5
TRIO_RECUR_FRAC = 0.05
SPAN_FLANK = 3  # a read "spans" a repeat run when it aligns this many bp past both ends
BIN_BP = 1_000_000
PILOT_FRACTION = 0.25
PARTITION_SALT = "r4-partition-v1"  # frozen: pilot and main draws share one partition


# ---------------------------------------------------------------------------
# Cache
# ---------------------------------------------------------------------------


def cache_dir() -> Path:
    """``$BAMCP_EVAL_CACHE``, else ``$XDG_CACHE_HOME/bamcp-eval``, else ``~/.cache/bamcp-eval``."""
    env = os.environ.get("BAMCP_EVAL_CACHE")
    if env:
        d = Path(env)
    else:
        base = os.environ.get("XDG_CACHE_HOME") or str(Path.home() / ".cache")
        d = Path(base) / "bamcp-eval"
    d.mkdir(parents=True, exist_ok=True)
    return d


def is_remote(path: str | Path) -> bool:
    return str(path).startswith(("http://", "https://", "ftp://"))


def cached(url: str | Path, name: str | None = None) -> Path:
    """A local copy of ``url`` (downloaded once into :func:`cache_dir`); local paths pass
    through unchanged. The name is prefixed with a short URL hash, so two releases that
    share a basename never collide."""
    if not is_remote(url):
        return Path(url)
    url = str(url)
    digest = hashlib.sha256(url.encode()).hexdigest()[:10]
    dest = cache_dir() / f"{digest}_{name or url.rsplit('/', 1)[-1]}"
    if not dest.exists():
        tmp = dest.with_name(f"{dest.name}.{os.getpid()}.part")  # concurrent runs never share
        for attempt in range(DOWNLOAD_RETRIES):
            try:
                _download(url, tmp)
                break
            except OSError as exc:
                tmp.unlink(missing_ok=True)
                missing = isinstance(exc, urllib.error.HTTPError) and exc.code == 404
                if missing or attempt == DOWNLOAD_RETRIES - 1:
                    raise
        tmp.replace(dest)
    return dest


DOWNLOAD_RETRIES = 3


def _download(url: str, dest: Path) -> None:
    """Fetch ``url`` to ``dest``; a short body (fewer bytes than Content-Length, as when a
    proxy drops the connection) raises instead of leaving a silently truncated file."""
    with urllib.request.urlopen(url, timeout=600) as resp, dest.open("wb") as fh:  # noqa: S310
        shutil.copyfileobj(resp, fh)
        want = resp.headers.get("Content-Length") if hasattr(resp, "headers") else None
    got = dest.stat().st_size
    if want is not None and int(want) != got:
        raise OSError(f"short download of {url}: {got} of {want} bytes")


def open_bam(path: str, index: str | Path | None = None) -> Any:
    """A pysam AlignmentFile; a remote BAM's ``.bai`` is fetched into the cache rather than
    the working directory (htslib's default)."""
    import pysam

    if index is None and is_remote(path):
        index = cached(path + ".bai")
    return pysam.AlignmentFile(path, index_filename=str(index) if index else None)


def fetch_vcf(holder: list[Any], path: str | Path, *args: Any, tries: int = 4) -> list[Any]:
    """``list(vf.fetch(*args))`` for the VariantFile in ``holder[0]``, re-opening ``path``
    into ``holder`` and retrying when a remote read drops (``OSError``)."""
    import time

    for attempt in range(tries):
        try:
            return list(holder[0].fetch(*args))
        except OSError:
            if attempt == tries - 1 or not is_remote(path):
                raise
            time.sleep(2 * (attempt + 1))
            holder[0] = open_vcf(path)
    raise AssertionError("unreachable")  # pragma: no cover


def open_vcf(path: str | Path) -> Any:
    """A pysam VariantFile, its remote ``.tbi`` cached like :func:`open_bam`'s index."""
    import pysam

    p = str(path)
    if is_remote(p):
        return pysam.VariantFile(p, index_filename=str(cached(p + ".tbi")))
    return pysam.VariantFile(p)


def load_bed_multi(path: str | Path, chroms: Iterable[str] | None = None) -> dict[str, Intervals]:
    """One :class:`Intervals` per chromosome (all chromosomes when ``chroms`` is None)."""
    p = cached(path)
    want = set(chroms) if chroms is not None else None
    ivs: dict[str, list[tuple[int, int]]] = {}
    opener: Any = gzip.open if p.suffix == ".gz" else open
    with opener(p, "rt") as fh:
        for line in fh:
            if not line.strip() or line.startswith(("#", "track", "browser")):
                continue
            f = line.split("\t", 3)
            if want is None or f[0] in want:
                ivs.setdefault(f[0], []).append((int(f[1]), int(f[2])))
    return {c: Intervals(v) for c, v in ivs.items()}


class BedSet:
    """Merged intervals for several chromosomes; missing chromosomes are empty."""

    def __init__(self, by_chrom: dict[str, Intervals]):
        self.by_chrom = by_chrom
        self._empty = Intervals([])

    @classmethod
    def load(cls, path: str | Path, chroms: Iterable[str] | None = None) -> BedSet:
        return cls(load_bed_multi(path, chroms))

    def get(self, chrom: str) -> Intervals:
        return self.by_chrom.get(chrom, self._empty)

    def contains(self, chrom: str, pos0: int) -> bool:
        return self.get(chrom).contains(pos0)

    def overlaps(self, chrom: str, s: int, e: int) -> bool:
        return self.get(chrom).overlaps(s, e)

    def covers(self, chrom: str, s: int, e: int) -> bool:
        return self.get(chrom).covers(s, e)

    def total_bp(self) -> int:
        return sum(e - s for iv in self.by_chrom.values() for s, e in iv.pairs())


# ---------------------------------------------------------------------------
# Truth records and sets
# ---------------------------------------------------------------------------


@dataclass(frozen=True)
class TruthRecord:
    """One VCF record. ``pos`` is the 1-based VCF POS; spans are 0-based half-open."""

    chrom: str
    pos: int
    ref: str
    alts: tuple[str, ...]
    gt: tuple[int | None, ...]
    filter: str
    source: str
    info: tuple[tuple[str, Any], ...] = ()

    @property
    def start0(self) -> int:
        return self.pos - 1

    @property
    def end0(self) -> int:
        """End of the full reference span (an insertion's span is its anchor base)."""
        return self.start0 + max(len(self.ref), 1)

    @property
    def is_pass(self) -> bool:
        return self.filter in ("PASS", ".", "")

    @property
    def biallelic(self) -> bool:
        return len(self.alts) == 1

    @property
    def kind(self) -> str:
        """``SNV``, ``ins``, ``del``, ``mnp``, ``complex`` or ``multi`` (several ALTs)."""
        if len(self.alts) != 1:
            return "multi"
        return variant_kind(self.ref, self.alts[0])

    @property
    def zygosity(self) -> str:
        calls = [a for a in self.gt if a is not None]
        if len(calls) != len(self.gt) or not calls:
            return "unknown"
        if all(a == 0 for a in calls):
            return "hom_ref"
        if len(set(calls)) == 1:
            return "hom_alt"
        return "het" if 0 in calls else "multi"

    def overlaps(self, s0: int, e0: int, bp: int = 0) -> bool:
        return self.start0 < e0 + bp and self.end0 > s0 - bp


def variant_kind(ref: str, alt: str) -> str:
    if len(ref) == len(alt) == 1:
        return "SNV"
    if len(ref) == len(alt):
        return "mnp"
    if len(ref) > len(alt) == 1 and ref[0] == alt[0]:
        return "del"
    if len(alt) > len(ref) == 1 and ref[0] == alt[0]:
        return "ins"
    return "complex"


def _record(rec: Any, source: str) -> TruthRecord:
    filt = ";".join(rec.filter.keys()) or "."
    gt: tuple[int | None, ...] = ()
    if rec.samples:
        try:
            gt = tuple(next(iter(rec.samples.values()))["GT"] or ())
        except KeyError:
            gt = ()
    return TruthRecord(rec.chrom, rec.pos, rec.ref, tuple(rec.alts or ()), gt, filt, source)


class TruthSet:
    """A truth VCF (+ optional BED).

    Records keep every FILTER value. :meth:`positives` yields PASS records inside the BED;
    :meth:`records_near` (isolation, negative clearance) counts records of any FILTER,
    inside or outside the BED. Chromosomes in ``chroms`` are loaded whole; for an indexed
    VCF any other chromosome is answered by region queries (so a genome-wide set such as
    v5.0q is never held in memory for a handful of mosaic sites), and loaded whole only if
    :meth:`records` / :meth:`positives` ask for it. An unindexed VCF (e.g. the Strelka2
    release file) is read once. The BED is always loaded for every chromosome.
    """

    REGION_PAD = 500  # longest small-variant span a region query must still see

    def __init__(
        self,
        name: str,
        records: Iterable[TruthRecord],
        bed: BedSet | None = None,
        *,
        vcf_handle: Any = None,
        vcf_path: str | Path | None = None,
    ):
        self.name = name
        self.bed = bed
        self._vf = vcf_handle
        self._vcf_path = vcf_path
        self.by_chrom: dict[str, list[TruthRecord]] = {}
        self._starts: dict[str, list[int]] = {}
        self._max_len: dict[str, int] = {}
        grouped: dict[str, list[TruthRecord]] = {}
        for r in records:
            grouped.setdefault(r.chrom, []).append(r)
        for c, recs in grouped.items():
            self._index(c, recs)

    def _index(self, chrom: str, recs: list[TruthRecord]) -> None:
        recs.sort(key=lambda r: (r.start0, r.end0, r.alts))
        self.by_chrom[chrom] = recs
        self._starts[chrom] = [r.start0 for r in recs]
        self._max_len[chrom] = max((r.end0 - r.start0 for r in recs), default=1)

    @classmethod
    def load(
        cls,
        name: str,
        vcf: str | Path,
        bed: str | Path | None = None,
        *,
        chroms: Iterable[str] | None = None,
    ) -> TruthSet:
        """``chroms=None`` loads every record."""
        want = list(chroms) if chroms is not None else None
        records: list[TruthRecord] = []
        src: str | Path = vcf
        if is_remote(vcf):
            try:
                cached(str(vcf) + ".tbi")
            except OSError:  # no server-side index: read the whole (small) file locally
                src = cached(vcf)
        vf = open_vcf(src)
        handle = None
        try:
            if want is None:
                raise ValueError("whole file")
            contigs = set(vf.header.contigs)
            holder = [vf]
            for c in want:
                if c in contigs:
                    records.extend(_record(r, name) for r in fetch_vcf(holder, src, c))
            handle = vf = holder[0]
        except ValueError:  # unindexed (or all chromosomes wanted): one pass
            records = [_record(r, name) for r in vf if want is None or r.chrom in want]
            vf.close()
        ts = cls(
            name,
            records,
            BedSet.load(bed) if bed is not None else None,
            vcf_handle=handle,
            vcf_path=src,
        )
        for c in want or ():
            if handle is not None and c not in ts.by_chrom:
                ts._index(c, [])  # loaded and empty, not "query by region"
        return ts

    def chroms(self) -> list[str]:
        return sorted(self.by_chrom)

    def _ensure(self, chrom: str) -> None:
        if chrom in self.by_chrom or self._vf is None:
            return
        recs = []
        if chrom in set(self._vf.header.contigs):
            recs = [_record(r, self.name) for r in self._region(chrom)]
        self._index(chrom, recs)

    def _region(self, *args: Any) -> list[Any]:
        holder = [self._vf]
        got = fetch_vcf(holder, self._vcf_path or "", *args)
        self._vf = holder[0]
        return got

    def records(self, chrom: str) -> list[TruthRecord]:
        self._ensure(chrom)
        return self.by_chrom.get(chrom, [])

    def covers(self, chrom: str, s0: int, e0: int | None = None) -> bool:
        """The BED covers [s0, e0) (one base when ``e0`` is None). No BED: False."""
        if self.bed is None:
            return False
        return self.bed.covers(chrom, s0, s0 + 1 if e0 is None else e0)

    def records_near(
        self, chrom: str, s0: int, e0: int, bp: int, *, exclude: TruthRecord | None = None
    ) -> list[TruthRecord]:
        """Records of any FILTER whose full reference span comes within ``bp`` of [s0, e0)."""
        recs = self.by_chrom.get(chrom)
        if recs is None and self._vf is not None:
            if chrom not in set(self._vf.header.contigs):
                return []
            lo = max(0, s0 - bp - self.REGION_PAD)
            got = [_record(r, self.name) for r in self._region(chrom, lo, e0 + bp)]
            return [r for r in got if r.overlaps(s0, e0, bp) and r != exclude]
        if not recs:
            return []
        starts = self._starts[chrom]
        lo = bisect.bisect_left(starts, s0 - bp - self._max_len[chrom])
        hi = bisect.bisect_left(starts, e0 + bp)
        return [r for r in recs[lo:hi] if r.overlaps(s0, e0, bp) and r is not exclude]

    def positives(self, chrom: str) -> Iterator[TruthRecord]:
        """PASS records whose whole span lies inside the BED (all PASS records if no BED)."""
        for r in self.records(chrom):
            if r.is_pass and (self.bed is None or self.covers(chrom, r.start0, r.end0)):
                yield r

    def asserts_reference(self, chrom: str, s0: int, e0: int, bp: int = 0) -> bool:
        """The set says the sample is reference over [s0, e0): BED covers it, no record near."""
        return self.covers(chrom, s0, e0) and not self.records_near(chrom, s0, e0, bp)

    def find(self, chrom: str, pos: int, ref: str, alt: str) -> TruthRecord | None:
        for r in self.records_near(chrom, pos - 1, pos, 0):
            if r.pos == pos and r.ref == ref and alt in r.alts:
                return r
        return None


# ---------------------------------------------------------------------------
# Normalisation and concordance
# ---------------------------------------------------------------------------

Fetch = Callable[[str, int, int], str]


def left_normalize(fetch: Fetch, chrom: str, pos: int, ref: str, alt: str) -> tuple[int, str, str]:
    """``bcftools norm``-style left-normalisation of one biallelic record (1-based ``pos``).

    Trims shared trailing bases (extending left from the reference when an allele would
    become empty), then shared leading bases while both alleles keep at least one base.
    ``fetch(chrom, start0, end0)`` returns reference bases (any case).
    """
    ref, alt = ref.upper(), alt.upper()
    if ref == alt:
        return pos, ref, alt
    while True:
        changed = False
        if ref and alt and ref[-1] == alt[-1]:
            ref, alt = ref[:-1], alt[:-1]
            changed = True
        if not ref or not alt:
            if pos <= 1:
                break
            pos -= 1
            base = fetch(chrom, pos - 1, pos).upper()
            ref, alt = base + ref, base + alt
            changed = True
        if not changed:
            break
    while len(ref) > 1 and len(alt) > 1 and ref[0] == alt[0]:
        ref, alt = ref[1:], alt[1:]
        pos += 1
    return pos, ref, alt


def normalized_key(fetch: Fetch, rec: TruthRecord, alt: str | None = None) -> tuple[Any, ...]:
    a = alt if alt is not None else rec.alts[0]
    return (rec.chrom, *left_normalize(fetch, rec.chrom, rec.pos, rec.ref, a))


def concordant(rec: TruthRecord, other: TruthSet, fetch: Fetch, *, window: int = 50) -> bool | None:
    """Do ``rec`` (a germline call) and ``other`` agree?

    ``None`` when ``other`` does not assert a genotype over the record's span (its BED does
    not cover it): no comparison is made. Otherwise True only if ``other`` has a record with
    the same left-normalised allele **and** zygosity; a different allele, a different
    zygosity, or no record at all (``other`` asserts reference) is discordant. Mosaic-class
    records must not be passed here (v4.2.1/v5.0q cannot see < 30% VAF).
    """
    if not other.covers(rec.chrom, rec.start0, rec.end0):
        return None
    if not rec.alts:
        return False
    want = normalized_key(fetch, rec)
    for o in other.records_near(rec.chrom, rec.start0, rec.end0, window):
        for i, a in enumerate(o.alts):
            if normalized_key(fetch, o, a) != want:
                continue
            z = o.zygosity
            if z == "multi":  # 1/2: compare presence only of this allele's copy number
                n = sum(1 for g in o.gt if g == i + 1)
                z = "hom_alt" if n == 2 else "het"
            return z == rec.zygosity
    return False


def reference_concordant(other: TruthSet, chrom: str, s0: int, e0: int, bp: int) -> bool | None:
    """For a negative over [s0, e0): ``None`` if ``other`` does not cover it, else whether
    ``other`` also says reference (no record of any FILTER within ``bp``)."""
    if not other.covers(chrom, s0, e0):
        return None
    return not other.records_near(chrom, s0, e0, bp)


# ---------------------------------------------------------------------------
# Repeat runs (reference only)
# ---------------------------------------------------------------------------


def repeat_run(seq: str, i: int, max_unit: int = 3) -> tuple[int, int, str]:
    """The longest tandem run (unit 1..``max_unit``) of ``seq`` that contains index ``i``.

    Returns ``(start, end, unit)`` within ``seq`` (half-open). A run needs at least two
    copies of its unit; a homopolymer of any length counts. Ties prefer the shorter unit.
    """
    seq = seq.upper()
    best = (i, i + 1, seq[i] if 0 <= i < len(seq) else "")
    for k in range(1, max_unit + 1):
        for phase in range(k):
            s = i - phase
            if s < 0 or s + k > len(seq):
                continue
            unit = seq[s : s + k]
            if "N" in unit or (k > 1 and len(set(unit)) == 1):
                continue  # a homopolymer is found at k = 1
            a = s
            while a - k >= 0 and seq[a - k : a] == unit:
                a -= k
            b = s + k
            while b + k <= len(seq) and seq[b : b + k] == unit:
                b += k
            # A partial trailing copy still belongs to the run (e.g. ACACA).
            while b < len(seq) and seq[b] == unit[(b - a) % k]:
                b += 1
            if (b - a) >= 2 * k and b - a > best[1] - best[0] and a <= i < b:
                best = (a, b, seq[a : a + k] if k > 1 else unit)
    return best


def is_tr_run(length: int, unit: str, hp_min: int, str_min: int = 20) -> bool:
    """A homopolymer >= ``hp_min`` bp, or a di/tri STR >= ``str_min`` bp (design 3.3)."""
    if len(unit) == 1:
        return length >= hp_min
    return len(unit) in (2, 3) and length >= str_min


# ---------------------------------------------------------------------------
# Read-level counting shared by the checks
# ---------------------------------------------------------------------------


def base_counts(bam: Any, chrom: str, pos0: int, *, min_mapq: int = 0) -> Counter[str]:
    """Bases (and ``-`` for a deletion) at ``pos0`` over primary, non-duplicate reads."""
    c: Counter[str] = Counter()
    for col in bam.pileup(
        chrom,
        pos0,
        pos0 + 1,
        truncate=True,
        min_base_quality=0,
        stepper="nofilter",
        ignore_orphans=False,
        ignore_overlaps=False,
        max_depth=100_000,
    ):
        for pr in col.pileups:
            a = pr.alignment
            if not _usable(a) or a.mapping_quality < min_mapq:
                continue
            if pr.is_del or pr.is_refskip or pr.query_position is None:
                c["-"] += 1
                continue
            q = a.query_sequence
            c[(q[pr.query_position] if q else "N").upper()] += 1
    return c


def read_run_delta(read: Any, s0: int, e0: int) -> int | None:
    """Net indel length (insertions positive) a read carries inside the run [s0, e0),
    with one base of tolerance at each edge; None if the read does not span the run with
    ``SPAN_FLANK`` aligned bases on each side."""
    if read.reference_start > s0 - SPAN_FLANK or (read.reference_end or 0) < e0 + SPAN_FLANK:
        return None
    delta = 0
    ref = read.reference_start
    for op, n in read.cigartuples or ():
        if op in (0, 7, 8):  # M, =, X
            ref += n
        elif op == 1:  # insertion sits before reference base `ref`
            if s0 - 1 <= ref <= e0 + 1:
                delta += n
        elif op in (2, 3):  # deletion / skip
            if s0 - 1 <= ref <= e0:
                delta -= n
            ref += n
    return delta


def spanning_hist(bam: Any, chrom: str, s0: int, e0: int, *, min_mapq: int = 0) -> Counter[int]:
    """Histogram of :func:`read_run_delta` over reads spanning the run."""
    hist: Counter[int] = Counter()
    for r in bam.fetch(chrom, s0, e0):
        if not _usable(r) or r.mapping_quality < min_mapq:
            continue
        d = read_run_delta(r, s0, e0)
        if d is not None:
            hist[d] += 1
    return hist


def unit_fraction(hist: Counter[int], unit_len: int, exclude: Iterable[int] = (0,)) -> float:
    """Fraction of spanning reads at +/-1 repeat unit, not counting ``exclude`` deltas."""
    n = sum(hist.values())
    if not n:
        return 0.0
    skip = set(exclude)
    return sum(hist[d] for d in (-unit_len, unit_len) if d not in skip) / n


# ---------------------------------------------------------------------------
# Orthogonal (HiFi) check
# ---------------------------------------------------------------------------


def _row_allele(row: dict[str, Any]) -> tuple[str, str]:
    """(kind, allele) of the candidate in a manifest row: ``("snv", "G")`` or
    ``("tr", "<delta>")`` for a repeat-length allele."""
    ctx = row.get("context") or {}
    if ctx.get("run") is not None and ctx.get("delta") is not None:
        return "tr", str(ctx["delta"])
    lab = row.get("label") or {}
    allele = ctx.get("candidate_allele") or lab.get("alt_allele") or row.get("vcf_alt") or ""
    return "snv", str(allele).upper()


def orthogonal_check(row: dict[str, Any], bam: Any) -> dict[str, Any]:
    """HiFi check of one manifest row (design 3.0); returns ``{"pass": bool, ...}``.

    - SNV positives: alt reads at MAPQ >= 20 >= 3 and fraction >= 0.03 (``low_vaf``),
      0.25 (``het``) or 0.80 (``hom_alt``).
    - Repeat-indel positives: the spanning-read length histogram's modal length(s) sit at
      the truth genotype (het: the two top lengths are {0, delta}; hom: the top is delta
      and reference-length reads are fewer than half as many).
    - Negatives: the candidate allele on <= 1 HiFi read and <= 2%; for repeat-length
      candidates (stratum 4), the HiFi +/-1-unit fraction must also be <= half the
      Illumina fraction (``features["stutter_frac"]``).
    """
    kind, allele = _row_allele(row)
    pos = row.get("polarity") == "pos"
    chrom, site0 = row["chrom"], int(row["site0"])
    out: dict[str, Any] = {"source": "hifi"}
    if kind == "tr":
        ctx = row["context"]
        s0, e0 = ctx["run"]
        delta = int(ctx["delta"])
        unit = len(ctx.get("unit") or "A")
        hist = spanning_hist(bam, chrom, s0, e0, min_mapq=HIFI_MIN_MAPQ)
        n = sum(hist.values())
        out.update(hist={str(k): v for k, v in sorted(hist.items())}, depth=n)
        if n < HIFI_MIN_DEPTH:
            return {**out, "pass": False, "reason": "hifi_low_depth"}
        top = [d for d, _ in hist.most_common()]
        if pos:
            zyg = (row.get("label") or {}).get("zygosity")
            if zyg == "hom_alt":  # the truth length is the mode and reference is not a co-mode
                ok = top[0] == delta and hist[0] < HIFI_HOM_REF_RATIO * hist[delta]
            else:
                ok = set(top[:2]) == {0, delta} if len(top) >= 2 else False
            return {**out, "pass": ok, "reason": None if ok else "hifi_mode_mismatch"}
        frac = hist[delta] / n
        ill = float((row.get("features") or {}).get("stutter_frac") or 0.0)
        hifi_unit = unit_fraction(hist, unit)
        out.update(alt_reads=hist[delta], alt_frac=frac, unit_frac=hifi_unit)
        ok = hist[delta] <= HIFI_NEG_MAX_READS or frac <= HIFI_NEG_MAX_FRAC
        ok = ok and hifi_unit <= HIFI_STUTTER_RATIO * ill
        return {**out, "pass": ok, "reason": None if ok else "hifi_carries_allele"}

    counts = base_counts(bam, chrom, site0, min_mapq=HIFI_MIN_MAPQ)
    n = sum(counts.values())
    alt = counts.get(allele, 0)
    frac = alt / n if n else 0.0
    out.update(depth=n, alt_reads=alt, alt_frac=round(frac, 4))
    if n < HIFI_MIN_DEPTH:
        return {**out, "pass": False, "reason": "hifi_low_depth"}
    if pos:
        zyg = (row.get("label") or {}).get("zygosity", "het")
        floor = HIFI_ALT_FRAC.get(zyg, HIFI_ALT_FRAC["het"])
        ok = alt >= HIFI_MIN_ALT and frac >= floor
        return {**out, "pass": ok, "reason": None if ok else "hifi_alt_absent"}
    ok = alt <= HIFI_NEG_MAX_READS and frac <= HIFI_NEG_MAX_FRAC
    return {**out, "pass": ok, "reason": None if ok else "hifi_carries_allele"}


# ---------------------------------------------------------------------------
# Trio recurrence
# ---------------------------------------------------------------------------


def trio_recurrence(
    row: dict[str, Any],
    parents: dict[str, Any],
    parent_truth: dict[str, TruthSet] | None = None,
) -> dict[str, Any]:
    """The candidate allele's fraction in each parent's 60x Illumina BAM.

    ``recurrent`` is True when at least one *applicable* parent carries it at >= 5%: a
    parent is applicable where its own v4.2.1 truth is reference at the site, or its BED
    does not cover it (a parent that truly carries the variant says nothing about
    artifacts). Without ``parent_truth`` every parent counts and ``parent_truth`` is
    reported as ``"unchecked"``.
    """
    kind, allele = _row_allele(row)
    chrom, site0 = row["chrom"], int(row["site0"])
    s0, e0 = row["context"]["run"] if kind == "tr" else (site0, site0 + 1)
    out: dict[str, Any] = {"parents": {}}
    recurrent = False
    for name, bam in parents.items():
        applicable: bool | str = "unchecked"
        if parent_truth and name in parent_truth:
            ts = parent_truth[name]
            near = ts.records_near(chrom, s0, e0, 0)
            applicable = not near  # reference (covered or not) => applicable
        if kind == "tr":
            hist = spanning_hist(bam, chrom, s0, e0)
            n = sum(hist.values())
            frac = hist[int(allele)] / n if n else 0.0
        else:
            c = base_counts(bam, chrom, site0)
            n = sum(c.values())
            frac = c.get(allele, 0) / n if n else 0.0
        out["parents"][name] = {"depth": n, "frac": round(frac, 4), "applicable": applicable}
        if applicable is not False and frac >= TRIO_RECUR_FRAC:
            recurrent = True
    out["recurrent"] = recurrent
    out["parent_truth"] = "checked" if parent_truth else "unchecked"
    return out


# ---------------------------------------------------------------------------
# Clusters and partition
# ---------------------------------------------------------------------------


@dataclass
class ClusterIndex:
    """Connected components of the ``genomicSuperDups`` pairing.

    Nodes are merged segdup intervals per chromosome (overlapping segments are one node);
    every table row links a segment's node to its partner's. A site outside every segment
    falls in a 1 Mb bin (``bin:<chrom>:<Mb>``).
    """

    nodes: dict[str, Intervals] = field(default_factory=dict)
    comp: dict[tuple[str, int], str] = field(default_factory=dict)  # (chrom, node idx) -> id
    sizes: Counter[str] = field(default_factory=Counter)

    def cluster_id(self, chrom: str, pos0: int) -> str:
        iv = self.nodes.get(chrom)
        if iv is not None and iv.contains(pos0):
            i = bisect.bisect_right(iv.starts, pos0) - 1
            return self.comp[(chrom, i)]
        return f"bin:{chrom}:{pos0 // BIN_BP}"


def cluster_index(path: str | Path | None = None) -> ClusterIndex:
    """Build the cluster index from a ``genomicSuperDups.txt.gz`` (default: the UCSC hg38
    table, downloaded once into the cache)."""
    p = cached(path or SUPERDUPS_URL)
    pairs: list[tuple[str, int, int, str, int, int]] = []
    raw: dict[str, list[tuple[int, int]]] = {}
    opener: Any = gzip.open if p.suffix == ".gz" else open
    with opener(p, "rt") as fh:
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) < 10:
                continue
            c1, s1, e1 = f[1], int(f[2]), int(f[3])
            c2, s2, e2 = f[7], int(f[8]), int(f[9])
            pairs.append((c1, s1, e1, c2, s2, e2))
            raw.setdefault(c1, []).append((s1, e1))
            raw.setdefault(c2, []).append((s2, e2))
    nodes = {c: Intervals(v) for c, v in raw.items()}
    keys = [(c, i) for c, iv in nodes.items() for i in range(len(iv))]
    parent = {k: k for k in keys}

    def find(k: tuple[str, int]) -> tuple[str, int]:
        while parent[k] != k:
            parent[k] = parent[parent[k]]
            k = parent[k]
        return k

    def node(c: str, s: int) -> tuple[str, int]:
        return (c, bisect.bisect_right(nodes[c].starts, s) - 1)

    for c1, s1, _e1, c2, s2, _e2 in pairs:
        a, b = find(node(c1, s1)), find(node(c2, s2))
        if a != b:
            parent[max(a, b, key=_node_order)] = min(a, b, key=_node_order)
    idx = ClusterIndex(nodes=nodes)
    for k in keys:
        root = find(k)
        cid = f"sd:{root[0]}:{nodes[root[0]].starts[root[1]]}"
        idx.comp[k] = cid
        idx.sizes[cid] += 1
    return idx


def _chrom_order(c: str) -> tuple[int, str]:
    tail = c[3:] if c.startswith("chr") else c
    return (int(tail), "") if tail.isdigit() else (1000, tail)


def _node_order(k: tuple[str, int]) -> tuple[tuple[int, str], int]:
    return (_chrom_order(k[0]), k[1])


def partition(
    cluster: str, *, salt: str = PARTITION_SALT, pilot_fraction: float = PILOT_FRACTION
) -> str:
    """``"pilot"`` or ``"main"`` from a seeded hash of the cluster ID (never of the site)."""
    h = int.from_bytes(hashlib.sha256(f"{salt}:{cluster}".encode()).digest()[:8], "big")
    return "pilot" if h / 2**64 < pilot_fraction else "main"


# ---------------------------------------------------------------------------
# Contamination
# ---------------------------------------------------------------------------


def estimate_contamination(
    bam: Any,
    truth: TruthSet,
    chroms: Iterable[str],
    *,
    n_sites: int = 2000,
    seed: int = 1,
    exclude: BedSet | None = None,
) -> dict[str, Any]:
    """Cross-sample contamination from the reference fraction at hom-alt truth SNVs.

    At a site where HG002 is hom-alt, reference reads come from sequencing error or from a
    contaminant carrying the reference allele. With the per-base error to one specific base
    estimated from the third/fourth bases, the excess reference fraction is the
    contamination share times the contaminant's reference-allele frequency; assuming that
    frequency is about 0.5 at common SNPs, ``contamination ~ 2 x excess``. Sites in
    ``exclude`` (e.g. segdups) are skipped.
    """
    import random

    rng = random.Random(f"contamination:{seed}")  # noqa: S311
    pool = [
        r
        for c in chroms
        for r in truth.positives(c)
        if r.kind == "SNV"
        and r.zygosity == "hom_alt"
        and not (exclude and exclude.contains(r.chrom, r.start0))
    ]
    rng.shuffle(pool)
    ref_n = alt_n = other_n = 0
    used = 0
    for r in pool[:n_sites]:
        c = base_counts(bam, r.chrom, r.start0, min_mapq=20)
        d = sum(v for k, v in c.items() if k != "-")
        if d < 10:
            continue
        used += 1
        ref_n += c.get(r.ref.upper(), 0)
        alt_n += c.get(r.alts[0].upper(), 0)
        other_n += d - c.get(r.ref.upper(), 0) - c.get(r.alts[0].upper(), 0)
    total = ref_n + alt_n + other_n
    ref_frac = ref_n / total if total else 0.0
    per_base_err = other_n / total / 2 if total else 0.0
    excess = max(0.0, ref_frac - per_base_err)
    return {
        "sites": used,
        "ref_fraction": round(ref_frac, 5),
        "per_base_error": round(per_base_err, 5),
        "contamination_est": round(2 * excess, 5),
        "assumption": "contaminant reference-allele frequency ~0.5 at hom-alt common SNPs",
    }


__all__ = [
    "AUTOSOMES",
    "BedSet",
    "ClusterIndex",
    "HG002_BAM",
    "TruthRecord",
    "TruthSet",
    "base_counts",
    "cache_dir",
    "cached",
    "cluster_index",
    "concordant",
    "estimate_contamination",
    "is_tr_run",
    "left_normalize",
    "load_bed_multi",
    "open_bam",
    "open_vcf",
    "orthogonal_check",
    "partition",
    "reference_concordant",
    "repeat_run",
    "spanning_hist",
    "trio_recurrence",
    "unit_fraction",
]
