"""Deterministic synthetic BAM for the visual-vs-text eval (``vt_eval.bam``).

Everything is derived from one integer seed: a non-periodic reference (several contigs),
one planted site per item at a seeded random slot, and ~30x of 100-150 bp single-end reads
with jittered starts, a realistic base-quality profile, ~0.3% random substitution errors,
hashed read names and both orientations. The same seed always produces byte-identical
BAM/BAI/FASTA/FAI files (``tests/unit/eval/visual`` asserts it).

The item plan (:func:`build_plan`) is pure and cheap, so the stimulus builder and the
labels writer re-derive it from the seed instead of persisting site metadata next to
the fixture — nothing label-bearing ever sits beside the BAM.

Simplifications (recorded in ``docs/eval/visual-preregistration.md``): single-end reads
only; no soft-clip breakpoint or large-deletion-dip cases; no misalignment cluster.
"""

from __future__ import annotations

import hashlib
import os
import random
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

from .cases import CASE_TYPES, CaseSpec

BAM_NAME = "vt_eval.bam"
REF_NAME = "vt_eval_ref.fa"
STAMP_NAME = "vt_eval.stamp"

N_CONTIGS = 3
SLOT_BP = 3000
CONTIG_PAD = 1000
HALF_WINDOW = 350  # reads are simulated within site +/- this many bp
READ_LEN_MIN = 100
READ_LEN_MAX = 150
BASE_ERROR_RATE = 0.003
DEFAULT_DEPTH = 30
_BASES = "ACGT"


def _derive_rng(seed: int, *parts: object) -> random.Random:
    # str seeds hash through sha512 inside random.Random — stable across processes,
    # unlike hash(), which is salted per interpreter.
    return random.Random(":".join(str(p) for p in (seed, *parts)))  # noqa: S311


def opaque_id(*parts: object, n: int = 8) -> str:
    return hashlib.sha256(":".join(str(p) for p in parts).encode()).hexdigest()[:n]


@dataclass
class Item:
    """One planted site. ``sites`` are 0-based; labels are 1-based."""

    item_id: str
    case_type: str
    instance: int
    contig: str
    sites: list[int]
    depth: int
    hap_edits: list[dict[str, dict[int, Any]]] = field(default_factory=list)
    params: dict[str, Any] = field(default_factory=dict)
    label: dict[str, Any] = field(default_factory=dict)

    @property
    def site(self) -> int:
        return self.sites[0]

    @property
    def center(self) -> int:
        return (self.sites[0] + self.sites[-1]) // 2

    def positions_1based(self) -> list[int]:
        return [s + 1 for s in self.sites]


@dataclass
class Plan:
    seed: int
    contigs: dict[str, list[str]]  # name -> mutable base list
    items: list[Item]

    def reference(self, contig: str) -> str:
        return "".join(self.contigs[contig])


# ---------------------------------------------------------------------------
# Reference
# ---------------------------------------------------------------------------


def _random_reference(rng: random.Random, length: int) -> list[str]:
    """GC ~41%, runs capped at 3 so only planted homopolymers exist."""
    weights = [0.295, 0.205, 0.205, 0.295]
    seq: list[str] = []
    while len(seq) < length:
        b = rng.choices(_BASES, weights)[0]
        if len(seq) >= 3 and seq[-1] == seq[-2] == seq[-3] == b:
            continue
        seq.append(b)
    return seq


def _other_base(rng: random.Random, ref: str, exclude: tuple[str, ...] = ()) -> str:
    choices = [b for b in _BASES if b != ref and b not in exclude]
    return rng.choice(choices)


def _ensure_unique_neighbors(seq: list[str], pos: int, rng: random.Random) -> None:
    """Make seq[pos] differ from both neighbours (a deletion there is unambiguous)."""
    if seq[pos] in (seq[pos - 1], seq[pos + 1]):
        seq[pos] = _other_base(rng, seq[pos - 1], exclude=(seq[pos + 1],))


# ---------------------------------------------------------------------------
# Plan
# ---------------------------------------------------------------------------


def build_plan(seed: int, instances: int = 3) -> Plan:
    """Deterministic reference + item plan for ``seed``."""
    rng = _derive_rng(seed, "plan")
    n_items = len(CASE_TYPES) * instances
    slots_per_contig = -(-n_items // N_CONTIGS)
    contig_len = 2 * CONTIG_PAD + slots_per_contig * SLOT_BP
    names = [f"chr{i + 1}" for i in range(N_CONTIGS)]
    contigs = {
        name: _random_reference(_derive_rng(seed, "ref", name), contig_len) for name in names
    }

    slots = [(c, k) for c in names for k in range(slots_per_contig)]
    rng.shuffle(slots)

    items: list[Item] = []
    idx = 0
    for spec in CASE_TYPES:
        for inst in range(instances):
            contig, k = slots[idx]
            idx += 1
            irng = _derive_rng(seed, "item", spec.name, inst)
            site = CONTIG_PAD + k * SLOT_BP + 600 + irng.randint(0, SLOT_BP - 1200)
            item = Item(
                item_id="i" + opaque_id(seed, "item", spec.name, inst),
                case_type=spec.name,
                instance=inst,
                contig=contig,
                sites=[site],
                depth=spec.depth,
            )
            _plant(item, spec, contigs[contig], irng)
            items.append(item)
    return Plan(seed=seed, contigs=contigs, items=items)


def _plant(item: Item, spec: CaseSpec, seq: list[str], rng: random.Random) -> None:
    """Edit the reference around the site and fill haplotype edits + labels."""
    s = item.site
    h0: dict[str, dict[int, Any]] = {"snv": {}, "del": {}, "ins": {}}
    h1: dict[str, dict[int, Any]] = {"snv": {}, "del": {}, "ins": {}}
    kind = spec.name
    alt_label = ""

    if kind == "homopolymer_del":
        flank = _other_base(rng, "A")
        seq[s - 1] = flank
        for i in range(9):
            seq[s + i] = "A"
        seq[s + 9] = _other_base(rng, "A")
        item.params["frac"] = rng.uniform(0.12, 0.18)
    elif kind == "real_del":
        _ensure_unique_neighbors(seq, s, rng)
        h1["del"][s] = 1
        alt_label = seq[s]
    elif kind == "real_ins":
        ins = "".join(rng.choice(_BASES) for _ in range(rng.randint(2, 3)))
        # Make the insertion unambiguous: it must not extend either flank.
        while ins[0] == seq[s + 1] or ins[-1] == seq[s] or len(set(ins)) == 1:
            ins = "".join(rng.choice(_BASES) for _ in range(len(ins)))
        h1["ins"][s] = ins
        alt_label = ins
    elif kind == "mnv":
        a1 = _other_base(rng, seq[s])
        a2 = _other_base(rng, seq[s + 1])
        h1["snv"][s] = a1
        h1["snv"][s + 1] = a2
        alt_label = a1 + a2
    elif kind == "multi_allelic":
        a1 = _other_base(rng, seq[s])
        a2 = _other_base(rng, seq[s], exclude=(a1,))
        h1["snv"][s] = a1
        h0["snv"][s] = a2
        alt_label = ",".join(sorted([a1, a2]))
    elif kind in ("phased_cis", "phased_trans"):
        s2 = s + rng.randint(25, 45)
        item.sites = [s, s2]
        a1 = _other_base(rng, seq[s])
        a2 = _other_base(rng, seq[s2])
        h1["snv"][s] = a1
        (h1 if kind == "phased_cis" else h0)["snv"][s2] = a2
        alt_label = ",".join(sorted([a1, a2]))
    elif kind == "hom_alt_snv":
        a = _other_base(rng, seq[s])
        h0["snv"][s] = a
        h1["snv"][s] = a
        alt_label = a
    elif kind in ("het_snv", "low_vaf_snv", "low_depth_het"):
        a = _other_base(rng, seq[s])
        h1["snv"][s] = a
        alt_label = a
    elif kind in (
        "low_baseq_alt",
        "strand_bias",
        "read_end_bias",
        "low_mapq_paralog",
        "pcr_duplicate",
    ):
        item.params["alt"] = _other_base(rng, seq[s])
        if kind == "low_mapq_paralog":
            offs = rng.sample([o for o in range(-40, 41) if abs(o) > 4], 3)
            item.params["paralog"] = {s + o: _other_base(rng, seq[s + o]) for o in offs}
    # clean_negative / scattered_errors: nothing planted.

    if kind == "het_snv":
        item.params["frac"] = rng.uniform(0.42, 0.58)
    elif kind == "low_vaf_snv":
        item.params["frac"] = rng.uniform(0.10, 0.14)

    item.hap_edits = [h0, h1]
    item.label = spec.label(alt_label)


# ---------------------------------------------------------------------------
# Read simulation
# ---------------------------------------------------------------------------


@dataclass
class SimRead:
    start: int
    length: int
    reverse: bool
    hap: int
    snv: dict[int, str] = field(default_factory=dict)
    dels: dict[int, int] = field(default_factory=dict)
    ins: dict[int, str] = field(default_factory=dict)
    mapq: int = 60
    bq_override: dict[int, int] = field(default_factory=dict)
    clip3: str = ""
    error_rate: float = BASE_ERROR_RATE
    error_window: tuple[int, int] | None = None
    copy_of: SimRead | None = None

    def covers(self, pos: int) -> bool:
        # Approximate aligned span (indels shift by a base or two; callers pad).
        return self.start + 2 <= pos <= self.start + self.length - 3


def _hap_reads(item: Item, rng: random.Random) -> list[SimRead]:
    lo = max(0, item.center - HALF_WINDOW)
    hi = item.center + HALF_WINDOW
    mean_len = (READ_LEN_MIN + READ_LEN_MAX) / 2
    n = round(item.depth * (hi - lo + mean_len) / mean_len)
    reads = []
    for _ in range(n):
        length = rng.randint(READ_LEN_MIN, READ_LEN_MAX)
        start = rng.randint(max(0, lo - length + 20), hi - 20)
        reads.append(SimRead(start, length, rng.random() < 0.5, int(rng.random() < 0.5)))
    return reads


def _choose(reads: list[SimRead], k: int, rng: random.Random) -> list[SimRead]:
    """Pick ``k`` reads spread evenly over both orientations."""
    fwd = [r for r in reads if not r.reverse]
    rev = [r for r in reads if r.reverse]
    k_fwd = min(len(fwd), round(k * len(fwd) / max(1, len(reads))))
    k_rev = min(len(rev), k - k_fwd)
    return rng.sample(fwd, k_fwd) + rng.sample(rev, k_rev)


def simulate_item(item: Item, ref: str, seed: int) -> list[SimRead]:
    rng = _derive_rng(seed, "reads", item.item_id)
    s = item.site
    kind = item.case_type
    reads = _hap_reads(item, rng)

    if kind == "low_depth_het":
        reads = [r for r in reads if not (r.start - 5 <= s <= r.start + r.length + 5)]
        for i in range(3):
            length = rng.randint(READ_LEN_MIN, READ_LEN_MAX)
            start = s - rng.randint(15, length - 15)
            reads.append(SimRead(start, length, i == 1, int(i != 2)))

    covering = [r for r in reads if r.covers(s)]
    if "frac" in item.params and kind != "homopolymer_del":
        k = max(3, round(item.params["frac"] * len(covering)))
        chosen = {id(r) for r in _choose(covering, k, rng)}
        for r in covering:
            r.hap = int(id(r) in chosen)

    for r in reads:
        edits = item.hap_edits[r.hap]
        r.snv.update(edits["snv"])
        r.dels.update(edits["del"])
        r.ins.update(edits["ins"])

    alt: str = item.params.get("alt", "")
    if kind == "homopolymer_del":
        spanning = [r for r in reads if r.start + 2 <= s and s + 12 <= r.start + r.length - 3]
        k = max(3, round(item.params["frac"] * len(spanning)))
        for r in _choose(spanning, k, rng):
            r.dels[s] = 1
    elif kind == "low_baseq_alt":
        for r in _choose(covering, round(rng.uniform(0.25, 0.35) * len(covering)), rng):
            r.snv[s] = alt
            r.bq_override[s] = rng.randint(4, 12)
    elif kind == "strand_bias":
        fwd = [r for r in covering if not r.reverse]
        for r in rng.sample(fwd, round(rng.uniform(0.35, 0.5) * len(fwd))):
            r.snv[s] = alt
    elif kind == "read_end_bias":
        for i in range(rng.randint(8, 10)):
            length = rng.randint(READ_LEN_MIN, READ_LEN_MAX)
            d = rng.randint(0, 4)
            reverse = i % 2 == 1
            start = s - d if reverse else s + d - length + 1
            clip = ""
            if i % 3 != 2:
                clip = "".join(rng.choice(_BASES) for _ in range(rng.randint(6, 15)))
            reads.append(SimRead(start, length, reverse, 0, snv={s: alt}, clip3=clip, mapq=60))
    elif kind == "low_mapq_paralog":
        for i in range(rng.randint(10, 12)):
            length = rng.randint(READ_LEN_MIN, READ_LEN_MAX)
            start = s - rng.randint(20, length - 20)
            snv = {s: alt, **item.params["paralog"]}
            reads.append(SimRead(start, length, i % 2 == 1, 0, snv=snv, mapq=rng.randint(0, 3)))
    elif kind == "pcr_duplicate":
        length = rng.randint(READ_LEN_MIN + 10, READ_LEN_MAX)
        start = s - rng.randint(30, length - 30)
        reverse = rng.random() < 0.5
        proto = SimRead(start, length, reverse, 0, snv={s: alt})
        reads.append(proto)
        for _ in range(rng.randint(6, 8)):
            reads.append(SimRead(start, length, reverse, 0, snv={s: alt}, copy_of=proto))
    elif kind == "scattered_errors":
        for r in reads:
            r.error_rate = 0.03
            r.error_window = (s - 30, s + 30)
    return reads


def _quality_profile(rng: random.Random, n: int) -> list[int]:
    out = []
    for i in range(n):
        q = 37 - 12 * (i / n) ** 3 + rng.gauss(0, 2.5)
        if rng.random() < 0.02:
            q = rng.uniform(10, 20)
        out.append(int(min(41, max(2, round(q)))))
    return out


def _render_read(
    r: SimRead, ref: str, rng: random.Random
) -> tuple[int, str, list[int], list[tuple[int, int]]]:
    """Walk the reference from ``r.start`` → (aligned start, sequence, qualities, cigar)."""
    seq: list[str] = []
    is_ref: list[int | None] = []  # ref pos per query base (None for inserted)
    cigar: list[tuple[int, int]] = []

    def push(op: int, n: int) -> None:
        if cigar and cigar[-1][0] == op:
            cigar[-1] = (op, cigar[-1][1] + n)
        else:
            cigar.append((op, n))

    p = r.start
    while p in r.dels:  # never start on a deletion
        p += 1
    start = p
    while len(seq) < r.length:
        if p in r.dels and seq:
            push(2, r.dels[p])
            p += r.dels[p]
            continue
        seq.append(r.snv.get(p, ref[p]))
        is_ref.append(p)
        push(0, 1)
        if p in r.ins and len(seq) + len(r.ins[p]) < r.length:
            seq.extend(r.ins[p])
            is_ref.extend([None] * len(r.ins[p]))
            push(1, len(r.ins[p]))
        p += 1

    n = len(seq)
    quals = _quality_profile(rng, n)
    if r.reverse:
        quals.reverse()  # sequencing runs right-to-left on the reverse strand
    for i, rp in enumerate(is_ref):
        rate = r.error_rate
        if r.error_window and rp is not None and not r.error_window[0] <= rp <= r.error_window[1]:
            rate = BASE_ERROR_RATE
        if rng.random() < rate:
            seq[i] = _other_base(rng, seq[i])
            quals[i] = rng.randint(5, 20)
        if rp is not None and rp in r.bq_override:
            quals[i] = r.bq_override[rp]

    if r.clip3:
        clip_q = [rng.randint(8, 25) for _ in r.clip3]
        if r.reverse:
            seq = list(r.clip3) + seq
            quals = clip_q + quals
            cigar.insert(0, (4, len(r.clip3)))
        else:
            seq.extend(r.clip3)
            quals.extend(clip_q)
            cigar.append((4, len(r.clip3)))
    return start, "".join(seq), quals, cigar


def _read_name(seed: int, item_id: str, idx: int) -> str:
    return hashlib.sha1(f"{seed}:{item_id}:{idx}".encode()).hexdigest()[:16]  # noqa: S324


# ---------------------------------------------------------------------------
# Writing
# ---------------------------------------------------------------------------


def _write_fasta(plan: Plan, path: Path) -> None:
    import pysam

    with path.open("w") as f:
        for name in plan.contigs:
            seq = plan.reference(name)
            f.write(f">{name}\n")
            for i in range(0, len(seq), 80):
                f.write(seq[i : i + 80] + "\n")
    fai = path.with_name(path.name + ".fai")
    if fai.exists():
        fai.unlink()
    pysam.faidx(str(path))


def _write_bam(plan: Plan, path: Path) -> None:
    import pysam

    names = list(plan.contigs)
    header = {
        "HD": {"VN": "1.6", "SO": "coordinate"},
        "SQ": [{"SN": n, "LN": len(plan.contigs[n])} for n in names],
    }
    records: list[tuple[int, int, str, int, str, list[int], list[tuple[int, int]], int]] = []
    for item in plan.items:
        ref = plan.reference(item.contig)
        rng = _derive_rng(plan.seed, "render", item.item_id)
        rendered: dict[int, tuple[int, str, list[int], list[tuple[int, int]]]] = {}
        for idx, r in enumerate(simulate_item(item, ref, plan.seed)):
            if r.copy_of is not None and id(r.copy_of) in rendered:
                start, seq, quals, cigar = rendered[id(r.copy_of)]
                # Duplicates share bases (one molecule); quality calls differ per cluster.
                quals = [max(2, min(41, q + rng.randint(-2, 2))) for q in quals]
            else:
                start, seq, quals, cigar = _render_read(r, ref, rng)
                rendered[id(r)] = (start, seq, quals, cigar)
            records.append(
                (
                    names.index(item.contig),
                    start,
                    _read_name(plan.seed, item.item_id, idx),
                    16 if r.reverse else 0,
                    seq,
                    quals,
                    cigar,
                    r.mapq,
                )
            )
    records.sort(key=lambda t: (t[0], t[1], t[2]))
    with pysam.AlignmentFile(str(path), "wb", header=header) as out:
        for tid, start, name, flag, seq, quals, cigar, mapq in records:
            a = pysam.AlignedSegment(out.header)
            a.query_name = name
            a.query_sequence = seq
            a.flag = flag
            a.reference_id = tid
            a.reference_start = start
            a.mapping_quality = mapq
            a.cigartuples = cigar
            a.query_qualities = pysam.qualitystring_to_array("".join(chr(q + 33) for q in quals))
            out.write(a)
    bai = path.with_name(path.name + ".bai")
    if bai.exists():
        bai.unlink()
    pysam.index(str(path))


def create_visual_eval_fixtures(
    seed: int = 1, out_dir: str | os.PathLike[str] | None = None, instances: int = 3
) -> dict[str, Path]:
    """Write ``vt_eval.bam`` (+ ``.bai``) and ``vt_eval_ref.fa`` (+ ``.fai``).

    Returns the written paths. Labels are NOT written here (see :mod:`.stimuli`).
    """
    out = Path(out_dir) if out_dir is not None else default_fixture_dir()
    out.mkdir(parents=True, exist_ok=True)
    plan = build_plan(seed, instances)
    ref_path = out / REF_NAME
    bam_path = out / BAM_NAME
    _write_fasta(plan, ref_path)
    _write_bam(plan, bam_path)
    (out / STAMP_NAME).write_text(f"seed={seed} instances={instances}\n")
    return {
        "bam": bam_path,
        "bai": bam_path.with_name(BAM_NAME + ".bai"),
        "reference": ref_path,
        "fai": ref_path.with_name(REF_NAME + ".fai"),
    }


def default_fixture_dir() -> Path:
    # src/bamcp/eval/visual/fixtures.py → repo root/tests/fixtures
    return Path(__file__).resolve().parents[4] / "tests" / "fixtures"


def ensure_fixtures(
    seed: int, out_dir: str | os.PathLike[str] | None = None, instances: int = 3
) -> dict[str, Path]:
    """Create the fixtures unless a matching (seed, instances) set is already present."""
    out = Path(out_dir) if out_dir is not None else default_fixture_dir()
    stamp = out / STAMP_NAME
    paths = {
        "bam": out / BAM_NAME,
        "bai": out / (BAM_NAME + ".bai"),
        "reference": out / REF_NAME,
        "fai": out / (REF_NAME + ".fai"),
    }
    expected = f"seed={seed} instances={instances}\n"
    if all(p.exists() for p in paths.values()) and stamp.exists() and stamp.read_text() == expected:
        return paths
    return create_visual_eval_fixtures(seed, out, instances)
