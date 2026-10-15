"""Deterministic vt_eval fixture generator."""

from __future__ import annotations

import hashlib
from pathlib import Path

import pysam

from bamcp.eval.visual.cases import CASE_TYPES
from bamcp.eval.visual.fixtures import (
    BAM_NAME,
    REF_NAME,
    build_plan,
    create_visual_eval_fixtures,
    ensure_fixtures,
)


def _digest(paths: dict[str, Path]) -> dict[str, str]:
    return {k: hashlib.sha256(p.read_bytes()).hexdigest() for k, p in paths.items()}


def test_same_seed_identical_bytes(tmp_path: Path) -> None:
    a = create_visual_eval_fixtures(3, tmp_path / "a", instances=1)
    b = create_visual_eval_fixtures(3, tmp_path / "b", instances=1)
    assert set(a) == {"bam", "bai", "reference", "fai"}
    assert _digest(a) == _digest(b)


def test_different_seed_differs(tmp_path: Path) -> None:
    a = create_visual_eval_fixtures(3, tmp_path / "a", instances=1)
    b = create_visual_eval_fixtures(4, tmp_path / "b", instances=1)
    assert _digest(a)["bam"] != _digest(b)["bam"]
    assert _digest(a)["reference"] != _digest(b)["reference"]


def test_plan_is_pure_and_covers_every_case() -> None:
    p1, p2 = build_plan(5, 2), build_plan(5, 2)
    assert [(i.item_id, i.contig, i.sites, i.label) for i in p1.items] == [
        (i.item_id, i.contig, i.sites, i.label) for i in p2.items
    ]
    assert len(p1.items) == 2 * len(CASE_TYPES)
    assert {i.case_type for i in p1.items} == {c.name for c in CASE_TYPES}
    assert len({i.item_id for i in p1.items}) == len(p1.items)
    # Opaque ids carry no case-type text.
    assert all(i.case_type not in i.item_id for i in p1.items)


def test_reference_is_not_periodic_and_has_planted_homopolymer() -> None:
    plan = build_plan(1, 1)
    seq = plan.reference("chr1")
    assert seq[:10] * 5 != seq[:50]
    hp = next(i for i in plan.items if i.case_type == "homopolymer_del")
    ref = plan.reference(hp.contig)
    assert ref[hp.site : hp.site + 9] == "A" * 9
    assert ref[hp.site - 1] != "A" and ref[hp.site + 9] != "A"


def test_bam_properties(tmp_path: Path) -> None:
    paths = create_visual_eval_fixtures(2, tmp_path, instances=1)
    with pysam.AlignmentFile(str(paths["bam"])) as bam:
        reads = list(bam.fetch(until_eof=True))
    lengths = {r.query_alignment_length for r in reads}
    assert min(lengths) >= 95 and max(r.query_length for r in reads) <= 170
    assert any(r.is_reverse for r in reads) and any(not r.is_reverse for r in reads)
    names = [r.query_name for r in reads]
    assert len(set(names)) == len(names)
    assert all(len(n) == 16 and all(c in "0123456789abcdef" for c in n) for n in names)
    quals = [q for r in reads[:200] for q in r.query_qualities]
    assert len(set(quals)) > 10  # realistic spread, not a flat Q40
    with pysam.FastaFile(str(paths["reference"])) as fa:
        assert list(fa.references) == ["chr1", "chr2", "chr3"]


def test_ensure_fixtures_reuses_matching_stamp(tmp_path: Path) -> None:
    first = ensure_fixtures(9, tmp_path, instances=1)
    mtime = first["bam"].stat().st_mtime_ns
    again = ensure_fixtures(9, tmp_path, instances=1)
    assert again["bam"].stat().st_mtime_ns == mtime
    ensure_fixtures(10, tmp_path, instances=1)  # seed change → regenerated
    assert (tmp_path / BAM_NAME).exists() and (tmp_path / REF_NAME).exists()
    assert (tmp_path / "vt_eval.stamp").read_text().startswith("seed=10")


def test_script_wrapper(tmp_path: Path) -> None:
    import importlib.util

    script = Path(__file__).resolve().parents[3] / "create_visual_eval_fixtures.py"
    spec = importlib.util.spec_from_file_location("cvef", script)
    assert spec and spec.loader
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    assert mod.main(["--seed", "2", "--instances", "1", "--out-dir", str(tmp_path)]) == 0
    assert (tmp_path / BAM_NAME).exists()
