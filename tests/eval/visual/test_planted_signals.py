"""The default-seed vt_eval BAM carries each case type's intended evidence.

Ground truth for the visual-vs-text eval is only meaningful if every planted site really
shows its pattern; this checks the pattern from the BAM itself (not from the plan).
"""

from __future__ import annotations

from collections import defaultdict
from pathlib import Path

import pysam
import pytest

from bamcp.eval.visual.fixtures import build_plan, create_visual_eval_fixtures
from bamcp.eval.visual.stimuli import _read_observations

SEED = 1


@pytest.fixture(scope="module")
def observed(tmp_path_factory: pytest.TempPathFactory):
    out = tmp_path_factory.mktemp("vt_full")
    paths = create_visual_eval_fixtures(SEED, out)
    plan = build_plan(SEED)
    result = {}
    with pysam.AlignmentFile(str(paths["bam"])) as bam:
        for item in plan.items:
            ref = plan.reference(item.contig)
            per_site = defaultdict(list)
            for read in bam.fetch(item.contig, min(item.sites), max(item.sites) + 1):
                obs = _read_observations(read, set(item.sites))
                for s, o in obs.items():
                    per_site[s].append(o)
            result[item.item_id] = (item, ref, per_site)
    return result


def _alt(item, ref, per_site, site=None):
    s = item.site if site is None else site
    obs = per_site[s]
    alt = [o for o in obs if o.allele != ref[s]]
    return obs, alt


def _each(observed, case):
    items = [v for v in observed.values() if v[0].case_type == case]
    assert len(items) == 3
    return items


def test_real_snvs(observed) -> None:
    for item, ref, ps in _each(observed, "het_snv"):
        obs, alt = _alt(item, ref, ps)
        assert 0.3 <= len(alt) / len(obs) <= 0.7
        assert any(o.reverse for o in alt) and any(not o.reverse for o in alt)
    for item, ref, ps in _each(observed, "hom_alt_snv"):
        obs, alt = _alt(item, ref, ps)
        assert len(alt) / len(obs) >= 0.9
    for item, ref, ps in _each(observed, "low_vaf_snv"):
        obs, alt = _alt(item, ref, ps)
        assert 0.05 <= len(alt) / len(obs) <= 0.2 and len(alt) >= 3
    for item, ref, ps in _each(observed, "low_depth_het"):
        obs, alt = _alt(item, ref, ps)
        assert len(obs) <= 4 and len(alt) >= 1


def test_artifacts(observed) -> None:
    for item, ref, ps in _each(observed, "strand_bias"):
        _, alt = _alt(item, ref, ps)
        assert len(alt) >= 4 and all(not o.reverse for o in alt if o.allele == item.params["alt"])
    for item, ref, ps in _each(observed, "low_baseq_alt"):
        _, alt = _alt(item, ref, ps)
        planted = [o for o in alt if o.allele == item.params["alt"]]
        assert len(planted) >= 5 and sum(o.bq < 15 for o in planted) >= len(planted) - 1
    for item, ref, ps in _each(observed, "read_end_bias"):
        _, alt = _alt(item, ref, ps)
        planted = [o for o in alt if o.allele == item.params["alt"]]
        assert len(planted) >= 7 and all(o.dist3 <= 4 for o in planted)
    for item, ref, ps in _each(observed, "low_mapq_paralog"):
        obs, alt = _alt(item, ref, ps)
        planted = [o for o in alt if o.allele == item.params["alt"]]
        assert len(planted) >= 10 and all(o.mapq <= 3 for o in planted)
        assert all(o.mapq == 60 for o in obs if o.allele == ref[item.site])
    for item, ref, ps in _each(observed, "pcr_duplicate"):
        _, alt = _alt(item, ref, ps)
        planted = [o for o in alt if o.allele == item.params["alt"]]
        from collections import Counter

        assert Counter(o.span for o in planted).most_common(1)[0][1] >= 7
    for item, ref, ps in _each(observed, "homopolymer_del"):
        obs, alt = _alt(item, ref, ps)
        dels = [o for o in obs if o.allele == "-"]
        assert 0.08 <= len(dels) / len(obs) <= 0.25


def test_indels_and_multisite(observed) -> None:
    for item, ref, ps in _each(observed, "real_del"):
        obs, _ = _alt(item, ref, ps)
        assert 0.3 <= sum(o.allele == "-" for o in obs) / len(obs) <= 0.7
    for item, ref, ps in _each(observed, "real_ins"):
        obs, _ = _alt(item, ref, ps)
        assert 0.3 <= sum("+" in o.allele for o in obs) / len(obs) <= 0.7
    for item, ref, ps in _each(observed, "multi_allelic"):
        obs, alt = _alt(item, ref, ps)
        assert len({o.allele for o in alt if o.allele in "ACGT"}) >= 2
    for case in ("phased_cis", "phased_trans"):
        for item, ref, ps in _each(observed, case):
            s1, s2 = item.sites
            assert 25 <= s2 - s1 <= 45
            _, a1 = _alt(item, ref, ps, s1)
            _, a2 = _alt(item, ref, ps, s2)
            assert len(a1) >= 5 and len(a2) >= 5
    for item, ref, ps in _each(observed, "clean_negative"):
        obs, alt = _alt(item, ref, ps)
        assert len(alt) <= 2


def test_negatives_have_no_planted_variant() -> None:
    plan = build_plan(SEED)
    for it in plan.items:
        if it.label["variant_present"] == "no":
            assert not any(e["snv"] or e["del"] or e["ins"] for e in it.hap_edits)


def test_fixture_dir_untouched(tmp_path: Path) -> None:
    # Generating into a custom dir never writes into tests/fixtures.
    create_visual_eval_fixtures(SEED, tmp_path, instances=1)
    assert (tmp_path / "vt_eval.bam").exists()
