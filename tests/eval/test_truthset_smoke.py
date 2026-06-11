"""Ground-truth regression gate.

Runs the deterministic truth-set scorer over the synthetic_v1 dataset and
asserts the metric floors. No LLM, no network — this is what makes the eval
harness load-bearing in CI: a real regression in variant detection or artifact
scoring fails the build here.

The conftest session fixtures generate the comprehensive BAM/reference on first
use, so this reproduces on a fresh checkout.
"""

from __future__ import annotations

from pathlib import Path

import pytest

from bamcp.config import BAMCPConfig
from bamcp.eval.truthset import load_truthset, score_truthset

MANIFEST = Path("tests/eval/datasets/synthetic_v1/manifest.yaml")


@pytest.fixture(scope="module")
def _ensure_fixtures(comprehensive_bam_path, comprehensive_ref_fasta_path):
    """Trigger fixture generation (paths in the manifest are repo-relative)."""
    return comprehensive_bam_path, comprehensive_ref_fasta_path


@pytest.fixture(scope="module")
def report(_ensure_fixtures):
    import asyncio

    truthset = load_truthset(MANIFEST)
    config = BAMCPConfig(reference=truthset.reference)
    return asyncio.run(score_truthset(truthset, config))


@pytest.mark.integration
def test_truthset_meets_floors(report):
    passed, failures = report.meets_floors(
        min_recall=1.0,
        min_precision=0.95,
        min_artifact_recall=1.0,
        max_false_positives=0,
    )
    assert passed, "Truth-set floors not met:\n" + "\n".join(failures)


@pytest.mark.integration
def test_variant_detection_perfect_on_synthetic(report):
    # The synthetic fixture is fully controlled — detection should be exact.
    assert report.detection.recall == 1.0
    assert report.detection.precision == 1.0
    assert report.detection.tp == 8
    assert not report.missing_calls
    assert not report.spurious_calls


@pytest.mark.integration
def test_no_false_positives_in_negative_regions(report):
    assert report.false_positives == 0


@pytest.mark.integration
def test_artifact_types_surfaced(report):
    # Every known-artifact site surfaces its expected risk type.
    assert report.artifact_recall == 1.0
    found = {c.site: c.expected_type for c in report.artifact_checks if c.found}
    assert found == {
        "chr1:1050": "strand_bias",
        "chr1:1200": "low_mapq",
        "chr1:1802": "homopolymer",
    }


@pytest.mark.integration
def test_clean_control_discriminated(report):
    assert report.clean_discriminated is True


@pytest.mark.integration
def test_truthset_cli_exits_zero(_ensure_fixtures, capsys):
    # Exercises the full CLI path against real fixtures.
    from bamcp.eval.truthset import main

    rc = main(["--manifest", str(MANIFEST)])
    out = capsys.readouterr().out
    assert rc == 0
    assert "PASS" in out
    assert "variant detection" in out


@pytest.mark.integration
def test_truthset_cli_json_output(_ensure_fixtures, capsys):
    import json

    from bamcp.eval.truthset import main

    rc = main(["--manifest", str(MANIFEST), "--json"])
    payload = json.loads(capsys.readouterr().out)
    assert rc == 0
    assert payload["passed"] is True
    assert payload["detection"]["recall"] == 1.0
