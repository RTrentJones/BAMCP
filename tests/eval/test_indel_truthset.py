"""Indel-detection regression gate.

Runs the deterministic truth-set scorer over the indel_v1 dataset — a 3bp
deletion and a 2bp insertion planted in the comprehensive fixture (Group I) —
and asserts perfect recovery. Complements the SNV-only synthetic_v1 gate so the
pileup-based indel caller (detect_indels) is load-bearing in CI: a regression in
indel detection or deletion-inclusive depth fails the build here. No LLM, no
network.
"""

from __future__ import annotations

import asyncio
from pathlib import Path

import pytest

from bamcp.config import BAMCPConfig
from bamcp.eval.truthset import load_truthset, score_truthset

MANIFEST = Path("tests/eval/datasets/indel_v1/manifest.yaml")


@pytest.fixture(scope="module")
def _ensure_fixtures(comprehensive_bam_path, comprehensive_ref_fasta_path):
    """Trigger fixture generation (paths in the manifest are repo-relative)."""
    return comprehensive_bam_path, comprehensive_ref_fasta_path


@pytest.fixture(scope="module")
def report(_ensure_fixtures):
    truthset = load_truthset(MANIFEST)
    config = BAMCPConfig(reference=truthset.reference)
    return asyncio.run(score_truthset(truthset, config))


@pytest.mark.integration
def test_indel_truthset_meets_floors(report):
    passed, failures = report.meets_floors(
        min_recall=1.0,
        min_precision=1.0,
        max_false_positives=0,
    )
    assert passed, "Indel truth-set floors not met:\n" + "\n".join(failures)


@pytest.mark.integration
def test_both_indels_detected(report):
    # The planted deletion and insertion are recovered exactly, with nothing
    # spurious called in the region.
    assert report.detection.tp == 2
    assert report.detection.recall == 1.0
    assert report.detection.precision == 1.0
    assert not report.missing_calls
    assert not report.spurious_calls


@pytest.mark.integration
def test_no_false_positive_indels_in_negative_regions(report):
    assert report.false_positives == 0


@pytest.mark.integration
def test_indel_cli_exits_zero(_ensure_fixtures, capsys):
    from bamcp.eval.truthset import main

    rc = main(["--manifest", str(MANIFEST), "--min-precision", "1.0"])
    out = capsys.readouterr().out
    assert rc == 0
    assert "PASS" in out
    assert "tp=2" in out
