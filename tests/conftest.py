"""Shared test fixtures for BAMCP tests."""

import os
import pytest

FIXTURES_DIR = os.path.join(os.path.dirname(__file__), "fixtures")


@pytest.fixture
def fixtures_dir():
    """Path to the test fixtures directory."""
    return FIXTURES_DIR


@pytest.fixture
def small_bam_path():
    """Path to the small test BAM file."""
    return os.path.join(FIXTURES_DIR, "small.bam")


@pytest.fixture
def empty_bam_path():
    """Path to an empty BAM file (header only)."""
    return os.path.join(FIXTURES_DIR, "empty.bam")


@pytest.fixture
def ref_fasta_path():
    """Path to the reference FASTA file."""
    return os.path.join(FIXTURES_DIR, "ref.fa")
