"""Shared test fixtures for BAMCP tests."""

import os

import pytest

from bamcp.core import tools as _tools_module

FIXTURES_DIR = os.path.join(os.path.dirname(__file__), "fixtures")


@pytest.fixture(autouse=True)
def _reset_service_registry():
    """Drop per-config services (cache + clients + region cache) between tests."""
    yield
    _tools_module._services_registry.clear()


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


@pytest.fixture(scope="session")
def comprehensive_ref_fasta_path():
    """Path to the comprehensive reference FASTA. Generated on first use if absent.

    CI (and fresh checkouts) won't have the file until something asks for it,
    so we lazily create it here instead of relying on the e2e session fixture.
    """
    path = os.path.join(FIXTURES_DIR, "comprehensive_ref.fa")
    if not os.path.exists(path):
        from tests.create_fixtures import create_comprehensive_reference

        create_comprehensive_reference()
    return path


@pytest.fixture(scope="session")
def comprehensive_bam_path(comprehensive_ref_fasta_path):
    """Path to the comprehensive test BAM. Generated on first use if absent."""
    path = os.path.join(FIXTURES_DIR, "comprehensive.bam")
    if not os.path.exists(path):
        from tests.create_fixtures import create_comprehensive_bam

        create_comprehensive_bam(comprehensive_ref_fasta_path)
    return path
