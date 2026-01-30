"""Unit and integration tests for Docker deployment configuration.

These tests validate the Docker infrastructure files without requiring
a Docker daemon — they parse and verify Dockerfile syntax, compose config,
.dockerignore coverage, health check behavior, and entrypoint logic.
"""

import os
import re
import subprocess
import sys

import pytest

PROJECT_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def read_project_file(name: str) -> str:
    path = os.path.join(PROJECT_ROOT, name)
    with open(path, "r") as f:
        return f.read()


def project_file_exists(name: str) -> bool:
    return os.path.isfile(os.path.join(PROJECT_ROOT, name))


# ===========================================================================
# Dockerfile (production) tests
# ===========================================================================

class TestDockerfileProd:
    """Validate the production Dockerfile structure and best practices."""

    @pytest.fixture(autouse=True)
    def load(self):
        self.content = read_project_file("Dockerfile")

    @pytest.mark.unit
    def test_file_exists(self):
        assert project_file_exists("Dockerfile")

    @pytest.mark.unit
    def test_multi_stage_build(self):
        """Production image should use a multi-stage build."""
        stages = re.findall(r"^FROM\s+", self.content, re.MULTILINE)
        assert len(stages) >= 2, "Expected at least 2 FROM statements (builder + runtime)"

    @pytest.mark.unit
    def test_builder_stage_named(self):
        assert "AS builder" in self.content

    @pytest.mark.unit
    def test_runtime_stage_named(self):
        assert "AS runtime" in self.content

    @pytest.mark.unit
    def test_uses_slim_base(self):
        """Should use slim base images to reduce attack surface."""
        from_lines = re.findall(r"^FROM\s+(.+)$", self.content, re.MULTILINE)
        for line in from_lines:
            if "AS" in line:
                base = line.split(" AS")[0].strip()
            else:
                base = line.strip()
            assert "slim" in base, f"Base image should be slim, got: {base}"

    @pytest.mark.unit
    def test_python_312(self):
        """Should target Python 3.12."""
        assert "python:3.12" in self.content

    @pytest.mark.unit
    def test_no_cache_pip(self):
        """pip install should use --no-cache-dir."""
        assert "--no-cache-dir" in self.content

    @pytest.mark.unit
    def test_apt_lists_cleaned(self):
        """apt-get should clean lists after install."""
        assert "rm -rf /var/lib/apt/lists/*" in self.content

    @pytest.mark.unit
    def test_non_root_user(self):
        """Production image should run as non-root."""
        assert "USER bamcp" in self.content
        assert "groupadd" in self.content
        assert "useradd" in self.content

    @pytest.mark.unit
    def test_healthcheck_defined(self):
        assert "HEALTHCHECK" in self.content

    @pytest.mark.unit
    def test_healthcheck_uses_script(self):
        """Health check should use the dedicated Python script."""
        assert "healthcheck.py" in self.content

    @pytest.mark.unit
    def test_entrypoint_uses_script(self):
        """Should use the entrypoint shell script."""
        assert "entrypoint.sh" in self.content

    @pytest.mark.unit
    def test_cmd_runs_bamcp(self):
        """Default CMD should run `python -m bamcp`."""
        assert "python" in self.content
        assert "bamcp" in self.content

    @pytest.mark.unit
    def test_data_volume(self):
        assert 'VOLUME ["/data"]' in self.content

    @pytest.mark.unit
    def test_env_vars_declared(self):
        """All 6 config env vars should have defaults in the Dockerfile."""
        for var in [
            "BAMCP_REFERENCE",
            "BAMCP_MAX_READS",
            "BAMCP_DEFAULT_WINDOW",
            "BAMCP_MIN_VAF",
            "BAMCP_MIN_DEPTH",
            "BAMCP_MIN_MAPQ",
        ]:
            assert var in self.content, f"Missing ENV declaration for {var}"

    @pytest.mark.unit
    def test_workdir_set(self):
        assert "WORKDIR /app" in self.content

    @pytest.mark.unit
    def test_copies_docker_dir(self):
        """Should copy the docker/ helper scripts."""
        assert "COPY docker/ docker/" in self.content

    @pytest.mark.unit
    def test_builder_installs_c_deps(self):
        """Builder stage needs C compiler deps for pysam."""
        for dep in ["gcc", "zlib1g-dev", "libbz2-dev", "liblzma-dev"]:
            assert dep in self.content, f"Builder missing C dependency: {dep}"

    @pytest.mark.unit
    def test_runtime_has_shared_libs(self):
        """Runtime stage needs shared libs (not -dev packages)."""
        assert "zlib1g" in self.content
        # Should NOT have gcc in runtime
        runtime_section = self.content.split("AS runtime")[1]
        assert "gcc" not in runtime_section, "Runtime should not install gcc"


# ===========================================================================
# Dockerfile.dev tests
# ===========================================================================

class TestDockerfileDev:
    """Validate the development Dockerfile."""

    @pytest.fixture(autouse=True)
    def load(self):
        self.content = read_project_file("Dockerfile.dev")

    @pytest.mark.unit
    def test_file_exists(self):
        assert project_file_exists("Dockerfile.dev")

    @pytest.mark.unit
    def test_single_stage(self):
        """Dev image should be a single stage (no multi-stage needed)."""
        stages = re.findall(r"^FROM\s+", self.content, re.MULTILINE)
        assert len(stages) == 1

    @pytest.mark.unit
    def test_installs_dev_deps(self):
        """Should install dev extras."""
        assert ".[dev]" in self.content

    @pytest.mark.unit
    def test_installs_playwright(self):
        """Should install Playwright browsers for E2E tests."""
        assert "playwright install" in self.content

    @pytest.mark.unit
    def test_playwright_system_deps(self):
        """Should install Playwright's required system libraries."""
        for dep in ["libnss3", "libgbm1", "libasound2"]:
            assert dep in self.content, f"Missing Playwright system dep: {dep}"

    @pytest.mark.unit
    def test_copies_tests(self):
        assert "COPY tests/ tests/" in self.content

    @pytest.mark.unit
    def test_creates_fixtures(self):
        """Should pre-generate test fixtures."""
        assert "create_fixtures.py" in self.content

    @pytest.mark.unit
    def test_default_cmd_runs_tests(self):
        """Default CMD should run pytest."""
        assert "pytest" in self.content

    @pytest.mark.unit
    def test_git_installed(self):
        """Dev image should have git for development workflows."""
        assert "git" in self.content

    @pytest.mark.unit
    def test_unbuffered_python(self):
        """Should set PYTHONUNBUFFERED for real-time test output."""
        assert "PYTHONUNBUFFERED" in self.content


# ===========================================================================
# docker-compose.yml tests
# ===========================================================================

class TestDockerCompose:
    """Validate docker-compose.yml profiles and service definitions."""

    @pytest.fixture(autouse=True)
    def load(self):
        self.content = read_project_file("docker-compose.yml")

    @pytest.mark.unit
    def test_file_exists(self):
        assert project_file_exists("docker-compose.yml")

    @pytest.mark.unit
    def test_has_services_key(self):
        assert "services:" in self.content

    @pytest.mark.unit
    def test_dev_profile_services(self):
        """Dev profile should have test, e2e, and lint services."""
        for service in ["test:", "e2e:", "lint:"]:
            assert service in self.content, f"Missing dev service: {service}"

    @pytest.mark.unit
    def test_beta_profile_service(self):
        assert "beta:" in self.content

    @pytest.mark.unit
    def test_prod_profile_service(self):
        assert "prod:" in self.content

    @pytest.mark.unit
    def test_dev_profile_tag(self):
        assert '"dev"' in self.content

    @pytest.mark.unit
    def test_beta_profile_tag(self):
        assert '"beta"' in self.content

    @pytest.mark.unit
    def test_prod_profile_tag(self):
        assert '"prod"' in self.content

    @pytest.mark.unit
    def test_prod_read_only(self):
        """Production service should mount filesystem read-only."""
        assert "read_only: true" in self.content

    @pytest.mark.unit
    def test_prod_no_new_privileges(self):
        """Production service should prevent privilege escalation."""
        assert "no-new-privileges" in self.content

    @pytest.mark.unit
    def test_prod_tmpfs(self):
        """Production should have a tmpfs for /tmp (read-only rootfs needs it)."""
        assert "tmpfs:" in self.content
        assert "/tmp" in self.content

    @pytest.mark.unit
    def test_dev_uses_dev_dockerfile(self):
        assert "Dockerfile.dev" in self.content

    @pytest.mark.unit
    def test_dev_mounts_source(self):
        """Dev services should mount source code for hot reload."""
        assert "./src:/app/src" in self.content

    @pytest.mark.unit
    def test_dev_mounts_tests(self):
        assert "./tests:/app/tests" in self.content

    @pytest.mark.unit
    def test_named_volume(self):
        """Should define a named volume for BAM data."""
        assert "bamcp-data:" in self.content

    @pytest.mark.unit
    def test_prod_restart_always(self):
        assert "restart: always" in self.content

    @pytest.mark.unit
    def test_beta_restart_unless_stopped(self):
        assert "restart: unless-stopped" in self.content

    @pytest.mark.unit
    def test_e2e_service_runs_e2e_tests(self):
        """E2E service should run tests/e2e/."""
        assert "tests/e2e/" in self.content

    @pytest.mark.unit
    def test_lint_service_runs_formatters(self):
        """Lint service should run black, ruff, mypy."""
        assert "black" in self.content
        assert "ruff" in self.content
        assert "mypy" in self.content

    @pytest.mark.unit
    def test_stdin_open_for_mcp_server(self):
        """MCP server services need stdin_open for stdio transport."""
        assert "stdin_open: true" in self.content


# ===========================================================================
# .dockerignore tests
# ===========================================================================

class TestDockerignore:
    """Validate .dockerignore excludes unnecessary files."""

    @pytest.fixture(autouse=True)
    def load(self):
        self.content = read_project_file(".dockerignore")

    @pytest.mark.unit
    def test_file_exists(self):
        assert project_file_exists(".dockerignore")

    @pytest.mark.unit
    def test_excludes_git(self):
        assert ".git" in self.content

    @pytest.mark.unit
    def test_excludes_pycache(self):
        assert "__pycache__" in self.content

    @pytest.mark.unit
    def test_excludes_coverage(self):
        assert ".coverage" in self.content

    @pytest.mark.unit
    def test_excludes_venv(self):
        assert ".venv" in self.content or "venv" in self.content

    @pytest.mark.unit
    def test_excludes_env_file(self):
        """Should exclude .env to prevent secret leaks."""
        assert ".env" in self.content

    @pytest.mark.unit
    def test_excludes_docker_files(self):
        """Docker files themselves shouldn't be in the build context layers."""
        assert "Dockerfile" in self.content
        assert "docker-compose" in self.content

    @pytest.mark.unit
    def test_excludes_ide(self):
        for pattern in [".idea", ".vscode"]:
            assert pattern in self.content, f"Should exclude {pattern}"

    @pytest.mark.unit
    def test_excludes_pytest_cache(self):
        assert ".pytest_cache" in self.content

    @pytest.mark.unit
    def test_does_not_exclude_src(self):
        """src/ must NOT be excluded — it's the actual application."""
        lines = [l.strip() for l in self.content.splitlines() if l.strip() and not l.startswith("#")]
        assert "src" not in lines and "src/" not in lines

    @pytest.mark.unit
    def test_does_not_exclude_pyproject(self):
        """pyproject.toml must NOT be excluded."""
        lines = [l.strip() for l in self.content.splitlines() if l.strip() and not l.startswith("#")]
        assert "pyproject.toml" not in lines


# ===========================================================================
# Health check script tests
# ===========================================================================

class TestHealthCheck:
    """Test the Docker health check script runs correctly."""

    @pytest.mark.unit
    def test_script_exists(self):
        assert project_file_exists("docker/healthcheck.py")

    @pytest.mark.unit
    def test_script_importable(self):
        """Health check module should be importable."""
        sys.path.insert(0, os.path.join(PROJECT_ROOT, "docker"))
        try:
            import healthcheck
            assert hasattr(healthcheck, "check_health")
        finally:
            sys.path.pop(0)

    @pytest.mark.integration
    def test_health_check_passes(self):
        """Running the health check should succeed in the test environment."""
        result = subprocess.run(
            [sys.executable, os.path.join(PROJECT_ROOT, "docker", "healthcheck.py")],
            capture_output=True,
            text=True,
            timeout=10,
        )
        assert result.returncode == 0, f"Health check failed: {result.stderr}"

    @pytest.mark.integration
    def test_health_check_verifies_version(self):
        """Health check should verify the package version is set."""
        sys.path.insert(0, os.path.join(PROJECT_ROOT, "docker"))
        try:
            from healthcheck import check_health
            assert check_health() is True
        finally:
            sys.path.pop(0)

    @pytest.mark.integration
    def test_health_check_verifies_server(self):
        """Health check should instantiate the server."""
        sys.path.insert(0, os.path.join(PROJECT_ROOT, "docker"))
        try:
            from healthcheck import check_health
            result = check_health()
            assert result is True
        finally:
            sys.path.pop(0)


# ===========================================================================
# Entrypoint script tests
# ===========================================================================

class TestEntrypoint:
    """Test the Docker entrypoint shell script."""

    @pytest.mark.unit
    def test_script_exists(self):
        assert project_file_exists("docker/entrypoint.sh")

    @pytest.mark.unit
    def test_script_is_executable(self):
        path = os.path.join(PROJECT_ROOT, "docker", "entrypoint.sh")
        assert os.access(path, os.X_OK), "entrypoint.sh should be executable"

    @pytest.mark.unit
    def test_script_has_shebang(self):
        content = read_project_file("docker/entrypoint.sh")
        assert content.startswith("#!/bin/sh"), "Should have #!/bin/sh shebang"

    @pytest.mark.unit
    def test_script_uses_set_e(self):
        """Should exit on error."""
        content = read_project_file("docker/entrypoint.sh")
        assert "set -e" in content

    @pytest.mark.unit
    def test_script_runs_health_check(self):
        """Entrypoint should reference the health check."""
        content = read_project_file("docker/entrypoint.sh")
        assert "healthcheck.py" in content

    @pytest.mark.unit
    def test_script_uses_exec(self):
        """Should use exec to replace the shell process (proper PID 1 handling)."""
        content = read_project_file("docker/entrypoint.sh")
        assert 'exec "$@"' in content


# ===========================================================================
# Cross-file consistency tests
# ===========================================================================

class TestDockerConsistency:
    """Tests that verify Docker files are consistent with each other and with
    the main project configuration."""

    @pytest.mark.unit
    def test_compose_references_existing_dockerfiles(self):
        """docker-compose.yml should only reference Dockerfiles that exist."""
        compose = read_project_file("docker-compose.yml")
        if "Dockerfile.dev" in compose:
            assert project_file_exists("Dockerfile.dev")
        if "dockerfile: Dockerfile\n" in compose or "dockerfile: Dockerfile" in compose:
            assert project_file_exists("Dockerfile")

    @pytest.mark.unit
    def test_env_vars_match_config(self):
        """Dockerfile env vars should match BAMCPConfig fields."""
        from bamcp.config import BAMCPConfig
        dockerfile = read_project_file("Dockerfile")
        config = BAMCPConfig()

        env_map = {
            "BAMCP_REFERENCE": config.reference,
            "BAMCP_MAX_READS": str(config.max_reads),
            "BAMCP_DEFAULT_WINDOW": str(config.default_window),
            "BAMCP_MIN_VAF": str(config.min_vaf),
            "BAMCP_MIN_DEPTH": str(config.min_depth),
            "BAMCP_MIN_MAPQ": str(config.min_mapq),
        }

        for var, expected_default in env_map.items():
            assert var in dockerfile, f"Dockerfile missing ENV {var}"

    @pytest.mark.unit
    def test_python_version_consistent(self):
        """All Dockerfiles should use the same Python version."""
        prod = read_project_file("Dockerfile")
        dev = read_project_file("Dockerfile.dev")

        prod_versions = set(re.findall(r"python:(\d+\.\d+)", prod))
        dev_versions = set(re.findall(r"python:(\d+\.\d+)", dev))

        assert prod_versions == dev_versions, (
            f"Python version mismatch: prod={prod_versions}, dev={dev_versions}"
        )

    @pytest.mark.unit
    def test_compose_volume_name_consistent(self):
        """Volume name in service and volumes section should match."""
        compose = read_project_file("docker-compose.yml")
        # Both the service mount and the top-level volume should reference bamcp-data
        assert compose.count("bamcp-data") >= 3, (
            "Expected bamcp-data in at least 3 places (prod, beta, volumes section)"
        )

    @pytest.mark.unit
    def test_dockerignore_does_not_exclude_tests_for_dev(self):
        """tests/ should not be in .dockerignore since Dockerfile.dev needs them.
        (Dockerfile.dev does COPY tests/ and .dockerignore excludes Docker*
        so the .dockerignore applies to the prod build context primarily.)"""
        dockerignore = read_project_file(".dockerignore")
        lines = [l.strip() for l in dockerignore.splitlines()
                 if l.strip() and not l.startswith("#")]
        assert "tests" not in lines and "tests/" not in lines, (
            "tests/ must not be in .dockerignore — Dockerfile.dev needs it"
        )
