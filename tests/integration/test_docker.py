"""Unit and integration tests for Docker deployment configuration.

These tests validate the Docker infrastructure files without requiring
a Docker daemon — they parse and verify Dockerfile syntax, compose config,
.dockerignore coverage, health check behavior, and entrypoint logic.
"""

import http.server
import importlib.util
import os
import re
import socket
import subprocess
import sys
import threading

import pytest

PROJECT_ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------


def read_project_file(name: str) -> str:
    path = os.path.join(PROJECT_ROOT, name)
    with open(path) as f:
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
            base = line.split(" AS")[0].strip() if "AS" in line else line.strip()
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
        """Non-threshold config env vars should have defaults in the Dockerfile."""
        for var in [
            "BAMCP_REFERENCE",
            "BAMCP_MAX_READS",
            "BAMCP_DEFAULT_WINDOW",
        ]:
            assert var in self.content, f"Missing ENV declaration for {var}"

    @pytest.mark.unit
    def test_calling_thresholds_not_overridden(self):
        """The image must not pin calling thresholds: the code defaults are the source of truth.
        An ENV BAMCP_MIN_MAPQ="0" here once silently disabled the MAPQ filter in prod."""
        for var in [
            "BAMCP_MIN_VAF",
            "BAMCP_MIN_DEPTH",
            "BAMCP_MIN_MAPQ",
            "BAMCP_MIN_BASEQ",
            "BAMCP_MIN_ALT_COUNT",
        ]:
            assert f"ENV {var}" not in self.content, f"Dockerfile pins {var}"

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
        assert "/tmp" in self.content  # noqa: S108

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
        """Lint service should run ruff and mypy."""
        assert "ruff" in self.content
        assert "mypy" in self.content

    @pytest.mark.unit
    def test_http_services_expose_port(self):
        """Beta and prod services should expose port 8000."""
        assert '"8000:8000"' in self.content

    @pytest.mark.unit
    def test_prod_profile_is_locked_down(self):
        """Prod binds loopback, rejects local paths, and allow-lists remote fetch hosts."""
        import yaml

        prod = yaml.safe_load(self.content)["services"]["prod"]
        assert prod["ports"] == ["127.0.0.1:8000:8000"]
        env = dict(e.split("=", 1) for e in prod["environment"])
        assert env["BAMCP_ALLOW_REMOTE_FILES"] == "true"
        assert env["BAMCP_ALLOW_LOCAL_FILES"] == "false"
        assert env["BAMCP_ALLOWED_REMOTE_HOSTS"], "remote fetch without an allow-list"

    @pytest.mark.unit
    def test_prod_allow_list_matches_setup_cloudflared(self):
        """The compose prod allow-list and setup-cloudflared.sh's default must stay in sync."""
        import yaml

        prod = yaml.safe_load(self.content)["services"]["prod"]
        env = dict(e.split("=", 1) for e in prod["environment"])
        script = read_project_file("scripts/setup-cloudflared.sh")
        match = re.search(r'^DEFAULT_ALLOWED_REMOTE_HOSTS="([^"]+)"', script, re.M)
        assert match
        compose_hosts = set(env["BAMCP_ALLOWED_REMOTE_HOSTS"].split(","))
        assert compose_hosts == set(match.group(1).split(","))
        # Set membership of allow-list entries, not URL matching.
        assert {"ftp.ensembl.org"} <= compose_hosts


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
        lines = [
            line.strip()
            for line in self.content.splitlines()
            if line.strip() and not line.startswith("#")
        ]
        assert "src" not in lines and "src/" not in lines

    @pytest.mark.unit
    def test_does_not_exclude_pyproject(self):
        """pyproject.toml must NOT be excluded."""
        lines = [
            line.strip()
            for line in self.content.splitlines()
            if line.strip() and not line.startswith("#")
        ]
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
            import healthcheck  # type: ignore[import-not-found]

            assert hasattr(healthcheck, "check_health")
        finally:
            sys.path.pop(0)

    @pytest.mark.integration
    def test_health_check_passes(self):
        """Running the health check (stdio transport) should succeed in the test environment."""
        env = {**os.environ, "BAMCP_TRANSPORT": "stdio"}
        result = subprocess.run(  # noqa: S603
            [sys.executable, os.path.join(PROJECT_ROOT, "docker", "healthcheck.py")],
            capture_output=True,
            text=True,
            timeout=10,
            env=env,
        )
        assert result.returncode == 0, f"Health check failed: {result.stderr}"

    @pytest.mark.integration
    def test_preflight_ignores_http_transport(self):
        """--preflight (the entrypoint's mode) runs before the server listens: import check only."""
        env = {**os.environ, "BAMCP_TRANSPORT": "sse", "BAMCP_PORT": str(_free_port())}
        result = subprocess.run(  # noqa: S603
            [sys.executable, os.path.join(PROJECT_ROOT, "docker", "healthcheck.py"), "--preflight"],
            capture_output=True,
            text=True,
            timeout=10,
            env=env,
        )
        assert result.returncode == 0, f"Preflight failed: {result.stderr}"

    @pytest.mark.integration
    def test_health_check_verifies_server(self, monkeypatch):
        """The stdio health check should import bamcp and instantiate the server."""
        monkeypatch.setenv("BAMCP_TRANSPORT", "stdio")
        assert _load_healthcheck().check_health() is True


def _load_healthcheck():
    """Import docker/healthcheck.py by path (docker/ is not a package)."""
    spec = importlib.util.spec_from_file_location(
        "bamcp_docker_healthcheck", os.path.join(PROJECT_ROOT, "docker", "healthcheck.py")
    )
    assert spec and spec.loader
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def _free_port() -> int:
    with socket.socket() as s:
        s.bind(("127.0.0.1", 0))
        return s.getsockname()[1]


class _ProbeHandler(http.server.BaseHTTPRequestHandler):
    status = 200
    seen: list[tuple[str, str | None]] = []

    def do_GET(self):  # noqa: N802
        type(self).seen.append((self.path, self.headers.get("Host")))
        self.send_response(type(self).status)
        self.end_headers()
        self.wfile.write(b'{"sha": null}')

    def log_message(self, *args):
        pass


@pytest.fixture
def probe_server():
    """A loopback HTTP server standing in for BAMCP; yields (port, handler class)."""

    class Handler(_ProbeHandler):
        status = 200
        seen = []

    server = http.server.HTTPServer(("127.0.0.1", 0), Handler)
    thread = threading.Thread(target=server.serve_forever, daemon=True)
    thread.start()
    try:
        yield server.server_address[1], Handler
    finally:
        server.shutdown()
        server.server_close()


class TestHealthCheckHttpProbe:
    """The runtime probe (Dockerfile HEALTHCHECK) for HTTP transports hits GET /__version."""

    @pytest.fixture(autouse=True)
    def _env(self, monkeypatch):
        for var in ("BAMCP_PORT", "PORT", "BAMCP_TRUSTED_HOSTS"):
            monkeypatch.delenv(var, raising=False)
        monkeypatch.setenv("BAMCP_TRANSPORT", "sse")

    @pytest.mark.unit
    @pytest.mark.parametrize("transport", ["sse", "streamable-http"])
    def test_http_transport_probes_version(self, monkeypatch, probe_server, transport):
        port, handler = probe_server
        monkeypatch.setenv("BAMCP_TRANSPORT", transport)
        monkeypatch.setenv("BAMCP_PORT", str(port))
        hc = _load_healthcheck()
        monkeypatch.setattr(hc, "_check_imports", lambda: pytest.fail("stdio path used"))
        assert hc.check_health() is True
        assert handler.seen[0][0] == "/__version"

    @pytest.mark.unit
    def test_non_200_is_unhealthy(self, monkeypatch, probe_server):
        port, handler = probe_server
        handler.status = 503
        monkeypatch.setenv("BAMCP_PORT", str(port))
        assert _load_healthcheck().check_health() is False

    @pytest.mark.unit
    def test_nothing_listening_is_unhealthy(self, monkeypatch):
        monkeypatch.setenv("BAMCP_PORT", str(_free_port()))
        assert _load_healthcheck().check_health() is False

    @pytest.mark.unit
    def test_cloud_port_used_when_bamcp_port_unset(self, monkeypatch, probe_server):
        """Mirrors entrypoint.sh: PORT applies only when BAMCP_PORT is unset."""
        port, _ = probe_server
        monkeypatch.setenv("PORT", str(port))
        assert _load_healthcheck().check_health() is True

    @pytest.mark.unit
    def test_bamcp_port_wins_over_cloud_port(self, monkeypatch, probe_server):
        port, _ = probe_server
        monkeypatch.setenv("PORT", str(_free_port()))
        monkeypatch.setenv("BAMCP_PORT", str(port))
        assert _load_healthcheck().check_health() is True

    @pytest.mark.unit
    @pytest.mark.parametrize(
        ("trusted", "expected"),
        [
            ("bamcp.example.com,other", "bamcp.example.com"),
            ("*.example.com", "healthcheck.example.com"),
        ],
    )
    def test_sends_a_trusted_host_header(self, monkeypatch, probe_server, trusted, expected):
        """TrustedHostMiddleware would 400 a bare 127.0.0.1 Host when BAMCP_TRUSTED_HOSTS is set."""
        port, handler = probe_server
        monkeypatch.setenv("BAMCP_PORT", str(port))
        monkeypatch.setenv("BAMCP_TRUSTED_HOSTS", trusted)
        assert _load_healthcheck().check_health() is True
        assert handler.seen[0][1] == expected

    @pytest.mark.unit
    def test_stdio_uses_import_check(self, monkeypatch):
        monkeypatch.setenv("BAMCP_TRANSPORT", "stdio")
        hc = _load_healthcheck()
        monkeypatch.setattr(hc, "_check_imports", lambda: True)
        monkeypatch.setattr(hc, "_check_http", lambda: pytest.fail("HTTP probed for stdio"))
        assert hc.check_health() is True

    @pytest.mark.unit
    def test_entrypoint_runs_preflight_mode(self):
        assert "healthcheck.py --preflight" in read_project_file("docker/entrypoint.sh")


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
            "BAMCP_REFERENCE": config.reference or "",
            "BAMCP_MAX_READS": str(config.max_reads),
            "BAMCP_DEFAULT_WINDOW": str(config.default_window),
        }

        # Values, not just names: an image default that differs from the code default silently
        # changes behavior in prod (it once pinned BAMCP_MIN_MAPQ="0" against a code default of 20).
        for var, expected in env_map.items():
            match = re.search(rf'^ENV {var}="([^"]*)"', dockerfile, re.MULTILINE)
            assert match, f"Dockerfile missing ENV {var}"
            assert match.group(1) == expected, f"ENV {var}={match.group(1)!r} != code {expected!r}"

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
        lines = [
            line.strip()
            for line in dockerignore.splitlines()
            if line.strip() and not line.startswith("#")
        ]
        assert "tests" not in lines and "tests/" not in lines, (
            "tests/ must not be in .dockerignore — Dockerfile.dev needs it"
        )


def test_node_version_is_pinned_once():
    """CI and the image build the viewer with the same Node major (src/bamcp/static/.nvmrc)."""
    with open(os.path.join(PROJECT_ROOT, "src", "bamcp", "static", ".nvmrc")) as f:
        pinned = f.read().strip()
    with open(os.path.join(PROJECT_ROOT, "Dockerfile")) as f:
        image = re.search(r"^FROM node:(\d+)", f.read(), re.MULTILINE)
    assert image is not None and image.group(1) == pinned
    workflows = os.path.join(PROJECT_ROOT, ".github", "workflows")
    for name in os.listdir(workflows):
        with open(os.path.join(workflows, name)) as f:
            text = f.read()
        assert "node-version:" not in text, f"{name} pins Node inline; use node-version-file"
        if "make install" in text:
            assert "node-version-file: src/bamcp/static/.nvmrc" in text, name
