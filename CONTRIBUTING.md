# Contributing to BAMCP

Thank you for your interest in contributing to BAMCP! This document provides guidelines for contributing.

## Development Setup

```bash
git clone https://github.com/RTrentJones/BAMCP.git
cd BAMCP
python -m venv venv
source venv/bin/activate
make install                          # needs Node/npm: builds the viewer bundle
pip install -e ".[telemetry,eval]"    # extras the test suite and eval-smoke need
python tests/create_fixtures.py
```

### Pre-commit Hooks

```bash
pip install pre-commit
pre-commit install
# Run against all files
pre-commit run --all-files
```

## Running Tests

```bash
# Unit and integration tests
make test

# With HTML coverage report
make coverage

# E2E tests (requires Playwright)
make test-e2e
```

## Docker Development

```bash
# Build and run tests in Docker
make docker-test

# Build production and dev images
make docker-build
```

## Code Style

```bash
# Check for violations
make lint

# Auto-format
make format

# Type checking
make typecheck
```

## Makefile Targets

| Target | Description |
|--------|-------------|
| `make verify` | Every ci.yml gate in one command (`VERIFY_E2E=1` adds Playwright) |
| `make install` | Build the viewer, then install the package with dev dependencies |
| `make build-viewer` | Build the viewer bundle (`static/dist/viewer.html`) |
| `make test` | Run unit and integration tests |
| `make test-network` | Live-network tests only (excluded from the default suite) |
| `make test-e2e` | Run E2E Playwright tests |
| `make lint` | Check formatting and linting |
| `make format` | Auto-format code |
| `make typecheck` | Run mypy type checking |
| `make docker-build` | Build Docker images |
| `make docker-test` | Run tests in Docker |
| `make coverage` | Generate HTML coverage report |
| `make coverage-strict` | Coverage with the 85% CI bar |
| `make eval-smoke` | Deterministic ground-truth gate (required ship-gate) |
| `make render-viewer` | Render the viewer to a PNG (`.render-dev/last.png`) |
| `make clean` | Remove build artifacts |

## Pull Request Process

1. Fork the repository
2. Create a branch off `main` (`git checkout -b feat/amazing-feature` or `fix/<slug>`; there is no `develop` branch)
3. Write tests for your changes
4. Ensure all checks pass (`make verify` runs every ci.yml gate: lockfile, viewer typecheck/lint/test/build, ruff format + check, mypy, tests at the 85% coverage bar, eval-smoke)
5. Commit in conventional style (`git commit -m 'feat(viewer): amazing feature'`); the merged commit carries the PR number, e.g. `fix(viewer): distinct ruler tick labels at zoom (#53)`
6. Push to the branch (`git push origin feat/amazing-feature`)
7. Open a Pull Request against `main`

## Reporting Issues

Please use [GitHub Issues](https://github.com/RTrentJones/BAMCP/issues) to report bugs or request features.
