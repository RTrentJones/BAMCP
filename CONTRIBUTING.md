# Contributing to BAMCP

Thank you for your interest in contributing to BAMCP! This document provides guidelines for contributing.

## Development Setup

```bash
git clone https://github.com/RTrentJones/BAMCP.git
cd BAMCP
python -m venv venv
source venv/bin/activate
make install
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
| `make install` | Install package with dev dependencies |
| `make test` | Run unit and integration tests |
| `make test-e2e` | Run E2E Playwright tests |
| `make lint` | Check formatting and linting |
| `make format` | Auto-format code |
| `make typecheck` | Run mypy type checking |
| `make docker-build` | Build Docker images |
| `make docker-test` | Run tests in Docker |
| `make coverage` | Generate HTML coverage report |
| `make clean` | Remove build artifacts |

## Pull Request Process

1. Fork the repository
2. Create a feature branch (`git checkout -b feature/amazing-feature`)
3. Write tests for your changes
4. Ensure all checks pass (`make lint && make typecheck && make test`)
5. Commit your changes (`git commit -m 'Add amazing feature'`)
6. Push to the branch (`git push origin feature/amazing-feature`)
7. Open a Pull Request

## Reporting Issues

Please use [GitHub Issues](https://github.com/RTrentJones/BAMCP/issues) to report bugs or request features.
