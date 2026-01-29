# Contributing to BAMCP

Thank you for your interest in contributing to BAMCP! This document provides guidelines for contributing.

## Development Setup

```bash
git clone https://github.com/yourusername/bamcp.git
cd bamcp
python -m venv venv
source venv/bin/activate
pip install -e ".[dev]"
```

## Running Tests

```bash
# Unit and integration tests
pytest

# With coverage
pytest --cov=bamcp

# E2E tests (requires Playwright)
playwright install chromium
pytest tests/e2e/
```

## Code Style

```bash
black src tests
isort src tests
ruff check src tests
mypy src
```

## Pull Request Process

1. Fork the repository
2. Create a feature branch (`git checkout -b feature/amazing-feature`)
3. Write tests for your changes
4. Ensure all tests pass
5. Commit your changes (`git commit -m 'Add amazing feature'`)
6. Push to the branch (`git push origin feature/amazing-feature`)
7. Open a Pull Request

## Reporting Issues

Please use [GitHub Issues](https://github.com/yourusername/bamcp/issues) to report bugs or request features.
