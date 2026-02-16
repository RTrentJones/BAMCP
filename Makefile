.PHONY: install test test-e2e test-all lint format typecheck docker-build docker-test clean coverage build-viewer

build-viewer:
	cd src/bamcp/static && npm install && npm run build

install: build-viewer
	pip install -e ".[dev]"

test:
	python -m pytest tests/ --ignore=tests/e2e -v

test-e2e:
	playwright install chromium
	python -m pytest tests/e2e/ -v --no-cov

test-all: test test-e2e

lint:
	ruff format --check src tests
	ruff check src tests

format:
	ruff format src tests
	ruff check --fix src tests

typecheck:
	mypy src

docker-build:
	docker compose --profile prod build
	docker compose --profile dev build

docker-test:
	docker compose --profile dev run --rm test
	docker compose --profile dev run --rm lint

clean:
	rm -rf .pytest_cache .mypy_cache .ruff_cache htmlcov .coverage .coverage.* dist build *.egg-info
	rm -rf src/bamcp/static/dist src/bamcp/static/node_modules
	find . -type d -name __pycache__ -exec rm -rf {} + 2>/dev/null || true

coverage:
	python -m pytest tests/ --ignore=tests/e2e -v --cov=bamcp --cov-report=html
	@echo "Coverage report: htmlcov/index.html"
