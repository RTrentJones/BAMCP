.PHONY: install test test-e2e test-all lint format typecheck docker-build docker-test clean coverage coverage-strict build-viewer eval eval-cached eval-compare eval-dry eval-vision-setup render-viewer render-viewer-all-modes

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

coverage-strict:
	python -m pytest tests/ --ignore=tests/e2e -v --cov=bamcp --cov-fail-under=85 --cov-report=term-missing

# Eval harness — MARRVEL-MCP-compatible runner against the BAMCP server.
eval:
	python scripts/run_eval.py --output-dir .eval-results

eval-cached:
	python scripts/run_eval.py --cache --output-dir .eval-results

eval-compare:
	python scripts/run_eval.py --with-vanilla --with-rendering-comparison --output-dir .eval-results

eval-dry:
	python scripts/run_eval.py --dry-run --output-dir .eval-results

eval-vision-setup:
	python -m playwright install chromium

# Developer feedback loop — render the viewer to a PNG without a browser.
# Defaults render the comprehensive fixture; override BAM/REF/REGION/MODE/OUT
# to point at your own data. Pair with `Read` on the output path to verify.
BAM ?= tests/fixtures/comprehensive.bam
REF ?= tests/fixtures/comprehensive_ref.fa
REGION ?= chr1:1000-1300
MODE ?= expanded
OUT ?= .render-dev/last.png

render-viewer: build-viewer
	python scripts/render_viewer.py \
		--bam $(BAM) --reference $(REF) --region $(REGION) \
		--mode $(MODE) --out $(OUT)

render-viewer-all-modes: build-viewer
	@mkdir -p .render-dev
	@for mode in expanded dv-strips dv-composite; do \
		python scripts/render_viewer.py \
			--bam $(BAM) --reference $(REF) --region $(REGION) \
			--mode $$mode --out .render-dev/$$mode.png ; \
	done
