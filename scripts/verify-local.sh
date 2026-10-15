#!/usr/bin/env bash
# Local regression loop — every gate CI's ci.yml enforces, in one command.
#
#   ./scripts/verify-local.sh            # full loop
#   VERIFY_SKIP_VIEWER=1 ./scripts/...   # skip the npm gates (Python-only change)
#   VERIFY_E2E=1 ./scripts/...           # also run the playwright e2e suite
#
# Runs against THIS checkout's src/ (also from a git worktree), never a
# different editable install, and leaves the tracked fixtures untouched.
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$ROOT"
export PYTHONPATH="$ROOT/src${PYTHONPATH:+:$PYTHONPATH}"

step() { printf '\n\033[1m==> %s\033[0m\n' "$*"; }

step "import path"
python -c "import bamcp, sys; p = bamcp.__file__; print(p); sys.exit(0 if p.startswith('$ROOT/src/') else 1)" \
  || { echo "bamcp imports from outside this checkout" >&2; exit 1; }

step "lockfile (uv lock --check)"
if command -v uv >/dev/null; then uv lock --check; else echo "uv not installed — skipped"; fi

step "fixtures"
python tests/create_fixtures.py >/dev/null
# create_fixtures.py rewrites the tracked fixtures non-byte-identically; keep the
# committed copies so the loop never produces binary churn.
git ls-files -z tests/fixtures | xargs -0 git restore --source=HEAD --

if [[ "${VERIFY_SKIP_VIEWER:-0}" != "1" ]]; then
  step "viewer (npm ci, typecheck, lint, test, build)"
  (
    cd src/bamcp/static
    npm ci --no-audit --no-fund --loglevel=error
    npm run typecheck
    npm run lint
    npm test
    npm run build
  )
  python scripts/render_viewer.py --record-build-hash
elif [[ ! -f src/bamcp/static/dist/viewer.html ]]; then
  echo "VERIFY_SKIP_VIEWER=1 but dist/viewer.html is missing; test_resources will fail" >&2
fi

step "lint (ruff format --check + ruff check)"
ruff format --check src tests
ruff check src tests

step "typecheck (mypy)"
python -m mypy src

step "tests + coverage (85% CI bar)"
python -m pytest tests/ --ignore=tests/e2e -q --cov=bamcp --cov-fail-under=85 --cov-report=term:skip-covered

step "eval-smoke (deterministic ship-gate)"
python -m bamcp.eval.truthset --manifest tests/eval/datasets/synthetic_v1/manifest.yaml
python -m bamcp.eval.truthset --manifest tests/eval/datasets/indel_v1/manifest.yaml --min-precision 1.0

step "eval-visual-dry (visual-vs-text eval: fixtures, stimuli, leak-lint, grader self-tests)"
if [[ -z "${PLAYWRIGHT_CHROMIUM_EXECUTABLE:-}" && -x /opt/pw-browsers/chromium ]]; then
  export PLAYWRIGHT_CHROMIUM_EXECUTABLE=/opt/pw-browsers/chromium
fi
python -m bamcp.eval.visual dry --out .eval-results/visual-dry

if [[ "${VERIFY_E2E:-0}" == "1" ]]; then
  step "e2e (playwright, needs a built viewer + chromium)"
  if [[ -z "${PLAYWRIGHT_CHROMIUM_EXECUTABLE:-}" && -x /opt/pw-browsers/chromium ]]; then
    export PLAYWRIGHT_CHROMIUM_EXECUTABLE=/opt/pw-browsers/chromium
  fi
  python -m pytest tests/e2e/ -q --no-cov
fi

step "OK — all local gates green"
