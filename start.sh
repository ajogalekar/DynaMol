#!/usr/bin/env bash
set -euo pipefail
cd -- "$(dirname -- "${BASH_SOURCE[0]}")"
if [ "${1:-}" != "--dev" ] && [ -x .venv/bin/python ] && .venv/bin/python -c 'import urllib.request; assert b"DynaMol" in urllib.request.urlopen("http://127.0.0.1:8765", timeout=1).read()' 2>/dev/null; then
  printf '%s\n' 'DynaMol is already running → http://127.0.0.1:8765'
  if [ "${DYNAMOL_OPEN_BROWSER:-0}" = "1" ]; then
    .venv/bin/python -m webbrowser 'http://127.0.0.1:8765'
  fi
  exit 0
fi
if ! command -v uv >/dev/null 2>&1; then
  printf '%s\n' 'DynaMol needs uv: https://docs.astral.sh/uv/getting-started/installation/'
  exit 1
fi
if ! command -v npm >/dev/null 2>&1; then
  printf '%s\n' 'DynaMol needs Node.js 22 or newer: https://nodejs.org/'
  exit 1
fi
printf '\n%s\n' 'DynaMol · preparing your local workspace'
uv sync --frozen --extra dev --python 3.12
if [ ! -d frontend/node_modules ]; then
  npm --prefix frontend ci
fi
npm --prefix frontend run build
.venv/bin/python scripts/seed_demo.py
if [ "${1:-}" = "--dev" ]; then
  printf '\n%s\n' 'API: http://127.0.0.1:8765  ·  UI: http://127.0.0.1:5173'
  .venv/bin/python -m uvicorn backend.main:app --host 127.0.0.1 --port 8765 &
  dynamol_api_pid=$!
  trap 'kill "$dynamol_api_pid" 2>/dev/null || true' EXIT INT TERM
  npm --prefix frontend run dev
else
  printf '\n%s\n' 'Open DynaMol → http://127.0.0.1:8765' 'Press Ctrl+C to stop the app. Simulation workers are independent; reopen the app to follow a running job.'
  if [ "${DYNAMOL_OPEN_BROWSER:-0}" = "1" ]; then
    .venv/bin/python scripts/open_browser.py &
  fi
  exec .venv/bin/python -m uvicorn backend.main:app --host 127.0.0.1 --port 8765
fi
