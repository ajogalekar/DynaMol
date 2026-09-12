#!/usr/bin/env bash
# Source checkouts only. The complete Mac app includes this private runtime.
set -euo pipefail
cd -- "$(dirname -- "${BASH_SOURCE[0]}")/.."
manager=""
for candidate in mamba conda "$HOME/miniforge3/bin/mamba" "$HOME/miniforge3/bin/conda"; do
  if command -v "$candidate" >/dev/null 2>&1; then
    manager="$candidate"
    break
  fi
done
if [ -z "$manager" ]; then
  printf '%s\n' 'Install Miniforge first: https://github.com/conda-forge/miniforge' >&2
  exit 1
fi
operation=create
if [ -f .tools/promod3/conda-meta/history ]; then operation=install; fi
"$manager" "$operation" -y --prefix "$PWD/.tools/promod3" --override-channels \
  -c conda-forge -c bioconda 'python=3.12' 'promod3=3.6.0' \
  'openstructure=2.11.1' 'openmm=8.5.1'
"$manager" list --prefix "$PWD/.tools/promod3" --explicit > .tools/promod3-explicit.txt
"$PWD/.tools/promod3/bin/python" -I -B -c 'import promod3, ost, openmm; print("Loop runtime:", promod3.__version__, ost.__version__, openmm.__version__)'
