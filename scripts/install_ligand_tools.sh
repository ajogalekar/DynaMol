#!/usr/bin/env bash
# Install a private, reproducible CPU AmberTools runtime for ligand preparation.
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
"$manager" create -y --prefix "$PWD/.tools/ambertools" --override-channels -c conda-forge 'ambertools=24.8'
"$manager" list --prefix "$PWD/.tools/ambertools" --explicit > .tools/ambertools-explicit.txt
printf '%s\n' 'AmberTools is installed privately for DynaMol. No global MD environment was changed.'
printf '%s\n' 'Run uv sync --frozen to install the Python preparation dependencies.'
