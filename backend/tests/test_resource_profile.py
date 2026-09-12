"""Explicit atom-cap opt-in retains the independent trajectory-memory bound."""
import json
import os
from pathlib import Path
import subprocess
import sys

import pytest


ROOT = Path(__file__).resolve().parents[2]


@pytest.mark.parametrize('setting,expected', [(None, 100000), ('200000', 200000), ('125000', 125000)])
def test_bounded_atom_profile_is_shared_with_selection_validation(tmp_path, setting, expected):
    env = {**os.environ, 'DYNAMOL_DATA_DIR': str(tmp_path)}
    env.pop('DYNAMOL_MAX_ATOMS', None)
    if setting is not None:
        env['DYNAMOL_MAX_ATOMS'] = setting
    script = "from backend import config, workspaces; from backend.structural_analysis import StructuralAnalysisRequest; import json; print(json.dumps([config.MAX_ATOMS, config.MAX_COORD_BYTES, len(workspaces._indices(list(range(config.MAX_ATOMS)), config.MAX_ATOMS)), StructuralAnalysisRequest.model_json_schema()['properties']['atoms']['maxItems']]))"
    result = subprocess.run([sys.executable, '-c', script], cwd=ROOT, env=env, capture_output=True, text=True, timeout=30)
    assert result.returncode == 0, result.stderr
    assert json.loads(result.stdout) == [expected, 256 * 1024**2, expected, expected]


@pytest.mark.parametrize('setting', ['0', '200001', 'invalid'])
def test_invalid_or_unbounded_profile_cannot_start(tmp_path, setting):
    env = {**os.environ, 'DYNAMOL_DATA_DIR': str(tmp_path), 'DYNAMOL_MAX_ATOMS': setting}
    result = subprocess.run([sys.executable, '-c', 'from backend import config'], cwd=ROOT, env=env, capture_output=True, text=True, timeout=10)
    assert result.returncode != 0
