"""Materialize, but never launch, a full-adduct reference after parent gates pass.

Run using DynaMol's .venv (RDKit is required for the parent stereochemistry check).
The output is an ordinary, fully specified input for the isolated reference worker.
"""
from __future__ import annotations
import argparse
import hashlib
import importlib.util
import json
from pathlib import Path
import numpy as np


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def materialize(plan_path, output_path):
    if output_path.exists():
        raise ValueError('Existing materialized requests are retained, not overwritten')
    plan = json.loads(plan_path.read_text())
    if plan['status'] != 'prepared_pending_parent_geometry_and_launch_review':
        raise ValueError('Expected explicit pending research-reference plan')
    provider = Path(plan['provider']['worker'])
    if sha(provider) != plan['provider']['worker_sha256']:
        raise ValueError('Qualified reference worker changed')
    qualification_path = Path(plan['qualification']['result_path'])
    if sha(qualification_path) != plan['qualification']['result_sha256']:
        raise ValueError('Qualification report changed')
    qualification = json.loads(qualification_path.read_text())
    if qualification['status'] != 'passed_numerical_qualification' or qualification['worker_sha256'] != sha(provider):
        raise ValueError('Matching passed numerical qualification is required')
    parent = plan['parent']
    input_path, output_dir = Path(parent['input_path']), Path(parent['output_dir'])
    if sha(input_path) != parent['input_sha256']:
        raise ValueError('Pinned parent request changed')
    original = json.loads(input_path.read_text())
    for name, expected_sha in parent['cap_files_sha256'].items():
        if sha(Path(parent['cap_dir'])/name) != expected_sha:
            raise ValueError('Pinned parent chemical graph or atom map changed')
    result_path = output_dir/'result.json'
    result = json.loads(result_path.read_text())
    if result.get('input_sha256') != sha(input_path):
        raise ValueError('Actual parent result is not bound to the pinned input')
    if result.get('charge') != original['charge'] or result.get('spin_2S') != original['spin']:
        raise ValueError('Actual parent result changed the electronic state')
    arrays_path = output_dir/'arrays.npz'
    if sha(arrays_path) != result['arrays_sha256']:
        raise ValueError('Parent numerical arrays do not match accepted result')
    helper_path = Path(plan['parent_geometry_validator']['path'])
    if sha(helper_path) != plan['parent_geometry_validator']['sha256']:
        raise ValueError('Parent geometry validation implementation changed')
    spec = importlib.util.spec_from_file_location('parent_geometry_validator', helper_path)
    helper = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(helper)
    with np.load(arrays_path, allow_pickle=False) as arrays:
        _, geometry_report = helper.check_geometry(Path(parent['cap_dir']), original, result, arrays)
        xyz = np.asarray(arrays['coords_bohr']).tolist()
    if original['atom_ids'] != plan['atom_ids'] or original['elements'] != plan['elements']:
        raise ValueError('Planned parent identities differ from actual optimized model')
    request = dict(plan['request_settings'], atom_ids=plan['atom_ids'], elements=plan['elements'], coords_bohr=xyz)
    request['provenance'] = {
        'research_only': True, 'purpose': 'First dispersion-corrected full-parent fixed-geometry reference; no fitting or physical model acceptance.',
        'plan_path': str(plan_path.resolve()), 'plan_sha256': sha(plan_path),
        'parent_input_path': str(input_path), 'parent_input_sha256': sha(input_path),
        'parent_result_path': str(result_path), 'parent_result_sha256': sha(result_path),
        'parent_arrays_path': str(arrays_path), 'parent_arrays_sha256': sha(arrays_path),
        'parent_geometry_checks': geometry_report, 'qualified_worker_sha256': sha(provider),
        'qualification_result_sha256': sha(qualification_path),
        'parent_geometry_method': 'DF-RHF/6-31G* with five frozen temporary cap heavy atoms',
        'limits': 'One reference geometry cannot validate a torsion profile; gas-phase capped-adduct electronic model, declared state only.'}
    output_path.write_text(json.dumps(request, indent=2)+'\n')
    return request


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('plan', type=Path)
    parser.add_argument('output', type=Path)
    args = parser.parse_args()
    materialize(args.plan, args.output)
