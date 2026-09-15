"""Recompute the small gradient audit from frozen evidence; performs no QM."""
import ast
import hashlib
import json
from pathlib import Path

import numpy as np

ROOT = Path.home() / '.cache/dynamol-research/df-gradient-integration-candidate-v1'
OUT = Path.home() / '.cache/dynamol-research/reviews/df-gradient-integration-independent-v1'


def read(path):
    return json.loads(path.read_bytes())


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main():
    manifest = ROOT / 'artifact-manifest.json'
    assert sha(manifest) == '22c8a6a53eb9ac399b8a73fc8323cbede879487bb4ebeebaa1b3fb79a311ade5'
    for entry in read(manifest)['artifacts']:
        path = ROOT / entry['file']
        assert path.stat().st_size == entry['bytes'] and sha(path) == entry['sha256']

    baseline = (ROOT / 'provider/baseline_rhf.py').read_text()
    candidate = (ROOT / 'provider/candidate_rhf.py').read_text()
    added_import = 'from pinned_bridge import contract as _dynamol_tiled_contract\n'
    old_call = "lib.einsum('pij,qji->pq', rhok_oo[k], rhok_oo[l])"
    new_call = '_dynamol_tiled_contract(rhok_oo[k], rhok_oo[l])'
    assert candidate.count(added_import) == candidate.count(new_call) == 1
    restored = candidate.replace(added_import, '').replace(new_call, old_call)
    assert restored == baseline
    assert ast.dump(ast.parse(restored)) == ast.dump(ast.parse(baseline))

    plan = read(ROOT / 'predeclared-plan.json')
    results = []
    paths = sorted((ROOT / 'runs').glob('*/result.json'))
    paths += sorted((ROOT / 'parent-settings-supplement/runs').glob('*/result.json'))
    assert len(paths) == 7
    estimates = 0
    for path in paths:
        record = read(path)
        assert record['status'] == 'passed'
        assert record['accepted_chemistry_model'] is False
        assert record['large_tensor_resource_qualification'] is False
        with np.load(path.parent / 'arrays.npz', allow_pickle=False) as arrays:
            left = arrays['baseline_gradient']
            right = arrays['candidate_gradient']
            assert left.shape == right.shape == arrays['coords_bohr'].shape
            assert np.isfinite(left).all() and np.isfinite(right).all()
            difference = float(np.abs(left - right).max())
            assert difference == 0
            for name in ['mo_coeff', 'mo_energy', 'mo_occ']:
                assert hashlib.sha256(arrays[name].tobytes()).hexdigest() == record['state'][name + '_sha256']
            errors = []
            for item in record['finite_differences']:
                minus, plus = item['energies_minus_plus_hartree']
                derivative = (plus - minus) / (2 * item['step_bohr'])
                error = abs(derivative - right[item['atom'], item['axis']])
                assert derivative == item['central_difference_hartree_per_bohr']
                assert error == item['errors']['candidate'] == item['errors']['baseline']
                assert error <= plan['tolerances'][f"fd_max_abs_{record['family']}_hartree_per_bohr"]
                estimates += 1
                errors.append(float(error))
        assert len(record['contraction_calls']) == 1
        assert record['full_gradient_calls'] == 2
        assert record['self_reported_peak_rss_bytes'] < 1_000_000_000
        for trace in record['traces'].values():
            assert trace['aux_response_enabled'] is True and trace['native_unit'] == 'au'
            assert any(any(v != 0 for v in row['aux_value']) for row in trace['aux_extra_force_events'])
        results.append({'case': str(path.relative_to(ROOT)), 'gradient_max_abs_difference': difference,
                        'fd_estimates': len(errors), 'max_fd_error': max(errors, default=None),
                        'peak_rss_bytes': record['self_reported_peak_rss_bytes']})
    assert estimates == 24
    OUT.mkdir(parents=True, exist_ok=False)
    result = {'verdict': 'No blocker in the small native numerical integration',
              'manifest_sha256': sha(manifest), 'verified_artifacts': 175,
              'source_delta_text_and_ast_verified': True, 'cases': results,
              'finite_difference_estimates_independently_recomputed': estimates,
              'new_QM_runs': 0, 'large_molecule_resource_qualification': False,
              'physical_model_acceptance': False}
    (OUT / 'result.json').write_text(json.dumps(result, indent=2) + '\n')
    print(json.dumps(result, indent=2))


if __name__ == '__main__':
    main()
