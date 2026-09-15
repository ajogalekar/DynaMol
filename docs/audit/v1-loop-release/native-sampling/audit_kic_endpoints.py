"""Read-only audit of saved KIC coordinates; does not call closure or sampling."""
import hashlib
import json
from pathlib import Path

import numpy as np


ROOT = Path(__file__).resolve().parents[4]
CACHE = Path.home() / '.cache/dynamol-research/v1-loop-release'
OUTPUT = CACHE / 'kic-endpoint-audit-v1'
SOURCE = ROOT / 'docs/audit/loop-fallback/runs/root-final/8K5R/workspace/datasets/8899bea3019d447b/topology.pdb'
CONTROL = ROOT / 'docs/audit/preparation-fixtures/six_residues_intact.pdb'
RUNS = ('8k5r-native-sampling-v1', '8k5r-coherent-native-sampling-v1', '8k5r-fragment-pool-v1')


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def pdb_atoms(path):
    result = {}
    for line in path.read_text().splitlines():
        if line[:6].strip() not in ('ATOM', 'HETATM') or line[21] != 'A':
            continue
        if line[12:16].strip() not in ('N', 'CA', 'C', 'O'):
            continue
        key = ('A', line[22:26].strip(), line[26].strip(), line[17:20].strip(), line[12:16].strip())
        if key in result:
            raise ValueError(f'Duplicate exact identity in {path}: {key}')
        result[key] = np.array([float(line[30:38]), float(line[38:46]), float(line[46:54])])
    return result


def saved_atoms(rows):
    result = {tuple(a['identity']): np.array(a['xyz_angstrom']) for a in rows}
    if len(result) != len(rows):
        raise ValueError('Duplicate saved atom identity')
    return result


def rigid_fit(moving, target):
    # Row-vector Kabsch fit, independent of OST/ProMod3.
    left, _, right = np.linalg.svd((moving-moving.mean(0)).T @ (target-target.mean(0)))
    sign = np.eye(3)
    sign[2, 2] = np.linalg.det(left @ right)
    rotation = left @ sign @ right
    return rotation, target.mean(0) - moving.mean(0) @ rotation


def norm(v):
    return float(np.linalg.norm(v))


def main():
    paths = [SOURCE, CONTROL, Path(__file__)] + [CACHE/r/'native/result.json' for r in RUNS]
    paths += sorted(OUTPUT.glob('source-release*.cc'))
    runtime = Path.home() / 'Library/Application Support/DynaMol/Runtimes/promod3-dfd2a559551cd2f2'
    paths += [runtime/'include/promod3/loop/backbone.hh', runtime/'include/promod3/modelling/kic.hh']
    hashes = {str(p): sha(p) for p in paths}
    cases = []
    checks = []
    for run in RUNS:
        data = json.loads((CACHE/run/'native/result.json').read_text())
        for case in data['cases']:
            source_path = CONTROL if case['name'] == 'intact-control' else SOURCE
            source = pdb_atoms(source_path)
            saved_hashes = data['source_sha256']
            expected = next(v for k, v in saved_hashes.items() if Path(k).name == source_path.name)
            assert sha(source_path) == expected
            start, end = case['start'], case['end']
            summary = {'run': run, 'case': case['name'], 'source': str(source_path),
                       'start': start, 'end': end, 'solution_count': 0,
                       'max_n_endpoint_N_CA_C_A': 0., 'max_c_endpoint_N_CA_C_A': 0.,
                       'max_endpoint_O_A': 0., 'max_observed_internal_A': 0.,
                       'endpoint_identities': [list(k) for k in source if int(k[1]) in (start, end)]}
            for ti, trial in enumerate(case['trials']):
                if not trial.get('solutions'):
                    continue
                incoming = saved_atoms(trial['input_backbone_atoms'])
                for si, solution in enumerate(trial['solutions']):
                    actual = saved_atoms(solution['backbone_atoms'])
                    assert set(actual) == set(incoming)
                    assert len(actual) == 4*(end-start+1)
                    recorded = {tuple(r['identity']): r['displacement_angstrom'] for r in solution['observed_backbone_displacements']}
                    observed = {k for k in actual if int(k[1]) not in case['modeled']}
                    assert observed == set(recorded) and observed <= set(source)
                    residuals = {k: norm(actual[k]-source[k]) for k in observed}
                    row = {'run': run, 'case': case['name'], 'trial': ti, 'solution': si,
                           'maximum_saved_displacement_recalculation_error_A': max(abs(residuals[k]-recorded[k]) for k in observed),
                           'endpoint_coordinate_residuals_A': {':'.join(k): v for k, v in residuals.items() if int(k[1]) in (start, end)}}
                    for label, index in (('n', start), ('c', end)):
                        keys = [next(k for k in actual if int(k[1]) == index and k[4] == a) for a in ('N', 'CA', 'C')]
                        rotation, shift = rigid_fit(np.array([incoming[k] for k in keys]), np.array([source[k] for k in keys]))
                        row[label+'_endpoint_fit_prediction_error_A'] = max(norm(incoming[k] @ rotation + shift - actual[k]) for k in keys)
                        summary['max_'+label+'_endpoint_N_CA_C_A'] = max(summary['max_'+label+'_endpoint_N_CA_C_A'], *(residuals[k] for k in keys))
                        oxygen = keys[0][:-1] + ('O',)
                        summary['max_endpoint_O_A'] = max(summary['max_endpoint_O_A'], residuals[oxygen])
                        next_n = next((k for k in source if int(k[1]) == end+1 and k[4] == 'N'), None)
                        if label == 'c' and next_n is not None:
                            carbon, ca = actual[keys[2]], actual[keys[1]]
                            v1, v2 = carbon-ca, carbon-source[next_n]
                            direction = v1/norm(v1) + v2/norm(v2)
                            predicted_o = carbon + 1.230 * direction/norm(direction)
                        else:
                            predicted_o = incoming[oxygen] @ rotation + shift
                        row[label+'_oxygen_prediction_error_A'] = norm(predicted_o-actual[oxygen])
                    summary['max_observed_internal_A'] = max(summary['max_observed_internal_A'], *(v for k, v in residuals.items() if int(k[1]) not in (start, end)), 0.)
                    summary['solution_count'] += 1
                    checks.append(row)
            cases.append(summary)
    metrics = {k: max(row[k] for row in checks) for k in (
        'maximum_saved_displacement_recalculation_error_A', 'n_endpoint_fit_prediction_error_A',
        'c_endpoint_fit_prediction_error_A', 'n_oxygen_prediction_error_A', 'c_oxygen_prediction_error_A')}
    assert all(v < 0.0001 for v in metrics.values()), metrics
    assert hashes == {str(p): sha(p) for p in paths}
    result = {'scope': 'Read-only saved-coordinate/native-implementation audit; no new closure, sampling, refinement, MD, or acceptance.',
              'source_sha256': hashes, 'coordinate_units': 'angstrom', 'case_summaries': cases,
              'solution_count_including_repeated_controls': len(checks), 'metrics': metrics,
              'solution_checks': checks, 'identity_checks_passed': True,
              'physical_model_validated': False, 'app_ready': False}
    OUTPUT.mkdir(exist_ok=True)
    (OUTPUT/'result.json').write_text(json.dumps(result, indent=2, allow_nan=False)+'\n')
    compact = {k: v for k, v in result.items() if k != 'solution_checks'}
    compact['full_result_sha256'] = sha(OUTPUT/'result.json')
    compact['full_result_path'] = str(OUTPUT/'result.json')
    Path(__file__).with_name('kic-endpoint-audit-summary.json').write_text(json.dumps(compact, indent=2)+'\n')
    print(json.dumps({'metrics': metrics, 'solutions': len(checks), 'output': str(OUTPUT/'result.json')}))


if __name__ == '__main__':
    main()
