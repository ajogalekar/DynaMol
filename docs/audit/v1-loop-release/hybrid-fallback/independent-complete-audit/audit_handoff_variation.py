"""Read-only runner artifact contract and repeat-output comparison; no dynamics."""
from pathlib import Path
import hashlib
import json
import numpy as np
from openmm import app

BASE = Path('/Users/ashujo/.cache/dynamol-research/v1-loop-release')
RUN = BASE / 'accepted-coordinate-handoff-v1'
FIRST = BASE / 'outer-torsion-anchors-v1/04-conditioned/validation'
SECOND = RUN / 'search/torsion-source-ccd-2027/validation'
DATASET = Path('/Users/ashujo/Documents/Science/DynaMol/docs/audit/loop-fallback/runs/root-final/1UA2/workspace/datasets/4dcb5c25f5c14559')
CHECKS = frozenset(('identity', 'retained_environment', 'geometry', 'stereochemistry',
                    'backbone_reference', 'observed_displacement', 'serialized_output', 'parameter_integrity'))

def load(path):
    return json.loads(path.read_text())

def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()

plan = load(RUN/'plan.json')
search = load(RUN/'search/search.json')
summary = load(RUN/'summary.json')
progress = load(RUN/'progress.json')
frozen = {**plan['frozen_inputs'], **plan['implementation_sha256']}
for p in [RUN/'plan.json', RUN/'summary.json', RUN/'search/search.json', RUN/'progress.json']:
    frozen[str(p)] = sha(p)
assert all(sha(Path(p)) == digest for p, digest in frozen.items())
assert search['status'] == summary['status'] == 'accepted'
assert len(search['attempts']) == 5
assert [r['status'] for r in search['attempts']] == ['rejected']*4 + ['accepted']
for i, row in enumerate(search['attempts']):
    result_path = RUN/'search'/row['name']/'validation/result.json'
    result = load(result_path)
    assert set(row['checks']) == CHECKS and row['checks'] == result['checks']
    assert all(type(v) is bool for v in row['checks'].values())
    if i < 4:
        assert row['checks']['backbone_reference'] is False and 'artifact' not in row
    else:
        assert all(row['checks'][k] is True for k in CHECKS)
    for p, digest in {**result['source_sha256'], **result['artifacts_sha256']}.items():
        assert sha(Path(p)) == digest
        frozen[p] = digest
    frozen[str(result_path)] = sha(result_path)
    assert [p['stage'] for p in progress if p['attempt'] == row['name']] == ['generating', 'validating', 'complete']
artifact = search['accepted_artifact']
assert artifact == summary['accepted_artifact'] == search['attempts'][-1]['artifact']
path = Path(artifact['path']).resolve()
assert path == (SECOND/'refined.npz').resolve()
assert path.is_relative_to((RUN/'search/torsion-source-ccd-2027').resolve())
assert sha(path) == artifact['sha256'] == load(SECOND/'result.json')['artifacts_sha256'][str(path)]
assert len(progress) == 15

atoms = list(app.PDBxFile(str(DATASET/'loop-refinement-topology.cif')).topology.atoms())
groups = {'all': list(range(len(atoms))),
          'heavy': [a.index for a in atoms if a.element != app.element.hydrogen],
          'hydrogen': [a.index for a in atoms if a.element == app.element.hydrogen]}
comparison = {}
for filename in ['candidate-seed.npz', 'refined.npz']:
    one, two = np.load(FIRST/filename)['xyz_nm'], np.load(SECOND/filename)['xyz_nm']
    assert one.shape == two.shape == (len(atoms), 3)
    delta = np.linalg.norm(one-two, axis=1)*10
    comparison[filename] = {}
    for group, ix in groups.items():
        maximum = max(ix, key=lambda i: delta[i])
        a = atoms[maximum]
        comparison[filename][group] = {
            'coordinates_bitwise_equal': bool(np.array_equal(one[ix], two[ix])),
            'maximum_difference_A': float(delta[maximum]),
            'rms_difference_A': float(np.sqrt(np.mean(delta[ix]**2))),
            'maximum_atom': [a.residue.chain.id, a.residue.id, a.residue.name, a.name],
        }
    frozen[str(FIRST/filename)], frozen[str(SECOND/filename)] = sha(FIRST/filename), sha(SECOND/filename)
assert all(sha(Path(p)) == digest for p, digest in frozen.items())
report = {
    'verified_hashes': frozen,
    'implementation_sha256': sha(Path(__file__)),
    'accepted_artifact': artifact,
    'five_attempts_four_rejections_then_acceptance_verified': True,
    'actual_refined_npz_path_hash_and_attempt_directory_verified': True,
    'progress_event_count': len(progress),
    'coordinate_repeat_comparison': comparison,
    'endpoint_independently_audited_by': str(Path(__file__).with_name('handoff-audit-summary.json')),
    'scope': 'Known-outcome orchestration regression and separate saved-artifact audit; not a success-rate benchmark, bitwise reproducibility claim, full preparation run or release qualification.',
    'new_refinement_md_or_qm': False,
}
for path in [Path(__file__).with_name('handoff-variation.json'), BASE/'handoff-independent-v1/handoff-variation.json']:
    assert not path.exists()
    path.write_text(json.dumps(report, indent=2)+'\n')
print(json.dumps({'accepted_artifact': artifact, 'comparison': comparison}, indent=2))
