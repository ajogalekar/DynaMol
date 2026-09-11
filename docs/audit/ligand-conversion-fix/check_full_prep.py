"""Inspect the retained real 6DBK prep and solvent jobs, without rerunning MD.

Run from the repository root with .venv/bin/python. The local integration job
records identify the actual completed outputs; input/user data are read only.
"""
import hashlib
import json
import os
import sys
from collections import Counter
from pathlib import Path

ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT))
RUN = ROOT / 'data/integration-checks/ligand-conversion-fix'
DATA = RUN / 'full-prep'
os.environ['DYNAMOL_DATA_DIR'] = str(DATA)

from openmm import app
from backend.modified_residues import register_topology_definitions
from backend.prepared_system import load_prepared_forcefield
from backend.tests.test_ligand_identity import parameter_terms


def read(path):
    return json.loads(path.read_text())


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


register_topology_definitions()
results = []
for stage, filename in [('prepared', 'full-prep-job.json'), ('solvated', 'full-solvation-job.json')]:
    job = read(DATA / 'jobs' / read(RUN / filename)['id'] / 'status.json')
    assert job['status'] == 'completed', job
    folder = DATA / 'datasets' / job['dataset_id']
    metadata = read(folder / 'metadata.json')
    preparation = metadata['preparation']
    pdb = app.PDBFile(str(folder / 'prepared.pdb'))
    ff, _ = load_prepared_forcefield(folder, preparation)
    system = ff.createSystem(pdb.topology, nonbondedMethod=app.NoCutoff, constraints=None, rigidWater=False)
    ligands = []
    for record in preparation['ligand_parameters']['ligands']:
        key = record['key'].split(':')
        residue = next(r for r in pdb.topology.residues()
                       if [r.chain.id, r.id, (r.insertionCode or '').strip(), r.name] == key)
        mapping = {entry['prepared_name']: entry['native_index'] for entry in record['parameter_atom_map']}
        indices = [-(i + 1) for i in range(system.getNumParticles())]
        for atom in residue.atoms():
            indices[atom.index] = mapping[atom.name]
        actual = Counter({term: count for term, count in parameter_terms(system, indices).items()
                          if all(index >= 0 for index in term[1])})
        native_file = folder / 'ligands' / record['artifacts_directory'] / 'ligand.prmtop'
        native = app.AmberPrmtopFile(str(native_file)).createSystem(nonbondedMethod=app.NoCutoff, constraints=None)
        expected = parameter_terms(native)
        assert actual == expected, {'missing': str(expected - actual), 'extra': str(actual - expected)}
        ligands.append({'key': record['key'], 'atom_count': len(mapping),
                        'native_parameter_terms': sum(expected.values()), 'parameter_terms_match': True,
                        'native_prmtop_sha256': digest(native_file),
                        'conversion_validation': record['conversion_validation']})
    results.append({'stage': stage, 'job_id': job['id'], 'dataset_id': metadata['id'],
                    'elapsed_seconds': job['elapsed_seconds'], 'atoms': metadata['n_atoms'],
                    'prepared_pdb_sha256': digest(folder / 'prepared.pdb'),
                    'ligands': ligands,
                    'modified_residues': [{k: r[k] for k in ['residue', 'resid', 'template', 'assembled_charge_e']}
                                          for r in preparation['modified_residues']],
                    'retained_geometry_max_displacement_angstrom': preparation['preserved_bound_geometry']['max_displacement_angstrom']})

report = {'date': '2026-09-11', 'input': '6DBK, retained source dataset ece1af90b70c42c5',
          'method': 'Real background preparation and explicit TIP3P preview, then fresh PDB/factory reload and full ligand parameter comparison to native Amber in both assembled systems.',
          'limits': 'No dynamics or convergence claim. Parameter comparison uses six-decimal canonical terms; separate three-pose native energy/force tests retain the unchanged production tolerances.',
          'results': results, 'passed': True,
          'source_sha256': {str(path.relative_to(ROOT)): digest(path) for path in [
              ROOT / 'backend/forcefield_identity.py', ROOT / 'backend/ligands.py',
              ROOT / 'backend/prepared_system.py', ROOT / 'backend/tests/test_ligand_identity.py', Path(__file__)]}}
destination = Path(__file__).with_name('full-prep-results.json')
destination.write_text(json.dumps(report, indent=2) + '\n')
print(json.dumps({'passed': True, 'stages': [{k: r[k] for k in ['stage', 'job_id', 'atoms']} for r in results]}))
