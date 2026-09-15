"""Fresh full worker preparations on three retained experimental input files.

Inputs are copied into this audit's data root. Source metadata/coordinates are
never modified. These are representative workflow checks, not pose accuracy,
metal coordination, phosphate titration, or universal chemistry validation.
"""
from pathlib import Path
from collections import Counter
import argparse
import datetime as dt
import hashlib
import json
import os
import shutil
import sys
import time
import traceback

ROOT = Path(__file__).resolve().parents[4]
OUT = Path(__file__).resolve().parent
DATA = ROOT / 'data/integration-checks/final-release-2026-09-12/chemistry'
os.environ['DYNAMOL_DATA_DIR'] = str(DATA)
os.environ['DYNAMOL_CPU_THREADS'] = '2'
sys.path.insert(0, str(ROOT))

import numpy as np
import openmm as mm
from openmm import app, unit
from backend import config, jobs, monomers, preparation, sources, storage
from backend.ions import validate_ion_system
from backend.prepared_system import load_prepared_forcefield
from backend.tests.test_ligand_identity import parameter_terms

CASES = [
    {'id': '6DBK', 'source': 'data/integration-checks/ligand-conversion-fix/full-prep/datasets/ece1af90b70c42c5/source.cif',
     'chain_index': None, 'engine': 'openmm', 'actions': {},
     'coverage': 'Observed protein, organic bound ligand with aromatic/improper and stereochemical mapping; prior native conversion failure regression.'},
    {'id': '6A93', 'source': 'data/integration-checks/chemistry-repair-ui/datasets/63837c81bba04bcf/source.cif',
     'chain_index': 0, 'engine': 'openmm', 'actions': {'F:3004::1PE': 'repair', 'H:3006::1PE': 'repair'},
     'coverage': 'One observed protein chain, Zn2+ nonbonded model, complete and two incomplete 1PE molecules with explicit local completion.'},
    {'id': '1UA2', 'source': 'data/datasets/f491cb784d7f4fbf/source.cif',
     'chain_index': 0, 'engine': 'openmm', 'actions': {},
     'coverage': 'One observed kinase chain, TPO at fixed -2, ATP, and sequence-supported internal missing loop.'},
]

def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()

def wait(job, label, timeout=600):
    begin = time.monotonic()
    last = None
    while True:
        state = jobs.get_job(job['id'])
        if state['stage'] != last:
            print(label, state['status'], state['stage'], flush=True)
            last = state['stage']
        if state['status'] in {'completed', 'failed', 'cancelled', 'interrupted'}:
            return state
        if time.monotonic() - begin > timeout:
            jobs.cancel_job(job['id'])
            raise RuntimeError(f'{label} exceeded bounded {timeout} seconds; cancellation requested.')
        time.sleep(.5)

def inspect_inputs():
    rows = []
    inputs = OUT / 'complex-inputs'
    inputs.mkdir(exist_ok=True)
    for case in CASES:
        original = ROOT / case['source']
        target = inputs / f"{case['id']}.cif"
        shutil.copy2(original, target)
        source_hash = digest(original)
        meta = sources.import_structure(target, name=f"Audit {case['id']}", provenance={
            'audit_original_input': case['source'], 'input_sha256': source_hash,
            'pdb_accession': case['id'], 'source_url': f"https://files.rcsb.org/download/{case['id']}.cif",
            'retrieval': 'Retained previously fetched RCSB file; fresh import/preparation, not a fresh network fetch.'})
        original_id = meta['id']
        if case['chain_index'] is not None:
            meta = monomers.create_monomer(meta['id'], monomers.MonomerRequest(chain_index=case['chain_index'], keep_associated_molecules=True))
        settings = {'dataset_id': meta['id'], 'name': f"Audit {case['id']} prepared", 'engine': case['engine'],
                    'ph': 7, 'seed': 2026, 'build_missing_residues': True, 'add_missing_atoms': True,
                    'optimize_sidechains': True, 'ligand_actions': case['actions']}
        inspection = preparation.inspect_preparation(meta['id'], ligand_actions=case['actions'])
        row = {**case, 'original_dataset_id': original_id, 'selected_dataset_id': meta['id'],
               'input_sha256': source_hash, 'source_unchanged': digest(original) == source_hash,
               'settings': settings, 'inspection': inspection, 'monomer_selection': meta.get('monomer_selection')}
        try:
            preparation.validate_preparation(settings)
            row['ready'] = True
        except Exception as exc:
            row.update(ready=False, error=str(exc))
        rows.append(row)
        storage.atomic_json(OUT / 'complex-inspection.json', {'created_at': dt.datetime.now(dt.timezone.utc).isoformat(), 'rows': rows})
        print(case['id'], 'ready', row['ready'], row.get('error'), flush=True)
    return rows

def execute():
    inspection = json.loads((OUT / 'complex-inspection.json').read_text())
    report = {'started_at': dt.datetime.now(dt.timezone.utc).isoformat(),
              'scope': 'Fresh preparation and explicit water preview on representative retained experimental inputs; no dynamics/convergence or universal chemistry claim.',
              'driver_sha256': digest(__file__), 'rows': []}
    for case in inspection['rows']:
        row = {k: case[k] for k in ('id', 'input_sha256', 'source', 'settings', 'coverage', 'selected_dataset_id')}
        row['passed'] = False
        try:
            assert case['ready'], case.get('error')
            start = time.monotonic()
            state = wait(preparation.submit_preparation(case['settings']), case['id'])
            row['preparation_job'] = {k: state.get(k) for k in ('id', 'status', 'stage', 'elapsed_seconds', 'dataset_id', 'error')}
            assert state['status'] == 'completed', state.get('error')
            prepared = storage.get_dataset(state['dataset_id'])
            folder = storage.dataset_dir(prepared['id'])
            pdb = app.PDBFile(str(folder / 'prepared.pdb'))
            forcefield, _ = load_prepared_forcefield(folder, prepared['preparation'], solvent='explicit')
            system = forcefield.createSystem(pdb.topology, nonbondedMethod=app.NoCutoff, constraints=None)
            ions = validate_ion_system(pdb.topology, system, prepared['preparation'].get('ions', []))
            assert np.isfinite(pdb.positions.value_in_unit(unit.nanometer)).all()
            ligands = prepared['preparation'].get('ligand_parameters', {}).get('ligands', [])
            assert all(item['conversion_validation']['passed'] for item in ligands)
            parity = []
            for ligand in ligands:
                key = ligand['key'].split(':')
                residue = next(residue for residue in pdb.topology.residues()
                               if [residue.chain.id, residue.id, (residue.insertionCode or '').strip(), residue.name] == key)
                mapping = {entry['prepared_name']: entry['native_index'] for entry in ligand['parameter_atom_map']}
                indices = [-(i + 1) for i in range(system.getNumParticles())]
                for atom in residue.atoms():
                    indices[atom.index] = mapping[atom.name]
                actual = Counter({term: count for term, count in parameter_terms(system, indices).items()
                                  if all(index >= 0 for index in term[1])})
                native_file = folder / 'ligands' / ligand['artifacts_directory'] / 'ligand.prmtop'
                native = app.AmberPrmtopFile(str(native_file)).createSystem(nonbondedMethod=app.NoCutoff, constraints=None)
                expected = parameter_terms(native)
                assert actual == expected, {'ligand': ligand['key'], 'missing': str(expected-actual), 'extra': str(actual-expected)}
                parity.append({'key': ligand['key'], 'passed': True, 'native_parameter_terms': sum(expected.values()),
                               'native_prmtop_sha256': digest(native_file)})
            row.update(prepared_atoms=prepared['n_atoms'], prepared_pdb_sha256=digest(folder / 'prepared.pdb'),
                       ligands=ligands, native_assembled_ligand_parity=parity, ions=ions, modified_residues=prepared['preparation'].get('modified_residues', []),
                       repaired_loops=prepared['preparation'].get('rebuilt_residues', prepared['preparation'].get('missing_residues', [])),
                       preparation=prepared['preparation'])
            solvation = wait(preparation.submit_solvation(prepared['id'], {'padding_nm': 1, 'ph': 7, 'seed': 2026}), case['id'] + ' water')
            row['solvation_job'] = {k: solvation.get(k) for k in ('id', 'status', 'stage', 'elapsed_seconds', 'dataset_id', 'error')}
            assert solvation['status'] == 'completed', solvation.get('error')
            solvated = storage.get_dataset(solvation['dataset_id'])
            assert solvated['n_atoms'] > prepared['n_atoms']
            row.update(solvated_atoms=solvated['n_atoms'], elapsed_seconds=round(time.monotonic()-start, 2), passed=True,
                       source_unchanged=digest(ROOT / case['source']) == case['input_sha256'])
            assert row['source_unchanged']
        except Exception as exc:
            row.update(error=str(exc), traceback=traceback.format_exc())
        report['rows'].append(row)
        report['passed'] = all(item['passed'] for item in report['rows'])
        report['finished_at'] = dt.datetime.now(dt.timezone.utc).isoformat()
        storage.atomic_json(OUT / 'complex-preparation.json', report)
        print('RESULT', row['id'], row['passed'], row.get('error'), flush=True)
    if not report['passed']:
        raise SystemExit(1)

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--run', action='store_true')
    args = parser.parse_args()
    if args.run:
        execute()
    else:
        inspect_inputs()
