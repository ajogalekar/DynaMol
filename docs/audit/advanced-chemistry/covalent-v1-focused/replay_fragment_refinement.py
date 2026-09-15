"""Paired local refinement from immutable saved seeds; no native MD run."""
import argparse
import json
import os
from pathlib import Path
import sys
import time

import numpy as np
from openmm import app, unit

from run_fragment_comparison import bounded, digest, write


def coordinates(path):
    structure = app.PDBFile(str(path))
    keys = [(a.residue.chain.id, a.residue.id, (a.residue.insertionCode or '').strip(),
             a.residue.name, a.name, a.element.symbol) for a in structure.topology.atoms()]
    if len(set(keys)) != len(keys):
        raise ValueError('Ambiguous saved seed atoms')
    return keys, np.asarray(structure.positions.value_in_unit(unit.nanometer))


def run(base, prior, output, repo, rama_runtime):
    output.mkdir(exist_ok=False)
    scripts = Path(__file__).resolve().parent
    original = json.loads((prior / 'result.json').read_text())
    trials = [r['name'] for r in original['trials'] if r['completed']]
    if len(trials) != original['trial_count'] or not 1 <= len(trials) <= 10:
        raise ValueError('All bounded prior candidates must be available for paired replay')
    sources = [Path(__file__), scripts / 'run_fragment_comparison.py', scripts / 'check_context_seed.py',
               scripts / 'rama_reference.py', scripts / 'check_loop_environment.py', prior / 'result.json',
               base / '6jxt-protein-v4/before-refinement.pdb', base / '6jxt-protein-v4/retained-environment.pdb',
               repo / 'backend/loop_refinement.py', repo / 'backend/loop_geometry.py',
               repo / 'backend/preparation_worker.py']
    for name in trials:
        sources.extend(prior / name / p for p in ['context/output.json', 'context/input.json', 'refinement/coherent-context-seed.pdb'])
    hashes = {str(p): digest(p) for p in sources}
    write(output / 'plan.json', {'trials': trials, 'protocols': ['replay-original-coverage', 'preserve-context-peptides'],
          'source_sha256': hashes, 'app_ready': False, 'physical_model_validated': False,
          'scope': 'Identical saved coordinates and hydrogens per pair. Add only temporary observed-context peptide coverage, preserving source isomer basins and the existing strength schedule. All admission limits unchanged.'})
    snapshot = output / 'implementation-snapshot'
    snapshot.mkdir()
    for p in sources:
        if p.suffix == '.py':
            (snapshot / p.name).write_bytes(p.read_bytes())
    env = os.environ.copy()
    env.update(PYTHONPATH=str(rama_runtime), OMP_NUM_THREADS='2', OPENBLAS_NUM_THREADS='2',
               DYNAMOL_CPU_THREADS='2')
    records = []
    started = time.monotonic()
    for name in trials:
        for mode in ['replay-original-coverage', 'preserve-context-peptides']:
            if time.monotonic() - started > 1200:
                raise RuntimeError('Overall comparison budget exhausted')
            folder = output / name / mode
            folder.mkdir(parents=True)
            seed = prior / name / 'refinement/coherent-context-seed.pdb'
            write(output / 'progress.json', {'active': [name, mode], 'completed': records})
            try:
                command = [sys.executable, str(scripts / 'check_context_seed.py'), '--repo', str(repo),
                    '--original', str(base / '6jxt-protein-v4'), '--context', str(prior / name / 'context'),
                    '--output', str(folder / 'refinement'), '--environment', str(base / '6jxt-protein-v4/retained-environment.pdb'),
                    '--backbone-construction', '--preserve-backbone-conformation', '--seed-pdb', str(seed)]
                if mode == 'preserve-context-peptides':
                    command.append('--preserve-context-peptides')
                bounded(command, folder, 'refinement', env)
                input_keys, input_xyz = coordinates(seed)
                saved_keys, saved_xyz = coordinates(folder / 'refinement/coherent-context-seed.pdb')
                if saved_keys != input_keys or not np.array_equal(saved_xyz, input_xyz):
                    raise ValueError('Refinement replay changed initial coordinates or hydrogens')
                local = json.loads((folder / 'refinement/result.json').read_text())
                bounded([sys.executable, str(scripts / 'rama_reference.py'),
                    '--pdb', str(folder / 'refinement/refined-candidate.pdb'),
                    '--selection-result', str(folder / 'refinement/result.json'),
                    '--output', str(folder / 'refined-rama')], folder, 'refined-rama', env)
                rama = json.loads((folder / 'refined-rama/result.json').read_text())
                geometry = json.loads((folder / 'refinement/geometry.json').read_text())
                refinement = json.loads((folder / 'refinement/refinement.json').read_text())
                targets = refinement['peptide_construction_restraints']['targets']
                covered = {tuple(r['atoms']) for r in targets}
                passing = bool(local['passed_local_checks'] and rama['all_selected_scored_without_outliers'])
                environment_passed = False
                if passing:
                    bounded([sys.executable, str(scripts / 'check_loop_environment.py'), '--repo', str(repo),
                        '--original', str(base / '6jxt-protein-v4'), '--candidate', str(folder / 'refinement'),
                        '--environment', str(base / '6jxt-protein-v4/retained-environment.pdb'),
                        '--output', str(folder / 'environment-audit')], folder, 'environment-audit', env)
                    environment_passed = json.loads((folder / 'environment-audit/result.json').read_text())['passed']
                record = {'trial': name, 'mode': mode, 'completed': True,
                    'identical_saved_seed_including_hydrogens': True,
                    'local_and_saved_reference_passed': passing, 'saved_environment_audit_passed': environment_passed,
                    'result': str(folder / 'refinement/result.json'),
                    'geometry_accepted': local['geometry_accepted'], 'stereo_error': local['stereo_error'],
                    'maximum_observed_context_displacement_nm': local['maximum_observed_context_displacement_nm'],
                    'outliers': [r['residue'] for r in rama['outliers']],
                    'temporary_peptide_target_count': len(targets),
                    'observed_context_targets': [r for r in targets if not r['touches_modeled_residue']],
                    'failed_peptides': [{'residues': r['residues'], 'deviation_degrees': r['deviation_degrees'],
                        'covered_by_temporary_peptide_restraint': tuple(r['atoms']) in covered}
                        for r in geometry['omega_checks'] if not r['accepted']],
                    'app_ready': False, 'physical_model_validated': False}
            except Exception as exc:
                record = {'trial': name, 'mode': mode, 'completed': False, 'error': str(exc),
                          'local_and_saved_reference_passed': False, 'saved_environment_audit_passed': False,
                          'app_ready': False, 'physical_model_validated': False}
            records.append(record)
            write(output / 'progress.json', {'active': None, 'completed': records})
            print(json.dumps({k: v for k, v in record.items() if k != 'observed_context_targets'}), flush=True)
    if hashes != {name: digest(Path(name)) for name in hashes}:
        raise ValueError('A paired comparison source changed during execution')
    report = {'stage': 'paired_saved_fragment_refinement_complete', 'results': records,
              'trial_count': len(trials), 'refinement_count': len(records),
              'all_completed': all(r['completed'] for r in records), 'source_sha256': hashes,
              'passing_candidates': [r['result'] for r in records if r['saved_environment_audit_passed']],
              'elapsed_seconds': time.monotonic() - started,
              'app_ready': False, 'physical_model_validated': False,
              'scope': 'Local construction/PDB/environment screen. Any passing research candidate still requires independent native assembly and dynamics; no force-field accuracy or native loop claim.'}
    write(output / 'result.json', report)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    for name in ['base', 'prior', 'output', 'repo', 'rama-runtime']:
        parser.add_argument('--' + name, type=Path, required=True)
    args = parser.parse_args()
    run(args.base, args.prior, args.output, args.repo, args.rama_runtime)
