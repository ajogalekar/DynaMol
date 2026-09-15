"""Real preparation-only benchmark; each entry uses an isolated workspace.

No MD, solvent generation, fabricated loop masks or deletion of failed cases.
Run separate subsets concurrently, with at most three native jobs in total.
"""
from pathlib import Path
import argparse
import datetime as dt
import hashlib
import json
import os
import subprocess
import sys
import time
import traceback

ROOT = Path(__file__).resolve().parents[3]
OUT = Path(__file__).resolve().parent


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def save(path, obj):
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_suffix('.tmp')
    temporary.write_text(json.dumps(obj, indent=2, allow_nan=False) + '\n')
    temporary.replace(path)


def sources():
    return {str(p.relative_to(ROOT)): sha(p) for p in sorted((ROOT / 'backend').glob('*.py'))}


def verify_final_refinement(target, record, pdb, source_target):
    """Check saved refinement arrays, identities and retained source chemistry."""
    import numpy as np
    from openmm import app, unit
    from backend.loop_geometry import loop_geometry_report

    def residue_key(atom):
        r = atom.residue
        return (r.chain.id, r.id, (r.insertionCode or '').strip(), r.name)

    def atom_key(atom):
        return (*residue_key(atom), atom.name, atom.element.symbol if atom.element else None)

    loop = record['loop_construction']
    artifacts = loop['refinement_artifacts']
    for name in artifacts.values():
        assert sha(target / name) == loop['refinement_artifacts_sha256'][name], f'Refinement artifact changed: {name}'
    with np.load(target / artifacts['input']) as data:
        before = data['xyz_nm'].copy()
    with np.load(target / artifacts['final']) as data:
        after = data['xyz_nm'].copy()
    topology = app.PDBxFile(str(target / artifacts['topology'])).topology
    atoms = list(topology.atoms())
    identities = [atom_key(a) for a in atoms]
    assert len(set(identities)) == len(identities), 'Refinement atom identities are not unique.'
    assert identities == [atom_key(a) for a in pdb.topology.atoms()], 'Saved atom identity/order changed after refinement.'
    bonds = lambda top: {tuple(sorted((a.index, b.index))) for a, b in top.bonds()}
    assert bonds(topology) == bonds(pdb.topology), 'Saved bond identities changed after refinement.'
    assert before.shape == after.shape == (len(atoms), 3), 'Refinement arrays do not match the saved topology.'
    assert np.isfinite(before).all() and np.isfinite(after).all(), 'Refinement arrays contain non-finite coordinates.'
    modeled = {tuple(a['identity'][:4]) for a in loop['modeler']['worker_report']['atoms']}
    parameters = json.loads((target / artifacts['parameters']).read_text())
    assert modeled == {tuple(key) for key in parameters['modeled_residue_keys']}, 'Modeled identities changed at refinement entry.'
    relaxation = loop['refinement'].get('flank_relaxation')
    flank_rows = relaxation['residues'] if relaxation else []
    flanks = {tuple(key) for key in flank_rows}
    assert len(flanks) == len(flank_rows) and not modeled & flanks, 'Flank identities duplicate or overlap missing residues.'
    assert flanks == {tuple(key) for key in loop.get('observed_flanks_remodeled', [])}, 'Recorded flank inventories disagree.'
    residue_keys = {residue_key(a) for a in atoms}
    assert modeled | flanks <= residue_keys, 'A modeled or movable flank residue is absent.'
    heavy = [a for a in atoms if a.element != app.element.hydrogen]
    fixed = [a.index for a in heavy if residue_key(a) not in modeled | flanks]
    assert fixed, 'No independent fixed heavy-atom comparison is possible.'
    assert np.array_equal(before[fixed], after[fixed]), 'Final refinement moved a heavy atom outside loops and recorded flanks.'
    flank_atoms = [a for a in heavy if residue_key(a) in flanks]
    flank_displacements = {(*residue_key(a), a.name): float(np.linalg.norm(after[a.index] - before[a.index])) for a in flank_atoms}
    maximum = max(flank_displacements.values(), default=0.0)
    if relaxation:
        assert flank_atoms, 'Flank relaxation has no identifiable heavy atoms.'
        reported = {tuple(row['atom']): row['displacement_nm'] for row in relaxation['atoms']}
        assert len(reported) == len(relaxation['atoms']) and set(reported) == set(flank_displacements), 'Reported flank atom inventory changed.'
        assert all(abs(reported[key] - value) <= 1e-12 for key, value in flank_displacements.items()), 'Flank movement differs from the refinement-entry baseline.'
        assert abs(maximum - relaxation['maximum_displacement_nm']) <= 1e-12
        assert relaxation['maximum_allowed_displacement_nm'] == .1
        assert maximum <= .1 and relaxation['within_displacement_limit'], 'A flank heavy atom exceeded the 1 Å movement limit.'
    expanded = modeled | flanks
    geometry = loop_geometry_report(topology, after * unit.nanometer, expanded)
    assert geometry['accepted'], geometry['errors']
    saved_xyz = np.asarray(pdb.positions.value_in_unit(unit.nanometer))
    rounding_error = float(np.max(np.abs(saved_xyz - after)))
    # PDB writes Å coordinates with three decimals: at most 0.00005 nm per
    # component. This tolerance applies only to serialization, never NPZ masks.
    assert rounding_error <= .00005 + 1e-12, 'Saved PDB coordinates differ from final refinement beyond PDB rounding.'
    saved_geometry = loop_geometry_report(pdb.topology, pdb.positions, expanded)
    assert saved_geometry['accepted'], saved_geometry['errors']

    source = app.PDBFile(str(source_target / 'topology.pdb'))
    source_atoms = list(source.topology.atoms())
    source_metadata = json.loads((source_target / 'metadata.json').read_text())
    assert len(source_metadata['atoms']) == len(source_atoms)
    source_xyz = np.asarray(source.positions.value_in_unit(unit.nanometer))
    destination = {key: index for index, key in enumerate(identities)}
    retained = []
    for row in source_metadata['atoms']:
        if row['category'] not in {'ligands', 'ions'} or row['element'] in {'H', 'D'}:
            continue
        atom = source_atoms[row['index']]
        key = atom_key(atom)
        assert key in destination, f'Retained ligand/ion heavy atom was lost: {key}'
        final_index = destination[key]
        delta = float(np.max(np.abs(after[final_index] - source_xyz[atom.index])))
        # The source PDB and NPZ may use arithmetically different Å→nm
        # conversions. Permit roundoff only (1e-11 Å), not coordinate movement.
        assert delta <= 1e-12, f'Retained ligand/ion heavy atom moved: {key}'
        assert np.max(np.abs(saved_xyz[final_index] - source_xyz[atom.index])) <= .00005 + 1e-12
        retained.append({'identity': list(key), 'maximum_component_delta_nm': delta})
    evidence = {'baseline': artifacts['input'], 'final': artifacts['final'], 'topology': artifacts['topology'],
                'atom_count': len(atoms), 'fixed_heavy_atoms_checked': len(fixed),
                'fixed_heavy_coordinates_bitwise_equal': True,
                'flank_residues': [list(key) for key in sorted(flanks)],
                'flank_heavy_atoms_checked': len(flank_atoms), 'maximum_flank_displacement_nm': maximum,
                'flank_displacement_limit_nm': .1,
                'flank_displacement_baseline': 'Saved coordinates entering refinement, after recorded sidechain sampling.',
                'saved_pdb_maximum_component_rounding_nm': rounding_error,
                'saved_pdb_component_rounding_tolerance_nm': .00005,
                'retained_ligand_ion_heavy_atoms_checked': len(retained),
                'retained_ligand_ion_roundoff_tolerance_nm': 1e-12,
                'retained_ligand_ion_atoms': retained, 'expanded_geometry': geometry}
    return evidence, saved_geometry


def single(case, folder, timeout):
    data = folder / 'workspace'
    if data.exists():
        raise ValueError('Use a new run label; prior benchmark workspaces are retained.')
    data.mkdir(parents=True)
    os.environ.update(DYNAMOL_DATA_DIR=str(data), DYNAMOL_CPU_THREADS='2',
                      OPENMM_CPU_THREADS='2', OPENBLAS_NUM_THREADS='2', OMP_NUM_THREADS='2')
    sys.path.insert(0, str(ROOT))
    import numpy as np
    import openmm as mm
    from openmm import app, unit
    from backend import config, jobs, monomers, preparation, storage
    from backend import sources as inputs
    from backend.prepared_system import load_prepared_forcefield

    started = time.monotonic()
    source = (ROOT / case['input_path']).resolve()
    result = {'id': case['id'], 'passed': False, 'phase': 'import', 'case': case,
              'started_at': dt.datetime.now(dt.timezone.utc).isoformat(),
              'backend_source_start': sources(), 'driver_sha256': sha(__file__),
              'input_sha256': sha(source), 'workspace': str(data), 'checks': {}}
    job = None
    report = folder / 'result.json'
    save(report, result)
    try:
        assert result['input_sha256'] == case['input_sha256'], 'Frozen input hash changed.'
        meta = inputs.import_structure(source, name=f"Prep benchmark {case['id']}",
                                      provenance={'pdb_accession': case.get('pdb_id', case['id']),
                                                  'source_url': case['source_url'],
                                                  'benchmark_input_sha256': case['input_sha256']})
        result['original_dataset_id'] = meta['id']
        result['phase'] = 'monomer selection'
        save(report, result)
        if case.get('chain_index') is not None:
            meta = monomers.create_monomer(meta['id'], monomers.MonomerRequest(
                chain_index=case['chain_index'], keep_associated_molecules=True))
        result.update(selected_dataset_id=meta['id'], monomer_selection=meta.get('monomer_selection'))
        settings = {'ph': 7.0, 'seed': 42, 'build_missing_residues': True,
                    'add_missing_atoms': True, 'optimize_sidechains': True,
                    'remove_waters': True, 'remove_heterogens': False,
                    **case.get('settings', {}), 'dataset_id': meta['id'],
                    'name': f"Prep benchmark {case['id']}"}
        result.update(phase='inspection', settings=settings)
        save(report, result)
        inspection = preparation.inspect_preparation(meta['id'], settings['ph'],
                    settings.get('ligand_overrides'), settings.get('ligand_actions'))
        result['inspection'] = inspection
        gaps = [g for g in inspection['missing_residues'] if not g['terminal']]
        result['checks']['internal_gaps_detected'] = bool(gaps)
        expected_lengths = sorted(g['length'] for g in case.get('source_gap_inventory', []))
        actual_lengths = sorted(g['count'] for g in gaps)
        result['checks']['source_gap_lengths_preserved'] = actual_lengths == expected_lengths
        assert gaps, 'Actual source gaps were lost during import/selection.'
        assert actual_lengths == expected_lengths, f'Source/inspection gap lengths disagree: {expected_lengths} vs {actual_lengths}'
        result['phase'] = 'preparation eligibility'
        save(report, result)
        job = preparation.submit_preparation(settings)
        result.update(phase='native preparation', job_id=job['id'])
        save(report, result)
        last = None
        wait_start = time.monotonic()
        while True:
            state = jobs.get_job(job['id'])
            result['job'] = {k: state.get(k) for k in ('id', 'status', 'stage', 'error', 'dataset_id', 'elapsed_seconds')}
            if state['stage'] != last:
                print(case['id'], state['status'], state['stage'], flush=True)
                last = state['stage']
                save(report, result)
            if state['status'] in {'completed', 'failed', 'cancelled', 'interrupted'}:
                break
            if time.monotonic() - wait_start > timeout:
                jobs.cancel_job(job['id'])
                raise TimeoutError(f'Preparation exceeded {timeout}s; job cancellation requested.')
            time.sleep(.5)
        assert state['status'] == 'completed', state.get('error', state['status'])
        result['phase'] = 'output verification'
        prepared = storage.get_dataset(state['dataset_id'])
        target = storage.dataset_dir(prepared['id'])
        record = prepared['preparation']
        result.update(prepared_dataset_id=prepared['id'], prepared_atoms=prepared['n_atoms'], preparation=record)
        result['checks']['prepared_dataset_saved'] = True
        assert record['simulation_ready']
        pdb = app.PDBFile(str(target / 'prepared.pdb'))
        xyz = np.asarray(pdb.positions.value_in_unit(unit.nanometer))
        assert np.isfinite(xyz).all()
        result['checks']['finite_saved_coordinates'] = True
        loop = record['loop_construction']
        assert loop['accepted'] and loop['final_geometry']['accepted']
        result['checks']['full_complex_loop_geometry_passed'] = True
        assert loop['observed_coordinates_preserved_exactly']
        before = np.load(target / 'loop-construction-before.npz')['xyz_nm']
        after = np.load(target / 'loop-construction-after.npz')['xyz_nm']
        context = app.PDBFile(str(target / 'loop-construction.pdb'))
        modeled = {tuple(a['identity'][:4]) for a in loop['modeler']['worker_report']['atoms']}
        indices = [a.index for a in context.topology.atoms()
                   if (a.residue.chain.id, a.residue.id, (a.residue.insertionCode or '').strip(), a.residue.name) not in modeled]
        assert before.shape == after.shape
        assert np.array_equal(before[indices], after[indices])
        result['checks']['independent_observed_atom_transplant_unchanged'] = True
        evidence, geometry = verify_final_refinement(target, record, pdb, storage.dataset_dir(meta['id']))
        result['final_refinement_verification'] = evidence
        result['saved_pdb_loop_geometry'] = geometry
        result['checks'].update(refinement_artifact_hashes_match=True, final_refinement_identities_and_bonds_preserved=True,
                                final_refinement_fixed_heavy_atoms_unchanged=True, final_refinement_flank_displacements_bounded=True,
                                final_refinement_expanded_geometry_passed=True, serialized_pdb_matches_final_refinement=True,
                                retained_ligand_ion_heavy_identities_and_coordinates_preserved=True)
        result['checks']['serialized_pdb_geometry_passed'] = True
        forcefield, _ = load_prepared_forcefield(target, record, solvent='explicit')
        system = forcefield.createSystem(pdb.topology, nonbondedMethod=app.NoCutoff, constraints=None)
        assert system.getNumParticles() == len(xyz)
        result['checks']['complete_forcefield_reconstructed'] = True
        integrator = mm.VerletIntegrator(.001 * unit.picosecond)
        context = mm.Context(system, integrator, mm.Platform.getPlatformByName('CPU'), {'Threads': '2'})
        context.setPositions(pdb.positions)
        static_state = context.getState(getEnergy=True, getForces=True)
        energy = static_state.getPotentialEnergy().value_in_unit(unit.kilojoule_per_mole)
        forces = np.asarray(static_state.getForces(asNumpy=True).value_in_unit(unit.kilojoule_per_mole / unit.nanometer))
        assert np.isfinite(energy)
        assert forces.shape == xyz.shape and np.isfinite(forces).all(), 'Full-system static forces are non-finite or incomplete.'
        del context, integrator
        result['static_potential_energy_kj_mol'] = energy
        result['static_force_components'] = {'maximum_absolute_kj_mol_nm': float(np.max(np.abs(forces))),
                                              'rms_kj_mol_nm': float(np.sqrt(np.mean(forces ** 2))),
                                              'scope': 'Finite derivatives at saved coordinates; no dynamics, stability or convergence claim.'}
        result['checks']['finite_static_forcefield_energy_no_dynamics'] = True
        result['checks']['finite_static_forcefield_forces_no_dynamics'] = True
        ligands = record.get('ligand_parameters', {}).get('ligands', [])
        expected_ligands = {l['key'] for l in inspection.get('ligands', []) if not l.get('removed')}
        assert {l['key'] for l in ligands} == expected_ligands, 'Retained ligand parameter inventory changed.'
        result['checks']['all_inspected_ligands_parameterized'] = True
        assert all(l['conversion_validation']['passed'] for l in ligands)
        result['checks']['ligand_conversion_validation_passed'] = True
        assert not record.get('removed_ligand_residues')
        result['checks']['no_ligand_removal'] = True
        result['passed'] = True
        result['phase'] = 'complete'
    except Exception as exc:
        result.update(error=str(exc), traceback=traceback.format_exc())
    finally:
        if job:
            state = jobs.get_job(job['id'])
            if state['status'] in {'queued', 'running', 'cancelling'}:
                jobs.cancel_job(job['id'])
                for _ in range(30):
                    if jobs.get_job(job['id'])['status'] not in {'queued', 'running', 'cancelling'}:
                        break
                    time.sleep(.5)
            native = config.JOBS_DIR / job['id']
            result['job_artifacts'] = {str(p.relative_to(data)): {'sha256': sha(p), 'bytes': p.stat().st_size}
                                       for p in sorted(native.rglob('*')) if p.is_file()}
        result.update(elapsed_seconds=round(time.monotonic() - started, 2),
                      source_unchanged=sha(source) == case['input_sha256'], backend_source_end=sources())
        result['backend_source_unchanged_during_case'] = result['backend_source_start'] == result['backend_source_end']
        save(report, result)
        print('RESULT', case['id'], result['passed'], result.get('error', ''), flush=True)
    return result['passed']


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--manifest', type=Path, default=OUT / 'benchmark-manifest.json')
    parser.add_argument('--ids', required=True)
    parser.add_argument('--run-label', required=True)
    parser.add_argument('--timeout', type=int, default=1200)
    parser.add_argument('--single', action='store_true')
    args = parser.parse_args()
    if not all(c.isalnum() or c in '-_' for c in args.run_label):
        raise SystemExit('Use a simple unique run label.')
    selected = set(args.ids.split(','))
    cases = [c for c in json.loads(args.manifest.read_text())['cases'] if c['id'] in selected]
    assert len(cases) == len(selected), 'Unknown case ID.'
    run = OUT / 'runs' / args.run_label
    if args.single:
        assert len(cases) == 1
        raise SystemExit(0 if single(cases[0], run / cases[0]['id'], args.timeout) else 1)
    statuses = []
    for case in cases:
        folder = run / case['id']
        folder.mkdir(parents=True, exist_ok=True)
        log = folder / 'driver.log'
        with log.open('a') as stream:
            command = [sys.executable, '-B', str(Path(__file__).resolve()), '--manifest', str(args.manifest.resolve()),
                       '--ids', case['id'], '--run-label', args.run_label, '--timeout', str(args.timeout), '--single']
            process = subprocess.run(command, cwd=ROOT, stdout=stream, stderr=subprocess.STDOUT)
        statuses.append({'id': case['id'], 'returncode': process.returncode, 'report': str(folder / 'result.json')})
        save(run / ('batch-' + '-'.join(args.ids.split(',')) + '.json'), statuses)
        print(json.dumps(statuses[-1]), flush=True)


if __name__ == '__main__':
    main()
