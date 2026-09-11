#!/usr/bin/env python3
"""Run a bounded real worker continuity check on a 1UA2-derived TPO fragment.

This intentionally truncated 13-residue fixture tests software preservation. It
is not a prepared full 1UA2 protein, an equilibrated trajectory, or a model of
1UA2 biological function. Its source provenance does not inherit full SEQRES.
"""
from __future__ import annotations

import argparse
import csv
import datetime as dt
import hashlib
import json
import os
from pathlib import Path
import sys
import time

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def run():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--data-dir', type=Path, default=ROOT / 'data/integration-checks/modified-residue-native')
    parser.add_argument('--source', type=Path, default=ROOT / 'data/datasets/f491cb784d7f4fbf/topology.pdb')
    parser.add_argument('--report', type=Path, default=ROOT / 'docs/audit/modified-residue-native-run.json')
    args = parser.parse_args()
    os.environ['DYNAMOL_DATA_DIR'] = str(args.data_dir.resolve())
    os.environ['DYNAMOL_CPU_THREADS'] = '2'
    import mdtraj as md
    import numpy as np
    import openmm as mm
    from openmm import app, unit
    from backend import config, jobs, preparation, storage
    from backend.complex_topology import residue_key, subset
    from backend.modified_residues import register_topology_definitions, modified_stereochemistry_report
    from backend.models import SimulationConfig
    from backend.prepared_system import load_prepared_forcefield

    report = {'passed': False, 'started_at': dt.datetime.now(dt.timezone.utc).isoformat(),
              'scope': 'Real preparation, solvation and 2-step OpenMM software continuity check on an explicitly truncated 13-residue TPO-containing 1UA2 fragment. This is not full 1UA2 preparation or scientific convergence validation.',
              'source': {'path': str(args.source.resolve()), 'sha256': sha(args.source), 'selection': 'chain A, residues 164–176 inclusive; full-source sequence deliberately not inherited'},
              'data_dir': str(config.DATA_ROOT), 'jobs': {}, 'checks': {}, 'versions': {'openmm': mm.__version__, 'mdtraj': md.__version__}}
    args.report.parent.mkdir(parents=True, exist_ok=True)
    def save():
        args.report.write_text(json.dumps(report, indent=2, allow_nan=False) + '\n')
    def wait(job, operation, timeout=360):
        report['jobs'][operation] = {'id': job['id']}
        started = time.monotonic()
        last = None
        while True:
            status = jobs.get_job(job['id'])
            if status['stage'] != last:
                print(f"{operation}: {status['status']} · {status['stage']}", flush=True)
                last = status['stage']
            if status['status'] in {'completed', 'failed', 'cancelled', 'interrupted'}:
                report['jobs'][operation].update(status=status['status'], elapsed_seconds=status['elapsed_seconds'], dataset_id=status.get('dataset_id'), error=status.get('error'))
                save()
                if status['status'] != 'completed':
                    raise RuntimeError(f"{operation} failed: {status.get('error')}")
                return storage.get_dataset(status['dataset_id']), config.JOBS_DIR / job['id']
            if time.monotonic() - started > timeout:
                jobs.cancel_job(job['id'])
                raise RuntimeError(f'{operation} exceeded the bounded {timeout}-second check; cancellation requested.')
            time.sleep(.25)
    def atom_key(atom):
        return (*residue_key(atom.residue), atom.name)
    def pdb_for(meta):
        register_topology_definitions()
        return app.PDBFile(str(storage.dataset_dir(meta['id']) / 'prepared.pdb'))
    def tpo_charge(structure, system):
        atoms = [atom for atom in structure.topology.atoms() if atom.residue.name == 'TPO']
        nonbonded = next(force for force in system.getForces() if isinstance(force, mm.NonbondedForce))
        return {atom.name: float(nonbonded.getParticleParameters(atom.index)[0].value_in_unit(unit.elementary_charge)) for atom in atoms}
    def charge_for(meta):
        structure = pdb_for(meta)
        forcefield, _ = load_prepared_forcefield(storage.dataset_dir(meta['id']), meta['preparation'])
        system = forcefield.createSystem(structure.topology, nonbondedMethod=app.NoCutoff, constraints=app.HBonds)
        return tpo_charge(structure, system)
    def snapshots(folder):
        return {path.name: sha(path) for path in sorted((folder / 'residue-parameters').iterdir()) if path.is_file()}
    def atom_positions(structure):
        coords = np.asarray(structure.positions.value_in_unit(unit.angstrom))
        return {atom_key(atom): coords[atom.index] for atom in structure.topology.atoms()}

    try:
        register_topology_definitions()
        original = app.PDBFile(str(args.source))
        keys = {residue_key(residue) for residue in original.topology.residues() if residue.chain.id == 'A' and 164 <= int(residue.id) <= 176}
        assert len(keys) == 13
        fragment = subset(original.topology, original.positions, keys)
        fixture_path = config.DATA_ROOT / 'fixture-A164-176.pdb'
        with fixture_path.open('w') as handle:
            app.PDBFile.writeFile(fragment.topology, fragment.positions, handle, keepIds=True)
        traj = md.Trajectory(np.asarray(fragment.positions.value_in_unit(unit.nanometer))[None], md.Topology.from_openmm(fragment.topology))
        fixture = storage.save_dataset(traj, 'Audit fixture · 1UA2 A164–176 (truncated)', 'Analytical continuity fixture', report['scope'], provenance={'audit_source': report['source'], 'no_full_source_sequence_inheritance': True})
        report['fixture'] = {'dataset_id': fixture['id'], 'residues': fixture['n_residues'], 'atoms': fixture['n_atoms'], 'pdb_sha256': sha(fixture_path)}
        original_tpo_atoms = [atom for atom in fragment.topology.atoms() if atom.residue.name == 'TPO' and atom.element != app.element.hydrogen]
        original_positions = atom_positions(fragment)
        tpo_keys = {atom_key(atom) for atom in original_tpo_atoms}
        report['fixture']['tpo_heavy_atoms'] = len(tpo_keys)
        assert len(tpo_keys) == 11
        inspection = preparation.inspect_preparation(fixture['id'])
        assert len(inspection['modified_residues']) == 1 and not inspection['ligands'] and not inspection['blockers'] and not inspection['gaps']
        assert inspection['missing_residues'] == []
        report['checks']['isolated_fragment_has_no_inherited_full_sequence'] = not inspection['has_sequence'] and not inspection['missing_residues']
        save()
        prepared, prep_job = wait(preparation.submit_preparation({'dataset_id': fixture['id'], 'name': 'Audit TPO fragment · prepared', 'ph': 7, 'remove_heterogens': True, 'optimize_sidechains': True, 'seed': 2026}), 'preparation')
        prepared_pdb = pdb_for(prepared)
        final_positions = atom_positions(prepared_pdb)
        assert tpo_keys.issubset(final_positions)
        delta = max(float(np.linalg.norm(final_positions[key] - original_positions[key])) for key in tpo_keys)
        assert delta <= .0018
        report['checks']['retained_all_TPO_heavy_atoms_despite_remove_heterogens'] = True
        report['checks']['TPO_heavy_coordinate_max_displacement_angstrom'] = delta
        tpo = next(residue for residue in prepared_pdb.topology.residues() if residue.name == 'TPO')
        tpo_bonds = [(atom_key(a), atom_key(b)) for a, b in prepared_pdb.topology.bonds() if a.residue == tpo or b.residue == tpo]
        assert sum((a[:4] == residue_key(tpo)) != (b[:4] == residue_key(tpo)) for a, b in tpo_bonds) == 2
        assert sum((a[-1] == 'P' or b[-1] == 'P') for a, b in tpo_bonds) == 4
        report['TPO_bonds_after_preparation'] = tpo_bonds
        report['checks']['TPO_peptide_links_and_four_phosphate_bonds_retained'] = True
        stereo = modified_stereochemistry_report(prepared_pdb.topology, prepared_pdb.positions)
        assert stereo['checked_centers'] == 2 and not stereo['violations']
        report['modified_stereochemistry'] = stereo
        hydrogens = [atom.name for atom in tpo.atoms() if atom.element == app.element.hydrogen]
        assert hydrogens
        report['TPO_added_hydrogens'] = hydrogens
        reference_charges = charge_for(prepared)
        assert abs(sum(reference_charges.values()) + 2) < 1e-5
        report['TPO_per_atom_charges_e'] = reference_charges
        report['checks']['TPO_total_charge_e'] = sum(reference_charges.values())
        reference_snapshot = snapshots(prep_job)
        assert reference_snapshot == snapshots(storage.dataset_dir(prepared['id']))
        report['residue_parameter_snapshot_sha256'] = reference_snapshot
        save()
        solvated, solvent_job = wait(preparation.submit_solvation(prepared['id'], {'padding_nm': 1, 'ph': 7, 'seed': 2026}), 'solvation')
        solvated_pdb = pdb_for(solvated)
        solvent_positions = atom_positions(solvated_pdb)
        # Every prepared solute atom, including H and terminal groups, stays at its selected position.
        assert all(key in solvent_positions and np.allclose(position, solvent_positions[key], atol=1e-8, rtol=0) for key, position in final_positions.items())
        report['checks']['solvation_preserves_all_prepared_solute_coordinates'] = True
        report['solvated_atoms'] = solvated['n_atoms']
        assert charge_for(solvated) == reference_charges
        for folder in [solvent_job, storage.dataset_dir(solvated['id'])]:
            assert snapshots(folder) == reference_snapshot
        save()
        trajectory, md_job = wait(jobs.submit_job(SimulationConfig(dataset_id=solvated['id'], engine='openmm', name='Audit TPO fragment · 2-step MD', duration_ps=.004, timestep_fs=2, report_interval=1, solvent='explicit', equilibration_steps=2, minimize=True, seed=2026)), 'dynamics')
        system = mm.XmlSerializer.deserialize((md_job / 'system.xml').read_text())
        run_pdb = app.PDBFile(str(md_job / 'prepared.pdb'))
        assert tpo_charge(run_pdb, system) == reference_charges
        assert charge_for(trajectory) == reference_charges
        for folder in [md_job, storage.dataset_dir(trajectory['id'])]:
            assert snapshots(folder) == reference_snapshot
        assert trajectory['preparation']['modified_residue_parameters'] == prepared['preparation']['modified_residue_parameters']
        physical = storage.load_physical(trajectory['id'])
        assert np.isfinite(physical.xyz).all() and physical.n_frames == 3
        with (md_job / 'energies.csv').open() as handle:
            energies = [{key: float(value) for key, value in row.items()} for row in csv.DictReader(handle)]
        assert np.isfinite([[row['potential_kj_mol'], row['kinetic_kj_mol']] for row in energies]).all()
        report['checks'].update(exact_per_atom_charges_and_parameter_snapshot_across_all_stages=True, finite_native_coordinates_and_energies=True, original_full_source_unchanged=sha(args.source) == report['source']['sha256'])
        report['trajectory'] = {'dataset_id': trajectory['id'], 'atoms': trajectory['n_atoms'], 'frames': trajectory['n_frames'], 'times_ps': trajectory['times_ps'], 'energies': energies}
        report['passed'] = True
        print('PASS: real TPO preparation → solvation → 2-step OpenMM continuity.', flush=True)
    except Exception as exc:
        report['error'] = f'{type(exc).__name__}: {exc}'
        raise
    finally:
        report['finished_at'] = dt.datetime.now(dt.timezone.utc).isoformat()
        save()


if __name__ == '__main__':
    run()
