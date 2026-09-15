"""Bounded paired full-complex replay of saved loop refinement, research only.

Compare the existing protocol with the previously defined temporary phi/psi
protection, preserving the exact native parameter bundle and starting H atoms.
"""
import argparse
import hashlib
import json
import os
from pathlib import Path
import sys

import numpy as np
import openmm as mm
from openmm import app, unit
from pdbfixer import PDBFixer

ROOT = Path(__file__).resolve().parents[3]
REFERENCE = ROOT / 'docs/audit/advanced-chemistry/covalent-v1-focused'
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(REFERENCE))
from backend.loop_refinement import _minimize_attempt
from backend.loop_geometry import loop_geometry_report
from backend.prepared_system import load_prepared_forcefield
from backend.preparation_worker import stereochemistry_report, require_valid_stereochemistry
from rama_reference import rama_report, residue_key
from run_fragment_comparison import bounded


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def save(path, data):
    path.write_text(json.dumps(data, indent=2, allow_nan=False) + '\n')


def atom_key(atom):
    return (*residue_key(atom.residue), atom.name, atom.element.symbol if atom.element else None)


def child(dataset, output, protocol):
    output.mkdir(exist_ok=False)
    parameters = json.loads((dataset / 'loop-refinement-parameters.json').read_text())
    topology = app.PDBxFile(str(dataset / 'loop-refinement-topology.cif')).topology
    atoms = list(topology.atoms())
    identities = [atom_key(a) for a in atoms]
    if len(set(identities)) != len(identities):
        raise ValueError('Duplicate saved atom identities')
    with np.load(dataset / 'loop-refinement-input.npz') as archive:
        initial = archive['xyz_nm'].copy()
    if initial.shape != (len(atoms), 3) or not np.isfinite(initial).all():
        raise ValueError('Invalid saved refinement input')
    with (output / 'exact-initial.pdb').open('w') as handle:
        app.PDBFile.writeFile(topology, initial * unit.nanometer, handle, keepIds=True)
    np.savez(output / 'exact-initial.npz', xyz_nm=initial)
    modeled = set(map(tuple, parameters['modeled_residue_keys']))
    affected = set(modeled)
    for a, b in topology.bonds():
        if a.residue is b.residue or {a.name, b.name} != {'C', 'N'}:
            continue
        pair = {residue_key(a.residue), residue_key(b.residue)}
        if pair & modeled:
            affected |= pair
    seed_reference = rama_report(topology, initial * unit.nanometer, affected)
    if seed_reference['unavailable']:
        raise ValueError('Every modeled and affected boundary residue must have a native reference score')
    save(output / 'seed-reference.json', seed_reference)
    forcefield, files = load_prepared_forcefield(dataset, parameters['preparation'], solvent=parameters['solvent'])
    def fresh_system():
        return forcefield.createSystem(topology, nonbondedMethod=app.NoCutoff, constraints=None, rigidWater=False)
    native_xml = mm.XmlSerializer.serialize(fresh_system())
    (output / 'native-system-before.xml').write_text(native_xml)
    targets = []
    for row in seed_reference['rows']:
        if row['classification'] == 'outlier':
            continue
        for kind in ['phi', 'psi']:
            targets.append({'residue': row['residue'], 'kind': kind, 'atoms': row[kind + '_atoms'],
                            'target_radians': float(np.deg2rad(row[kind + '_degrees'])),
                            'construction_strength_kj_mol': 1000.0, 'seed_score': row['score']})
    class ProtectedForceField:
        def createSystem(self, input_topology, **kwargs):
            system = forcefield.createSystem(input_topology, **kwargs)
            tether = mm.CustomTorsionForce('k*(1-cos(theta-target))')
            tether.addPerTorsionParameter('k')
            tether.addPerTorsionParameter('target')
            for row in targets:
                tether.addTorsion(*row['atoms'], [row['construction_strength_kj_mol'], row['target_radians']])
            system.addForce(tether)
            return system
    positions, refinement = _minimize_attempt(topology, initial * unit.nanometer,
        ProtectedForceField() if protocol == 'protect-allowed' else forcefield, modeled, max_iterations=1000)
    refined = np.asarray(positions.value_in_unit(unit.nanometer))
    np.savez(output / 'refined.npz', xyz_nm=refined)
    with (output / 'refined.pdb').open('w') as handle:
        app.PDBFile.writeFile(topology, positions, handle, keepIds=True)
    geometry = loop_geometry_report(topology, positions, modeled)
    templates = PDBFixer(filename=str(dataset / 'prepared.pdb')).templates
    stereo = stereochemistry_report(topology, positions, templates)
    stereo_error = None
    try:
        require_valid_stereochemistry(stereo, 'Paired research refinement')
    except ValueError as exc:
        stereo_error = str(exc)
    reference = rama_report(topology, positions, affected)
    saved = app.PDBFile(str(output / 'refined.pdb'))
    if [atom_key(a) for a in saved.topology.atoms()] != identities:
        raise ValueError('Saved refined PDB atom identity/order changed')
    roundtrip = rama_report(saved.topology, saved.positions, affected)
    original_classes = {tuple(r['residue']): r['classification'] for r in reference['rows']}
    if original_classes != {tuple(r['residue']): r['classification'] for r in roundtrip['rows']}:
        raise ValueError('Saved refined PDB changes an affected reference classification')
    fixed = [a.index for a in atoms if residue_key(a.residue) not in modeled and a.element != app.element.hydrogen]
    if not np.array_equal(initial[fixed], refined[fixed]):
        raise ValueError('An observed heavy atom moved')
    native_after = mm.XmlSerializer.serialize(fresh_system())
    if native_after != native_xml:
        raise ValueError('Temporary construction changed a fresh native physical System')
    refinement['extra_phi_psi_construction'] = {'enabled': protocol == 'protect-allowed',
        'targets': targets if protocol == 'protect-allowed' else [],
        'unprotected_seed_outliers': seed_reference['outliers'],
        'scope': 'Previously defined temporary construction tethers; no physical torsion fitting or acceptance cutoff change. Endpoint/phase energies include all temporary forces.'}
    for name, value in [('refinement', refinement), ('geometry', geometry), ('stereo', stereo),
                        ('reference', reference), ('saved-reference', roundtrip)]:
        save(output / (name + '.json'), value)
    result = {'protocol': protocol, 'initial_xyz_sha256': hashlib.sha256(initial.tobytes()).hexdigest(),
              'atom_identity_sha256': hashlib.sha256(json.dumps(identities).encode()).hexdigest(),
              'atom_count': len(atoms), 'modeled_residues': sorted(modeled), 'affected_residues': sorted(affected),
              'parameter_files': files, 'native_system_sha256': hashlib.sha256(native_xml.encode()).hexdigest(),
              'fresh_native_system_unchanged': True, 'construction_forces_exported': False,
              'fixed_heavy_atom_count': len(fixed), 'fixed_heavy_coordinates_bitwise_equal': True,
              'geometry_pass': geometry['accepted'], 'stereo_error': stereo_error,
              'modeled_outliers': [r['residue'] for r in reference['outliers'] if tuple(r['residue']) in modeled],
              'boundary_outliers': [r['residue'] for r in reference['outliers'] if tuple(r['residue']) not in modeled],
              'backbone_reference_pass': reference['all_selected_scored_without_outliers'],
              'passed_recorded_checks': geometry['accepted'] and stereo_error is None and reference['all_selected_scored_without_outliers'],
              'elapsed_refinement_seconds': refinement['elapsed_seconds'], 'app_ready': False,
              'physical_model_validated': False,
              'scope': 'One saved complete-complex refinement replay. No new preparation, dynamics, package or release admission.'}
    save(output / 'result.json', result)
    print(json.dumps(result, indent=2))


def main(dataset, output):
    output.mkdir(parents=True, exist_ok=False)
    sources = [Path(__file__), REFERENCE / 'rama_reference.py', REFERENCE / 'run_fragment_comparison.py',
               ROOT / 'backend/loop_refinement.py', ROOT / 'backend/loop_geometry.py',
               ROOT / 'backend/prepared_system.py', ROOT / 'backend/preparation_worker.py']
    sources += [p for p in dataset.rglob('*') if p.is_file() and p.suffix in {'.json', '.xml', '.npz', '.pdb', '.cif'}]
    hashes = {str(p): digest(p) for p in sources}
    save(output / 'plan.json', {'dataset': str(dataset), 'source_sha256': hashes,
        'protocols': ['original', 'protect-allowed'], 'iteration_budget_per_protocol': 1000,
        'purpose': 'Paired 12-residue loop replay with identical saved H atoms and prepared full-complex parameters. No threshold, observed-motion limit, or physical parameter change.'})
    snapshots = output / 'implementation'; snapshots.mkdir()
    for p in sources:
        if p.suffix == '.py':
            (snapshots / p.name).write_bytes(p.read_bytes())
    rows = []
    for protocol in ['original', 'protect-allowed']:
        try:
            supervisor = bounded([sys.executable, str(Path(__file__).resolve()), '--dataset', str(dataset),
                '--output', str(output / protocol), '--protocol', protocol], output, protocol, os.environ.copy())
            result = json.loads((output / protocol / 'result.json').read_text())
            result['supervisor'] = supervisor
            rows.append(result)
        except Exception as exc:
            rows.append({'protocol': protocol, 'error': str(exc), 'passed_recorded_checks': False})
        save(output / 'progress.json', {'protocols': rows})
    paired = not any('error' in r for r in rows)
    if paired:
        for field in ['initial_xyz_sha256', 'atom_identity_sha256', 'native_system_sha256']:
            if rows[0][field] != rows[1][field]:
                raise ValueError('Paired protocols did not use the same ' + field)
    if hashes != {name: digest(Path(name)) for name in hashes}:
        raise ValueError('A source changed during the paired experiment')
    save(output / 'result.json', {'protocols': rows, 'identical_start_and_physical_parameters_verified': paired,
        'source_sha256': hashes, 'source_files_unchanged': True, 'app_ready': False,
        'scope': 'Bounded research comparison; even a recorded pass does not establish native loop accuracy or MD validity.'})


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--dataset', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--protocol', choices=['original', 'protect-allowed'])
    args = parser.parse_args()
    if args.protocol:
        child(args.dataset.resolve(), args.output.resolve(), args.protocol)
    else:
        main(args.dataset.resolve(), args.output.resolve())
