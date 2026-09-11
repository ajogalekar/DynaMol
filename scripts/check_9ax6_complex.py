#!/usr/bin/env python3
"""Independently audit the retained 9AX6 complex and native ligand parameters.

Only --correct-annotations mutates datasets: it preserves each old annotation
and transparently repairs first-copy provenance, never coordinates or parameters.
"""
from __future__ import annotations

import argparse
import copy
import datetime as dt
import hashlib
import json
from pathlib import Path
import sys

import numpy as np
import openmm as mm
from openmm import app, unit

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from backend.prepared_system import load_prepared_forcefield


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def write_json(path, value):
    Path(path).write_text(json.dumps(value, indent=2, allow_nan=False) + '\n')


def key(atom):
    r = atom.residue
    return (r.chain.id, r.id, r.insertionCode.strip(), r.name, atom.name)


def native_hashes(folder):
    return {str(p.relative_to(folder)): sha(p) for p in sorted(folder.rglob('*')) if p.is_file()
            and p.name not in {'provenance.json', 'provenance-before-copy-fix.json', 'provenance-copy-correction.json'}
            and not p.name.startswith('residue-')}


def correct_annotations(folder, ligands):
    corrections = []
    for name in sorted({ligand['artifacts_directory'] for ligand in ligands}):
        entries = [ligand for ligand in ligands if ligand['artifacts_directory'] == name]
        artifact = folder / 'ligands' / name
        note_path = artifact / 'provenance-copy-correction.json'
        if note_path.exists():
            note = json.loads(note_path.read_text())
            assert sha(artifact / 'provenance.json') == note['corrected_provenance_sha256']
            assert native_hashes(artifact) == note['unchanged_parameter_and_native_files']
            corrections.append(note)
            continue
        before = native_hashes(artifact)
        original = artifact / 'provenance.json'
        preserved = artifact / 'provenance-before-copy-fix.json'
        assert not preserved.exists(), 'Refusing to overwrite preserved historical annotation.'
        preserved.write_bytes(original.read_bytes())
        original_record = json.loads(original.read_text())
        source = entries[0]
        source_key = source['key'].split(':')
        per_copy = []
        for entry in entries:
            corrected = copy.deepcopy(entry)
            corrected['parameterization_source_residue_key'] = source_key
            corrected['parameter_reuse_note'] = 'Identical ordered chemical states share parameters computed from the first bound copy; each copy retains its own heavy-atom coordinates.'
            corrected['annotation_correction'] = 'Added after completion to repair the folder annotation overwritten by the final identical copy. Native artifacts, coordinates, and XML were not changed.'
            instance_id = hashlib.sha256(entry['key'].encode()).hexdigest()[:12]
            instance_file = artifact / f'residue-{instance_id}.json'
            write_json(instance_file, corrected)
            per_copy.append({'key': entry['key'], 'path': str(instance_file.relative_to(ROOT)), 'sha256': sha(instance_file)})
            if entry is source:
                write_json(original, corrected)
        assert native_hashes(artifact) == before
        note = {'corrected_at': dt.datetime.now(dt.timezone.utc).isoformat(),
                'scope': 'Annotation correction only; byte identity is asserted for parameters/native outputs, not the corrected provenance annotation.',
                'artifact_directory': str(artifact.relative_to(ROOT)),
                'original_annotation_key': original_record.get('key'), 'correct_source_key': source['key'],
                'preserved_original_path': str(preserved.relative_to(ROOT)), 'preserved_original_sha256': sha(preserved),
                'corrected_provenance_sha256': sha(original), 'per_residue_annotations': per_copy,
                'unchanged_parameter_and_native_files': before}
        write_json(note_path, note)
        corrections.append(note)
    return corrections


def compare_coordinates(source, prepared, keys):
    original = {key(a): np.asarray(source.positions[a.index].value_in_unit(unit.angstrom)) for a in source.topology.atoms()}
    final = {key(a): np.asarray(prepared.positions[a.index].value_in_unit(unit.angstrom)) for a in prepared.topology.atoms()}
    assert all(k in final for k in keys), 'A protected heavy atom was lost.'
    displacement = {':'.join(k): float(np.linalg.norm(final[k] - original[k])) for k in keys}
    maximum = max(displacement.values(), default=0)
    assert maximum <= .0018, 'Protected atom moved beyond PDB decimal-rounding tolerance.'
    return {'count': len(keys), 'maximum_displacement_angstrom': maximum, 'displacements_angstrom': displacement}


def independently_check_template(xml_path):
    folder = xml_path.parent
    native_topology = app.AmberPrmtopFile(str(folder / 'ligand.prmtop'))
    positions = np.asarray(app.AmberInpcrdFile(str(folder / 'ligand.inpcrd')).positions.value_in_unit(unit.nanometer))
    systems = [native_topology.createSystem(nonbondedMethod=app.NoCutoff, constraints=None),
               app.ForceField(str(xml_path)).createSystem(native_topology.topology, nonbondedMethod=app.NoCutoff, constraints=None)]
    comparisons = []
    for trial in range(3):
        coordinates = positions if trial == 0 else positions + np.random.default_rng(2026 + trial).normal(0, .0003, positions.shape)
        values = []
        for system in systems:
            integrator = mm.VerletIntegrator(.001)
            context = mm.Context(system, integrator, mm.Platform.getPlatformByName('Reference'))
            context.setPositions(coordinates)
            state = context.getState(getEnergy=True, getForces=True)
            values.append((state.getPotentialEnergy().value_in_unit(unit.kilojoule_per_mole),
                           np.asarray(state.getForces(asNumpy=True).value_in_unit(unit.kilojoule_per_mole / unit.nanometer))))
            del context, integrator
        delta_energy = abs(values[1][0] - values[0][0])
        delta_force = float(np.max(np.abs(values[1][1] - values[0][1])))
        assert np.isfinite([delta_energy, delta_force]).all() and delta_energy < 1e-4 and delta_force < 1e-3
        comparisons.append({'pose': 'bound' if trial == 0 else f'perturbed-{trial}', 'energy_difference_kj_mol': delta_energy, 'maximum_force_difference_kj_mol_nm': delta_force})
    return comparisons


def audit(source_id, prepared_id, job_id, correction):
    source_folder = ROOT / 'data' / 'datasets' / source_id
    prepared_folder = ROOT / 'data' / 'datasets' / prepared_id
    metadata = json.loads((prepared_folder / 'metadata.json').read_text())
    preparation = metadata['preparation']
    ligands = preparation['ligand_parameters']['ligands']
    report = {'checked_at': dt.datetime.now(dt.timezone.utc).isoformat(), 'source_dataset': source_id,
              'prepared_dataset': prepared_id, 'preparation_job': job_id, 'script_sha256': sha(__file__),
              'source_topology_sha256': sha(source_folder / 'topology.pdb'),
              'prepared_pdb_sha256': sha(prepared_folder / 'prepared.pdb'),
              'source_version_limitation': 'The application source was edited while this development job ran. Its recorded implementation hash identifies file bytes observed at completion, not a guaranteed immutable snapshot of loaded Python modules. This audit independently rechecks preserved molecular artifacts with the current audit script; it does not retroactively claim a frozen build.'}
    if correction:
        report['annotation_corrections'] = correct_annotations(ROOT / 'data' / 'jobs' / job_id, ligands) + correct_annotations(prepared_folder, ligands)
    else:
        report['annotation_corrections'] = [json.loads(p.read_text()) for folder in [ROOT / 'data' / 'jobs' / job_id, prepared_folder]
                                            for p in sorted((folder / 'ligands').glob('*/provenance-copy-correction.json'))]
    source = app.PDBFile(str(source_folder / 'topology.pdb'))
    prepared = app.PDBFile(str(prepared_folder / 'prepared.pdb'))
    ligand_names = {'GNP', 'A1A'}
    ligand_keys = [key(a) for a in source.topology.atoms() if a.residue.name in ligand_names and a.element != app.element.hydrogen]
    metal_keys = [key(a) for a in source.topology.atoms() if a.residue.name == 'MG']
    assert len(ligand_keys) == 180 and len(metal_keys) == 2
    prepared_ligand_count = sum(a.residue.name in ligand_names and a.element != app.element.hydrogen for a in prepared.topology.atoms())
    assert prepared_ligand_count == 180
    source_xyz = {key(a): np.asarray(source.positions[a.index].value_in_unit(unit.angstrom)) for a in source.topology.atoms()}
    prepared_xyz = {key(a): np.asarray(prepared.positions[a.index].value_in_unit(unit.angstrom)) for a in prepared.topology.atoms()}
    contacts = []
    for metal in metal_keys:
        for donor in source.topology.atoms():
            if donor.element.symbol not in {'N', 'O', 'S'}:
                continue
            donor_key = key(donor)
            distance = float(np.linalg.norm(source_xyz[metal] - source_xyz[donor_key]))
            if distance <= 2.8:
                new_distance = float(np.linalg.norm(prepared_xyz[metal] - prepared_xyz[donor_key]))
                contacts.append({'metal': list(metal), 'donor': list(donor_key), 'original_distance_angstrom': distance,
                                 'prepared_distance_angstrom': new_distance, 'absolute_delta_angstrom': abs(distance - new_distance)})
    assert len(contacts) == 12
    water_keys = sorted({tuple(c['donor']) for c in contacts if c['donor'][3] == 'HOH'})
    assert len(water_keys) == 4
    assert sum(a.residue.name == 'MG' for a in prepared.topology.atoms()) == 2
    assert sum(a.residue.name == 'HOH' and a.element.symbol == 'O' for a in prepared.topology.atoms()) == 4
    protected_keys = [key(a) for a in source.topology.atoms() if a.residue.chain.id in {'A', 'B'} and (a.residue.id, a.residue.name) in {('17', 'SER'), ('35', 'THR')} and a.element != app.element.hydrogen]
    assert len(protected_keys) == 26
    report['ligand_heavy_atoms'] = compare_coordinates(source, prepared, ligand_keys)
    report['magnesium_atoms'] = compare_coordinates(source, prepared, metal_keys)
    report['coordinating_water_oxygens'] = compare_coordinates(source, prepared, water_keys)
    report['protected_ser17_thr35_heavy_atoms'] = compare_coordinates(source, prepared, protected_keys)
    report['metal_contacts'] = contacts
    report['maximum_contact_distance_change_angstrom'] = max(c['absolute_delta_angstrom'] for c in contacts)
    assert report['maximum_contact_distance_change_angstrom'] <= .0036
    ff, files = load_prepared_forcefield(prepared_folder, preparation, solvent='explicit')
    system = ff.createSystem(prepared.topology, nonbondedMethod=app.NoCutoff, constraints=None)
    nb = next(force for force in system.getForces() if isinstance(force, mm.NonbondedForce))
    residues = {':'.join((r.chain.id, r.id, r.insertionCode.strip(), r.name)): r for r in prepared.topology.residues()}
    charge_checks = []
    for ligand in ligands:
        atoms = list(residues[ligand['key']].atoms())
        native = app.AmberPrmtopFile(str(prepared_folder / 'ligands' / ligand['artifacts_directory'] / 'ligand.prmtop')).createSystem(nonbondedMethod=app.NoCutoff, constraints=None)
        native_nb = next(force for force in native.getForces() if isinstance(force, mm.NonbondedForce))
        actual = [nb.getParticleParameters(a.index)[0].value_in_unit(unit.elementary_charge) for a in atoms]
        assigned = [native_nb.getParticleParameters(i)[0].value_in_unit(unit.elementary_charge) for i in range(native.getNumParticles())]
        assert len(atoms) == len(assigned) == len(ligand['charges_e'])
        maximum = float(np.max(np.abs(np.asarray(actual) - assigned)))
        assert maximum < 1e-8
        assert np.max(np.abs(np.asarray(actual) - ligand['charges_e'])) < 1e-8
        assert abs(sum(actual) - ligand['formal_charge']) < 1e-6
        charge_checks.append({'key': ligand['key'], 'component_id': ligand['component_id'], 'formal_charge': ligand['formal_charge'],
                              'assembled_charge_sum_e': sum(actual), 'native_charge_sum_e': sum(assigned), 'maximum_atom_charge_difference_e': maximum,
                              'atoms': [{'name': atom.name, 'assembled_index': atom.index, 'assembled_charge_e': a, 'native_charge_e': n} for atom, a, n in zip(atoms, actual, assigned)]})
    report['assembled_system_particles'] = system.getNumParticles()
    report['assembled_forcefield_files'] = files
    report['ligand_charge_checks'] = charge_checks
    report['conversion_validation_reports'] = []
    for entry in preparation['ligand_parameters']['files']:
        xml_path = prepared_folder / entry['path']
        assert sha(xml_path) == entry['sha256']
        path = xml_path.parent / 'conversion-validation.json'
        validation = json.loads(path.read_text())
        assert validation['passed'] and len(validation['conformations']) == 3
        assert all(abs(c['energy_delta_kj_mol']) < 1e-4 and c['maximum_force_delta_kj_mol_nm'] < 1e-3 for c in validation['conformations'])
        report['conversion_validation_reports'].append({'path': str(path.relative_to(ROOT)), 'sha256': sha(path), 'validation': validation, 'independent_recalculation': independently_check_template(xml_path)})
    report['passed'] = True
    return report


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--source', default='9fb26d98962c453b')
    parser.add_argument('--prepared', default='4cf5bcedf38847c1')
    parser.add_argument('--job', default='57083a2415754858')
    parser.add_argument('--correct-annotations', action='store_true')
    parser.add_argument('--output', type=Path, default=ROOT / 'docs' / 'audit' / 'complex-preparation-checks.json')
    arguments = parser.parse_args()
    result = audit(arguments.source, arguments.prepared, arguments.job, arguments.correct_annotations)
    write_json(arguments.output, result)
    print(json.dumps({'passed': result['passed'], 'ligand_heavy_atoms': result['ligand_heavy_atoms']['count'],
                      'metal_contacts': len(result['metal_contacts']), 'maximum_ligand_displacement_angstrom': result['ligand_heavy_atoms']['maximum_displacement_angstrom'],
                      'maximum_contact_change_angstrom': result['maximum_contact_distance_change_angstrom'], 'report': str(arguments.output)}, indent=2))
