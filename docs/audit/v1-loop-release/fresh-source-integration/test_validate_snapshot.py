"""Private wrapper integrity tests; refinement is always replaced with a stub.

These synthetic fixtures establish byte/topology binding only. They are not
physical model validation and never construct a Context or run minimization/MD.
"""
import hashlib
import importlib.util
import json
from pathlib import Path
from types import SimpleNamespace

import numpy as np
import openmm as mm
from openmm import app, unit
import pytest

from backend.loop_search import REQUIRED_CHECKS
from backend.loop_snapshot import create_snapshot, verify_accepted_bundle, _inventory


SPEC = importlib.util.spec_from_file_location(
    'private_validate_snapshot', Path(__file__).with_name('validate_snapshot.py'))
validator = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(validator)


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


@pytest.fixture
def frozen_case(tmp_path, monkeypatch):
    root = tmp_path.resolve()
    inputs = root / 'inputs'; inputs.mkdir()
    old = app.Topology()
    chain = old.addChain('L')
    residue = old.addResidue('LIG', chain, '99')
    old.addAtom('CX', app.element.carbon, residue)
    top = app.Topology()
    chain = top.addChain('L')
    residue = top.addResidue('LIG', chain, '99')
    top.addAtom('CX', app.element.carbon, residue)
    residue = top.addResidue('GLY', top.addChain('A'), '2', 'B')
    atoms = [top.addAtom(name, element, residue) for name, element in
             [('N', app.element.nitrogen), ('CA', app.element.carbon),
              ('C', app.element.carbon), ('O', app.element.oxygen)]]
    for a, b in zip(atoms, atoms[1:]):
        top.addBond(a, b, type=app.Single, order=1)
    xyz = np.asarray([[0., 0., 0.], [1., 1., 1.], [1.1, 1., 1.],
                      [1.2, 1.1, 1.], [1.2, 1.2, 1.1]])
    with (inputs / 'prepared-topology.cif').open('w') as handle:
        app.PDBxFile.writeFile(top, xyz * unit.nanometer, handle, keepIds=True)
    system = mm.System()
    for atom in top.atoms():
        system.addParticle(atom.element.mass)
    system.addForce(mm.HarmonicBondForce())
    xml = mm.XmlSerializer.serialize(system)
    (inputs / 'native-system.xml').write_text(xml)
    roots = ['amber14/protein.ff14SB.xml', 'amber14/tip3p.xml']
    state = {'solvent': 'explicit', 'parameter_state': {}, 'forcefield_files': roots}
    (inputs / 'parameter-state.json').write_text(json.dumps(state))
    (inputs / 'loop-model-context.pdb').write_text('Synthetic immutable context bytes.\n')
    (inputs / 'loop-model-input.json').write_text('{}\n')
    parameters = {}
    for i, name in enumerate(roots):
        path = inputs / f'root-{i}.xml'; path.write_text('<ForceField/>\n')
        parameters[name] = path
    snap = create_snapshot(
        root / 'snapshot', source_topology=old,
        source_positions=xyz[:1] * unit.nanometer,
        topology=top, positions=xyz * unit.nanometer,
        modeled_keys=[('A', '2', 'B', 'GLY')],
        requested_modeled_keys=[('A', '2', 'B', 'GLY')], source_to_prepared=[0],
        preparation_provenance={'purpose': 'synthetic static integrity fixture'},
        parameter_files=parameters,
        source_files={p.name: p for p in inputs.iterdir() if not p.name.startswith('root-')})
    proposal = root / 'proposal.json'
    proposal.write_text(json.dumps({'source_sha256': sha(inputs / 'loop-model-context.pdb'),
                                    'marker': 'original'}))
    request = {'snapshot': str(snap.manifest_path.parent), 'snapshot_sha256': snap.sha256,
               'proposal': str(proposal), 'proposal_sha256': sha(proposal),
               'output': str(root / 'output')}
    monkeypatch.setattr(validator, 'load_prepared_forcefield',
                        lambda *args, **kwargs: (SimpleNamespace(createSystem=lambda *a, **k: system), roots))
    monkeypatch.setattr(validator, 'new_reference_keys', lambda *args: set())

    def static_stub(topology, positions, forcefield, modeled, candidate, output, **kwargs):
        output.mkdir()
        np.savez(output / 'refined.npz', xyz_nm=xyz)
        np.savez(output / 'candidate-seed.npz', xyz_nm=xyz)
        with (output / 'topology.cif').open('w') as handle:
            app.PDBxFile.writeFile(topology, positions, handle, keepIds=True)
        with (output / 'refined.pdb').open('w') as handle:
            app.PDBFile.writeFile(topology, positions, handle, keepIds=True)
        (output / 'native-system.xml').write_text(xml)
        result = {'accepted': True, 'checks': dict.fromkeys(REQUIRED_CHECKS, True),
                  'observed_context_residue_keys': [],
                  'artifacts_sha256': {p.name: sha(p) for p in output.iterdir()}}
        (output / 'result.json').write_text(json.dumps(result))
        return result

    monkeypatch.setattr(validator, 'validate_complete_loop', static_stub)
    return SimpleNamespace(root=root, inputs=inputs, snapshot=snap, top=top,
                           xyz=xyz, proposal=proposal, request=request, stub=static_stub)


def test_static_fixture_binds_original_candidate_and_physical_system(frozen_case):
    case = frozen_case
    result = validator.run(case.request)
    assert result['status'] == 'static_checks_passed'
    assert result['native_parameter_check']['all_atom_indexed_parameters_equal']
    assert result['full_preparation_published'] is False
    bundle = verify_accepted_bundle(result['accepted_bundle'],
                                    expected_sha256=result['accepted_bundle_sha256'])
    provenance = {f.name: f for f in bundle.artifacts['provenance']}
    assert provenance['proposal.json'].sha256 == case.request['proposal_sha256']


def test_topology_preserves_author_identities_elements_and_typed_bonds(frozen_case):
    case = frozen_case
    top, error = validator.topology_from_snapshot(case.snapshot, case.inputs / 'prepared-topology.cif')
    assert _inventory(top) == _inventory(case.top)
    assert error < 1e-12


def test_topology_preserves_periodic_cell(frozen_case):
    case = frozen_case
    box = np.diag([3., 4., 5.]) * unit.nanometer
    case.top.setPeriodicBoxVectors(box)
    path = case.root / 'boxed.cif'
    with path.open('w') as handle:
        app.PDBxFile.writeFile(case.top, case.xyz * unit.nanometer, handle, keepIds=True)
    top, _ = validator.topology_from_snapshot(case.snapshot, path)
    assert top.getPeriodicBoxVectors() is not None
    np.testing.assert_allclose(top.getPeriodicBoxVectors().value_in_unit(unit.nanometer),
                               box.value_in_unit(unit.nanometer), atol=1e-12)


def test_changed_candidate_cannot_be_bundled_with_original_hash(frozen_case, monkeypatch):
    case = frozen_case

    def mutate(*args, **kwargs):
        result = case.stub(*args, **kwargs)
        case.proposal.write_text(json.dumps({'source_sha256': '0' * 64, 'marker': 'replacement'}))
        return result

    monkeypatch.setattr(validator, 'validate_complete_loop', mutate)
    # Rejecting a changed source or preserving the original verified buffer are
    # both valid policies. Claiming the original hash for replacement bytes is not.
    try:
        result = validator.run(case.request)
    except ValueError:
        return
    assert result['status'] == 'static_checks_passed'
    assert sha(case.root / 'output' / 'proposal.json') == case.request['proposal_sha256']


def test_changed_validated_artifact_is_rejected_before_bundle(frozen_case, monkeypatch):
    case = frozen_case

    def mutate(*args, **kwargs):
        result = case.stub(*args, **kwargs)
        # Change a ligand coordinate after the validator recorded its output
        # hashes, preserving count and finite coordinate shape.
        xyz = case.xyz.copy(); xyz[0, 0] += 0.2
        np.savez(args[5] / 'refined.npz', xyz_nm=xyz)
        return result

    monkeypatch.setattr(validator, 'validate_complete_loop', mutate)
    with pytest.raises(ValueError):
        validator.run(case.request)


def test_bundle_capture_rejects_late_artifact_change(frozen_case, monkeypatch):
    case = frozen_case
    original = validator.create_accepted_bundle

    def mutate_then_capture(*args, **kwargs):
        path = kwargs['artifacts']['topology']['native-system.xml']
        path.write_text('<changed-after-validator-hash-check/>\n')
        return original(*args, **kwargs)

    monkeypatch.setattr(validator, 'create_accepted_bundle', mutate_then_capture)
    with pytest.raises(ValueError, match='validator-approved artifact bytes'):
        validator.run(case.request)


def test_proposal_from_other_source_is_rejected(frozen_case):
    case = frozen_case
    case.proposal.write_text(json.dumps({'source_sha256': '0' * 64}))
    case.request['proposal_sha256'] = sha(case.proposal)
    with pytest.raises(ValueError, match='different observed source context'):
        validator.run(case.request)


def test_symlink_proposal_is_rejected(frozen_case):
    case = frozen_case
    link = case.root / 'proposal-link.json'; link.symlink_to(case.proposal)
    case.request['proposal'] = str(link)
    with pytest.raises(ValueError, match='Symlink'):
        validator.run(case.request)


def test_nonmatching_native_system_is_rejected(frozen_case, monkeypatch):
    case = frozen_case
    system = mm.System()
    for atom in case.top.atoms():
        system.addParticle(atom.element.mass)
    system.addForce(mm.HarmonicBondForce())
    system.setParticleMass(0, 13.0 * unit.dalton)
    names = ['amber14/protein.ff14SB.xml', 'amber14/tip3p.xml']
    monkeypatch.setattr(validator, 'load_prepared_forcefield',
                        lambda *args, **kwargs: (SimpleNamespace(createSystem=lambda *a, **k: system), names))
    with pytest.raises(ValueError, match='atom-indexed physical System'):
        validator.run(case.request)


@pytest.mark.parametrize('payload', ['{"x": 1, "x": 2}', '{"x": NaN}'])
def test_json_rejects_duplicate_and_nonfinite_values(tmp_path, payload):
    path = tmp_path / 'bad.json'; path.write_text(payload)
    with pytest.raises(ValueError):
        validator.read_json(path)
