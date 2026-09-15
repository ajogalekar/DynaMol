"""Private fresh-snapshot qualification child; never publishes an app dataset.

Run only under backend.loop_process supervision. The snapshot, proposal and
request are immutable inputs. A successful numerical result is still a static
model check, not evidence for the native conformation or a released workflow.
"""
import hashlib
import json
from pathlib import Path
import sys

import numpy as np
import openmm as mm
from openmm import app, unit
from openmm.app import topology as topology_types

from backend.loop_candidate_validation import (atom_key, physical_parameter_comparison,
                                                validate_complete_loop)
from backend.loop_context import defining_indices
from backend.loop_snapshot import (load_snapshot, create_accepted_bundle,
                                   verify_accepted_bundle, _inventory, _read, _json)
from backend.prepared_system import load_prepared_forcefield
from backend.residue_identity import STANDARD_PROTEINS, residue_key


def read_json(path):
    def unique(pairs):
        result = {}
        for key, value in pairs:
            if key in result:
                raise ValueError('Duplicate request field')
            result[key] = value
        return result
    return json.loads(Path(path).read_text(), object_pairs_hook=unique,
                      parse_constant=lambda x: (_ for _ in ()).throw(ValueError('Nonfinite JSON')))


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def write_json(path, value):
    data = (json.dumps(value, indent=2, allow_nan=False) + '\n').encode()
    with Path(path).open('xb') as handle:
        handle.write(data)
    return hashlib.sha256(data).hexdigest()


def topology_from_snapshot(snapshot, common_cif):
    """Preserve CIF chain segmentation and the exact captured bonded inventory."""
    from backend.modified_residues import register_topology_definitions
    register_topology_definitions()
    common = app.PDBxFile(str(common_cif))
    ordered = [(atom_key(a), a.element.symbol) for a in common.topology.atoms()]
    expected = [(a.identity, a.element) for a in snapshot.prepared.atoms]
    if ordered != expected:
        raise ValueError('Common CIF does not reproduce the captured ordered atom inventory')
    xyz = np.asarray(common.positions.value_in_unit(unit.nanometer))
    if xyz.shape != snapshot.prepared_xyz_nm.shape or not np.isfinite(xyz).all():
        raise ValueError('Common CIF coordinates are malformed')
    # CIF stores finite decimal coordinates; the separately frozen NPY is exact.
    error = float(np.max(np.linalg.norm(xyz - snapshot.prepared_xyz_nm, axis=1)))
    if error > .00009:
        raise ValueError('Common CIF differs from the captured coordinates beyond serialization rounding')
    top = app.Topology()
    atoms = []
    for old_chain in common.topology.chains():
        chain = top.addChain(old_chain.id)
        for old_residue in old_chain.residues():
            residue = top.addResidue(old_residue.name, chain, old_residue.id, old_residue.insertionCode)
            for old in old_residue.atoms():
                atoms.append(top.addAtom(old.name, old.element, residue, old.id,
                                         formalCharge=old.formalCharge))
    types = {name: getattr(topology_types, name) for name in ('Single', 'Double', 'Triple', 'Aromatic', 'Amide')}
    for i, j, kind, order in snapshot.prepared.bonds:
        if kind is not None and kind not in types:
            raise ValueError('Unknown captured bond type')
        top.addBond(atoms[i], atoms[j], type=types.get(kind), order=order)
    top.setPeriodicBoxVectors(common.topology.getPeriodicBoxVectors())
    inventory = _inventory(top)
    if tuple(tuple(row) for row in inventory['bonds']) != snapshot.prepared.bonds:
        raise ValueError('Reconstructed native bonded inventory differs from the snapshot')
    return top, error


def new_reference_keys(snapshot, topology):
    """Find every newly defined phi/psi/omega support using explicit raw mapping."""
    from backend.loop_backbone_reference import rama_report
    retained = {i for i in snapshot.source_to_prepared if i is not None}
    new = {i for i, a in enumerate(snapshot.prepared.atoms)
           if i not in retained and a.atomic_number != 1 and a.identity[-1] in {'N', 'CA', 'C'}}
    report = rama_report(topology, snapshot.prepared_xyz_nm * unit.nanometer)
    keys = {tuple(row['residue']) for row in report['rows']
            if defining_indices(topology, row) & new}
    return keys


def source_movement(snapshot, final, flanks):
    mapped = [(old, new) for old, new in enumerate(snapshot.source_to_prepared)
              if new is not None and snapshot.source.atoms[old].atomic_number != 1]
    rows = [{'source': snapshot.source.atoms[old].identity,
             'prepared': snapshot.prepared.atoms[new].identity,
             'common_displacement_angstrom': float(np.linalg.norm(snapshot.prepared_xyz_nm[new] - snapshot.source_xyz_nm[old]) * 10),
             'final_displacement_angstrom': float(np.linalg.norm(final[new] - snapshot.source_xyz_nm[old]) * 10)}
            for old, new in mapped]
    scoped = [row for row in rows if tuple(row['prepared'][:4]) in flanks]
    maximum = max((row['final_displacement_angstrom'] for row in scoped), default=0.)
    return {'mapped_heavy_atoms': len(rows), 'rows': rows,
            'maximum_original_source_flank_displacement_angstrom': maximum,
            'original_source_flank_cap_passed': maximum <= 1.0,
            'scope': 'Original-source movement is separate from the common-preparation baseline; the same 1 A final flank cap applies to both.'}


def run(request):
    required = {'snapshot', 'snapshot_sha256', 'proposal', 'proposal_sha256', 'output'}
    if not required <= request.keys() or set(request) - required - {'preserve_seed_torsions'}:
        raise ValueError('Exact snapshot, proposal, hashes and output are required')
    if type(request.get('preserve_seed_torsions', False)) is not bool:
        raise ValueError('Seed preservation requires an explicit boolean')
    output = Path(request['output']).resolve()
    output.mkdir(parents=True, exist_ok=False)
    snapshot = load_snapshot(request['snapshot'], expected_sha256=request['snapshot_sha256'])
    sources = {f.name: f for f in snapshot.source_files}
    required_sources = {'prepared-topology.cif', 'native-system.xml', 'parameter-state.json',
                        'loop-model-context.pdb', 'loop-model-input.json'}
    if not required_sources <= sources.keys():
        raise ValueError('Fresh snapshot omits required common-preparation dependencies')
    # One bounded read binds the validated object and the retained generator
    # evidence to exactly the requested bytes, including symlink/race checks.
    proposal_path = Path(request['proposal']).absolute()
    proposal_bytes = _read(proposal_path, request['proposal_sha256'], limit=32 * 1024**2)
    proposal = _json(proposal_bytes)
    source_hash = sources['loop-model-context.pdb'].sha256
    if proposal.get('source_sha256') != source_hash:
        raise ValueError('Native proposal was built against a different observed source context')
    state = read_json(sources['parameter-state.json'].path)
    solvent = state['solvent']
    root_names = ['amber14/protein.ff14SB.xml', 'amber14/tip3p.xml' if solvent == 'explicit' else 'implicit/gbn2.xml']
    parameters = {f.name: f.path for f in snapshot.parameter_files}
    roots = {name: parameters[name] for name in root_names}
    top, cif_error = topology_from_snapshot(snapshot, sources['prepared-topology.cif'].path)
    ff, names = load_prepared_forcefield(snapshot.manifest_path.parent / 'parameters',
                                       state['parameter_state'], solvent=solvent, base_parameter_paths=roots)
    if names != state['forcefield_files']:
        raise ValueError('Pinned force-field root selection differs from common preparation')
    native = mm.XmlSerializer.serialize(ff.createSystem(top, nonbondedMethod=app.NoCutoff,
                                                        constraints=None, rigidWater=False))
    parameter_check = physical_parameter_comparison(sources['native-system.xml'].path.read_text(), native)
    if not parameter_check['all_atom_indexed_parameters_equal']:
        raise ValueError('Pinned parameters did not reproduce the captured atom-indexed physical System')
    templates = {Path(name).stem: app.PDBFile(str(record.path)) for name, record in sources.items()
                 if name.startswith('templates/') and name.endswith('.pdb')}
    required_templates = {r.name for r in top.residues() if r.name in STANDARD_PROTEINS and r.name != 'GLY'}
    if not required_templates <= templates.keys():
        raise ValueError('Snapshot is missing standard stereochemistry reference templates')
    newly_defined = new_reference_keys(snapshot, top)
    result = validate_complete_loop(top, snapshot.prepared_xyz_nm * unit.nanometer, ff,
                                    snapshot.modeled_keys, proposal, output / 'validation',
                                    source_sha256=source_hash, templates=templates,
                                    protected_residues=state.get('protected_residues', ()),
                                    newly_defined_backbone_residue_keys=newly_defined,
                                    preserve_seed_torsions=request.get('preserve_seed_torsions', False))
    hashes = result['artifacts_sha256']
    for name, expected in hashes.items():
        if Path(name).name != name:
            raise ValueError('Validator artifact name escapes its private output')
        _read(output / 'validation' / name, expected)
    result_bytes = _read(output / 'validation' / 'result.json')
    if _json(result_bytes) != _json(json.dumps(result, allow_nan=False)):
        raise ValueError('Saved validator decision differs from the returned result')
    result_hash = hashlib.sha256(result_bytes).hexdigest()
    import io
    with np.load(io.BytesIO(_read(output / 'validation' / 'refined.npz', hashes['refined.npz'])), allow_pickle=False) as archive:
        final = archive['xyz_nm']
    movement = source_movement(snapshot, final, set(map(tuple, result['observed_context_residue_keys'])))
    movement_hash = write_json(output / 'original-source-movement.json', movement)
    binding = {'snapshot_sha256': snapshot.sha256, 'proposal_sha256': request['proposal_sha256'],
               'original_source_context_sha256': source_hash, 'common_cif_maximum_rounding_nm': cif_error,
               'native_parameter_check': parameter_check, 'newly_defined_reference_keys': sorted(newly_defined),
               'all_required_checks': result['checks'], 'loop_candidate_accepted': result['accepted'],
               'original_source_flank_cap_passed': movement['original_source_flank_cap_passed'],
               'preserve_seed_torsions': request.get('preserve_seed_torsions', False),
               'app_default_enabled': False, 'full_preparation_published': False, 'physical_model_validated': False}
    binding_hash = write_json(output / 'binding.json', binding)
    # Bind the exact supplied candidate bytes; never regenerate JSON as evidence
    # of the original generator output hash.
    with (output / 'proposal.json').open('xb') as handle:
        handle.write(proposal_bytes)
    if not result['accepted'] or not movement['original_source_flank_cap_passed']:
        return {'status': 'rejected', **binding}
    artifacts = {
        'coordinates': {'refined.npz': output / 'validation' / 'refined.npz'},
        'topology': {name: output / 'validation' / name for name in ['topology.cif', 'refined.pdb', 'native-system.xml']},
        'parameters': parameters,
        'provenance': {p.name: p for p in (output / 'validation').glob('*.json')},
    }
    artifacts['provenance'].update({p.name: p for p in [output / 'binding.json', output / 'proposal.json',
                                                      output / 'original-source-movement.json',
                                                      output / 'validation' / 'candidate-seed.npz']})
    bundle = create_accepted_bundle(output / 'accepted-bundle', snapshot=snapshot, artifacts=artifacts,
                                    checks=result['checks'])
    verified = verify_accepted_bundle(bundle.manifest_path, expected_sha256=bundle.sha256)
    expected_outputs = {**hashes, 'result.json': result_hash,
                        'proposal.json': request['proposal_sha256'],
                        'binding.json': binding_hash,
                        'original-source-movement.json': movement_hash}
    for role in ('coordinates', 'topology', 'provenance'):
        for record in verified.artifacts[role]:
            if record.name not in expected_outputs or record.sha256 != expected_outputs[record.name]:
                raise ValueError('Captured bundle differs from the validator-approved artifact bytes')
    return {'status': 'static_checks_passed', **binding,
            'accepted_bundle': str(verified.manifest_path), 'accepted_bundle_sha256': verified.sha256}


if __name__ == '__main__':
    request = read_json(sys.argv[1])
    try:
        result = run(request)
    except Exception as exc:
        result = {'status': 'failed', 'error_type': type(exc).__name__, 'error': str(exc),
                  'full_preparation_published': False, 'app_default_enabled': False}
        write_json(sys.argv[2], result)
        raise
    write_json(sys.argv[2], result)
    raise SystemExit(0 if result['status'] == 'static_checks_passed' else 2)
