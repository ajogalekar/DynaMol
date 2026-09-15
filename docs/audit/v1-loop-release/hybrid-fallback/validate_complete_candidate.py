"""Private complete-complex candidate check against archived native parameters.

No parameter fitting or app admission. Candidate backbones are proposals;
hydrogen identities/states and every retained molecule come from the archive.
"""
import argparse
import hashlib
import json
from pathlib import Path
import random
import sys

import numpy as np
import openmm as mm
from openmm import app, unit
from pdbfixer import PDBFixer

ROOT = Path(__file__).resolve().parents[4]
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / 'docs/audit/advanced-chemistry/covalent-v1-focused'))
from backend.loop_search import REQUIRED_CHECKS
from backend.complex_topology import metal_environment
from backend.loop_geometry import loop_geometry_report, _dihedral
from backend.loop_refinement import _minimize_attempt
from backend.prepared_system import load_prepared_forcefield
from backend.preparation_worker import stereochemistry_report, require_valid_stereochemistry
from backend.residue_identity import STANDARD_PROTEINS, residue_key
from rama_reference import rama_report
from context_integrity import (outer_torsion_anchor_indices, affected_reference_keys,
                               outside_observed_torsions, AnchoredConstructionForceField)


def save(path, value):
    path.write_text(json.dumps(value, indent=2, allow_nan=False) + '\n')


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def key(atom):
    return (*residue_key(atom.residue), atom.name)


def bond_keys(topology):
    return {tuple(sorted((key(a), key(b)))) for a, b in topology.bonds()}


def observed_peptides(topology, before, after, modeled, flanks):
    rows = []
    for a, b in topology.bonds():
        if a.residue is b.residue or {a.name, b.name} != {'C', 'N'}:
            continue
        pair = {residue_key(a.residue), residue_key(b.residue)}
        if pair & modeled or not pair & flanks:
            continue
        left, right = (a, b) if a.name == 'C' else (b, a)
        names = [{x.name: x.index for x in r.atoms()} for r in [left.residue, right.residue]]
        indices = [names[0]['CA'], left.index, right.index, names[1]['CA']]
        original, final = _dihedral(before[indices]), _dihedral(after[indices])
        valid = original is not None and final is not None
        valid = valid and abs(abs(original) - 90) > 1e-6 and (abs(original) < 90) == (abs(final) < 90)
        rows.append({'residues': sorted(pair), 'source_omega_degrees': original,
                     'final_omega_degrees': final, 'same_observed_basin': bool(valid)})
    return rows


def frame(named, xyz):
    origin = xyz[named['CA']]
    x = xyz[named['C']] - origin
    x /= np.linalg.norm(x)
    y = xyz[named['N']] - origin
    y -= np.dot(y, x) * x
    y /= np.linalg.norm(y)
    rotation = np.column_stack((x, y, np.cross(x, y)))
    if not np.isfinite(rotation).all() or np.linalg.det(rotation) < .999:
        raise ValueError('Degenerate residue frame')
    return origin, rotation


def candidate_rows(candidate):
    if candidate.get('archive_baseline'):
        return [], [], set()
    if 'generated_atoms' in candidate:
        generated = candidate['generated_atoms']
        context = candidate.get('research_context_atoms', [])
        selected = {tuple(r['identity'][:4]) for r in generated}
    else:
        selected = set(map(tuple, candidate['modeled_residue_keys']))
        generated = candidate.get('modeled_heavy_atoms') or [r for r in candidate['backbone_atoms'] if tuple(r['identity'][:4]) in selected]
        context = candidate.get('observed_context_atoms') or [r for r in candidate['backbone_atoms'] if tuple(r['identity'][:4]) not in selected]
    return generated, context, selected


def run(dataset, candidate_path, output, observed_context_policy='mobile_flanks'):
    output.mkdir(parents=True, exist_ok=False)
    inputs = [candidate_path, dataset / 'loop-refinement-parameters.json',
              dataset / 'loop-refinement-topology.cif', dataset / 'loop-refinement-input.npz']
    hashes = {str(p): sha(p) for p in inputs}
    candidate = json.loads(candidate_path.read_text())
    parameters = json.loads(inputs[1].read_text())
    structure = app.PDBxFile(str(inputs[2]))
    topology = structure.topology
    initial = np.load(inputs[3])['xyz_nm'].copy()
    atoms = list(topology.atoms())
    index = {key(a): a.index for a in atoms}
    if len(index) != len(atoms):
        raise ValueError('Archived atom identity is ambiguous')
    modeled = set(map(tuple, parameters['modeled_residue_keys']))
    scaffold_binding = {'archive_baseline': bool(candidate.get('archive_baseline'))}
    if not candidate.get('archive_baseline'):
        source_path = Path(candidate.get('source') or candidate['immutable_scaffold'])
        if sha(source_path) != candidate['source_sha256']:
            raise ValueError('Candidate scaffold hash differs from generation provenance')
        hashes[str(source_path)] = sha(source_path)
        source = app.PDBFile(str(source_path))
        source_xyz = np.asarray(source.positions.value_in_unit(unit.nanometer))
        source_rows = [a for a in source.topology.atoms()
                       if a.residue.name in STANDARD_PROTEINS and a.name in {'N', 'CA', 'C', 'O'}]
        source_keys = {key(a) for a in source_rows}
        expected_keys = {key(a) for a in atoms if a.residue.name in STANDARD_PROTEINS
                         and residue_key(a.residue) not in modeled and a.name in {'N', 'CA', 'C', 'O'}}
        if len(source_keys) != len(source_rows) or source_keys != expected_keys:
            raise ValueError('Candidate scaffold does not bind to the complete archived protein backbone')
        displacement = max(np.linalg.norm(source_xyz[a.index] - initial[index[key(a)]])
                           for a in source_rows)
        if displacement > .0002:
            raise ValueError('Candidate scaffold uses a different observed backbone frame')
        scaffold_binding = {'source': str(source_path), 'source_sha256': candidate['source_sha256'],
                            'matched_observed_backbone_atoms': len(source_rows),
                            'maximum_source_rounding_difference_angstrom': float(displacement * 10),
                            'retained_molecules_from': str(dataset)}
    generated, context_rows, candidate_selected = candidate_rows(candidate)
    if not candidate_selected <= modeled:
        raise ValueError('Candidate sequence differs from the requested missing sequence')
    rows = generated + context_rows
    if len({tuple(r['identity']) for r in rows}) != len(rows):
        raise ValueError('Duplicate candidate atom identity')
    updated = initial.copy()
    explicit = set()
    flanks = set()
    for row in rows:
        identity = tuple(row['identity'])
        if identity not in index:
            raise ValueError('Unknown candidate atom: ' + str(identity))
        atom = atoms[index[identity]]
        if atom.element == app.element.hydrogen:
            continue
        if 'element' in row and row['element'] != atom.element.symbol:
            raise ValueError('Candidate atom element changed')
        if identity[:4] not in modeled:
            if identity[:4] not in {tuple(r['identity'][:4]) for r in context_rows} or atom.residue.name not in STANDARD_PROTEINS:
                raise ValueError('Candidate would change an unapproved retained residue')
            flanks.add(identity[:4])
        xyz = row.get('xyz_nm')
        if xyz is None:
            xyz = np.asarray(row['xyz_angstrom'], dtype=float) / 10
        xyz = np.asarray(xyz, dtype=float)
        if xyz.shape != (3,) or not np.isfinite(xyz).all():
            raise ValueError('Invalid candidate coordinate')
        updated[atom.index] = xyz
        explicit.add(identity)
    for selected in candidate_selected:
        if not all((*selected, name) in explicit for name in ('N', 'CA', 'C', 'O')):
            raise ValueError('Candidate omitted a modeled backbone atom')
    if observed_context_policy not in {'mobile_flanks', 'fixed_observed', 'preserve_outer_torsions'}:
        raise ValueError('Unknown observed-context refinement policy')
    proposed_flanks = set(flanks)
    extra_fixed = set()
    if observed_context_policy == 'fixed_observed':
        # Keep the experimental scaffold exact. This can initially strain a
        # transplanted join: record it and require all final checks after the
        # existing bounded Cartesian refinement. Restoration is never a pass.
        for atom in atoms:
            if residue_key(atom.residue) not in modeled:
                updated[atom.index] = initial[atom.index]
        flanks.clear()
    elif observed_context_policy == 'preserve_outer_torsions':
        extra_fixed = outer_torsion_anchor_indices(topology, modeled, flanks)
        for index_to_fix in extra_fixed:
            updated[index_to_fix] = initial[index_to_fix]
    # Fill only unprovided canonical sidechain atoms by preserving the seed's
    # internal sidechain frame. This is initialization, not rotamer prediction.
    transferred = []
    for residue in topology.residues():
        if residue_key(residue) not in candidate_selected | flanks:
            continue
        named = {a.name: a.index for a in residue.atoms()}
        old_origin, old_frame = frame(named, initial)
        new_origin, new_frame = frame(named, updated)
        for atom in residue.atoms():
            if atom.element == app.element.hydrogen or key(atom) in explicit:
                continue
            if atom.name in {'N', 'CA', 'C', 'O', 'OXT'}:
                raise ValueError('Partial context backbone is not supported')
            updated[atom.index] = new_origin + new_frame @ old_frame.T @ (initial[atom.index] - old_origin)
            transferred.append(list(key(atom)))
    # No observed chemical crosslink or metal-binding residue may be moved.
    metal_contacts = metal_environment(topology, initial * unit.nanometer)
    if flanks & metal_contacts['protected_keys']:
        raise ValueError('Candidate context touches an observed metal-binding residue')
    eligible = set()
    for a, b in topology.bonds():
        if a.residue is b.residue or {a.name, b.name} != {'C', 'N'}:
            continue
        akey, bkey = residue_key(a.residue), residue_key(b.residue)
        if akey in modeled and bkey not in modeled:
            eligible.add(bkey)
        if bkey in modeled and akey not in modeled:
            eligible.add(akey)
    context_scope_valid = flanks <= eligible
    for a, b in topology.bonds():
        if a.residue is not b.residue and {a.name, b.name} != {'C', 'N'}:
            if {residue_key(a.residue), residue_key(b.residue)} & flanks:
                raise ValueError('Candidate context touches a chemical crosslink')
    # Preserve exact protonation by recording every archived H-parent pair.
    neighbors = {a.index: [] for a in atoms}
    for a, b in topology.bonds():
        neighbors[a.index].append(b); neighbors[b.index].append(a)
    variants = []
    for residue in topology.residues():
        hydrogens = []
        for atom in residue.atoms():
            if atom.element != app.element.hydrogen:
                continue
            parents = neighbors[atom.index]
            if len(parents) != 1 or parents[0].residue is not residue:
                raise ValueError('Cannot preserve archived hydrogen parent identity')
            hydrogens.append((atom.name, parents[0].name))
        variants.append(hydrogens)
    forcefield, files = load_prepared_forcefield(dataset, parameters['preparation'], solvent=parameters['solvent'])
    # Record current base files separately: the old archive did not hash its
    # protein/water distribution, so replay consistency is not historical proof.
    from xml.etree import ElementTree
    pending = [dataset / p if (dataset / p).exists() else Path(app.__file__).parent / 'data' / p for p in files]
    parameter_hashes = {}
    while pending:
        path = pending.pop().resolve()
        if str(path) in parameter_hashes:
            continue
        parameter_hashes[str(path)] = sha(path)
        for include in ElementTree.parse(path).getroot().findall('Include'):
            pending.append(path.parent / include.attrib['file'])
    hashes.update(parameter_hashes)
    def system():
        return forcefield.createSystem(topology, nonbondedMethod=app.NoCutoff, constraints=None, rigidWater=False)
    native_xml = mm.XmlSerializer.serialize(system())
    (output / 'native-system.xml').write_text(native_xml)
    modeller = app.Modeller(topology, updated * unit.nanometer)
    regenerated_keys = modeled | flanks
    if not candidate.get('archive_baseline'):
        modeller.delete([a for a in modeller.topology.atoms() if a.element == app.element.hydrogen and residue_key(a.residue) in regenerated_keys])
        random.seed(914)
        platform = mm.Platform.getPlatformByName('CPU')
        platform.setPropertyDefaultValue('Threads', '2')
        modeller.addHydrogens(forcefield, variants=variants, platform=platform)
        rebuilt = {key(a): a for a in modeller.topology.atoms()}
        if set(rebuilt) != set(index):
            raise ValueError('Hydrogen reconstruction changed the archived atom/state inventory')
        before_bonds = {tuple(sorted((key(a), key(b)))) for a, b in topology.bonds()}
        after_bonds = {tuple(sorted((key(a), key(b)))) for a, b in modeller.topology.bonds()}
        if before_bonds != after_bonds:
            raise ValueError('Hydrogen reconstruction changed connectivity')
        rebuilt_xyz = np.asarray(modeller.positions.value_in_unit(unit.nanometer))
        seeded = np.asarray([rebuilt_xyz[rebuilt[key(a)].index] for a in atoms])
        fixed = [a.index for a in atoms if a.element != app.element.hydrogen or residue_key(a.residue) not in regenerated_keys]
        if not np.array_equal(seeded[fixed], updated[fixed]):
            raise ValueError('Hydrogen reconstruction moved a retained atom')
    else:
        seeded = updated
    np.savez(output / 'candidate-seed.npz', xyz_nm=seeded)
    seed_geometry = loop_geometry_report(topology, seeded * unit.nanometer, modeled | flanks)
    save(output / 'seed-geometry.json', seed_geometry)
    reference = seeded.copy()
    for atom in atoms:
        if atom.element != app.element.hydrogen:
            reference[atom.index] = initial[atom.index]
    construction_ff = AnchoredConstructionForceField(forcefield, extra_fixed) if extra_fixed else forcefield
    positions, refinement = _minimize_attempt(topology, seeded * unit.nanometer, construction_ff, modeled,
        mobile_flank_keys=flanks, reference_positions=reference * unit.nanometer,
        preserve_context_peptides=bool(flanks), max_iterations=1000)
    final = np.asarray(positions.value_in_unit(unit.nanometer))
    np.savez(output / 'refined.npz', xyz_nm=final)
    with (output / 'refined.pdb').open('w') as handle:
        app.PDBFile.writeFile(topology, positions, handle, keepIds=True)
    selected = modeled | flanks
    affected = set(selected)
    for a, b in topology.bonds():
        if a.residue is not b.residue and {a.name, b.name} == {'C', 'N'}:
            pair = {residue_key(a.residue), residue_key(b.residue)}
            if pair & selected:
                affected |= pair
    geometry = loop_geometry_report(topology, positions, selected)
    templates = PDBFixer(filename=str(dataset / 'prepared.pdb')).templates
    stereo = stereochemistry_report(topology, positions, templates)
    stereo_error = None
    try:
        require_valid_stereochemistry(stereo, 'Complete candidate')
    except ValueError as exc:
        stereo_error = str(exc)
    preliminary_rama = rama_report(topology, positions, affected)
    affected = affected_reference_keys(topology, initial, final, modeled, preliminary_rama)
    rama = rama_report(topology, positions, affected)
    outer_torsions = outside_observed_torsions(topology, initial, final, modeled, flanks)
    extra_fixed_exact = np.array_equal(initial[sorted(extra_fixed)], final[sorted(extra_fixed)])
    if not extra_fixed_exact:
        raise ValueError('A fixed observed torsion-support atom moved during construction')
    if observed_context_policy == 'preserve_outer_torsions' and not outer_torsions['all_defining_coordinates_exactly_preserved']:
        raise ValueError('The observed outer-torsion support contract was not preserved')
    peptide_checks = observed_peptides(topology, initial, final, modeled, flanks)
    peptide_preserved = all(row['same_observed_basin'] for row in peptide_checks)
    saved = app.PDBFile(str(output / 'refined.pdb'))
    saved_rama = rama_report(saved.topology, saved.positions, affected)
    saved_atoms = list(saved.topology.atoms())
    saved_xyz = np.asarray(saved.positions.value_in_unit(unit.nanometer))
    serialization_ok = [(key(a), a.element.symbol) for a in saved_atoms] == [(key(a), a.element.symbol) for a in atoms]
    serialization_ok &= bond_keys(saved.topology) == bond_keys(topology)
    serialization_ok &= saved_xyz.shape == final.shape and bool(np.max(np.linalg.norm(saved_xyz - final, axis=1)) <= .00009)
    serialization_ok &= {tuple(r['residue']): r['classification'] for r in rama['rows']} == {tuple(r['residue']): r['classification'] for r in saved_rama['rows']}
    saved_geometry = loop_geometry_report(saved.topology, saved.positions, selected)
    saved_stereo = stereochemistry_report(saved.topology, saved.positions, templates)
    try:
        require_valid_stereochemistry(saved_stereo, 'Serialized complete candidate')
    except ValueError:
        serialization_ok = False
    serialization_ok &= saved_geometry['accepted']
    saved_peptides = observed_peptides(topology, initial, saved_xyz, modeled, flanks) if len(saved_atoms) == len(atoms) else []
    serialization_ok &= len(saved_peptides) == len(peptide_checks) and all(row['same_observed_basin'] for row in saved_peptides)
    unchanged_indices = [a.index for a in atoms if a.element != app.element.hydrogen and residue_key(a.residue) not in selected]
    unchanged = bool(np.array_equal(initial[unchanged_indices], final[unchanged_indices]))
    moved_indices = [a.index for a in atoms if a.element != app.element.hydrogen and residue_key(a.residue) in flanks]
    max_displacement = float(max(np.linalg.norm(final[moved_indices] - initial[moved_indices], axis=1), default=0))
    native_unchanged = native_xml == mm.XmlSerializer.serialize(system())
    checks = {'identity': True, 'retained_environment': unchanged and not geometry['gross_collisions'],
              'geometry': geometry['accepted'] and peptide_preserved, 'stereochemistry': stereo_error is None,
              'backbone_reference': rama['all_selected_scored_without_outliers'],
              'observed_displacement': max_displacement <= .1 and context_scope_valid,
              'serialized_output': bool(serialization_ok), 'parameter_integrity': native_unchanged}
    assert set(checks) == REQUIRED_CHECKS
    reasons = list(geometry['errors'])
    if stereo_error: reasons.append(stereo_error)
    if not peptide_preserved: reasons.append('An observed peptide changed its cis/trans basin.')
    if rama['outliers']: reasons.append('Affected backbone reference contains outliers.')
    if not context_scope_valid: reasons.append('Proposed context extends beyond immediate eligible flanks.')
    if max_displacement > .1: reasons.append('Observed heavy-atom displacement exceeds1 Å.')
    if not serialization_ok: reasons.append('Saved-file identity/reference mismatch.')
    for name, data in [('geometry', geometry), ('stereo', stereo), ('reference', rama), ('saved-reference', saved_rama),
                       ('saved-geometry', saved_geometry), ('saved-stereo', saved_stereo),
                       ('observed-peptides', peptide_checks), ('saved-observed-peptides', saved_peptides),
                       ('outside-observed-torsions', outer_torsions), ('preliminary-reference', preliminary_rama), ('refinement', refinement)]:
        save(output / (name + '.json'), data)
    if hashes != {name: sha(Path(name)) for name in hashes}:
        raise ValueError('A frozen source changed during candidate validation')
    result = {'checks': checks, 'reasons': reasons, 'accepted_by_recorded_checks': all(checks.values()),
              'modeled_residue_keys': sorted(modeled), 'candidate_modeled_residue_keys': sorted(candidate_selected),
              'observed_context_residue_keys': sorted(flanks), 'atom_count': len(atoms),
              'proposed_observed_context_residue_keys': sorted(proposed_flanks),
              'observed_context_policy': observed_context_policy,
              'additional_fixed_observed_atoms': [list(key(atoms[i])) for i in sorted(extra_fixed)],
              'additional_fixed_observed_atoms_exactly_preserved': bool(extra_fixed_exact),
              'reference_selection': 'All new definitions and every actually changed phi/psi/preceding-omega definition. Fully unchanged source-supported external rows are excluded by exact defining-coordinate equality, never by outlier classification.',
              'seed_geometry_accepted': seed_geometry['accepted'],
              'maximum_observed_heavy_displacement_angstrom': max_displacement * 10,
              'sidechain_frame_initialization_atoms': transferred,
              'hydrogen_states': 'Exact archived H names and parent bonds; only changed/modelled-residue H regenerated.',
              'native_system_sha256': hashlib.sha256(native_xml.encode()).hexdigest(),
              'scaffold_binding': scaffold_binding,
              'source_sha256': hashes, 'parameter_files': files, 'app_ready': False,
              'parameter_integrity_scope': 'Archived ligand/modified-residue checksum checks, frozen current base XML files, and fresh native-system replay with no construction-force leakage. Historical protein/water file equivalence is unverified.',
              'physical_model_validated': False, 'new_md_or_qm': False,
              'artifacts_sha256': {str(output / name): sha(output / name) for name in ['refined.npz', 'refined.pdb', 'native-system.xml']},
              'topology_file': str(inputs[2]),
              'scope': 'Private archived complete-complex refinement and saved-file screening; no full preparation rerun or release admission.'}
    save(output / 'result.json', result)
    print(json.dumps({k: result[k] for k in ['checks', 'reasons', 'accepted_by_recorded_checks', 'maximum_observed_heavy_displacement_angstrom']}), flush=True)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    for name in ['dataset', 'candidate', 'output']:
        parser.add_argument('--' + name, required=True, type=Path)
    parser.add_argument('--observed-context-policy', choices=['mobile_flanks', 'fixed_observed', 'preserve_outer_torsions'], default='mobile_flanks')
    args = parser.parse_args()
    run(args.dataset.resolve(), args.candidate.resolve(), args.output.resolve(), args.observed_context_policy)
