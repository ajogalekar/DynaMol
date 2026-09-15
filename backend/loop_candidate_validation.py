"""Complete static loop checks, isolated from default preparation until qualified.

The caller supplies a verified fresh preparation snapshot and its force field.
This module does not fetch structures, choose chemistry, fit parameters, publish
datasets or turn a generator score into an acceptance decision.
"""
from collections import Counter
from copy import deepcopy
import hashlib
import json
from pathlib import Path
import random
from xml.etree import ElementTree as ET

import numpy as np
import openmm as mm
from openmm import app, unit

from .loop_context import (admit_loop_proposal, affected_reference_keys,
                           outside_observed_torsions, AnchoredConstructionForceField)
from .loop_geometry import loop_geometry_report, _dihedral
from .loop_refinement import _minimize_attempt
from .loop_search import REQUIRED_CHECKS
from .residue_identity import STANDARD_PROTEINS, residue_key


def atom_key(atom):
    return (*residue_key(atom.residue), atom.name)


def bond_keys(topology):
    return {tuple(sorted((atom_key(a), atom_key(b)))) for a, b in topology.bonds()}


def physical_parameter_comparison(first_xml, second_xml):
    """Compare undirected harmonic bonds, preserving parameters/multiplicity."""
    first, second = ET.fromstring(first_xml), ET.fromstring(second_xml)
    forces1, forces2 = first.find('Forces'), second.find('Forces')
    if forces1 is None or forces2 is None:
        raise ValueError('Physical System serialization is incomplete.')
    hb1 = [f for f in forces1 if f.attrib.get('type') == 'HarmonicBondForce']
    hb2 = [f for f in forces2 if f.attrib.get('type') == 'HarmonicBondForce']
    if len(hb1) != 1 or len(hb2) != 1:
        raise ValueError('Expected exactly one physical harmonic-bond force.')
    canonical = lambda e: ET.canonicalize(ET.tostring(e, encoding='unicode'), strip_text=True)
    b1, b2 = hb1[0].find('Bonds'), hb2[0].find('Bonds')
    if b1 is None or b2 is None:
        raise ValueError('Physical harmonic-bond terms are missing.')
    def canonical_bond(term):
        term = deepcopy(term)
        if term.tag != 'Bond':
            raise ValueError('Unexpected child in physical harmonic-bond terms.')
        try:
            endpoints = sorted(int(term.attrib[key]) for key in ('p1', 'p2'))
        except (KeyError, TypeError, ValueError) as exc:
            raise ValueError('Malformed physical harmonic-bond particle indices.') from exc
        if endpoints[0] < 0:
            raise ValueError('Negative physical harmonic-bond particle index.')
        # HarmonicBondForce depends on the distance between the endpoints;
        # reversing their order leaves this term exactly the same. No other
        # force (in particular angle/torsion orientation) is normalized here.
        term.set('p1', str(endpoints[0])); term.set('p2', str(endpoints[1]))
        return canonical(term)
    equal_bonds = Counter(canonical_bond(term) for term in b1) == Counter(canonical_bond(term) for term in b2)
    count = len(b1)
    # Only term order and undirected endpoints are ignored; retain attributes and
    # any other child data in the exact remainder comparison.
    b1[:] = []; b2[:] = []
    equal_rest = canonical(first) == canonical(second)
    return {'harmonic_bond_parameter_multisets_equal': equal_bonds, 'bond_terms': count,
            'other_system_fields_and_terms_exactly_equal': equal_rest,
            'all_atom_indexed_parameters_equal': equal_bonds and equal_rest}


def _observed_peptides(topology, before, after, modeled, flanks):
    rows = []
    for a, b in topology.bonds():
        if a.residue is b.residue or {a.name, b.name} != {'C', 'N'}:
            continue
        pair = {residue_key(a.residue), residue_key(b.residue)}
        if pair & modeled or not pair & flanks:
            continue
        left, right = (a, b) if a.name == 'C' else (b, a)
        names = [{x.name: x.index for x in r.atoms()} for r in (left.residue, right.residue)]
        indices = [names[0]['CA'], left.index, right.index, names[1]['CA']]
        initial, final = _dihedral(before[indices]), _dihedral(after[indices])
        valid = initial is not None and final is not None
        valid = valid and abs(abs(initial) - 90) > 1e-6 and (abs(initial) < 90) == (abs(final) < 90)
        rows.append({'residues': sorted(pair), 'source_omega_degrees': initial,
                     'final_omega_degrees': final, 'same_observed_basin': bool(valid)})
    return rows


def _rebuild_hydrogens(topology, reference, seeded, forcefield, regenerated):
    atoms = list(topology.atoms())
    index = {atom_key(a): a.index for a in atoms}
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
                raise ValueError('Cannot preserve the prepared hydrogen parent identity.')
            hydrogens.append((atom.name, parents[0].name))
        variants.append(hydrogens)
    modeller = app.Modeller(topology, seeded * unit.nanometer)
    modeller.delete([a for a in modeller.topology.atoms()
                     if a.element == app.element.hydrogen and residue_key(a.residue) in regenerated])
    random.seed(914)
    platform = mm.Platform.getPlatformByName('CPU')
    platform.setPropertyDefaultValue('Threads', '2')
    modeller.addHydrogens(forcefield, variants=variants, platform=platform)
    rebuilt = {atom_key(a): a for a in modeller.topology.atoms()}
    if len(rebuilt) != len(atoms) or set(rebuilt) != set(index):
        raise ValueError('Hydrogen reconstruction changed the prepared atom/state inventory.')
    if bond_keys(topology) != bond_keys(modeller.topology):
        raise ValueError('Hydrogen reconstruction changed prepared connectivity.')
    xyz = np.asarray(modeller.positions.value_in_unit(unit.nanometer))
    result = np.asarray([xyz[rebuilt[atom_key(a)].index] for a in atoms])
    fixed = [a.index for a in atoms if a.element != app.element.hydrogen
             or residue_key(a.residue) not in regenerated]
    if not np.array_equal(result[fixed], seeded[fixed]):
        raise ValueError('Hydrogen reconstruction moved a retained atom.')
    return result


def validate_complete_loop(topology, reference_positions, forcefield, modeled_keys,
                           proposal, output, *, source_sha256, templates,
                           protected_residues=(), newly_defined_backbone_residue_keys=(),
                           preserve_seed_torsions=False):
    """Refine one complete proposal and return its eight-check result record.

    Run in a supervised child: native hydrogen building/minimization can block.
    ``reference_positions`` is the verified common-preparation reference, never
    a previous failed proposal. The caller binds the original source, explicit
    state decisions, all parameter bytes and final output bundle separately.
    """
    from .loop_backbone_reference import rama_report
    from .preparation_worker import stereochemistry_report, require_valid_stereochemistry
    if type(preserve_seed_torsions) is not bool:
        raise ValueError('Seed torsion preservation requires an explicit boolean')
    output = Path(output).resolve()
    output.mkdir(parents=True, exist_ok=False)
    atoms = list(topology.atoms())
    reference = np.asarray(reference_positions.value_in_unit(unit.nanometer), dtype=float)
    if reference.shape != (len(atoms), 3) or not np.isfinite(reference).all():
        raise ValueError('The complete preparation reference is malformed.')
    modeled = set(map(tuple, modeled_keys))
    allowed = set()
    for a, b in topology.bonds():
        if a.residue is b.residue or {a.name, b.name} != {'C', 'N'}:
            continue
        for first, second in ((a, b), (b, a)):
            if residue_key(first.residue) in modeled and residue_key(second.residue) not in modeled:
                if second.residue.name in STANDARD_PROTEINS:
                    allowed.add(residue_key(second.residue))
    admitted = admit_loop_proposal(proposal, topology, reference, modeled, allowed,
                                  source_sha256=source_sha256,
                                  protected_residues=protected_residues)
    flanks = set(admitted.mobile_flank_keys)
    extra_fixed = set(admitted.extra_fixed_indices)
    seeded = _rebuild_hydrogens(topology, reference, admitted.seeded_xyz_nm,
                                forcefield, modeled | flanks)
    seed_targets = None
    if preserve_seed_torsions:
        seed_reference = rama_report(topology, seeded * unit.nanometer, modeled | flanks)
        (output / 'seed-backbone-reference.json').write_text(json.dumps(seed_reference, indent=2, allow_nan=False) + '\n')
        if not seed_reference['all_selected_scored_without_outliers']:
            raise ValueError('Torsion preservation requires every modeled and flank seed definition to pass the backbone reference')
        seed_targets = [{'residue': row['residue'], 'kind': kind, 'atoms': row[kind + '_atoms'],
                         'target_radians': float(np.deg2rad(row[kind + '_degrees']))}
                        for row in seed_reference['rows'] for kind in ('phi', 'psi')]

    def native_system(top=topology):
        return forcefield.createSystem(top, nonbondedMethod=app.NoCutoff,
                                       constraints=None, rigidWater=False)

    native_xml = mm.XmlSerializer.serialize(native_system())
    (output / 'native-system.xml').write_text(native_xml)
    np.savez(output / 'candidate-seed.npz', xyz_nm=seeded)
    construction_ff = AnchoredConstructionForceField(forcefield, extra_fixed) if extra_fixed else forcefield
    restraint_reference = seeded.copy()
    heavy = [a.index for a in atoms if a.element != app.element.hydrogen]
    restraint_reference[heavy] = reference[heavy]
    positions, refinement = _minimize_attempt(
        topology, seeded * unit.nanometer, construction_ff, modeled,
        mobile_flank_keys=flanks, reference_positions=restraint_reference * unit.nanometer,
        preserve_context_peptides=bool(flanks), max_iterations=1000,
        **({'seed_torsion_targets': seed_targets} if preserve_seed_torsions else {}))
    final = np.asarray(positions.value_in_unit(unit.nanometer))
    np.savez(output / 'refined.npz', xyz_nm=final)
    with (output / 'refined.pdb').open('w') as stream:
        app.PDBFile.writeFile(topology, positions, stream, keepIds=True)
    with (output / 'topology.cif').open('w') as stream:
        app.PDBxFile.writeFile(topology, positions, stream, keepIds=True)
    selected = modeled | flanks
    affected = set(selected)
    for a, b in topology.bonds():
        if a.residue is not b.residue and {a.name, b.name} == {'C', 'N'}:
            pair = {residue_key(a.residue), residue_key(b.residue)}
            if pair & selected:
                affected |= pair
    new_backbone = set(map(tuple, newly_defined_backbone_residue_keys))
    affected |= new_backbone
    neighboring_scope = set(affected)
    preliminary = rama_report(topology, positions, neighboring_scope)
    affected = affected_reference_keys(topology, reference, final, modeled, preliminary) | new_backbone
    rama = rama_report(topology, positions, affected)
    outer = outside_observed_torsions(topology, reference, final, modeled, flanks)
    if not np.array_equal(reference[sorted(extra_fixed)], final[sorted(extra_fixed)]):
        raise ValueError('A fixed outer torsion-support atom moved during construction.')
    if not outer['all_defining_coordinates_exactly_preserved']:
        raise ValueError('An observed external backbone torsion definition changed.')
    geometry = loop_geometry_report(topology, positions, selected)
    stereo = stereochemistry_report(topology, positions, templates)
    stereo_error = None
    try:
        require_valid_stereochemistry(stereo, 'Complete loop candidate')
    except ValueError as exc:
        stereo_error = str(exc)
    peptides = _observed_peptides(topology, reference, final, modeled, flanks)
    saved = app.PDBFile(str(output / 'refined.pdb'))
    saved_xyz = np.asarray(saved.positions.value_in_unit(unit.nanometer))
    saved_atoms = list(saved.topology.atoms())
    identity_saved = [(atom_key(a), a.element.symbol) for a in saved_atoms] == [
        (atom_key(a), a.element.symbol) for a in atoms]
    identity_saved &= bond_keys(saved.topology) == bond_keys(topology)
    if not identity_saved:
        raise ValueError('Serialized complete topology changed atom identities or bonds.')
    saved_preliminary = rama_report(saved.topology, saved.positions, neighboring_scope)
    saved_affected = affected_reference_keys(saved.topology, reference, saved_xyz, modeled, saved_preliminary) | new_backbone
    saved_rama = rama_report(saved.topology, saved.positions, saved_affected)
    saved_outer = outside_observed_torsions(saved.topology, reference, saved_xyz, modeled, flanks)
    saved_geometry = loop_geometry_report(saved.topology, saved.positions, selected)
    saved_stereo = stereochemistry_report(saved.topology, saved.positions, templates)
    saved_peptides = _observed_peptides(topology, reference, saved_xyz, modeled, flanks)
    saved_valid = saved_xyz.shape == final.shape and bool(np.max(np.linalg.norm(saved_xyz - final, axis=1)) <= .00009)
    saved_valid &= saved_geometry['accepted'] and saved_rama['all_selected_scored_without_outliers']
    saved_valid &= saved_outer['all_defining_coordinates_exactly_preserved']
    saved_valid &= np.array_equal(reference[sorted(extra_fixed)], saved_xyz[sorted(extra_fixed)])
    saved_valid &= {tuple(r['residue']): r['classification'] for r in rama['rows']} == {
        tuple(r['residue']): r['classification'] for r in saved_rama['rows']}
    saved_valid &= len(saved_peptides) == len(peptides) and all(r['same_observed_basin'] for r in saved_peptides)
    try:
        require_valid_stereochemistry(saved_stereo, 'Saved complete loop candidate')
    except ValueError:
        saved_valid = False
    parameter_reload = physical_parameter_comparison(native_xml, mm.XmlSerializer.serialize(native_system(saved.topology)))
    fixed = [a.index for a in atoms if a.element != app.element.hydrogen and residue_key(a.residue) not in selected]
    moved = [a.index for a in atoms if a.element != app.element.hydrogen and residue_key(a.residue) in flanks]
    unchanged = bool(np.array_equal(reference[fixed], final[fixed]))
    maximum = float(max(np.linalg.norm(final[moved] - reference[moved], axis=1), default=0))
    saved_maximum = float(max(np.linalg.norm(saved_xyz[moved] - reference[moved], axis=1), default=0))
    saved_valid &= saved_maximum <= .1
    checks = {
        'identity': True, 'retained_environment': unchanged and not geometry['gross_collisions'],
        'geometry': geometry['accepted'] and all(r['same_observed_basin'] for r in peptides),
        'stereochemistry': stereo_error is None,
        'backbone_reference': rama['all_selected_scored_without_outliers'],
        'observed_displacement': maximum <= .1,
        'serialized_output': bool(saved_valid),
        'parameter_integrity': native_xml == mm.XmlSerializer.serialize(native_system())
                               and parameter_reload['all_atom_indexed_parameters_equal'],
    }
    if set(checks) != REQUIRED_CHECKS:
        raise ValueError('Complete loop validation omitted a required check.')
    reports = {'geometry': geometry, 'stereochemistry': stereo, 'reference': rama,
               'saved-geometry': saved_geometry, 'saved-stereochemistry': saved_stereo,
               'saved-reference': saved_rama, 'refinement': refinement,
               'outside-observed-torsions': outer, 'parameter-roundtrip': parameter_reload,
               'saved-outside-observed-torsions': saved_outer,
               'observed-peptides': peptides, 'saved-observed-peptides': saved_peptides,
               'proposal-admission': dict(admitted.provenance)}
    for name, report in reports.items():
        (output / (name + '.json')).write_text(json.dumps(report, indent=2, allow_nan=False) + '\n')
    reasons = [name for name, passed in checks.items() if not passed]
    report = {'checks': checks, 'reasons': reasons, 'accepted': not reasons,
              'modeled_residue_keys': sorted(modeled), 'observed_context_residue_keys': sorted(flanks),
              'maximum_observed_heavy_displacement_angstrom': maximum * 10,
              'maximum_saved_observed_heavy_displacement_angstrom': saved_maximum * 10,
              'newly_defined_backbone_residue_keys': sorted(new_backbone),
              'displacement_baseline': 'Verified common preparation reference; original-source movement is bound and reported separately by the caller.',
              'source_sha256': source_sha256, 'atom_count': len(atoms),
              'artifacts_sha256': {p.name: hashlib.sha256(p.read_bytes()).hexdigest()
                                   for p in output.iterdir() if p.is_file()},
              'scope': 'Complete static checks only; caller must bind source/state/parameter snapshot and verify all sidecars before publication.',
              'app_ready': False, 'physical_model_validated': False}
    (output / 'result.json').write_text(json.dumps(report, indent=2, allow_nan=False) + '\n')
    return report
