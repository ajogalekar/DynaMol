"""Research-only coherent ProMod3 context seed under unchanged final gates.

The native context may start outside the displacement cap. Acceptance still
requires the original one-angstrom cap and independent geometry/stereo checks.
"""
import argparse
import json
import os
from pathlib import Path
import sys

import numpy as np
import openmm as mm
from openmm import app, unit
from pdbfixer import PDBFixer

from prepare_adduct import digest


def run(repo, original, context, output, backbone_construction=False, environment=None,
        preserve_backbone_conformation=False, preserve_context_peptides=False, seed_pdb=None):
    output.mkdir(exist_ok=False)
    os.environ['DYNAMOL_CPU_THREADS'] = '2'
    sys.path.insert(0, str(repo))
    from backend.loop_refinement import _minimize_attempt
    from backend.loop_geometry import _RADII, loop_geometry_report
    from backend.preparation_worker import atom_key, stereochemistry_report
    from backend.residue_identity import residue_key

    evidence = json.loads((context/'output.json').read_text())
    config = json.loads((context/'input.json').read_text())
    selected = {tuple([r['chain'], r['resid'], r['insertion_code'], r['residue']])
        for chain in config['chains'] for r in chain['residues'] if not r['observed']}
    source = PDBFixer(filename=str(original/'before-refinement.pdb'))
    source_paths = [original/'before-refinement.pdb', context/'output.json', context/'input.json']
    if seed_pdb is not None:
        source_paths.append(seed_pdb)
    external = []
    external_records = []
    if environment is not None:
        source_paths.append(environment)
        retained = app.PDBFile(str(environment))
        for atom in retained.topology.atoms():
            if atom.element is None:
                raise ValueError('Retained environment contains an unspecified element')
            if atom.element.symbol in {'H', 'D'}:
                continue
            xyz = np.asarray(retained.positions[atom.index].value_in_unit(unit.nanometer))
            radius = _RADII.get(atom.element.symbol, .170)
            external.append((xyz, radius))
            external_records.append({'identity': atom_key(atom), 'xyz_nm': xyz.tolist(),
                                     'steric_radius_nm': radius})
        if len({tuple(row['identity']) for row in external_records}) != len(external_records):
            raise ValueError('Ambiguous retained environment atom identity')
    source_hashes = {str(p): digest(p) for p in source_paths}
    original_atoms = {atom_key(a): np.array(source.positions[a.index].value_in_unit(unit.nanometer))
        for a in source.topology.atoms() if a.element != app.element.hydrogen}
    modeller = app.Modeller(source.topology, source.positions)
    modeller.delete([a for a in modeller.topology.atoms() if a.element == app.element.hydrogen])
    atoms = list(modeller.topology.atoms())
    indices = {atom_key(a): a.index for a in atoms}
    if len(indices) != len(atoms):
        raise ValueError('Ambiguous source atom identity')
    if set(indices) & {tuple(row['identity']) for row in external_records}:
        raise ValueError('Retained environment duplicates protein atom identities')
    seeded = np.asarray(modeller.positions.value_in_unit(unit.nanometer)).copy()
    flanks = set()
    changed = []
    for row in evidence['research_context_atoms']:
        key = tuple(row['identity'])
        if key not in indices:
            raise ValueError('Native context introduced an unknown observed atom')
        if np.max(abs(np.asarray(row['original_xyz_nm'])-original_atoms[key])) > .00011:
            raise ValueError('Native context is not bound to the original observed coordinates')
        if row['displacement_nm'] <= 1e-7:
            continue
        if key[1] == '797':
            raise ValueError('Native context would move the covalent attachment residue')
        seeded[indices[key]] = row['xyz_nm']
        flanks.add(key[:4])
        changed.append(row)
    for row in evidence['atoms']:
        key = tuple(row['identity'])
        if key[:4] not in selected or key not in indices:
            raise ValueError('Native missing-loop identity mismatch')
        seeded[indices[key]] = row['xyz_nm']
    modeller.positions = seeded*unit.nanometer
    ff = app.ForceField('amber14/protein.ff14SB.xml', 'amber14/tip3p.xml')
    platform = mm.Platform.getPlatformByName('CPU')
    platform.setPropertyDefaultValue('Threads', '2')
    if seed_pdb is None:
        modeller.addHydrogens(ff, pH=7.4, platform=platform)
    else:
        saved = app.PDBFile(str(seed_pdb))
        saved_heavy = {atom_key(a): a for a in saved.topology.atoms() if a.element != app.element.hydrogen}
        if saved_heavy.keys() != indices.keys():
            raise ValueError('Saved replay seed changed the expected heavy-atom inventory')
        saved_xyz = np.asarray(saved.positions.value_in_unit(unit.nanometer))
        for key, atom in saved_heavy.items():
            if atom.element != atoms[indices[key]].element or np.max(abs(saved_xyz[atom.index]-seeded[indices[key]])) > .00011:
                raise ValueError('Saved replay seed differs from the recorded native context beyond PDB precision')
        if len({atom_key(a) for a in saved.topology.atoms()}) != saved.topology.getNumAtoms():
            raise ValueError('Saved replay seed has ambiguous atom identities')
        modeller = app.Modeller(saved.topology, saved.positions)
    reference = np.asarray(modeller.positions.value_in_unit(unit.nanometer)).copy()
    for atom in modeller.topology.atoms():
        if atom.element != app.element.hydrogen:
            reference[atom.index] = original_atoms[atom_key(atom)]
    before = loop_geometry_report(modeller.topology, modeller.positions, selected|flanks,
                                  environment_positions=external)
    with (output/'coherent-context-seed.pdb').open('w') as handle:
        app.PDBFile.writeFile(modeller.topology, modeller.positions, handle, keepIds=True)
    (output/'seed-geometry.json').write_text(json.dumps(before, indent=2)+'\n')
    (output/'context-displacements.json').write_text(json.dumps(changed, indent=2)+'\n')
    construction_angles = []
    conformation_targets = []
    seed_rama = None
    refinement_ff = ff
    if backbone_construction:
        class ConstructionForceField:
            def createSystem(self, topology, **kwargs):
                system = ff.createSystem(topology, **kwargs)
                extra = mm.CustomAngleForce('0.5*k*(theta-target)^2')
                extra.addPerAngleParameter('k')
                extra.addPerAngleParameter('target')
                atoms = list(topology.atoms())
                for force in system.getForces():
                    if not isinstance(force, mm.HarmonicAngleForce):
                        continue
                    for i in range(force.getNumAngles()):
                        a, b, c, target, strength = force.getAngleParameters(i)
                        triple = [atoms[j] for j in [a, b, c]]
                        if not all(atom.name in {'N', 'CA', 'C'} for atom in triple):
                            continue
                        if not any(residue_key(atom.residue) in selected|flanks for atom in triple):
                            continue
                        k = 4*strength.value_in_unit(unit.kilojoule_per_mole/unit.radian**2)
                        theta = target.value_in_unit(unit.radian)
                        extra.addAngle(a, b, c, [k, theta])
                        construction_angles.append({'atoms': [a, b, c],
                            'target_radians_from_native_forcefield': theta,
                            'additional_strength_kj_mol_rad2': k})
                system.addForce(extra)
                return system
        refinement_ff = ConstructionForceField()
    if preserve_backbone_conformation:
        from rama_reference import rama_report
        seed_rama = rama_report(modeller.topology, modeller.positions, selected|flanks)
        (output/'seed-ramachandran.json').write_text(json.dumps(seed_rama, indent=2)+'\n')
        if seed_rama['unavailable']:
            raise ValueError('Conformation protection requires complete native Ramachandran scores')
        # Protect only statistically allowed seed conformations. An existing
        # seed outlier is left free to improve and must pass the final screen.
        for row in seed_rama['rows']:
            if row['classification'] == 'outlier':
                continue
            for kind in ['phi', 'psi']:
                conformation_targets.append({'residue': row['residue'], 'kind': kind,
                    'atoms': row[kind+'_atoms'], 'target_radians': np.deg2rad(row[kind+'_degrees']),
                    'seed_rama_score': row['score'], 'construction_strength_kj_mol': 1000.0})
        base_conformation_ff = refinement_ff
        class ConformationConstructionForceField:
            def createSystem(self, topology, **kwargs):
                system = base_conformation_ff.createSystem(topology, **kwargs)
                tether = mm.CustomTorsionForce('k*(1-cos(theta-target))')
                tether.addPerTorsionParameter('k')
                tether.addPerTorsionParameter('target')
                for row in conformation_targets:
                    tether.addTorsion(*row['atoms'], [row['construction_strength_kj_mol'], row['target_radians']])
                system.addForce(tether)
                return system
        refinement_ff = ConformationConstructionForceField()
    environment_pair_count = 0
    if external:
        base_refinement_ff = refinement_ff
        class EnvironmentConstructionForceField:
            def createSystem(self, topology, **kwargs):
                nonlocal environment_pair_count
                system = base_refinement_ff.createSystem(topology, **kwargs)
                # Fixed geometric obstacles during construction only. These
                # radii and forces are NOT ligand/ion force-field parameters.
                obstacle = mm.CustomCompoundBondForce(1,
                    '0.5*k*max(0,rmin-r)^2; '
                    'r=sqrt((x1-x0)^2+(y1-y0)^2+(z1-z0)^2+1e-12)')
                for name in ['k', 'rmin', 'x0', 'y0', 'z0']:
                    obstacle.addPerBondParameter(name)
                for atom in topology.atoms():
                    if atom.element is None:
                        raise ValueError('Construction particle has an unspecified element')
                    if atom.element.symbol in {'H', 'D'} or residue_key(atom.residue) not in selected|flanks:
                        continue
                    radius = _RADII.get(atom.element.symbol, .170)
                    for xyz, other_radius in external:
                        # Include every environment pair, with no stale neighbor
                        # list when a modeled atom moves during minimization.
                        obstacle.addBond([atom.index], [100000.0, .9*(radius+other_radius), *xyz])
                environment_pair_count = obstacle.getNumBonds()
                system.addForce(obstacle)
                return system
        refinement_ff = EnvironmentConstructionForceField()
    positions, refinement = _minimize_attempt(modeller.topology, modeller.positions, refinement_ff, selected,
        mobile_flank_keys=flanks, reference_positions=reference*unit.nanometer, max_iterations=1000,
        preserve_context_peptides=preserve_context_peptides)
    geometry = loop_geometry_report(modeller.topology, positions, selected|flanks,
                                    environment_positions=external)
    stereo = stereochemistry_report(modeller.topology, positions, source.templates)
    final_rama = None
    if preserve_backbone_conformation:
        final_rama = rama_report(modeller.topology, positions, selected|flanks)
        (output/'ramachandran.json').write_text(json.dumps(final_rama, indent=2)+'\n')
    with (output/'refined-candidate.pdb').open('w') as handle:
        app.PDBFile.writeFile(modeller.topology, positions, handle, keepIds=True)
    (output/'geometry.json').write_text(json.dumps(geometry, indent=2)+'\n')
    (output/'stereochemistry.json').write_text(json.dumps(stereo, indent=2)+'\n')
    refinement['flank_relaxation']['scope'] = 'Research-only exact native changed-context residues; original one-angstrom cap unchanged. App integration not performed.'
    refinement['backbone_angle_construction'] = {'enabled': backbone_construction,
        'additional_weight_relative_to_native': 4 if backbone_construction else 0,
        'angles': construction_angles, 'exported_to_simulation': False,
        'scope': 'Temporary local geometry construction with exact native equilibrium angles; original MD force-field parameters and independent acceptance limits are unchanged.'}
    refinement['retained_environment_construction'] = {'enabled': bool(external),
        'heavy_atoms': external_records, 'pair_count': environment_pair_count,
        'radial_strength_kj_mol_nm2': 100000.0, 'target_radius_sum_factor': .9,
        'coordinates_fixed': True, 'exported_to_simulation': False,
        'scope': 'Temporary steric obstacles for all retained heavy atoms against modeled/context heavy atoms. No electrostatics or environment force-field accuracy is implied; independent complete-complex mechanics and dynamics are still required.'}
    refinement['backbone_conformation_protection'] = {'enabled': preserve_backbone_conformation,
        'targets': conformation_targets, 'exported_to_simulation': False,
        'unprotected_seed_outliers': seed_rama['outliers'] if seed_rama else [],
        'scope': 'Research construction tethers retain allowed native-fragment phi/psi neighborhoods during anchor/environment refinement. Strength matches the existing final peptide construction strength; it is not a physical torsion parameter or a fitted Ramachandran energy. Independent final native CCTBX scoring is additional to all original gates.'}
    (output/'refinement.json').write_text(json.dumps(refinement, indent=2)+'\n')
    # Reuse the actual app stereo admission rather than interpreting its report.
    from backend.preparation_worker import require_valid_stereochemistry
    stereo_error = None
    try:
        require_valid_stereochemistry(stereo, 'Context-seeded research candidate')
    except Exception as exc:
        stereo_error = str(exc)
    displacement = refinement['flank_relaxation']['maximum_displacement_nm']
    report = {'stage': 'coherent_context_seed_checked', 'app_ready': False,
        'physical_model_validated': False, 'full_environment_validated': False,
        'geometry_accepted': geometry['accepted'], 'stereo_error': stereo_error,
        'maximum_observed_context_displacement_nm': displacement,
        'unchanged_displacement_limit_nm': .1,
        'passed_local_checks': (geometry['accepted'] and stereo_error is None and displacement <= .1
                               and (final_rama is None or final_rama['all_selected_scored_without_outliers'])),
        'ramachandran_reference_checked': final_rama is not None,
        'ramachandran_no_selected_outliers': final_rama['all_selected_scored_without_outliers'] if final_rama else None,
        'modeled_residues': sorted(selected), 'native_changed_context_residues': sorted(flanks),
        'seed_geometry_error_count': len(before['errors']), 'final_geometry_errors': geometry['errors'],
        'retained_environment_geometry_checked': bool(external),
        'replayed_seed_pdb': str(seed_pdb) if seed_pdb is not None else None,
        'source_sha256': source_hashes}
    if source_hashes != {str(p): digest(p) for p in source_paths}:
        raise ValueError('Source changed during construction')
    (output/'result.json').write_text(json.dumps(report, indent=2)+'\n')
    print(json.dumps(report, indent=2))


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    for name in ['repo', 'original', 'context', 'output']:
        parser.add_argument('--'+name, type=Path, required=True)
    parser.add_argument('--backbone-construction', action='store_true')
    parser.add_argument('--environment', type=Path)
    parser.add_argument('--preserve-backbone-conformation', action='store_true')
    parser.add_argument('--preserve-context-peptides', action='store_true')
    parser.add_argument('--seed-pdb', type=Path)
    args = parser.parse_args()
    try:
        run(args.repo, args.original, args.context, args.output, args.backbone_construction, args.environment,
            args.preserve_backbone_conformation, args.preserve_context_peptides, args.seed_pdb)
    except Exception as exc:
        if args.output.is_dir():
            (args.output/'failure.json').write_text(json.dumps({'type': type(exc).__name__, 'error': str(exc)}, indent=2)+'\n')
        raise
