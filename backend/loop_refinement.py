"""Bounded loop relaxation, with a recorded restrained-flank fallback.

This creates starting geometry before solvation/dynamics. A finite, lower
energy is not acceptance: independent stereo and geometry gates must pass.
"""
from __future__ import annotations

import time
import numpy as np
import openmm as mm
from openmm import app, unit

from . import config
from .residue_identity import residue_key


def _validated_seed_torsion_targets(topology, xyz, modeled_keys, flank_keys, targets):
    """Validate caller-qualified phi/psi targets, without loading reference tables.

    The caller must first qualify every modeled/flank seed with its independent
    backbone reference. This checks exact topology, scope and seed angles; it
    does not itself certify Ramachandran quality or a native conformation.
    """
    from .loop_geometry import _dihedral
    from .residue_identity import STANDARD_PROTEINS, backbone

    selected = set(modeled_keys) | set(flank_keys)
    if not isinstance(targets, (list, tuple)) or not targets or len(targets) != 2 * len(selected):
        raise ValueError('Seed torsion targets require exactly phi and psi for every modeled/flank residue.')
    atoms = list(topology.atoms())
    xyz = np.asarray(xyz, dtype=float)
    if xyz.shape != (len(atoms), 3) or not np.isfinite(xyz).all():
        raise ValueError('Seed torsion targets require finite coordinates matching the topology.')
    residues = list(topology.residues())
    by_key = {}
    for residue in residues:
        key = residue_key(residue)
        if key in by_key:
            raise ValueError('Seed torsion targets require unambiguous residue identities.')
        by_key[key] = residue
    if not selected or not selected <= by_key.keys() or set(modeled_keys) & set(flank_keys):
        raise ValueError('Seed torsion target scope must contain exact, separate modeled and flank residues.')
    if any(by_key[key].name not in STANDARD_PROTEINS or backbone(by_key[key]) is None for key in selected):
        raise ValueError('Seed torsion targets support only complete canonical protein backbones.')
    names = {}
    for residue in residues:
        named = {}
        for atom in residue.atoms():
            if atom.name in named:
                raise ValueError('Seed torsion targets require unambiguous atom identities.')
            named[atom.name] = atom
        names[residue] = named
    previous, following, bonds = {}, {}, set()
    for a, b in topology.bonds():
        bonds.add(frozenset((a.index, b.index)))
        if a.residue is b.residue:
            continue
        peptide = (a.residue.chain is b.residue.chain and {a.name, b.name} == {'C', 'N'}
                   and all(atom.element is not None and atom.element.symbol == atom.name for atom in (a, b))
                   and all(residue.name in STANDARD_PROTEINS and backbone(residue) is not None
                           for residue in (a.residue, b.residue)))
        if not peptide:
            if {residue_key(a.residue), residue_key(b.residue)} & selected:
                raise ValueError('Seed torsion targets cannot include crosslinked or nonprotein connections.')
            continue
        c, n = (a, b) if a.name == 'C' else (b, a)
        if c.residue in following or n.residue in previous:
            raise ValueError('Seed torsion targets require unambiguous peptide connectivity.')
        following[c.residue], previous[n.residue] = n.residue, c.residue

    expected = {}
    for key in selected:
        residue = by_key[key]
        if residue not in previous or residue not in following:
            raise ValueError('Seed torsion targets require complete peptide neighborhoods.')
        named = names[residue]
        expected[key, 'phi'] = [names[previous[residue]]['C'].index,
                                *[named[name].index for name in ('N', 'CA', 'C')]]
        expected[key, 'psi'] = [*[named[name].index for name in ('N', 'CA', 'C')],
                                names[following[residue]]['N'].index]
    checked, seen = [], set()
    for row in targets:
        if not isinstance(row, dict) or set(row) != {'residue', 'kind', 'atoms', 'target_radians'}:
            raise ValueError('Seed torsion target fields must be residue, kind, atoms and target_radians.')
        key, kind, indices, target = row['residue'], row['kind'], row['atoms'], row['target_radians']
        if (not isinstance(key, (list, tuple)) or len(key) != 4
                or any(not isinstance(value, str) for value in key)):
            raise ValueError('Seed torsion targets require exact residue identity strings.')
        key = tuple(key)
        if not isinstance(kind, str) or (key, kind) not in expected:
            raise ValueError('Seed torsion target lies outside the modeled/flank phi/psi scope.')
        if (key, kind) in seen:
            raise ValueError('Duplicate seed torsion target.')
        if (not isinstance(indices, (list, tuple)) or len(indices) != 4
                or any(type(i) is not int or not 0 <= i < len(atoms) for i in indices)):
            raise ValueError('Seed torsion targets require four valid integer atom indices.')
        if list(indices) != expected[key, kind] or any(
                frozenset(pair) not in bonds for pair in zip(indices, indices[1:])):
            raise ValueError('Seed torsion atom indices must match the exact bonded peptide phi/psi definition.')
        if (type(target) not in (int, float) or not np.isfinite(target)
                or abs(target) > np.pi + 1e-12):
            raise ValueError('Seed torsion targets require finite canonical radians in [-pi, pi].')
        initial = _dihedral(xyz[list(indices)])
        if initial is None:
            raise ValueError('Seed torsion targets cannot preserve a degenerate peptide angle.')
        difference = (np.deg2rad(initial) - target + np.pi) % (2 * np.pi) - np.pi
        if abs(difference) > 1e-8:
            raise ValueError('Seed torsion target differs from the supplied initial peptide angle.')
        seen.add((key, kind))
        checked.append({'residue': list(key), 'kind': kind, 'atoms': list(indices),
                        'target_radians': float(target), 'initial_degrees': initial})
    if seen != set(expected):
        raise ValueError('Seed torsion targets must cover every modeled/flank phi and psi exactly once.')
    return checked


def _minimize_attempt(topology, positions, forcefield, modeled_keys, *, max_iterations=1000,
                      mobile_flank_keys=(), reference_positions=None, target_overrides=None,
                      preserve_context_peptides=False, seed_torsion_targets=None):
    """Return (positions, report); move loops, H and explicitly allowed flanks.

    Caller must have removed input hydrogens and regenerated them during prep.
    Optional seed_torsion_targets must have passed the caller's independent
    backbone reference; exact scope/topology/seed-angle checks are repeated here.
    """
    if not 1 <= max_iterations <= 1000:
        raise ValueError("Loop refinement is limited to 1–1000 minimization iterations.")
    started = time.monotonic()
    keys = {tuple(key) for key in modeled_keys}
    flank_keys={tuple(key) for key in mobile_flank_keys}
    atoms = list(topology.atoms())
    xyz = np.asarray(positions.value_in_unit(unit.nanometer), dtype=float)
    if xyz.shape != (len(atoms), 3) or not np.isfinite(xyz).all():
        raise ValueError("Loop refinement needs finite coordinates matching its topology.")
    reference=np.asarray(reference_positions.value_in_unit(unit.nanometer),dtype=float) if reference_positions is not None else xyz
    if preserve_context_peptides:
        if reference_positions is None or reference.shape != xyz.shape or not np.isfinite(reference).all():
            raise ValueError('Observed context peptide preservation requires finite original reference coordinates.')
        found_flanks={residue_key(residue) for residue in topology.residues() if residue_key(residue) in flank_keys}
        if not flank_keys or found_flanks!=flank_keys or keys&flank_keys:
            raise ValueError('Observed context peptide preservation requires exact, separate mobile flank identities.')
    selected = [residue for residue in topology.residues() if residue_key(residue) in keys]
    found = [residue_key(residue) for residue in selected]
    if not keys or set(found) != keys or len(found) != len(set(found)):
        raise ValueError("Loop refinement requires exact modeled residue identities.")
    if any(not list(residue.atoms()) for residue in selected):
        raise ValueError("Every modeled loop residue must contain atoms.")
    seed_targets = None if seed_torsion_targets is None else _validated_seed_torsion_targets(
        topology, xyz, keys, flank_keys, seed_torsion_targets)
    # Recreate without constraints: constraints cannot couple massless observed
    # heavy atoms to mobile hydrogens. The same parameterized ForceField is used.
    system = forcefield.createSystem(topology, nonbondedMethod=app.NoCutoff,
                                    constraints=None, rigidWater=False)
    if seed_targets is not None:
        seed_restraints = mm.CustomTorsionForce('k*(1-cos(theta-target))')
        seed_restraints.setName('QualifiedSeedPhiPsiPreservation')
        seed_restraints.addPerTorsionParameter('k')
        seed_restraints.addPerTorsionParameter('target')
        for row in seed_targets:
            seed_restraints.addTorsion(*row['atoms'], [200.0, row['target_radians']])
        system.addForce(seed_restraints)
    # Temporary construction restraints enforce peptide planarity while the
    # loop adapts to fixed crystal anchors. New non-Pro links start trans;
    # Pro links retain the candidate's cis/trans basin. These restraints are
    # confined to this temporary System, never exported to the MD system.
    from .loop_geometry import _dihedral
    from .modified_residues import PARENT_RESIDUES
    peptide_restraints = mm.CustomTorsionForce("k*(1-cos(theta-target))")
    peptide_restraints.addPerTorsionParameter("k")
    peptide_restraints.addPerTorsionParameter("target")
    peptide_targets = []
    for a, b in topology.bonds():
        if a.residue == b.residue or {a.name, b.name} != {"C", "N"}:
            continue
        left, right = (a, b) if a.name == "C" else (b, a)
        pair_keys={residue_key(left.residue), residue_key(right.residue)}
        touches_modeled=bool(pair_keys & keys)
        if not touches_modeled and not (preserve_context_peptides and pair_keys & flank_keys):
            continue
        before = {atom.name: atom for atom in left.residue.atoms()}
        after = {atom.name: atom for atom in right.residue.atoms()}
        if "CA" not in before or "CA" not in after:
            raise ValueError("Modeled peptide restraints need both alpha-carbon anchors.")
        indices = [before["CA"].index, left.index, right.index, after["CA"].index]
        omega = _dihedral(xyz[indices])
        if omega is None:
            raise ValueError("Cannot initialize a peptide restraint from a degenerate loop.")
        is_proline = PARENT_RESIDUES.get(right.residue.name, right.residue.name) == "PRO"
        target = 0.0 if is_proline and abs(omega) < 90 else np.pi
        reference_omega=None
        if not touches_modeled:
            # Both residues were observed. Preserve the source cis/trans basin
            # even if the proposed context fragment changed it. This includes
            # observed non-Pro cis peptides; no new isomer state is guessed.
            reference_omega=_dihedral(reference[indices])
            if reference_omega is None or abs(abs(reference_omega)-90)<1e-6:
                raise ValueError('Observed context peptide has an undefined or ambiguous source basin.')
            target=0.0 if abs(reference_omega)<90 else np.pi
        if target_overrides is not None and touches_modeled:
            target=target_overrides.get(tuple(indices),target)
        peptide_restraints.addTorsion(*indices, [200.0, target])
        peptide_targets.append({"atoms": indices, "target_degrees": float(np.degrees(target)),
                                "initial_omega_degrees": omega,
                                "modeled_bond_before_residue": list(residue_key(right.residue)),
                                "touches_modeled_residue": touches_modeled,
                                "target_origin": 'modeled-link policy' if touches_modeled else 'original observed peptide basin',
                                "reference_omega_degrees": reference_omega})
    if peptide_targets:
        system.addForce(peptide_restraints)
    fixed, mobile = [], []
    for atom in atoms:
        if residue_key(atom.residue) in keys|flank_keys or atom.element == app.element.hydrogen:
            mobile.append(atom.index)
        else:
            system.setParticleMass(atom.index, 0 * unit.dalton)
            fixed.append(atom.index)
    if not mobile:
        raise ValueError("Loop refinement requires mobile modeled atoms.")
    flank_atoms=[atom.index for atom in atoms if residue_key(atom.residue) in flank_keys and atom.element!=app.element.hydrogen]
    chiral_targets=[]
    if flank_atoms:
        restraint=mm.CustomExternalForce(
            '0.5*k*r2+0.5*k_bound*max(0,sqrt(r2+1e-12)-r_bound)^2; '
            'r2=(x-x0)^2+(y-y0)^2+(z-z0)^2')
        restraint.addGlobalParameter('k',10000.0)
        # Guide minimization away from the endpoint displacement limit. A
        # harmonic anchor alone can be overpowered by strained covalent terms.
        # This soft boundary is a construction force, not coordinate clipping;
        # the independent 1 Å acceptance cap remains authoritative.
        restraint.addGlobalParameter('k_bound',1e6)
        restraint.addGlobalParameter('r_bound',0.075)
        for name in ('x0','y0','z0'):
            restraint.addPerParticleParameter(name)
        for index in flank_atoms:
            restraint.addParticle(index,reference[index].tolist())
        system.addForce(restraint)
        # The initial standard-residue stereochemistry was validated by prep.
        # Maintain those signs while flanks move, rather than allowing a low
        # energy but inverted residue to cross a planar intermediate.
        chirality=mm.CustomCompoundBondForce(4,
            '0.5*k*min(0,s*v-vmin)^2; '
            'v=(x2-x1)*((y3-y1)*(z4-z1)-(z3-z1)*(y4-y1))'
            '-(y2-y1)*((x3-x1)*(z4-z1)-(z3-z1)*(x4-x1))'
            '+(z2-z1)*((x3-x1)*(y4-y1)-(y3-y1)*(x4-x1))')
        for name in ('k','s','vmin'):
            chirality.addPerBondParameter(name)
        from .residue_identity import STANDARD_PROTEINS
        for residue in topology.residues():
            if residue_key(residue) not in keys|flank_keys or residue.name not in STANDARD_PROTEINS or residue.name=='GLY':
                continue
            named={atom.name:atom.index for atom in residue.atoms()}
            definitions=[('CA',('N','C','CB'))]
            if residue.name=='ILE':definitions.append(('CB',('CA','CG1','CG2')))
            if residue.name=='THR':definitions.append(('CB',('CA','OG1','CG2')))
            for center,neighbors in definitions:
                if not all(name in named for name in (center,*neighbors)):
                    continue
                indices=[named[name] for name in (center,*neighbors)]
                volume=float(np.linalg.det([reference[index]-reference[indices[0]] for index in indices[1:]]))
                if abs(volume)<1e-4:
                    raise ValueError('Cannot restrain a near-planar source stereocenter during local flank relaxation.')
                minimum=max(2e-4,.25*abs(volume))
                chirality.addBond(indices,[1e9,float(np.sign(volume)),minimum])
                chiral_targets.append({'residue':list(residue_key(residue)),'center':center,
                    'reference_signed_volume_nm3':volume,'minimum_oriented_volume_nm3':minimum})
        system.addForce(chirality)
    platform = mm.Platform.getPlatformByName("CPU")
    integrator = mm.VerletIntegrator(.001 * unit.picosecond)
    properties={"Threads":str(config.CPU_THREADS),"DeterministicForces":"true"}
    context = mm.Context(system, integrator, platform, properties)
    try:
        context.setPositions(positions)
        before = context.getState(getEnergy=True).getPotentialEnergy().value_in_unit(unit.kilojoule_per_mole)
        if not np.isfinite(before):
            raise ValueError("Loop refinement starts with non-finite energy.")
        from .loop_geometry import loop_geometry_report
        # Relieve collisions before tightening peptide planes. Applying the
        # strong construction force to a colliding initial fragment can drive
        # other bonds or chirality into a bad basin. Keep the total budget.
        first_budget=max(1,max_iterations//2)
        budgets=[first_budget,max_iterations-first_budget]
        phases=[];strength=200.0
        for phase,budget in enumerate(budgets):
            if budget==0:
                continue
            if phase and peptide_targets and not phases[-1]['geometry_accepted']:
                strength=1000.0
                for index in range(peptide_restraints.getNumTorsions()):
                    a,b,c,d,parameters=peptide_restraints.getTorsionParameters(index)
                    peptide_restraints.setTorsionParameters(index,a,b,c,d,[strength,parameters[1]])
                peptide_restraints.updateParametersInContext(context)
            phase_before=context.getState(getEnergy=True).getPotentialEnergy().value_in_unit(unit.kilojoule_per_mole)
            try:
                mm.LocalEnergyMinimizer.minimize(context,10 * unit.kilojoule_per_mole / unit.nanometer,budget)
            except mm.OpenMMException as exc:
                # The native minimizer reports a NaN coordinate itself; keep the
                # rejection explicit instead of surfacing the raw engine text.
                if 'nan' not in str(exc).lower():
                    raise
                raise ValueError('Loop refinement produced non-finite coordinates for the reconstructed loop; '
                                 'the candidate is rejected and no prepared model was accepted. '
                                 'Supply a repaired structure or try a different seed.') from exc
            phase_state=context.getState(getPositions=True,getEnergy=True)
            phase_xyz=np.asarray(phase_state.getPositions(asNumpy=True).value_in_unit(unit.nanometer))
            phase_after=phase_state.getPotentialEnergy().value_in_unit(unit.kilojoule_per_mole)
            if not np.isfinite(phase_after) or not np.isfinite(phase_xyz).all():
                raise ValueError('Loop refinement produced non-finite coordinates or energy.')
            if not np.array_equal(phase_xyz[fixed],xyz[fixed]):
                raise ValueError('Loop refinement moved an observed heavy atom; the candidate is rejected.')
            if phase_after>phase_before+max(1e-3,abs(phase_before)*1e-8):
                raise ValueError('A loop refinement stage increased its potential energy; the candidate is rejected.')
            geometry=loop_geometry_report(topology,phase_state.getPositions(),keys|flank_keys)
            phases.append({'phase':'collision relief' if phase==0 else 'peptide geometry refinement',
                           'iteration_budget':budget,'peptide_strength_kj_mol':strength,
                           'energy_before_kj_mol':float(phase_before),'energy_after_kj_mol':float(phase_after),
                           'geometry_accepted':geometry['accepted'],'geometry_error_count':len(geometry['errors'])})
        state = context.getState(getPositions=True, getEnergy=True, getForces=True)
        after = state.getPotentialEnergy().value_in_unit(unit.kilojoule_per_mole)
        refined = np.asarray(state.getPositions(asNumpy=True).value_in_unit(unit.nanometer))
        forces = np.asarray(state.getForces(asNumpy=True).value_in_unit(unit.kilojoule_per_mole / unit.nanometer))
        unchanged = bool(np.array_equal(refined[fixed], xyz[fixed]))
        if not (np.isfinite(after) and np.isfinite(refined).all() and np.isfinite(forces).all()):
            raise ValueError("Loop refinement produced non-finite coordinates, forces or energy.")
        if not unchanged:
            raise ValueError("Loop refinement moved an observed heavy atom; the candidate is rejected.")
        # Report comparable endpoint energies under the FINAL temporary
        # Hamiltonian; changing construction strength changes the objective.
        if strength!=200.0:
            context.setPositions(positions)
            before=context.getState(getEnergy=True).getPotentialEnergy().value_in_unit(unit.kilojoule_per_mole)
        if after > before + max(1e-3, abs(before) * 1e-8):
            raise ValueError("Bounded loop refinement increased the potential energy; the candidate is rejected.")
        report = {"method": "OpenMM LocalEnergyMinimizer with the prepared force field and temporary peptide construction restraints, NoCutoff, no constraints; all non-loop heavy atoms fixed by zero mass. Loop atoms and regenerated hydrogens mobile.",
                  "max_iterations": max_iterations, "tolerance_kj_mol_nm": 10,
                  "platform_properties":properties,
                  "fixed_heavy_atoms": len(fixed), "mobile_atoms": len(mobile),
                  "fixed_heavy_coordinates_preserved_exactly": unchanged,
                  "energy_before_kj_mol": float(before), "energy_after_kj_mol": float(after),
                  "energy_comparison": "Both endpoint energies use the final construction strength; per-phase energies use that phase's strength.",
                  "phases": phases,
                  "energy_includes_construction_restraints": bool(peptide_targets or seed_targets),
                  "mobile_force_rms_kj_mol_nm": float(np.sqrt(np.mean(forces[mobile] ** 2))),
                  "elapsed_seconds": time.monotonic() - started,
                  "peptide_construction_restraints": {"initial_strength_kj_mol": 200, "strength_kj_mol": strength, "targets": peptide_targets,
                      "observed_context_peptides_preserved": bool(preserve_context_peptides),
                      "exported_to_simulation": False,
                      "assumption": "New non-Pro peptide bonds initialized trans; Pro bonds retain the candidate cis/trans basin. Native isomer states of missing residues are unknown."},
                  "scope": "Local starting-geometry relaxation before periodic solvent equilibration; no convergence or native-loop claim. Stereo and geometric acceptance are checked separately."}
        if seed_targets is not None:
            report['initial_torsion_preservation'] = {
                'strength_kj_mol': 200.0, 'targets': seed_targets,
                'seed_target_angles_verified': True, 'native_reference_check_performed_here': False,
                'caller_reference_qualification_required': True, 'exported_to_simulation': False,
                'scope': 'Temporary circular phi/psi prior around caller-qualified seed angles. '
                         'Construction energies/forces include this prior; no independent validation, '
                         'native-conformation or unbiased-equilibrium claim.'}
        if flank_atoms:
            report['method']='OpenMM local minimization with loop atoms, regenerated hydrogens and explicitly listed restrained flanks mobile. Other heavy atoms fixed by zero mass; temporary construction forces are not exported.'
            moved=[{'atom':[ *residue_key(atoms[index].residue),atoms[index].name],
                    'displacement_nm':float(np.linalg.norm(refined[index]-reference[index]))} for index in flank_atoms]
            report['flank_relaxation']={'residues':[list(key) for key in sorted(flank_keys)],
                'atoms':moved,'maximum_displacement_nm':max(row['displacement_nm'] for row in moved),
                'maximum_allowed_displacement_nm':0.1,'positional_restraint_k_kj_mol_nm2':10000,
                'displacement_boundary_start_nm':0.075,'displacement_boundary_k_kj_mol_nm2':1e6,
                'displacement_reference':'Coordinates entering loop refinement, after recorded sidechain sampling; backbone coordinates still equal the input backbone.',
                'chirality_construction':{'strength_kj_mol_nm6':1e9,'targets':chiral_targets,
                    'scope':'Flat-bottom signed-volume barriers from previously validated standard-residue centers; independent final stereochemistry checks remain required.', 'exported_to_simulation':False},
                'exported_to_simulation':False,'scope':'At most one adjacent standard residue at each end of a failed loop; original atom identities retained.'}
        return refined * unit.nanometer, report
    finally:
        del context, integrator


def _failed_loop_flanks(topology, modeled_keys, geometry, protected_residues):
    """Immediate neighbors of failed stretches, excluding chemical crosslinks."""
    from .residue_identity import STANDARD_PROTEINS
    keys={tuple(key) for key in modeled_keys};protected={tuple(key) for key in protected_residues}
    atoms=list(topology.atoms());residues={residue_key(r):r for r in topology.residues()}
    neighbors={key:set() for key in residues}
    for a,b in topology.bonds():
        ka,kb=residue_key(a.residue),residue_key(b.residue)
        if ka==kb:
            continue
        if {a.name,b.name}=={'C','N'} and a.residue.chain==b.residue.chain:
            neighbors[ka].add(kb);neighbors[kb].add(ka)
        else:
            protected.update((ka,kb))
    failed_atoms=set()
    for field in ('bond_checks','angle_checks','omega_checks'):
        for check in geometry.get(field,[]):
            if not check.get('accepted',False):
                failed_atoms.update(check.get('atoms',[]))
    for pair in geometry.get('gross_collisions',[]):
        failed_atoms.update(pair.get('atoms',[]) if isinstance(pair,dict) else pair)
    failed_keys={residue_key(atoms[index].residue) for index in failed_atoms if 0<=index<len(atoms)}
    remaining=set(keys);flanks=set()
    while remaining:
        first=remaining.pop();segment={first};pending=[first]
        while pending:
            for key in neighbors[pending.pop()]&remaining:
                remaining.remove(key);segment.add(key);pending.append(key)
        adjacent=set().union(*(neighbors[key] for key in segment))-keys
        if not failed_keys&segment:
            continue
        flanks.update(key for key in adjacent if key not in protected and residues[key].name in STANDARD_PROTEINS)
    return flanks


def minimize_modeled_loops(topology, positions, forcefield, modeled_keys, *, max_iterations=1000, protected_residues=()):
    """Try fixed observed atoms first, then bounded restrained local flanks.

    Each attempt is capped at 1,000 iterations. The fallback may move only
    immediate, uncrosslinked standard flanks, at most 1 Å; it is provisional
    local remodeling, recorded separately from the original fixed-atom path.
    """
    from .loop_geometry import loop_geometry_report
    refined,report=_minimize_attempt(topology,positions,forcefield,modeled_keys,max_iterations=max_iterations)
    geometry=loop_geometry_report(topology,refined,modeled_keys)
    report['total_iteration_budget']=max_iterations
    if geometry['accepted']:
        return refined,report
    flanks=_failed_loop_flanks(topology,modeled_keys,geometry,protected_residues)
    if not flanks:
        return refined,report
    first_report=report
    targets={tuple(row['atoms']):np.deg2rad(row['target_degrees']) for row in report['peptide_construction_restraints']['targets']}
    refined,report=_minimize_attempt(topology,refined,forcefield,modeled_keys,max_iterations=max_iterations,
                                    mobile_flank_keys=flanks,reference_positions=positions,target_overrides=targets)
    initial_omegas={tuple(row['atoms']):row['initial_omega_degrees'] for row in first_report['peptide_construction_restraints']['targets']}
    for row in report['peptide_construction_restraints']['targets']:
        row['fallback_start_omega_degrees']=row['initial_omega_degrees']
        row['initial_omega_degrees']=initial_omegas[tuple(row['atoms'])]
    report['fixed_anchor_attempt']=first_report
    report['fixed_anchor_geometry']=geometry
    report['total_iteration_budget']=2*max_iterations
    report['elapsed_seconds']+=first_report['elapsed_seconds']
    report['flank_relaxation']['geometry']=loop_geometry_report(topology,refined,{tuple(key) for key in modeled_keys}|flanks)
    report['flank_relaxation']['within_displacement_limit']=report['flank_relaxation']['maximum_displacement_nm']<=0.1
    return refined,report
