"""Independent CCTBX/coordinate screen of one completed native KIC comparison."""
import argparse
import hashlib
import json
from pathlib import Path
import sys

import numpy as np
from openmm import app, unit

ROOT = Path(__file__).resolve().parents[4]
sys.path.insert(0, str(ROOT/'docs/audit/advanced-chemistry/covalent-v1-focused'))
from rama_reference import rama_report, residue_key, dihedral64


def sha(p):
    return hashlib.sha256(p.read_bytes()).hexdigest()


def key(a):
    return (*residue_key(a.residue), a.name)


def angle(xyz):
    left, right = xyz[0]-xyz[1], xyz[2]-xyz[1]
    return float(np.degrees(np.arccos(np.clip(np.dot(left, right)/(np.linalg.norm(left)*np.linalg.norm(right)), -1, 1))))


def geometry(rows):
    lookup = {(int(r['identity'][1]), r['identity'][4]): np.asarray(r['xyz_angstrom']) for r in rows}
    numbers = sorted({n for n, a in lookup})
    bonds, angles, omegas = [], [], []
    def bond(members, bounds):
        value = float(np.linalg.norm(lookup[members[0]]-lookup[members[1]]))
        bonds.append({'atoms': members, 'angstrom': value, 'accepted': bounds[0] <= value <= bounds[1]})
    def triple(members):
        value = angle(np.asarray([lookup[k] for k in members]))
        angles.append({'atoms': members, 'degrees': value, 'accepted': 80 <= value <= 150})
    for n in numbers:
        for a, b, bounds in [('N', 'CA', (1.25, 1.70)), ('CA', 'C', (1.30, 1.75)), ('C', 'O', (1.05, 1.45))]:
            bond([(n, a), (n, b)], bounds)
        triple([(n, 'N'), (n, 'CA'), (n, 'C')])
        triple([(n, 'CA'), (n, 'C'), (n, 'O')])
    for n, m in zip(numbers, numbers[1:]):
        bond([(n, 'C'), (m, 'N')], (1.15, 1.60))
        for members in [[(n, 'CA'), (n, 'C'), (m, 'N')], [(n, 'O'), (n, 'C'), (m, 'N')], [(n, 'C'), (m, 'N'), (m, 'CA')]]:
            triple(members)
        value = dihedral64(np.asarray([lookup[k] for k in [(n, 'CA'), (n, 'C'), (m, 'N'), (m, 'CA')]]))
        departure = min(abs(value), 180-abs(value))
        omegas.append({'first': n, 'second': m, 'degrees': value, 'departure_degrees': departure, 'accepted': departure <= 35})
    return {'bond_checks': bonds, 'angle_checks': angles, 'omega_checks': omegas,
            'accepted': all(r['accepted'] for r in bonds+angles+omegas),
            'scope': 'Backbone-only gross checks; no sidechain chirality, contact, rotamer, environment or force-field validation.'}


def save_backbone(path, topology, xyz, first, last):
    new = app.Topology()
    chain = new.addChain('A')
    coords = []
    names = []
    for r in topology.residues():
        if r.chain.id != 'A' or not first <= int(r.id) <= last:
            continue
        added = new.addResidue(r.name, chain, r.id, r.insertionCode)
        lookup = {}
        for a in r.atoms():
            if a.name in {'N', 'CA', 'C', 'O'}:
                lookup[a.name] = new.addAtom(a.name, a.element, added)
                coords.append(xyz[a.index])
        for a, b in [('N', 'CA'), ('CA', 'C'), ('C', 'O')]:
            new.addBond(lookup[a], lookup[b])
        if names:
            new.addBond(names[-1]['C'], lookup['N'])
        names.append(lookup)
    with path.open('w') as stream:
        app.PDBFile.writeFile(new, np.asarray(coords)*unit.nanometer, stream, keepIds=True)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--run', type=Path, required=True)
    parser.add_argument('--control', type=Path, required=True)
    parser.add_argument('--seed', type=Path)
    parser.add_argument('--source', type=Path)
    args = parser.parse_args()
    output = args.run/'reference-screen'
    output.mkdir(exist_ok=False)
    native_path = args.run/'native/result.json'
    seed_path = args.seed or args.run/'initial-refinement-seed.pdb'
    files = [native_path, seed_path, args.control, Path(__file__)]
    if args.source is not None:
        files.append(args.source)
    hashes = {str(p): sha(p) for p in files}
    (output/'implementation.py').write_bytes(Path(__file__).read_bytes())
    native = json.loads(native_path.read_text())
    reports = []
    for case in native['cases']:
        control = case['name'] == 'intact-control'
        seed = app.PDBFile(str(args.control if control else seed_path))
        atoms = list(seed.topology.atoms())
        lookup = {key(a): a.index for a in atoms}
        if len(lookup) != len(atoms):
            raise ValueError('Ambiguous exact atom identity')
        entry = np.asarray(seed.positions.value_in_unit(unit.nanometer))
        # Coherent research seeds can move observed context. The comparison
        # baseline and atoms outside the candidate span must remain the exact
        # original observed backbone, not those temporary native coordinates.
        if args.source is not None and not control:
            source = app.PDBFile(str(args.source))
            source_xyz = np.asarray(source.positions.value_in_unit(unit.nanometer))
            for atom in source.topology.atoms():
                if atom.name in {'N', 'CA', 'C', 'O'} and key(atom) in lookup:
                    entry[lookup[key(atom)]] = source_xyz[atom.index]
        entry_report = rama_report(seed.topology, entry*unit.nanometer)
        case_dir = output/case['name']
        case_dir.mkdir()
        candidates = []
        for ti, trial in enumerate(case['trials']):
            input_geometry = geometry(trial.get('input_backbone_atoms', case['initial_backbone_atoms']))
            for si, sol in enumerate(trial['solutions']):
                xyz = entry.copy()
                for row in sol['backbone_atoms']:
                    xyz[lookup[tuple(row['identity'])]] = np.asarray(row['xyz_angstrom'])/10
                modeled = {residue_key(r) for r in seed.topology.residues() if r.chain.id == 'A' and int(r.id) in case['modeled']}
                # 1e-4 Å distinguishes meaningful float-coordinate change from
                # native float32 rotation/closure roundoff, not quality acceptance.
                changed = set(np.where(np.linalg.norm(xyz-entry, axis=1)*10 > 1e-4)[0])
                affected = set(modeled)
                for row in entry_report['rows']:
                    indices = row['phi_atoms']+row['psi_atoms']
                    if any(residue_key(atoms[i].residue) in modeled or i in changed for i in indices):
                        affected.add(tuple(row['residue']))
                if control:
                    affected = {tuple(r['residue']) for r in entry_report['rows']}
                report = rama_report(seed.topology, xyz*unit.nanometer, affected)
                gross = geometry(sol['backbone_atoms'])
                bond_delta = max(abs(a['angstrom']-b['angstrom']) for a, b in zip(gross['bond_checks'], input_geometry['bond_checks']))
                angle_delta = max(abs(a['degrees']-b['degrees']) for a, b in zip(gross['angle_checks'], input_geometry['angle_checks']))
                candidate = {'trial': ti, 'solution': si, 'pivots': trial['pivots'], 'draw_index': trial['draw_index'],
                             'selected_residues': sorted(affected),
                             'maximum_observed_N_CA_C_O_displacement_angstrom': sol['maximum_observed_N_CA_C_O_displacement_angstrom'],
                             'observed_cap_pass': sol['within_original_1A_observed_backbone_cap'],
                             'observed_omega_basins_preserved': sol['observed_omega_basins_preserved'],
                             'rama': report, 'gross_backbone_geometry': gross,
                             'maximum_bond_length_change_from_native_trial_input_angstrom': bond_delta,
                             'maximum_backbone_angle_change_from_native_trial_input_degrees': angle_delta}
                candidate['plausible_backbone_screen_pass'] = bool(sol['within_original_1A_observed_backbone_cap']
                    and sol['observed_omega_basins_preserved'] and report['all_selected_scored_without_outliers'] and gross['accepted'])
                if candidate['plausible_backbone_screen_pass']:
                    path = case_dir/f'trial-{ti}-solution-{si}-BACKBONE-ONLY.pdb'
                    save_backbone(path, seed.topology, xyz, case['start'], case['end'])
                    candidate['backbone_only_pdb'] = str(path)
                    candidate['backbone_only_pdb_sha256'] = sha(path)
                candidates.append(candidate)
        result = {'name': case['name'], 'candidate_count': len(candidates),
                  'within_observed_cap_and_omega': sum(c['observed_cap_pass'] and c['observed_omega_basins_preserved'] for c in candidates),
                  'cap_omega_and_rama_pass': sum(c['observed_cap_pass'] and c['observed_omega_basins_preserved'] and c['rama']['all_selected_scored_without_outliers'] for c in candidates),
                  'plausible_backbone_screen_pass_count': sum(c['plausible_backbone_screen_pass'] for c in candidates),
                  'candidates': candidates}
        (case_dir/'result.json').write_text(json.dumps(result, indent=2)+'\n')
        print(json.dumps({k:v for k,v in result.items() if k!='candidates'}),flush=True)
        reports.append(result)
    if hashes != {str(p):sha(p) for p in files}:
        raise ValueError('Frozen inputs changed during independent screening')
    result = {'source_sha256':hashes,'cases': reports,
              'native_coordinate_movement_roundoff_tolerance_angstrom':1e-4,
              'app_ready':False,'physical_model_validated':False,
              'sidechains_or_full_atom_environment_screened':False,
              'scope':'Backbone feasibility comparison only; no replacement for sidechain rebuilding, chirality/contact checks, complete force-field refinement or final admission.'}
    (output/'result.json').write_text(json.dumps(result,indent=2)+'\n')


if __name__=='__main__':
    main()
