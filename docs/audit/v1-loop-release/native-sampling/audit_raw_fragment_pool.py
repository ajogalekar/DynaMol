"""Independent raw-fragment screen, evaluated after the closure plan is frozen."""
import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
from openmm import app, unit
from screen_kic_backbones import geometry, rama_report, residue_key, dihedral64, key


def sha(p):
    return hashlib.sha256(p.read_bytes()).hexdigest()


def main():
    parser = argparse.ArgumentParser()
    for name in ['run','source','base-seed']:
        parser.add_argument('--'+name,type=Path,required=True)
    args = parser.parse_args()
    output = args.run/'raw-reference-screen'
    output.mkdir(exist_ok=False)
    pool_path = args.run/'native/raw-pool.json'
    frozen_plan = args.run/'native/frozen-closure-plan.json'
    inputs = [pool_path,frozen_plan,args.source,args.base_seed,Path(__file__)]
    hashes = {str(p):sha(p) for p in inputs}
    pool = json.loads(pool_path.read_text())
    base = app.PDBFile(str(args.base_seed))
    atoms = list(base.topology.atoms())
    lookup = {key(a):a.index for a in atoms}
    if len(lookup)!=len(atoms):
        raise ValueError('Duplicate exact backbone identity')
    original = app.PDBFile(str(args.source))
    source_xyz = np.asarray(original.positions.value_in_unit(unit.nanometer))
    source = {key(a):source_xyz[a.index] for a in original.topology.atoms() if a.name in {'N','CA','C','O'}}
    baseline = np.asarray(base.positions.value_in_unit(unit.nanometer))
    for identity, xyz in source.items():
        if identity in lookup:
            baseline[lookup[identity]] = xyz
    original_report = rama_report(base.topology, baseline*unit.nanometer)
    modeled = {('A','177','','ALA'),('A','178','','LYS'),('A','179','','ASN')}
    reports = []
    for candidate in pool['candidates']:
        xyz = baseline.copy()
        by_number = {}
        displacements = []
        for row in candidate['backbone_atoms']:
            identity = tuple(row['identity'])
            xyz[lookup[identity]] = np.asarray(row['xyz_angstrom'])/10
            by_number.setdefault(int(identity[1]),{})[identity[4]] = identity
            if identity[:4] not in modeled:
                displacements.append({'identity':identity,'angstrom':float(np.linalg.norm(xyz[lookup[identity]]-source[identity])*10)})
        maximum = max(d['angstrom'] for d in displacements)
        if abs(maximum-candidate['raw_maximum_observed_backbone_displacement_angstrom']) > 1e-4:
            raise ValueError('Native/raw independent complete-anchor displacement differs beyond floating-point allowance')
        changed = set(np.where(np.linalg.norm(xyz-baseline,axis=1)*10>1e-4)[0])
        selected = set(modeled)
        for row in original_report['rows']:
            indices = row['phi_atoms']+row['psi_atoms']
            if any(i in changed or residue_key(atoms[i].residue) in modeled for i in indices):
                selected.add(tuple(row['residue']))
        rama = rama_report(base.topology,xyz*unit.nanometer,selected)
        omega = []
        numbers = sorted(by_number)
        for first,last in zip(numbers,numbers[1:]):
            identifiers = [by_number[first]['CA'],by_number[first]['C'],by_number[last]['N'],by_number[last]['CA']]
            if any(ident[:4] in modeled for ident in identifiers):
                continue
            original_value = dihedral64(np.asarray([source[i] for i in identifiers]))
            value = dihedral64(xyz[[lookup[i] for i in identifiers]])
            original_basin = 'cis' if abs(original_value)<90 else 'trans'
            candidate_basin = 'cis' if abs(value)<90 else 'trans'
            omega.append({'first':first,'last':last,'source_degrees':original_value,
                          'candidate_degrees':value,'source_basin':original_basin,
                          'candidate_basin':candidate_basin,'preserved':original_basin==candidate_basin})
        gross = geometry(candidate['backbone_atoms'])
        report = {'candidate_id':candidate['id'],'source_span':[candidate['source_start'],candidate['source_end']],
                  'maximum_observed_N_CA_C_O_displacement_angstrom':maximum,'within_original_cap':maximum<=1,
                  'observed_omega_checks':omega,'observed_omega_basins_preserved':all(r['preserved'] for r in omega),
                  'gross_backbone_geometry':gross,'rama':rama,'selected_residues':sorted(selected)}
        report['backbone_screen_pass'] = maximum<=1 and report['observed_omega_basins_preserved'] and gross['accepted'] and rama['all_selected_scored_without_outliers']
        reports.append(report)
    if hashes != {str(p):sha(p) for p in inputs}:
        raise ValueError('Frozen inputs or plan changed during raw reference screen')
    (output/'implementation.py').write_bytes(Path(__file__).read_bytes())
    result = {'source_sha256':hashes,'candidates':reports,'candidate_count':len(reports),
              'within_original_cap_count':sum(r['within_original_cap'] for r in reports),
              'gross_backbone_geometry_pass_count':sum(r['gross_backbone_geometry']['accepted'] for r in reports),
              'rama_pass_count':sum(r['rama']['all_selected_scored_without_outliers'] for r in reports),
              'original_observed_omega_basins_preserved_count':sum(r['observed_omega_basins_preserved'] for r in reports),
              'combined_backbone_pass_count':sum(r['backbone_screen_pass'] for r in reports),
              'scope':'Raw native fragment diagnostics; no KIC or source-restored join has been introduced in these coordinates. No full-atom model acceptance.'}
    (output/'result.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({k:v for k,v in result.items() if k not in {'candidates','source_sha256'}},indent=2))


if __name__=='__main__':
    main()
