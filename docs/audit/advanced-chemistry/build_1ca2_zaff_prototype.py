"""Reproducible full 1CA2 native ZAFF preparation specimen, outside the app.

Explicit Zn-HHH-water (center 6), not an inferred protonation state. Retain all
source protein, Zn and crystal water heavy atoms. No dynamics or model-accuracy
claim is made by construction or finite mechanical checks.
"""
from __future__ import annotations
import argparse
import hashlib
import json
import os
from pathlib import Path
import random
import shutil
import subprocess
import sys

import numpy as np
from openmm import app, unit, Platform
import parmed
from pdbfixer import PDBFixer

BASE = Path(__file__).resolve().parent
ROOT = BASE.parents[2]
sys.path.insert(0, str(ROOT))
from backend.bonded_site_validation import validate_bonded_site, check_finite_forces
from export_native_parameter_manifest import export


def identity(atom):
    r = atom.residue
    return f"{r.chain.id}:{r.id}:{r.insertionCode.strip()}:{r.name}:{atom.name}"


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--output', type=Path, default=BASE/'prototype-1ca2-zaff-water')
    args = parser.parse_args()
    output = args.output.resolve(); output.mkdir(exist_ok=False)
    source = BASE/'inputs/1CA2.cif'
    fixer = PDBFixer(filename=str(source))
    before = np.asarray(fixer.positions.value_in_unit(unit.angstrom))
    original = {identity(a): before[a.index].copy() for a in fixer.topology.atoms()}
    if len(original) != len(before): raise ValueError('Nonunique original identities')
    fixer.findMissingResidues(); fixer.findMissingAtoms()
    missing_terminal = [(r.id, r.name, names) for r,names in fixer.missingTerminals.items()]
    if fixer.missingResidues or fixer.missingAtoms or missing_terminal != [('260','PHE',['OXT'])]:
        raise ValueError('Frozen 1CA2 expected repair inventory changed')
    fixer.addMissingAtoms(seed=2026)
    modeller = app.Modeller(fixer.topology, fixer.positions)
    variants = []
    for r in modeller.topology.residues():
        variants.append({'94':'HID','96':'HID','119':'HIE'}.get(r.id) if r.chain.id=='A' else None)
    random.seed(2026)
    actual_variants = modeller.addHydrogens(forcefield=None, pH=7.0, variants=variants,
        platform=Platform.getPlatformByName('Reference'))
    xyz = np.asarray(modeller.positions.value_in_unit(unit.angstrom))
    source_map = {}
    before_rename = {}
    for a in modeller.topology.atoms():
        ident = identity(a); before_rename[a.index] = ident
        if ident in original and not np.allclose(xyz[a.index], original[ident], atol=1e-12,rtol=0):
            raise ValueError('Observed source heavy atom moved during OXT/H placement')
    observed = {before_rename[a.index] for a in modeller.topology.atoms() if before_rename[a.index] in original}
    if observed != set(original): raise ValueError('Observed atom identity lost')
    res = {(r.chain.id,r.id):r for r in modeller.topology.residues()}
    atoms = lambda r:{a.name:a for a in r.atoms()}
    context_distances = {}
    for hid in ('94','96'):
        for one,two in [('O','N'),('N','O')]:
            a,b = atoms(res['A',hid])[one],atoms(res['A','119'])[two]
            context_distances[f'His{hid}:{one}-His119:{two}'] = float(np.linalg.norm(xyz[a.index]-xyz[b.index]))
    if not (context_distances['His94:O-His119:N'] < 3.2 and context_distances['His96:O-His119:N'] > 5):
        raise ValueError('Published distinct-HID context mapping not reproduced')
    assignments = {('A','94'):'HD4',('A','96'):'HD5',('A','119'):'HE2',('B','262'):'ZN6',('C','263'):'WT1'}
    histidines = []
    for r,variant in zip(modeller.topology.residues(),actual_variants,strict=True):
        if r.name == 'HIS':
            histidines.append({'source_chain':r.chain.id,'source_resid':r.id,'variant':variant,
                'selection':'explicit metal donor' if (r.chain.id,r.id) in assignments else 'OpenMM pH7 geometric default; not pKa calculation'})
            r.name = variant
        if (r.chain.id,r.id) in assignments: r.name = assignments[r.chain.id,r.id]
    # OpenMM's N-terminal H/H2/H3 and native Amber H1/H2/H3 are aliases.
    # Rename the explicitly identified added H, never add a fourth N hydrogen.
    first = next(modeller.topology.residues())
    first_atoms = atoms(first)
    if set(('H','H2','H3')) <= set(first_atoms) and 'H1' not in first_atoms:
        first_atoms['H'].name = 'H1'
    for r in modeller.topology.residues():
        for a in r.atoms():
            source_map[(r.index,a.name)] = {'source_or_added_id':before_rename[a.index],
                'role':'observed source heavy atom' if before_rename[a.index] in original else 'added hydrogen' if a.element==app.element.hydrogen else 'added terminal OXT',
                'prepared_atom_name':a.name,'prepared_residue_name':r.name,'prepared_xyz_angstrom':xyz[a.index].tolist()}
        # Native LEaP input and output use the exact sequential residue map.
        r.id = str(r.index+1)
    with (output/'prepared-input.pdb').open('w') as stream:
        app.PDBFile.writeFile(modeller.topology,modeller.positions,stream,keepIds=True)
    # A bound solvent template has PREP main-chain connection metadata. Explicit
    # TER boundaries prevent LEaP from attempting connections to the next water.
    lines = (output/'prepared-input.pdb').read_text().splitlines()
    separated = []
    for i,line in enumerate(lines):
        separated.append(line)
        if line.startswith(('ATOM  ','HETATM')) and line[17:20].strip() in {'WT1','ZN6','HOH'}:
            following = lines[i+1] if i+1<len(lines) else ''
            if following.startswith(('ATOM  ','HETATM')) and following[21:27]!=line[21:27]:
                separated.append('TER')
    (output/'prepared-input.pdb').write_text('\n'.join(separated)+'\n')
    for name in ('ZAFF.prep','ZAFF.frcmod'): shutil.copy2(BASE/'reference-zaff'/name,output/name)
    zn = res['B','262']; donors = [(res['A','94'],'NE2'),(res['A','96'],'NE2'),(res['A','119'],'ND1'),(res['C','263'],'O')]
    leap = ['source leaprc.protein.ff14SB','source leaprc.water.tip3p',
        'addAtomTypes { { "ZN" "Zn" "sp3" } { "N5" "N" "sp2" } { "N6" "N" "sp2" } { "N7" "N" "sp2" } { "O1" "O" "sp3" } }',
        'loadamberprep ZAFF.prep','loadamberparams ZAFF.frcmod','mol = loadpdb prepared-input.pdb']
    leap += [f'bond mol.{zn.index+1}.ZN mol.{r.index+1}.{name}' for r,name in donors]
    leap += ['check mol','savepdb mol prepared-native.pdb','saveamberparm mol prepared.prmtop prepared.inpcrd','quit']
    (output/'tleap.in').write_text('\n'.join(leap)+'\n')
    amber=ROOT/'.tools/ambertools';env=dict(os.environ,AMBERHOME=str(amber),OMP_NUM_THREADS='2',OPENBLAS_NUM_THREADS='2')
    with (output/'tleap.log').open('w') as stream:
        proc=subprocess.run([str(amber/'bin/tleap'),'-s','-f','tleap.in'],cwd=output,env=env,
            stdout=stream,stderr=subprocess.STDOUT,timeout=120)
    if proc.returncode or not (output/'prepared.prmtop').exists(): raise ValueError('Native LEaP failed; preserved log')
    native=parmed.load_file(str(output/'prepared.prmtop'),xyz=str(output/'prepared.inpcrd'))
    mapping=[]; matched=set(); max_rounding=0.
    for a in native.atoms:
        record=source_map.get((a.residue.idx,a.name))
        if record is None: raise ValueError(f'Unmapped native atom {a.residue.idx}:{a.name}')
        record=dict(record,native_index=a.idx)
        delta=float(np.max(np.abs(native.coordinates[a.idx]-record['prepared_xyz_angstrom'])))
        max_rounding=max(max_rounding,delta)
        if delta>.0005001: raise ValueError('Native coordinate change exceeds PDB serialization rounding')
        if record['role']=='observed source heavy atom': matched.add(record['source_or_added_id'])
        mapping.append(record)
    if matched!=set(original) or len(mapping)!=len(source_map): raise ValueError('Final native atom inventory differs')
    (output/'source-atom-map.json').write_text(json.dumps({'atom_map':mapping},indent=2)+'\n')
    manifest=export(output/'prepared.prmtop',output/'native-parameters.json',output/'source-atom-map.json')
    manifest['model_id']='1CA2-published-ZAFF-center6-ff14SB-TIP3P'
    state={'status':'explicit published water-state benchmark model; not an experimental protonation assignment',
        'oxidation_state':'Zn(II)','electronic_spin':0,'bound_solvent':'neutral water, two hydrogens; center6',
        'alternative':'Separate center7 hydroxide: HD6/HD7/HE3/OH1/ZN7, one H and whole-site charge +1; not substituted automatically',
        'chosen_state_rationale':'Retains deposited HOH263 identity using the published Zn-HHH-water parameter family; source CIF has no crystallization pH field. Protonation uncertainty remains.',
        'forcefield':'ff14SB + official ZAFF.prep/frcmod + TIP3P; unmodified original ZAFF Zn Lennard-Jones parameters',
        'histidines':histidines}
    manifest['chemical_state']=state
    (output/'native-parameters.json').write_text(json.dumps(manifest,indent=2)+'\n')
    loaded=app.AmberPrmtopFile(str(output/'prepared.prmtop'))
    coords=app.AmberInpcrdFile(str(output/'prepared.inpcrd')).positions
    system=loaded.createSystem(nonbondedMethod=app.NoCutoff,constraints=None,rigidWater=False)
    mechanical=validate_bonded_site(loaded.topology,system,manifest)
    (output/'mechanical-validation.json').write_text(json.dumps(mechanical,indent=2)+'\n')
    finite=check_finite_forces(system,coords.value_in_unit(unit.nanometer))
    expected_pairs={tuple(sorted([next(a.idx for a in native.atoms if a.residue.idx==zn.index and a.name=='ZN'),
        next(a.idx for a in native.atoms if a.residue.idx==r.index and a.name==name)])) for r,name in donors}
    actual_pairs={tuple(sorted([b.atom1.idx,b.atom2.idx])) for b in native.bonds if b.atom1.atomic_number==30 or b.atom2.atomic_number==30}
    if actual_pairs != expected_pairs: raise ValueError('Native Zn donor bond graph differs')
    site_residues=set(assignments.values());site_charge=sum(a.charge for a in native.atoms if a.residue.name in site_residues)
    if abs(site_charge-2)>1e-4: raise ValueError('Whole-site fitted charge differs from explicit +2 state')
    report={'status':'native model constructed; mechanical checks recorded separately', 'source_sha256':hashlib.sha256(source.read_bytes()).hexdigest(),
        'script_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),'parameters':json.loads((BASE/'reference-zaff/sources.json').read_text()),
        'model_state':state,'published_template_context_angstrom':context_distances,'template_assignments':{':'.join(k):v for k,v in assignments.items()},
        'observed_heavy_atoms_retained':len(matched),'source_observed_heavy_atoms':len(original),
        'maximum_native_serialization_error_angstrom':max_rounding,'native_atoms':len(native.atoms),'native_residues':len(native.residues),
        'native_Zn_bonds':[list(x) for x in sorted(actual_pairs)],'site_partial_charge_sum_e':site_charge,'total_system_charge_e':sum(a.charge for a in native.atoms),
        'mechanical_terms_match':mechanical['accepted'],'finite_configuration_check':finite,'dynamics':'not run',
        'interpretation':'Construction and same-parameter mechanical implementation checks only. No new site fitting, experimental validation, or app integration.'}
    (output/'preparation-report.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps({k:report[k] for k in ('native_atoms','observed_heavy_atoms_retained','mechanical_terms_match','finite_configuration_check')},indent=2))


if __name__=='__main__': main()
