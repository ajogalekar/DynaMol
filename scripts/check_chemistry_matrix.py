"""Bounded representative chemistry checks; no coverage or convergence inference.

Uses an isolated data root. Writes incremental machine-readable outcomes, exact
inputs and native logs. Run with DYNAMOL_DATA_DIR pointing below data/integration-checks.
"""
from pathlib import Path
import datetime as dt
import hashlib
import importlib.metadata
import json
import os
import random
import sys
import traceback

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
import numpy as np
import openmm as mm
from openmm import app, unit
from backend import config, sources, storage
from backend import ligands
from backend.ions import inspect_ions
from backend.preparation_worker import protonation_inventory
from backend.complex_topology import subset, residue_key

OUT = ROOT / 'docs/audit/release-readiness'
WORK = config.DATA_ROOT / 'chemistry-matrix'
REPORT = OUT / 'chemistry-matrix.json'
SEED = 2026


def finite_run(topology, positions, forcefield):
    system = forcefield.createSystem(topology, nonbondedMethod=app.NoCutoff, constraints=None)
    integrator = mm.LangevinMiddleIntegrator(100, 1, .0002)
    integrator.setRandomNumberSeed(SEED)
    context = mm.Context(system, integrator, mm.Platform.getPlatformByName('Reference'))
    context.setPositions(positions)
    mm.LocalEnergyMinimizer.minimize(context, 10, 100)
    context.setVelocitiesToTemperature(100, SEED)
    integrator.step(5)
    state = context.getState(getPositions=True, getEnergy=True)
    energy = state.getPotentialEnergy().value_in_unit(unit.kilojoule_per_mole)
    xyz = np.asarray(state.getPositions(asNumpy=True).value_in_unit(unit.nanometer))
    assert np.isfinite(xyz).all() and np.isfinite(energy)
    charges = [f.getParticleParameters(i)[0].value_in_unit(unit.elementary_charge)
               for f in system.getForces() if isinstance(f, mm.NonbondedForce)
               for i in range(system.getNumParticles())]
    del context, integrator
    return {'atoms': system.getNumParticles(), 'charge_e': sum(charges), 'charges_e': charges,
            'finite_energy_kj_mol': energy, 'steps': 5, 'time_ps': .001,
            'note': 'Reference-platform, nonperiodic vacuum numerical continuity only; not aqueous or converged MD.'}


def ccd_dataset(component):
    block, reference = ligands._ccd(component)
    table = block.get_mmcif_category('_chem_comp_atom.')
    topology = app.Topology(); residue = topology.addResidue(component, topology.addChain('A'), '1')
    xyz=[]
    for i, name in enumerate(table['atom_id']):
        if table['type_symbol'][i] in {'H', 'D'}: continue
        element=app.element.Element.getBySymbol(table['type_symbol'][i].title())
        topology.addAtom(name, element, residue)
        xyz.append([float(table[f'pdbx_model_Cartn_{axis}_ideal'][i]) for axis in 'xyz'])
    folder=WORK/'inputs'; folder.mkdir(parents=True, exist_ok=True)
    path=folder/f'{component}.pdb'
    with path.open('w') as f: app.PDBFile.writeFile(topology, np.asarray(xyz)*unit.angstrom, f, keepIds=True)
    import shutil
    saved=OUT/'chemistry-inputs'; saved.mkdir(parents=True, exist_ok=True)
    shutil.copy2(config.DATA_ROOT/'chemistry/ccd'/f'{component}.cif', saved/f'{component}.cif')
    return sources.import_structure(path, name=f'CCD {component} isolated ideal-geometry fixture',
                                    provenance={'audit_source':reference,'audit_geometry':'CCD ideal coordinates; not bound experimental coordinates'})


def run_ligand(label, smiles=None):
    meta=sources.create_smiles(smiles, name=label, seed=SEED) if smiles else ccd_dataset(label)
    from backend.preparation import exact_input_path
    before=app.PDBFile(str(exact_input_path(meta['id'])))
    original_xyz=np.asarray(before.positions.value_in_unit(unit.nanometer))
    inspected=ligands.inspect_ligands(meta['id'], 7)
    assert len(inspected)==1 and inspected[0]['error'] is None, inspected
    prepared=ligands.prepare_ligands(meta['id'], 7, SEED, WORK/label,
        on_progress=lambda s: print(label,s,flush=True))[0]
    xyz=np.asarray(prepared.positions.value_in_unit(unit.nanometer))
    maximum=max(float(np.linalg.norm(xyz[r['ligand_index']]-original_xyz[r['original_index']])) for r in prepared.atom_map)
    assert maximum < 1e-6
    result=finite_run(prepared.topology,prepared.positions,app.ForceField(str(prepared.ffxml)))
    assert abs(result['charge_e']-inspected[0]['formal_charge']) < 1e-4
    result.pop('charges_e')
    return {'id':label,'kind':'organic ligand/cofactor','passed':True,'dataset_id':meta['id'],
            'reference':inspected[0]['reference'],'selected_smiles':inspected[0]['selected_smiles'],
            'formal_charge':inspected[0]['formal_charge'],'protonation_method':inspected[0]['protonation_method'],
            'warnings':inspected[0]['warnings'],'heavy_atoms':inspected[0]['heavy_atoms'],
            'heavy_coordinate_max_displacement_nm':maximum,'conversion_validation':prepared.provenance['conversion_validation'],
            'native_runtime':prepared.provenance['native_runtime'],'runtime':result,
            'artifacts':str(prepared.ffxml.parent.relative_to(ROOT))}


def standard_inputs():
    variants=['ALA','ARG','ASN','ASP','CYS','GLN','GLU','GLY','HID','HIE','HIP','ILE','LEU','LYS','MET','PHE','PRO','SER','THR','TRP','TYR','VAL','ASH','GLH','LYN','CYM']
    folder=WORK/'protein-native'; folder.mkdir(parents=True,exist_ok=True)
    lines=['source leaprc.protein.ff14SB']
    for variant in variants:
        lines += [f'p = sequence {{ ACE {variant} NME }}',f'saveamberparm p {variant}.prmtop {variant}.inpcrd',f'savepdb p {variant}.pdb']
    lines += ['p = sequence { NALA CGLY }','saveamberparm p termini.prmtop termini.inpcrd','savepdb p termini.pdb',
              'p = sequence { ACE CYX GLY CYX NME }','bond p.2.SG p.4.SG','saveamberparm p disulfide.prmtop disulfide.inpcrd','savepdb p disulfide.pdb','quit']
    (folder/'leap.in').write_text('\n'.join(lines)+'\n')
    amber=ligands.ambertools_home()
    env=dict(os.environ,AMBERHOME=str(amber),PATH=str(amber/'bin')+os.pathsep+os.environ.get('PATH',''))
    ligands._run([str(amber/'bin/tleap'),'-f','leap.in'],folder,env,None,None,timeout=60)
    return folder, variants+['termini','disulfide']


def run_standard(folder, name):
    native=app.AmberPrmtopFile(str(folder/f'{name}.prmtop'))
    pos=app.AmberInpcrdFile(str(folder/f'{name}.inpcrd')).positions
    ff=app.ForceField('amber14/protein.ff14SB.xml','amber14/tip3p.xml')
    result=finite_run(native.topology,pos,ff)
    charges=result.pop('charges_e')
    native_system=native.createSystem(nonbondedMethod=app.NoCutoff, constraints=None)
    native_charges=[f.getParticleParameters(i)[0].value_in_unit(unit.elementary_charge) for f in native_system.getForces() if isinstance(f,mm.NonbondedForce) for i in range(native_system.getNumParticles())]
    np.testing.assert_allclose(charges,native_charges,atol=1e-7,rtol=0)
    result['native_per_atom_charges_match']=True
    states=protonation_inventory(native.topology, [], charges)
    # Exercise DynaMol's hydrogen removal followed by OpenMM pH assignment on
    # the canonical PDB topology. Disulfide topology explicitly retains S-S.
    ph=2 if name in {'ASH','GLH','HIP'} else 12 if name in {'LYN','CYM'} else 7
    pdb_text=(folder/f'{name}.pdb').read_text()
    import io
    top=app.PDBFile(io.StringIO(pdb_text))
    from backend.preparation import canonicalize_protonation_aliases
    canonicalize_protonation_aliases(top.topology)
    if name=='disulfide':
        sg=[a for a in top.topology.atoms() if a.name=='SG']
        if not any({a,b}==set(sg) for a,b in top.topology.bonds()): top.topology.addBond(*sg)
    model=subset(top.topology,top.positions,{residue_key(r) for r in top.topology.residues()},remove_hydrogens=True)
    random.seed(SEED)
    variants=model.addHydrogens(ff,pH=ph,platform=mm.Platform.getPlatformByName('Reference'))
    regenerated=ff.createSystem(model.topology,nonbondedMethod=app.NoCutoff)
    rech=[f.getParticleParameters(i)[0].value_in_unit(unit.elementary_charge) for f in regenerated.getForces() if isinstance(f,mm.NonbondedForce) for i in range(regenerated.getNumParticles())]
    reassigned=protonation_inventory(model.topology,variants,rech)
    return {'id':name,'kind':'protein-template','passed':True,'native_template_states':states,
            'pH_reassignment':{'ph':ph,'states':reassigned,'charge_e':sum(rech)},'runtime':result}


def run_ion(name,symbol,charge):
    top=app.Topology(); residue=top.addResidue(name,top.addChain('I'),'1'); top.addAtom(name,app.element.Element.getBySymbol(symbol),residue)
    rows=inspect_ions(top); assert len(rows)==1 and rows[0]['supported']
    result=finite_run(top,np.zeros((1,3))*unit.nanometer,app.ForceField('amber14/tip3p.xml'))
    assert abs(result['charge_e']-charge)<1e-8
    result.pop('charges_e')
    return {'id':name,'kind':'ion','passed':True,'inspection':rows[0],'runtime':result}


def main():
    if config.DATA_ROOT == ROOT/'data': raise SystemExit('Use isolated DYNAMOL_DATA_DIR to protect the user library.')
    WORK.mkdir(parents=True,exist_ok=True); OUT.mkdir(parents=True,exist_ok=True)
    report={'created_at':dt.datetime.now(dt.timezone.utc).isoformat(),'seed':SEED,'scope':'Representative software/numerical checks, not universal chemistry or biological accuracy. Native artifacts retained in isolated data root.','rows':[],
            'versions':{name:importlib.metadata.version(name) for name in ['openmm','rdkit','pdbfixer','mdtraj','numpy','parmed','gemmi','dimorphite-dl']},
            'script_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest()}
    def run(label, fn):
        print('CHECK',label,flush=True)
        try: row=fn()
        except Exception as exc: row={'id':label,'passed':False,'error':str(exc),'traceback':traceback.format_exc()}
        report['rows'].append(row); report['passed']=all(r['passed'] for r in report['rows'])
        storage.atomic_json(REPORT,report)
        print('RESULT',label,row['passed'],row.get('error',''),flush=True)
    native, names=standard_inputs()
    report['native_source']={'runtime':ligands._amber_versions(),
        'files':{str(path.relative_to(ROOT)):hashlib.sha256(path.read_bytes()).hexdigest() for path in sorted(native.iterdir()) if path.suffix in {'.pdb','.prmtop','.inpcrd','.in'}}}
    for name in names: run(name,lambda name=name:run_standard(native,name))
    for name,symbol,charge in [('NA','Na',1),('CL','Cl',-1),('K','K',1),('MG','Mg',2),('CA','Ca',2)]:
        run(name,lambda n=name,s=symbol,c=charge:run_ion(n,s,c))
    for label,smiles in [('ethanol','CCO'),('acetate','CC(=O)O'),('methylammonium','CN'),('halogen-aromatic','Fc1cc(Cl)cc(Br)c1'),('phosphate','COP(=O)(O)O'),('chiral-amide','CC(=O)N[C@@H](C)CN(C)C')]:
        run(label,lambda l=label,s=smiles:run_ligand(l,s))
    for name in ['ATP','ADP','NAD','FAD']: run(name,lambda name=name:run_ligand(name))
    report['finished_at']=dt.datetime.now(dt.timezone.utc).isoformat()
    storage.atomic_json(REPORT,report)
    print(json.dumps({'passed':report['passed'],'cases':len(report['rows'])}))

if __name__=='__main__': main()
