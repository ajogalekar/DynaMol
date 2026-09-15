"""Independent Sander/OpenMM static energy+force checks of an Amber topology.

An engineering parity test, not a QM/experimental accuracy or MD stability test.
The Sander subprocess runs in the bundled AmberTools Python runtime.
"""
from __future__ import annotations
import argparse
import hashlib
import json
import subprocess
import sys
from pathlib import Path
import numpy as np


def sander_values(request, output):
    import sander
    d=json.loads(request.read_text());options=sander.gas_input()
    # The native API selects igb=6 for vacuum: forcing igb=0 switches to its
    # non-GB/Ewald setup path and is not the documented gas_input default.
    if options.igb!=6:raise ValueError('Unexpected native vacuum configuration.')
    options.ntb=0;options.cut=999.0
    if len(d['poses']) != 1:
        raise ValueError('Each native process must evaluate exactly one pose; CMAP setup is not safely reusable.')
    result=[]
    for pose in d['poses']:
        with sander.setup(d['prmtop'],np.asarray(pose['angstrom']).ravel(),None,options):
            energy,forces=sander.energy_forces()
            components={name:float(getattr(energy,name))*4.184 for name in dir(energy)
                        if not name.startswith('_') and isinstance(getattr(energy,name),(float,int))}
            result.append({'name':pose['name'],'energy_kj_mol':float(energy.tot)*4.184,
                           'components_kj_mol':components,
                           'forces_kj_mol_nm':(np.asarray(forces).reshape(-1,3)*41.84).tolist()})
    output.write_text(json.dumps(result)+'\n')


def isolated_sander_values(prmtop, poses, folder, runtime, helper=None):
    """Keep native setup/cleanup state in one bounded subprocess per pose."""
    helper=Path(helper or __file__).resolve();results=[]
    for index,pose in enumerate(poses):
        request=folder/f'sander-pose-{index:02d}-input.json'
        output=folder/f'sander-pose-{index:02d}-native.json'
        if request.exists() or output.exists():
            raise FileExistsError('Preserve prior native evidence; select a new report directory.')
        request.write_text(json.dumps({'prmtop':str(prmtop),'poses':[pose]})+'\n')
        run=subprocess.run([str(runtime),str(helper),str(request),'--sander','--output',str(output)],
                           capture_output=True,text=True,timeout=120)
        (folder/f'sander-pose-{index:02d}.log').write_text(run.stdout+'\n'+run.stderr)
        if run.returncode:
            raise RuntimeError(f'Independent Sander pose {index} failed; native log preserved.')
        values=json.loads(output.read_text())
        if len(values)!=1 or values[0]['name']!=pose['name']:
            raise ValueError('Native result is not the requested pose.')
        results.extend(values)
    return results


def main():
    parser=argparse.ArgumentParser();parser.add_argument('folder',type=Path)
    parser.add_argument('--sander',action='store_true');parser.add_argument('--output',type=Path)
    parser.add_argument('--stem',default='ligand')
    parser.add_argument('--report-dir',type=Path)
    args=parser.parse_args()
    if args.sander:
        sander_values(args.folder,args.output);return
    import openmm as mm
    from openmm import app,unit
    folder=args.folder.resolve();root=Path(__file__).resolve().parents[3]
    prmtop=folder/(args.stem+'.prmtop');inpcrd=folder/(args.stem+'.inpcrd')
    native=app.AmberPrmtopFile(str(prmtop));positions=app.AmberInpcrdFile(str(inpcrd)).positions
    xyz=np.asarray(positions.value_in_unit(unit.angstrom));poses=[{'name':'bound','angstrom':xyz.tolist()}]
    for seed in (2027,2028):
        poses.append({'name':f'perturb-{seed}','angstrom':(xyz+np.random.default_rng(seed).normal(0,.003,xyz.shape)).tolist()})
    report_folder=args.report_dir.resolve() if args.report_dir else folder
    report_folder.mkdir(parents=True,exist_ok=True)
    request=report_folder/'sander-parity-input.json';output=report_folder/'sander-parity-native.json'
    if request.exists() or output.exists():raise FileExistsError('Select a fresh report directory to preserve prior evidence.')
    request.write_text(json.dumps({'prmtop':str(prmtop),'poses':poses})+'\n')
    reference=isolated_sander_values(prmtop,poses,report_folder,root/'.tools/ambertools/bin/python')
    output.write_text(json.dumps(reference)+'\n')
    system=native.createSystem(nonbondedMethod=app.NoCutoff,constraints=None,rigidWater=False)
    # Measure OpenMM's electrostatic constant with a separate two-particle
    # system. Amber's stored charges use its legacy18.2223 scale. Do not infer
    # a correction from the molecule being tested or change production files.
    probe=mm.System();probe.addParticle(1);probe.addParticle(1);force=mm.NonbondedForce()
    force.addParticle(1,1,0);force.addParticle(1,1,0);probe.addForce(force)
    probe_integrator=mm.VerletIntegrator(.001);probe_context=mm.Context(probe,probe_integrator,mm.Platform.getPlatformByName('Reference'))
    probe_context.setPositions([[0,0,0],[1,0,0]])
    openmm_coulomb=probe_context.getState(getEnergy=True).getPotentialEnergy().value_in_unit(unit.kilojoule_per_mole)
    del probe_context,probe_integrator
    amber_coulomb=18.2223**2*4.184/10
    ratio=amber_coulomb/openmm_coulomb
    normalized=mm.XmlSerializer.deserialize(mm.XmlSerializer.serialize(system))
    nb=next(f for f in normalized.getForces() if isinstance(f,mm.NonbondedForce))
    for i in range(nb.getNumParticles()):
        charge,sigma,epsilon=nb.getParticleParameters(i);nb.setParticleParameters(i,charge*np.sqrt(ratio),sigma,epsilon)
    for i in range(nb.getNumExceptions()):
        a,b,charge_product,sigma,epsilon=nb.getExceptionParameters(i)
        nb.setExceptionParameters(i,a,b,charge_product*ratio,sigma,epsilon)
    integrator=mm.VerletIntegrator(.001);context=mm.Context(system,integrator,mm.Platform.getPlatformByName('Reference'))
    other_integrator=mm.VerletIntegrator(.001);other_context=mm.Context(normalized,other_integrator,mm.Platform.getPlatformByName('Reference'))
    records=[]
    for pose,expected in zip(poses,reference):
        context.setPositions(np.asarray(pose['angstrom'])*unit.angstrom)
        state=context.getState(getEnergy=True,getForces=True)
        energy=float(state.getPotentialEnergy().value_in_unit(unit.kilojoule_per_mole))
        forces=np.asarray(state.getForces(asNumpy=True).value_in_unit(unit.kilojoule_per_mole/unit.nanometer))
        de=energy-expected['energy_kj_mol'];df=float(np.max(np.abs(forces-np.asarray(expected['forces_kj_mol_nm']))))
        other_context.setPositions(np.asarray(pose['angstrom'])*unit.angstrom)
        normalized_state=other_context.getState(getEnergy=True,getForces=True)
        normalized_energy=normalized_state.getPotentialEnergy().value_in_unit(unit.kilojoule_per_mole)
        normalized_forces=np.asarray(normalized_state.getForces(asNumpy=True).value_in_unit(unit.kilojoule_per_mole/unit.nanometer))
        normalized_de=normalized_energy-expected['energy_kj_mol'];normalized_df=float(np.max(np.abs(normalized_forces-np.asarray(expected['forces_kj_mol_nm']))))
        records.append({'pose':pose['name'],'sander_energy_kj_mol':expected['energy_kj_mol'],
                        'openmm_energy_kj_mol':energy,'energy_difference_kj_mol':de,
                        'maximum_force_difference_kj_mol_nm':df,
                        'raw_constants_agree_within_tolerance':bool(np.isfinite([energy,de,df]).all() and abs(de)<1e-3 and df<1e-2),
                        'matched_constant_energy_difference_kj_mol':float(normalized_de),
                        'matched_constant_maximum_force_difference_kj_mol_nm':normalized_df,
                        'passed':bool(np.isfinite([normalized_energy,normalized_de,normalized_df]).all() and abs(normalized_de)<1e-3 and normalized_df<1e-2)})
    del context,integrator,other_context,other_integrator
    report={'method':'Independent AmberTools Sander API vacuum (igb6,ntb0,cut999Å) versus OpenMM Reference AmberPrmtopFile NoCutoff, no constraints; bound coordinates and two0.003Å perturbations.',
            'scope':'Checks engine implementation and unit conversion of the same classical Hamiltonian, not its chemical accuracy.',
            'energy_tolerance_kj_mol':1e-3,'force_tolerance_kj_mol_nm':1e-2,
            'electrostatic_convention':{'amber_legacy_charge_scale':18.2223,
                'amber_coulomb_kj_mol_nm_e2':amber_coulomb,'openmm_coulomb_kj_mol_nm_e2':openmm_coulomb,
                'energy_scale_ratio':ratio,'method':'Diagnostic copy only: particle charges and1–4 products scaled to match Amber legacy electrostatic convention. Both raw and convention-matched errors recorded; saved production parameters untouched.'},
            'input_sha256':{p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in (prmtop,inpcrd)},
            'conformations':records,'passed':all(x['passed'] for x in records)}
    (report_folder/'sander-parity.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(report,indent=2))
    if not report['passed']:raise ValueError('Sander/OpenMM parity failed.')


if __name__=='__main__':main()
