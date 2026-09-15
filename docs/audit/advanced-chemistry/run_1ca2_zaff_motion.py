"""Bounded, explicitly unrestrained water-box mechanics check of native 1CA2.

This does not establish thermodynamic convergence or chemical model accuracy.
Two velocities/noise seeds use the same minimized full native model.
"""
from __future__ import annotations
import argparse
import hashlib
import json
import os
from pathlib import Path
import shutil
import subprocess
import time

import numpy as np
import openmm as mm
from openmm import app, unit
import parmed

BASE=Path(__file__).resolve().parent
ROOT=BASE.parents[2]


def main():
    parser=argparse.ArgumentParser();parser.add_argument('--input',type=Path,default=BASE/'prototype-1ca2-zaff-water-v3')
    parser.add_argument('--output',type=Path,default=BASE/'prototype-1ca2-zaff-motion')
    parser.add_argument('--steps',type=int,default=5000)
    args=parser.parse_args();source=args.input.resolve();output=args.output.resolve();output.mkdir(exist_ok=False)
    if not json.loads((source/'mechanical-validation.json').read_text())['accepted']:
        raise ValueError('Native dry parameter comparison did not pass')
    if not 1000<=args.steps<=10000:raise ValueError('This bounded audit permits 1–10 ps per seed')
    for name in ('ZAFF.prep','ZAFF.frcmod','prepared-input.pdb'):shutil.copy2(source/name,output/name)
    setup=(source/'tleap.in').read_text().split('check mol')[0]
    setup+='solvatebox mol TIP3PBOX 10.0\ncheck mol\nsavepdb mol solvated.pdb\nsaveamberparm mol solvated.prmtop solvated.inpcrd\nquit\n'
    (output/'tleap.in').write_text(setup)
    amber=ROOT/'.tools/ambertools';env=dict(os.environ,AMBERHOME=str(amber),OMP_NUM_THREADS='2',OPENBLAS_NUM_THREADS='2')
    with (output/'tleap.log').open('w') as stream:
        proc=subprocess.run([str(amber/'bin/tleap'),'-s','-f','tleap.in'],cwd=output,env=env,stdout=stream,stderr=subprocess.STDOUT,timeout=120)
    if proc.returncode or not (output/'solvated.prmtop').exists():raise ValueError('Native solvation failed; see retained log')
    dry=parmed.load_file(str(source/'prepared.prmtop'),xyz=str(source/'prepared.inpcrd'))
    native=parmed.load_file(str(output/'solvated.prmtop'),xyz=str(output/'solvated.inpcrd'))
    if [(a.name,a.residue.name,a.atomic_number) for a in native.atoms[:len(dry.atoms)]] != [(a.name,a.residue.name,a.atomic_number) for a in dry.atoms]:
        raise ValueError('Solvation changed the original atom inventory/order')
    shift=native.coordinates[:len(dry.atoms)]-dry.coordinates
    delta=np.max(np.abs(shift-shift[0]))
    if delta>1e-6:raise ValueError('Solvation changed retained coordinates beyond a common box translation')
    top=app.AmberPrmtopFile(str(output/'solvated.prmtop'));inp=app.AmberInpcrdFile(str(output/'solvated.inpcrd'))
    system=top.createSystem(nonbondedMethod=app.PME,nonbondedCutoff=1.0*unit.nanometer,constraints=app.HBonds,
        rigidWater=True,ewaldErrorTolerance=1e-5)
    if any(f.__class__.__name__ not in {'HarmonicBondForce','HarmonicAngleForce','PeriodicTorsionForce','NonbondedForce','CMMotionRemover'} for f in system.getForces()):
        raise ValueError('Unexpected force class; no positional restraint is allowed in this specimen')
    zn=next(a for a in native.atoms if a.residue.name=='ZN6')
    site=[{'a':b.atom1.idx,'b':b.atom2.idx,'donor':b.atom2.name if b.atom1.idx==zn.idx else b.atom1.name,
        'residue':b.atom2.residue.name if b.atom1.idx==zn.idx else b.atom1.residue.name,'equilibrium_angstrom':b.type.req}
        for b in native.bonds if zn in (b.atom1,b.atom2)]
    if len(site)!=4:raise ValueError('Expected four physical Zn–donor bonds')
    ca=[a.idx for a in native.atoms if a.name=='CA' and a.residue.idx<256]
    platform=mm.Platform.getPlatformByName('CPU');properties={'Threads':'2','DeterministicForces':'true'}
    started=time.monotonic();records=[];minimized=None
    report={'model':'1CA2 exact published ZAFF center6, ff14SB/TIP3P; explicit neutral bound water',
        'scope':'Short classical mechanics exercise; not experimental accuracy, catalytic validity, or converged sampling',
        'source_input':str(source),'original_atoms_retained':len(dry.atoms),'total_atoms':len(native.atoms),
        'total_residues':len(native.residues),'translation_angstrom':shift[0].tolist(),'translation_residual_angstrom':float(delta),
        'temperature_kelvin':300,'time_step_ps':.001,'steps_per_seed':args.steps,'duration_ps_per_seed':args.steps*.001,
        'solvent_box_padding_angstrom':10,'restraints':False,'constraints':'H bonds and rigid solvent; Zn donor heavy bonds remain dynamic',
        'platform':'CPU','platform_properties':properties,'seeds':[2026,2027],'site_bonds':site,'runs':records,'status':'running',
        'source_hashes':{n:hashlib.sha256((output/n).read_bytes()).hexdigest() for n in ('solvated.prmtop','solvated.inpcrd')}}
    def save():
        report['elapsed_seconds']=time.monotonic()-started
        (output/'motion-report.json').write_text(json.dumps(report,indent=2)+'\n')
    save()
    for seed in report['seeds']:
        print(f'seed {seed}: initialize',flush=True)
        integrator=mm.LangevinMiddleIntegrator(300*unit.kelvin,1/unit.picosecond,.001*unit.picosecond)
        integrator.setRandomNumberSeed(seed)
        simulation=app.Simulation(top.topology,system,integrator,platform,properties)
        simulation.context.setPeriodicBoxVectors(*inp.boxVectors)
        simulation.context.setPositions(inp.positions if minimized is None else minimized)
        if minimized is None:
            print('minimize complete unrestrained system (at most 1,000 iterations)',flush=True)
            simulation.minimizeEnergy(tolerance=10*unit.kilojoule_per_mole/unit.nanometer,maxIterations=1000)
            minimized=simulation.context.getState(getPositions=True).getPositions()
            with (output/'minimized.pdb').open('w') as stream:app.PDBFile.writeFile(top.topology,minimized,stream)
        simulation.context.setVelocitiesToTemperature(300*unit.kelvin,seed)
        simulation.reporters.append(app.DCDReporter(str(output/f'seed-{seed}.dcd'),100))
        simulation.reporters.append(app.StateDataReporter(str(output/f'seed-{seed}.csv'),100,step=True,time=True,
            potentialEnergy=True,kineticEnergy=True,temperature=True,speed=True))
        frames=[];metrics=[];reference=np.asarray(minimized.value_in_unit(unit.angstrom))[ca]
        for step in range(0,args.steps+1,100):
            if step:simulation.step(100)
            state=simulation.context.getState(getPositions=True,getEnergy=True,getForces=True)
            xyz=np.asarray(state.getPositions(asNumpy=True).value_in_unit(unit.angstrom))
            energy=float(state.getPotentialEnergy().value_in_unit(unit.kilojoule_per_mole))
            force=np.asarray(state.getForces(asNumpy=True).value_in_unit(unit.kilojoule_per_mole/unit.nanometer))
            if not np.isfinite(xyz).all() or not np.isfinite(force).all() or not np.isfinite(energy):
                report['status']='nonfinite_failure';save();raise ValueError('Nonfinite full system during motion')
            distances=[float(np.linalg.norm(xyz[b['a']]-xyz[b['b']])) for b in site]
            p=xyz[ca]-xyz[ca].mean(0);q=reference-reference.mean(0);u,_,v=np.linalg.svd(p.T@q)
            sign=np.eye(3);sign[-1,-1]=np.linalg.det(u@v);rotation=u@sign@v
            rmsd=float(np.sqrt(np.mean(np.sum((p@rotation-q)**2,axis=1))))
            metrics.append({'time_ps':step*.001,'potential_energy_kj_mol':energy,
                'maximum_force_kj_mol_nm':float(np.linalg.norm(force,axis=1).max()),'protein_CA_aligned_rmsd_angstrom':rmsd,
                'Zn_donor_distances_angstrom':distances})
            frames.append(xyz.astype(np.float32))
            if step%1000==0: print(f'seed {seed}: {step*.001:.1f} ps; Zn distances '+','.join(f'{d:.3f}' for d in distances),flush=True)
            if time.monotonic()-started>1800:raise TimeoutError('Bounded motion audit exceeded 30 minutes')
        np.savez_compressed(output/f'seed-{seed}-frames.npz',positions_angstrom=np.asarray(frames),times_ps=np.arange(len(frames))*.1)
        distances=np.asarray([m['Zn_donor_distances_angstrom'] for m in metrics])
        bounds=[[1.6,2.8] if b['donor']=='O' else [1.6,2.6] for b in site]
        geometry=all(np.all((distances[:,i]>=lo)&(distances[:,i]<=hi)) for i,(lo,hi) in enumerate(bounds))
        records.append({'seed':seed,'finite_all_saved_frames':True,'frame_count':len(frames),'screening_bounds_angstrom':bounds,
            'bounds_meaning':'Declared broad geometry alarms for this short Zn model exercise; not confidence intervals or scientific validation',
            'Zn_geometry_within_screening_bounds':bool(geometry),'Zn_distance_min_angstrom':distances.min(0).tolist(),
            'Zn_distance_max_angstrom':distances.max(0).tolist(),'Zn_distance_mean_angstrom':distances.mean(0).tolist(),
            'maximum_CA_rmsd_angstrom':max(m['protein_CA_aligned_rmsd_angstrom'] for m in metrics),'metrics':metrics})
        save();del simulation;del integrator
    report['status']='completed';report['short_motion_checks_passed']=all(r['finite_all_saved_frames'] and r['Zn_geometry_within_screening_bounds'] for r in records);save()
    print(json.dumps({k:report[k] for k in ('status','short_motion_checks_passed','elapsed_seconds','total_atoms')},indent=2))


if __name__=='__main__':main()
