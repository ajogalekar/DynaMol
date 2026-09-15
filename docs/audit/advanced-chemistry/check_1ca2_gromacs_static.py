"""Actual mixed-precision GROMACS energies/forces versus both source readers.

All native terms are checked separately; numerical agreement is not a chemical
model validation. Use identical saved coordinates and retain raw errors.
"""
from pathlib import Path
import hashlib
import json
import subprocess
import time
import numpy as np
import mdtraj
import parmed
import openmm as mm
from openmm import app,unit,Vec3

BASE=Path(__file__).resolve().parent
ROOT=BASE.parents[2]
folder=BASE/'prototype-1ca2-gromacs-static'
source=BASE/'prototype-1ca2-zaff-water-v3'
exe=ROOT/'.gromacs/bin/gmx'
native=parmed.load_file(str(source/'prepared.prmtop'),xyz=str(source/'prepared.inpcrd'))
n=len(native.atoms); native.box=[200,200,200,90,90,90]
original=native.coordinates*.1+5.
relaxed=np.load(BASE/'prototype-1ca2-zaff-motion/seed-2026-frames.npz')['positions_angstrom'][0,:n].astype(float)*.1
poses=[('raw-centered',original),('relaxed',relaxed)]+[(f'relaxed-perturb-{seed}',relaxed+np.random.default_rng(seed).normal(0,.0003,relaxed.shape)) for seed in (2027,2028)]
systems=[]
for label,top in [('native_amber',app.AmberPrmtopFile(str(source/'prepared.prmtop'))),('gromacs_export',app.GromacsTopFile(str(folder/'system.top'),
        periodicBoxVectors=(Vec3(20,0,0),Vec3(0,20,0),Vec3(0,0,20))*unit.nanometer,defines={'FLEXIBLE':''}))]:
    sys=top.createSystem(nonbondedMethod=app.NoCutoff,constraints=None,rigidWater=False)
    integ=mm.VerletIntegrator(.001)
    systems.append((label,sys,integ,mm.Context(sys,integ,mm.Platform.getPlatformByName('Reference'))))
records=[];started=time.monotonic()
for label,coordinates in poses:
    out=folder/label;out.mkdir(exist_ok=False)
    native.coordinates=coordinates*10;native.save(str(out/'input.gro'),precision=7)
    def run(args,input=None):
        proc=subprocess.run([str(exe)]+args,cwd=out,input=input,capture_output=True,text=True,timeout=90)
        (out/(args[0]+'.log')).write_text(proc.stdout+'\n'+proc.stderr)
        if proc.returncode:raise RuntimeError(f'{label}: {args[0]} failed; original output retained')
    run(['grompp','-f','../static.mdp','-c','input.gro','-p','../system.top','-o','run.tpr'])
    run(['mdrun','-s','run.tpr','-deffnm','run','-ntmpi','1','-ntomp','2','-nb','cpu'])
    run(['energy','-f','run.edr','-o','potential.xvg','-xvg','none'],'Potential\n0\n')
    energy=float(np.loadtxt(out/'potential.xvg')[1])
    with mdtraj.formats.TRRTrajectoryFile(str(out/'run.trr')) as f:
        arrays=f._read(1,None,get_forces=True)
    xyz=arrays[0][0].astype(float);forces=arrays[-1][0].astype(float)
    if np.max(np.abs(xyz-coordinates))>2e-6:raise ValueError('Native static coordinates changed beyond float32 serialization')
    np.savez(out/'native-arrays.npz',coords_nm=xyz,forces_kj_mol_nm=forces,energy_kj_mol=energy)
    record={'pose':label,'gromacs_energy_kj_mol':energy,'native_precision':'mixed, recorded positions and forces float32',
        'maximum_coordinate_serialization_error_nm':float(np.max(np.abs(xyz-coordinates))),'readers':[]}
    raw=[]
    for reader,_,_,ctx in systems:
        ctx.setPositions(xyz*unit.nanometer)
        state=ctx.getState(getEnergy=True,getForces=True)
        e=state.getPotentialEnergy().value_in_unit(unit.kilojoule_per_mole)
        f=np.asarray(state.getForces(asNumpy=True).value_in_unit(unit.kilojoule_per_mole/unit.nanometer))
        if not np.isfinite([e,energy]).all() or not np.isfinite(f).all() or not np.isfinite(forces).all():raise ValueError('Nonfinite native/static comparison')
        delta=f-forces;maxref=float(np.max(np.abs(f)));rmsref=float(np.sqrt(np.mean(f*f)))
        record['readers'].append({'reader':reader,'energy_kj_mol':e,'energy_difference_kj_mol':float(e-energy),
            'energy_relative_difference':abs(e-energy)/max(1,abs(e)),
            'maximum_force_difference_kj_mol_nm':float(np.max(np.abs(delta))),
            'rms_force_difference_kj_mol_nm':float(np.sqrt(np.mean(delta*delta))),
            'force_maximum_reference_kj_mol_nm':maxref,'force_rms_reference_kj_mol_nm':rmsref,
            'force_maximum_difference_normalized_to_maximum':float(np.max(np.abs(delta)))/max(1,maxref),
            'force_rms_relative_error':float(np.sqrt(np.mean(delta*delta)))/max(1,rmsref)})
        raw.append((e,f))
    record['double_precision_source_vs_export']={'energy_difference_kj_mol':float(raw[1][0]-raw[0][0]),
        'maximum_force_difference_kj_mol_nm':float(np.max(np.abs(raw[1][1]-raw[0][1])))}
    records.append(record)
    print(label,json.dumps(record['readers'][0]),flush=True)
report={'scope':'Native GROMACS mixed-precision evaluation, original Amber/OpenMM Reference and exported GROMACS/OpenMM Reference on identical saved coordinates',
    'source_model':'Complete 1CA2 published ZAFF6 + ff14SB/TIP3P including all crystal waters',
    'forcefield_terms_or_charges_modified':False,'numerical_setup':'20nm box,9nm unshifted cutoffs encompass every pair without interacting with periodic images; no constraints and FLEXIBLE waters',
    'time_integration':'none; actual nsteps0 static evaluations','parameter_mechanics_report':'gromacs-reader-native-mechanics-remapped.json',
    'interpretation':'Raw numerical errors are reported; double-precision parameter conversion and mixed-precision native arithmetic are distinct checks. No tolerance adjusted in the application.',
    'conformations':records,'elapsed_seconds':time.monotonic()-started,'hashes':{n:hashlib.sha256((folder/n).read_bytes()).hexdigest() for n in ['system.top','static.mdp']}}
(folder/'native-static-matrix.json').write_text(json.dumps(report,indent=2)+'\n')
