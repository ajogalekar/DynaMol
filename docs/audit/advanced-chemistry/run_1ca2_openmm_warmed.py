"""Unrestrained NVT comparison from an explicitly reviewed warmed endpoint."""
from pathlib import Path
import argparse
import hashlib
import json
import time
import numpy as np
import openmm as mm
from openmm import app,unit
import parmed

BASE=Path(__file__).resolve().parent
p=argparse.ArgumentParser();p.add_argument('--input',type=Path,required=True);p.add_argument('--output',type=Path,required=True)
p.add_argument('--duration-ps',type=int,default=20);p.add_argument('--seed',type=int,default=3131)
a=p.parse_args();source=a.input.resolve();out=a.output.resolve()
if not 10<=a.duration_ps<=50:raise ValueError('Warmed comparison is bounded to10–50ps')
review=json.loads((source/'comparison-readiness.json').read_text())
if review.get('ready_for_short_comparison') is not True:raise ValueError('Endpoint trends have not been explicitly reviewed')
out.mkdir(exist_ok=False)
prmtop=BASE/'prototype-1ca2-zaff-motion/solvated.prmtop';native=parmed.load_file(str(prmtop))
end=np.load(source/'endpoint.npz');raw=end['coords_nm'].astype(float);box=end['box_nm'].astype(float);inverse=np.linalg.inv(box)
if len(native.atoms)!=len(raw) or not np.isfinite(raw).all() or np.linalg.det(box)<=0:raise ValueError('Invalid native endpoint inventory/box')
adj=[[] for _ in native.atoms]
for bond in native.bonds:
    i,j=bond.atom1.idx,bond.atom2.idx;adj[i].append(j);adj[j].append(i)
whole=np.full_like(raw,np.nan);components=0
def minimum_image(v):return v-np.rint(v@inverse)@box
for anchor in range(len(raw)):
    if np.isfinite(whole[anchor]).all():continue
    components+=1;whole[anchor]=raw[anchor];queue=[anchor]
    for i in queue:
        for j in adj[i]:
            expected=whole[i]+minimum_image(raw[j]-raw[i])
            if np.isnan(whole[j]).any():whole[j]=expected;queue.append(j)
            elif not np.allclose(whole[j],expected,atol=1e-6,rtol=0):raise ValueError('Inconsistent periodic bonded cycle')
images=(whole-raw)@inverse
if not np.allclose(images,np.rint(images),atol=1e-6,rtol=0):raise ValueError('Coordinate reconstruction was not an integer-image change')
np.savez_compressed(out/'initial-image-mapping.npz',raw_coords_nm=raw,whole_coords_nm=whole,box_nm=box,image_indices=np.rint(images).astype(int))
top=app.AmberPrmtopFile(str(prmtop),periodicBoxVectors=tuple(mm.Vec3(*v) for v in box)*unit.nanometer)
sys=top.createSystem(nonbondedMethod=app.PME,nonbondedCutoff=1*unit.nanometer,constraints=app.HBonds,rigidWater=True,ewaldErrorTolerance=1e-5)
dof=3*sum(sys.getParticleMass(i)>0*unit.dalton for i in range(sys.getNumParticles()))-sys.getNumConstraints()-3
volume=float(np.linalg.det(box));density=sum(atom.mass for atom in native.atoms)*.00166053906660/volume
integ=mm.LangevinMiddleIntegrator(300*unit.kelvin,1/unit.picosecond,.001*unit.picosecond);integ.setRandomNumberSeed(a.seed)
integ.setConstraintTolerance(1e-6)
(out/'system.xml').write_text(mm.XmlSerializer.serialize(sys))
sim=app.Simulation(top.topology,sys,integ,mm.Platform.getPlatformByName('CPU'),{'Threads':'2','DeterministicForces':'true'})
sim.context.setPositions(whole*unit.nanometer);sim.context.applyConstraints(1e-6)
constrained=np.array(sim.context.getState(getPositions=True).getPositions(asNumpy=True).value_in_unit(unit.nanometer))
projection=float(np.max(np.linalg.norm(constrained-whole,axis=1)))
sim.context.setVelocitiesToTemperature(300*unit.kelvin,a.seed)
sim.reporters.append(app.DCDReporter(str(out/'comparison.dcd'),1000));sim.reporters.append(app.StateDataReporter(str(out/'thermodynamics.csv'),100,
    step=True,time=True,potentialEnergy=True,kineticEnergy=True,temperature=True,volume=True,density=True))
zn=next(atom for atom in native.atoms if atom.residue.name=='ZN6')
bonds=[(bond.atom1.idx,bond.atom2.idx) for bond in native.bonds if zn in (bond.atom1,bond.atom2)]
if len(bonds)!=4:raise ValueError('Four native Zn bonds required')
angles=[(angle.atom1.idx,angle.atom2.idx,angle.atom3.idx) for angle in native.angles if angle.atom2.idx==zn.idx]
if len(angles)!=6:raise ValueError('Six native Zn-centered angles required')
ca=[atom.idx for atom in native.atoms if atom.name=='CA' and atom.residue.idx<256]
reference=whole[ca];q=reference-reference.mean(0);frames=[];metrics=[];started=time.monotonic()
report={'status':'running','scope':'Short unrestrained warmed NVT numerical/coordination comparison; not catalytic or thermodynamic-convergence validation',
    'duration_ps':a.duration_ps,'source_endpoint':str(source),'source_endpoint_sha256':hashlib.sha256((source/'endpoint.npz').read_bytes()).hexdigest(),
    'temperature_target_kelvin':300,'mobile_dof':dof,'volume_nm3':volume,'density_g_ml':density,
    'maximum_initial_constraint_projection_nm':projection,
    'image_mapping':'Whole molecules reconstructed by exact native bond graph; periodic cycle and integer-image checks passed','molecular_components':components,
    'restraints':False,'initial_velocities':'New Maxwell velocities and Langevin noise; independently seeded, not matched microscopic trajectories','seed':a.seed,
    'constraints':'Hydrogen bonds and rigid ordinary water; relative tolerance1e-6; Zn-heavy bonds remain dynamic',
    'native_Zn_bonds':bonds,'native_Zn_centered_angles':angles,
    'system_xml_sha256':hashlib.sha256((out/'system.xml').read_bytes()).hexdigest(),
    'native_prmtop_sha256':hashlib.sha256(prmtop.read_bytes()).hexdigest()}
def save():
    report['elapsed_seconds']=time.monotonic()-started;(out/'comparison-report.json').write_text(json.dumps(report,indent=2)+'\n')
save()
for step in range(0,a.duration_ps*1000+1,100):
    if step:sim.step(100)
    st=sim.context.getState(getPositions=True,getEnergy=True,getForces=True)
    xyz=np.asarray(st.getPositions(asNumpy=True).value_in_unit(unit.nanometer));f=st.getForces(asNumpy=True).value_in_unit(unit.kilojoule_per_mole/unit.nanometer)
    potential=st.getPotentialEnergy().value_in_unit(unit.kilojoule_per_mole);kinetic=st.getKineticEnergy().value_in_unit(unit.kilojoule_per_mole);temp=2*kinetic/(dof*.00831446261815324)
    if not np.isfinite(xyz).all() or not np.isfinite(f).all() or not np.isfinite([potential,kinetic,temp]).all():raise ValueError('Nonfinite complete system')
    distances=[float(np.linalg.norm(xyz[i]-xyz[j])*10) for i,j in bonds]
    angle_values=[]
    for i,j,k in angles:
        u=xyz[i]-xyz[j];v=xyz[k]-xyz[j]
        angle_values.append(float(np.degrees(np.arccos(np.clip(np.dot(u,v)/np.linalg.norm(u)/np.linalg.norm(v),-1,1)))))
    pc=xyz[ca]-xyz[ca].mean(0);u,_,v=np.linalg.svd(pc.T@q);sign=np.eye(3);sign[-1,-1]=np.linalg.det(u@v);rmsd=float(np.sqrt(np.mean(np.sum((pc@u@sign@v-q)**2,axis=1)))*10)
    metrics.append({'time_ps':step*.001,'potential_kj_mol':float(potential),'kinetic_kj_mol':float(kinetic),'temperature_kelvin':float(temp),
        'Zn_distances_angstrom':distances,'Zn_angles_degrees':angle_values,'protein_CA_aligned_rmsd_angstrom':rmsd,
        'maximum_force_kj_mol_nm':float(np.linalg.norm(f,axis=1).max()),'rms_force_kj_mol_nm':float(np.sqrt(np.mean(np.sum(f*f,axis=1))))})
    if step%1000==0:
        frames.append(xyz.astype(np.float32));report['last_sample']=metrics[-1];save()
        print(f'{step*.001:.1f}ps T={temp:.1f}K Zn='+','.join(f'{d:.3f}' for d in distances),flush=True)
    if time.monotonic()-started>3600:raise TimeoutError('Warmed comparison exceeded one hour')
d=np.array([m['Zn_distances_angstrom'] for m in metrics]);bounds=[(1.6,2.8) if any(native.atoms[j].name=='O' for j in pair) else (1.6,2.6) for pair in bonds]
ad=np.array([m['Zn_angles_degrees'] for m in metrics]);last=[m for m in metrics if m['time_ps']>=a.duration_ps-10]
thermal={}
for key in ['potential_kj_mol','temperature_kelvin']:
    values=np.array([m[key] for m in last]);times=np.array([m['time_ps'] for m in last])
    thermal[key]={'mean':float(values.mean()),'standard_deviation':float(values.std()),'slope_per_ps':float(np.polyfit(times-times[0],values,1)[0]),
        'first5ps_mean':float(values[times<times[0]+5].mean()),'last5ps_mean':float(values[times>=times[0]+5].mean())}
report.update(status='completed',finite_every_saved_force=True,metrics=metrics,temperature_mean_kelvin=float(np.mean([m['temperature_kelvin'] for m in metrics])),
    Zn_min_angstrom=d.min(0).tolist(),Zn_max_angstrom=d.max(0).tolist(),Zn_mean_angstrom=d.mean(0).tolist(),
    Zn_angle_min_degrees=ad.min(0).tolist(),Zn_angle_max_degrees=ad.max(0).tolist(),Zn_angle_mean_degrees=ad.mean(0).tolist(),last10ps_metrics=thermal,
    broad_geometry_alarms_clear=all(np.all((d[:,i]>=lo)&(d[:,i]<=hi)) for i,(lo,hi) in enumerate(bounds)),
    maximum_CA_rmsd_angstrom=max(m['protein_CA_aligned_rmsd_angstrom'] for m in metrics));save()
np.savez_compressed(out/'frames.npz',positions_nm=np.array(frames),times_ps=np.arange(len(frames)),box_nm=box)
print(json.dumps({k:report[k] for k in ['status','temperature_mean_kelvin','maximum_CA_rmsd_angstrom','elapsed_seconds']},indent=2))
