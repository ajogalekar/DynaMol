"""One bounded warm-up or comparison stage; always preserve previous stages."""
from pathlib import Path
import argparse
import hashlib
import json
import re
import shutil
import subprocess
import time
import numpy as np
import mdtraj
import parmed

BASE=Path(__file__).resolve().parent;ROOT=BASE.parents[2]
p=argparse.ArgumentParser();p.add_argument('--input',type=Path,required=True);p.add_argument('--input-stem',default='motion')
p.add_argument('--output',type=Path,required=True);p.add_argument('--duration-ps',type=int,required=True)
p.add_argument('--ensemble',choices=['npt','nvt'],required=True);p.add_argument('--seed',type=int,default=3030)
a=p.parse_args();source=a.input.resolve();out=a.output.resolve();out.mkdir(exist_ok=False)
if not 10<=a.duration_ps<=100:raise ValueError('A stage is bounded to10–100ps')
for src,dst in [('system.top','system.top'),(a.input_stem+'.gro','input.gro'),(a.input_stem+'.cpt','input.cpt')]:shutil.copy2(source/src,out/dst)
native=parmed.load_file(str(out/'system.top'));mass=sum(atom.mass for atom in native.atoms)
if not any(x.residue.name=='ZN6' for x in native.atoms):raise ValueError('Expected explicit ZAFF Zn model absent')
mdp={}
for line in (BASE/'prototype-1ca2-gromacs-motion-v2/motion.mdp').read_text().splitlines():
    if '=' in line:
        k,v=line.split('=',1);mdp[k.strip()]=v.strip()
mdp.update({'integrator':'md','nsteps':str(a.duration_ps*1000),'nstxout':'1000','nstfout':'1000','nstlog':'1000',
    'tcoupl':'v-rescale','tau-t':'1.0','nsttcouple':'10','continuation':'yes','gen-vel':'no','ld-seed':str(a.seed),
    'pcoupl':'C-rescale' if a.ensemble=='npt' else 'no'})
mdp.pop('gen-temp',None);mdp.pop('gen-seed',None)
if a.ensemble=='npt':mdp.update({'pcoupltype':'isotropic','ref-p':'1.0','tau-p':'2.0','compressibility':'4.5e-5','nstpcouple':'10'})
(out/'stage.mdp').write_text('\n'.join(f'{k} = {v}' for k,v in mdp.items())+'\n')
report={'status':'running','stage_duration_ps':a.duration_ps,'ensemble':a.ensemble,'target_temperature_kelvin':300,
    'target_pressure_bar':1 if a.ensemble=='npt' else None,'restraints':False,
    'thermostat':'native stochastic velocity rescaling; tau1ps','barostat':'native stochastic cell rescaling; tau2ps' if a.ensemble=='npt' else None,
    'mobile_dof':'read from native grompp after constraints','initial_source':str(source),'initial_checkpoint':a.input_stem+'.cpt',
    'model':'unchanged native ff14SB/ZAFF6/TIP3P','atoms':len(native.atoms),'mass_da':mass,
    'meaning':'Thermal/density warm-up or comparison specimen. Completion alone is not convergence or model validation.',
    'source_hashes':{n:hashlib.sha256((out/n).read_bytes()).hexdigest() for n in ['system.top','input.gro','input.cpt','stage.mdp']}}
(out/'stage-status.json').write_text(json.dumps(report,indent=2)+'\n');started=time.monotonic();exe=ROOT/'.gromacs/bin/gmx'
def run(args,timeout,input=None):
    proc=subprocess.run([str(exe)]+args,cwd=out,input=input,text=True,stdout=(out/(args[0]+'.log')).open('w'),stderr=subprocess.STDOUT,timeout=timeout)
    if proc.returncode:raise RuntimeError(args[0]+' failed; original output retained')
run(['grompp','-f','stage.mdp','-c','input.gro','-t','input.cpt','-p','system.top','-o','stage.tpr'],90)
print(f'Native {a.ensemble.upper()} {a.duration_ps}ps stage prepared; CPU2; no position restraints',flush=True)
run(['mdrun','-s','stage.tpr','-deffnm','stage','-ntmpi','1','-ntomp','2','-nb','cpu'],3600)
if 'LINCS WARNING' in (out/'stage.log').read_text():raise ValueError('Constraint warning during stage')
fields=['Potential','Kinetic-En.','Temperature','Pressure']+(['Volume','Density'] if a.ensemble=='npt' else [])
run(['energy','-f','stage.edr','-o','thermodynamics.xvg','-xvg','none'],30,'\n'.join(fields)+'\n0\n')
x=np.loadtxt(out/'thermodynamics.xvg');dof=float(re.search(r'group System is ([0-9.]+)',(out/'grompp.log').read_text()).group(1))
temperature=2*x[:,2]/(dof*.00831446261815324)
if not np.isfinite(x).all() or np.max(np.abs(temperature-x[:,3]))>.001:raise ValueError('Energy/temperature or constrained DOF check failed')
with mdtraj.formats.TRRTrajectoryFile(str(out/'stage.trr')) as stream:frames=stream._read(10000,None,get_forces=True)
if not np.isfinite(frames[0]).all() or not np.isfinite(frames[-1]).all():raise ValueError('Nonfinite saved coordinate/force')
volume=np.abs(np.linalg.det(frames[3].astype(float)));density=mass*.00166053906660/volume
np.savez_compressed(out/'endpoint.npz',coords_nm=frames[0][-1],box_nm=frames[3][-1],forces_kj_mol_nm=frames[-1][-1],time_ps=frames[1][-1])
metrics={};last=x[:,0]>=x[-1,0]-10
for col,key in [(1,'potential_energy_kj_mol'),(3,'temperature_kelvin'),(4,'pressure_bar')]+([(5,'volume_nm3'),(6,'density_kg_m3')] if a.ensemble=='npt' else []):
    y=x[last,col];t=x[last,0];slope=np.polyfit(t-t[0],y,1)[0]
    metrics[key]={'mean':float(y.mean()),'standard_deviation':float(y.std()),'slope_per_ps':float(slope),
        'first5ps_mean':float(y[t<t[0]+5].mean()),'last5ps_mean':float(y[t>=t[0]+5].mean())}
report.update(status='completed',elapsed_seconds=time.monotonic()-started,mobile_dof=dof,
    temperature_maximum_recalculation_error_kelvin=float(np.max(np.abs(temperature-x[:,3]))),
    columns=['time_ps','potential_kj_mol','kinetic_kj_mol','temperature_kelvin','pressure_bar']+(['volume_nm3','density_kg_m3'] if a.ensemble=='npt' else []),
    thermodynamic_samples=x.tolist(),last10ps_metrics=metrics,box_samples_nm=frames[3].tolist(),
    box_sample_times_ps=frames[1].tolist(),calculated_density_g_ml=density.tolist(),
    finite_all_saved_forces=True,frame_count=len(frames[0]),equilibrated_acceptance='not automatic; assess retained trends before next stage')
(out/'stage-status.json').write_text(json.dumps(report,indent=2)+'\n')
print(json.dumps({'status':report['status'],'elapsed_seconds':report['elapsed_seconds'],'last10ps_metrics':metrics},indent=2))
