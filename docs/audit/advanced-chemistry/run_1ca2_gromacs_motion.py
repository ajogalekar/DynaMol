"""Bounded GROMACS motion of the unchanged native ZAFF water-box model."""
from pathlib import Path
import argparse
import hashlib
import json
import subprocess
import time
import numpy as np
import mdtraj
import parmed

BASE=Path(__file__).resolve().parent
ROOT=BASE.parents[2]
source=BASE/'prototype-1ca2-zaff-motion'
static=BASE/'prototype-1ca2-gromacs-static'
parser=argparse.ArgumentParser()
parser.add_argument('--output',type=Path,default=BASE/'prototype-1ca2-gromacs-motion')
output=parser.parse_args().output.resolve()
if not json.loads((static/'gromacs-reader-native-mechanics-null-terms-reviewed.json').read_text())['accepted']:
    raise ValueError('Exact exported nonzero parameters and native exceptions must pass before dynamics')
matrix=json.loads((static/'native-static-matrix.json').read_text())
if any(abs(r['double_precision_source_vs_export']['energy_difference_kj_mol'])>.001 or r['double_precision_source_vs_export']['maximum_force_difference_kj_mol_nm']>.01 for r in matrix['conformations']):
    raise ValueError('Independent source-vs-export force/energy comparison failed')
output.mkdir(exist_ok=False)
native=parmed.load_file(str(source/'solvated.prmtop'),xyz=str(source/'solvated.inpcrd'))
coords=np.load(source/'seed-2026-frames.npz')['positions_angstrom'][0].astype(float)
if coords.shape!=(len(native.atoms),3):raise ValueError('Minimized full-system atom inventory differs')
native.coordinates=coords
native.save(str(output/'system.top'));native.save(str(output/'input.gro'),precision=7)
readback=parmed.load_file(str(output/'system.top'))
if [(a.name,a.residue.name,a.atomic_number) for a in readback.atoms] != [(a.name,a.residue.name,a.atomic_number) for a in native.atoms]:
    raise ValueError('Exported full-system atom identity/order differs')
zn=next(a for a in native.atoms if a.residue.name=='ZN6')
bonds=[{'atoms':[b.atom1.idx,b.atom2.idx],'native_equilibrium_angstrom':b.type.req,
        'donor':b.atom2.name if b.atom1.idx==zn.idx else b.atom1.name} for b in native.bonds if zn in (b.atom1,b.atom2)]
if len(bonds)!=4:raise ValueError('Four Zn–donor bonds required')
(output/'motion.mdp').write_text('''integrator = sd
nsteps = 5000
dt = 0.001
nstxout = 100
nstvout = 0
nstfout = 100
nstenergy = 100
nstlog = 100
cutoff-scheme = Verlet
nstlist = 10
verlet-buffer-tolerance = 0.005
coulombtype = PME
rcoulomb = 1.0
ewald-rtol = 0.00001
pme-order = 4
fourierspacing = 0.12
vdwtype = Cut-off
vdw-modifier = None
rvdw = 1.0
DispCorr = EnerPres
pbc = xyz
constraints = h-bonds
constraint-algorithm = lincs
lincs-order = 6
lincs-iter = 2
tc-grps = System
tau-t = 1.0
ref-t = 300
pcoupl = no
comm-mode = Linear
nstcomm = 100
gen-vel = yes
gen-temp = 300
gen-seed = 2029
ld-seed = 2030
continuation = no
''')
started=time.monotonic();exe=ROOT/'.gromacs/bin/gmx'
def run(args,timeout):
    with (output/(args[0]+'.log')).open('w') as stream:
        proc=subprocess.run([str(exe)]+args,cwd=output,stdout=stream,stderr=subprocess.STDOUT,timeout=timeout)
    if proc.returncode:raise RuntimeError(args[0]+' failed; retained native log')
run(['grompp','-f','motion.mdp','-c','input.gro','-p','system.top','-o','motion.tpr'],90)
print('Native GROMACS 34,414-atom input passed; starting unrestrained 5ps',flush=True)
run(['mdrun','-s','motion.tpr','-deffnm','motion','-ntmpi','1','-ntomp','2','-nb','cpu'],1800)
if 'LINCS WARNING' in (output/'motion.log').read_text():raise ValueError('Native LINCS warnings during motion')
with mdtraj.formats.TRRTrajectoryFile(str(output/'motion.trr')) as stream:
    arrays=stream._read(10000,None,get_forces=True)
xyz=arrays[0].astype(float)*10;times=arrays[1].astype(float);boxes=arrays[3].astype(float)*10;forces=arrays[-1]
if len(xyz)!=51 or abs(times[-1]-5)>.0001 or not np.isfinite(xyz).all() or not np.isfinite(forces).all():
    raise ValueError('Expected 51 finite native frames through 5ps')
ca=[a.idx for a in native.atoms if a.name=='CA' and a.residue.idx<256]
metrics=[]
for frame,t,box,force in zip(xyz,times,boxes,forces,strict=True):
    inv=np.linalg.inv(box)
    def minimum_image(delta):return delta-np.rint(delta@inv)@box
    distances=[float(np.linalg.norm(minimum_image(frame[r['atoms'][0]]-frame[r['atoms'][1]]))) for r in bonds]
    unwrapped=coords[ca]+minimum_image(frame[ca]-coords[ca])
    p=unwrapped-unwrapped.mean(0);q=coords[ca]-coords[ca].mean(0)
    u,_,v=np.linalg.svd(p.T@q);sign=np.eye(3);sign[-1,-1]=np.linalg.det(u@v)
    rmsd=float(np.sqrt(np.mean(np.sum((p@u@sign@v-q)**2,axis=1))))
    metrics.append({'time_ps':float(t),'Zn_distances_angstrom':distances,'protein_CA_aligned_rmsd_angstrom':rmsd,
        'maximum_force_kj_mol_nm':float(np.linalg.norm(force,axis=1).max())})
d=np.array([r['Zn_distances_angstrom'] for r in metrics]);bounds=[[1.6,2.8] if r['donor']=='O' else [1.6,2.6] for r in bonds]
accepted=all(np.all((d[:,i]>=lo)&(d[:,i]<=hi)) for i,(lo,hi) in enumerate(bounds))
report={'status':'completed','short_motion_checks_passed':bool(accepted),'model':'Unmodified native ff14SB/ZAFF6/TIP3P; retained Zn/water/protein model',
    'engine':'GROMACS2025.4 mixed precision, CPU2','duration_ps':5,'frame_count':len(xyz),'atom_count':len(native.atoms),
    'restraints':False,'constraints':'H bonds via LINCS; ordinary waters via SETTLE; metal-donor heavy bonds dynamic',
    'initial_state':'The same complete minimized structure used for the two OpenMM seeds; new native velocities/noise',
    'periodic_coordinates':'Zn bonds use minimum images; Cα coordinates unwrapped relative to the shared initial structure before rigid alignment',
    'native_site_bonds':bonds,'screening_bounds_angstrom':bounds,'screening_meaning':'Broad geometry alarms only, not model accuracy or thermodynamic convergence',
    'Zn_min_angstrom':d.min(0).tolist(),'Zn_max_angstrom':d.max(0).tolist(),'Zn_mean_angstrom':d.mean(0).tolist(),
    'maximum_CA_rmsd_angstrom':max(r['protein_CA_aligned_rmsd_angstrom'] for r in metrics),'finite_every_saved_force':True,
    'elapsed_seconds':time.monotonic()-started,'metrics':metrics,'hashes':{n:hashlib.sha256((output/n).read_bytes()).hexdigest() for n in ['system.top','input.gro','motion.mdp','motion.tpr']}}
(output/'motion-report.json').write_text(json.dumps(report,indent=2)+'\n')
print(json.dumps({k:report[k] for k in ['status','short_motion_checks_passed','duration_ps','maximum_CA_rmsd_angstrom','elapsed_seconds']},indent=2))
