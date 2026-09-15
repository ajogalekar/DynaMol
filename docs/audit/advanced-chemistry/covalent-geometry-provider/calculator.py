"""Bounded research GFN2-xTB energy/gradient interface; no optimizer or charges.

Call evaluate(request) only inside execution_guard(root), with an external
process-group timeout. Results are in Hartree and Hartree/Bohr, not forces.
"""
from contextlib import contextmanager
import ctypes
import fcntl
import hashlib
import importlib.metadata
import json
import os
from pathlib import Path
import shutil
import shlex
import subprocess
import sys

HERE=Path(__file__).resolve().parent
RUNTIME=Path('/Users/ashujo/.cache/dynamol-runtimes/covalent-geometry-v1')
ELEMENTS={'H':1,'C':6,'N':7,'O':8,'F':9,'S':16}
_ACTIVE=None


def is_quantum_worker(command):
    """Classify the executed Python script, not a test/supervisor substring."""
    try:parts=shlex.split(command)
    except ValueError:return False
    if not parts or 'python' not in Path(parts[0]).name:return False
    index=1
    while index<len(parts) and parts[index].startswith('-'):
        flag=parts[index]
        if flag in ('-c','-m','-'):return False
        index+=2 if flag in ('-W','-X') else 1
    if index>=len(parts):return False
    script=Path(parts[index]).name
    if script in ('qm_worker.py','covalent_reference_worker.py','qualify_calculator.py'):return True
    return script=='geometric_adapter.py' and '--worker' in parts[index+1:]


def sha(path):
    h=hashlib.sha256()
    with Path(path).open('rb') as stream:
        for block in iter(lambda:stream.read(1024**2),b''):h.update(block)
    return h.hexdigest()


def verify_runtime():
    manifest=json.loads((HERE/'runtime-manifest.json').read_text())
    if Path(sys.prefix).resolve()!=RUNTIME.resolve():raise ValueError('Use the isolated geometry runtime')
    for package,version in manifest['packages'].items():
        if importlib.metadata.version(package)!=version:raise ValueError(f'Runtime package changed: {package}')
    for item in manifest['files']:
        if sha(item['path'])!=item['sha256']:raise ValueError('Runtime file changed: '+item['path'])
    libs=RUNTIME/'lib/python3.12/site-packages/tblite/.dylibs'
    omp=ctypes.CDLL(str(libs/'libgomp.1.dylib'));omp.omp_set_num_threads(1);omp.omp_get_max_threads.restype=ctypes.c_int
    blas=ctypes.CDLL(str(libs/'libopenblasp-r0.3.33.dylib'));blas.openblas_set_num_threads(1);blas.openblas_get_num_threads.restype=ctypes.c_int
    observed={'tblite_openmp_max_threads':omp.omp_get_max_threads(),'tblite_openblas_threads':blas.openblas_get_num_threads()}
    if set(observed.values())!={1}:raise ValueError('Native thread limit was not applied')
    return {'manifest_sha256':sha(HERE/'runtime-manifest.json'),'files_checked':len(manifest['files']),
        'packages':manifest['packages'],'native_threads':observed,'runtime':str(RUNTIME)}


@contextmanager
def execution_guard(root):
    global _ACTIVE
    if _ACTIVE is not None:raise ValueError('Nested geometry execution guards are not supported')
    root=Path(root).resolve();lock_path=root/'build/advanced-chemistry/QM-LAUNCH.lock'
    with lock_path.open('a') as lock:
        fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
        active=[]
        for row in subprocess.check_output(['ps','-axo','pid=,command='],text=True).splitlines():
            parts=row.strip().split(None,1)
            if len(parts)==2 and int(parts[0])!=os.getpid() and is_quantum_worker(parts[1]):
                active.append({'pid':int(parts[0]),'command':parts[1]})
        if len(active)>=3:raise ValueError('No free quantum execution slot')
        if shutil.disk_usage(root).free<2*1024**3:raise ValueError('Small geometry qualification requires 2 GiB free')
        for key in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','VECLIB_MAXIMUM_THREADS','NUMEXPR_NUM_THREADS']:os.environ[key]='1'
        runtime=verify_runtime()
        _ACTIVE={**runtime,'pid':os.getpid(),'active_other_quantum_workers':active,'lock':str(lock_path)}
        try:yield dict(_ACTIVE)
        finally:_ACTIVE=None


def validate(request):
    import numpy as np
    allowed={'schema_version','method','atom_ids','elements','charge','spin','coords_bohr','accuracy','max_iter','temperature_hartree'}
    if set(request)-allowed:raise ValueError('Unknown geometry-calculator setting')
    if request.get('schema_version')!=1 or request.get('method')!='GFN2-xTB':raise ValueError('Explicit schema and GFN2-xTB method required')
    ids=request.get('atom_ids');elements=request.get('elements')
    if not isinstance(ids,list) or not 1<=len(ids)<=32 or any(not isinstance(i,str) or not i for i in ids) or len(ids)!=len(set(ids)):
        raise ValueError('One to 32 unique atom IDs required for this qualification')
    if not isinstance(elements,list) or len(elements)!=len(ids) or any(x not in ELEMENTS for x in elements):raise ValueError('Explicit H/C/N/O/F/S elements required')
    if type(request.get('charge')) is not int or type(request.get('spin')) is not int or request['spin']!=0:raise ValueError('Explicit integer charge and closed-shell spin 0 required')
    electrons=sum(ELEMENTS[x] for x in elements)-request['charge']
    if electrons<=0 or electrons%2:raise ValueError('Closed-shell state requires a positive even electron count')
    xyz=np.asarray(request.get('coords_bohr'),dtype=float)
    if xyz.shape!=(len(ids),3) or not np.isfinite(xyz).all():raise ValueError('Finite Nx3 coordinates in Bohr required')
    if len(xyz)>1:
        distances=np.linalg.norm(xyz[:,None]-xyz[None,:],axis=2);np.fill_diagonal(distances,np.inf)
        if distances.min()<.3:raise ValueError('Coincident or implausibly close atoms')
    accuracy=request.get('accuracy',.001);max_iter=request.get('max_iter',250);temperature=request.get('temperature_hartree',.00095)
    if type(accuracy) not in (int,float) or not 1e-5<=accuracy<=1:raise ValueError('Invalid SCC accuracy')
    if type(max_iter) is not int or not 1<=max_iter<=250:raise ValueError('Invalid SCC iteration limit')
    if type(temperature) not in (int,float) or not 0<=temperature<=.001:raise ValueError('Invalid electronic temperature (Hartree)')
    return xyz,{'accuracy':accuracy,'max_iter':max_iter,'temperature_hartree':temperature}


def evaluate(request):
    if _ACTIVE is None:raise RuntimeError('Use execution_guard(root) and an external timeout')
    import numpy as np
    from tblite.interface import Calculator
    from tblite.library import get_version
    xyz,settings=validate(request)
    calc=Calculator('GFN2-xTB',np.array([ELEMENTS[x] for x in request['elements']]),xyz,charge=request['charge'],uhf=0)
    for key,value in [('accuracy',settings['accuracy']),('max-iter',settings['max_iter']),('temperature',settings['temperature_hartree']),('verbosity',0)]:calc.set(key,value)
    native=calc.singlepoint()  # Native nonconvergence raises; no last-iterate acceptance.
    energy=float(native.get('energy'));gradient=np.asarray(native.get('gradient'),dtype=float)
    if not np.isfinite(energy) or gradient.shape!=xyz.shape or not np.isfinite(gradient).all():raise ValueError('Nonfinite/malformed native energy or gradient')
    return {'accepted':True,'status':'numerically_complete_research_geometry_reference','physical_acceptance':False,
        'method':'GFN2-xTB','native_method_scope':'Self-consistent GFN2 Hamiltonian with its native D4 dispersion; gas phase, no external restraints or penalties.',
        'atom_ids':list(request['atom_ids']),'elements':list(request['elements']),'charge':request['charge'],'spin':0,
        'coords_bohr':xyz.tolist(),'energy_hartree':energy,'gradient_hartree_per_bohr':gradient.tolist(),
        'settings':settings,'units':{'coordinates':'Bohr','energy':'Hartree','gradient':'Hartree/Bohr; positive dE/dx'},
        'native_version':list(get_version()),'runtime':dict(_ACTIVE),'calculator_sha256':sha(__file__),
        'limitations':'Numerical implementation only. No geometry/stereochemistry/chemical-state or force-field accuracy inference; no charges exported.'}
