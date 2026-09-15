"""Durable, single-claim research DFT continuation. No launch on import/preflight.

The pinned materializer owns chemistry checks. This controller binds the actual
parent and completed RESP chain, verifies files, and owns external resource and
process lifetime limits. It never substitutes charges or accepts a force field.
"""
from __future__ import annotations
import argparse
from contextlib import contextmanager
import fcntl
import hashlib
import json
import math
import os
from pathlib import Path
import shutil
import signal
import subprocess
import sys
import time

GIB=1024**3
class Cancelled(Exception): pass


def sha(path):
    with Path(path).open('rb') as stream:
        h=hashlib.sha256()
        for chunk in iter(lambda:stream.read(1024**2),b''):h.update(chunk)
    return h.hexdigest()


def read(path):return json.loads(Path(path).read_text())
def write(path,value,exclusive=False):
    if exclusive:
        with Path(path).open('x') as stream:json.dump(value,stream,indent=2);stream.write('\n')
    else:
        temp=Path(path).with_suffix('.json.tmp');temp.write_text(json.dumps(value,indent=2)+'\n');temp.replace(path)

def check_file(path,expected,observed=None):
    actual=sha(path)
    if actual!=expected:raise ValueError(f'Pinned file changed: {path}')
    if observed is not None:observed[str(Path(path).resolve())]=actual
    return actual


def verify_static(config):
    """Every explicit manifest file hash is checked, including native libraries."""
    checked={}
    for item in config['pins']:check_file(item['path'],item['sha256'],checked)
    plan=read(config['plan_path']);provider=plan['provider']
    required=[config['plan_path'],config['materializer'],config['materializer_runtime'],
              provider['runtime'],config['source_manifest'],config['shared_runtime_manifest']]
    if any(str(Path(p).resolve()) not in checked for p in required):
        raise ValueError('Controller configuration omits a mandatory executable/manifest pin')
    for path_key,hash_key in [('worker','worker_sha256'),('runtime_manifest','runtime_manifest_sha256')]:
        check_file(provider[path_key],provider[hash_key],checked)
    check_file(plan['qualification']['result_path'],plan['qualification']['result_sha256'],checked)
    check_file(plan['parent_geometry_validator']['path'],plan['parent_geometry_validator']['sha256'],checked)
    check_file(plan['parent']['input_path'],plan['parent']['input_sha256'],checked)
    for name,digest in plan['parent']['cap_files_sha256'].items():check_file(Path(plan['parent']['cap_dir'])/name,digest,checked)
    runtime=read(provider['runtime_manifest'])
    if Path(runtime['runtime']).resolve()!=Path(provider['runtime']).resolve():raise ValueError('Runtime manifest names a different executable')
    for item in runtime['files']:check_file(item['path'],item['sha256'],checked)
    sources_path=Path(config['source_manifest']);sources=read(sources_path)
    for item in sources['sources']:check_file(sources_path.parent/item['filename'],item['sha256'],checked)
    shared_path=Path(config['shared_runtime_manifest']);shared=read(shared_path)
    for relative,digest in shared['native_libraries'].items():check_file(shared_path.parent/relative,digest,checked)
    check_file(shared['openmp_library']['path'],shared['openmp_library']['sha256'],checked)
    qualification=read(plan['qualification']['result_path'])
    if qualification['status']!='passed_numerical_qualification' or qualification['worker_sha256']!=provider['worker_sha256']:
        raise ValueError('Provider does not have matching passed numerical qualification')
    for fixture in qualification['tests']['small_molecules']:
        rp=Path(fixture['result_path']);check_file(rp,fixture['result_sha256'],checked);r=read(rp)
        check_file(rp.parent/r['arrays_file'],r['arrays_sha256'],checked)
    settings=plan['request_settings']
    if settings['method']!='B3LYP-D3BJ/Psi4-DZVP' or settings['threads']!=2 or settings['max_memory_mb']!=8000:
        raise ValueError('Unexpected reviewed parent method/resource settings')
    if not 0<settings['max_wall_seconds']<=14400 or not 0<config['wait_seconds']<=28800:
        raise ValueError('Controller budget exceeds reviewed limit')
    if not 0<config['poll_seconds']<=30:raise ValueError('Invalid monitoring interval')
    return plan,checked


def dependencies(config,plan,observed):
    """Pending is allowed, but changed observed artifacts or wrong lineage fail."""
    for path,digest in list(observed.items()):check_file(path,digest)
    parent=plan['parent'];optimizer=Path(parent['output_dir']);respdir=Path(config['resp_dir'])
    check_file(parent['input_path'],parent['input_sha256'])
    if (optimizer/'input.json').exists():check_file(optimizer/'input.json',parent['input_sha256'])
    progress_path=respdir/'progress.json'
    if not progress_path.exists():raise ValueError('The existing RESP continuation is required')
    progress=read(progress_path)
    for r in [progress]+([read(respdir/'result.json')] if (respdir/'result.json').exists() else []):
        if r.get('optimizer_input_sha256')!=parent['input_sha256'] or Path(r.get('optimizer','')).resolve()!=optimizer.resolve():
            raise ValueError('RESP continuation refers to a different parent')
        if r.get('script_sha256')!=config['resp_controller_sha256']:raise ValueError('Existing RESP controller identity changed')
        if r.get('cap_graph_sha256')!=parent['cap_files_sha256']['capped-adduct.json'] or r.get('cap_sdf_sha256')!=parent['cap_files_sha256']['capped-adduct.sdf']:
            raise ValueError('RESP chemical graph changed')
        if r.get('status') not in ('running','research_candidate'):raise ValueError('RESP continuation failed or has an unreviewed state')
    result_path=optimizer/'result.json'
    if not result_path.exists():
        for path in [optimizer/'progress.json',optimizer.parent/'progress.json']:
            if path.exists() and read(path).get('status') in ('failed','timed_out','cancelled'):raise ValueError('Parent optimization failed')
        return 'waiting_for_parent'
    result=read(result_path)
    if not result.get('accepted') or not result.get('optimization',{}).get('converged'):
        raise ValueError('Actual parent did not converge')
    if result.get('input_sha256')!=parent['input_sha256']:raise ValueError('Parent result input binding differs')
    observed[str(result_path.resolve())]=sha(result_path)
    check_file(optimizer/'arrays.npz',result['arrays_sha256'],observed)
    resp_path=respdir/'result.json'
    if not resp_path.exists():return 'waiting_for_resp'
    resp=read(resp_path)
    if resp.get('status')!='research_candidate' or progress.get('status')!='research_candidate':
        raise ValueError('A completed RESP research candidate is required')
    stage=resp['stages']['parent']
    if stage['result_sha256']!=sha(result_path) or stage['arrays_sha256']!=result['arrays_sha256']:
        raise ValueError('RESP result belongs to a different optimized geometry')
    esp_stage=resp['stages']['exact_esp_result'];esp_dir=Path(esp_stage['path'])
    if esp_dir.resolve()!=(respdir/'exact-esp').resolve():raise ValueError('Unexpected existing ESP output location')
    check_file(esp_dir/'result.json',esp_stage['result_sha256'],observed)
    esp=read(esp_dir/'result.json');esp_input=read(respdir/'exact-esp-input.json')
    if not esp.get('accepted') or esp['input_sha256']!=sha(respdir/'exact-esp-input.json'):
        raise ValueError('ESP is not bound to its input')
    if esp_input['initial_checkpoint']['result_sha256']!=sha(result_path) or Path(esp_input['initial_checkpoint']['result_path']).resolve()!=result_path.resolve():
        raise ValueError('ESP checkpoint differs from the accepted parent')
    check_file(esp_dir/'arrays.npz',esp_stage['arrays_sha256'],observed)
    if esp['arrays_sha256']!=esp_stage['arrays_sha256']:raise ValueError('ESP array hash mismatch')
    fit=resp['stages']['resp']['fit']
    if not fit.get('constraint_checks_passed'):raise ValueError('RESP constraint checks did not pass')
    check_file(respdir/'hf-resp.esp',fit['esp_sha256'],observed)
    check_file(fit['resp_executable'],fit['resp_executable_sha256'],observed)
    check_file(config['resp_fitter'],fit['fitter_source_sha256'],observed)
    if resp['stages']['resp']['fitter_sha256']!=fit['fitter_source_sha256']:raise ValueError('RESP fitter lineage differs')
    observed[str(resp_path.resolve())]=sha(resp_path)
    observed[str((respdir/'exact-esp-input.json').resolve())]=sha(respdir/'exact-esp-input.json')
    return 'ready'


def resource_snapshot(root):
    active=[]
    for row in subprocess.check_output(['ps','-axo','pid=,command='],text=True).splitlines():
        parts=row.strip().split(None,1)
        if len(parts)!=2 or int(parts[0])==os.getpid():continue
        command=parts[1]
        if 'python' in Path(command.split()[0]).name and ('qm_worker.py' in command or 'covalent_reference_worker.py' in command):
            active.append({'pid':int(parts[0]),'command':command})
    available=True;path=Path(root)/'build/advanced-chemistry/QM-LAUNCH.lock'
    with path.open('a') as lock:
        try:fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
        except BlockingIOError:available=False
        # Release before provider launch: the provider owns this lock throughout.
    return {'active_qm':active,'free_disk_bytes':shutil.disk_usage(root).free,'shared_lock_available':available}


def terminate_group(child,grace=10):
    def alive_group():
        rows=subprocess.check_output(['ps','-axo','pgid=,stat='],text=True).splitlines()
        return any(len(parts:=row.split())==2 and int(parts[0])==child.pid and not parts[1].startswith('Z') for row in rows)
    try:os.killpg(child.pid,signal.SIGTERM)
    except ProcessLookupError:return
    except PermissionError:
        if child.poll() is not None and not alive_group():return
        raise
    deadline=time.monotonic()+grace
    while time.monotonic()<deadline:
        child.poll()
        # macOS can return EPERM for signal-0 during group teardown. Inspect
        # executable group members instead, without hiding a live child.
        if not alive_group():return
        time.sleep(.05)
    try:os.killpg(child.pid,signal.SIGKILL)
    except ProcessLookupError:pass
    child.wait()


@contextmanager
def uninterrupted_cleanup():
    """Ignore repeat cancellations while the owned process group is reaped.

    Block both signals atomically while installing/restoring handlers; changing
    just one handler at a time would leave a second-signal interruption window.
    This is cleanup only: the ordinary controller cancellation handler remains
    active during work and monitoring.
    """
    signals={signal.SIGINT,signal.SIGTERM}
    old_mask=signal.pthread_sigmask(signal.SIG_BLOCK,signals)
    handlers={sig:signal.signal(sig,signal.SIG_IGN) for sig in signals}
    signal.pthread_sigmask(signal.SIG_SETMASK,old_mask)
    try:
        yield
    finally:
        signal.pthread_sigmask(signal.SIG_BLOCK,signals)
        for sig,handler in handlers.items():signal.signal(sig,handler)
        signal.pthread_sigmask(signal.SIG_SETMASK,old_mask)


def run_child(command,output,timeout,poll,root,on_start,on_tick,reserve=8*GIB):
    """A separate process group contains the native worker and descendants."""
    started=time.monotonic();child=None
    with Path(output).open('x') as log:
        try:
            env=dict(os.environ,OMP_NUM_THREADS='2',OPENBLAS_NUM_THREADS='2',MKL_NUM_THREADS='2',
                     VECLIB_MAXIMUM_THREADS='2',NUMEXPR_NUM_THREADS='2')
            child=subprocess.Popen(command,stdout=log,stderr=subprocess.STDOUT,start_new_session=True,env=env)
            on_start(child.pid)
            while child.poll() is None:
                if time.monotonic()-started>=timeout:raise TimeoutError('Outer native process-group deadline expired')
                if shutil.disk_usage(root).free<reserve:raise RuntimeError('Running free disk reserve fell below 8 GiB')
                on_tick()
                time.sleep(min(poll,max(.01,timeout-(time.monotonic()-started))))
            if child.returncode:raise RuntimeError(f'Child exited with return code {child.returncode}; output retained')
            return {'pid':child.pid,'returncode':child.returncode,'elapsed_seconds':time.monotonic()-started}
        finally:
            if child is not None:
                with uninterrupted_cleanup():
                    terminate_group(child)
                    child.wait()


def validate_reference(plan,request_path,calculation,expected_request_sha256):
    # Use the digest captured from the materializer's first completed bytes,
    # never a self-consistent replacement request encountered at completion.
    request_bytes=Path(request_path).read_bytes()
    if hashlib.sha256(request_bytes).hexdigest()!=expected_request_sha256:
        raise ValueError('Materialized request changed after its original snapshot')
    request=json.loads(request_bytes);result=read(calculation/'result.json')
    if result.get('accepted') is not True or result.get('status')!='numerically_complete_research_reference':raise ValueError('Reference provider did not succeed')
    if result['source']['worker_sha256']!=plan['provider']['worker_sha256'] or result['source']['input_sha256']!=expected_request_sha256:raise ValueError('Reference source/input identity mismatch')
    check_file(calculation/'input.json',expected_request_sha256)
    check_file(calculation/'arrays.npz',result['arrays_sha256'])
    for field in ('atom_ids','elements','charge','spin'):
        if result[field]!=request[field]:raise ValueError('Reference atom/state identity differs')
    if result['method']['name']!=request['method'] or result['method']['density_fitting'] or result['method']['libxc_id']!=402:raise ValueError('Reference method differs')
    if not all(math.isfinite(result[k]) for k in ['energy_hartree','electronic_energy_hartree','dispersion_energy_hartree']):raise ValueError('Nonfinite reference energy')
    return {'result_path':str(calculation/'result.json'),'result_sha256':sha(calculation/'result.json'),'arrays_sha256':result['arrays_sha256'],
            'energy_hartree':result['energy_hartree'],'method':result['method'],'physical_acceptance':False}


@contextmanager
def claim(config,plan,output):
    identity={'parent_input_sha256':plan['parent']['input_sha256'],'parent_output_dir':str(Path(plan['parent']['output_dir']).resolve()),'method':plan['request_settings']['method']}
    key=hashlib.sha256(json.dumps(identity,sort_keys=True).encode()).hexdigest()
    folder=Path(config['claims_dir']);folder.mkdir(parents=True,exist_ok=True)
    with (folder/f'{key}.lock').open('a') as lock:
        fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
        # Claim survives failure/interruption/completion; a restart requires review.
        write(folder/f'{key}.json',{'identity':identity,'pid':os.getpid(),'output':str(output),'created_unix':time.time()},exclusive=True)
        yield key


def execute(config_path,output,preflight=False):
    config=read(config_path);plan,pins=verify_static(config);observed={}
    if preflight:return {'launched':False,'stage':dependencies(config,plan,observed),'verified_files':len(pins),'resource':resource_snapshot(config['root'])}
    output=Path(output).resolve()
    if output.exists():raise FileExistsError('Use a fresh output directory; no prior calculation is overwritten')
    with claim(config,plan,output) as key:
        output.mkdir(parents=True,exist_ok=False)
        state={'schema_version':1,'status':'running','stage':'checking_dependencies','controller_pid':os.getpid(),'child_pid':None,
            'started_unix':time.time(),'wait_deadline_unix':time.time()+config['wait_seconds'],'claim':key,'physical_acceptance':False,
            'scope':'B3LYP-D3BJ/DZVP single point at accepted DF-RHF/6-31G* capped parent; existing HF ESP/RESP retained. No torsion or force-field validation.',
            'config_sha256':sha(config_path),'script_sha256':sha(__file__),'verified_static_files':pins,'observed_dependencies':observed}
        write(output/'controller-config.json',config,exclusive=True);write(output/'controller-source.py.json',{'path':str(Path(__file__).resolve()),'sha256':sha(__file__)},exclusive=True)
        shutil.copy2(__file__,output/'controller-source.py')
        handlers={}
        def cancel(signum,frame):raise Cancelled(f'Controller cancellation signal {signum}')
        for sig in (signal.SIGINT,signal.SIGTERM):handlers[sig]=signal.signal(sig,cancel)
        def update(**kw):
            state.update(kw);state['updated_unix']=time.time();write(output/'progress.json',state)
        def record_child(pid,command,stage):
            state.setdefault('children',[]).append({'pid':pid,'process_group_id':pid,
                'stage':stage,'command':command,'started_unix':time.time()})
            update(child_pid=pid)
        def verify():
            check_file(config_path,state['config_sha256']);verify_static(config)
            return dependencies(config,plan,observed)
        try:
            while True:
                if time.time()>=state['wait_deadline_unix']:raise TimeoutError('Bounded eight-hour dependency/resource wait expired')
                stage=verify();resource=resource_snapshot(config['root'])
                update(stage=stage,resources=resource)
                if time.time()>=state['wait_deadline_unix']:raise TimeoutError('Bounded eight-hour dependency/resource wait expired')
                if stage=='ready' and len(resource['active_qm'])<3 and resource['free_disk_bytes']>=25*GIB and resource['shared_lock_available']:break
                time.sleep(min(config['poll_seconds'],max(.01,state['wait_deadline_unix']-time.time())))
            update(stage='materializing_accepted_parent')
            request_path=output/'reference-input.json'
            materializer_command=[config['materializer_runtime'],config['materializer'],config['plan_path'],str(request_path)]
            materialized=run_child(materializer_command,
                output/'materializer.log',120,config['poll_seconds'],config['root'],
                lambda pid:record_child(pid,materializer_command,'materialization'),lambda:None)
            request_bytes=request_path.read_bytes()
            request_sha256=hashlib.sha256(request_bytes).hexdigest()
            request=json.loads(request_bytes)
            update(request_sha256=request_sha256)
            if request['provenance']['parent_result_sha256']!=observed[str((Path(plan['parent']['output_dir'])/'result.json').resolve())]:raise ValueError('Materialized parent changed')
            if request['provenance']['qualified_worker_sha256']!=plan['provider']['worker_sha256']:raise ValueError('Materialized provider changed')
            if verify()!='ready':raise ValueError('Dependencies changed after materialization')
            resource=resource_snapshot(config['root'])
            if len(resource['active_qm'])>=3 or resource['free_disk_bytes']<25*GIB or not resource['shared_lock_available']:raise RuntimeError('Resources changed before provider launch; no retry or overwritten output')
            update(stage='reference_single_point',materializer=materialized,outer_timeout_seconds=plan['request_settings']['max_wall_seconds']+120)
            # Intentionally no shared QM lock here. The qualified provider owns it.
            native_command=[plan['provider']['runtime'],plan['provider']['worker'],str(request_path),str(output/'calculation')]
            check_file(request_path,request_sha256)
            run=run_child(native_command,
                output/'native-controller.log',plan['request_settings']['max_wall_seconds']+120,config['poll_seconds'],config['root'],
                lambda pid:record_child(pid,native_command,'reference_single_point'),lambda:update(child_monitor_unix=time.time()))
            verify();check_file(request_path,request_sha256)
            reference=validate_reference(plan,request_path,output/'calculation',request_sha256)
            update(status='numerically_complete_research_reference',stage='awaiting_varied_geometries_and_physical_validation',child_pid=None,run=run,reference=reference)
        except (Exception,KeyboardInterrupt) as error:
            update(status='cancelled' if isinstance(error,(Cancelled,KeyboardInterrupt)) else 'failed',stage='stopped',child_pid=None,error_type=type(error).__name__,error=str(error))
        finally:
            for sig,handler in handlers.items():signal.signal(sig,handler)
            state['elapsed_seconds']=time.time()-state['started_unix']
            write(output/'progress.json',state)
            write(output/'result.json',state,exclusive=True)
        return state


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('config',type=Path);parser.add_argument('--output',type=Path)
    parser.add_argument('--preflight',action='store_true')
    args=parser.parse_args()
    if not args.preflight and args.output is None:parser.error('--output is required unless --preflight')
    result=execute(args.config,args.output,args.preflight);print(json.dumps(result,indent=2))
    if result.get('status') in ('failed','cancelled'):raise SystemExit(1)
