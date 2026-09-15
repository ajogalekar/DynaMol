"""Run the reviewed 1MNC research request once, with an independent deadline.

This controller reports numerical completion only. The separately frozen
endpoint policy must be assessed before any further model or Hessian work.
"""
from __future__ import annotations
import fcntl
import hashlib
import json
import os
from pathlib import Path
import shlex
import shutil
import signal
import subprocess
import time

ROOT=Path(__file__).resolve().parents[3]
OUT=ROOT/'build/advanced-chemistry/metal/qm-1mnc-anchored-optimization-v1'
RUNTIME=ROOT/'.tools/qm-parallel/bin/python'
INPUT_SHA='3b6fc0c51200556bd78bf884760b77c7feeeb079d614ceee253480e9887a5cc5'
POLICY_SHA='61c06fa5679fbaafd7db06afb225b57b021985b7976fa3e7aca8a94da3ebd93f'
WORKER_SHA='ddb67a1f8e9be6ed05c3c95e9864e20b79f0233604e27926e30e5ffc09295bc1'
GIB=1024**3


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def write(path,value):
    temporary=path.with_suffix(path.suffix+'.tmp')
    temporary.write_text(json.dumps(value,indent=2)+'\n')
    temporary.replace(path)


def verify_manifest(path,base):
    manifest=json.loads(path.read_text())
    for name,record in manifest.items():
        target=base/name
        if target.stat().st_size!=record['bytes'] or sha(target)!=record['sha256']:
            raise ValueError('Changed or unavailable pinned artifact: '+str(target))
    return len(manifest)


def active_workers():
    found=[]
    for line in subprocess.check_output(['ps','-axo','pid=,command='],text=True).splitlines():
        fields=line.strip().split(None,1)
        if len(fields)!=2:continue
        words=shlex.split(fields[1])
        if not words or 'python' not in Path(words[0]).name.lower():continue
        if not any(Path(word).name in ('qm_worker.py','covalent_reference_worker.py') for word in words[1:]):continue
        request_path=Path(words[-2])
        try:density_fit=bool(json.loads(request_path.read_text()).get('density_fit',False))
        except (OSError,ValueError):density_fit=True
        found.append({'pid':int(fields[0]),'command':fields[1],'density_fit':density_fit})
    return found


def terminate(child):
    if child is None or child.poll() is not None:return
    os.killpg(child.pid,signal.SIGTERM)
    try:child.wait(timeout=10)
    except subprocess.TimeoutExpired:
        os.killpg(child.pid,signal.SIGKILL)
        child.wait(timeout=10)


def run():
    lock=(OUT/'controller.lock').open('a')
    fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
    if (OUT/'controller-progress.json').exists() or (OUT/'optimization').exists():
        raise ValueError('This reviewed request has already been claimed; inspect retained evidence instead of repeating it')
    started=time.time()
    state={'status':'preflight','controller_pid':os.getpid(),'started_unix':started,
        'controller_sha256':sha(__file__),'scope':'Constrained-cluster numerical optimization only; not a physical model or full-system stationary point.',
        'model_accuracy_validated':False,'full_preparation_ready':False,
        'native_wall_seconds':21600,'outer_wall_seconds':21720,
        'resource_wait_deadline_unix':started+14400,
        'root_review':'Exact three-carbon hypothesis and frozen diagnostic policy reviewed before launch; continuing user-authorized local research.'}
    child=None
    def update(**values):
        state.update(values);state['elapsed_seconds']=time.time()-started
        write(OUT/'controller-progress.json',state)
    def interrupted(signum,frame):
        raise InterruptedError('Controller received signal '+str(signum))
    signal.signal(signal.SIGTERM,interrupted)
    signal.signal(signal.SIGINT,interrupted)
    update()
    try:
        for name,expected in [('input.json',INPUT_SHA),('optimization-policy.json',POLICY_SHA),('qm_worker.py',WORKER_SHA)]:
            if sha(OUT/name)!=expected:raise ValueError('Reviewed source changed: '+name)
        artifact_count=verify_manifest(OUT/'artifact-manifest.json',OUT)
        source_count=verify_manifest(OUT/'source-manifest.json',ROOT)
        (OUT/'controller-source.py').write_bytes(Path(__file__).read_bytes())
        update(status='waiting_for_resources',input_sha256=INPUT_SHA,worker_sha256=WORKER_SHA,
            policy_sha256=POLICY_SHA,artifact_manifest_sha256=sha(OUT/'artifact-manifest.json'),
            source_manifest_sha256=sha(OUT/'source-manifest.json'),
            verified_artifact_count=artifact_count,verified_source_count=source_count,
            runtime=str(RUNTIME),runtime_executable_sha256=sha(RUNTIME))
        with (ROOT/'build/advanced-chemistry/QM-LAUNCH.lock').open('a') as launch_lock:
            while child is None:
                if time.time()>state['resource_wait_deadline_unix']:
                    raise TimeoutError('Resource wait expired before any quantum launch')
                try:fcntl.flock(launch_lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
                except BlockingIOError:
                    update(resource_hold='Shared quantum launch lock is held')
                    time.sleep(15);continue
                try:
                    workers=active_workers()
                    if any(str(OUT) in row['command'] for row in workers):
                        raise ValueError('A worker already owns this exact request')
                    count_df=sum(row['density_fit'] for row in workers)
                    free=shutil.disk_usage(OUT).free
                    minimum=(35 if count_df else 22)*GIB
                    update(active_qm_workers=workers,free_disk_bytes=free,minimum_start_free_bytes=minimum)
                    if len(workers)<3 and count_df<2 and free>=minimum:
                        with (OUT/'controller-native.log').open('x') as log:
                            child=subprocess.Popen([str(RUNTIME),str(OUT/'qm_worker.py'),str(OUT/'input.json'),str(OUT/'optimization')],
                                cwd=ROOT,stdout=log,stderr=subprocess.STDOUT,start_new_session=True)
                        update(status='optimization_running',worker_pid=child.pid,
                            qm_started_unix=time.time(),outer_deadline_unix=time.time()+21720)
                finally:fcntl.flock(launch_lock,fcntl.LOCK_UN)
                if child is None:time.sleep(15)
        while child.poll() is None:
            if time.time()>state['outer_deadline_unix']:
                raise TimeoutError('Independent six-hour optimization deadline reached')
            free=shutil.disk_usage(OUT).free
            update(free_disk_bytes=free)
            if free<8*GIB:raise RuntimeError('Stopped own child to preserve the 8 GiB disk reserve')
            time.sleep(15)
        result_path=OUT/'optimization/result.json'
        result=json.loads(result_path.read_text())
        if child.returncode or not result.get('accepted') or not result.get('optimization',{}).get('converged'):
            raise ValueError('Optimization did not converge successfully; original outputs retained')
        if result['input_sha256']!=INPUT_SHA or result['worker_sha256']!=WORKER_SHA:
            raise ValueError('Numerical result changed pinned input or worker')
        if sha(OUT/'optimization/arrays.npz')!=result['arrays_sha256']:
            raise ValueError('Numerical output hash mismatch')
        update(status='numerically_complete_pending_endpoint_review',result_sha256=sha(result_path),
            endpoint_policy_assessed=False,chemical_accuracy_validated=False)
    except Exception as error:
        terminate(child)
        update(status='failed',error_type=type(error).__name__,error=str(error),
            worker_returncode=child.returncode if child else None)
    write(OUT/'controller-result.json',state)
    return state


if __name__=='__main__':
    final=run()
    print(json.dumps(final,indent=2))
    raise SystemExit(1 if final['status']=='failed' else 0)
