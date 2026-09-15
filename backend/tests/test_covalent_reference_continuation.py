"""Bounded controller integration tests use tiny local Python stubs, never QM."""
import hashlib
import importlib.util
import json
import os
from pathlib import Path
import subprocess
import sys
import time
from types import SimpleNamespace

import pytest

PATH=Path(__file__).resolve().parents[2]/'docs/audit/advanced-chemistry/continue_covalent_reference.py'
spec=importlib.util.spec_from_file_location('reference_continuation_audit',PATH)
c=importlib.util.module_from_spec(spec);spec.loader.exec_module(c)

def j(path,value):path.parent.mkdir(parents=True,exist_ok=True);path.write_text(json.dumps(value,indent=2)+'\n');return path

def h(path):return hashlib.sha256(path.read_bytes()).hexdigest()

@pytest.fixture
def prepared(tmp_path,monkeypatch):
    root=tmp_path;parent=root/'parent';resp=root/'resp';cap=root/'cap';shared=root/'shared'
    input_path=j(root/'parent-input.json',{'atom_ids':['A'],'elements':['H'],'charge':0,'spin':0})
    parent.mkdir();(parent/'input.json').write_bytes(input_path.read_bytes());(parent/'arrays.npz').write_bytes(b'parent arrays')
    result=j(parent/'result.json',{'accepted':True,'optimization':{'converged':True},'input_sha256':h(input_path),'arrays_sha256':h(parent/'arrays.npz')})
    for name in ['capped-adduct.json','capped-adduct.sdf']:j(cap/name,{'fixture':name})
    fitter=j(root/'fitter.py',{'fixture':'fitter'});native=j(root/'native-resp',{'fixture':'native'})
    resp.mkdir();(resp/'hf-resp.esp').write_bytes(b'esp fixture')
    espi=j(resp/'exact-esp-input.json',{'initial_checkpoint':{'result_path':str(result),'result_sha256':h(result)}})
    espdir=resp/'exact-esp';espdir.mkdir();(espdir/'arrays.npz').write_bytes(b'esp arrays')
    espr=j(espdir/'result.json',{'accepted':True,'input_sha256':h(espi),'arrays_sha256':h(espdir/'arrays.npz')})
    fit={'constraint_checks_passed':True,'esp_sha256':h(resp/'hf-resp.esp'),'resp_executable':str(native),'resp_executable_sha256':h(native),'fitter_source_sha256':h(fitter)}
    rr={'status':'research_candidate','script_sha256':'resident-controller','optimizer_input_sha256':h(input_path),'optimizer':str(parent),
        'cap_graph_sha256':h(cap/'capped-adduct.json'),'cap_sdf_sha256':h(cap/'capped-adduct.sdf'),
        'stages':{'parent':{'result_sha256':h(result),'arrays_sha256':h(parent/'arrays.npz')},
        'exact_esp_result':{'path':str(espdir),'result_sha256':h(espr),'arrays_sha256':h(espdir/'arrays.npz')},
        'resp':{'fit':fit,'fitter_sha256':h(fitter)}}}
    j(resp/'progress.json',rr);j(resp/'result.json',rr)
    materializer=root/'materializer.py';materializer.write_text('''import json,sys,hashlib
from pathlib import Path
p=json.loads(Path(sys.argv[1]).read_text());r=p['request_settings'];r.update(atom_ids=p['atom_ids'],elements=p['elements'],coords_bohr=[[0,0,0]])
r['provenance']={'parent_result_sha256':hashlib.sha256((Path(p['parent']['output_dir'])/'result.json').read_bytes()).hexdigest(),'qualified_worker_sha256':p['provider']['worker_sha256']}
Path(sys.argv[2]).write_text(json.dumps(r,indent=2)+'\\n')
''')
    worker=root/'covalent_reference_worker.py';worker.write_text('''import json,sys,hashlib,time,signal,subprocess
from pathlib import Path
if (Path(__file__).parent/'fail-worker').exists():raise SystemExit(17)
if (Path(__file__).parent/'ignore-termination').exists():
 signal.signal(signal.SIGTERM,signal.SIG_IGN)
 ready=Path(__file__).parent/'grandchild-ready'
 p=subprocess.Popen([sys.executable,'-c','import signal,time,sys;from pathlib import Path;signal.signal(signal.SIGTERM,signal.SIG_IGN);Path(sys.argv[1]).write_text("ready");time.sleep(60)',str(ready)])
 (Path(__file__).parent/'grandchild.pid').write_text(str(p.pid))
if (Path(__file__).parent/'hang-worker').exists():time.sleep(60)
p=Path(sys.argv[1]);o=Path(sys.argv[2]);o.mkdir();r=json.loads(p.read_text())
if (Path(__file__).parent/'mutate-request-before-read').exists():
 r['coords_bohr']=[[9.,0.,0.]];p.write_text(json.dumps(r,indent=2)+'\\n')
(o/'input.json').write_bytes(p.read_bytes());(o/'arrays.npz').write_bytes(b'stub arrays')
d={k:r[k] for k in ['atom_ids','elements','charge','spin']};d.update(accepted=True,status='numerically_complete_research_reference',source={'worker_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),'input_sha256':hashlib.sha256(p.read_bytes()).hexdigest()},arrays_sha256=hashlib.sha256((o/'arrays.npz').read_bytes()).hexdigest(),energy_hartree=0.,electronic_energy_hartree=0.,dispersion_energy_hartree=0.,method={'name':r['method'],'density_fitting':False,'libxc_id':402})
(o/'result.json').write_text(json.dumps(d))
''')
    runtime=j(root/'runtime.json',{'runtime':sys.executable,'files':[]})
    source=j(root/'sources.json',{'sources':[]})
    shared.mkdir();lib=j(shared/'lib',{'fixture':'lib'});omp=j(shared/'omp',{'fixture':'omp'})
    sharedmanifest=j(shared/'manifest.json',{'native_libraries':{'lib':h(lib)},'openmp_library':{'path':str(omp),'sha256':h(omp)}})
    qualification=j(root/'qualification.json',{'status':'passed_numerical_qualification','worker_sha256':h(worker),'tests':{'small_molecules':[]}})
    validator=j(root/'validator.py',{'fixture':'validator'})
    plan=j(root/'plan.json',{'atom_ids':['A'],'elements':['H'],'parent':{'input_path':str(input_path),'input_sha256':h(input_path),'output_dir':str(parent),'cap_dir':str(cap),'cap_files_sha256':{n:h(cap/n) for n in ['capped-adduct.json','capped-adduct.sdf']}},
      'provider':{'worker':str(worker),'worker_sha256':h(worker),'runtime':sys.executable,'runtime_manifest':str(runtime),'runtime_manifest_sha256':h(runtime)},
      'qualification':{'result_path':str(qualification),'result_sha256':h(qualification)},'parent_geometry_validator':{'path':str(validator),'sha256':h(validator)},
      'request_settings':{'method':'B3LYP-D3BJ/Psi4-DZVP','threads':2,'max_memory_mb':8000,'max_wall_seconds':1,'charge':0,'spin':0}})
    pins=[plan,materializer,Path(sys.executable),source,sharedmanifest]
    config=j(root/'config.json',{'root':str(root),'plan_path':str(plan),'materializer':str(materializer),'materializer_runtime':sys.executable,
       'source_manifest':str(source),'shared_runtime_manifest':str(sharedmanifest),'pins':[{'path':str(p),'sha256':h(p)} for p in pins],
       'resp_dir':str(resp),'resp_fitter':str(fitter),'resp_controller_sha256':'resident-controller','claims_dir':str(root/'claims'),
       'wait_seconds':1,'poll_seconds':.01})
    monkeypatch.setattr(c,'resource_snapshot',lambda root:{'active_qm':[],'free_disk_bytes':30*c.GIB,'shared_lock_available':True})
    return SimpleNamespace(root=root,parent=parent,resp=resp,config=config,plan=plan)


def test_full_stub_pipeline_and_persistent_duplicate_claim(prepared):
    f=prepared;out=f.root/'run';r=c.execute(f.config,out)
    assert r['status']=='numerically_complete_research_reference'
    assert r['physical_acceptance'] is False and r['child_pid'] is None
    assert (out/'controller-source.py').is_file() and (out/'calculation/result.json').is_file()
    with pytest.raises(FileExistsError):c.execute(f.config,f.root/'duplicate')
    assert not (f.root/'duplicate').exists()


def test_stub_provider_failure_preserves_request_and_logs(prepared):
    f=prepared;(f.root/'fail-worker').write_text('fail');r=c.execute(f.config,f.root/'failed')
    assert r['status']=='failed' and '17' in r['error']
    assert (f.root/'failed/reference-input.json').exists()
    assert (f.root/'failed/native-controller.log').exists()


def test_wait_has_a_real_deadline_without_materializing(prepared):
    f=prepared;(f.parent/'result.json').unlink();(f.resp/'result.json').unlink()
    p=c.read(f.resp/'progress.json');p['status']='running';j(f.resp/'progress.json',p)
    cfg=c.read(f.config);cfg['wait_seconds']=.05;j(f.config,cfg)
    result=c.execute(f.config,f.root/'waiting')
    assert result['error_type']=='TimeoutError'
    assert not (f.root/'waiting/reference-input.json').exists()


def test_readiness_observed_after_wait_deadline_cannot_launch(prepared,monkeypatch):
    f=prepared;cfg=c.read(f.config);cfg['wait_seconds']=.05;j(f.config,cfg)
    def late_ready(*args):
        time.sleep(.07)
        return 'ready'
    monkeypatch.setattr(c,'dependencies',late_ready)
    result=c.execute(f.config,f.root/'late')
    assert result['error_type']=='TimeoutError'
    assert not (f.root/'late/materializer.log').exists()


def test_materialized_request_change_before_launch_is_rejected(prepared,monkeypatch):
    f=prepared;out=f.root/'prelaunch-mutation';calls=0
    def resources(root):
        nonlocal calls
        calls+=1
        if calls==2:
            path=out/'reference-input.json';request=c.read(path)
            request['coords_bohr']=[[9.,0.,0.]];j(path,request)
        return {'active_qm':[],'free_disk_bytes':30*c.GIB,'shared_lock_available':True}
    monkeypatch.setattr(c,'resource_snapshot',resources)
    result=c.execute(f.config,out)
    assert result['status']=='failed' and 'Pinned file changed' in result['error']
    assert result['request_sha256']!=h(out/'reference-input.json')
    assert not (out/'calculation').exists()


def test_self_consistent_changed_request_during_native_startup_cannot_pass(prepared):
    f=prepared;(f.root/'mutate-request-before-read').write_text('test mutation')
    out=f.root/'startup-mutation';result=c.execute(f.config,out)
    assert (out/'calculation/result.json').exists()
    assert result['status']=='failed' and 'Pinned file changed' in result['error']
    assert result['request_sha256']!=h(out/'reference-input.json')


def test_provider_input_must_match_original_digest_even_if_request_is_restored(prepared):
    f=prepared;out=f.root/'restored';result=c.execute(f.config,out)
    assert result['status']=='numerically_complete_research_reference'
    changed=c.read(out/'calculation/input.json');changed['coords_bohr']=[[9.,0.,0.]]
    j(out/'calculation/input.json',changed)
    provider_result=c.read(out/'calculation/result.json')
    provider_result['source']['input_sha256']=h(out/'calculation/input.json')
    j(out/'calculation/result.json',provider_result)
    assert h(out/'reference-input.json')==result['request_sha256']
    with pytest.raises(ValueError,match='source/input identity'):
        c.validate_reference(c.read(f.plan),out/'reference-input.json',out/'calculation',result['request_sha256'])


@pytest.mark.parametrize('mutation',['resp_input','parent_not_converged','wrong_parent_result','failed_resp','esp_checkpoint'])
def test_dependency_chain_mismatches_are_rejected(prepared,mutation):
    f=prepared;cfg=c.read(f.config);plan=c.read(f.plan)
    if mutation=='resp_input':
        p=c.read(f.resp/'progress.json');p['optimizer_input_sha256']='wrong';j(f.resp/'progress.json',p)
    elif mutation=='parent_not_converged':
        p=c.read(f.parent/'result.json');p['optimization']['converged']=False;j(f.parent/'result.json',p)
    elif mutation=='wrong_parent_result':
        p=c.read(f.resp/'result.json');p['stages']['parent']['result_sha256']='wrong';j(f.resp/'result.json',p)
    elif mutation=='failed_resp':
        p=c.read(f.resp/'progress.json');p['status']='failed';j(f.resp/'progress.json',p)
    else:
        p=c.read(f.resp/'exact-esp-input.json');p['initial_checkpoint']['result_sha256']='wrong';j(f.resp/'exact-esp-input.json',p)
    with pytest.raises(ValueError):c.dependencies(cfg,plan,{})


def test_observed_parent_cannot_change_while_resp_is_pending(prepared):
    f=prepared;(f.resp/'result.json').unlink();cfg=c.read(f.config);plan=c.read(f.plan);observed={}
    assert c.dependencies(cfg,plan,observed)=='waiting_for_resp'
    p=c.read(f.parent/'result.json');p['extra']='altered';j(f.parent/'result.json',p)
    with pytest.raises(ValueError,match='Pinned file changed'):c.dependencies(cfg,plan,observed)


def test_pinned_materializer_change_is_rejected_before_claim(prepared):
    f=prepared;(f.root/'materializer.py').write_text('changed')
    with pytest.raises(ValueError,match='Pinned file changed'):c.execute(f.config,f.root/'run')
    assert not (f.root/'claims').exists()


def test_actual_timeout_kills_worker_and_descendant(tmp_path):
    script=tmp_path/'hanging.py';pidfile=tmp_path/'grandchild.pid'
    script.write_text('import subprocess,sys,time\nfrom pathlib import Path\np=subprocess.Popen([sys.executable,"-c","import time;time.sleep(60)"])\nPath(sys.argv[1]).write_text(str(p.pid))\ntime.sleep(60)\n')
    pids=[]
    with pytest.raises(TimeoutError):
        c.run_child([sys.executable,str(script),str(pidfile)],tmp_path/'log',.25,.01,tmp_path,pids.append,lambda:None,reserve=0)
    assert len(pids)==1
    for pid in [pids[0],int(pidfile.read_text())]:
        # A reparented zombie has stopped executing even before init reaps it.
        state=subprocess.run(['ps','-p',str(pid),'-o','stat='],capture_output=True,text=True).stdout.strip()
        assert not state or state.startswith('Z')


def test_runtime_disk_reserve_stops_child(tmp_path,monkeypatch):
    pids=[];monkeypatch.setattr(c.shutil,'disk_usage',lambda root:SimpleNamespace(free=1))
    with pytest.raises(RuntimeError,match='reserve'):
        c.run_child([sys.executable,'-c','import time;time.sleep(60)'],tmp_path/'log',10,.01,tmp_path,pids.append,lambda:None)
    assert len(pids)==1


def test_real_controller_cancellation_cleans_child_and_holds_live_claim(prepared):
    f=prepared;(f.root/'hang-worker').write_text('hang')
    bootstrap=f.root/'bootstrap.py'
    bootstrap.write_text('''import importlib.util,sys
from pathlib import Path
s=importlib.util.spec_from_file_location('controller',sys.argv[1]);c=importlib.util.module_from_spec(s);s.loader.exec_module(c)
c.resource_snapshot=lambda root: {'active_qm': [], 'free_disk_bytes': 30*c.GIB, 'shared_lock_available': True}
c.execute(Path(sys.argv[2]),Path(sys.argv[3]))
''')
    out=f.root/'cancelled';p=subprocess.Popen([sys.executable,str(bootstrap),str(PATH),str(f.config),str(out)],start_new_session=True)
    child=None
    try:
        deadline=time.monotonic()+5
        while time.monotonic()<deadline:
            if (out/'progress.json').exists():
                progress=c.read(out/'progress.json')
                if progress['stage']=='reference_single_point' and progress['child_pid']:
                    child=progress['child_pid'];break
            if p.poll() is not None:pytest.fail('Controller exited before test cancellation')
            time.sleep(.01)
        assert child is not None
        with pytest.raises(BlockingIOError):c.execute(f.config,f.root/'concurrent-duplicate')
        p.terminate();p.wait(timeout=5)
        assert c.read(out/'result.json')['status']=='cancelled'
        state=subprocess.run(['ps','-p',str(child),'-o','stat='],capture_output=True,text=True).stdout.strip()
        assert not state or state.startswith('Z')
    finally:
        if p.poll() is None:p.kill();p.wait()


def test_repeated_term_int_cannot_interrupt_cleanup_of_ignoring_descendants(prepared):
    f=prepared;(f.root/'hang-worker').write_text('hang');(f.root/'ignore-termination').write_text('ignore')
    bootstrap=f.root/'repeated-cancellation-bootstrap.py'
    bootstrap.write_text('''import importlib.util,sys
from pathlib import Path
s=importlib.util.spec_from_file_location('controller',sys.argv[1]);c=importlib.util.module_from_spec(s);s.loader.exec_module(c)
c.resource_snapshot=lambda root: {'active_qm': [], 'free_disk_bytes': 30*c.GIB, 'shared_lock_available': True}
terminate=c.terminate_group
c.terminate_group=lambda child: terminate(child,grace=.5)
c.execute(Path(sys.argv[2]),Path(sys.argv[3]))
''')
    out=f.root/'repeat-cancelled';p=subprocess.Popen([sys.executable,str(bootstrap),str(PATH),str(f.config),str(out)],start_new_session=True)
    child=None;grandchild=None
    try:
        deadline=time.monotonic()+5
        while time.monotonic()<deadline:
            if (out/'progress.json').exists() and (f.root/'grandchild-ready').exists():
                progress=c.read(out/'progress.json')
                if progress['stage']=='reference_single_point' and progress['child_pid']:
                    child=progress['child_pid'];grandchild=int((f.root/'grandchild.pid').read_text());break
            if p.poll() is not None:pytest.fail('Controller exited before repeated cancellation')
            time.sleep(.01)
        assert child is not None and grandchild is not None
        p.send_signal(c.signal.SIGTERM);time.sleep(.1)
        p.send_signal(c.signal.SIGINT);time.sleep(.1)
        p.send_signal(c.signal.SIGTERM)
        p.wait(timeout=5)
        assert c.read(out/'result.json')['status']=='cancelled'
        for pid in [child,grandchild]:
            state=subprocess.run(['ps','-p',str(pid),'-o','stat='],capture_output=True,text=True).stdout.strip()
            assert not state or state.startswith('Z')
    finally:
        if p.poll() is None:p.kill();p.wait()
        # A failed test must not leave its deliberate ignoring descendants.
        if child is not None:
            try:os.killpg(child,c.signal.SIGKILL)
            except ProcessLookupError:pass
