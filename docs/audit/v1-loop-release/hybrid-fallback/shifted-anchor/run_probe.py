"""Freeze all nine saved inputs before one native closure per input."""
import hashlib,json,os,signal,subprocess,time
from pathlib import Path
import psutil
ROOT=Path('/Users/ashujo/Documents/Science/DynaMol')
CACHE=Path('/Users/ashujo/.cache/dynamol-research/v1-loop-release')
OUT=CACHE/'shifted-anchor-native-v1';PRIOR=CACHE/'boundary-conditioned-native-v1'
NATIVE=Path('/Users/ashujo/Library/Application Support/DynaMol/Runtimes/promod3-dfd2a559551cd2f2/bin/python')
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def save(p,d):Path(p).write_text(json.dumps(d,indent=2,allow_nan=False)+'\n')
def supervise(command,folder):
    env={k:v for k,v in os.environ.items() if not k.startswith(('PYTHON','DYLD_','OPENMM','PM3_'))}
    env.update(OMP_NUM_THREADS='2',OPENBLAS_NUM_THREADS='2',PM3_OPENMM_CPU_THREADS='2')
    started=time.monotonic();peak=0;reason=None
    with (folder/'worker.log').open('w') as log:
        proc=subprocess.Popen(command,stdout=log,stderr=subprocess.STDOUT,env=env,start_new_session=True);monitor=psutil.Process(proc.pid)
        try:
            while proc.poll() is None:
                try:peak=max(peak,sum(p.memory_info().rss for p in [monitor,*monitor.children(recursive=True)] if p.is_running()))
                except psutil.NoSuchProcess:pass
                if time.monotonic()-started>120 or peak>4*1024**3:
                    reason='wall_time' if time.monotonic()-started>120 else 'memory';os.killpg(proc.pid,signal.SIGTERM)
                    try:proc.wait(timeout=3)
                    except subprocess.TimeoutExpired:os.killpg(proc.pid,signal.SIGKILL)
                    break
                time.sleep(.1)
        except BaseException:
            if proc.poll() is None:os.killpg(proc.pid,signal.SIGKILL)
            proc.wait();raise
        code=proc.wait()
    result={'returncode':code,'stop_reason':reason,'elapsed_seconds':time.monotonic()-started,'peak_sampled_rss_bytes':peak,
       'limits':{'seconds':120,'rss_bytes':4*1024**3,'threads':2},'command':command}
    save(folder/'supervisor.json',result);return result

def main():
    OUT.mkdir(exist_ok=False);impl=OUT/'implementation';impl.mkdir()
    for p in [Path(__file__),Path(__file__).with_name('native_probe.py')]:
        (impl/p.name).write_bytes(p.read_bytes())
    helper=impl/'frozen_native_helper.py';helper.write_bytes((PRIOR/'implementation/native_mc_adapter.py').read_bytes())
    prior=json.loads((PRIOR/'plan.json').read_text())
    assert all(sha(p)==h for k in ['source_sha256','implementation_sha256'] for p,h in prior[k].items())
    cases=[c for c in prior['cases'] if c['case'] in ['1UA2','masked-control']];attempts=[]
    files=[PRIOR/'plan.json',*impl.iterdir()]
    for c in cases:
        for protocol in prior['protocols']:
            for seed in c['seeds']:
                a,b=c['spans'][0];p=PRIOR/c['case']/protocol['name']/'native'/f'{a}-{b}-seed-{seed}'/'candidate.json'
                assert p.exists()
                attempts.append({'name':protocol['name']+'-seed-'+str(seed),'case':c['case'],'source_protocol':protocol['name'],
                  'seed':seed,'input_candidate':str(p),'input_candidate_sha256':sha(p)})
                files.append(p)
        files.extend([Path(c['source']),Path(c['input']),Path(json.loads(Path(c['input']).read_text())['input_pdb'])])
        if c.get('heldout_control_source'):files.append(Path(c['heldout_control_source']))
    assert len(attempts)==9 and len({a['input_candidate'] for a in attempts})==9
    plan={'cases':cases,'attempts':attempts,'helper':str(helper),'frozen_sha256':{str(p):sha(p) for p in files},
      'selection':'All six original-span 1UA2 and three masked-control candidates from the previous fixed comparison; no selection by quality outcomes.',
      'maximum_native_closure_calls':9,'maximum_closure_calls_per_attempt':1,'native_maximum_iterations_per_call':1000,
      'native_c_stem_N_CA_C_rmsd_cutoff_A':.1,'limits_per_case':{'seconds':120,'rss_bytes':4*1024**3,'threads':2},
      'intervention':'Fit first observed N/CA/C; coherent sequential=True first-psi rotation to source carbonyl orientation; use first modeled residue as virtual native CCD anchor with coherent preceding residue; close suffix once and recombine without a transform.',
      'invariants':'No MC resampling, original source identities/environment, original 1 Angstrom observed N/CA/C/O cap, independent CCTBX and geometry gates. No Cartesian O restoration or added forces.',
      'screen':'All produced output backbones independently evaluated after frozen native run; full-complex validation is separate.',
      'app_ready':False,'physical_model_validated':False,'md_or_qm':False}
    save(OUT/'plan.json',plan);reports=[]
    for case in cases:
        folder=OUT/case['case'];folder.mkdir()
        command=[str(NATIVE),'-I','-B',str(impl/'native_probe.py'),'--plan',str(OUT/'plan.json'),'--case',case['case'],'--output',str(folder/'native')]
        result=supervise(command,folder);reports.append(dict(case=case['case'],**result));save(OUT/'supervised-batches.json',reports);print(json.dumps(reports[-1]),flush=True)
    assert all(sha(p)==h for p,h in plan['frozen_sha256'].items())
if __name__=='__main__':main()
