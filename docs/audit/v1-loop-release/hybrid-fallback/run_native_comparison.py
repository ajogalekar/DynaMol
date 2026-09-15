"""Freeze bounded native comparison inputs, supervise independent mode batches."""
import hashlib
import json
import os
from pathlib import Path
import signal
import subprocess
import sys
import time

import psutil

ROOT=Path(__file__).resolve().parents[4]
CACHE=Path.home()/'.cache/dynamol-research/v1-loop-release'
OUTPUT=CACHE/'hybrid-native-mc-v1'
NATIVE=Path.home()/'Library/Application Support/DynaMol/Runtimes/promod3-dfd2a559551cd2f2/bin/python'


def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def save(p,d):p.write_text(json.dumps(d,indent=2,allow_nan=False)+'\n')


def supervise(command,folder):
    env={k:v for k,v in os.environ.items() if not k.startswith(('PYTHON','DYLD_','OPENMM','PM3_'))}
    env.update(OMP_NUM_THREADS='2',OPENBLAS_NUM_THREADS='2',PM3_OPENMM_CPU_THREADS='2')
    start=time.monotonic();peak=0;reason=None
    with (folder/'worker.log').open('w') as log:
        proc=subprocess.Popen(command,stdout=log,stderr=subprocess.STDOUT,env=env,start_new_session=True)
        monitor=psutil.Process(proc.pid)
        try:
            while proc.poll() is None:
                try:peak=max(peak,sum(p.memory_info().rss for p in [monitor,*monitor.children(recursive=True)] if p.is_running()))
                except psutil.NoSuchProcess:pass
                if time.monotonic()-start>180 or peak>4*1024**3:
                    reason='wall_time' if time.monotonic()-start>180 else 'memory'
                    os.killpg(proc.pid,signal.SIGTERM)
                    try:proc.wait(timeout=3)
                    except subprocess.TimeoutExpired:os.killpg(proc.pid,signal.SIGKILL)
                    break
                time.sleep(.2)
        except BaseException:
            if proc.poll() is None:os.killpg(proc.pid,signal.SIGKILL)
            proc.wait();raise
        code=proc.wait()
    report={'returncode':code,'stop_reason':reason,'elapsed_seconds':time.monotonic()-start,
            'peak_sampled_rss_bytes':peak,'limits':{'seconds':180,'rss_bytes':4*1024**3,'threads':2},'command':command}
    save(folder/'supervisor.json',report)
    return report


def prepare():
    OUTPUT.mkdir(parents=True,exist_ok=False)
    impl=OUTPUT/'implementation';impl.mkdir()
    for p in [Path(__file__),Path(__file__).with_name('native_mc_adapter.py')]:
        (impl/p.name).write_bytes(p.read_bytes())
    cases=[]
    rootfinal=ROOT/'docs/audit/loop-fallback/runs/root-final'
    for name,source_id,modeled,spans,fragment_length in [
        ('8K5R','8899bea3019d447b',range(177,180),[[176,180],[175,181]],3),
        ('1UA2','3a661cfc670048c5',range(44,56),[[43,56],[42,57]],9),
    ]:
        folder=OUTPUT/name;folder.mkdir()
        prior=CACHE/'carbonyl-ranking-v1'/name
        config=json.loads((prior/'input.json').read_text())
        config.pop('research_include_carbonyl_anchors',None)
        config['input_pdb']=str(folder/'context.pdb')
        (folder/'context.pdb').write_bytes((prior/'context.pdb').read_bytes())
        save(folder/'input.json',config)
        source=rootfinal/name/'workspace/datasets'/source_id/'topology.pdb'
        (folder/'source.pdb').write_bytes(source.read_bytes())
        keys=[[c['chain_id'],r['resid'],r['insertion_code'],r['residue']] for c in config['chains'] for r in c['residues']
              if c['chain_id']=='A' and int(r['resid']) in modeled]
        assert len(keys)==len(modeled) and all(not r['observed'] for c in config['chains'] for r in c['residues'] if int(r['resid']) in modeled)
        cases.append({'case':name,'source':str(folder/'source.pdb'),'original_source':str(source),
                      'input':str(folder/'input.json'),'modeled_residue_keys':keys,'spans':spans,'fragment_length':fragment_length})
    control=ROOT/'docs/audit/preparation-fixtures/six_residues_intact.pdb'
    folder=OUTPUT/'masked-control';folder.mkdir()
    lines=control.read_text().splitlines(keepends=True)
    observed=[l for l in lines if not(l.startswith(('ATOM  ','HETATM')) and l[21]=='A' and int(l[22:26]) in (2,3,4))]
    (folder/'source.pdb').write_text(''.join(observed))
    (folder/'context.pdb').write_text(''.join(observed))
    letters=dict(zip('MET GLN ILE PHE VAL LYS'.split(),'MQIFVK'));rows=[]
    for number,name in enumerate(letters,1):rows.append({'chain':'A','resid':str(number),'insertion_code':'','residue':name,'one_letter':letters[name],'observed':number not in (2,3,4)})
    save(folder/'input.json',{'input_pdb':str(folder/'context.pdb'),'chains':[{'chain_id':'A','residues':rows}], 'max_res_extension':0,'requested_seed':2026})
    cases.append({'case':'masked-control','source':str(folder/'source.pdb'),'input':str(folder/'input.json'),
                  'modeled_residue_keys':[['A',str(n),'',rows[n-1]['residue']] for n in (2,3,4)],'spans':[[1,5]],'fragment_length':3,
                  'heldout_control_source':str(control),'heldout_control_sha256':sha(control),
                  'control_limit':'Masked coordinates excluded from generator input/scoring. Database may contain homologous or identical fragments; not a held-out-database accuracy benchmark.'})
    files=list(OUTPUT.glob('*/source.pdb'))+list(OUTPUT.glob('*/input.json'))+list(OUTPUT.glob('*/context.pdb'))+list(impl.glob('*.py'))
    plan={'cases':cases,'modes':['torsion_mc','fragment_mc'],'seeds':[2026,2027],'steps':5000,
          'source_sha256':{str(p):sha(p) for p in files},'maximum_candidates':20,
          'contexts':'Original stems and symmetric one-residue extension each side (2 observed internal residues); control original only.',
          'source_observed_N_CA_C_O_cap_A':1.0,'limits_per_case_mode':{'seconds':180,'rss_bytes':4*1024**3,'threads':2},
          'native_components':'Existing native MC reduced/cb_packing/clash weights and cooling/DirtyCCDCloser; no force-field parameter changes.',
          'selection':'Exactly one result per predeclared seed/context/mode. No adaptive resampling after quality outcomes.',
          'local_validation':'Independent CCTBX and original gross backbone/omega/cap checks; sidechains/full complex remain root-owned validation.',
          'app_ready':False,'physical_model_validated':False}
    save(OUTPUT/'plan.json',plan)
    return plan


def main():
    plan=prepare();reports=[]
    for case in plan['cases']:
        for mode in plan['modes']:
            folder=OUTPUT/case['case']/mode;folder.mkdir()
            command=[str(NATIVE),'-I','-B',str(OUTPUT/'implementation/native_mc_adapter.py'),'--plan',str(OUTPUT/'plan.json'),
                     '--case',case['case'],'--mode',mode,'--output',str(folder/'native')]
            result=supervise(command,folder);reports.append({'case':case['case'],'mode':mode,**result})
            save(OUTPUT/'supervised-batches.json',reports);print(json.dumps(reports[-1]),flush=True)
    assert all(sha(Path(p))==h for p,h in plan['source_sha256'].items())


if __name__=='__main__':main()
