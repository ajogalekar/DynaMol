"""Orchestration controls without fragment generation or additional native modeling."""
import copy
import json
import os
from pathlib import Path
import signal
import subprocess

import numpy as np
from openmm import app, unit
import pytest

from backend import loop_modeling as lm


def scaffold():
    top = app.Topology()
    chain = top.addChain('A')
    for i in range(1, 4):
        residue = top.addResidue('ALA', chain, str(i))
        for name, symbol in [('N','N'),('CA','C'),('C','C'),('O','O'),('CB','C'),('H','H')]:
            top.addAtom(name, app.element.Element.getBySymbol(symbol), residue)
    # A retained water can reuse a protein chain ID and must remain outside the alignment.
    water = top.addResidue('HOH', top.addChain('A'), '7')
    top.addAtom('O', app.element.oxygen, water)
    xyz = np.arange(top.getNumAtoms()*3).reshape(-1,3)/100
    observed = {lm._key(a) for a in top.atoms() if a.residue.id != '2' and a.name != 'CB'}
    return top, xyz, observed


@pytest.fixture
def runtime(tmp_path, monkeypatch):
    prefix = tmp_path/'runtime'
    (prefix/'conda-meta').mkdir(parents=True)
    (prefix/'bin').mkdir()
    python = prefix/'bin/python'
    python.write_text('private python fixture')
    python.chmod(0o700)
    (prefix/'lib/plugins').mkdir(parents=True)
    (prefix/'share/promod3').mkdir(parents=True)
    (prefix/'share/promod3/fragment.dat').write_bytes(b'finite fixture database')
    for name, version in lm._PINS.items():
        (prefix/f'conda-meta/{name}-{version}-fixture.json').write_text(json.dumps(dict(name=name,version=version,build='fixture')))
    monkeypatch.setenv('DYNAMOL_PROMOD3', str(prefix))
    return prefix


def fake_candidate(top):
    return {'status':'candidate','versions':{'promod3':'3.6.0','ost':'2.11.1'},'atoms':[{'identity':list(lm._key(a)), 'xyz_nm':[a.index+1,2,3]}
        for a in top.atoms() if a.residue.id=='2' and a.element != app.element.hydrogen]}


def test_exact_transplant_and_context_excludes_unknown_loop_and_hydrogens(runtime,tmp_path,monkeypatch):
    top, xyz, observed = scaffold()
    before = xyz.copy()
    env_top=app.Topology(); r=env_top.addResidue('NA',env_top.addChain('A'),'99')
    env_top.addAtom('NA',app.element.sodium,r)
    environment=app.Modeller(env_top,np.array([[1.25,2.5,3.75]])*unit.nanometer)
    commands=[]
    def run(command, folder, runtime, check_cancel):
        commands.append(command)
        request=json.loads(Path(command[-2]).read_text())
        assert [r['observed'] for r in request['chains'][0]['residues']]==[True,False,True]
        context=app.PDBFile(request['input_pdb'])
        residues=list(context.topology.residues())
        assert [(r.name,r.id) for r in residues]==[('ALA','1'),('ALA','3'),('HOH','1'),('NA','1')]
        assert len({c.id for c in context.topology.chains()})==3
        assert all(a.element != app.element.hydrogen for a in context.topology.atoms())
        # The observed residues' repaired CB atoms remain part of the scoring context.
        assert sum(a.name=='CB' for a in context.topology.atoms())==2
        Path(command[-1]).write_text(json.dumps(fake_candidate(top)))
        return 0
    monkeypatch.setattr(lm,'_run',run)
    new, report=lm.generate_loop_model(top,xyz*unit.nanometer,observed,tmp_path/'attempt',2026,environment=environment)
    after=new.value_in_unit(unit.nanometer)
    moved=[a.index for a in top.atoms() if a.residue.id=='2' and a.element!=app.element.hydrogen]
    fixed=[a.index for a in top.atoms() if a.index not in moved]
    np.testing.assert_array_equal(xyz,before)
    np.testing.assert_array_equal(after[fixed],before[fixed])
    assert not np.array_equal(after[moved],before[moved])
    assert commands[0][1:3]==['-I','-B']
    assert report['status']=='candidate'
    assert report['seed_applied'] is False
    assert {m['source'] for m in report['environment_context_mapping']}=={'scaffold','environment'}
    persisted=json.loads((tmp_path/'attempt/loop-model-provenance.json').read_text())
    assert persisted==report
    assert persisted['files']['loop-model-context.pdb']['sha256']
    assert persisted['databases'][0]['sha256']


@pytest.mark.parametrize('change,message',[
    (lambda r:r['atoms'].pop(),'missing required'),
    (lambda r:r['atoms'].append(copy.deepcopy(r['atoms'][0])),'duplicate or unexpected'),
    (lambda r:r['atoms'][0]['identity'].__setitem__(1,'1'),'duplicate or unexpected'),
    (lambda r:r['atoms'][0]['identity'].__setitem__(4,'H'),'duplicate or unexpected'),
    (lambda r:r['atoms'][0].update(xyz_nm=[float('nan'),1,2]),'finite triples'),
    (lambda r:r['atoms'][0].update(xyz_nm=[1,2]),'finite triples'),
    (lambda r:r.update(status='failed',error='No fragments'),'No fragments'),
])
def test_reject_candidate_inventory_without_mutating_input(change,message):
    top, xyz, _=scaffold();before=xyz.copy();result=fake_candidate(top)
    target={lm._key(a):a.index for a in top.atoms() if a.residue.id=='2' and a.element!=app.element.hydrogen}
    change(result)
    with pytest.raises(ValueError,match=message):lm._transplant(result,target,xyz)
    np.testing.assert_array_equal(xyz,before)


def test_failed_worker_preserves_failure_and_input_evidence(runtime,tmp_path,monkeypatch):
    top,xyz,observed=scaffold()
    def failed(command,*args):
        Path(command[-1]).write_text(json.dumps({'status':'failed','error':'database cannot close gap'}))
        return 1
    monkeypatch.setattr(lm,'_run',failed)
    with pytest.raises(ValueError,match='cannot close gap'):
        lm.generate_loop_model(top,xyz,observed,tmp_path/'attempt',2026)
    saved=json.loads((tmp_path/'attempt/loop-model-provenance.json').read_text())
    assert saved['status']=='failed' and saved['files']['loop-model-output.json']['sha256']


def test_invalid_runtime_versions_fail_read_only(runtime):
    marker=runtime/'conda-meta/openmm-8.5.1-fixture.json'
    marker.write_text(json.dumps(dict(name='openmm',version='8.6.1')))
    before=marker.read_bytes()
    status=lm.loop_runtime_status()
    assert status['available'] is False and '8.5.1' in status['error']
    assert marker.read_bytes()==before


def test_missing_runtime_is_actionable_and_read_only(tmp_path,monkeypatch):
    path=tmp_path/'not-installed'
    monkeypatch.setenv('DYNAMOL_PROMOD3',str(path))
    assert lm.loop_runtime_status()['available'] is False
    assert not path.exists()


@pytest.mark.parametrize('mode',['duplicate_protein_chain','duplicate_resid','unsupported_missing','noncanonical_id'])
def test_ambiguous_input_fails_before_native_call(runtime,tmp_path,monkeypatch,mode):
    top,xyz,observed=scaffold()
    if mode=='unsupported_missing':list(top.residues())[1].name='TPO'
    elif mode=='noncanonical_id':list(top.residues())[1].id='0002'
    elif mode=='duplicate_resid':
        list(top.residues())[1].id='1';list(top.residues())[1].name='GLY'
    else:
        extra=top.addResidue('ALA',top.addChain('A'),'9')
        top.addAtom('N',app.element.nitrogen,extra);xyz=np.vstack([xyz,[1,2,3]])
    monkeypatch.setattr(lm,'_run',lambda *a:pytest.fail('must reject before native launch'))
    with pytest.raises(ValueError):lm.generate_loop_model(top,xyz,observed,tmp_path/'attempt',2026)


def test_worker_version_mismatch_rejects_candidate(runtime,tmp_path,monkeypatch):
    top,xyz,observed=scaffold()
    def run(command,*args):
        result=fake_candidate(top);result['versions']['promod3']='different'
        Path(command[-1]).write_text(json.dumps(result));return 0
    monkeypatch.setattr(lm,'_run',run)
    with pytest.raises(ValueError,match='unexpected ProMod3'):
        lm.generate_loop_model(top,xyz,observed,tmp_path/'attempt',2026)


class Process:
    pid=424242
    returncode=None
    def __init__(self,waits):self.waits=iter(waits);self.reaped=False
    def poll(self):return None
    def wait(self,timeout):
        value=next(self.waits)
        if value=='timeout':raise subprocess.TimeoutExpired('fixture',timeout)
        self.reaped=True;self.returncode=value;return value


@pytest.mark.parametrize('mode',['cancel','timeout'])
def test_cancel_and_timeout_terminate_group_and_reap(tmp_path,monkeypatch,mode):
    process=Process(['timeout',-9]);signals=[];calls=[]
    def popen(command,**kwargs):calls.append(kwargs);return process
    monkeypatch.setattr(lm.subprocess,'Popen',popen)
    monkeypatch.setattr(lm.os,'killpg',lambda pid,sig:signals.append((pid,sig)))
    monkeypatch.setattr(lm.time,'sleep',lambda _:None)
    monkeypatch.setenv('PYTHONPATH','untrusted')
    monkeypatch.setenv('DYLD_LIBRARY_PATH','untrusted')
    monkeypatch.setenv('OPENMM_PLUGIN_DIR','untrusted')
    checks=0
    def check():
        nonlocal checks
        checks+=1
        if mode=='cancel' and checks>1:raise RuntimeError('cancelled fixture')
    clock=iter([0,181]);monkeypatch.setattr(lm.time,'monotonic',lambda:next(clock))
    with pytest.raises(RuntimeError if mode=='cancel' else TimeoutError):
        lm._run(['private-python'],tmp_path,{'prefix':'/private/runtime'},check)
    assert signals==[(process.pid,signal.SIGTERM),(process.pid,signal.SIGKILL)]
    assert process.reaped and calls[0]['start_new_session'] is True
    env=calls[0]['env']
    assert 'PYTHONPATH' not in env and 'DYLD_LIBRARY_PATH' not in env
    assert env['OPENMM_PLUGIN_DIR']=='/private/runtime/lib/plugins'
    assert env['OPENMM_CPU_THREADS']==env['PM3_OPENMM_CPU_THREADS']=='2'


def test_cancel_before_launch_never_spawns(tmp_path,monkeypatch):
    monkeypatch.setattr(lm.subprocess,'Popen',lambda *a,**k:pytest.fail('must not spawn'))
    def cancel():raise RuntimeError('cancelled')
    with pytest.raises(RuntimeError,match='cancelled'):
        lm._run(['private-python'],tmp_path,{'prefix':'unused'},cancel)
