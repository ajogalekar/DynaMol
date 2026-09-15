"""The audit helper must isolate native CMAP lifecycle and preserve evidence."""
import importlib.util
import json
from pathlib import Path
from types import SimpleNamespace

import pytest

MODULE=Path(__file__).resolve().parents[2]/'docs/audit/advanced-chemistry/native_sander_parity.py'
spec=importlib.util.spec_from_file_location('native_sander_parity_audit',MODULE)
helper=importlib.util.module_from_spec(spec)
spec.loader.exec_module(helper)


def test_every_pose_uses_a_separate_single_pose_native_invocation(tmp_path,monkeypatch):
    seen=[]
    poses=[{'name':f'pose-{i}','angstrom':[[i,0,0]]} for i in range(3)]
    def run(command,**kwargs):
        request=json.loads(Path(command[2]).read_text())
        assert len(request['poses'])==1
        assert kwargs['timeout']==120
        seen.append(request)
        Path(command[5]).write_text(json.dumps([{'name':request['poses'][0]['name'],'energy_kj_mol':0}]))
        return SimpleNamespace(returncode=0,stdout='',stderr='')
    monkeypatch.setattr(helper.subprocess,'run',run)
    actual=helper.isolated_sander_values(Path('/fixture.prmtop'),poses,tmp_path,Path('/native/python'))
    assert [r['poses'][0] for r in seen]==poses
    assert [x['name'] for x in actual]==[p['name'] for p in poses]
    assert len(list(tmp_path.glob('*-input.json')))==3
    assert len(list(tmp_path.glob('*.log')))==3


def test_prior_evidence_prevents_native_rerun(tmp_path,monkeypatch):
    prior=tmp_path/'sander-pose-00-input.json';prior.write_text('preserve this failed request')
    monkeypatch.setattr(helper.subprocess,'run',lambda *a,**k:pytest.fail('Must not overwrite or rerun'))
    with pytest.raises(FileExistsError):
        helper.isolated_sander_values(Path('/fixture.prmtop'),[{'name':'bound'}],tmp_path,Path('/native/python'))
    assert prior.read_text()=='preserve this failed request'


def test_wrong_native_pose_is_rejected(tmp_path,monkeypatch):
    def run(command,**kwargs):
        Path(command[5]).write_text(json.dumps([{'name':'a different pose'}]))
        return SimpleNamespace(returncode=0,stdout='',stderr='')
    monkeypatch.setattr(helper.subprocess,'run',run)
    with pytest.raises(ValueError,match='not the requested pose'):
        helper.isolated_sander_values(Path('/fixture.prmtop'),[{'name':'bound'}],tmp_path,Path('/native/python'))
    assert (tmp_path/'sander-pose-00-native.json').exists()
