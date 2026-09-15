"""Protocol safeguards plus opt-in calculations in the isolated QM runtime."""
import importlib.util
import hashlib
import json
import os
from pathlib import Path
import subprocess
import shutil

import numpy as np
import pytest

ROOT=Path(__file__).resolve().parents[2]
WORKER=ROOT/'backend'/'qm_worker.py'
spec=importlib.util.spec_from_file_location('isolated_qm_worker',WORKER)
qm=importlib.util.module_from_spec(spec)
spec.loader.exec_module(qm)
NATIVE=os.environ.get('DYNAMOL_QM_TESTS')=='1'
PYTHON=Path(os.environ.get('DYNAMOL_QM_PYTHON',ROOT/'.tools/qm/bin/python'))
native=pytest.mark.skipif(not NATIVE or not PYTHON.exists(),reason='Set DYNAMOL_QM_TESTS=1 to run private PySCF numerical validation')


def water():
    return {'schema_version':1,'atom_ids':['water/O','water/H1','water/H2'],
            'elements':['O','H','H'],'coords_angstrom':[[0,0,0],[0.75,0,0.59],[-0.76,0,0.60]],
            'charge':0,'spin':0,'method':'RHF','basis':'sto-3g','operations':['gradient']}


def peroxide():
    return {'schema_version':1,'atom_ids':['H-left','O-left','O-right','H-right'],
        'elements':['H','O','O','H'],'coords_angstrom':[[-0.30,0.92,0],[0,0,0],[1.45,0,0],[1.75,0,0.92]],
        'charge':0,'spin':0,'method':'RHF','basis':'sto-3g','operations':['gradient']}


def dihedral_optimization(target=60):
    return {'enabled':True,'max_steps':60,'freeze_atom_ids':['O-left','O-right'],
        'dihedral_constraints':[{'atom_ids':peroxide()['atom_ids'],'target_degrees':target}]}


@pytest.mark.parametrize('patch,message',[
    ({'atom_ids':['same','same','other']},'unique'),
    ({'coords_bohr':[[0,0,0]]*3},'exactly one'),
    ({'spin':1},'closed-shell'),
    ({'charge':True},'integer'),
    ({'coords_angstrom':[[0,0,0],[0,float('nan'),1],[0,1,0]]},'finite'),
    ({'method':'RKS','functional':'B3LYP'},'Ambiguous'),
    ({'method':'RKS'},'explicit functional'),
    ({'method':'RHF','functional':'B3LYPG'},'not applicable'),
    ({'operations':['esp']},'esp_points_bohr'),
    ({'threads':3},'threads'),
    ({'max_memory_mb':9000},'max_memory'),
    ({'ecp':{'Zn':'lanl2dz'}},'absent elements'),
    ({'optimization':{'enabled':True,'freeze_atom_ids':['not-present']}},'freeze_atom_ids'),
    ({'optimization':{'enabled':False,'max_steps':2}},'disabled'),
    ({'optimize':True},'Unknown'),
    ({'initial_checkpoint':{}},'exactly'),
    ({'initial_checkpoint':{'result_path':'relative','result_sha256':'0'*64,'checkpoint_path':'/tmp/scf.chk','checkpoint_sha256':'0'*64}},'absolute'),
    ({'initial_checkpoint':{'result_path':'/tmp/result.json','result_sha256':'invalid','checkpoint_path':'/tmp/scf.chk','checkpoint_sha256':'0'*64}},'SHA-256'),
])
def test_invalid_model_or_ignored_options_are_rejected(patch,message):
    with pytest.raises(ValueError,match=message):
        qm.validate_input(water()|patch)


def test_units_are_explicit_and_rks_gaussian_variant_is_preserved():
    d=water();d['coords_bohr']=d.pop('coords_angstrom')
    result=qm.validate_input(d|{'method':'RKS','functional':'B3LYPG'})
    assert 'coords_bohr' in result and 'coords_angstrom' not in result
    assert result['functional']=='B3LYPG'


def test_existing_output_is_never_overwritten(tmp_path):
    out=tmp_path/'evidence';out.mkdir();p=out/'result.json';p.write_text('original')
    with pytest.raises(FileExistsError):qm.run(tmp_path/'absent.json',out)
    assert p.read_text()=='original'


def test_invalid_input_writes_failed_result_without_arrays(tmp_path):
    p=tmp_path/'input.json';p.write_text(json.dumps(water()|{'spin':1}))
    result=qm.run(p,tmp_path/'failed')
    assert result['accepted'] is False and result['status']=='failed'
    assert result['phase']=='input' and 'closed-shell' in result['error']['message']
    assert not (tmp_path/'failed'/'arrays.npz').exists()


def test_hessian_atom_and_cartesian_axes_are_distinct():
    native=np.zeros((2,2,3,3));native[1,0,2,1]=7.3
    converted=qm.flatten_hessian(native,2)
    assert converted[5,1]==7.3
    assert np.count_nonzero(converted)==1
    with pytest.raises(ValueError):qm.flatten_hessian(np.zeros((6,6)),2)


@pytest.mark.parametrize('constraint,message',[
    ({'atom_ids':['H-left','O-left','O-right','O-right'],'target_degrees':60},'distinct'),
    ({'atom_ids':['H-left','O-left','O-right','missing'],'target_degrees':60},'existing'),
    ({'atom_ids':peroxide()['atom_ids'],'target_degrees':float('nan')},'finite'),
    ({'atom_ids':peroxide()['atom_ids'],'target_degrees':60,'force_constant':100},'exactly'),
])
def test_dihedral_malformed_targets_are_rejected(constraint,message):
    with pytest.raises(ValueError,match=message):
        qm.validate_input(peroxide()|{'optimization':dihedral_optimization()|{'dihedral_constraints':[constraint]}})


def test_conflicting_reversed_constraints_and_cartesian_freezes_are_rejected():
    opt=dihedral_optimization()
    opt['dihedral_constraints'].append({'atom_ids':peroxide()['atom_ids'][::-1],'target_degrees':-60})
    with pytest.raises(ValueError,match='Duplicate or conflicting'):qm.validate_input(peroxide()|{'optimization':opt})
    data=peroxide();data['atom_ids']+=['He'];data['elements']+=['He'];data['coords_angstrom']+=[[8,8,8]]
    opt=dihedral_optimization();opt['freeze_atom_ids']=peroxide()['atom_ids']
    with pytest.raises(ValueError,match='conflicting or redundant'):qm.validate_input(data|{'optimization':opt})


def test_dihedral_singularity_and_periodic_final_target_checks():
    data=peroxide();data['coords_angstrom'][0]=[-0.97,0,0]
    with pytest.raises(ValueError,match='Singular dihedral'):qm.validate_input(data|{'optimization':dihedral_optimization()})
    data['coords_angstrom'][0]=[-0.97,1e-7,0]
    with pytest.raises(ValueError,match='nearly linear'):qm.validate_input(data|{'optimization':dihedral_optimization()})
    coords=np.array([[0,1,0],[0,0,0],[1,0,0],[1,-1,0]])
    constraints=[{'atom_ids':['A','B','C','D'],'target_degrees':-180}]
    result=qm._dihedral_rows(coords,['A','B','C','D'],constraints,0.1)[0]
    assert result['target_satisfied'] and result['absolute_periodic_error_degrees']==0
    constraints[0]['target_degrees']=150
    result=qm._dihedral_rows(coords,['A','B','C','D'],constraints,0.1)[0]
    assert not result['target_satisfied'] and result['absolute_periodic_error_degrees']==30


def test_runtime_manifest_is_copied_as_declared_provenance(tmp_path,monkeypatch):
    prefix=tmp_path/'runtime';prefix.mkdir();out=tmp_path/'output';out.mkdir()
    monkeypatch.setattr(qm.sys,'prefix',str(prefix))
    assert qm.capture_runtime_manifest(out)=={'present':False}
    raw=b'{"declared_build":"test-only"}\n';(prefix/'dynamol-qm-runtime.json').write_bytes(raw)
    result=qm.capture_runtime_manifest(out)
    assert result['sha256']==hashlib.sha256(raw).hexdigest()
    assert (out/'runtime-manifest.json').read_bytes()==raw
    assert 'not verified' in result['meaning']


@pytest.fixture
def native_run(tmp_path,request):
    base=Path(os.environ['DYNAMOL_QM_AUDIT_DIR'])/request.node.name if 'DYNAMOL_QM_AUDIT_DIR' in os.environ else tmp_path
    base.mkdir(parents=True,exist_ok=True)
    def call(name,data,accepted=True):
        inp=base/f'{name}.json';inp.write_text(json.dumps(data,indent=2,allow_nan=False)+'\n')
        out=base/name
        process=subprocess.run([str(PYTHON),str(WORKER),str(inp),str(out)],capture_output=True,text=True,timeout=120)
        (base/f'{name}-process.log').write_text(process.stdout+'\n'+process.stderr)
        report=json.loads((out/'result.json').read_text())
        assert report['accepted'] is accepted,report
        assert process.returncode==(0 if accepted else 1)
        arrays=dict(np.load(out/'arrays.npz',allow_pickle=False)) if accepted else {}
        return report,arrays
    call.base=base
    return call


@native
def test_native_water_derivatives_match_full_finite_differences(native_run):
    data=water()|{'operations':['gradient','hessian','esp'],
                  'esp_points_bohr':[[3,2,4],[-4,3,2],[0,4,-3]]}
    report,arrays=native_run('central',data)
    step=1e-4
    coords=arrays['coords_bohr']
    gradient_fd=np.empty(9);hessian_fd=np.empty((9,9))
    for j in range(9):
        results=[]
        for sign,label in [(1,'plus'),(-1,'minus')]:
            shifted=coords.copy();shifted.flat[j]+=sign*step
            trial={k:v for k,v in water().items() if k!='coords_angstrom'}
            results.append(native_run(f'axis{j}-{label}',trial|{'coords_bohr':shifted.tolist()}))
        (plus,gplus),(minus,gminus)=results
        gradient_fd[j]=(plus['energy_hartree']-minus['energy_hartree'])/(2*step)
        hessian_fd[:,j]=(gplus['gradient_hartree_per_bohr']-gminus['gradient_hartree_per_bohr']).reshape(-1)/(2*step)
    gradient_error=float(np.max(np.abs(gradient_fd-arrays['gradient_hartree_per_bohr'].reshape(-1))))
    hessian_error=float(np.max(np.abs(hessian_fd-arrays['hessian_hartree_per_bohr2'])))
    assert gradient_error<2e-6
    assert hessian_error<2e-5
    assert report['scf']['converged'] and report['hessian']['response_solver_completed']
    assert report['esp']['passed']
    assert np.allclose(arrays['coords_bohr']*report['units']['bohr_to_angstrom'],arrays['coords_angstrom'],atol=1e-13)
    (native_run.base/'finite-difference-check.json').write_text(json.dumps({
        'step_bohr':step,'gradient_max_error_hartree_per_bohr':gradient_error,
        'hessian_max_error_hartree_per_bohr2':hessian_error,
        'gradient_tolerance':2e-6,'hessian_tolerance':2e-5,
        'all_nine_atom_major_coordinates_checked':True},indent=2)+'\n')


@native
def test_native_esp_sign_far_field_and_ecp_core_charge(native_run):
    # Be2+ has a closed 1s shell and an unambiguous positive 2/r asymptote.
    points=[[0,0,30],[0,0,60],[0,0,120]]
    ion={'atom_ids':['Be2+'],'elements':['Be'],'coords_bohr':[[0,0,0]],
         'charge':2,'spin':0,'method':'RHF','basis':'sto-3g','operations':['esp'],
         'esp_points_bohr':points}
    report,arrays=native_run('positive-ion',ion)
    assert np.max(np.abs(arrays['esp_hartree_per_e']*np.array([30,60,120])-2))<1e-8
    # Zn LANL2DZ removes its core electrons; using Z=30 as the ESP point
    # nuclear charge would produce the wrong far-field charge by 18 electrons.
    ecp=ion|{'atom_ids':['Zn2+'],'elements':['Zn'],'basis':'lanl2dz','ecp':{'Zn':'lanl2dz'}}
    report,arrays=native_run('ecp-ion',ecp)
    assert arrays['atomic_numbers'][0]==30
    assert arrays['effective_nuclear_charges'][0]<30
    assert abs(arrays['effective_nuclear_charges'].sum()-report['explicit_electron_count']-2)<1e-12
    assert np.max(np.abs(arrays['esp_hartree_per_e']*np.array([30,60,120])-2))<1e-7


@native
def test_native_rks_b3lypg_all_requested_properties(native_run):
    data={'atom_ids':['H1','H2'],'elements':['H','H'],'coords_bohr':[[0,0,0],[0,0,1.4]],
        'charge':0,'spin':0,'method':'RKS','functional':'B3LYPG','grid_level':3,
        'basis':'sto-3g','operations':['gradient','hessian','esp'],
        'esp_points_bohr':[[4,0,3],[-3,2,0]]}
    report,arrays=native_run('rks',data)
    assert report['method']['functional_requested']=='B3LYPG'
    assert arrays['hessian_hartree_per_bohr2'].shape==(6,6)
    assert report['esp']['max_absolute_error_hartree_per_e']<1e-8
    assert np.isfinite(arrays['gradient_hartree_per_bohr']).all()


@native
def test_native_optimization_preserves_frozen_cap_coordinates(native_run):
    data=water()|{'optimization':{'enabled':True,'max_steps':40,'freeze_atom_ids':['water/O']}}
    report,arrays=native_run('optimized',data)
    assert report['optimization']['converged']
    assert report['optimization']['max_frozen_displacement_bohr']<1e-5
    assert np.linalg.norm(arrays['coords_angstrom'][0])<1e-5
    assert np.linalg.norm(arrays['coords_angstrom']-np.array(data['coords_angstrom']))>0.01
    assert np.linalg.norm(arrays['gradient_hartree_per_bohr'][1:])<0.001


@native
def test_native_fixed_dihedrals_reach_signed_periodic_targets_and_honor_freezes(native_run):
    for label,target in [('positive',60),('periodic',300)]:
        data=peroxide()|{'optimization':dihedral_optimization(target)}
        result,arrays=native_run(label,data)
        final=result['optimization']['final_dihedral_constraints'][0]
        assert result['optimization']['converged'] and result['optimization']['dihedral_constraints_satisfied']
        assert final['target_degrees']==target and final['absolute_periodic_error_degrees']<0.1
        assert abs(qm.dihedral_degrees(arrays['coords_bohr'])-(60 if target==60 else -60))<0.1
        np.testing.assert_allclose(arrays['coords_angstrom'][[1,2]],np.array(data['coords_angstrom'])[[1,2]],atol=6e-6,rtol=0)
        assert result['optimization']['max_frozen_displacement_bohr']<1e-5
        constraint=(native_run.base/label/'constraints.txt').read_text()
        assert '$freeze\nxyz 2,3\n$set\ndihedral 1 2 3 4 ' in constraint
        assert result['requested_threads']==2 and 1<=result['actual_pyscf_threads']<=2
        assert result['method']['actual_pyscf_threads']==result['actual_pyscf_threads']


@native
def test_native_singular_dihedral_request_fails_before_optimization(native_run):
    data=peroxide();data['coords_angstrom'][0]=[-0.97,0,0]
    result,arrays=native_run('singular',data|{'optimization':dihedral_optimization()},accepted=False)
    assert result['phase']=='input' and 'Singular dihedral' in result['error']['message'] and not arrays


@native
def test_optimizer_convergence_flag_cannot_substitute_for_actual_final_dihedral(native_run):
    data=peroxide()|{'optimization':dihedral_optimization()}
    input_path=native_run.base/'input.json';input_path.write_text(json.dumps(data))
    output=native_run.base/'false-convergence'
    # A deliberately incorrect optimizer stub returns the unchanged90° input
    # while claiming convergence. The worker must independently reject it.
    script="""import importlib.util,sys
from pyscf.geomopt import geometric_solver
spec=importlib.util.spec_from_file_location('qm_worker_under_test',sys.argv[1])
qm=importlib.util.module_from_spec(spec);spec.loader.exec_module(qm)
geometric_solver.kernel=lambda mean_field,**kwargs:(True,mean_field.mol.copy())
qm.run(sys.argv[2],sys.argv[3])
"""
    subprocess.run([str(PYTHON),'-c',script,str(WORKER),str(input_path),str(output)],
        check=True,capture_output=True,text=True,timeout=30)
    result=json.loads((output/'result.json').read_text())
    assert not result['accepted'] and result['optimization']['converged']
    assert result['optimization']['dihedral_constraints_satisfied'] is False
    assert 'Final actual dihedral' in result['error']['message']
    assert not (output/'arrays.npz').exists()


@native
def test_native_nonconvergence_and_singular_esp_are_never_accepted(native_run):
    report,arrays=native_run('scf-limit',water()|{'max_scf_cycles':1},accepted=False)
    assert report['scf']['converged'] is False and not arrays
    report,_=native_run('singular-esp',water()|{'operations':['esp'],'esp_points_bohr':[[0,0,0]]},accepted=False)
    assert 'singular' in report['error']['message']
    report,_=native_run('optimization-limit',water()|{'optimization':{'enabled':True,'max_steps':1}},accepted=False)
    assert report['optimization']['converged'] is False
    report,_=native_run('unresolved-ecp',water()|{'ecp':{'H':'lanl2dz'}},accepted=False)
    assert 'ECP was not resolved' in report['error']['message']


def checkpoint_request(directory):
    return {key:str((directory/name).resolve()) for key,name in [('result_path','result.json'),('checkpoint_path','scf.chk')]} | {
        key:hashlib.sha256((directory/name).read_bytes()).hexdigest() for key,name in [('result_sha256','result.json'),('checkpoint_sha256','scf.chk')]}


@native
def test_native_df_checkpoint_initializes_independent_exact_scf_and_optimization(native_run):
    data=water()|{'operations':['gradient','esp'],'esp_points_bohr':[[3,2,4],[-4,3,2]]}
    source,source_arrays=native_run('df-source',data|{'density_fit':True})
    source_directory=native_run.base/'df-source'
    request=checkpoint_request(source_directory)
    cold,cold_arrays=native_run('cold-exact',data)
    warm,warm_arrays=native_run('warm-exact',data|{'initial_checkpoint':request})
    assert warm['scf']['converged'] and warm['initial_checkpoint']['accepted_as_initial_guess']
    assert warm['initial_checkpoint']['source_density_fit'] is True
    assert warm['initial_checkpoint']['requested_density_fit'] is False
    assert abs(warm['energy_hartree']-cold['energy_hartree'])<1e-10
    np.testing.assert_allclose(warm_arrays['gradient_hartree_per_bohr'],cold_arrays['gradient_hartree_per_bohr'],atol=1e-8,rtol=0)
    np.testing.assert_allclose(warm_arrays['esp_hartree_per_e'],cold_arrays['esp_hartree_per_e'],atol=1e-9,rtol=0)
    assert checkpoint_request(source_directory)==request  # Never overwrite the source.
    assert (native_run.base/'warm-exact'/'initial-checkpoint'/'scf.chk').read_bytes()==(source_directory/'scf.chk').read_bytes()
    optimized,arrays=native_run('warm-optimization',water()|{'initial_checkpoint':request,
        'optimization':{'enabled':True,'max_steps':40,'freeze_atom_ids':['water/O']}})
    assert optimized['optimization']['converged']
    assert optimized['initial_checkpoint']['requested_method_initial_scf_converged']
    assert optimized['optimization']['max_frozen_displacement_bohr']<1e-5


@native
def test_native_checkpoint_rejects_changed_identity_geometry_state_basis_and_hash(native_run):
    native_run('source',water()|{'operations':['energy']})
    request=checkpoint_request(native_run.base/'source')
    shifted=water();shifted['coords_angstrom'][1][0]+=0.001
    trials=[('geometry',shifted,'geometry differs'),
        ('ids',water()|{'atom_ids':['water/O','water/H2','water/H1']},'atom order'),
        ('charge',water()|{'charge':2},'charge'),
        ('basis',water()|{'basis':'6-31g'},'AO count'),
        ('method',water()|{'method':'RKS','functional':'B3LYPG'},'method differs')]
    for name,data,message in trials:
        result,_=native_run(name,data|{'initial_checkpoint':request},accepted=False)
        assert message in result['error']['message']
    result,_=native_run('hash',water()|{'initial_checkpoint':request|{'checkpoint_sha256':'0'*64}},accepted=False)
    assert 'hash mismatch' in result['error']['message']


@native
def test_native_checkpoint_checks_energy_and_orbital_density_not_only_file_hashes(native_run):
    native_run('source',water()|{'operations':['energy']})
    source=native_run.base/'source'
    def clone(name,mutation=None):
        folder=native_run.base/name;folder.mkdir()
        for file in ['result.json','arrays.npz','resolved-method.json','scf.chk']:shutil.copy2(source/file,folder/file)
        result=json.loads((folder/'result.json').read_text())
        # Earlier accepted worker artifacts did not bind their checkpoint hash;
        # explicit input hashes plus molecule/density/energy checks cover them.
        result.pop('checkpoint_sha256');result.pop('checkpoint_file')
        (folder/'result.json').write_text(json.dumps(result))
        if mutation:
            code="import h5py,sys\nwith h5py.File(sys.argv[1], 'r+') as f:\n    "+mutation
            subprocess.run([str(PYTHON),'-c',code,str(folder/'scf.chk')],check=True,capture_output=True,text=True,timeout=30)
        return checkpoint_request(folder)
    legacy=clone('legacy-artifacts')
    result,_=native_run('legacy-restart',water()|{'initial_checkpoint':legacy})
    assert result['initial_checkpoint']['accepted_as_initial_guess']
    for name,mutation,message in [
        ('energy',"f['scf/e_tot'][...]=float(f['scf/e_tot'][()])+1.0",'energy does not match'),
        ('occupancy',"f['scf/mo_occ'][0]=1.0",'invalid closed-shell'),
        ('coefficients',"f['scf/mo_coeff'][0,0]+=0.1",'not orthonormal')]:
        request=clone(name+'-artifacts',mutation)
        result,_=native_run(name+'-restart',water()|{'initial_checkpoint':request},accepted=False)
        assert message in result['error']['message']


@native
def test_native_rks_checkpoint_requires_same_functional_and_grid(native_run):
    data={'atom_ids':['H1','H2'],'elements':['H','H'],'coords_bohr':[[0,0,0],[0,0,1.4]],
        'charge':0,'spin':0,'method':'RKS','functional':'B3LYPG','grid_level':3,'basis':'sto-3g','operations':['energy']}
    source,_=native_run('source',data)
    request=checkpoint_request(native_run.base/'source')
    warm,_=native_run('same',data|{'initial_checkpoint':request})
    assert abs(warm['energy_hartree']-source['energy_hartree'])<1e-10
    for name,patch in [('functional',{'functional':'PBE'}),('grid',{'grid_level':4})]:
        result,_=native_run(name,data|patch|{'initial_checkpoint':request},accepted=False)
        assert 'resolved method differs' in result['error']['message']
