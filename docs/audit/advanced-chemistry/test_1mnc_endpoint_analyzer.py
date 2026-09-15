"""Synthetic numerical/geometry failure tests. These do not calculate QM."""
from __future__ import annotations

import copy
import hashlib
import importlib.util
import json
from pathlib import Path

import numpy as np
import pytest

SPEC=importlib.util.spec_from_file_location('mmp_endpoint',Path(__file__).with_name('analyze_1mnc_anchored_endpoint.py'))
m=importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(m)


@pytest.fixture(scope='module')
def ctx():
    return m.load_context()


def analyze(ctx,changes=None,gradient=None):
    xyz=np.array(ctx['request']['coords_angstrom'])
    if changes:
        for atom,position in changes.items():xyz[ctx['index'][atom]]=position
    gradient=np.zeros((97,3)) if gradient is None else gradient
    return m.chemistry_report(ctx,xyz,gradient)


def codes(report):
    return {r['code'] for r in report['alarms']}


def test_synthetic_zero_gradient_control_is_not_a_chemical_model(ctx):
    result=analyze(ctx)
    assert result['endpoint_checks_passed']
    assert not result['model_accuracy_validated'] and not result['full_preparation_ready']
    assert len(result['coordination'])==5 and len(result['zinc_centered_angles'])==10
    assert len(result['atomic_gradient_rows'])==97
    assert len(result['hydrogen_parent_checks'])==50
    assert len(result['covalent_bond_distances'])==96
    assert len(result['proper_ligand_heavy_torsions'])==36
    assert len(result['source_frame_displacements'])==97
    assert result['omitted_protein_contacts']['omitted_heavy_atoms']==1212
    assert result['omitted_protein_contacts']['contacts_within_cutoff_at_either_geometry']


def test_large_cap_reaction_force_is_visible_and_not_free_gradient(ctx):
    gradient=np.zeros((97,3));gradient[1]=[1.,-2.,3.]
    result=analyze(ctx,gradient=gradient)
    assert result['endpoint_checks_passed']
    assert result['gradient_statistics']['free']['max_atomic_vector']==0
    assert result['gradient_statistics']['frozen_caps']['max_atomic_vector']==pytest.approx(np.sqrt(14))
    assert result['atomic_gradient_rows'][1]['gradient_hartree_per_bohr']==[1.,-2.,3.]
    assert result['atomic_gradient_rows'][1]['force_hartree_per_bohr']==[-1.,2.,-3.]
    assert not result['model_accuracy_validated']


def test_one_large_free_force_is_not_hidden_by_rms(ctx):
    gradient=np.zeros((97,3));gradient[ctx['index']['E:280::PLH:C11'],0]=.00046
    result=analyze(ctx,gradient=gradient)
    assert 'free_gradient_max_exceeded' in codes(result)
    assert 'free_gradient_rms_exceeded' not in codes(result)


def test_many_small_free_forces_fail_rms_even_below_max(ctx):
    gradient=np.zeros((97,3));gradient[:,0]=.00031
    result=analyze(ctx,gradient=gradient)
    assert 'free_gradient_rms_exceeded' in codes(result)
    assert 'free_gradient_max_exceeded' not in codes(result)


def test_cap_movement_cannot_be_hidden_by_rigid_alignment(ctx):
    atom=ctx['request']['optimization']['freeze_atom_ids'][0]
    point=np.array(ctx['request']['coords_angstrom'][ctx['index'][atom]])+[.001,0,0]
    assert 'frozen_cap_displacement_exceeded' in codes(analyze(ctx,{atom:point}))


def test_named_zinc_donor_loss_is_detected(ctx):
    atom='E:280::PLH:O1'
    point=np.array(ctx['request']['coords_angstrom'][ctx['index'][atom]])+[0,0,3.]
    assert 'named_coordination_distance_outside_bounds' in codes(analyze(ctx,{atom:point}))


def test_new_ligand_oxygen_near_zinc_is_detected(ctx):
    zn=np.array(ctx['request']['coords_angstrom'][ctx['index']['B:281::ZN:ZN']])
    assert 'undeclared_N_O_near_zinc' in codes(analyze(ctx,{'E:280::PLH:O4':zn+[2.,0,0]}))


def test_zinc_angle_change_is_detected_with_preserved_distance(ctx):
    zn=np.array(ctx['request']['coords_angstrom'][ctx['index']['B:281::ZN:ZN']])
    donor=np.array(ctx['request']['coords_angstrom'][ctx['index']['A:218::HID:NE2']])
    result=analyze(ctx,{'A:218::HID:NE2':2*zn-donor})
    assert 'coordination_angle_change' in codes(result)
    target=next(r for r in result['coordination'] if r['atom_ids'][1]=='A:218::HID:NE2')
    assert abs(target['change_angstrom'])<1e-12


def test_hydroxamate_proton_transfer_is_not_hidden_by_name(ctx):
    oxygen=np.array(ctx['request']['coords_angstrom'][ctx['index']['E:280::PLH:O1']])
    result=analyze(ctx,{'E:280::PLH:H24':oxygen+[0,0,.98]})
    assert 'ambiguous_or_transferred_hydrogen' in codes(result)
    assert 'unprotonated_donor_has_close_hydrogen' in codes(result)


def test_stereo_inversion_erases_inherited_source_tags(ctx):
    xyz=np.array(ctx['request']['coords_angstrom'])
    a,b='E:280::PLH:C2','E:280::PLH:C4'
    result=analyze(ctx,{a:xyz[ctx['index'][b]].copy(),b:xyz[ctx['index'][a]].copy()})
    row=next(r for r in result['ligand_stereochemistry'] if r['atom_id']=='E:280::PLH:C3')
    assert not row['sign_preserved'] and row['final_graph_based_CIP']=='S'
    assert 'ligand_stereochemistry_changed_or_ambiguous' in codes(result)


def test_nonbonded_collision_and_ligand_drift_are_visible(ctx):
    xyz=np.array(ctx['request']['coords_angstrom'])
    result=analyze(ctx,{'E:280::PLH:C16':xyz[ctx['index']['A:218::HID:CG']]+[0,0,.1]})
    assert {'nonbonded_heavy_collision','retained_heavy_max_displacement_exceeded'} <= codes(result)


def test_whole_ligand_translation_is_not_removed_by_ligand_alignment(ctx):
    xyz=np.array(ctx['request']['coords_angstrom'])
    result=analyze(ctx,{a:xyz[ctx['index'][a]]+[2.1,0,0] for a in ctx['ligand_ids']})
    assert result['group_displacements']['whole_PLH_heavy']['source_frame_rmsd_angstrom']==pytest.approx(2.1)
    assert {'retained_heavy_source_frame_rmsd_exceeded','retained_heavy_max_displacement_exceeded'} <= codes(result)


def test_distal_ligand_collision_with_omitted_protein_is_detected(ctx):
    atom='E:280::PLH:C16'
    omitted=ctx['omitted_protein'][0]
    result=analyze(ctx,{atom:omitted['coords_angstrom']})
    assert 'omitted_protein_heavy_collision' in codes(result)


@pytest.mark.parametrize('bad',['nan_coordinates','nan_gradient','missing_atom'])
def test_nonfinite_or_wrong_sized_arrays_cannot_pass(ctx,bad):
    xyz=np.array(ctx['request']['coords_angstrom']);g=np.zeros((97,3))
    if bad=='nan_coordinates':xyz[0,0]=np.nan
    if bad=='nan_gradient':g[0,0]=np.nan
    if bad=='missing_atom':xyz=xyz[:-1]
    with pytest.raises(m.EvidenceError):m.chemistry_report(ctx,xyz,g)


def dump(path,data):
    path.write_text(json.dumps(data,indent=2)+'\n')


def synthetic_endpoint(ctx,tmp_path):
    """Explicitly fake worker output for schema checks; no native checkpoint."""
    folder=tmp_path/'SYNTHETIC-NEVER-QM';folder.mkdir()
    df=Path(ctx['request']['initial_checkpoint']['result_path']).parent
    result=copy.deepcopy(m.load(df/'result.json'))
    result.update(input_sha256=m.INPUT_HASH,energy_hartree=-1.0,density_electron_count=372.0)
    result['scf'].update(conv_tol_hartree=1e-10,conv_tol_grad=1e-7)
    settings=m.load(ctx['base']/'validated-input.json')['optimization']
    indices=[ctx['index'][a] for a in settings['freeze_atom_ids']]
    constraint=('$freeze\nxyz '+','.join(str(i+1) for i in indices)+'\n').encode()
    (folder/'constraints.txt').write_bytes(constraint)
    result['optimization']={'requested':True,'converged':True,'settings':settings,
        'constraint_atom_indices_zero_based':indices,'constraint_file_sha256':hashlib.sha256(constraint).hexdigest(),
        'dihedral_constraints_satisfied':True,'final_dihedral_constraints':[],
        'steps':[{'cycle':1,'scf_converged':True,'energy_hartree':-1.0}]}
    result['resources']['max_wall_seconds']=21600
    (folder/'input.json').write_bytes((ctx['base']/'input.json').read_bytes())
    (folder/'scf.chk').write_bytes(b'SYNTHETIC TEST BYTES: not an HDF5 checkpoint; never use for QM')
    (folder/'resolved-method.json').write_bytes((df/'resolved-method.json').read_bytes())
    (folder/'runtime-manifest.json').write_bytes((df/'runtime-manifest.json').read_bytes())
    xyz=np.array(ctx['request']['coords_angstrom']);numbers=np.array([m.Z[e] for e in ctx['elements']])
    np.savez(folder/'arrays.npz',coords_angstrom=xyz,coords_bohr=xyz/m.BOHR,gradient_hartree_per_bohr=np.zeros((97,3)),
             atomic_numbers=numbers,effective_nuclear_charges=numbers.astype(float))
    for name,key in [('arrays.npz','arrays_sha256'),('scf.chk','checkpoint_sha256'),('resolved-method.json','resolved_method_file_sha256')]:
        result[key]=m.sha(folder/name)
    dump(folder/'result.json',result)
    return folder,result


def test_valid_synthetic_artifact_envelope(ctx,tmp_path):
    folder,_=synthetic_endpoint(ctx,tmp_path)
    endpoint=m.read_endpoint(ctx,folder)
    assert endpoint['arrays']['coords_angstrom'].shape==(97,3)
    assert endpoint['result']['energy_hartree']==-1.0


@pytest.mark.parametrize('field,change',[
    ('atom_ids',lambda r:r['atom_ids'].__setitem__(0,'WRONG-ID')),
    ('charge',lambda r:r.__setitem__('charge',0)),
    ('spin',lambda r:r.__setitem__('spin_2S',2)),
    ('convergence',lambda r:r['optimization'].__setitem__('converged',False)),
    ('hidden_constraint',lambda r:r['optimization']['settings']['freeze_atom_ids'].append('B:281::ZN:ZN')),
    ('input_hash',lambda r:r.__setitem__('input_sha256','0'*64)),
    ('density_count',lambda r:r.__setitem__('density_electron_count',371.0)),
    ('unconverged_step',lambda r:r['optimization']['steps'][0].__setitem__('scf_converged',False)),
])
def test_false_metadata_cannot_become_an_endpoint(ctx,tmp_path,field,change):
    folder,result=synthetic_endpoint(ctx,tmp_path);change(result);dump(folder/'result.json',result)
    with pytest.raises(m.EvidenceError):m.read_endpoint(ctx,folder)


def test_altered_file_fails_even_when_result_says_completed(ctx,tmp_path):
    folder,_=synthetic_endpoint(ctx,tmp_path)
    with (folder/'arrays.npz').open('ab') as f:f.write(b'TAMPERED')
    with pytest.raises(m.EvidenceError,match='checksum'):m.read_endpoint(ctx,folder)


def test_different_resolved_method_fails_even_with_updated_file_hash(ctx,tmp_path):
    folder,result=synthetic_endpoint(ctx,tmp_path);method=m.load(folder/'resolved-method.json')
    method['functional_resolved_sha256']='0'*64;dump(folder/'resolved-method.json',method)
    result['resolved_method_file_sha256']=m.sha(folder/'resolved-method.json');dump(folder/'result.json',result)
    with pytest.raises(m.EvidenceError,match='method changed'):m.read_endpoint(ctx,folder)


@pytest.mark.parametrize('mutation',['nonfinite_gradient','wrong_atomic_number','units','unrequested_hessian'])
def test_consistently_rehashed_invalid_arrays_are_rejected(ctx,tmp_path,mutation):
    folder,result=synthetic_endpoint(ctx,tmp_path)
    with np.load(folder/'arrays.npz') as payload:arrays={k:payload[k].copy() for k in payload.files}
    if mutation=='nonfinite_gradient':arrays['gradient_hartree_per_bohr'][0,0]=np.nan
    if mutation=='wrong_atomic_number':arrays['atomic_numbers'][0]=6
    if mutation=='units':arrays['coords_bohr']*=10
    if mutation=='unrequested_hessian':arrays['hessian_hartree_per_bohr2']=np.zeros((291,291))
    np.savez(folder/'arrays.npz',**arrays)
    result['arrays_sha256']=m.sha(folder/'arrays.npz');dump(folder/'result.json',result)
    with pytest.raises(m.EvidenceError):m.read_endpoint(ctx,folder)


def test_extra_native_constraint_is_rejected_despite_correct_requested_settings(ctx,tmp_path):
    folder,result=synthetic_endpoint(ctx,tmp_path)
    (folder/'constraints.txt').write_text('$freeze\nxyz 2,17,32,97\n')
    result['optimization']['constraint_file_sha256']=m.sha(folder/'constraints.txt');dump(folder/'result.json',result)
    with pytest.raises(m.EvidenceError,match='constraint file'):m.read_endpoint(ctx,folder)


def test_changed_declared_gradient_units_are_rejected(ctx,tmp_path):
    folder,result=synthetic_endpoint(ctx,tmp_path)
    result['units']['gradient_hartree_per_bohr']='kcal/mol/angstrom';dump(folder/'result.json',result)
    with pytest.raises(m.EvidenceError,match='units differ'):m.read_endpoint(ctx,folder)


def test_torsion_sign_and_wrapping_against_analytic_right_angle():
    points=np.array([[1.,0,0],[0,0,0],[0,0,1],[0,1,1]])
    assert m.torsion(points)==pytest.approx(90.)
    points[-1]=[0,-1,1]
    assert m.torsion(points)==pytest.approx(-90.)


def test_absent_endpoint_returns_without_reading_result(ctx,tmp_path,monkeypatch):
    def forbidden(*a,**k):raise AssertionError('No result should be read')
    monkeypatch.setattr(m,'load',forbidden)
    assert m.read_endpoint(ctx,tmp_path/'absent') is None


def test_endpoint_alarm_prevents_any_conventional_request(ctx,tmp_path):
    folder,_=synthetic_endpoint(ctx,tmp_path);endpoint=m.read_endpoint(ctx,folder)
    g=np.ones((97,3));analysis=m.chemistry_report(ctx,endpoint['arrays']['coords_angstrom'],g)
    with pytest.raises(m.EvidenceError,match='all endpoint checks'):m.conventional_request(ctx,endpoint,analysis)


def test_conventional_request_preserves_actual_bohr_geometry_and_state(ctx,tmp_path):
    folder,_=synthetic_endpoint(ctx,tmp_path);endpoint=m.read_endpoint(ctx,folder)
    analysis=m.chemistry_report(ctx,endpoint['arrays']['coords_angstrom'],endpoint['arrays']['gradient_hartree_per_bohr'])
    request=m.conventional_request(ctx,endpoint,analysis)
    assert 'coords_angstrom' not in request and 'auxbasis' not in request
    assert np.array_equal(request['coords_bohr'],endpoint['arrays']['coords_bohr'])
    assert request['charge']==1 and request['spin']==0
    assert request['density_fit'] is False and request['optimization']=={'enabled':False}
    assert request['operations']==['gradient']
    assert request['initial_checkpoint']['result_sha256']==m.sha(folder/'result.json')
