"""No-QM lineage and rejection tests; no synthetic accepted endpoint is written."""
import copy
import importlib.util
import inspect
from pathlib import Path
import numpy as np
import pytest

spec=importlib.util.spec_from_file_location('recovery_endpoint',Path(__file__).with_name('analyze_1mnc_recovery_endpoint_v1.py'))
m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m)


@pytest.fixture(scope='module')
def ctx():
    return m.core.load_context()


@pytest.fixture
def inputs(ctx):
    return copy.deepcopy(ctx['request']),m.load(m.RECOVERED,m.PINS[m.RECOVERED]),m.load(m.RECOVERY/'prepared/input.json',m.RECOVERY_INPUT_HASH)


def test_genuine_prepared_lineage_preserves_original_context(ctx):
    before=copy.deepcopy(ctx['request'])
    spec,lineage=m.load_lineage(ctx)
    assert ctx['request']==before
    assert lineage['free_atom_count']==94 and lineage['fixed_atom_count']==3
    assert lineage['original_policy_sha256']==m.core.POLICY_HASH
    assert spec['pins'][spec['input']]==m.RECOVERY_INPUT_HASH
    assert lineage['maximum_anchor_restore_angstrom']<6e-11


@pytest.mark.parametrize('change',[
    'cycle','accepted','status','source_hash','log_hash','source_ids','source_elements','source_charge','source_spin',
    'source_free_coordinate','source_cap_coordinate','charge','spin','ids','elements','extra_anchor',
    'removed_atom','proton_coordinate','free_coordinate','memory','basis','functional','method','grid',
    'checkpoint','extra_operation','gradient_threshold','scf_threshold','arbitrary_extra_key'])
def test_wrong_lineage_state_or_unreviewed_change_rejected(inputs,change):
    original,recovered,request=inputs
    if change=='cycle':recovered['completed_cycle']=9
    elif change=='accepted':recovered['accepted']=True
    elif change=='status':recovered['status']='converged'
    elif change=='source_hash':recovered['source_input_sha256']='0'*64
    elif change=='log_hash':recovered['source_native_log_sha256']='0'*64
    elif change=='source_ids':recovered['atom_ids'][0]='different'
    elif change=='source_elements':recovered['elements'][0]='N'
    elif change=='source_charge':recovered['charge']=0
    elif change=='source_spin':recovered['spin']=2
    elif change=='source_free_coordinate':recovered['coords_bohr'][0][0]+=.001
    elif change=='source_cap_coordinate':recovered['coords_bohr'][1][0]+=.001
    elif change=='charge':request['charge']=0
    elif change=='spin':request['spin']=2
    elif change=='ids':request['atom_ids'][0]='different'
    elif change=='elements':request['elements'][0]='N'
    elif change=='extra_anchor':request['optimization']['freeze_atom_ids'].append(request['atom_ids'][0])
    elif change=='removed_atom':request['atom_ids'].pop()
    elif change=='proton_coordinate':request['coords_angstrom'][request['elements'].index('H')][0]+=.01
    elif change=='free_coordinate':request['coords_angstrom'][0][0]+=.01
    elif change=='memory':request['max_memory_mb']=5000
    elif change=='basis':request['basis']='sto-3g'
    elif change=='functional':request['functional']='PBE'
    elif change=='method':request['method']='RHF'
    elif change=='grid':request['grid_level']=3
    elif change=='checkpoint':request['initial_checkpoint']=original['initial_checkpoint']
    elif change=='extra_operation':request['operations'].append('hessian')
    elif change=='gradient_threshold':request['optimization']['convergence_gmax']*=2
    elif change=='scf_threshold':request['scf_conv_tol']*=100
    elif change=='arbitrary_extra_key':request['guess_oxidation_state']=True
    with pytest.raises(m.EvidenceError):m.validate_lineage(original,recovered,request)


def test_actual_recovery_geometry_uses_original_source_frame(ctx,inputs):
    original,recovered,request=inputs
    xyz=np.asarray(request['coords_angstrom'])
    report=m.core.chemistry_report(ctx,xyz,np.asarray(recovered['gradient_hartree_per_bohr']))
    assert {a['code'] for a in report['alarms']}=={'free_gradient_rms_exceeded','free_gradient_max_exceeded'}
    assert not report['endpoint_checks_passed'] and not report['model_accuracy_validated']
    assert len(report['coordination'])==5 and len(report['zinc_centered_angles'])==10
    assert report['omitted_protein_contacts']['omitted_heavy_atoms']==1212
    assert report['omitted_protein_contacts']['source_sha256']==ctx['source_sha256']
    assert report['group_displacements']['whole_PLH_heavy']['source_frame_rmsd_angstrom']>0
    for row in report['source_frame_displacements']:
        index=ctx['index'][row['atom_id']]
        np.testing.assert_array_equal(row['vector_angstrom'],xyz[index]-np.asarray(original['coords_angstrom'][index]))


def test_native_reader_preserves_every_original_gate_except_three_declared_changes():
    original=inspect.getsource(m.core.read_endpoint)
    expected=original.replace('def read_endpoint(ctx, folder):','def read_native_endpoint(ctx, folder):').replace('INPUT_HASH','RECOVERY_INPUT_HASH')
    expected=expected.replace("runtime = load(ctx['base']/'source-manifest.json',SOURCE_MANIFEST_HASH)[str(ROOT/'.tools/qm-parallel/dynamol-qm-runtime.json')]['sha256']",'runtime = RECOVERY_RUNTIME_HASH')
    expected=expected.replace("get('max_memory_mb') == 8000","get('max_memory_mb') == 2000")
    assert inspect.getsource(m.read_native_endpoint)==expected
    assert m.sha(m.AUDIT/'analyze_1mnc_anchored_endpoint.py')==m.ORIGINAL_ANALYZER_HASH


def test_no_terminal_controller_does_not_read_native_evidence(tmp_path,monkeypatch):
    def forbidden(*args,**kwargs):raise AssertionError('Live native evidence was read')
    monkeypatch.setattr(m,'load',forbidden)
    assert m.read_recovery_endpoint({}, {'output':str(tmp_path)}) is None


@pytest.mark.parametrize('extra',['hessian.npy','esp.npz','initial-checkpoint','result-alternative.json','symlink','missing'])
def test_extra_missing_or_symlink_native_artifact_rejected(tmp_path,extra):
    for name in m.NATIVE_FILES:(tmp_path/name).write_text('synthetic rejection inventory only')
    if extra=='missing':(tmp_path/'arrays.npz').unlink()
    elif extra=='symlink':
        (tmp_path/'arrays.npz').unlink();(tmp_path/'arrays.npz').symlink_to(tmp_path/'input.json')
    elif extra=='initial-checkpoint':(tmp_path/extra).mkdir()
    else:(tmp_path/extra).write_text('unrequested')
    with pytest.raises(m.EvidenceError):m.validate_native_inventory(tmp_path)


@pytest.mark.parametrize('change',['failed','running','child','spec','controller','result','physical','simulation',
                                    'memory','rss_limit','threads','wall','peak','missing_peak',
                                    'missing_child','schema_bool','schema_float','physical_zero','simulation_zero',
                                    'missing_physical','missing_simulation','memory_float','rss_float',
                                    'threads_float','wall_float'])
def test_terminal_controller_mismatches_rejected(change):
    # No endpoint arrays/result file or successful endpoint is constructed.
    c={'schema_version':1,'status':'numerically_complete_pending_geometry_review','child_pid':None,
       'spec_sha256':m.SPEC_HASH,'controller_sha256':m.CONTROLLER_HASH,'result_sha256':'a'*64,
       'physical_acceptance':False,'simulation_ready':False,'peak_sampled_rss_bytes':7_000_000_000,
       'limits':{'native_memory_mb':2000,'max_rss_bytes':8_000_000_000,'threads':2,'outer_wall_seconds':21720}}
    if change=='failed':c['status']='failed'
    elif change=='running':c['status']='running'
    elif change=='child':c['child_pid']=123
    elif change=='spec':c['spec_sha256']='0'*64
    elif change=='controller':c['controller_sha256']='0'*64
    elif change=='result':c['result_sha256']='0'*64
    elif change=='physical':c['physical_acceptance']=True
    elif change=='simulation':c['simulation_ready']=True
    elif change=='memory':c['limits']['native_memory_mb']=5000
    elif change=='rss_limit':c['limits']['max_rss_bytes']=16_000_000_000
    elif change=='threads':c['limits']['threads']=4
    elif change=='wall':c['limits']['outer_wall_seconds']=43200
    elif change=='peak':c['peak_sampled_rss_bytes']=8_000_000_001
    elif change=='missing_peak':c.pop('peak_sampled_rss_bytes')
    elif change=='missing_child':c.pop('child_pid')
    elif change=='schema_bool':c['schema_version']=True
    elif change=='schema_float':c['schema_version']=1.0
    elif change=='physical_zero':c['physical_acceptance']=0
    elif change=='simulation_zero':c['simulation_ready']=0
    elif change=='missing_physical':c.pop('physical_acceptance')
    elif change=='missing_simulation':c.pop('simulation_ready')
    elif change=='memory_float':c['limits']['native_memory_mb']=2000.0
    elif change=='rss_float':c['limits']['max_rss_bytes']=8_000_000_000.0
    elif change=='threads_float':c['limits']['threads']=2.0
    elif change=='wall_float':c['limits']['outer_wall_seconds']=21720.0
    with pytest.raises(m.EvidenceError):m.validate_controller(c,'a'*64)


@pytest.mark.parametrize('key',['coords_bohr','coords_angstrom','gradient_hartree_per_bohr',
                              'atomic_numbers','effective_nuclear_charges'])
@pytest.mark.parametrize('imaginary',[0.0,1e-10])
def test_complex_array_rejected_before_any_float_coercion(key,imaginary):
    array=np.ones((97,3) if key not in ['atomic_numbers','effective_nuclear_charges'] else (97,),dtype=complex)
    array+=imaginary*1j
    with pytest.raises(m.EvidenceError):m.validate_real_arrays({key:array})


def test_real_array_dtype_guard_accepts_worker_storage_types_without_chemical_claim():
    m.validate_real_arrays({'coordinates':np.zeros((97,3),dtype=np.float64),
                            'atomic_numbers':np.ones(97,dtype=np.int64)})
