"""Synthetic software-contract tests only: no QM, fitted charges or model acceptance."""
import copy
import importlib.util
import json
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch

MODULE=Path(__file__).with_name('materialize_6oim_recovery_esp.py')
spec=importlib.util.spec_from_file_location('recovery_materializer',MODULE)
m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m)


def controller():
    return dict(schema_version=1,status='numerically_complete_pending_geometry_review',child_pid=None,
        physical_acceptance=False,simulation_ready=False,spec_sha256=m.SPEC_SHA,controller_sha256=m.CONTROLLER_SHA,
        limits=dict(native_memory_mb=2000,max_rss_bytes=8_000_000_000,threads=2,outer_wall_seconds=21720),
        peak_sampled_rss_bytes=1)


def synthetic():
    # Two invented atom labels and values deliberately differ from the biological source.
    parent=dict(schema_version=1,atom_ids=['fixture:H1','fixture:H2'],elements=['H','H'],charge=0,spin=0,
        method='RHF',basis='6-31g*',density_fit=True,auxbasis='def2-universal-jkfit',
        coords_angstrom=[[0.,0.,0.],[0.,0.,1.]],operations=['gradient'],threads=2,max_memory_mb=2000,
        max_wall_seconds=21600,optimization={'enabled':True})
    coords=[[0.,0.,0.],[0.,0.,1.]];points=[[3.,0.,0.],[0.,3.,0.],[0.,0.,3.]]
    request=m.build_cold_request(parent,coords,points)
    method=dict(method='RHF',basis_requested='6-31g*',density_fit=False,basis_sha256='b'*64,
        ecp_sha256='e'*64,spherical_basis=True,symmetry=False,functional_requested=None,grid_level=None,ecp_requested=None)
    binding=dict(input_sha256='a'*64,worker_sha256=m.WORKER_SHA,runtime_sha256=m.RUNTIME_SHA,
                 arrays_sha256='c'*64,resolved_method_file_sha256='d'*64)
    result=dict(schema_version=1,accepted=True,status='completed',atom_ids=parent['atom_ids'][:],
        elements=parent['elements'][:],charge=0,spin_2S=0,method=method,energy_hartree=-1.,
        scf=dict(converged=True,max_cycles=100,cycles=5,conv_tol_hartree=1e-10,conv_tol_grad=1e-7),
        units=dict(energy='hartree',coords_bohr='bohr',esp_points_bohr='bohr',esp_hartree_per_e='hartree/e',bohr_to_angstrom=.52917721092),
        arrays_file='arrays.npz',resolved_method_file='resolved-method.json',runtime_manifest={'sha256':m.RUNTIME_SHA},
        **{k:v for k,v in binding.items() if k!='runtime_sha256'})
    arrays=dict(coords_bohr=copy.deepcopy(coords),esp_points_bohr=copy.deepcopy(points),esp_hartree_per_e=[.1,.2,.3])
    return request,result,parent,method,coords,points,arrays,m.COLD,binding


class ColdContracts(unittest.TestCase):
    def test_synthetic_contract_passes_without_acceptance_return(self):
        self.assertIsNone(m.validate_cold_esp(*synthetic()))

    def test_build_cold_exact_geometry_and_original_chemistry(self):
        request,_,parent,_,coords,points,*_=synthetic()
        before=copy.deepcopy(parent)
        self.assertFalse(request['density_fit'])
        self.assertEqual(request['coords_bohr'],coords);self.assertEqual(request['esp_points_bohr'],points)
        self.assertEqual(request['method'],'RHF');self.assertEqual(request['basis'],'6-31g*')
        self.assertEqual(request['operations'],['esp']);self.assertEqual(parent,before)
        for key in ['initial_checkpoint','optimization','coords_angstrom','auxbasis']:self.assertNotIn(key,request)
        request['coords_bohr'][0][0]=99
        self.assertEqual(coords[0][0],0.)

    def test_any_checkpoint_key_rejects_even_null_or_empty(self):
        for location in [0,1]:
            for value in [None,{},False,{'result_path':'old-result.json'}]:
                with self.subTest(location=location,value=value):
                    args=list(synthetic());args[location]['initial_checkpoint']=value
                    with self.assertRaisesRegex(ValueError,'checkpoint'):m.validate_cold_esp(*args)

    def test_unspecified_or_checkpoint_policy_rejects(self):
        for policy in [None,'','reuse_if_available','verified_checkpoint']:
            args=list(synthetic());args[7]=policy
            with self.subTest(policy=policy),self.assertRaises(ValueError):m.validate_cold_esp(*args)

    def test_request_modification_rejects(self):
        mutations=[('charge',False),('spin',0.0),('density_fit',0),('schema_version',True),('threads',2.0),
            ('method','RKS'),('basis','6-31g'),('max_memory_mb',8000),('optimization',{}),
            ('atom_ids',['fixture:H2','fixture:H1']),('elements',['H','He'])]
        for key,value in mutations:
            args=list(synthetic());args[0][key]=value
            with self.subTest(key=key),self.assertRaises(ValueError):m.validate_cold_esp(*args)

    def test_request_geometry_or_grid_mismatch_rejects(self):
        for key in ['coords_bohr','esp_points_bohr']:
            args=list(synthetic());args[0][key][0][0]+=.001
            with self.subTest(key=key),self.assertRaises(ValueError):m.validate_cold_esp(*args)

    def test_nonfinite_or_boolean_declared_geometry_rejects(self):
        for value in [float('nan'),float('inf'),True]:
            args=list(synthetic());args[4][0][0]=value;args[0]['coords_bohr']=copy.deepcopy(args[4])
            with self.subTest(value=value),self.assertRaises(ValueError):m.validate_cold_esp(*args)

    def test_result_state_completion_and_schema_reject(self):
        for key,value in [('accepted',False),('accepted',1),('status','running'),('status','failed'),
                          ('charge',False),('spin_2S',0.0),('schema_version',True),
                          ('atom_ids',['fixture:H2','fixture:H1']),('elements',['H','He'])]:
            args=list(synthetic());args[1][key]=value
            with self.subTest(key=key,value=value),self.assertRaises(ValueError):m.validate_cold_esp(*args)

    def test_scf_settings_or_nonconvergence_reject(self):
        for key,value in [('converged',False),('converged',1),('cycles',0),('cycles',101),('max_cycles',100.),
                          ('conv_tol_hartree',1e-6),('conv_tol_grad',1e-4)]:
            args=list(synthetic());args[1]['scf'][key]=value
            with self.subTest(key=key),self.assertRaises(ValueError):m.validate_cold_esp(*args)

    def test_resolved_method_modifications_reject(self):
        for key,value in [('method','RKS'),('density_fit',True),('density_fit',0),('basis_sha256','0'*64),
                          ('ecp_sha256','0'*64),('basis_requested','sto-3g'),('ecp_requested',{}),
                          ('spherical_basis',1),('symmetry',0),('functional_requested','B3LYP'),
                          ('grid_level',3),('auxbasis_sha256','a'*64)]:
            args=list(synthetic());args[1]['method']=copy.deepcopy(args[1]['method']);args[1]['method'][key]=value
            with self.subTest(key=key),self.assertRaises(ValueError):m.validate_cold_esp(*args)

    def test_immutable_binding_mismatch_reject(self):
        for key in ['input_sha256','worker_sha256','arrays_sha256','resolved_method_file_sha256']:
            args=list(synthetic());args[1][key]='0'*64
            with self.subTest(key=key),self.assertRaises(ValueError):m.validate_cold_esp(*args)
        args=list(synthetic());args[1]['runtime_manifest']['sha256']='0'*64
        with self.assertRaises(ValueError):m.validate_cold_esp(*args)

    def test_unreviewed_binding_worker_runtime_reject(self):
        for key in ['worker_sha256','runtime_sha256']:
            args=list(synthetic());args[8][key]='0'*64
            if key=='runtime_sha256':args[1]['runtime_manifest']['sha256']='0'*64
            else:args[1][key]='0'*64
            with self.subTest(key=key),self.assertRaises(ValueError):m.validate_cold_esp(*args)

    def test_missing_binding_or_bad_paths_reject(self):
        args=list(synthetic());args[8].pop('input_sha256')
        with self.assertRaises(ValueError):m.validate_cold_esp(*args)
        for key in ['arrays_file','resolved_method_file']:
            args=list(synthetic());args[1][key]='../other'
            with self.subTest(key=key),self.assertRaises(ValueError):m.validate_cold_esp(*args)

    def test_wrong_units_nonfinite_energy_or_esp_reject(self):
        args=list(synthetic());args[1]['units']['esp_hartree_per_e']='kcal/mol'
        with self.assertRaises(ValueError):m.validate_cold_esp(*args)
        for value in [float('nan'),float('inf'),True]:
            args=list(synthetic());args[1]['energy_hartree']=value
            with self.subTest(value=value),self.assertRaises(ValueError):m.validate_cold_esp(*args)
            args=list(synthetic());args[6]['esp_hartree_per_e'][0]=value
            with self.subTest(value=value),self.assertRaises(ValueError):m.validate_cold_esp(*args)

    def test_result_geometry_grid_or_length_mismatch_reject(self):
        for key in ['coords_bohr','esp_points_bohr']:
            args=list(synthetic());args[6][key][0][0]+=.000001
            with self.subTest(key=key),self.assertRaises(ValueError):m.validate_cold_esp(*args)
        args=list(synthetic());args[6]['esp_hartree_per_e'].pop()
        with self.assertRaises(ValueError):m.validate_cold_esp(*args)

    def test_boolean_result_coordinate_rejects_even_equal_to_zero(self):
        for key in ['coords_bohr','esp_points_bohr']:
            args=list(synthetic());args[6][key][0][1]=False
            with self.subTest(key=key),self.assertRaises(ValueError):m.validate_cold_esp(*args)


class AdmissionTests(unittest.TestCase):
    def setUp(self):
        scratch=m.CACHE/'software-test-scratch';scratch.mkdir(parents=True,exist_ok=True)
        self.temp=tempfile.TemporaryDirectory(dir=scratch);self.addCleanup(self.temp.cleanup)
        self.root=Path(self.temp.name)

    def test_controller_reference_metadata_only(self):m.validate_controller(controller())

    def test_running_failed_or_missing_child_pid_reject(self):
        for value in ['running','failed','queued','stopped']:
            data=controller();data['status']=value
            with self.subTest(status=value),self.assertRaises(ValueError):m.validate_controller(data)
        data=controller();del data['child_pid']
        with self.assertRaises(ValueError):m.validate_controller(data)

    def test_controller_strict_types_and_hashes_reject(self):
        for key,value in [('schema_version',True),('physical_acceptance',0),('simulation_ready',0),
                          ('child_pid',123),('spec_sha256','wrong'),('controller_sha256','wrong'),
                          ('peak_sampled_rss_bytes',0),('peak_sampled_rss_bytes',8_000_000_001)]:
            data=controller();data[key]=value
            with self.subTest(key=key),self.assertRaises(ValueError):m.validate_controller(data)
        data=controller();data['limits']['threads']=2.0
        with self.assertRaises(ValueError):m.validate_controller(data)

    def test_missing_terminal_emits_only_refusal_without_native_reads(self):
        self.early_refusal(None,'terminal parent absent')

    def test_failed_terminal_emits_only_refusal_without_native_reads(self):
        data=controller();data['status']='failed';self.early_refusal(data,'running, failed or nonterminal')

    def test_running_terminal_emits_only_refusal_without_native_reads(self):
        data=controller();data['status']='running';self.early_refusal(data,'running, failed or nonterminal')

    def early_refusal(self,data,expected):
        parent=self.root/'parent';parent.mkdir()
        if data is not None:(parent/'controller-result.json').write_text(json.dumps(data))
        plan=self.root/'plan.json';plan.write_text(json.dumps({'parent':str(parent),'checkpoint_policy':m.COLD}))
        output=self.root/'out'
        with patch.object(m,'PARENT',parent),patch.object(m,'CACHE',self.root),\
             patch.object(m,'validate_preparation_runtime',side_effect=AssertionError('No downstream reads allowed')):
            result=m.materialize(plan,output)
        self.assertEqual(result['status'],'refused');self.assertIn(expected,result['error'])
        self.assertFalse(result['request_emitted']);self.assertFalse(result['QM_launched']);self.assertFalse(result['RESP_launched'])
        self.assertEqual({p.name for p in output.iterdir()},{'result.json'})

    def test_pending_runtime_requires_explicit_independent_admission(self):
        with self.assertRaisesRegex(ValueError,'pending'):m.validate_preparation_runtime({'root_reviewed':False})

    def test_omitted_frozen_source_pin_refuses_before_native_inputs(self):
        parent=self.root/'parent';parent.mkdir();(parent/'controller-result.json').write_text(json.dumps(controller()))
        plan=self.root/'plan.json';plan.write_text(json.dumps({'parent':str(parent),'checkpoint_policy':m.COLD,
            'sources':str(m.CACHE/'prepared/sources'),'source_pins':{}}))
        with patch.object(m,'PARENT',parent),patch.object(m,'validate_preparation_runtime',return_value={}):
            result=m.materialize(plan,self.root/'out')
        self.assertEqual(result['status'],'refused');self.assertIn('frozen chemistry sources',result['error'])
        self.assertFalse(result['request_emitted']);self.assertEqual({p.name for p in (self.root/'out').iterdir()},{'result.json'})

    def test_publication_cannot_overwrite(self):
        path=self.root/'record.json';m.publish(path,{'a':1})
        with self.assertRaises(FileExistsError):m.publish(path,{'a':2})
        self.assertEqual(json.loads(path.read_text()),{'a':1})
        self.assertEqual(len(list(self.root.glob('.*.tmp'))),1)

    def test_interrupted_publication_retains_temporary(self):
        path=self.root/'record.json'
        with patch.object(m.os,'link',side_effect=OSError('synthetic ENOSPC')):
            with self.assertRaises(OSError):m.publish(path,{'synthetic':True})
        self.assertFalse(path.exists());self.assertEqual(len(list(self.root.glob('.*.tmp'))),1)

    def test_pinned_source_alteration_and_symlink_reject(self):
        path=self.root/'source';path.write_text('original');pin=m.sha(path)
        path.write_text('changed')
        with self.assertRaises(ValueError):m.raw(path,pin)
        link=self.root/'link';link.symlink_to(path)
        with self.assertRaises(ValueError):m.raw(link)

    def test_empty_runtime_package_file_requires_its_exact_empty_pin(self):
        path=self.root/'empty';path.write_bytes(b'')
        self.assertEqual(m.raw(path,m.hashlib.sha256(b'').hexdigest()),b'')
        with self.assertRaises(ValueError):m.raw(path)

    def test_output_outside_stable_cache_reject(self):
        with patch.object(m,'CACHE',self.root/'allowed'):
            with self.assertRaises(ValueError):m.materialize(self.root/'plan',self.root/'outside')
        self.assertFalse((self.root/'outside').exists())


if __name__=='__main__':unittest.main(verbosity=2)
