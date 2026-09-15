"""Native-file-format software fixtures only; no SCF, ESP, RESP or biological result."""
import copy
import ast
import importlib.util
import json
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch
import numpy as np

BASE=Path(__file__).parent
s=importlib.util.spec_from_file_location('reader',BASE/'read_6oim_cold_resp_evidence.py');r=importlib.util.module_from_spec(s);s.loader.exec_module(r)
s=importlib.util.spec_from_file_location('fixtures',BASE/'test_materialize_6oim_recovery_esp.py');fixtures=importlib.util.module_from_spec(s);s.loader.exec_module(fixtures)


class NativeFormatTests(unittest.TestCase):
    def setUp(self):
        self.m=r.core.materializer_module();scratch=self.m.CACHE/'software-test-scratch';scratch.mkdir(exist_ok=True)
        self.temp=tempfile.TemporaryDirectory(dir=scratch);self.addCleanup(self.temp.cleanup);self.root=Path(self.temp.name)

    def fixture(self):
        request,base_result,parent,_,coords,points,_,_,_=fixtures.synthetic()
        xyz=np.array(coords);nuclei=np.array([1,1],dtype=np.int64)
        basis={'H':[[0,[1.,1.]]]}  # Explicitly invented metadata basis for two labelled fixture atoms.
        method=dict(method='RHF',basis_requested='6-31g*',basis_resolved=basis,ecp_requested=None,ecp_resolved={},
                    basis_sha256=r.digest(basis),ecp_sha256=r.digest({}),spherical_basis=True,symmetry=False,
                    density_fit=True,functional_requested=None,grid_level=None,requested_threads=2,actual_pyscf_threads=2,
                    auxbasis_resolved={'synthetic':True},auxbasis_sha256=r.digest({'synthetic':True}))
        context={'materializer':self.m,'source':self.m.CACHE/'prepared/sources','request':parent,'method':method,
                 'result':{'explicit_electron_count':2,'nao':2},'arrays':{'coords_bohr':xyz,'atomic_numbers':nuclei,
                 'effective_nuclear_charges':nuclei.copy()}}
        folder=self.root/'synthetic-native';folder.mkdir()
        raw=(json.dumps(request,indent=2)+'\n').encode();(folder/'input.json').write_bytes(raw)
        materialized={'request':request,'files':{'exact-esp-input.json':raw}}
        resolved=r.worker_module(context).validate_input(request)
        (folder/'resolved-input.json').write_text(json.dumps(resolved)+'\n')
        (folder/'native.log').write_text('SYNTHETIC SOFTWARE FORMAT FIXTURE. No native calculation performed.\n')
        runtime=self.m.CACHE.parent/'6oim-local-recovery-v2-memory2000/run/optimization/runtime-manifest.json'
        # Use only the installed declared manifest; do not inspect the running native output.
        runtime=Path('/Users/ashujo/.cache/dynamol-runtimes/qm-parallel-local-v1/dynamol-qm-runtime.json')
        (folder/'runtime-manifest.json').write_bytes(r.data(runtime,self.m.RUNTIME_SHA))
        (folder/'progress.json').write_text(json.dumps({'phase':'completed','status':'completed','accepted':True})+'\n')
        native_method={k:v for k,v in method.items() if k not in ('auxbasis_resolved','auxbasis_sha256')};native_method['density_fit']=False
        (folder/'resolved-method.json').write_text(json.dumps(native_method)+'\n')
        arrays={'coords_bohr':xyz,'coords_angstrom':xyz*.52917721092,'atomic_numbers':nuclei,'effective_nuclear_charges':nuclei.copy(),
                'esp_points_bohr':np.array(points),'esp_hartree_per_e':np.array([.1,.2,.3])}
        np.savez_compressed(folder/'arrays.npz',**arrays)
        (folder/'scf.chk').write_bytes(b'SYNTHETIC UNUSABLE CHECKPOINT. No density or reuse permitted.\n')
        result=copy.deepcopy(base_result)
        result.update(phase='completed',method={k:v for k,v in native_method.items() if not k.endswith('_resolved')},
            input_sha256=r.sha(folder/'input.json'),arrays_sha256=r.sha(folder/'arrays.npz'),
            resolved_method_file_sha256=r.sha(folder/'resolved-method.json'),checkpoint_file='scf.chk',checkpoint_sha256=r.sha(folder/'scf.chk'),
            optimization={'requested':False,'converged':None,'settings':resolved['optimization'],'steps':[]},
            explicit_electron_count=2,density_electron_count=2.,nao=2,requested_threads=2,actual_pyscf_threads=2,
            arrays={k:{'shape':list(v.shape),'dtype':str(v.dtype)} for k,v in arrays.items()},
            resources={'threads':2,'max_memory_mb':2000,'max_wall_seconds':14400,'requested_threads':2,'actual_pyscf_threads':2,
                'memory_limit_kind':'PySCF allocation setting, not an operating-system hard limit'},
            esp={'method':'Batched int3c2e/fakemol_for_charges; independent int1e_rinv','passed':True,'batch_size':64,
                 'checked_point_indices':[0,1,2],'max_absolute_error_hartree_per_e':0.,'tolerance_hartree_per_e':1e-8})
        result['units'].update(coords_angstrom='angstrom',effective_nuclear_charges='e')
        (folder/'result.json').write_text(json.dumps(result)+'\n')
        return folder,materialized,context

    def mutate_result(self,folder,fun):
        result=r.load(folder/'result.json');fun(result);(folder/'result.json').write_text(json.dumps(result)+'\n')

    def test_tiny_synthetic_native_format(self):
        folder,materialized,context=self.fixture();evidence=r.read_native_esp(folder,materialized,context)
        self.assertFalse(evidence['physical_acceptance']);self.assertFalse(evidence['checkpoint_reused'])

    def test_resources_agree_with_actual_pinned_producer_statements(self):
        folder,materialized,context=self.fixture()
        tree=ast.parse(r.data(context['source']/'qm_worker.py',self.m.WORKER_SHA))
        run=next(x for x in tree.body if isinstance(x,ast.FunctionDef) and x.name=='run')
        calculate=next(x for x in tree.body if isinstance(x,ast.FunctionDef) and x.name=='calculate')
        assignment=next(x for x in ast.walk(run) if isinstance(x,ast.Assign) and
                        any(ast.unparse(t)=="report['resources']" for t in x.targets))
        update=next(x for x in calculate.body if isinstance(x,ast.Expr) and isinstance(x.value,ast.Call) and
                    ast.unparse(x.value.func)=="report['resources'].update")
        namespace={'report':{},'data':materialized['request'],'actual_threads':2}
        # Execute only these two exact dictionary bookkeeping statements; no native imports or call.
        exec(compile(ast.Module(body=[assignment,update],type_ignores=[]),'<pinned producer resources>','exec'),namespace)
        result=r.load(folder/'result.json');self.assertEqual(r.digest(result['resources']),r.digest(namespace['report']['resources']))
        r.read_native_esp(folder,materialized,context)
        del result['resources']['memory_limit_kind'];(folder/'result.json').write_text(json.dumps(result)+'\n')
        with self.assertRaisesRegex(ValueError,'resource metadata'):r.read_native_esp(folder,materialized,context)

    def test_checkpoint_key_or_failed_scf_result_rejects(self):
        for key,value in [('initial_checkpoint',None),('accepted',False),('phase','scf')]:
            folder,materialized,context=self.fixture();self.mutate_result(folder,lambda d:d.update({key:value}))
            with self.subTest(key=key),self.assertRaises(ValueError):r.read_native_esp(folder,materialized,context)
            for p in folder.iterdir():p.unlink()
            folder.rmdir()

    def test_nuclear_state_and_ao_contradictions_reject(self):
        for key,value in [('explicit_electron_count',2.),('explicit_electron_count',4),('nao',3),('density_electron_count',1.5),
                          ('density_electron_count',True),('spin_2S',0.)]:
            folder,materialized,context=self.fixture();self.mutate_result(folder,lambda d:d.update({key:value}))
            with self.subTest(key=key),self.assertRaises(ValueError):r.read_native_esp(folder,materialized,context)
            for p in folder.iterdir():p.unlink()
            folder.rmdir()

    def test_unreported_array_nuclei_change_rejects_even_with_new_hash(self):
        folder,materialized,context=self.fixture()
        with np.load(folder/'arrays.npz') as a:values={k:a[k].copy() for k in a.files}
        values['atomic_numbers'][0]=2;np.savez_compressed(folder/'arrays.npz',**values)
        self.mutate_result(folder,lambda d:d.update(arrays_sha256=r.sha(folder/'arrays.npz')))
        with self.assertRaisesRegex(ValueError,'nuclear identity'):r.read_native_esp(folder,materialized,context)

    def test_contradictory_optimization_metadata_rejects(self):
        folder,materialized,context=self.fixture();self.mutate_result(folder,lambda d:d['optimization'].update(requested=0))
        with self.assertRaisesRegex(ValueError,'optimization'):r.read_native_esp(folder,materialized,context)

    def test_wrong_progress_type_or_missing_artifact_rejects(self):
        folder,materialized,context=self.fixture();(folder/'progress.json').write_text(json.dumps({'phase':'completed','status':'completed','accepted':1}))
        with self.assertRaises(ValueError):r.read_native_esp(folder,materialized,context)
        (folder/'scf.chk').unlink()
        with self.assertRaises(ValueError):r.read_native_esp(folder,materialized,context)

    def test_esp_crosscheck_invalid_indices_threshold_or_failure_rejects(self):
        for key,value in [('passed',1),('passed',False),('checked_point_indices',[0,2]),('checked_point_indices',[0.,1,2]),
                          ('batch_size',64.),('max_absolute_error_hartree_per_e',1e-6),('max_absolute_error_hartree_per_e',-1.),
                          ('max_absolute_error_hartree_per_e',float('nan')),('tolerance_hartree_per_e',1e-6)]:
            folder,materialized,context=self.fixture();self.mutate_result(folder,lambda d:d['esp'].update({key:value}))
            with self.subTest(key=key,value=value),self.assertRaises(ValueError):r.read_native_esp(folder,materialized,context)
            for p in folder.iterdir():p.unlink()
            folder.rmdir()

    def test_unreviewed_extra_resolved_method_key_rejects_even_rehashed(self):
        folder,materialized,context=self.fixture();method=r.load(folder/'resolved-method.json');method['solvation']='extra'
        (folder/'resolved-method.json').write_text(json.dumps(method)+'\n')
        self.mutate_result(folder,lambda d:d.update(resolved_method_file_sha256=r.sha(folder/'resolved-method.json'),
            method={k:v for k,v in method.items() if not k.endswith('_resolved')}))
        with self.assertRaisesRegex(ValueError,'Complete resolved'):r.read_native_esp(folder,materialized,context)

    def test_wrong_arrays_or_resource_metadata_rejects(self):
        for field,change in [('arrays',lambda d:d['coords_bohr'].update(dtype='float32')),
                             ('resources',lambda d:d.update(threads=2.))]:
            folder,materialized,context=self.fixture();self.mutate_result(folder,lambda d:change(d[field]))
            with self.subTest(field=field),self.assertRaises(ValueError):r.read_native_esp(folder,materialized,context)
            for p in folder.iterdir():p.unlink()
            folder.rmdir()

    def test_exact_native_esp_format_and_precision(self):
        folder,materialized,context=self.fixture();esp=r.read_native_esp(folder,materialized,context)
        arrays=esp['arrays'];text=r.native_esp_text(arrays['coords_bohr'],arrays['esp_points_bohr'],arrays['esp_hartree_per_e'])
        r.check_native_esp_roundtrip(text,esp)
        self.assertEqual(len(text.decode().splitlines()[1]),65)
        self.assertEqual(len(text.decode().splitlines()[3]),65)
        with self.assertRaises(ValueError):r.check_native_esp_roundtrip(text.replace(b'1.0000000E-01',b'1.1000000E-01'),esp)

    def test_whole_synthetic_launch_capture_and_mutation(self):
        native,materialized,context=self.fixture();job=self.root/'synthetic-job';job.mkdir();native.rename(job/'native')
        (job/'scratch').mkdir();(job/'materialized').mkdir()
        files={name:b'SYNTHETIC SOFTWARE CONTAINER ONLY\n' for name in r.core.MATERIALIZED_FILES}
        files['exact-esp-input.json']=materialized['files']['exact-esp-input.json']
        for name,blob in files.items():(job/'materialized'/name).write_bytes(blob)
        materialized.update(files=files,hashes=r.core.hashes(files))
        spec=r.expected_launch_spec(materialized,context)
        (job/'launch-spec.json').write_text(json.dumps(spec)+'\n')
        (job/'input.json').write_bytes(files['exact-esp-input.json'])
        (job/'controller-source.py').write_bytes(r.data(r.LAUNCHER,r.LAUNCHER_SHA))
        (job/'qm_worker.py').write_bytes(r.data(context['source']/'qm_worker.py',self.m.WORKER_SHA))
        (job/'controller-native.log').write_text('SYNTHETIC: no launch.\n')
        receipt,_,_=ReceiptTests().fixture()
        receipt.update({k:spec[k] for k in ('core_sha256','materializer_sha256','controller_sha256','worker_sha256',
            'runtime_manifest_sha256','runtime_inventory_sha256','limits','materialized_hashes','input_sha256')})
        pins=r.core.hashes(r.core.snapshot(job/'native',r.core.ESP_NATIVE_FILES))
        receipt.update(spec_sha256=r.sha(job/'launch-spec.json'),native_hashes=pins,result_sha256=pins['result.json'])
        for name in ('controller-progress.json','controller-result.json'):(job/name).write_text(json.dumps(receipt)+'\n')
        captured=r.read_launch_receipt(job,materialized,context);self.assertEqual(captured['native_hashes'],pins)
        spec['limits']['threads']=2.
        (job/'launch-spec.json').write_text(json.dumps(spec)+'\n')
        with self.assertRaisesRegex(ValueError,'Complete launch'):r.read_launch_receipt(job,materialized,context)

    def test_original_charge_gate_function_is_exact(self):
        context={'materializer':self.m,'source':self.m.CACHE/'prepared/sources'}
        function=r.original_charge_validator(context)
        self.assertEqual(function.__name__,'validate_charges')
        constraints={'atom_count':2,'formal_charge':0}
        with self.assertRaises(ValueError):function({'atom_map':[],'formal_charge':0},constraints,constraints,{},[])


class ReceiptTests(unittest.TestCase):
    def fixture(self):
        m=r.core.materializer_module();context={'materializer':m}
        materialized={'hashes':{name:'1'*64 for name in r.core.MATERIALIZED_FILES}}
        spec=r.expected_launch_spec(materialized,context);pins={name:'2'*64 for name in r.core.ESP_NATIVE_FILES}
        receipt={k:spec[k] for k in ('core_sha256','materializer_sha256','controller_sha256','worker_sha256',
            'runtime_manifest_sha256','runtime_inventory_sha256','limits','materialized_hashes','input_sha256')}
        receipt.update(schema_version=1,status='numerically_complete_pending_esp_review',child_pid=None,controller_pid=42,
            physical_acceptance=False,simulation_ready=False,exit_code=0,cleanup_state='verified_no_live_group',
            started_unix=1.,qm_started_unix=2.,completed_unix=4.,updated_unix=5.,
            peak_sampled_process_rss_bytes=10,peak_sampled_group_rss_bytes=20,
            spec_sha256='3'*64,native_hashes=pins,result_sha256=pins['result.json'])
        return receipt,spec,pins

    def check(self,receipt,spec,pins):r.validate_launch_receipt(receipt,copy.deepcopy(receipt),spec,'3'*64,pins)

    def test_success_receipt_format(self):self.check(*self.fixture())

    def test_final_launcher_review_and_enablement_pins(self):
        result=r.launcher_admission();self.assertEqual(result['controller_sha256'],r.LAUNCHER_SHA)

    def test_failed_running_missing_child_and_type_fields_reject(self):
        for key,value in [('schema_version',True),('schema_version',1.),('status','running'),('status','failed'),
                          ('child_pid',123),('physical_acceptance',0),('simulation_ready',0),('exit_code',False),
                          ('cleanup_state','unknown_requires_attention'),('controller_pid',True)]:
            receipt,spec,pins=self.fixture();receipt[key]=value
            with self.subTest(key=key,value=value),self.assertRaises(ValueError):self.check(receipt,spec,pins)
        receipt,spec,pins=self.fixture();del receipt['child_pid']
        with self.assertRaises(ValueError):self.check(receipt,spec,pins)

    def test_source_input_output_and_limits_mutation_reject(self):
        for key in ('core_sha256','materializer_sha256','controller_sha256','worker_sha256','runtime_manifest_sha256',
                    'runtime_inventory_sha256','input_sha256','spec_sha256','result_sha256'):
            receipt,spec,pins=self.fixture();receipt[key]='9'*64
            with self.subTest(key=key),self.assertRaises(ValueError):self.check(receipt,spec,pins)
        receipt,spec,pins=self.fixture();receipt['limits']=dict(receipt['limits'],threads=2.)
        with self.assertRaises(ValueError):self.check(receipt,spec,pins)

    def test_missing_final_progress_time_resources_and_extra_claim_reject(self):
        receipt,spec,pins=self.fixture()
        with self.assertRaises(ValueError):r.validate_launch_receipt(receipt,{},spec,'3'*64,pins)
        for key,value in [('started_unix',True),('qm_started_unix',6.),('completed_unix',float('nan')),
            ('peak_sampled_process_rss_bytes',10.),('peak_sampled_group_rss_bytes',8_000_000_001),('accepted',True)]:
            receipt,spec,pins=self.fixture();receipt[key]=value
            with self.subTest(key=key),self.assertRaises(ValueError):self.check(receipt,spec,pins)

    def test_actual_parent_refusal_precedes_any_downstream_read(self):
        with patch.object(r.core,'revalidate_actual_parent',side_effect=ValueError('unfinished parent')) as parent, \
             patch.object(r,'launcher_admission',side_effect=AssertionError('downstream read')), \
             patch.object(r.core,'verify_materialized_request',side_effect=AssertionError('downstream read')):
            with self.assertRaisesRegex(ValueError,'unfinished parent'):r.read_completed_esp('plan','job')
            parent.assert_called_once_with('plan')

    def test_sample_must_fit_peaks_disk_and_time(self):
        receipt,spec,pins=self.fixture();sample={'maximum_process_rss_bytes':10,'process_group_rss_bytes':20,
            'sample_unix':3.,'free_disk_bytes':r.LAUNCH_LIMITS['running_disk_bytes']};receipt['resource_sample']=sample
        self.check(receipt,spec,pins)
        for key,value in [('maximum_process_rss_bytes',11),('process_group_rss_bytes',False),
                          ('free_disk_bytes',0),('sample_unix',1.)]:
            altered=copy.deepcopy(receipt);altered['resource_sample'][key]=value
            with self.subTest(key=key),self.assertRaises(ValueError):self.check(altered,spec,pins)


class ExistingSyntheticRESPTests(unittest.TestCase):
    """Parse only already-computed synthetic ESP fits, never fit or admit a model."""
    @classmethod
    def setUpClass(cls):
        cls.m=r.core.materializer_module();cls.source=cls.m.CACHE/'prepared/sources'
        cls.plan=r.load(cls.m.CACHE/'prepared/plan-v3-resp-admitted.json',r.core.MATERIALIZER_PLAN_SHA)
        cls.root=Path('/Users/ashujo/.cache/dynamol-runtimes/resp-heap-local-v1/local-qualification')
        manifest=r.load(cls.root/'artifact-manifest.json','f519f2266c1de508994292303192b14c612b9610e079f72d894879fcae7cf189')
        cls.pins={x['path']:x['sha256'] for x in manifest['artifacts']}
        cls.relative='native-parity-v3/synthetic94/'

    def fixture_bytes(self,name):
        relative=self.relative+name;return r.data(self.root/relative,self.pins[relative])

    def stage(self,number):
        name='stage'+str(number)
        files={suffix:self.fixture_bytes(name+'.'+suffix) for suffix in ('out','punch','qout','esout')}
        files['log']=self.fixture_bytes(name+'.stdout')+self.fixture_bytes(name+'.stderr')
        input_=self.fixture_bytes(name+'.in');initial=self.fixture_bytes('canonical.qin' if number==1 else 'stage1.qout')
        nuclei=[int(x.split()[0]) for x in input_.decode().splitlines()[5:99]]
        return files,input_,initial,nuclei,self.fixture_bytes('esp.dat')

    def test_both_existing_synthetic_stages_parse_but_never_admit_execution(self):
        for stage in (1,2):
            parsed=r.parse_native_resp_stage(*self.stage(stage))
            self.assertEqual(len(parsed['charges']),94);self.assertFalse(parsed['execution_provenance_admitted'])
            self.assertFalse(parsed['physical_acceptance'])

    def test_qout_only_open_error_and_truncation_reject(self):
        for mode in ('qout-only','open-error','truncated-punch','truncated-esout','nonfinite-qout','punch-mismatch'):
            files,input_,initial,nuclei,esp=self.stage(2)
            if mode=='qout-only':files={'qout':files['qout']}
            elif mode=='open-error':files['log']+=b'Unit 5 Error on OPEN\n'
            elif mode=='truncated-punch':files['punch']=files['punch'][:1200]
            elif mode=='truncated-esout':files['esout']=files['esout'][:1200]
            elif mode=='nonfinite-qout':files['qout']=files['qout'].replace(b'-0.366200',b'      nan',1)
            else:files['punch']=files['punch'].replace(b'-0.366200',b'-0.466200',1)
            with self.subTest(mode=mode),self.assertRaises((ValueError,IndexError)):r.parse_native_resp_stage(files,input_,initial,nuclei,esp)

    def test_native_identity_initial_charge_and_statistics_mutation_reject(self):
        for mode in ('nuclei','initial','statistics','ivary','restraint-weight'):
            files,input_,initial,nuclei,esp=self.stage(2)
            if mode=='nuclei':nuclei[0]=8
            elif mode=='initial':initial=initial.replace(b'-0.366200',b'-0.466200',1)
            elif mode=='restraint-weight':files['punch']=files['punch'].replace(b'0.001000',b'0.002000',1)
            elif mode=='ivary':input_=input_.replace(b'    6   -1',b'    6    0',1)
            else:files['out']=files['out'].replace(b'Statistics of the fitting:',b'INCOMPLETE').split(b'The residual')[0]
            with self.subTest(mode=mode),self.assertRaises(ValueError):r.parse_native_resp_stage(files,input_,initial,nuclei,esp)

    def test_joint_candidate_keeps_all_values_but_is_explicitly_unadmitted(self):
        source=r.load(self.source/'cap/capped-adduct.json');constraints=r.load(self.source/'cap/resp/constraints.json')
        context={'source':self.source,'materializer':self.m,'plan':self.plan,
                 'request':{'atom_ids':[f'{x["index"]:03d}:{x["name"]}' for x in source['atom_map']]}}
        raw=self.fixture_bytes('esp.dat');lines=raw.decode().splitlines();n=int(lines[0].split()[0])
        xyz=np.asarray([[float(x) for x in line.split()] for line in lines[1:1+n]])
        points=np.asarray([[float(x) for x in line.split()] for line in lines[1+n:]])
        esp={'arrays':{'coords_bohr':xyz,'esp_points_bohr':points[:,1:],'esp_hartree_per_e':points[:,0]}}
        stage1=r.parse_native_resp_stage(*self.stage(1));stage2=r.parse_native_resp_stage(*self.stage(2))
        candidate=r.validate_joint_charge_candidate(context,esp,stage1,stage2,raw,source,constraints)
        self.assertFalse(candidate['accepted']);self.assertFalse(candidate['execution_provenance_admitted'])
        self.assertFalse(candidate['physical_acceptance']);self.assertEqual(len(candidate['retained_atoms']),82)
        for row in candidate['retained_atoms']:
            i=context['request']['atom_ids'].index(row['atom_id']);self.assertEqual(row['charge_e'],stage2['charges'][i])
        for key in ('esp_sha256','initial_charge_sha256','input_sha256'):
            changed=copy.deepcopy(stage2);changed[key]='0'*64
            with self.subTest(key=key),self.assertRaises(ValueError):
                r.validate_joint_charge_candidate(context,esp,stage1,changed,raw,source,constraints)

    def test_same_native_synthetic_charges_retain_exact_values_and_neutral_cap(self):
        # Explicitly synthetic input to the pure charge gate, not an actual parent/fit receipt.
        source=r.load(self.source/'cap/capped-adduct.json');constraints=r.load(self.source/'cap/resp/constraints.json')
        context={'source':self.source,'materializer':self.m}
        parsed=r.parse_native_resp_stage(*self.stage(2));charges=parsed['charges']
        values=r.original_charge_validator(context)(source,constraints,constraints,
            {'charges':charges,'constraint_checks_passed':True},charges)
        self.assertEqual(values['charges_e'],charges);self.assertEqual(values['canonical_fixed_count'],18)
        self.assertEqual(len(values['removed_cap_indices']),12)
        self.assertLess(abs(values['retained_charge_sum_e']),1e-4)
        for index in (constraints['frozen_atoms'][0]['index'],constraints['stage2_methyl_equivalence_groups'][0][1]):
            changed=charges.copy();changed[index]+=.01
            with self.assertRaises(ValueError):r.original_charge_validator(context)(source,constraints,constraints,
                {'charges':changed,'constraint_checks_passed':True},changed)


if __name__=='__main__':unittest.main(verbosity=2)

