"""Native-file-format software fixtures only; no SCF, ESP, RESP or biological result."""
import copy
import importlib.util
import json
from pathlib import Path
import tempfile
import unittest
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
            resources={'threads':2,'max_memory_mb':2000,'max_wall_seconds':14400,'requested_threads':2,'actual_pyscf_threads':2},
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

    def test_original_charge_gate_function_is_exact(self):
        context={'materializer':self.m,'source':self.m.CACHE/'prepared/sources'}
        function=r.original_charge_validator(context)
        self.assertEqual(function.__name__,'validate_charges')
        with self.assertRaises(ValueError):function({'atom_map':[]},{},{},{},[])


if __name__=='__main__':unittest.main(verbosity=2)
