"""Tiny synthetic file-contract tests; never a genuine parent result or native run."""
import copy
import importlib.util
import io
import json
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch
import numpy as np

P=Path(__file__).with_name('recovery_cold_evidence_core.py')
s=importlib.util.spec_from_file_location('cold_core',P);c=importlib.util.module_from_spec(s);s.loader.exec_module(c)
T=Path(__file__).with_name('test_materialize_6oim_recovery_esp.py')
s=importlib.util.spec_from_file_location('synthetic_contracts',T);fixture=importlib.util.module_from_spec(s);s.loader.exec_module(fixture)


class CoreTests(unittest.TestCase):
    def setUp(self):
        self.m=c.materializer_module();scratch=self.m.CACHE/'software-test-scratch';scratch.mkdir(exist_ok=True)
        self.temp=tempfile.TemporaryDirectory(dir=scratch);self.addCleanup(self.temp.cleanup);self.root=Path(self.temp.name)

    def context(self):
        esp,_,parent,method,coords,_,_,_,_=fixture.synthetic()
        source=self.m.CACHE/'prepared/sources';manifest=c.load(source/'source-manifest.json')
        return dict(materializer=self.m,source=source,request=parent,method=method,arrays={'coords_bohr':np.array(coords)},
            esp_request=esp,before={'result.json':'a'*64},controller={'synthetic_controller':True},
            geometry={'scope':'invented two-atom serialization fixture, not accepted biology','passed':False},
            qualification={'inventory_sha256':'synthetic-not-runtime-admission'},resp_provider={'provider_admitted':False},
            plan={'preparation_runtime_admission':{'root_reviewed':False,'synthetic':True},
                  'source_pins':{k:v['sha256'] for k,v in manifest.items()}})

    def materialized(self):
        context=self.context();folder=self.root/'synthetic';folder.mkdir()
        for name,value in c.expected_materialized_payloads(context).items():(folder/name).write_text(json.dumps(value)+'\n')
        result=c.expected_materialized_report(folder,context)
        (folder/'result.json').write_text(json.dumps(result)+'\n');return folder,context

    def mutate(self,folder,file,function):
        obj=c.load(folder/file);function(obj);(folder/file).write_text(json.dumps(obj)+'\n')

    def test_synthetic_complete_payload_binding(self):
        folder,context=self.materialized();result=c.verify_materialized_request(folder,context)
        self.assertEqual(set(result['hashes']),c.MATERIALIZED_FILES)

    def test_full_downstream_contract_mutations_reject(self):
        for key,value in [('accepted',True),('accepted',0),('fit_inputs',{}),('resp_provider',{'provider_admitted':True}),
                          ('cold_branch_handoff_implemented',True)]:
            folder,context=self.materialized()
            self.mutate(folder,'downstream-contract.json',lambda d:d.update({key:value}))
            with self.subTest(key=key),self.assertRaisesRegex(ValueError,'Complete materialized payload'):c.verify_materialized_request(folder,context)
            for p in folder.iterdir():p.unlink()
            folder.rmdir()

    def test_full_parent_binding_mutations_reject(self):
        mutations=[('source_manifest_sha256','wrong'),('preparation_runtime_admission',{}),
                   ('preparation_runtime_inventory_sha256','wrong'),('accepted',0),('checkpoint_geometry_validated',0),
                   ('extra_unreviewed_key',True)]
        for key,value in mutations:
            folder,context=self.materialized();self.mutate(folder,'parent-binding.json',lambda d:d.update({key:value}))
            with self.subTest(key=key),self.assertRaisesRegex(ValueError,'Complete materialized payload'):c.verify_materialized_request(folder,context)
            for p in folder.iterdir():p.unlink()
            folder.rmdir()

    def test_missing_terminal_success_refuses_even_if_input_exists(self):
        folder,context=self.materialized();(folder/'result.json').unlink()
        with self.assertRaisesRegex(ValueError,'inventory'):c.verify_materialized_request(folder,context)

    def test_failed_terminal_status_refuses(self):
        folder,context=self.materialized();self.mutate(folder,'result.json',lambda d:d.update(status='refused',request_emitted=False))
        with self.assertRaises(ValueError):c.verify_materialized_request(folder,context)

    def test_request_changed_with_updated_local_input_hash_still_refuses(self):
        folder,context=self.materialized();self.mutate(folder,'exact-esp-input.json',lambda d:d['coords_bohr'][0].__setitem__(0,False))
        self.mutate(folder,'result.json',lambda d:d.update(input_sha256=c.sha(folder/'exact-esp-input.json')))
        with self.assertRaisesRegex(ValueError,'Complete materialized payload'):c.verify_materialized_request(folder,context)

    def test_geometry_extra_field_refuses(self):
        folder,context=self.materialized();self.mutate(folder,'geometry.json',lambda d:d.update(accepted=True))
        with self.assertRaises(ValueError):c.verify_materialized_request(folder,context)

    def test_source_plan_hash_change_refuses(self):
        folder,context=self.materialized();self.mutate(folder,'result.json',lambda d:d.update(plan_sha256='0'*64))
        with self.assertRaises(ValueError):c.verify_materialized_request(folder,context)

    def test_complete_terminal_extra_or_contradictory_remaining_claim_refuses(self):
        for key,value in [('extra_acceptance_claim',True),('remaining',[])]:
            folder,context=self.materialized();self.mutate(folder,'result.json',lambda d:d.update({key:value}))
            with self.subTest(key=key),self.assertRaisesRegex(ValueError,'Complete materializer terminal report'):c.verify_materialized_request(folder,context)
            for p in folder.iterdir():p.unlink()
            folder.rmdir()

    def test_extra_file_and_symlink_refuse(self):
        folder,context=self.materialized();(folder/'untracked').write_text('x')
        with self.assertRaises(ValueError):c.verify_materialized_request(folder,context)
        (folder/'untracked').unlink();original=(folder/'geometry.json').read_bytes();(folder/'geometry.json').unlink()
        outside=self.root/'geometry-target';outside.write_bytes(original);(folder/'geometry.json').symlink_to(outside)
        with self.assertRaises(ValueError):c.verify_materialized_request(folder,context)

    def test_missing_parent_refuses_without_writes_or_downstream_calls(self):
        parent=self.root/'missing-parent'
        self.m.PARENT=parent
        read_parent=c.readonly_parent_function(self.m)
        with self.assertRaisesRegex(ValueError,'terminal parent absent'):
            read_parent({'parent':str(parent),'checkpoint_policy':self.m.COLD})
        self.assertFalse(parent.exists())
        self.assertNotIn('publish',read_parent.__code__.co_names)

    def test_source_hash_change_refuses_before_compilation(self):
        changed=self.root/'source.py';changed.write_text('raise AssertionError("must not execute")')
        self.m.__file__=str(changed)
        with self.assertRaisesRegex(ValueError,'digest mismatch'):c.readonly_parent_function(self.m)

    def test_small_real_array_container(self):
        buf=io.BytesIO();np.savez(buf,x=np.array([[1.,2.,3.]]));result=c.real_arrays(buf.getvalue(),{'x':(1,3)})
        self.assertEqual(result['x'].tolist(),[[1.,2.,3.]])

    def test_arrays_reject_complex_boolean_nonfinite_and_wrong_shape(self):
        for value in [np.array([[1+0j,2,3]]),np.array([[True,False,True]]),np.array([[np.nan,2,3]]),np.zeros((2,3))]:
            buf=io.BytesIO();np.savez(buf,x=value)
            with self.subTest(dtype=str(value.dtype),shape=value.shape),self.assertRaises(ValueError):c.real_arrays(buf.getvalue(),{'x':(1,3)})

    def test_arrays_reject_extra_entries_and_bounded_inflation(self):
        buf=io.BytesIO();np.savez(buf,x=np.ones((1,3)),y=np.ones((1,3)))
        with self.assertRaises(ValueError):c.real_arrays(buf.getvalue(),{'x':(1,3)})
        with patch.object(c,'MAX_ARTIFACT_BYTES',1):
            with self.assertRaisesRegex(ValueError,'Oversized'):c.real_arrays(buf.getvalue(),{'x':(1,3),'y':(1,3)})


if __name__=='__main__':unittest.main(verbosity=2)
