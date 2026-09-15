"""Synthetic candidate callback integration; stdlib only, no NumPy/PySCF/QM."""
import ast
import copy
import hashlib
import importlib.util
import json
import math
from pathlib import Path
import tempfile
from types import SimpleNamespace
import unittest
from unittest.mock import patch

HERE=Path(__file__).resolve().parent
spec=importlib.util.spec_from_file_location('candidate_worker',HERE/'qm_worker.py')
worker=importlib.util.module_from_spec(spec);spec.loader.exec_module(worker)


class FakeMolecule:
    def __init__(self,coords):
        self.xyz=copy.deepcopy(coords)
        self.symbols=['O','H','H'];self.natm=3;self.charge=0;self.spin=0
        self._basis={'O':[[0,[1.0,1.0]]],'H':[[0,[1.0,1.0]]]}
        self._ecp={};self.cart=False;self.symmetry=False

    def atom_coords(self,unit):
        assert unit=='Bohr'
        return copy.deepcopy(self.xyz)

    def atom_symbol(self,index):return self.symbols[index]


class FakeSCF:
    def __init__(self,mol):
        self.mol=mol;self.with_df=SimpleNamespace(auxbasis='SYNTHETIC-AUX')
        self.max_cycle=100;self.conv_tol=1e-10;self.conv_tol_grad=1e-7;self.max_memory=128
        self.e_tot=-123.456;self.converged=True


class FakeNumpy:
    @staticmethod
    def asarray(x):return x

    @staticmethod
    def isfinite(x):return SimpleNamespace(all=lambda:all(math.isfinite(v) for row in x for v in row))

    linalg=SimpleNamespace(norm=lambda x:math.sqrt(sum(v*v for row in x for v in row)))


def existing_callback(record,data,report,progress):
    """Execute the actual nested candidate callback, replacing only NumPy I/O."""
    tree=ast.parse((HERE/'qm_worker.py').read_text())
    calculate=next(x for x in tree.body if isinstance(x,ast.FunctionDef) and x.name=='calculate')
    callback=next(x for x in ast.walk(calculate) if isinstance(x,ast.FunctionDef) and x.name=='callback')
    namespace={'record_evaluation':record,'np':FakeNumpy,'math':math,'dihedrals':[],
        'report':report,'progress':progress,'opt':data['optimization']}
    exec(compile(ast.fix_missing_locations(ast.Module(body=[callback],type_ignores=[])),
                 str(HERE/'qm_worker.py'),'exec'),namespace)
    return namespace['callback']


class IntegrationTests(unittest.TestCase):
    def setUp(self):
        self.temp=tempfile.TemporaryDirectory();self.addCleanup(self.temp.cleanup)
        self.out=Path(self.temp.name)
        self.bohr=.52917721092
        self.initial=[[0.,0.,0.],[2.,0.,0.],[0.,2.,0.]]
        self.input={'schema_version':1,'atom_ids':['SYNTHETIC:O','SYNTHETIC:H1','SYNTHETIC:H2'],
            'elements':['O','H','H'],'coords_angstrom':[[x*self.bohr for x in row] for row in self.initial],
            'charge':0,'spin':0,'method':'RHF','basis':'SYNTHETIC-NOT-EVALUATED',
            'density_fit':True,'auxbasis':'SYNTHETIC-AUX','operations':['energy','gradient'],
            'threads':1,'max_memory_mb':128,'max_wall_seconds':60,
            'optimization':{'enabled':True}}
        self.data=worker.validate_input(self.input)
        self.report={'input_sha256':hashlib.sha256(json.dumps(self.input).encode()).hexdigest(),
            'worker_sha256':hashlib.sha256((HERE/'qm_worker.py').read_bytes()).hexdigest(),
            'accepted':False,'optimization':{'steps':[]}}
        self.mean_field=FakeSCF(FakeMolecule(self.initial))
        self.scanner=SimpleNamespace(base=copy.deepcopy(self.mean_field),unit='au',
            de=[[.1,-.2,.3],[-.1,.2,-.3],[0.,0.,0.]],converged=True)
        self.scanner.mol=self.scanner.base.mol
        self.env={'mol':copy.deepcopy(self.scanner.mol),'g_scanner':self.scanner,
            'self':SimpleNamespace(cycle=1),'coords':copy.deepcopy(self.initial),
            'energy':self.scanner.base.e_tot,'gradients':copy.deepcopy(self.scanner.de)}
        self.record=worker._create_evaluation_journal(self.data,self.data['elements'],self.initial,
            self.mean_field,self.out,self.report,self.bohr)

    def read_record(self):
        return json.loads((self.out/'evaluation-journal/evaluation-000001.json').read_bytes())['payload']

    def manifest(self):
        return json.loads((self.out/'evaluation-journal/manifest.json').read_bytes())['payload']

    def closure(self,name):
        return dict(zip(self.record.__code__.co_freevars,[x.cell_contents for x in self.record.__closure__]))[name]

    def test_actual_callback_captures_current_geometry_gradient_energy_and_progress_reference(self):
        moved=[[.01,.02,.03],[2.1,0.,0.],[0.,1.9,0.]]
        self.env['coords']=copy.deepcopy(moved)
        self.env['mol'].xyz=copy.deepcopy(moved);self.scanner.mol.xyz=copy.deepcopy(moved)
        events=[]
        callback=existing_callback(self.record,self.data,self.report,lambda *args,**kw:events.append((args,kw)))
        callback(self.env)
        record=self.read_record()
        self.assertEqual(record['coords_bohr'],moved)
        self.assertNotEqual(record['coords_bohr'],self.initial)
        self.assertEqual(record['gradient_hartree_per_bohr'],self.env['gradients'])
        self.assertEqual(record['energy_hartree'],self.env['energy'])
        self.assertEqual(record['identity']['atom_ids'],self.input['atom_ids'])
        self.assertEqual(record['units'],{'coordinates':'bohr','gradient':'hartree/bohr','energy':'hartree'})
        self.assertEqual(record['status'],'UNCONVERGED');self.assertFalse(record['accepted'])
        self.assertFalse(record['checkpoint_reuse_authorized']);self.assertFalse(self.report['accepted'])
        step=self.report['optimization']['steps'][0]
        reference=step['evaluation_record']
        self.assertEqual(reference['file'],'evaluation-journal/evaluation-000001.json')
        self.assertEqual(reference['file_sha256'],hashlib.sha256((self.out/reference['file']).read_bytes()).hexdigest())
        self.assertEqual(events[0][1]['evaluation_record'],reference)

    def test_input_conversion_and_current_method_fingerprints_are_bound(self):
        self.record(self.env)
        manifest=self.manifest();context=manifest['callback_contract'];record=self.read_record()
        self.assertEqual(context['input_coordinate_unit'],'angstrom')
        self.assertEqual(context['bohr_to_angstrom'],self.bohr)
        self.assertEqual(context['conversion_source'],'pyscf.lib.param.BOHR')
        self.assertEqual(context['original_resolved_input_sha256'],worker._hash(self.data))
        self.assertEqual(context['native_initial_geometry_bohr_sha256'],worker._hash(self.initial))
        self.assertEqual(context['native_initial_method']['basis_sha256'],worker._hash(self.scanner.mol._basis))
        self.assertEqual(context['native_initial_method']['auxbasis'],'SYNTHETIC-AUX')
        self.assertEqual(context['native_initial_method_sha256'],worker._hash(context['native_initial_method']))
        self.assertEqual(record['requested_method_fingerprint_sha256'],manifest['requested_method_fingerprint_sha256'])
        self.assertEqual(record['manifest_sha256'],worker._hash(manifest))

    def test_requested_bohr_input_is_bound_without_conversion(self):
        data=copy.deepcopy(self.data);data['coords_bohr']=data.pop('coords_angstrom')
        data['coords_bohr']=copy.deepcopy(self.initial)
        output=self.out/'bohr';output.mkdir()
        record=worker._create_evaluation_journal(data,data['elements'],self.initial,
            self.mean_field,output,self.report,self.bohr)
        record(self.env)
        manifest=json.loads((output/'evaluation-journal/manifest.json').read_bytes())['payload']
        self.assertEqual(manifest['callback_contract']['input_coordinate_unit'],'bohr')
        self.assertEqual(manifest['resolved_request_sha256'],worker._hash(data))

    def test_native_element_normalization_has_both_source_and_native_fingerprints(self):
        data=copy.deepcopy(self.data);data['elements']=['o','h','h']
        output=self.out/'normalized';output.mkdir()
        record=worker._create_evaluation_journal(data,['O','H','H'],self.initial,
            self.mean_field,output,self.report,self.bohr)
        record(self.env)
        manifest=json.loads((output/'evaluation-journal/manifest.json').read_bytes())['payload']
        self.assertEqual(manifest['callback_contract']['original_input_elements'],['o','h','h'])
        self.assertEqual(manifest['callback_contract']['canonical_native_elements'],['O','H','H'])
        self.assertEqual(manifest['callback_contract']['original_resolved_input_sha256'],worker._hash(data))
        self.assertNotEqual(manifest['resolved_request_sha256'],worker._hash(data))

    def test_false_scf_flag_skips_journal_and_preserves_existing_callback_status(self):
        self.scanner.converged=False;self.scanner.base.converged=False
        callback=existing_callback(self.record,self.data,self.report,lambda *args,**kw:None)
        callback(self.env)
        self.assertFalse(list((self.out/'evaluation-journal').glob('evaluation-*.json')))
        self.assertFalse(self.report['optimization']['steps'][0]['scf_converged'])
        self.assertNotIn('evaluation_record',self.report['optimization']['steps'][0])
        self.assertFalse(self.report['accepted'])

    def test_truthy_nonboolean_scf_flag_is_refused(self):
        for value in (1,'true',None):
            with self.subTest(value=value),self.assertRaises(ValueError):
                self.scanner.converged=value;self.record(self.env)

    def test_disagreeing_scanner_base_convergence_is_refused(self):
        self.scanner.base.converged=False
        with self.assertRaises(ValueError):self.record(self.env)

    def test_stale_geometry_in_each_native_molecule_is_refused(self):
        for name in ('engine','gradient_scanner','scf_scanner'):
            with self.subTest(name=name):
                env=copy.deepcopy(self.env)
                if name=='engine':env['mol'].xyz[0][0]=.1
                elif name=='gradient_scanner':
                    env['g_scanner'].mol=copy.deepcopy(env['g_scanner'].mol)
                    env['g_scanner'].mol.xyz[0][0]=.1
                else:env['g_scanner'].base.mol.xyz[0][0]=.1
                with self.assertRaisesRegex(ValueError,'geometry differs'):self.record(env)

    def test_stale_energy_and_gradient_are_refused(self):
        self.env['energy']+=.1
        with self.assertRaisesRegex(ValueError,'energy differs'):self.record(self.env)
        self.env['energy']=self.scanner.base.e_tot
        self.env['gradients'][0][0]+=.1
        with self.assertRaisesRegex(ValueError,'gradient differs'):self.record(self.env)

    def test_changed_method_settings_basis_and_auxbasis_are_refused(self):
        changes=[lambda m:setattr(m,'conv_tol',1e-6),lambda m:setattr(m,'max_cycle',50),
            lambda m:setattr(m,'xc','OTHER'),lambda m:setattr(m,'with_df',None),
            lambda m:setattr(m.with_df,'auxbasis','OTHER'),
            lambda m:m.mol._basis.update(O=[[0,[2.,1.]]])]
        for change in changes:
            env=copy.deepcopy(self.env);change(env['g_scanner'].base)
            with self.subTest(change=change),self.assertRaisesRegex(ValueError,'method or resolved'):
                self.record(env)

    def test_changed_native_method_family_is_refused(self):
        self.scanner.base=SimpleNamespace(**self.scanner.base.__dict__)
        with self.assertRaisesRegex(ValueError,'method family'):self.record(self.env)

    def test_current_rks_functional_and_grid_are_bound(self):
        data=copy.deepcopy(self.data);data.update(method='RKS',functional='B3LYPG',grid_level=4)
        mean_field=copy.deepcopy(self.mean_field);mean_field.xc='B3LYPG';mean_field.grids=SimpleNamespace(level=4)
        output=self.out/'rks';output.mkdir()
        record=worker._create_evaluation_journal(data,data['elements'],self.initial,mean_field,
            output,self.report,self.bohr)
        env=copy.deepcopy(self.env);env['g_scanner'].base=copy.deepcopy(mean_field)
        record(env)
        manifest=json.loads((output/'evaluation-journal/manifest.json').read_bytes())['payload']
        self.assertEqual(manifest['callback_contract']['native_initial_method']['xc'],'B3LYPG')
        self.assertEqual(manifest['callback_contract']['native_initial_method']['grid_level'],4)
        env['self'].cycle=2;env['g_scanner'].base.grids.level=3
        with self.assertRaisesRegex(ValueError,'method or resolved'):record(env)

    def test_native_identity_state_and_unit_mismatch_are_refused(self):
        for name in ('element','charge','spin','unit'):
            env=copy.deepcopy(self.env)
            if name=='element':env['mol'].symbols=['H','O','H']
            elif name=='charge':env['mol'].charge=False
            elif name=='spin':env['mol'].spin=0.0
            else:env['g_scanner'].unit='Eh/Ang'
            with self.subTest(name=name),self.assertRaises(ValueError):self.record(env)

    def test_gradient_atom_subset_or_reordering_is_refused(self):
        for indices in ([0,1],[0,2,1],[False,1,2],[0,1.0,2]):
            self.scanner.atmlst=indices
            with self.subTest(indices=indices),self.assertRaisesRegex(ValueError,'gradient atom indices'):
                self.record(self.env)
        self.scanner.atmlst=[0,1,2]
        self.record(self.env)

    def test_nonfinite_callback_values_are_refused(self):
        for name in ('coords','energy','gradients'):
            env=copy.deepcopy(self.env)
            if name=='energy':env[name]=float('nan')
            else:env[name][0][0]=float('nan')
            with self.subTest(name=name),self.assertRaises(ValueError):self.record(env)

    def test_mutable_request_cannot_change_frozen_declared_method(self):
        self.data['basis']='MUTATED-AFTER-CREATE'
        self.record(self.env)
        self.assertEqual(self.manifest()['requested_method_and_settings']['basis'],'SYNTHETIC-NOT-EVALUATED')

    def test_persistence_failure_propagates_from_actual_callback_without_progress_claim(self):
        callback=existing_callback(self.record,self.data,self.report,lambda *args,**kw:self.fail('progress emitted'))
        with patch.object(worker.os,'link',side_effect=OSError('synthetic disk failure')):
            with self.assertRaisesRegex(OSError,'disk failure'):callback(self.env)
        self.assertEqual(self.report['optimization']['steps'],[])
        self.assertFalse(self.report['accepted'])
        self.assertFalse(list((self.out/'evaluation-journal').glob('evaluation-*.json')))

    def test_failure_after_publication_leaves_only_unaccepted_recovery_record(self):
        module=self.closure('journal_module')
        callback=existing_callback(self.record,self.data,self.report,lambda *args,**kw:self.fail('progress emitted'))
        with patch.object(module,'_sync_directory',side_effect=OSError('synthetic sync failure')):
            with self.assertRaisesRegex(OSError,'sync failure'):callback(self.env)
        self.assertEqual(self.report['optimization']['steps'],[])
        self.assertFalse(self.read_record()['accepted'])
        self.assertEqual(self.read_record()['status'],'UNCONVERGED')

    def test_worker_run_marks_persistence_failure_failed_and_removes_result_arrays(self):
        source=self.out/'input.json';source.write_text(json.dumps(self.input))
        output=self.out/'failed-worker'

        def calculate(data,out,report):
            report['phase']='optimization';report['optimization']={'steps':[]}
            record=worker._create_evaluation_journal(data,data['elements'],self.initial,
                self.mean_field,out,report,self.bohr)
            (out/'arrays.npz').write_bytes(b'SYNTHETIC-PENDING-NOT-A-NUMPY-ARRAY')
            callback=existing_callback(record,data,report,lambda *args,**kw:None)
            with patch.object(worker.os,'link',side_effect=OSError('synthetic journal publication failure')):
                callback(self.env)

        with patch.object(worker,'calculate',calculate):result=worker.run(source,output)
        self.assertEqual(result['status'],'failed');self.assertFalse(result['accepted'])
        self.assertIn('journal publication failure',result['error']['message'])
        self.assertFalse((output/'arrays.npz').exists())
        self.assertEqual(json.loads((output/'progress.json').read_text())['status'],'failed')
        self.assertTrue(list((output/'evaluation-journal').glob('.evaluation-*.tmp')))

    def test_journal_hook_is_before_progress_and_native_assertion_remains_enabled(self):
        tree=ast.parse((HERE/'qm_worker.py').read_text())
        calculate=next(x for x in tree.body if isinstance(x,ast.FunctionDef) and x.name=='calculate')
        callback=next(x for x in ast.walk(calculate) if isinstance(x,ast.FunctionDef) and x.name=='callback')
        self.assertEqual(ast.unparse(callback.body[0]),'evaluation_record = record_evaluation(env)')
        calls=[x for x in ast.walk(calculate) if isinstance(x,ast.Call)
               and ast.unparse(x.func)=='geometric_solver.kernel']
        self.assertEqual(len(calls),1)
        self.assertTrue(next(k.value.value for k in calls[0].keywords if k.arg=='assert_convergence'))


if __name__=='__main__':unittest.main()
