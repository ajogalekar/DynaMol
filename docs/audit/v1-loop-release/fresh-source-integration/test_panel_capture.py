"""Pure-Python contract checks; these never execute preparation or native tools."""
import ast
import importlib.util
import json
import os
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch

HERE = Path(__file__).parent

def module():
    spec = importlib.util.spec_from_file_location('panel_capture_contract', HERE / 'panel_capture.py')
    value = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(value)
    return value


def node(path, name):
    return next(n for n in ast.parse(path.read_text()).body if isinstance(n, (ast.FunctionDef, ast.ClassDef)) and n.name == name)


class PanelCaptureContract(unittest.TestCase):
    def setUp(self):
        self.directory = tempfile.TemporaryDirectory()
        self.addCleanup(self.directory.cleanup)
        self.root = Path(self.directory.name)
        self.m = module()
        self.value = {'case': 'TEST-2', 'root': str(self.root), 'dataset_id': 'test_raw',
                      'expected_raw_atoms': 19, 'seed': 47}

    def configure(self, **changes):
        path = self.root / 'config.json'
        path.write_text(json.dumps({**self.value, **changes}))
        self.m.configure(path)

    def test_explicit_case_values_reach_job_and_frozen_source_paths(self):
        self.configure()
        self.assertEqual((self.m.CASE, self.m.DATASET, self.m.EXPECTED_RAW_ATOMS, self.m.SEED), ('TEST-2', 'test_raw', 19, 47))
        self.assertEqual(self.m.JOB, 'fresh-test-2-capture')
        self.assertEqual(self.m.REPO, self.root / 'implementation/source')

    def test_environment_configuration(self):
        with patch.dict(os.environ, {'DYNAMOL_CAPTURE_CASE': 'Q2', 'DYNAMOL_CAPTURE_ROOT': str(self.root),
                                     'DYNAMOL_CAPTURE_DATASET': 'd2', 'DYNAMOL_CAPTURE_RAW_ATOMS': '7',
                                     'DYNAMOL_CAPTURE_SEED': '9'}, clear=True):
            self.m.configure()
        self.assertEqual((self.m.DATASET, self.m.EXPECTED_RAW_ATOMS, self.m.SEED), ('d2', 7, 9))

    def test_invalid_configuration_rejected_before_native_entry(self):
        for changes in [{'expected_raw_atoms': True}, {'expected_raw_atoms': 0}, {'seed': -1},
                        {'seed': 2**32}, {'dataset_id': '../escape'}, {'case': 'A/B'},
                        {'root': 'relative'}, {'unexpected': True}]:
            with self.subTest(changes=changes), self.assertRaises(ValueError):
                self.configure(**changes)

    def test_all_raw_evidence_and_actual_selected_sequence_are_bound(self):
        self.configure()
        raw = self.root / 'workspace/datasets/test_raw'
        records = []
        for name in ['topology.pdb', 'sequence-source.pdb', 'original-source/first.cif', 'original-source/second.pdb']:
            path = raw / name
            path.parent.mkdir(parents=True, exist_ok=True)
            path.write_text(name)
            records.append({'source': str(path), 'copy': str(path), 'sha256': self.m.digest(path)})
        (self.root / 'raw-source-provenance.json').write_text(json.dumps({'files': records}))
        files = self.m.raw_snapshot_files(raw, raw / 'sequence-source.pdb')
        self.assertEqual(set(files), {'raw-dataset/' + Path(r['copy']).relative_to(raw).as_posix() for r in records} | {'sequence-evidence.pdb'})
        self.assertEqual(files['sequence-evidence.pdb'], raw / 'sequence-source.pdb')
        with self.assertRaisesRegex(ValueError, 'actual sequence evidence'):
            self.m.raw_snapshot_files(raw, None)
        with self.assertRaisesRegex(ValueError, 'one of the pinned'):
            self.m.raw_snapshot_files(raw, self.root / 'external.cif')
        (raw / 'unlisted.cif').write_text('unlisted')
        with self.assertRaisesRegex(ValueError, 'enumerate every'):
            self.m.raw_snapshot_files(raw, raw / 'sequence-source.pdb')

    def test_frozen_source_copies_include_bundled_data_and_stay_independent(self):
        repo = self.root / 'repo'
        for name in ['backend/__init__.py', 'backend/example.py', 'backend/data/modified_residues/TPO.cif', 'pyproject.toml']:
            path = repo / name
            path.parent.mkdir(parents=True, exist_ok=True)
            path.write_text(name)
        self.configure(source_repo=str(repo))
        self.m.freeze_sources()
        records = json.loads((self.root / 'frozen-source-provenance.json').read_text())['files']
        self.assertEqual(len(records), 4)
        self.m.verify(records, 'copy')
        original = repo / 'backend/example.py'
        original.write_text('later repository edit')
        self.assertEqual((self.m.REPO / 'backend/example.py').read_text(), 'backend/example.py')
        with self.assertRaisesRegex(ValueError, 'new attempt root'):
            self.m.freeze_sources()

    def test_owned_tree_and_supervisor_lifecycle_ast_unchanged(self):
        original, panel = HERE / 'capture_preparation.py', HERE / 'panel_capture.py'
        self.assertEqual(ast.dump(node(original, 'OwnedTree')), ast.dump(node(panel, 'OwnedTree')))
        def lifecycle(path):
            body = node(path, 'launch').body
            start = next(i for i, item in enumerate(body) if isinstance(item, ast.Assign) and any(isinstance(t, ast.Name) and t.id == 'tree' for t in item.targets))
            return ast.dump(ast.Module(body=body[start:], type_ignores=[]))
        self.assertEqual(lifecycle(original), lifecycle(panel))
        for name in ['parameter_map', 'source_mapping', 'backbone_references']:
            self.assertEqual(ast.dump(node(original, name)), ast.dump(node(panel, name)))

    def test_refinement_and_publication_interception_ast_unchanged(self):
        def capture_class(path):
            return next(n for n in ast.walk(node(path, 'child')) if isinstance(n, ast.ClassDef) and n.name == 'CaptureWorker')
        self.assertEqual(ast.dump(capture_class(HERE / 'capture_preparation.py')), ast.dump(capture_class(HERE / 'panel_capture.py')))


if __name__ == '__main__':
    unittest.main()
