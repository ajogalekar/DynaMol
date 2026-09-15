"""Synthetic integrity/lifecycle checks; no molecular or quantum result is made."""
import hashlib
import importlib.util
import json
import os
from pathlib import Path
import shutil
import sys
import tempfile
import time
import types
import unittest
from unittest.mock import patch

SOURCE = Path(__file__).with_name('launch_cold_esp.py')
spec = importlib.util.spec_from_file_location('esp_launcher_candidate', SOURCE)
launch = importlib.util.module_from_spec(spec)
spec.loader.exec_module(launch)


class LaunchChecks(unittest.TestCase):
    def setUp(self):
        self.temporary = tempfile.TemporaryDirectory(prefix='dynamol-launch-contract-')
        self.folder = Path(self.temporary.name).resolve()

    def tearDown(self):
        self.temporary.cleanup()

    def test_unadmitted_cli_never_loads_or_spawns(self):
        with patch.object(launch, 'CACHE', self.folder), patch.object(launch, 'LAUNCH_ADMITTED', False), patch.object(launch, 'module', side_effect=AssertionError('must not load')) as loader:
            result = launch.run(self.folder/'absent-materialization', self.folder/'job')
            loader.assert_not_called()
        self.assertEqual(result['status'], 'failed')
        self.assertIsNone(result['child_pid'])
        self.assertFalse(result['physical_acceptance'])
        self.assertFalse(result['simulation_ready'])
        self.assertEqual(result['cleanup_state'], 'no_child_launched')
        self.assertEqual({x.name for x in (self.folder/'job').iterdir()}, {'controller-progress.json', 'controller-result.json'})

    def test_failed_directory_cannot_be_reused(self):
        with patch.object(launch, 'CACHE', self.folder), patch.object(launch, 'LAUNCH_ADMITTED', False):
            launch.run(self.folder/'absent', self.folder/'job')
            before = launch.raw(self.folder/'job/controller-result.json')
            with self.assertRaisesRegex(ValueError, 'new directory'):
                launch.run(self.folder/'absent', self.folder/'job')
            self.assertEqual(before, launch.raw(self.folder/'job/controller-result.json'))

    def test_output_must_be_local_recovery_descendant(self):
        with patch.object(launch, 'CACHE', self.folder/'allowed'):
            with self.assertRaises(ValueError):
                launch.run(self.folder/'absent', self.folder/'elsewhere')
        self.assertFalse((self.folder/'elsewhere').exists())

    def test_pin_rejects_changed_or_symlinked_bytes(self):
        p = self.folder/'source'; p.write_bytes(b'original')
        expected = launch.sha(p)
        p.write_bytes(b'changed')
        with self.assertRaises(ValueError): launch.raw(p, expected)
        link = self.folder/'link'; link.symlink_to(p)
        with self.assertRaises(ValueError): launch.raw(link)

    def test_inventory_rejects_missing_extra_and_symlink_artifacts(self):
        p = self.folder/'artifact'; p.write_bytes(b'synthetic')
        self.assertEqual(set(launch.inventory(self.folder, {'artifact'})), {'artifact'})
        for names in (set(), {'artifact', 'missing'}):
            with self.subTest(names=names), self.assertRaises(ValueError): launch.inventory(self.folder, names)
        p.unlink(); p.symlink_to(SOURCE)
        with self.assertRaises(ValueError): launch.inventory(self.folder, {'artifact'})

    def state(self):
        return dict(child_pid=123, peak_sampled_process_rss_bytes=20, peak_sampled_group_rss_bytes=30)

    def test_cleanup_failure_retains_live_pid(self):
        details = launch.failure_details(self.state(), RuntimeError('cleanup failed'), types.SimpleNamespace(_group_is_alive=lambda _: True))
        self.assertEqual(details['child_pid'], 123)
        self.assertEqual(details['cleanup_state'], 'live_group_requires_attention')

    def test_cleanup_uncertainty_retains_pid(self):
        def failed(_): raise RuntimeError('ps unavailable')
        details = launch.failure_details(self.state(), RuntimeError('failed'), types.SimpleNamespace(_group_is_alive=failed))
        self.assertEqual(details['child_pid'], 123)
        self.assertEqual(details['cleanup_state'], 'unknown_requires_attention')

    def test_cleanup_verified_absence_clears_pid(self):
        details = launch.failure_details(self.state(), RuntimeError('failed'), types.SimpleNamespace(_group_is_alive=lambda _: False))
        self.assertIsNone(details['child_pid'])
        self.assertEqual(details['cleanup_state'], 'verified_no_live_group')

    def test_actual_memory_crossing_is_retained(self):
        sample = {'maximum_process_rss_bytes': 8_000_000_001, 'process_group_rss_bytes': 8_000_000_010}
        details = launch.failure_details(self.state(), MemoryError('Owned quantum process exceeded sampled RSS limit: '+json.dumps(sample)), types.SimpleNamespace(_group_is_alive=lambda _: False))
        self.assertEqual(details['error_observed_resource_sample'], sample)
        self.assertEqual(details['peak_sampled_process_rss_bytes'], 8_000_000_001)
        self.assertEqual(details['peak_sampled_group_rss_bytes'], 8_000_000_010)

    def test_unparsed_memory_failure_remains_recorded(self):
        details = launch.failure_details(self.state(), MemoryError('unstructured failure'), types.SimpleNamespace(_group_is_alive=lambda _: False))
        self.assertEqual(details['error'], 'unstructured failure')
        self.assertNotIn('error_observed_resource_sample', details)

    def runtime_fixture(self):
        prefix = self.folder/'runtime'; prefix.mkdir()
        manifest = prefix/'dynamol-qm-runtime.json'; manifest.write_text('{"scope":"synthetic software fixture"}')
        source = prefix/'safe.py'; source.write_text('pass\n')
        evidence = prefix/'local-recovery-evidence'; evidence.mkdir()
        inv = {'runtime': str(prefix), 'files': [{'path': str(p), 'relative_path': p.name, 'sha256': launch.sha(p), 'bytes': p.stat().st_size} for p in (manifest, source)], 'symlinks': []}
        p = evidence/'runtime-files.json'; p.write_text(json.dumps(inv))
        return prefix, p, launch.sha(manifest)

    def test_runtime_rejects_unlisted_executable_payloads(self):
        prefix, inv, manifest_sha = self.runtime_fixture()
        with patch.multiple(launch, QM_PREFIX=prefix, INVENTORY=inv, INVENTORY_SHA=launch.sha(inv), RUNTIME_SHA=manifest_sha):
            launch.verify_runtime()
            for name in ('sitecustomize.py', 'injected.pth', 'new.dylib'):
                with self.subTest(name=name):
                    p = prefix/name; p.write_text('unreviewed')
                    with self.assertRaisesRegex(ValueError, 'payload inventory'): launch.verify_runtime()
                    p.unlink()

    def test_runtime_only_excludes_declared_bytecode_and_evidence(self):
        prefix, inv, manifest_sha = self.runtime_fixture()
        cache = prefix/'__pycache__'; cache.mkdir(); bytecode = cache/'safe.cpython-312.pyc'; bytecode.write_bytes(b'excluded, never executed')
        (prefix/'local-recovery-evidence/note.json').write_text('{}')
        with patch.multiple(launch, QM_PREFIX=prefix, INVENTORY=inv, INVENTORY_SHA=launch.sha(inv), RUNTIME_SHA=manifest_sha):
            launch.verify_runtime()
            bytecode.unlink(); bytecode.symlink_to(SOURCE)
            with self.assertRaises(ValueError): launch.verify_runtime()

    def test_runtime_rejects_changed_pinned_source(self):
        prefix, inv, manifest_sha = self.runtime_fixture()
        (prefix/'safe.py').write_text('raise RuntimeError()\n')
        with patch.multiple(launch, QM_PREFIX=prefix, INVENTORY=inv, INVENTORY_SHA=launch.sha(inv), RUNTIME_SHA=manifest_sha):
            with self.assertRaisesRegex(ValueError, 'Pinned bytes'): launch.verify_runtime()

    def test_native_nonzero_and_boolean_exit_reject_before_reads(self):
        for value in (1, True, False, 0.0, None):
            with self.subTest(value=repr(value)), self.assertRaisesRegex(ValueError, 'child failed'):
                launch.verify_native(self.folder/'absent', value, 'synthetic')

    def test_native_zero_exit_without_actual_outputs_rejects(self):
        with self.assertRaisesRegex(ValueError, 'inventory'):
            launch.verify_native(self.folder, 0, 'synthetic')

    def test_real_child_isolated_flags_and_cleanup(self):
        supervisor = launch.module(launch.SUPERVISOR, launch.SUPERVISOR_SHA, 'reviewed_test_supervisor')
        prefix = self.folder/'fresh-bytecode'; prefix.mkdir(); started = []; samples = []
        command = [sys.executable, '-I', '-B', '-X', 'pycache_prefix='+str(prefix), '-c',
                   'import json,sys; print(json.dumps({"isolated":sys.flags.isolated,"no_bytecode":sys.dont_write_bytecode,"cache":sys.pycache_prefix}))']
        code = supervisor.supervise(command, self.folder, 5, dict(os.environ), started.append, samples.append,
                                    128*1024**2, 0)
        self.assertEqual(code, 0)
        result = json.loads((self.folder/'controller-native.log').read_text())
        self.assertEqual(result, {'isolated': 1, 'no_bytecode': True, 'cache': str(prefix)})
        self.assertFalse(supervisor._group_is_alive(started[0]))
        self.assertFalse(list(prefix.rglob('*.pyc')))

    def test_real_child_deadline_leaves_no_live_group(self):
        supervisor = launch.module(launch.SUPERVISOR, launch.SUPERVISOR_SHA, 'reviewed_timeout_supervisor')
        started = []
        with self.assertRaises(TimeoutError):
            supervisor.supervise([sys.executable, '-I', '-B', '-c', 'import time;time.sleep(60)'],
                                 self.folder, .25, dict(os.environ), started.append, lambda _: None,
                                 128*1024**2, 0)
        self.assertEqual(len(started), 1)
        self.assertFalse(supervisor._group_is_alive(started[0]))


if __name__ == '__main__':
    unittest.main(verbosity=2)
