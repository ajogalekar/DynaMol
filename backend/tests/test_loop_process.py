import json
import os
from pathlib import Path
import sys
import time

import pytest

from backend import loop_process as lp


def launch(tmp_path, code, **kwargs):
    return lp.supervise_loop_child([sys.executable, '-c', code], tmp_path,
                                   environment=dict(os.environ), **kwargs)


def test_nonzero_code_is_returned_without_becoming_quality_acceptance(tmp_path):
    assert launch(tmp_path, 'raise SystemExit(2)') == 2
    report = json.loads((tmp_path / 'supervisor.json').read_text())
    assert report['returncode'] == 2
    assert report['owned_group_cleanup_completed']


def test_timeout_cleans_owned_group(tmp_path):
    with pytest.raises(TimeoutError):
        launch(tmp_path, 'import time; time.sleep(60)', timeout_seconds=.1)
    report = json.loads((tmp_path / 'supervisor.json').read_text())
    assert not lp._group_is_alive(report['pid'])
    assert report['returncode'] is not None


def test_cancel_during_native_work_cleans_group(tmp_path):
    started = time.monotonic()
    def cancel():
        if time.monotonic() - started > .1:
            raise lp.Cancelled('Test cancellation')
    with pytest.raises(lp.Cancelled):
        launch(tmp_path, 'import time; time.sleep(60)', check_cancel=cancel)
    report = json.loads((tmp_path / 'supervisor.json').read_text())
    assert report['status'] == 'cancelled'
    assert not lp._group_is_alive(report['pid'])


def test_leader_exit_does_not_leave_background_descendant(tmp_path):
    code = "import subprocess,sys; subprocess.Popen([sys.executable,'-c','import time; time.sleep(60)'])"
    assert launch(tmp_path, code) == 0
    report = json.loads((tmp_path / 'supervisor.json').read_text())
    assert not lp._group_is_alive(report['pid'])


def test_memory_failure_stops_without_reusing_a_stale_sample(tmp_path, monkeypatch):
    monkeypatch.setattr(lp, '_sample_group', lambda pid: {'group_id': pid, 'members': [], 'rss_bytes': 2})
    with pytest.raises(MemoryError):
        launch(tmp_path, 'import time; time.sleep(60)', max_rss_bytes=1)
    report = json.loads((tmp_path / 'supervisor.json').read_text())
    assert report['error_type'] == 'MemoryError'
    assert not lp._group_is_alive(report['pid'])


def test_slow_monitor_cannot_turn_a_late_exit_into_success(tmp_path, monkeypatch):
    original = lp._sample_group
    def delayed(pid):
        time.sleep(.15)
        return original(pid)
    monkeypatch.setattr(lp, '_sample_group', delayed)
    with pytest.raises(TimeoutError):
        launch(tmp_path, 'pass', timeout_seconds=.1)


def test_lifecycle_matches_the_reviewed_source():
    import ast
    root = Path(__file__).resolve().parents[2]
    source = ast.parse((root / 'docs/audit/advanced-chemistry/local_qm_supervisor.py').read_text())
    candidate = ast.parse(Path(lp.__file__).read_text())
    names = {'Cancelled', '_cancellation_scope', '_uninterrupted_cleanup', '_group_is_alive', '_terminate_owned_group'}
    functions = lambda tree: {node.name: ast.dump(node, include_attributes=False)
                              for node in tree.body if isinstance(node, (ast.FunctionDef, ast.ClassDef)) and node.name in names}
    assert functions(source) == functions(candidate)
