"""Bounded native loop subprocesses with owned-group cleanup.

Lifecycle functions are extracted unchanged from the reviewed local supervisor.
Resource limits are monitored bounds, never an acceptance result for a model.
"""
from __future__ import annotations
from contextlib import contextmanager
import json
import os
from pathlib import Path
import signal
import subprocess
import threading
import time


class Cancelled(Exception):
    """A supervisor cancellation, recorded only after its group is cleaned up."""

@contextmanager
def _cancellation_scope():
    """Latch cancellation without raising asynchronously at any boundary.

    Ordinary code checks the flag while polling. Signals therefore cannot abort
    Popen assignment, a finally block, or entry into the cleanup context itself.
    The child does not inherit a blocked signal mask.
    """
    if threading.current_thread() is not threading.main_thread():
        raise RuntimeError('The native process supervisor must run on the main thread')
    signals = {signal.SIGINT, signal.SIGTERM}
    state = {'requested': None}
    def cancel(signum, frame):
        if state['requested'] is None:
            state['requested'] = signum
    old_mask = signal.pthread_sigmask(signal.SIG_BLOCK, signals)
    handlers = {sig: signal.signal(sig, cancel) for sig in signals}
    try:
        signal.pthread_sigmask(signal.SIG_SETMASK, old_mask)
        yield state
    finally:
        signal.pthread_sigmask(signal.SIG_BLOCK, signals)
        for sig, handler in handlers.items():
            signal.signal(sig, handler)
        signal.pthread_sigmask(signal.SIG_SETMASK, old_mask)

@contextmanager
def _uninterrupted_cleanup():
    """Repeated cancellation cannot interrupt termination/reaping of our group."""
    signals = {signal.SIGINT, signal.SIGTERM}
    old_mask = signal.pthread_sigmask(signal.SIG_BLOCK, signals)
    handlers = {sig: signal.signal(sig, signal.SIG_IGN) for sig in signals}
    try:
        signal.pthread_sigmask(signal.SIG_SETMASK, old_mask)
        yield
    finally:
        signal.pthread_sigmask(signal.SIG_BLOCK, signals)
        for sig, handler in handlers.items():
            signal.signal(sig, handler)
        signal.pthread_sigmask(signal.SIG_SETMASK, old_mask)

def _group_is_alive(group_id):
    # A leader may already be reaped while its descendants remain. macOS can
    # return EPERM for killpg(..., 0) during teardown; inspect live group members.
    rows = subprocess.check_output(['ps', '-axo', 'pgid=,stat='], text=True,
                                   timeout=3).splitlines()
    return any(len(parts := row.split()) == 2 and int(parts[0]) == group_id and
               not parts[1].startswith('Z') for row in rows)

def _terminate_owned_group(child, grace_seconds=1.0):
    """Terminate descendants even if their leader exited; always reap the leader.

    macOS does not expose Linux subreapers: grandchildren are reaped by their
    parent or launchd. We verify there are no live executable group members.
    """
    try:
        try:
            os.killpg(child.pid, signal.SIGTERM)
        except ProcessLookupError:
            return
        except PermissionError:
            if not _group_is_alive(child.pid):
                return
            raise
        deadline = time.monotonic()+grace_seconds
        while time.monotonic() < deadline:
            child.poll()
            if not _group_is_alive(child.pid):
                return
            time.sleep(.025)
        try:
            os.killpg(child.pid, signal.SIGKILL)
        except ProcessLookupError:
            pass
        except PermissionError:
            if _group_is_alive(child.pid):
                raise
        deadline = time.monotonic()+3
        while _group_is_alive(child.pid):
            if time.monotonic() >= deadline:
                raise RuntimeError('Owned process group still has live members after SIGKILL')
            child.poll()
            time.sleep(.025)
    finally:
        child.wait(timeout=3)


def _sample_group(group_id):
    rows = subprocess.check_output(['ps', '-axo', 'pid=,pgid=,rss=,stat='],
                                   text=True, timeout=3).splitlines()
    members = []
    for row in rows:
        fields = row.split()
        if len(fields) != 4:
            raise ValueError('Malformed native resource sample.')
        pid, pgid, rss = map(int, fields[:3])
        if pgid == group_id and not fields[3].startswith('Z'):
            members.append({'pid': pid, 'rss_bytes': rss * 1024})
    return {'group_id': group_id, 'members': members,
            'rss_bytes': sum(row['rss_bytes'] for row in members)}


def supervise_loop_child(command, folder, *, environment, check_cancel=lambda: None,
                         timeout_seconds=180, max_rss_bytes=4 * 1024**3,
                         on_sample=lambda sample: None):
    """Return the actual exit code; always clean and reap our process group.

    A caller maps exit codes to explicit native result classes. Cancellation,
    deadline, memory and monitoring failures are never candidate rejection.
    Both generation and complete validation use this same bounded lifecycle.
    """
    if not 0 < timeout_seconds <= 1800 or not 0 < max_rss_bytes <= 8 * 1024**3:
        raise ValueError('Native loop work requires bounded time and memory.')
    if not command or not all(isinstance(arg, str) and arg for arg in command):
        raise ValueError('Native loop command must contain explicit arguments.')
    folder = Path(folder).resolve()
    folder.mkdir(parents=True, exist_ok=True)
    report_path = folder / 'supervisor.json'
    if report_path.exists():
        raise ValueError('Native attempt already has a supervisor report.')
    started = time.monotonic()
    child = None
    report = {'command': command, 'timeout_seconds': timeout_seconds,
              'max_rss_bytes': max_rss_bytes, 'samples': [], 'status': 'running'}
    try:
        check_cancel()
        with _cancellation_scope() as cancellation:
            with (folder / 'worker.log').open('x') as log:
                try:
                    child = subprocess.Popen(command, cwd=folder, stdout=log,
                                             stderr=subprocess.STDOUT,
                                             start_new_session=True, env=environment)
                    report['pid'] = child.pid
                    next_sample = started
                    while True:
                        check_cancel()
                        if cancellation['requested'] is not None:
                            raise Cancelled('Cancellation signal ' + str(cancellation['requested']))
                        now = time.monotonic()
                        if now - started >= timeout_seconds:
                            raise TimeoutError('Native loop attempt exceeded its deadline.')
                        if now >= next_sample:
                            sample = _sample_group(child.pid)
                            sample['elapsed_seconds'] = now - started
                            report['samples'].append(sample)
                            on_sample(sample)
                            if sample['rss_bytes'] > max_rss_bytes:
                                raise MemoryError('Native loop attempt exceeded its memory limit.')
                            next_sample = now + 1
                        remaining = timeout_seconds - (time.monotonic() - started)
                        if remaining <= 0:
                            raise TimeoutError('Native loop attempt exceeded its deadline.')
                        try:
                            code = child.wait(timeout=min(.1, remaining))
                            break
                        except subprocess.TimeoutExpired:
                            pass
                finally:
                    if child is not None:
                        with _uninterrupted_cleanup():
                            _terminate_owned_group(child)
                check_cancel()
                if cancellation['requested'] is not None:
                    raise Cancelled('Cancellation during native cleanup.')
                report.update(status='exited', returncode=code, owned_group_cleanup_completed=True)
                return code
    except BaseException as exc:
        report.update(status='cancelled' if isinstance(exc, Cancelled) else 'failed',
                      error_type=type(exc).__name__, error=str(exc))
        raise
    finally:
        report['elapsed_seconds'] = time.monotonic() - started
        report['peak_sampled_rss_bytes'] = max((r['rss_bytes'] for r in report['samples']), default=None)
        if child is not None:
            report['returncode'] = child.poll()
        temporary = report_path.with_suffix('.json.tmp')
        temporary.write_text(json.dumps(report, indent=2, allow_nan=False) + '\n')
        temporary.replace(report_path)
