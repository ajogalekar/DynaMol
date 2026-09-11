#!/usr/bin/env python3
"""Verify packaged startup, private worker launches, reuse and restart."""
from __future__ import annotations
import argparse
import json
import os
from pathlib import Path
import subprocess
import time
import urllib.request


def fetch(url, payload=None, headers=None):
    data = None if payload is None else json.dumps(payload).encode()
    request = urllib.request.Request(url, data=data, headers={'Content-Type': 'application/json', **(headers or {})})
    with urllib.request.urlopen(request, timeout=20) as response:
        raw = response.read()
        return json.loads(raw) if raw else None


def stop(home, process=None):
    state_path = home / 'service.json'
    if not state_path.exists():
        return
    state = json.loads(state_path.read_text())
    fetch(state['control_url'], {}, {'X-DynaMol-Key': state['control_token']})
    for _ in range(100):
        if not state_path.exists():
            if process is not None and process.poll() is not None:
                return
            try:
                os.kill(state['launcher_pid'], 0)
            except ProcessLookupError:
                return
        time.sleep(0.1)
    raise RuntimeError('Packaged service did not shut down cleanly.')


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('resources', type=Path)
    parser.add_argument('home', type=Path)
    args = parser.parse_args()
    resources, home = args.resources.resolve(), args.home.resolve()
    home.mkdir(parents=True, exist_ok=True)
    stop(home)
    executable = resources/'python/bin/python3.12'
    command = [str(executable), '-I', '-B', str(resources/'bootstrap.py'), '--no-browser', '--home', str(home)]
    env = {'PATH': '/usr/bin:/bin', 'HOME': str(Path.home())}
    children = []
    def start():
        child = subprocess.Popen(command, env=env, stdout=subprocess.PIPE, stderr=(home/'app-check.stderr.log').open('a'), text=True)
        children.append(child)
        first = json.loads(child.stdout.readline())
        progress_url = first['progress_url']
        with urllib.request.urlopen(progress_url, timeout=5) as response:
            html = response.read().decode()
            assert 'wheel' in html and 'Starting DynaMol' in html
        url = child.stdout.readline().strip()
        assert url.startswith('http://127.0.0.1:')
        state = json.loads((home/'service.json').read_text())
        assert state['url'] == url
        return child, state
    try:
        child, state = start()
        url = state['url']
        health = fetch(url+'api/health')
        assert {item['id'] for item in health['engines'] if item['available']} == {'openmm','gromacs'}
        library = fetch(url+'api/datasets')
        existing_jobs = fetch(url+'api/jobs')
        generated = {item.get('dataset_id') for item in existing_jobs if item.get('name', '').startswith('Packaged ')}
        assert {item['id'] for item in library} - generated == {'demo','ubiquitin-start'}
        reuse = subprocess.run(command, env=env, capture_output=True, text=True, timeout=15)
        assert reuse.returncode == 0 and reuse.stdout.strip() == url
        assert json.loads((home/'service.json').read_text())['pid'] == state['pid']
        assert fetch(url+'api/datasets/demo')['n_frames'] > 1
        jobs = []
        for engine in ('openmm','gromacs'):
            previous = next((item for item in existing_jobs if item.get('name') == f'Packaged {engine} worker check' and item['status'] == 'completed'), None)
            job = previous or fetch(url+'api/jobs', {'dataset_id': 'ubiquitin-start', 'engine': engine, 'name': f'Packaged {engine} worker check', 'duration_ps': 0.02, 'report_interval': 2, 'equilibration_steps': 0, 'minimize': True, 'solvent': 'implicit' if engine=='openmm' else 'explicit'})
            for _ in range(900):
                status = fetch(url+f"api/jobs/{job['id']}")
                if status['status'] not in ('running','queued'):
                    break
                time.sleep(0.25)
            assert status['status'] == 'completed', status
            dataset = fetch(url+f"api/datasets/{status['dataset_id']}")
            assert dataset['n_frames'] >= 2
            folder = home/'Workspace/jobs'/job['id']
            provenance = json.loads((folder/'provenance.json').read_text())
            # Both workers run through the private interpreter; their saved
            # source hash must match the exact code copied into this app.
            import hashlib
            expected = hashlib.sha256((resources/'app/backend/worker.py').read_bytes()).hexdigest()
            assert provenance['worker_source_sha256'] == expected
            jobs.append({'engine': engine, 'job_id': job['id'], 'status': status['status'], 'dataset_id': dataset['id'], 'atoms': dataset['n_atoms'], 'frames': dataset['n_frames'], 'worker_source_sha256': expected})
        before_restart = {item['id'] for item in fetch(url+'api/datasets')}
        stop(home, child)
        assert child.wait(timeout=15) == 0
        restarted, new_state = start()
        assert before_restart == {item['id'] for item in fetch(new_state['url']+'api/datasets')}
        assert len(fetch(new_state['url']+'api/jobs')) == 2
        stop(home, restarted)
        assert restarted.wait(timeout=15) == 0
        report = {'passed': True, 'resources': str(resources), 'home': str(home), 'build_id': state['build_id'], 'health': health, 'initial_datasets': sorted(item['id'] for item in library), 'startup_progress_page': True, 'same_service_reused': True, 'shutdown_and_restart_preserved_data': True, 'sanitized_path': env['PATH'], 'jobs': jobs, 'limitations': 'Local relocated app with paths containing spaces on the development Mac; no separate clean-machine, notarization, or Gatekeeper validation.'}
        (home/'app-validation.json').write_text(json.dumps(report, indent=2)+'\n')
        print(json.dumps({'passed': True, 'report': str(home/'app-validation.json')},indent=2))
    finally:
        for child in children:
            if child.poll() is None:
                child.terminate()
                child.wait(timeout=15)


if __name__ == '__main__':
    main()
