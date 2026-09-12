#!/usr/bin/env python3
"""Offline first-run setup and a private local service; stdlib-only entry point."""
from __future__ import annotations
import argparse
import fcntl
import hashlib
import json
import os
from pathlib import Path
import shutil
import signal
import secrets
import socket
import subprocess
import sys
import tarfile
import threading
import time
import urllib.request
from http.server import BaseHTTPRequestHandler, ThreadingHTTPServer

RESOURCES = Path(__file__).resolve().parent
APP = RESOURCES / 'app'
PYTHON = RESOURCES / 'python' / 'bin' / 'python3.12'
DEFAULT_HOME = Path.home() / 'Library' / 'Application Support' / 'DynaMol'
CONTROL_TOKEN = secrets.token_urlsafe(32)
STATE = {'phase': 'Starting DynaMol', 'detail': 'Preparing your private workspace.', 'ready': False, 'error': None, 'url': None}

PAGE = '''<!doctype html><html><head><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1"><title>Opening DynaMol</title><style>body{margin:0;background:#0b121b;color:#e3f0ed;font:16px -apple-system,BlinkMacSystemFont,sans-serif;display:grid;place-items:center;min-height:100vh}.card{width:min(460px,80vw);padding:42px;border:1px solid #294139;border-radius:22px;background:linear-gradient(130deg,#142a29,#121e2a)}h1{font-size:30px;letter-spacing:-1px}p{color:#b7c8cc;line-height:1.65}.wheel{border:3px solid #28483f;border-top-color:#83e5c5;width:28px;height:28px;border-radius:50%;animation:spin 1s linear infinite}@keyframes spin{to{transform:rotate(360deg)}}small{color:#91aaa4}pre{white-space:pre-wrap;color:#f4bc97;font-size:13px}</style></head><body><main class="card"><div class="wheel" id="wheel"></div><h1>Dyna<span style="color:#83e5c5">Mol</span></h1><h3 id="phase">Starting DynaMol</h3><p id="detail">Preparing your private workspace.</p><small>Your molecules stay on this computer. First launch unpacks the included engines; no downloads or installers are needed.</small><pre id="error"></pre></main><script>async function poll(){try{const s=await(await fetch('/status',{cache:'no-store'})).json();document.querySelector('#phase').textContent=s.phase;document.querySelector('#detail').textContent=s.detail;if(s.error){document.querySelector('#error').textContent=s.error;document.querySelector('#wheel').style.display='none';return}if(s.ready&&s.url){location.replace(s.url);return}}catch(e){}setTimeout(poll,600)}poll()</script></body></html>'''


def atomic_json(path: Path, value: object):
    temporary = path.with_suffix('.tmp')
    temporary.write_text(json.dumps(value, indent=2) + '\n')
    temporary.replace(path)


def digest(path: Path) -> str:
    with path.open('rb') as stream:
        return hashlib.file_digest(stream, 'sha256').hexdigest()


def private_environment(data: Path, engines: dict[str, Path]) -> dict[str, str]:
    env = {key: value for key, value in os.environ.items() if key in {'HOME', 'USER', 'LOGNAME', 'TMPDIR', 'LANG', 'LC_ALL', 'SHELL'}}
    env.update(PATH='/usr/bin:/bin:/usr/sbin:/sbin', PYTHONNOUSERSITE='1', PYTHONDONTWRITEBYTECODE='1', PYTHONUNBUFFERED='1', DYNAMOL_DATA_DIR=str(data), DYNAMOL_GMX=str(engines['gromacs'] / 'bin' / 'gmx'), DYNAMOL_AMBERTOOLS=str(engines['ambertools']), OPENMM_CPU_THREADS='2', OMP_NUM_THREADS='2', OPENBLAS_NUM_THREADS='2', MKL_NUM_THREADS='2')
    return env


def ensure_engines(home: Path, manifest: dict, log) -> dict[str, Path]:
    roots: dict[str, Path] = {}
    for name, record in manifest['engines'].items():
        archive = RESOURCES / 'engines' / record['archive']
        target = home / 'Runtimes' / f"{name}-{record['sha256'][:16]}"
        marker = target / '.dynamol-complete.json'
        if marker.exists():
            saved = json.loads(marker.read_text())
            if saved.get('sha256') == record['sha256'] and saved.get('prefix') == str(target):
                roots[name] = target
                continue
        STATE.update(phase=f"Preparing {record['label']}", detail='Unpacking the engine included with this application. This happens once.')
        if digest(archive) != record['sha256']:
            raise RuntimeError(f"The bundled {name} archive failed its integrity check. Download a fresh copy of DynaMol.")
        # conda-unpack fixes absolute paths once. Extract to its FINAL prefix,
        # protected by the launch lock, and mark success only after relocation.
        if target.exists():
            shutil.rmtree(target)
        target.mkdir(parents=True)
        with tarfile.open(archive, 'r:gz') as packed:
            packed.extractall(target, filter='data')
        unpack = target / 'bin' / 'conda-unpack'
        # -I ignores PYTHONDONTWRITEBYTECODE. Explicit -B keeps first-run
        # engine relocation from adding caches inside the sealed app bundle.
        subprocess.run([str(PYTHON), '-I', '-B', str(unpack)], check=True, stdout=log, stderr=log, env={**os.environ, 'PATH': '/usr/bin:/bin:/usr/sbin:/sbin'}, timeout=300)
        atomic_json(marker, {'sha256': record['sha256'], 'prefix': str(target)})
        roots[name] = target
    return roots


def live_state(path: Path) -> dict | None:
    try:
        state = json.loads(path.read_text())
        os.kill(int(state['pid']), 0)
        with urllib.request.urlopen(state['url'] + 'api/health', timeout=1) as response:
            if response.status == 200:
                return state
    except (ValueError, OSError, KeyError, urllib.error.URLError):
        pass
    return None


class ProgressHandler(BaseHTTPRequestHandler):
    def do_GET(self):
        body = json.dumps(STATE).encode() if self.path == '/status' else PAGE.encode()
        self.send_response(200)
        self.send_header('Content-Type', 'application/json' if self.path == '/status' else 'text/html; charset=utf-8')
        self.send_header('Cache-Control', 'no-store')
        self.send_header('Content-Length', str(len(body)))
        self.end_headers()
        self.wfile.write(body)
    def do_POST(self):
        if self.path != '/shutdown' or not secrets.compare_digest(self.headers.get('X-DynaMol-Key', ''), CONTROL_TOKEN):
            self.send_error(403)
            return
        self.send_response(204)
        self.end_headers()
        threading.Timer(0.1, lambda: os.kill(os.getpid(), signal.SIGTERM)).start()
    def log_message(self, *args):
        pass


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument('--home', type=Path, default=DEFAULT_HOME, help=argparse.SUPPRESS)
    parser.add_argument('--no-browser', action='store_true', help=argparse.SUPPRESS)
    parser.add_argument('--prepare-only', action='store_true', help=argparse.SUPPRESS)
    args = parser.parse_args()
    def terminate(signum, frame):
        raise KeyboardInterrupt
    signal.signal(signal.SIGTERM, terminate)
    home = args.home.expanduser().resolve()
    home.mkdir(parents=True, exist_ok=True)
    state_path = home / 'service.json'
    existing = live_state(state_path)
    if existing:
        if not args.no_browser:
            subprocess.run(['/usr/bin/open', existing['url']], check=False)
        print(existing['url'], flush=True)
        return 0
    lock = (home / 'launcher.lock').open('a+')
    try:
        fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
    except BlockingIOError:
        for _ in range(30):
            existing = live_state(state_path)
            if existing:
                if not args.no_browser:
                    subprocess.run(['/usr/bin/open', existing['url']], check=False)
                print(existing['url'], flush=True)
                return 0
            try:
                fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
                break
            except BlockingIOError:
                time.sleep(0.1)
        else:
            startup = home / 'startup.json'
            if startup.exists() and not args.no_browser:
                address = json.loads(startup.read_text()).get('url', '')
                if address.startswith('http://127.0.0.1:'):
                    subprocess.run(['/usr/bin/open', address], check=False)
            print('DynaMol is already starting. Please wait for its browser window.', flush=True)
            return 0
    progress = ThreadingHTTPServer(('127.0.0.1', 0), ProgressHandler)
    threading.Thread(target=progress.serve_forever, daemon=True).start()
    progress_url = f'http://127.0.0.1:{progress.server_port}/'
    atomic_json(home / 'startup.json', {'url': progress_url, 'pid': os.getpid()})
    print(json.dumps({'phase': 'starting', 'progress_url': progress_url}), flush=True)
    if not args.no_browser:
        subprocess.run(['/usr/bin/open', progress_url], check=False)
    logs = home / 'Logs'
    logs.mkdir(exist_ok=True)
    child = None
    try:
        manifest = json.loads((RESOURCES / 'manifest.json').read_text())
        with (logs / 'launcher.log').open('a', buffering=1) as log:
            engines = ensure_engines(home, manifest, log)
            env = private_environment(home / 'Workspace', engines)
            for name in ('demo', 'ubiquitin-start'):
                source, target = APP / 'examples' / name, home / 'Workspace' / 'datasets' / name
                if not target.exists():
                    target.parent.mkdir(parents=True, exist_ok=True)
                    shutil.copytree(source, target)
            atomic_json(home / 'runtime.json', {'python': str(PYTHON), 'engines': {key: str(value) for key, value in engines.items()}, 'environment': env | {'HOME': str(Path.home())}, 'manifest': manifest['build_id']})
            if args.prepare_only:
                print(json.dumps({'prepared': True, 'home': str(home), 'engines': {key: str(value) for key, value in engines.items()}}), flush=True)
                return 0
            STATE.update(phase='Opening your molecular workspace', detail='Starting the local service with the bundled engines.')
            listener = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
            listener.bind(('127.0.0.1', 0))
            listener.listen(128)
            url = f'http://127.0.0.1:{listener.getsockname()[1]}/'
            env['DYNAMOL_ORIGIN'] = url.rstrip('/')
            child = subprocess.Popen([str(PYTHON), '-B', '-m', 'uvicorn', 'backend.main:app', '--fd', str(listener.fileno()), '--log-level', 'info'], cwd=APP, env=env, pass_fds=(listener.fileno(),), stdout=log, stderr=log)
            listener.close()
            for _ in range(120):
                if child.poll() is not None:
                    raise RuntimeError(f'The local service exited. Details: {logs / "launcher.log"}')
                try:
                    with urllib.request.urlopen(url + 'api/health', timeout=1) as response:
                        if response.status == 200:
                            break
                except OSError:
                    time.sleep(0.25)
            else:
                raise RuntimeError(f'The local service did not start. Details: {logs / "launcher.log"}')
            atomic_json(state_path, {'url': url, 'pid': child.pid, 'launcher_pid': os.getpid(), 'control_url': progress_url + 'shutdown', 'control_token': CONTROL_TOKEN, 'build_id': manifest['build_id']})
            STATE.update(ready=True, url=url)
            print(url, flush=True)
            return child.wait()
    except KeyboardInterrupt:
        return 0
    except Exception as error:
        STATE.update(phase='DynaMol could not open', detail='Your existing molecules and jobs have been kept.', error=str(error))
        print(str(error), file=sys.stderr, flush=True)
        if not args.no_browser:
            time.sleep(120)
        return 1
    finally:
        if child and child.poll() is None:
            child.terminate()
            child.wait(timeout=15)
        if child:
            state_path.unlink(missing_ok=True)
        progress.shutdown()
        (home / 'startup.json').unlink(missing_ok=True)
        lock.close()


if __name__ == '__main__':
    raise SystemExit(main())
