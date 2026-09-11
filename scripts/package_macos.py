#!/usr/bin/env python3
"""Build the local, offline Apple Silicon application. No publishing/signing service."""
from __future__ import annotations
import argparse
import datetime as dt
import hashlib
import importlib.metadata
import json
import os
import re
from pathlib import Path
import platform
import plistlib
import shutil
import subprocess
import sys

ROOT = Path(__file__).resolve().parents[1]
STAGE = ROOT / 'build' / 'packaging'
RELEASES = ROOT / 'build' / 'releases'
VERSION = '0.1.0'
ENGINES = {'ambertools': (ROOT / '.tools' / 'ambertools', '24.8', 'AmberTools'), 'gromacs': (ROOT / '.gromacs', '2025.4', 'GROMACS')}


def sha(path: Path) -> str:
    with path.open('rb') as stream:
        return hashlib.file_digest(stream, 'sha256').hexdigest()


def command(*args: str | Path, **kwargs):
    return subprocess.run([str(value) for value in args], check=True, **kwargs)


def ignore_copy(directory, names):
    return [name for name in names if name in {'__pycache__', '.DS_Store', '_virtualenv.pth', '_virtualenv.py'} or name.endswith(('.pyc', '.pyo'))]


def prepare_runtime() -> dict:
    if platform.system() != 'Darwin' or platform.machine() != 'arm64':
        raise SystemExit('Build this target on Apple Silicon macOS.')
    STAGE.mkdir(parents=True, exist_ok=True)
    python_source = Path(command(ROOT / '.venv' / 'bin' / 'python', '-c', 'import sys; print(sys.base_prefix)', capture_output=True, text=True).stdout.strip())
    destination = STAGE / 'python'
    if not (destination / '.dynamol-runtime.json').exists():
        if destination.exists():
            shutil.rmtree(destination)
        shutil.copytree(python_source, destination, symlinks=True, ignore=ignore_copy)
        packages = destination / 'lib' / 'python3.12' / 'site-packages'
        shutil.copytree(ROOT / '.venv' / 'lib' / 'python3.12' / 'site-packages', packages, dirs_exist_ok=True, symlinks=True, ignore=ignore_copy)
        (destination / '.dynamol-runtime.json').write_text(json.dumps({'source_python_build': (python_source / 'BUILD').read_text().strip(), 'uv_lock_sha256': sha(ROOT / 'uv.lock')}, indent=2) + '\n')
    else:
        prior = json.loads((destination / '.dynamol-runtime.json').read_text())
        if prior['uv_lock_sha256'] != sha(ROOT / 'uv.lock'):
            raise SystemExit('uv.lock changed: remove build/packaging/python and prepare the runtime again.')
    # A locally built wheel can advertise the host OS even when its native
    # extension has an older deployment floor. The Mach-O load command is the
    # relevant floor for this bundled interpreter (no recipient pip install).
    extension = destination / 'lib/python3.12/site-packages/parmed/amber/_rdparm.cpython-312-darwin.so'
    commands = command('/usr/bin/otool', '-l', extension, capture_output=True, text=True).stdout
    versions = re.findall(r'\bminos (\d+)\.(\d+)', commands)
    if not versions or max(int(major) for major, _ in versions) > 14:
        raise SystemExit('Packaged ParmEd extension does not establish macOS 14 compatibility. Rebuild it with the required deployment target.')
    marker = destination / '.dynamol-runtime.json'
    details = json.loads(marker.read_text())
    details['native_deployment_floor'] = {'parmed': '.'.join(max(versions, key=lambda pair: tuple(map(int, pair))))}
    if 'packaged_overrides' in details:
        details['packaged_overrides']['parmed']['reason'] = 'Rebuilt pinned sdist in the private build environment; Mach-O target is independently checked. Source app environment unchanged.'
    marker.write_text(json.dumps(details, indent=2) + '\n')
    engines = {}
    packer = ROOT / '.tools' / 'packaging' / 'bin' / 'conda-pack'
    for key, (prefix, version, label) in ENGINES.items():
        archive = STAGE / f'{key}-{version}-macos-arm64.tar.gz'
        if not archive.exists():
            if not packer.exists():
                raise SystemExit('Build tool missing. Create .tools/packaging with conda-pack==0.8.1 and setuptools==80.9.0.')
            command(packer, '--prefix', prefix, '--output', archive, '--n-threads', '2', '--compress-level', '4')
        engines[key] = {'archive': archive.name, 'version': version, 'label': label, 'sha256': sha(archive)}
    return {'engines': engines, 'python': json.loads((destination / '.dynamol-runtime.json').read_text())}


def notices(destination: Path):
    destination.mkdir(parents=True, exist_ok=True)
    inventory = {'native': [], 'python': [], 'coverage': 'Bundled package license texts and recipes/source URLs are included. Corresponding source archives for all copyleft binaries have not been assembled or legally reviewed; this is a local prototype, not a public release compliance certification.'}
    for family, (prefix, _, _) in ENGINES.items():
        for metadata in sorted((prefix / 'conda-meta').glob('*.json')):
            package = json.loads(metadata.read_text())
            key = f"{package['name']}-{package['version']}-{package['build']}"
            target = destination / family / key
            target.mkdir(parents=True, exist_ok=True)
            cache = Path(package.get('extracted_package_dir') or package.get('link', {}).get('source', ''))
            found = []
            for name in ('licenses', 'recipe'):
                source = cache / 'info' / name
                if cache.is_absolute() and source.is_dir():
                    shutil.copytree(source, target / name, symlinks=False, ignore=ignore_copy)
                    found.append(name)
            installed_licenses = []
            for relative in package.get('files', []):
                if any(word in relative.lower() for word in ('license', 'copying', 'copyright')):
                    installed = prefix / relative
                    if installed.is_file():
                        copied = target / 'installed-license-files' / relative
                        copied.parent.mkdir(parents=True, exist_ok=True)
                        shutil.copy2(installed, copied)
                        installed_licenses.append(relative)
            if installed_licenses:
                found.append('installed-license-files')
            record = {'name': package['name'], 'version': package['version'], 'build': package['build'], 'license': package.get('license'), 'url': package.get('url'), 'sha256': package.get('sha256'), 'materials': found, 'installed_license_files': installed_licenses}
            (target / 'package.json').write_text(json.dumps(record, indent=2) + '\n')
            inventory['native'].append(record | {'engine': family})
    packages = STAGE / 'python' / 'lib' / 'python3.12' / 'site-packages'
    for distribution in sorted(importlib.metadata.distributions(path=[str(packages)]), key=lambda dist: dist.metadata.get('Name', '')):
        meta = distribution.metadata
        inventory['python'].append({'name': meta.get('Name'), 'version': distribution.version, 'license_expression': meta.get('License-Expression'), 'license': meta.get('License'), 'home_page': meta.get('Home-page'), 'project_urls': meta.get_all('Project-URL'), 'license_files': meta.get_all('License-File')})
    (destination / 'DEPENDENCY_INVENTORY.json').write_text(json.dumps(inventory, indent=2) + '\n')
    shutil.copy2(ROOT / 'THIRD_PARTY.md', destination / 'THIRD_PARTY.md')
    shutil.copy2(ROOT / 'LICENSE', destination / 'DynaMol-LICENSE')
    shutil.copy2(ROOT / 'packaging' / 'README.md', destination / 'PACKAGING-README.md')


def build() -> Path:
    runtime = prepare_runtime()
    if not (ROOT / 'frontend' / 'dist' / 'index.html').exists():
        raise SystemExit('Build the frontend before packaging.')
    RELEASES.mkdir(parents=True, exist_ok=True)
    bundle = RELEASES / 'DynaMol.app'
    if bundle.exists():
        shutil.rmtree(bundle)
    resources = bundle / 'Contents' / 'Resources'
    macos = bundle / 'Contents' / 'MacOS'
    resources.mkdir(parents=True)
    macos.mkdir()
    shutil.copytree(STAGE / 'python', resources / 'python', symlinks=True, ignore=ignore_copy)
    (resources / 'engines').mkdir()
    for record in runtime['engines'].values():
        shutil.copy2(STAGE / record['archive'], resources / 'engines' / record['archive'])
    app = resources / 'app'
    app.mkdir()
    for folder in ('backend', 'examples'):
        shutil.copytree(ROOT / folder, app / folder, symlinks=False, ignore=ignore_copy)
    shutil.copytree(ROOT / 'frontend' / 'dist', app / 'frontend' / 'dist')
    # Include the cited audit evidence and screenshots alongside the guides.
    shutil.copytree(ROOT / 'docs', app / 'docs', ignore=ignore_copy)
    for name in ('pyproject.toml', 'uv.lock', 'LICENSE', 'README.md', 'THIRD_PARTY.md'):
        shutil.copy2(ROOT / name, app / name)
    # The app source is included alongside the binaries; user datasets/jobs are
    # deliberately never inputs to this builder.
    shutil.copytree(ROOT / 'frontend' / 'src', app / 'frontend-source' / 'src', ignore=ignore_copy)
    shutil.copytree(ROOT / 'frontend' / 'tests', app / 'frontend-source' / 'tests', ignore=ignore_copy)
    for name in ('package.json', 'package-lock.json', 'vite.config.ts', 'tsconfig.json', 'tsconfig.app.json', 'tsconfig.node.json', 'index.html', 'playwright.config.ts', 'TESTING.md'):
        if (ROOT / 'frontend' / name).exists():
            shutil.copy2(ROOT / 'frontend' / name, app / 'frontend-source' / name)
    shutil.copytree(ROOT / 'scripts', app / 'scripts', ignore=ignore_copy)
    shutil.copytree(ROOT / 'packaging', app / 'packaging', ignore=ignore_copy)
    shutil.copy2(ROOT / 'packaging' / 'macos' / 'bootstrap.py', resources / 'bootstrap.py')
    notices(resources / 'Notices')
    source_hashes = {str(path.relative_to(app)): sha(path) for path in sorted(app.rglob('*')) if path.is_file()}
    build_id = hashlib.sha256(json.dumps(source_hashes, sort_keys=True).encode()).hexdigest()[:16]
    manifest = {'application': 'DynaMol', 'version': VERSION, 'target': 'macos-arm64', 'built_at': dt.datetime.now(dt.timezone.utc).isoformat(), 'build_id': build_id, **runtime, 'source_hashes': source_hashes, 'signing': 'ad-hoc launcher only; not Apple Developer signed or notarized'}
    (resources / 'manifest.json').write_text(json.dumps(manifest, indent=2) + '\n')
    info = {'CFBundleIdentifier': 'org.dynamol.desktop', 'CFBundleName': 'DynaMol', 'CFBundleDisplayName': 'DynaMol', 'CFBundleExecutable': 'DynaMol', 'CFBundlePackageType': 'APPL', 'CFBundleShortVersionString': VERSION, 'CFBundleVersion': '1', 'LSMinimumSystemVersion': '14.0', 'LSUIElement': True, 'NSHighResolutionCapable': True, 'CFBundleIconFile': 'DynaMol.icns'}
    (bundle / 'Contents' / 'Info.plist').write_bytes(plistlib.dumps(info))
    iconset = STAGE / 'DynaMol.iconset'
    command('/usr/bin/xcrun', 'swift', ROOT / 'packaging' / 'macos' / 'DrawIcon.swift', iconset)
    command('/usr/bin/iconutil', '-c', 'icns', iconset, '-o', resources / 'DynaMol.icns')
    command('/usr/bin/xcrun', 'swiftc', ROOT / 'packaging' / 'macos' / 'Launcher.swift', '-o', STAGE / 'DynaMolLauncher', '-framework', 'AppKit', '-target', 'arm64-apple-macosx14.0')
    command('/usr/bin/codesign', '--force', '--sign', '-', STAGE / 'DynaMolLauncher')
    shutil.copyfile(STAGE / 'DynaMolLauncher', macos / 'DynaMol')
    (macos / 'DynaMol').chmod(0o755)
    archive = RELEASES / f'DynaMol-{VERSION}-macos-arm64.zip'
    archive.unlink(missing_ok=True)
    command('/usr/bin/ditto', '-c', '-k', '--sequesterRsrc', '--keepParent', bundle, archive)
    checksum = sha(archive)
    archive.with_suffix('.zip.sha256').write_text(f'{checksum}  {archive.name}\n')
    print(json.dumps({'bundle': str(bundle), 'archive': str(archive), 'sha256': checksum, 'build_id': build_id}, indent=2))
    return bundle


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--prepare-runtimes', action='store_true')
    args = parser.parse_args()
    if args.prepare_runtimes:
        print(json.dumps(prepare_runtime(), indent=2))
    else:
        build()
