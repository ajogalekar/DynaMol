#!/usr/bin/env python3
"""Create and verify a drag-to-Applications DMG with built-in macOS tools and a build-only layout writer.

The input app is never modified. The ds_store/mac_alias libraries write installer metadata without controlling
Finder. No application is installed in /Applications or run.
"""
from __future__ import annotations

import argparse
from contextlib import contextmanager
from datetime import datetime, timezone
import errno
import hashlib
import importlib.metadata
import json
import os
from pathlib import Path
import platform
import plistlib
import re
import shutil
import stat
import subprocess
import tempfile
import time

ROOT = Path(__file__).resolve().parents[1]
HDIUTIL = '/usr/bin/hdiutil'


def run(*args: str | Path, **kwargs) -> subprocess.CompletedProcess:
    return subprocess.run([str(arg) for arg in args], check=True, **kwargs)


def sha(path: Path) -> str:
    with path.open('rb') as stream:
        return hashlib.file_digest(stream, 'sha256').hexdigest()


def verify_app_signature(app: Path) -> None:
    """Reject an incomplete or modified app before distributing its bytes.

    Ad-hoc integrity is separate from Developer ID trust and notarization.
    Signing a standalone launcher does not seal its containing app bundle.
    """
    result = subprocess.run(
        ['/usr/bin/codesign', '--verify', '--deep', '--strict', '--verbose=2', str(app)],
        capture_output=True, text=True,
    )
    if result.returncode:
        raise ValueError(f'App signature integrity failed at {app}: {result.stderr.strip()}')


def inventory(app: Path) -> dict:
    """Include every regular file, directory mode, and literal symlink target."""
    result = {}
    for directory, dirs, files in os.walk(app, followlinks=False):
        for name in sorted(dirs + files):
            path = Path(directory) / name
            relative = str(path.relative_to(app))
            info = path.lstat()
            record = {'mode': stat.S_IMODE(info.st_mode)}
            if path.is_symlink():
                record.update(type='symlink', target=os.readlink(path))
            elif path.is_file():
                record.update(type='file', bytes=info.st_size, sha256=sha(path))
            elif path.is_dir():
                record.update(type='directory')
            else:
                raise ValueError(f'Unsupported special file inside app: {relative}')
            result[relative] = record
    return dict(sorted(result.items()))


def inventory_hash(records: dict) -> str:
    return hashlib.sha256(json.dumps(records, sort_keys=True, separators=(',', ':')).encode()).hexdigest()


def assert_inventory(app: Path, expected: dict) -> None:
    actual = inventory(app)
    if actual != expected:
        changed = sorted(key for key in expected.keys() | actual.keys() if expected.get(key) != actual.get(key))
        raise ValueError(f'App content/mode mismatch at {app}: {changed[:10]}')


@contextmanager
def mounted(image: Path, location: Path, readonly: bool):
    args = [HDIUTIL, 'attach', image, '-nobrowse', '-noautoopen', '-mountpoint', location, '-plist']
    if readonly:
        args.append('-readonly')
    result = plistlib.loads(run(*args, capture_output=True).stdout)
    devices = [entry['dev-entry'] for entry in result['system-entities'] if entry.get('dev-entry')]
    mounted_paths = [entry['mount-point'] for entry in result['system-entities'] if entry.get('mount-point')]
    if str(location) not in mounted_paths:
        if devices:
            run(HDIUTIL, 'detach', devices[0])
        raise RuntimeError('Disk image did not mount at its isolated requested path')
    try:
        yield result
    finally:
        # Retry ordinary detach; never detach an unrelated device or silently
        # force an installer volume containing a user's open document.
        for attempt in range(4):
            detached = subprocess.run([HDIUTIL, 'detach', devices[0]], capture_output=True, text=True)
            if detached.returncode == 0:
                break
            time.sleep(1)
        else:
            raise RuntimeError(f'Could not detach {devices[0]} at {location}: {detached.stderr.strip()}')


def write_layout(folder: Path, app_name: str) -> dict:
    try:
        from ds_store import DSStore
    except ImportError as error:
        raise SystemExit('Install build-only ds_store==1.3.3 mac_alias==2.2.3 into .tools/packaging, then run this script with that Python; or explicitly use --plain-layout.') from error
    window = {'WindowBounds': '{{360, 180}, {680, 380}}', 'ShowToolbar': False,
              'ShowStatusBar': False, 'ShowSidebar': False, 'ShowPathbar': False,
              'ShowTabView': False, 'ContainerShowSidebar': False, 'SidebarWidth': 0}
    icons = {'viewOptionsVersion': 1, 'arrangeBy': 'none', 'iconSize': 104.0,
             'textSize': 14.0, 'gridSpacing': 100.0, 'gridOffsetX': 0.0,
             'gridOffsetY': 0.0, 'labelOnBottom': True, 'showIconPreview': True,
             'showItemInfo': False, 'backgroundType': 0,
             'scrollPositionX': 0.0, 'scrollPositionY': 0.0}
    positions = {app_name: (175, 170), 'Applications': (495, 170)}
    with DSStore.open(str(folder / '.DS_Store'), 'w+') as store:
        store['.']['vSrn'] = ('long', 1)
        store['.']['bwsp'] = window
        store['.']['icvp'] = icons
        store['.']['icvl'] = ('type', b'icnv')
        for name, position in positions.items():
            store[name]['Iloc'] = position
    with DSStore.open(str(folder / '.DS_Store'), 'r') as store:
        assert store['.']['bwsp'] == window
        assert store['.']['icvp'] == icons
        assert all(store[name]['Iloc'] == position for name, position in positions.items())
    return {'method': 'direct DS_Store metadata, no Finder automation',
            'ds_store': importlib.metadata.version('ds_store'),
            'mac_alias': importlib.metadata.version('mac_alias'),
            'window': window, 'icons': icons, 'positions': positions,
            'metadata_sha256': sha(folder / '.DS_Store')}


def create(args: argparse.Namespace) -> dict:
    if platform.system() != 'Darwin':
        raise SystemExit('Run this builder on macOS; hdiutil and ditto are built in.')
    if args.app is None:
        signature_record = ROOT / 'build/releases/DynaMol.app.signature.json'
        if not signature_record.exists():
            raise SystemExit('Build the app with package_macos.py first, or provide --app.')
        app = Path(json.loads(signature_record.read_text())['bundle']).resolve()
    else:
        app = args.app.resolve()
    if not app.is_dir() or app.suffix != '.app':
        raise SystemExit('--app must be an existing self-contained .app bundle')
    verify_app_signature(app)
    manifest_path = app / 'Contents/Resources/manifest.json'
    manifest = json.loads(manifest_path.read_text())
    version = manifest['version']
    if not re.fullmatch(r'\d+\.\d+\.\d+(?:[-+][A-Za-z0-9.-]+)?', version):
        raise SystemExit('The app manifest must specify a release version.')
    output = (args.output or ROOT / f'build/releases/DynaMol-{version}-macos-arm64.dmg').resolve()
    if output.suffix.lower() != '.dmg' or output.is_relative_to(app):
        raise SystemExit('--output must be a .dmg outside the input app')
    if output.exists():
        raise SystemExit(f'Output already exists; choose another path to preserve it: {output}')
    build_id = manifest['build_id']
    check_dir = (args.validation_dir or Path(tempfile.mkdtemp(prefix=f'DynaMol DMG Validation {build_id} '))).resolve()
    if check_dir.is_relative_to(app):
        raise SystemExit('Validation directory must be outside the input app')
    installed = check_dir / 'Applications' / app.name
    if installed.exists():
        raise SystemExit(f'Isolated install-copy destination already exists: {installed}')
    output.parent.mkdir(parents=True, exist_ok=True)
    check_dir.mkdir(parents=True, exist_ok=True)
    expected = inventory(app)
    app_inventory = check_dir / 'app-inventory.json'
    app_inventory.write_text(json.dumps(expected, indent=2) + '\n')
    started = time.monotonic()
    report = {'created_at': datetime.now(timezone.utc).isoformat(), 'build_id': build_id,
              'input_app': str(app), 'input_manifest_sha256': sha(manifest_path),
              'script_sha256': sha(Path(__file__)), 'app_inventory': str(app_inventory),
              'app_inventory_sha256': sha(app_inventory), 'app_entries_checked': len(expected),
              'app_tree_sha256': inventory_hash(expected), 'volume_name': args.volume_name,
              'finder_layout': False, 'unsigned': True, 'filesystem': 'APFS',
              'checks': {'input_app_signature_integrity': True}}
    # Keep app copies away from document-sync/Finder surfaces that can attach
    # presentation metadata during signing or copying. Only the finished disk
    # image is published to the requested output directory.
    with tempfile.TemporaryDirectory(prefix='DynaMol DMG ') as temporary:
        temp = Path(temporary).resolve()
        staging = temp / 'Staging'
        staging.mkdir()
        print('Copying the validated app into an isolated installer stage...', flush=True)
        run('/usr/bin/ditto', '--rsrc', '--extattr', app, staging / app.name)
        (staging / 'Applications').symlink_to('/Applications', target_is_directory=True)
        if not args.plain_layout:
            print('Writing the two-icon installer layout...', flush=True)
            report['finder_layout'] = write_layout(staging, app.name)
        assert_inventory(staging / app.name, expected)
        verify_app_signature(staging / app.name)
        compressed = temp / 'installer.dmg'
        print('Creating the compressed read-only disk image...', flush=True)
        run(HDIUTIL, 'create', '-srcfolder', staging, '-volname', args.volume_name,
            '-fs', 'APFS', '-format', 'UDZO', '-imagekey', 'zlib-level=9', '-nospotlight', compressed)
        run(HDIUTIL, 'verify', compressed, capture_output=True)
        image_info = plistlib.loads(run(HDIUTIL, 'imageinfo', '-plist', compressed, capture_output=True).stdout)
        if image_info.get('Format') != 'UDZO':
            raise RuntimeError('Final image is not the requested compressed read-only UDZO format')
        with mounted(compressed, temp / 'Read Only Mount', readonly=True):
            mount = temp / 'Read Only Mount'
            disk_info = plistlib.loads(run('/usr/sbin/diskutil', 'info', '-plist', mount, capture_output=True).stdout)
            if disk_info.get('WritableVolume') is not False or not (os.statvfs(mount).f_flag & os.ST_RDONLY):
                raise RuntimeError('Mounted installer is not read-only')
            try:
                (mount / '.dynamol-write-probe').write_text('must fail')
            except OSError as error:
                if error.errno != errno.EROFS:
                    raise
            else:
                raise RuntimeError('Mounted installer unexpectedly accepted a write')
            if os.readlink(mount / 'Applications') != '/Applications':
                raise RuntimeError('Installer Applications shortcut has the wrong target')
            visible = sorted(p.name for p in mount.iterdir() if not p.name.startswith('.'))
            if visible != sorted([app.name, 'Applications']):
                raise RuntimeError(f'Unexpected visible installer items: {visible}')
            assert_inventory(mount / app.name, expected)
            verify_app_signature(mount / app.name)
            installed.parent.mkdir(parents=True, exist_ok=True)
            print('Testing the drag-style copy into an isolated Applications folder...', flush=True)
            run('/usr/bin/ditto', '--rsrc', '--extattr', mount / app.name, installed)
            assert_inventory(installed, expected)
            verify_app_signature(installed)
            info = plistlib.loads((installed / 'Contents/Info.plist').read_bytes())
            executable = installed / 'Contents/MacOS' / info['CFBundleExecutable']
            if not os.access(executable, os.X_OK):
                raise RuntimeError('Copied app launcher lost its executable permission')
            report['checks'].update(compressed_udzo=True, hdiutil_verify=True,
                mounted_read_only=True, write_rejected_erofs=True, applications_shortcut=True,
                exactly_two_visible_items=True, mounted_app_matches=True,
                isolated_drag_copy_matches=True, copied_launcher_executable=True,
                staged_app_signature_integrity=True, mounted_app_signature_integrity=True,
                installed_app_signature_integrity=True)
        assert_inventory(app, expected)
        report['checks'].update(input_app_unchanged=True, validation_mount_detached=True)
        shutil.move(str(compressed), str(output))
    checksum = output.with_suffix('.dmg.sha256')
    checksum.write_text(f'{sha(output)}  {output.name}\n')
    report.update(passed=True, dmg=str(output), dmg_sha256=sha(output), dmg_bytes=output.stat().st_size,
                  isolated_installed_app=str(installed), elapsed_seconds=round(time.monotonic() - started, 2),
                  limitations=['Container and copy validation only; the script does not launch the app or run simulations.',
                               'No Developer ID signing, notarization, Gatekeeper bypass or /Applications mutation.',
                               'Finder layout is prepared only when not explicitly disabled.'])
    record = output.with_suffix('.dmg.validation.json')
    record.write_text(json.dumps(report, indent=2) + '\n')
    print(json.dumps({'dmg': str(output), 'sha256': report['dmg_sha256'], 'validation': str(record),
                      'isolated_installed_app': str(installed), 'passed': True}, indent=2), flush=True)
    return report


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--app', type=Path, help='Sealed app; defaults to the private path recorded by package_macos.py.')
    parser.add_argument('--output', type=Path, help='DMG destination; defaults to build/releases with the app manifest version.')
    parser.add_argument('--volume-name', default='Drag DynaMol to Applications')
    parser.add_argument('--validation-dir', type=Path)
    parser.add_argument('--plain-layout', action='store_true', help='Omit the build-only ds_store layout; records the plain fallback explicitly.')
    create(parser.parse_args())
