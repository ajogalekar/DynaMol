"""Native regression for the incomplete app signature shipped in 0.1.0."""
import argparse
import importlib.util
from pathlib import Path
import platform
import plistlib
import subprocess

import pytest


spec = importlib.util.spec_from_file_location('dynamol_dmg', Path(__file__).with_name('create_macos_dmg.py'))
dmg = importlib.util.module_from_spec(spec)
spec.loader.exec_module(dmg)
package_spec = importlib.util.spec_from_file_location('dynamol_package', Path(__file__).with_name('package_macos.py'))
package = importlib.util.module_from_spec(package_spec)
package_spec.loader.exec_module(package)

pytestmark = pytest.mark.skipif(platform.system() != 'Darwin', reason='macOS codesign regression')


def test_distribution_rejects_launcher_only_and_tampered_bundle(tmp_path):
    app = tmp_path / 'Signature Fixture.app'
    macos = app / 'Contents/MacOS'
    resources = app / 'Contents/Resources'
    macos.mkdir(parents=True)
    resources.mkdir()
    source = tmp_path / 'main.c'
    source.write_text('int main(void) { return 0; }\n')
    executable = macos / 'Fixture'
    subprocess.run(['/usr/bin/clang', str(source), '-o', str(executable)], check=True, capture_output=True)
    (app / 'Contents/Info.plist').write_bytes(plistlib.dumps({
        'CFBundleIdentifier': 'org.dynamol.signature-test',
        'CFBundleExecutable': 'Fixture',
        'CFBundlePackageType': 'APPL',
    }))
    marker = resources / 'sealed.txt'
    marker.write_text('original\n')

    # Reproduce the old recipe: sign an executable before it enters a bundle.
    bare = tmp_path / 'BareLauncher'
    bare.write_bytes(executable.read_bytes())
    bare.chmod(0o755)
    subprocess.run(['/usr/bin/codesign', '--force', '--sign', '-', str(bare)], check=True, capture_output=True)
    executable.write_bytes(bare.read_bytes())
    with pytest.raises(ValueError, match='signature integrity failed'):
        dmg.verify_app_signature(app)
    output = tmp_path / 'rejected.dmg'
    with pytest.raises(ValueError, match='signature integrity failed'):
        dmg.create(argparse.Namespace(app=app, output=output))
    assert not output.exists()

    # Use the actual packager: preserve valid nested code and quarantine, while
    # clearing only Finder presentation metadata that codesign disallows.
    inner = resources / 'runtime-tool'
    inner.write_bytes(bare.read_bytes())
    inner.chmod(0o755)
    original_inner = inner.read_bytes()
    subprocess.run(['/usr/bin/xattr', '-wx', 'com.apple.FinderInfo', (b'APPL' + bytes(28)).hex(), str(app)], check=True)
    quarantine = '0081;65000000;DynaMolSignatureRegression;'
    subprocess.run(['/usr/bin/xattr', '-w', 'com.apple.quarantine', quarantine, str(app)], check=True)
    signed = package.seal_application(app)
    assert signed['strict_deep_verification_passed']
    assert not signed['inner_signatures_repaired']
    assert inner.read_bytes() == original_inner
    assert 'com.apple.FinderInfo' not in subprocess.check_output(['/usr/bin/xattr', str(app)], text=True).splitlines()
    assert subprocess.check_output(['/usr/bin/xattr', '-p', 'com.apple.quarantine', str(app)], text=True).strip() == quarantine
    dmg.verify_app_signature(app)
    marker.write_text('changed after signing\n')
    with pytest.raises(ValueError, match='signature integrity failed'):
        dmg.verify_app_signature(app)
