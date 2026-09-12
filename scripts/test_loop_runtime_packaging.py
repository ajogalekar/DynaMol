"""Offline controls for the optional manifest runtime and its required ABI pins."""
import hashlib
import importlib.util
import io
import json
from pathlib import Path
import sys
import tarfile

import pytest


ROOT = Path(__file__).resolve().parents[1]


def load(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


bootstrap = load('loop_bootstrap_test', ROOT / 'packaging/macos/bootstrap.py')
package = load('loop_package_test', ROOT / 'scripts/package_macos.py')


def test_private_environment_uses_only_manifest_loop_runtime(tmp_path, monkeypatch):
    monkeypatch.setenv('DYNAMOL_PROMOD3', '/external/promod3')
    monkeypatch.setenv('PYTHONPATH', '/external/python')
    engines = {name: tmp_path / name for name in ('gromacs', 'ambertools', 'promod3')}
    env = bootstrap.private_environment(tmp_path / 'Workspace', engines)
    assert env['DYNAMOL_PROMOD3'] == str(engines['promod3'])
    assert 'PYTHONPATH' not in env
    assert env['PYTHONDONTWRITEBYTECODE'] == '1'
    del engines['promod3']
    assert 'DYNAMOL_PROMOD3' not in bootstrap.private_environment(tmp_path / 'Workspace', engines)


def test_generic_unpack_includes_third_runtime_and_reuses_relocated_marker(tmp_path, monkeypatch):
    resources = tmp_path / 'Resources'
    (resources / 'engines').mkdir(parents=True)
    monkeypatch.setattr(bootstrap, 'RESOURCES', resources)
    monkeypatch.setattr(bootstrap, 'PYTHON', Path(sys.executable))
    records = {}
    for name in ('gromacs', 'ambertools', 'promod3'):
        archive = resources / 'engines' / f'{name}.tar.gz'
        # A tiny relocation script exercises the real extraction/marker path,
        # without a native runtime download or computation.
        content = b'from pathlib import Path\nPath(__file__).parent.joinpath("relocated").write_text("once")\n'
        with tarfile.open(archive, 'w:gz') as packed:
            entry = tarfile.TarInfo('bin/conda-unpack')
            entry.size = len(content)
            packed.addfile(entry, io.BytesIO(content))
        records[name] = {'archive': archive.name, 'sha256': hashlib.sha256(archive.read_bytes()).hexdigest(), 'label': name}
    home = tmp_path / 'Private Home'
    with (tmp_path / 'unpack.log').open('w') as log:
        first = bootstrap.ensure_engines(home, {'engines': records}, log)
        assert set(first) == set(records)
        for name, prefix in first.items():
            assert (prefix / 'bin/relocated').read_text() == 'once'
            marker = json.loads((prefix / '.dynamol-complete.json').read_text())
            assert marker == {'sha256': records[name]['sha256'], 'prefix': str(prefix)}
            # A complete cache remains usable even if the source archive is
            # no longer present; repeat launch must not relocate it again.
            (resources / 'engines' / records[name]['archive']).unlink()
        assert bootstrap.ensure_engines(home, {'engines': records}, log) == first


def test_loop_runtime_rejects_missing_or_incompatible_abi_metadata(tmp_path):
    metadata = tmp_path / 'conda-meta'
    metadata.mkdir()
    with pytest.raises(SystemExit, match='promod3==3.6.0'):
        package.validate_loop_runtime_metadata(tmp_path, '3.6.0')
    for name, version in {'promod3': '3.6.0', 'openstructure': '2.11.1', 'openmm': '8.6.1'}.items():
        (metadata / f'{name}.json').write_text(json.dumps({'name': name, 'version': version}))
    with pytest.raises(SystemExit, match='openmm==8.5.1'):
        package.validate_loop_runtime_metadata(tmp_path, '3.6.0')
    (metadata / 'openmm.json').write_text(json.dumps({'name': 'openmm', 'version': '8.5.1'}))
    package.validate_loop_runtime_metadata(tmp_path, '3.6.0')


def test_new_runtime_supplemental_notice_materials_are_complete_and_hash_verified():
    versions = {'freetype': '2.14.3', 'gcc-runtime': '16.2.0', 'libsqlite': '3.53.4'}
    assert set(package.PROMOD3_NOTICE_FAMILIES.values()) == set(versions)
    for name, version in versions.items():
        folder = ROOT / 'packaging/licenses' / f'{name}-{version}'
        record = json.loads((folder / 'provenance.json').read_text())
        assert (record['name'], record['version']) == (name, version)
        assert record['sources']
        assert 'LICENSE' in record['files_sha256']
        for filename, expected in record['files_sha256'].items():
            assert Path(filename).name == filename
            assert hashlib.sha256((folder / filename).read_bytes()).hexdigest() == expected
    gcc = ROOT / 'packaging/licenses/gcc-runtime-16.2.0'
    assert 'GCC RUNTIME LIBRARY EXCEPTION' in (gcc / 'RUNTIME.LIBRARY.EXCEPTION').read_text()
    assert 'GNU GENERAL PUBLIC LICENSE' in (gcc / 'LICENSE').read_text()
