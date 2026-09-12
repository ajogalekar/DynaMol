#!/usr/bin/env python3
"""Collect source evidence for a frozen DynaMol bundle; never change the app ZIP.

Requires PyYAML (present in the DynaMol build environment). This collects
materials, not a legal certification or a claim of bit-identical rebuilding.
"""
from __future__ import annotations

import argparse
from concurrent.futures import ThreadPoolExecutor, as_completed
import csv
import datetime as dt
import hashlib
import json
from pathlib import Path
import re
import shutil
import sys
import tarfile
import time
import tomllib
import urllib.parse
import urllib.request

import yaml

ROOT = Path(__file__).resolve().parents[1]
DEFAULT_BUNDLE = ROOT / 'build/releases/DynaMol.app'
DEFAULT_OUTPUT = ROOT / 'build/releases/source-materials'
DIRECT_PYTHON = {'parmed', 'mdtraj', 'gemmi', 'certifi'}
PERMISSIVE_ALTERNATIVES = {'freetype', 'libfreetype', 'libfreetype6', 'libev', 'perl'}


def sha(path, algorithm='sha256'):
    with Path(path).open('rb') as stream:
        return hashlib.file_digest(stream, algorithm).hexdigest()


def dump(path, value):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_suffix(path.suffix + '.tmp')
    temporary.write_text(json.dumps(value, indent=2, sort_keys=True) + '\n')
    temporary.replace(path)


def sources_in(value):
    """Rendered conda/rattler recipes may place source under staging caches."""
    found = []
    if isinstance(value, dict):
        for key, child in value.items():
            if key == 'source':
                for item in child if isinstance(child, list) else [child]:
                    if isinstance(item, dict) and ('url' in item or 'git_url' in item):
                        found.append(item)
            elif isinstance(child, (dict, list)):
                found.extend(sources_in(child))
    elif isinstance(value, list):
        for child in value:
            found.extend(sources_in(child))
    return found


def pinned_native_git_source(source, name, engine, manifest, resources, target):
    """Resolve a recipe Git source without substituting a mutable release tag."""
    revision = str(source.get('git_rev') or '')
    authority = 'full Git commit in shipped rendered Conda recipe'
    if name == 'openmm':
        runtime = manifest['engines'][engine]
        archive = resources / 'engines' / runtime['archive']
        if sha(archive) != runtime['sha256']:
            raise ValueError('Native OpenMM runtime archive does not match its manifest')
        version_bytes = None
        # Read the exact bundled archive, not a development prefix. Stream to
        # the small version record without extracting or modifying the runtime.
        with tarfile.open(archive, 'r|gz') as packed:
            for member in packed:
                if re.fullmatch(r'lib/python[0-9.]+/site-packages/openmm/version\.py', member.name):
                    if not member.isfile() or member.size > 64 * 1024:
                        raise ValueError('Invalid bundled OpenMM version record')
                    version_bytes = packed.extractfile(member).read()
                    break
        if version_bytes is None:
            raise ValueError('Bundled native OpenMM version record is missing')
        match = re.search(r"git_revision = ['\"]([0-9a-f]{40})", version_bytes.decode())
        if not match or not re.fullmatch(r'[0-9a-f]{7,40}', revision) or not match.group(1).startswith(revision):
            raise ValueError('Native OpenMM embedded source revision conflicts with its recipe')
        revision = match.group(1)
        (target / 'shipped-version.py').write_bytes(version_bytes)
        authority = 'full Git revision embedded in manifest-verified native runtime archive; matches shipped recipe'
    if not re.fullmatch(r'[0-9a-f]{40}', revision):
        raise ValueError('Git source requires a full recorded commit, not a branch or tag')
    repository = re.fullmatch(r'https://github\.com/([A-Za-z0-9_.-]+/[A-Za-z0-9_.-]+?)(?:\.git)?', source['git_url'])
    if not repository:
        raise ValueError('Git source host requires an explicit pinned archive mapping')
    return ({**source, 'url': f'https://codeload.github.com/{repository.group(1)}/tar.gz/{revision}',
             'git_revision': revision, 'filename': f'{name}-{revision}.tar.gz'}, authority)


def plan(bundle, output):
    resources = bundle / 'Contents/Resources'
    notices = resources / 'Notices'
    manifest_path = resources / 'manifest.json'
    inventory_path = notices / 'DEPENDENCY_INVENTORY.json'
    manifest = json.loads(manifest_path.read_text())
    inventory = json.loads(inventory_path.read_text())
    lock_path = resources / 'app/uv.lock'
    locked = {p['name'].lower().replace('_', '-'): p
              for p in tomllib.loads(lock_path.read_text())['package']}
    result = {'schema_version': 1, 'build_id': manifest['build_id'],
              'created_utc': dt.datetime.now(dt.timezone.utc).isoformat(),
              'application_version': manifest['version'],
              'inputs': {'bundle_manifest_sha256': sha(manifest_path),
                         'dependency_inventory_sha256': sha(inventory_path),
                         'uv_lock_sha256': sha(lock_path),
                         'engines': manifest['engines']},
              'packages': [], 'artifacts': {}, 'excluded': [], 'gaps': [],
              'scope': 'Source/material collection for selected shipped copyleft components. '
                       'No legal certification or bit-identical rebuild claim. Alternate '
                       'licenses and compiler-runtime exceptions are recorded separately.'}

    def artifact(source, owner, authority):
        urls = source['url'] if isinstance(source['url'], list) else [source['url']]
        expected = source.get('sha256')
        expected512 = source.get('sha512')
        revision = source.get('git_revision')
        if not expected and not expected512 and not revision:
            raise ValueError(f'Unpinned source for {owner}: {urls}')
        if any('${{' in u or '{{' in u for u in urls):
            raise ValueError(f'Unrendered source URL for {owner}')
        key = expected or expected512 or revision
        filename = source.get('filename') or Path(urllib.parse.urlparse(urls[0]).path).name
        entry = result['artifacts'].setdefault(key, {
            'id': key, 'urls': urls, 'expected_sha256': expected,
            'expected_sha512': expected512,
            'git_revision': revision, 'filename': filename,
            'path': f'archives/{key[:12]}-{filename}',
            'authority': authority, 'status': 'pending', 'used_by': []})
        if owner not in entry['used_by']:
            entry['used_by'].append(owner)
        return key

    for package in inventory['native']:
        license_text = str(package.get('license') or '')
        if not re.search(r'\b(?:L?GPL|MPL)[-\d]', license_text):
            continue
        name, version, build = package['name'], package['version'], package['build']
        owner = f"native:{package['engine']}:{name}:{version}:{build}"
        if name in PERMISSIVE_ALTERNATIVES:
            result['excluded'].append({'package': owner, 'license': license_text,
                'reason': 'Metadata offers a non-GPL alternative; source not collected '
                          'solely because GPL appears in an OR expression. Notices remain in the app.'})
            continue
        source_dir = notices / package['engine'] / f'{name}-{version}-{build}'
        material_path = Path('materials/native') / package['engine'] / source_dir.name
        target = output / material_path
        shutil.copytree(source_dir, target, dirs_exist_ok=True)
        recipe_dir = target / 'recipe'
        recipe_path = recipe_dir / 'rendered_recipe.yaml'
        if not recipe_path.exists():
            recipe_path = recipe_dir / 'meta.yaml'
        record = {'id': owner, 'kind': 'native', **package,
                  'materials_path': str(material_path), 'source_artifacts': [],
                  'selection': 'compiler-runtime exception; collected conservatively'
                  if 'WITH GCC-exception' in license_text else 'explicit copyleft package license'}
        try:
            recipe = yaml.safe_load(recipe_path.read_text())
            specs = sources_in(recipe)
            if not specs:
                raise ValueError('No pinned source in rendered recipe')
            record['source_recipe'] = str(recipe_path.relative_to(output))
            record['source_specs'] = specs
            for spec in specs:
                authority = 'shipped rendered Conda recipe'
                if 'git_url' in spec:
                    spec, authority = pinned_native_git_source(spec, name, package['engine'], manifest, resources, target)
                    record['source_revision'] = spec['git_revision']
                record['source_artifacts'].append(artifact(spec, owner, authority))
                for patch in spec.get('patches', []):
                    if not isinstance(patch, str):
                        raise ValueError('Unrendered conditional patch')
                    candidates = [recipe_dir / patch, recipe_dir / 'parent' / patch]
                    if not any(p.is_file() for p in candidates):
                        result['gaps'].append({'package': owner, 'issue': f'Referenced recipe patch is missing: {patch}'})
            if name == 'ambertools':
                record['upstream_update_level'] = int(version.split('.')[1])
                updates_path = ROOT / 'packaging/source-materials' / f'ambertools-{version}-updates.json'
                if updates_path.is_file():
                    updates = json.loads(updates_path.read_text())
                    expected_numbers = list(range(1, record['upstream_update_level'] + 1))
                    if [u['number'] for u in updates['updates']] != expected_numbers:
                        raise ValueError('Incomplete AmberTools numbered update series')
                    record['upstream_updates'] = updates
                    for update in updates['updates']:
                        record['source_artifacts'].append(artifact(
                            {'url': update['url'], 'sha256': update['sha256'],
                             'filename': f"AmberTools24-update.{update['number']}"},
                            owner, 'official numbered upstream update, checksum frozen by source-materials metadata'))
                else:
                    result['gaps'].append({'package': owner, 'issue': 'Numbered AmberTools upstream update manifest is missing.', 'code': 'amber_updates'})
        except Exception as exc:
            record['error'] = str(exc)
            result['gaps'].append({'package': owner, 'issue': str(exc)})
        result['packages'].append(record)
    site = resources / 'python/lib/python3.12/site-packages'
    # The SciPy wheel's permissive project license does not identify the source
    # revision of every vendored compiler runtime. Preserve the actual evidence
    # and exception text; do not silently equate compiler and runtime versions.
    scipy = next((p for p in inventory['python'] if p['name'].lower() == 'scipy'), None)
    if scipy:
        target = output / 'materials/python' / f"scipy-{scipy['version']}-vendored-runtimes"
        target.mkdir(parents=True, exist_ok=True)
        dist = next(site.glob(f"scipy-{scipy['version']}.dist-info"))
        shutil.copy2(dist / 'LICENSE.txt', target / 'LICENSE.txt')
        shutil.copy2(site / 'scipy/__config__.py', target / 'shipped-build-config.py')
        compiler = re.search(r'"fortran":\s*\{.*?"version":\s*"([^"]+)"',
                             (site / 'scipy/__config__.py').read_text(), re.S)
        files = [p for p in (site / 'scipy/.dylibs').glob('*.dylib')
                 if any(term in p.name for term in ('gfortran', 'quadmath', 'gcc_s'))]
        result['conditional_components'] = [{
            'id': f"python:scipy:{scipy['version']}:vendored-compiler-runtimes",
            'materials_path': str(target.relative_to(output)),
            'files': [{'path': str(p.relative_to(resources)), 'sha256': sha(p),
                       'bytes': p.stat().st_size} for p in files],
            'reported_fortran_compiler_version': compiler.group(1) if compiler else None,
            'evidence': 'Shipped SciPy build configuration does not independently establish the exact source revision of '
                        'the three vendored runtime libraries.',
            'licenses': 'GCC runtime exception text and LGPL libquadmath notice retained verbatim in LICENSE.txt.',
            'status': 'runtime source provenance not independently established; '
                      'compiler-runtime exception applicability is not adjudicated by this collector'}]
        supplement_path = ROOT / 'packaging/source-materials/gcc-13.4.0-supplement.json'
        if compiler and compiler.group(1) == '13.4.0' and supplement_path.is_file():
            supplement = json.loads(supplement_path.read_text())
            component = result['conditional_components'][0]
            component['conservative_source_supplement'] = artifact(
                supplement['source'], component['id'], supplement['authority'])
            component['supplement_limit'] = supplement['limit']

    for package in inventory['python']:
        name = package['name'].lower().replace('_', '-')
        if name not in DIRECT_PYTHON and name != 'openmm':
            continue
        version = package['version']
        owner = f'python:{name}:{version}'
        record = {'id': owner, 'kind': 'python', 'name': name, 'version': version,
                  'license': package.get('license_expression') or package.get('license'),
                  'source_artifacts': [], 'binary_files': []}
        distributions = [p for p in site.glob('*.dist-info')
                         if p.name.split('-')[0].lower().replace('_', '-') == name]
        if len(distributions) != 1:
            raise ValueError(f'Cannot identify exactly one shipped Python distribution: {name}')
        dist = distributions[0]
        target = output / 'materials/python' / f'{name}-{version}'
        shutil.copytree(dist, target / dist.name, dirs_exist_ok=True)
        record['materials_path'] = str(target.relative_to(output))
        for path, _, _ in csv.reader((dist / 'RECORD').open()):
            source = (site / path).resolve()
            if source.is_relative_to(resources.resolve()) and source.is_file() and source.suffix in {'.so', '.dylib'}:
                record['binary_files'].append({'path': str(source.relative_to(resources)), 'sha256': sha(source), 'bytes': source.stat().st_size})
        if name == 'openmm':
            version_file = site / 'openmm/version.py'
            revision = re.search(r"git_revision = ['\"]([0-9a-f]{40})", version_file.read_text()).group(1)
            shutil.copy2(version_file, target / 'shipped-version.py')
            record['source_revision'] = revision
            record['license_scope'] = 'MIT core/CPU/application; OpenCL/HIP/CUDA platform implementations use LGPL. Collecting exact reported native source revision because OpenCL plugins are shipped.'
            spec = {'url': f'https://codeload.github.com/openmm/openmm/tar.gz/{revision}',
                    'git_revision': revision, 'filename': f'openmm-{revision}.tar.gz'}
            record['source_artifacts'].append(artifact(spec, owner, 'git revision embedded in shipped openmm/version.py'))
        else:
            lock = locked.get(name)
            if not lock or lock['version'] != version or not lock.get('sdist'):
                raise ValueError(f'No exact locked source distribution for {owner}')
            source = lock['sdist']
            spec = {'url': source['url'], 'sha256': source['hash'].removeprefix('sha256:')}
            record['source_artifacts'].append(artifact(spec, owner, 'uv.lock bundled with this exact app'))
            if name == 'parmed':
                record['local_build_note'] = manifest['python'].get('packaged_overrides', {}).get('parmed')
        result['packages'].append(record)
    for path in [manifest_path, inventory_path, lock_path]:
        shutil.copy2(path, output / path.name)
    return result


def download(entry, output):
    target = output / entry['path']
    target.parent.mkdir(parents=True, exist_ok=True)
    if target.is_file():
        actual = sha(target)
        expected = entry.get('expected_sha256') or entry.get('sha256')
        sha512_matches = (sha(target, 'sha512') == entry['expected_sha512']) if entry.get('expected_sha512') else None
        if ((expected and actual == expected and sha512_matches is not False)
                or (not expected and sha512_matches is True)):
            return {**entry, 'status': 'verified', 'sha256': actual, 'bytes': target.stat().st_size}
        if expected or entry.get('expected_sha512'):
            raise ValueError(f'Existing source file checksum mismatch: {target.name}')
    errors = []
    urls = list(entry['urls'])
    for url in list(urls):
        if url.startswith('http://'):
            urls.insert(0, 'https://' + url[7:])
    for url in dict.fromkeys(urls):
        partial = target.with_suffix(target.suffix + '.part')
        try:
            request = urllib.request.Request(url, headers={'User-Agent': 'DynaMol-source-materials/0.1'})
            started = time.monotonic()
            with urllib.request.urlopen(request, timeout=45) as response, partial.open('wb') as stream:
                size = 0
                while chunk := response.read(1024 * 1024):
                    size += len(chunk)
                    if size > 1024**3 or time.monotonic() - started > 900:
                        raise ValueError('Single source download exceeded 1 GiB / 15 minute bound')
                    stream.write(chunk)
                final_url = response.url
            actual = sha(partial)
            if entry.get('expected_sha256') and actual != entry['expected_sha256']:
                raise ValueError(f'SHA-256 mismatch: expected {entry["expected_sha256"]}, received {actual}')
            if entry.get('expected_sha512') and sha(partial, 'sha512') != entry['expected_sha512']:
                raise ValueError('SHA-512 mismatch against the recorded official upstream checksum')
            partial.replace(target)
            return {**entry, 'status': 'verified', 'sha256': actual, 'bytes': size,
                    'download_url': final_url,
                    'hash_basis': 'matched pre-recorded source hash' if entry.get('expected_sha256') or entry.get('expected_sha512') else 'archive pinned by embedded source commit; archive hash recorded on retrieval'}
        except Exception as exc:
            partial.unlink(missing_ok=True)
            errors.append({'url': url, 'error': str(exc)})
    return {**entry, 'status': 'unavailable', 'errors': errors}


def verify(index, output):
    failures = []
    checks = 0
    for entry in index['artifacts'].values():
        path = (output / entry['path']).resolve()
        expected = entry.get('sha256')
        checks += 1
        if (not path.is_relative_to(output.resolve()) or not expected
                or not path.is_file() or sha(path) != expected
                or (entry.get('expected_sha256') and expected != entry['expected_sha256'])
                or (entry.get('expected_sha512') and sha(path, 'sha512') != entry['expected_sha512'])):
            failures.append(entry['path'])
    for relative, entry in index.get('material_files', {}).items():
        path = (output / relative).resolve()
        checks += 1
        if not path.is_relative_to(output.resolve()) or not path.is_file() or sha(path) != entry['sha256']:
            failures.append(relative)
    return {'passed': not failures, 'checks': checks, 'failures': failures}


def collect_supporting_materials(output):
    """Keep instructions and the collector with the separate source companion."""
    metadata = ROOT / 'packaging/source-materials'
    shutil.copy2(metadata / 'README.md', output / 'README.md')
    target = output / 'materials/collector'
    target.mkdir(parents=True, exist_ok=True)
    shutil.copy2(Path(__file__), target / 'collect_release_sources.py')
    shutil.copy2(ROOT / 'scripts/test_release_sources.py', target / 'test_release_sources.py')
    for path in metadata.iterdir():
        if path.is_file() and path.name not in {'INDEX.json', 'VERIFICATION.json'}:
            shutil.copy2(path, target / path.name)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--bundle', type=Path, default=DEFAULT_BUNDLE)
    parser.add_argument('--output', type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument('--plan-only', action='store_true')
    parser.add_argument('--resume', action='store_true')
    parser.add_argument('--verify-only', action='store_true')
    parser.add_argument('--workers', type=int, default=3, choices=range(1, 5))
    args = parser.parse_args()
    output = args.output.resolve()
    output.mkdir(parents=True, exist_ok=True)
    index_path = output / 'INDEX.json'
    if args.verify_only:
        result = verify(json.loads(index_path.read_text()), output)
        print(json.dumps(result), flush=True)
        return 0 if result['passed'] else 1
    index = json.loads(index_path.read_text()) if args.resume else plan(args.bundle.resolve(), output)
    if args.resume and index['inputs']['bundle_manifest_sha256'] != sha(args.bundle / 'Contents/Resources/manifest.json'):
        raise ValueError('Resume bundle does not match the indexed release')
    dump(index_path, index)
    print(json.dumps({'build_id': index['build_id'], 'packages': len(index['packages']),
                      'archives': len(index['artifacts']), 'gaps': index['gaps']}), flush=True)
    if args.plan_only:
        return 0
    with ThreadPoolExecutor(max_workers=args.workers) as pool:
        futures = {pool.submit(download, entry, output): key for key, entry in index['artifacts'].items()}
        for future in as_completed(futures):
            key = futures[future]
            try:
                result = future.result()
            except Exception as exc:
                result = {**index['artifacts'][key], 'status': 'unavailable', 'error': str(exc)}
            index['artifacts'][key] = result
            dump(index_path, index)
            print(json.dumps({'file': result['filename'], 'status': result['status'],
                              'bytes': result.get('bytes'), 'errors': result.get('errors', result.get('error'))}), flush=True)
    index['download_summary'] = {'verified': sum(e['status'] == 'verified' for e in index['artifacts'].values()),
                                 'total': len(index['artifacts']),
                                 'bytes': sum(e.get('bytes', 0) for e in index['artifacts'].values())}
    collect_supporting_materials(output)
    index['material_files'] = {str(p.relative_to(output)): {'sha256': sha(p), 'bytes': p.stat().st_size}
                               for p in sorted((output / 'materials').rglob('*')) if p.is_file()}
    for name in ('README.md', 'manifest.json', 'DEPENDENCY_INVENTORY.json', 'uv.lock'):
        path = output / name
        index['material_files'][name] = {'sha256': sha(path), 'bytes': path.stat().st_size}
    dump(index_path, index)
    print(json.dumps(index['download_summary']), flush=True)
    return 0 if index['download_summary']['verified'] == index['download_summary']['total'] and not index['gaps'] else 1


if __name__ == '__main__':
    raise SystemExit(main())
