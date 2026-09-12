#!/usr/bin/env python3
"""Copy completed benchmark results into a local DynaMol library, preserving its open workspace."""
import argparse
import hashlib
import json
from pathlib import Path
import shutil
import tempfile
import urllib.request


def digest(path):
    with path.open('rb') as handle:
        return hashlib.file_digest(handle, 'sha256').hexdigest()


def inventory(folder):
    files = {}
    for path in sorted(folder.rglob('*')):
        if path.is_symlink():
            raise ValueError(f'Symlinks are not copied: {path}')
        if path.is_file():
            files[str(path.relative_to(folder))] = digest(path)
    return files


def api(base, path, payload=None):
    request = urllib.request.Request(base + path, data=None if payload is None else json.dumps(payload).encode(),
                                     headers={'Content-Type': 'application/json'})
    with urllib.request.urlopen(request, timeout=30) as response:
        return json.load(response)


def copy_verified(source, target, staging):
    expected = inventory(source)
    if target.exists():
        if inventory(target) != expected:
            raise ValueError(f'Existing library entry differs; nothing was overwritten: {target}')
    else:
        with tempfile.TemporaryDirectory(prefix='pressure-copy-', dir=staging) as temporary:
            copied = Path(temporary) / target.name
            shutil.copytree(source, copied)
            if inventory(copied) != expected:
                raise ValueError(f'Copy verification failed: {source}')
            target.parent.mkdir(parents=True, exist_ok=True)
            copied.rename(target)
    return {'source': str(source), 'target': str(target), 'sha256': expected}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--manifest', type=Path, required=True)
    parser.add_argument('--destination', type=Path, required=True)
    parser.add_argument('--base-url', default='http://127.0.0.1:8765')
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    cases = json.loads(args.manifest.read_text())['cases']
    destination = args.destination.resolve()
    current = destination / 'workspaces/current.json'
    before = digest(current) if current.exists() else None
    if any(j['status'] in {'running', 'queued', 'cancelling'} for j in api(args.base_url, '/api/jobs')):
        raise ValueError('A user computation is active; defer library delivery.')
    for case in cases:
        job = json.loads((Path(case['dataRoot']) / 'jobs' / case['jobId'] / 'status.json').read_text())
        if job['status'] != 'completed' or job['dataset_id'] != case['datasetId']:
            raise ValueError(f'Case is not a completed matching trajectory: {case["name"]}')
    records, copied = [], set()
    for case in cases:
        source = Path(case['dataRoot'])
        pending, ancestors = [case['datasetId']], set()
        while pending:
            dataset_id = pending.pop()
            if dataset_id in ancestors:
                continue
            ancestors.add(dataset_id)
            meta = json.loads((source / 'datasets' / dataset_id / 'metadata.json').read_text())
            if meta.get('parent_dataset_id'):
                pending.append(meta['parent_dataset_id'])
        folders = [('datasets', identifier) for identifier in sorted(ancestors)]
        for path in (source / 'jobs').glob('*/status.json'):
            status = json.loads(path.read_text())
            if status.get('status') == 'completed' and status.get('dataset_id') in ancestors:
                folders.append(('jobs', path.parent.name))
        for kind, identifier in folders:
            key = kind, identifier
            if key not in copied:
                records.append(copy_verified(source / kind / identifier, destination / kind / identifier, destination))
                copied.add(key)
        snapshot = api(args.base_url, f'/api/jobs/{case["jobId"]}/measurements')
        if snapshot['status'] != 'completed' or snapshot['output_dataset_id'] != case['datasetId'] or snapshot['errors']:
            raise ValueError('Published live measurements are not complete and valid.')
        measurements = []
        for index, value in enumerate(snapshot['measurements']):
            measurements.append({**value, 'atoms': value['output_atoms'], 'visible': True,
                                 'trackDuringRun': value.get('trackDuringRun') is not False,
                                 'color': value.get('color') or ['#087f8c', '#5465ca', '#c77828'][index % 3]})
        project_id = 'pressure-' + case['jobId']
        previous = next((p for p in api(args.base_url, '/api/projects') if p['id'] == project_id), None)
        if previous and previous['dataset_id'] != case['datasetId']:
            raise ValueError('Existing project ID belongs to another dataset; nothing was overwritten.')
        project = api(args.base_url, '/api/projects', {
            'id': project_id, 'name': f'Pressure test · {case["pdbId"]} · {case["name"]} · 100 ps',
            'state': {'version': 1, 'dataset_id': case['datasetId'], 'frame': 0,
                      'analysis_expanded': True, 'measurements': measurements,
                      'active_measurement': measurements[0]['id'] if measurements else None}})
        records.append({'project_id': project['id'], 'dataset_id': case['datasetId'], 'name': project['name']})
    after = digest(current) if current.exists() else None
    report = {'base_url': args.base_url, 'workspace_sha256_before': before, 'workspace_sha256_after': after,
              'current_workspace_unchanged': before == after, 'records': records}
    args.output.write_text(json.dumps(report, indent=2) + '\n')
    print(json.dumps({'projects': len(cases), 'copied_entries': len(copied), 'current_workspace_unchanged': before == after}))


if __name__ == '__main__':
    main()
