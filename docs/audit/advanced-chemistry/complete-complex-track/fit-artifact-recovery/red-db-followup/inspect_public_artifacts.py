"""Bounded read-only audit of observed public REDDB DMP project artifacts.

The two public search POSTs are preserved separately in the cache. This script
downloads only observed summary/archive links for their exact DMP hits. It never
executes downloaded scripts or assigns any parameter to a molecular system.
"""
from concurrent.futures import ThreadPoolExecutor
from datetime import datetime, timezone
from pathlib import Path
import hashlib
import html
import io
import json
import re
import tarfile
import urllib.request

CACHE = Path('/Users/ashujo/.cache/dynamol-research/complete-complex-track/fit-artifact-recovery-red-db-v1')
BASE = 'https://upjv.q4md-forcefieldtools.org/REDDB/'
OUT = Path(__file__).resolve().parent


def sha(data):
    return hashlib.sha256(data).hexdigest()


def clean(s):
    return ' '.join(html.unescape(re.sub('<[^>]*>', ' ', s)).split())


def fetch(item):
    project, url, filename = item
    destination = CACHE / project / filename
    destination.parent.mkdir(parents=True, exist_ok=True)
    meta_path = destination.with_name(filename + '.retrieval.json')
    if meta_path.exists():
        return json.loads(meta_path.read_text())
    metadata = {'project': project, 'url': url, 'path': str(destination),
                'requested_utc': datetime.now(timezone.utc).isoformat(),
                'tls_verification': 'default verified urllib context',
                'request_method': 'GET', 'source_link': 'lookup-2.html'}
    try:
        request = urllib.request.Request(url, headers={'User-Agent': 'DynaMol-public-source-audit/1.0'})
        with urllib.request.urlopen(request, timeout=25) as response:
            data = response.read(5_000_001)
            metadata.update(status=response.status, final_url=response.url,
                            content_type=response.headers.get('Content-Type'),
                            last_modified=response.headers.get('Last-Modified'),
                            bytes=len(data), sha256=sha(data))
        if len(data) > 5_000_000:
            raise ValueError('bounded download size exceeded')
        destination.write_bytes(data)
    except Exception as exc:
        metadata['error'] = repr(exc)
    meta_path.write_text(json.dumps(metadata, indent=2) + '\n')
    return metadata


def mol2_inventory(text):
    sections = {}
    current = None
    for line in text.splitlines():
        if line.startswith('@<TRIPOS>'):
            current = line[9:]
            sections[current] = []
        elif current and line.strip():
            sections[current].append(line)
    atoms = []
    for line in sections.get('ATOM', []):
        fields = line.split()
        atoms.append({'id': int(fields[0]), 'name': fields[1],
                      'coordinates_A': [float(v) for v in fields[2:5]],
                      'type': fields[5], 'residue_id': int(fields[6]),
                      'residue_name': fields[7], 'charge_e': float(fields[8])})
    bonds = [{'id': int(f[0]), 'a': int(f[1]), 'b': int(f[2]), 'type': f[3]}
             for f in (line.split() for line in sections.get('BOND', []))]
    return {'atoms': atoms, 'bonds': bonds, 'atom_count': len(atoms),
            'bond_count': len(bonds), 'total_charge_e': round(sum(a['charge_e'] for a in atoms), 8)}


def main():
    listing = (CACHE / 'lookup-2.html').read_text()
    tasks = []
    for number in range(11, 15):
        project = f'W-{number}'
        for relative, filename in [(f'projects/{project}/', 'summary.html'),
                                   (f'projects/{project}/{project}.tar.bz2', f'{project}.tar.bz2')]:
            assert relative in re.findall(r'href=[\"\x27]([^\"\x27]*)', listing)
            tasks.append((project, BASE + relative, filename))
    with ThreadPoolExecutor(max_workers=2) as pool:
        retrievals = list(pool.map(fetch, tasks))
    projects = []
    for number in range(11, 15):
        project = f'W-{number}'
        directory = CACHE / project
        record = {'project': project, 'summary_url': BASE + f'projects/{project}/',
                  'retrievals': [r for r in retrievals if r['project'] == project]}
        summary = directory / 'summary.html'
        if summary.exists():
            summary_text = clean(summary.read_text(errors='replace'))
            (directory / 'summary.txt').write_text(summary_text + '\n')
            record['summary_text_cache'] = str(directory / 'summary.txt')
            record['summary_text_sha256'] = sha((summary_text + '\n').encode())
        archive = directory / f'{project}.tar.bz2'
        if archive.exists():
            record['archive_members'] = []
            try:
                with tarfile.open(fileobj=io.BytesIO(archive.read_bytes()), mode='r:bz2') as tar:
                    for member in tar.getmembers():
                        if not member.isfile():
                            record['archive_members'].append({'name': member.name, 'regular_file': False})
                            continue
                        if member.size > 2_000_000:
                            raise ValueError('member exceeds bounded inspection size')
                        data = tar.extractfile(member).read()
                        item = {'name': member.name, 'bytes': len(data), 'sha256': sha(data),
                                'regular_file': True, 'archive_mtime_utc': datetime.fromtimestamp(member.mtime, timezone.utc).isoformat()}
                        record['archive_members'].append(item)
                        # Copy flat, allowlisted names for reading only; no archive extraction.
                        name = Path(member.name).name
                        if name in {'tripos1.mol2', 'mol1.pdb', 'input1.in', 'input2.in',
                                    'script1.ff', 'script2.ff', 'script3.ff', 'script4.ff', 'script5.ff',
                                    'summary.html', 'index.html'}:
                            target = directory / 'inspected-members' / name
                            target.parent.mkdir(exist_ok=True)
                            target.write_bytes(data)
                            item['inspected_copy'] = str(target)
                        if name == 'tripos1.mol2':
                            record['mol2'] = mol2_inventory(data.decode())
            except Exception as exc:
                record['archive_error'] = repr(exc)
        projects.append(record)
    result = {'schema': 'dynamol.reddb-public-dmp-artifacts.v1',
              'generated_utc': datetime.now(timezone.utc).isoformat(),
              'scope': 'source artifact recovery only; no assignment or model validation',
              'public_database_lookups_used': 2, 'public_database_lookup_limit': 3,
              'lookups': [json.loads((CACHE / f'lookup-{i}.json').read_text()) for i in (1, 2)],
              'projects': projects}
    (OUT / 'artifact-inventory.json').write_text(json.dumps(result, indent=2) + '\n')
    print(json.dumps([{'project': p['project'], 'summary_text_cache': p.get('summary_text_cache'),
                       'members': p.get('archive_members'), 'mol2': p.get('mol2'),
                       'error': p.get('archive_error')} for p in projects], indent=2))


if __name__ == '__main__':
    main()
