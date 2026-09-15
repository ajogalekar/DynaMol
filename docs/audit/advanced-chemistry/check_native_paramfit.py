"""Replay unmodified upstream synthetic fitting fixtures, preserving outputs.

These tests qualify a native executable, not parameters for a biological system.
The historical fixture name 'perfect fit' is not an accuracy assertion here.
"""
from __future__ import annotations

import hashlib
import json
import os
from pathlib import Path
import re
import shutil
import subprocess
import time

ROOT = Path(__file__).resolve().parents[3]
EVIDENCE = ROOT / 'build/advanced-chemistry/paramfit-native-audit-v1'
SOURCE = EVIDENCE / 'source/amber24_src/AmberTools/test/paramfit'
CACHE = Path.home() / '.cache/dynamol-native-audits/paramfit-v1'
BINARY = ROOT / '.tools/ambertools/bin/paramfit'


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def table(path):
    return [[float(value) for value in line.split()] for line in Path(path).read_text().splitlines()
            if line.strip() and not line.startswith('#')]


def run():
    CACHE.mkdir(parents=True, exist_ok=False)
    report = {'scope': 'Native implementation regression against upstream synthetic NMA fixtures, not quantum data or a DynaMol force-field model.',
        'executable': str(BINARY), 'executable_sha256': sha(BINARY),
        'runner_sha256': sha(__file__), 'source_manifest_sha256': sha(EVIDENCE / 'source-manifest.json'),
        'maximum_seconds_per_case': 180, 'threads': 1, 'cases': [], 'running': True}
    cases = [
        ('dihedral_least_squares_NMA', 'NMA.prmtop', 'mdcrd', 'amber_energy'),
        ('simplex_perfect_fit_NMA', 'prmtop', 'mdcrd_creation/mdcrd', 'mdcrd_creation/amber_energy.dat'),
    ]
    for name, topology, coordinates, energies in cases:
        target = CACHE / name
        shutil.copytree(SOURCE / name, target)
        arguments = [str(BINARY), '-i', 'Job_Control.in', '-p', topology,
                     '-c', coordinates, '-q', energies, '--random-seed', '5000']
        inputs = {p: sha(target / p) for p in ('Job_Control.in', topology, coordinates, energies)}
        record = {'name': name, 'directory': str(target), 'arguments': arguments,
                  'input_sha256': inputs, 'started_unix': time.time()}
        report['cases'].append(record)
        (EVIDENCE / 'progress.json').write_text(json.dumps(report, indent=2) + '\n')
        environment = dict(os.environ, OMP_NUM_THREADS='1', OPENBLAS_NUM_THREADS='1',
                           VECLIB_MAXIMUM_THREADS='1')
        try:
            with (target / 'native.log').open('w') as log:
                process = subprocess.run(arguments, cwd=target, env=environment,
                    stdin=subprocess.DEVNULL, stdout=log, stderr=subprocess.STDOUT, timeout=180)
            record['returncode'] = process.returncode
            record['input_unchanged'] = all(sha(target / p) == digest for p, digest in inputs.items())
            record['output_sha256'] = {p.name: sha(p) for p in (target / 'native.log', target / 'energy.dat', target / 'frcmod') if p.exists()}
            if process.returncode:
                raise RuntimeError('Native process failed; retained native.log contains details')
            values = table(target / 'energy.dat')
            if not values or any(len(row) != 3 for row in values):
                raise ValueError('Unexpected native energy-table format')
            residual = [row[1] - row[2] for row in values]
            record.update(structures=len(values),
                synthetic_target_energy_rmse_kcal_mol=(sum(v*v for v in residual) / len(residual))**0.5,
                synthetic_target_energy_max_error_kcal_mol=max(map(abs, residual)))
            reference = target / 'saved_output/energy.dat.saved'
            if reference.exists():
                expected = table(reference)
                if len(expected) != len(values):
                    raise ValueError('Upstream reference and native energy row counts differ')
                record['reference_fitted_energy_max_difference_kcal_mol'] = max(abs(a[1]-b[1]) for a,b in zip(values,expected))
                record['reference_target_energy_max_difference_kcal_mol'] = max(abs(a[2]-b[2]) for a,b in zip(values,expected))
            else:
                # The upstream simplex reference contains only the fitted column.
                expected = [row[0] for row in table(target / 'saved_output/energy.out.saved')]
                if len(expected) != len(values):
                    raise ValueError('Upstream fitted-energy row count differs')
                record['reference_fitted_energy_max_difference_kcal_mol'] = max(abs(row[1]-v) for row,v in zip(values,expected))
            if (target / 'frcmod').exists():
                actual = (target / 'frcmod').read_text()
                expected_text = (target / 'saved_output/frcmod.saved').read_text()
                record['frcmod_exact_reference_match'] = actual == expected_text
                pattern = r'(?<![A-Za-z0-9_])[+-]?(?:\d+\.\d*|\.\d+)(?:[eE][+-]?\d+)?'
                a, b = [float(x) for x in re.findall(pattern,actual)], [float(x) for x in re.findall(pattern,expected_text)]
                record['frcmod_numeric_reference_max_difference'] = max(abs(x-y) for x,y in zip(a,b)) if len(a)==len(b) and a else None
            record['status'] = 'completed_for_review'
        except Exception as error:
            record.update(status='failed', error={'type': type(error).__name__, 'message': str(error)})
        record['elapsed_seconds'] = time.time() - record['started_unix']
    report.update(running=False, executable_unchanged=sha(BINARY) == report['executable_sha256'],
                  biological_model_ready=False)
    (EVIDENCE / 'result.json').write_text(json.dumps(report, indent=2) + '\n')
    (EVIDENCE / 'progress.json').write_text(json.dumps(report, indent=2) + '\n')
    print(json.dumps(report, indent=2))


if __name__ == '__main__':
    run()
