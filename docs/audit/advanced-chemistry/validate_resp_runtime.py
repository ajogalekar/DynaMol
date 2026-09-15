"""Repeat the current native 1OKL RESP fit with the storage-only rebuild.

This compares executable behavior for identical inputs, not physical model
quality or agreement with the older tutorial's separately parameterized model.
"""
from pathlib import Path
import argparse
import hashlib
import json
import time

import numpy as np

from mcpb_array_adapter import run_native_resp


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--executable', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--repeats', type=int, default=10)
    args = parser.parse_args()
    base = Path(__file__).resolve().parent
    reference = base / 'reference-mcpb-1okl/native-regenerated'
    output = args.output.resolve()
    output.mkdir(parents=True, exist_ok=False)
    atom_ids = (reference / '1OKL_large.fingerprint').read_text().splitlines()
    native_charges = {step: np.asarray([float(s) for s in (reference / f'direct-resp{step}.chg').read_text().split()])
                      for step in (1, 2)}
    report = {'passed': False, 'executable': str(args.executable.resolve()),
              'executable_sha256': hashlib.sha256(args.executable.read_bytes()).hexdigest(),
              'reference': str(reference), 'repeats': [],
              'scope': 'Storage-only rebuild versus current native RESP on identical 113-atom two-stage inputs. Neither historical tutorial charge equivalence nor new-site physical accuracy is asserted.'}
    for repeat in range(args.repeats):
        started = time.monotonic()
        folder = output / f'repeat-{repeat + 1:02d}'
        row = {'repeat': repeat + 1, 'passed': False}
        try:
            fit = run_native_resp(executable=args.executable, workdir=folder,
                                  stage1_input=reference / 'resp1.in', stage2_input=reference / 'resp2.in',
                                  esp_file=reference / '1OKL_large_mk.esp', atom_ids=atom_ids, total_charge=1)
            comparisons = []
            for step in (1, 2):
                actual = np.asarray([float(s) for s in (folder / f'resp{step}.chg').read_text().split()])
                maximum = float(np.max(np.abs(actual - native_charges[step])))
                comparisons.append({'stage': step, 'atom_count': len(actual), 'maximum_charge_difference_e': maximum,
                                    'printed_charges_exactly_equal': bool(np.array_equal(actual, native_charges[step]))})
            row.update(fit=fit, comparisons=comparisons, passed=all(c['printed_charges_exactly_equal'] for c in comparisons))
        except Exception as exc:
            row['error'] = str(exc)
        row['elapsed_seconds'] = time.monotonic() - started
        report['repeats'].append(row)
        report['passed'] = len(report['repeats']) == args.repeats and all(r['passed'] for r in report['repeats'])
        (output / 'validation.json').write_text(json.dumps(report, indent=2) + '\n')
        print(repeat + 1, row['passed'], round(row['elapsed_seconds'], 2), row.get('error', ''), flush=True)
    if not report['passed']:
        raise SystemExit(1)


if __name__ == '__main__':
    main()
