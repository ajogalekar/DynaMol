"""Freeze and supervise the source-boundary conditioning intervention."""
import argparse
import hashlib
import json
from pathlib import Path
import time

from run_native_comparison import NATIVE, CACHE, supervise, save, sha


def run(output):
    output.mkdir(parents=True, exist_ok=False)
    implementation = output / 'implementation'
    implementation.mkdir()
    for path in [Path(__file__), Path(__file__).with_name('native_mc_adapter.py'), Path(__file__).with_name('run_native_comparison.py')]:
        (implementation / path.name).write_bytes(path.read_bytes())
    original = CACHE / 'hybrid-native-mc-v1/plan.json'
    old = json.loads(original.read_text())
    assert all(sha(Path(p)) == h for p, h in old['source_sha256'].items())
    cases = old['cases']
    for case in cases:
        case['spans'] = [case['spans'][0]]
        case['seeds'] = [2026] if case['case'] == 'masked-control' else [2026, 2027]
    protocols = [
        {'name': 'torsion-source-dirty', 'mode': 'torsion_mc', 'boundary_protocol': 'source_dirty'},
        {'name': 'torsion-source-ccd', 'mode': 'torsion_mc', 'boundary_protocol': 'source_ccd'},
        {'name': 'fragment-source-ccd', 'mode': 'fragment_mc', 'boundary_protocol': 'source_ccd'},
    ]
    plan = {'cases': cases, 'seeds': [2026, 2027], 'steps': 5000, 'protocols': protocols,
            'source_plan': str(original), 'source_plan_sha256': sha(original),
            'source_sha256': old['source_sha256'],
            'implementation_sha256': {str(p): sha(p) for p in implementation.iterdir()},
            'maximum_native_attempts': 15, 'overall_native_seconds': 900,
            'intervention': 'Supply observed flanking phi/psi and neighbor identities to the native torsion sampler, then separately test the native torsion-aware CCD closer. Compare prior same-input/seed native-default artifacts without overwriting or rerunning them.',
            'invariants': 'Native MC steps, scoring weights, cooling, seeds, original stems, source coordinates and independent gates remain unchanged. No stronger construction forces, source coordinate restoration or additional random retries.',
            'scope': 'Generation feasibility; native probabilities are not the independent CCTBX acceptance gate. Complete retained environment and final force-field refinement remain separate.',
            'app_ready': False, 'new_md_or_qm': False}
    save(output / 'plan.json', plan)
    rows = []
    started = time.monotonic()
    for case in cases:
        for protocol in protocols:
            if time.monotonic() - started > 900:
                save(output / 'overall-stop.json', {'reason': 'Frozen overall native wall budget reached', 'completed_batches': len(rows)})
                return
            folder = output / case['case'] / protocol['name']
            folder.mkdir(parents=True)
            command = [str(NATIVE), '-I', '-B', str(implementation / 'native_mc_adapter.py'),
                       '--plan', str(output / 'plan.json'), '--case', case['case'], '--mode', protocol['mode'],
                       '--boundary-protocol', protocol['boundary_protocol'], '--output', str(folder / 'native')]
            report = supervise(command, folder)
            rows.append({'case': case['case'], 'protocol': protocol['name'], **report})
            save(output / 'supervised-batches.json', rows)
            print(json.dumps(rows[-1]), flush=True)
    assert all(sha(Path(p)) == h for p, h in plan['source_sha256'].items())


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--output', type=Path, required=True)
    run(parser.parse_args().output.resolve())
