"""Bounded relocation-only regression of the qualified native RESP executable."""
import hashlib
import json
import math
import os
from pathlib import Path
import shutil
import signal
import subprocess
import time

ROOT = Path.home() / '.cache/dynamol-runtimes/resp-heap-local-v1'
EVIDENCE = ROOT / 'local-qualification'
RUN = EVIDENCE / 'native-parity-v3'


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def write(path, data):
    path.write_text(json.dumps(data, indent=2) + '\n')


def main():
    RUN.mkdir(exist_ok=False)
    plan = {'scope': 'Native RESP storage-only relocation regression',
            'source_binary_sha256': '1bde6a28988bfb9d08caee06e54d4f1b086808d154861bc9e0fd48cb3935ac3b',
            'candidate_sha256': '2d759392c83363f4294c01bc7c8602134e51b3a58f6b350044e1becce61048a3',
            'rss_limit_bytes': 3_000_000_000, 'per_process_seconds': 20,
            'campaign_seconds': 100, 'threads': 1,
            'native_workspace_bytes': 4 * 8000 * 8000 * 8,
            'native_file_arguments': 'Short relative names in a private local working directory; retained v2 demonstrates native absolute-path truncation and exit0 without fit results.',
            'resource_change_basis': 'Original source initializes four maxq-by-maxq double arrays, maxq=8000. The retained v1 1GB test budget was below its 2.048GB workspace alone. No binary, input or fit acceptance change.',
            'acceptance': 'All printed charges identical to retained outputs for both stages of the 113-atom reference and constrained synthetic 94-atom fixture; finite results, completed fitting statistics and all18frozen charges preserved.',
            'created_unix': time.time()}
    assert sha(ROOT / 'bin/resp') == plan['candidate_sha256']
    write(RUN / 'plan.json', plan)
    source = Path('/Users/ashujo/Documents/Science/DynaMol/build/advanced-chemistry/resp-runtime/heap-v2')
    sources = {}
    for name in ['resp.F', 'limits.h', 'resp_heap_workspace.f90']:
        original = source / name
        data = original.read_bytes()
        assert len(data) == original.stat().st_size and data
        (RUN / name).write_bytes(data)
        sources[name] = sha(RUN / name)
    assert b'call initialize_resp_heap(maxq)' in (RUN / 'resp.F').read_bytes()
    assert b'parameter (maxq   = 8000)' in (RUN / 'limits.h').read_bytes()
    write(RUN / 'source-pins.json', sources)
    cap_path = Path.home() / '.cache/dynamol-research/6oim-resp-recovery-v1/prepared/sources/cap/resp/constraints.json'
    assert sha(cap_path) == '7e19a5ae0c40cf3bec67fbce7b2cbf55b293722cf3212f16e0f53af0164ca1b7'
    constraints = json.loads(cap_path.read_bytes())
    (RUN / 'constraints.json').write_bytes(cap_path.read_bytes())
    started = time.monotonic()
    results = []
    try:
        for case, count in [('reference113', 113), ('synthetic94', 94)]:
            inputs = EVIDENCE / 'fixtures' / case
            output = RUN / case
            output.mkdir()
            for row in json.loads((inputs / 'manifest.json').read_bytes()):
                assert sha(inputs / row['name']) == row['sha256']
            for name in ['stage1.in', 'stage2.in']:
                (output / name).write_bytes((inputs / name).read_bytes())
            (output / 'esp.dat').write_bytes((inputs / ('input.esp' if count == 113 else 'esp.dat')).read_bytes())
            if count == 94:
                (output / 'canonical.qin').write_bytes((inputs / 'canonical.qin').read_bytes())
            for stage in (1, 2):
                prefix = output / f'stage{stage}'
                args = [str(ROOT / 'bin/resp'), '-O', '-i', f'stage{stage}.in',
                        '-o', f'stage{stage}.out', '-p', f'stage{stage}.punch',
                        '-t', f'stage{stage}.qout', '-s', f'stage{stage}.esout', '-e', 'esp.dat']
                if stage == 2:
                    args += ['-q', 'stage1.qout']
                elif count == 94:
                    args += ['-q', 'canonical.qin']
                env = {k: v for k, v in os.environ.items() if not k.startswith('DYLD_')}
                env.update(DYLD_PRINT_LIBRARIES='1', OMP_NUM_THREADS='1',
                           OPENBLAS_NUM_THREADS='1', VECLIB_MAXIMUM_THREADS='1')
                before = time.monotonic()
                peak = 0
                with prefix.with_suffix('.stdout').open('x') as stdout, prefix.with_suffix('.stderr').open('x') as stderr:
                    proc = subprocess.Popen(args, cwd=output, env=env, stdout=stdout, stderr=stderr, start_new_session=True)
                    write(prefix.with_suffix('.launch.json'), {'pid': proc.pid, 'args': args, 'started_unix': time.time()})
                    try:
                        with prefix.with_suffix('.memory.jsonl').open('x') as samples:
                            while proc.poll() is None:
                                observed = subprocess.run(['ps', '-p', str(proc.pid), '-o', 'rss='], capture_output=True, text=True, timeout=2).stdout.strip()
                                if observed:
                                    rss = int(observed) * 1024
                                    peak = max(peak, rss)
                                    samples.write(json.dumps({'elapsed_seconds': time.monotonic() - before, 'rss_bytes': rss}) + '\n')
                                    samples.flush()
                                    if rss > plan['rss_limit_bytes']:
                                        raise RuntimeError(f'Observed RSS {rss} exceeds limit')
                                if time.monotonic() - before > 20 or time.monotonic() - started > 100:
                                    raise TimeoutError('Native regression time bound exceeded')
                                if shutil.disk_usage(output).free < 8 * 1024**3:
                                    raise RuntimeError('Disk reserve crossed')
                                time.sleep(0.05)
                        assert proc.returncode == 0
                    except BaseException:
                        try:
                            os.killpg(proc.pid, signal.SIGKILL)
                        except ProcessLookupError:
                            pass
                        proc.wait(timeout=3)
                        raise
                assert 'Statistics of the fitting:' in prefix.with_suffix('.out').read_text()
                reference = inputs / (f'resp{stage}.chg' if count == 113 else f'stage{stage}.qout')
                observed = prefix.with_suffix('.qout').read_text().split()
                assert len(observed) == count and observed == reference.read_text().split()
                charges = [float(q.replace('D', 'E')) for q in observed]
                assert all(math.isfinite(q) for q in charges)
                loaded = prefix.with_suffix('.stderr').read_text()
                assert '/Documents/' not in loaded
                assert all(str(ROOT / 'lib' / n) in loaded for n in ['libgfortran.5.dylib', 'libquadmath.0.dylib', 'libgcc_s.1.1.dylib'])
                fixed = cap = None
                if count == 94:
                    assert len(constraints['frozen_atoms']) == 18
                    fixed = max(abs(charges[a['index']] - a['charge_e']) for a in constraints['frozen_atoms'])
                    cap = math.fsum(charges[i] for i in constraints['cap_indices'])
                    assert fixed == 0 and abs(cap) < 1e-12
                results.append({'case': case, 'stage': stage, 'pid': proc.pid, 'exit_code': proc.returncode,
                                'atoms': count, 'printed_charges_identical': True,
                                'candidate_qout_sha256': sha(prefix.with_suffix('.qout')),
                                'reference_qout_sha256': sha(reference), 'peak_sampled_rss_bytes': peak,
                                'wall_seconds': time.monotonic() - before,
                                'maximum_fixed_charge_error_e': fixed, 'cap_charge_e': cap,
                                'loaded_native_dependencies_local': True})
        assert sha(ROOT / 'bin/resp') == plan['candidate_sha256']
        result = {'passed': True, 'plan_sha256': sha(RUN / 'plan.json'), 'stages': results,
                  'wall_seconds': time.monotonic() - started, 'physical_acceptance': False,
                  'new_quantum_evaluations': 0}
        write(RUN / 'result.json', result)
        print(json.dumps(result, indent=2))
    except BaseException as error:
        write(RUN / 'failure.json', {'error': type(error).__name__, 'message': str(error),
                                    'completed_stages': results, 'physical_acceptance': False})
        raise


if __name__ == '__main__':
    main()
