"""Isolated geomeTRIC 1.1.1 / GFN2-xTB constrained geometry qualification.

This is a small-molecule numerical adapter, not an admitted full-adduct geometry
provider. It accepts at most 32 atoms and never reports physical acceptance. A
future full-parent lane requires a separately reviewed admission contract binding
the accepted upstream geometry, complete reacted atom graph, stereochemistry and
electronic state; changing this atom limit alone does not provide that contract.

The native entry point is run_native(), which supervises a separate process group
with a hard wall limit. The worker holds calculator.execution_guard throughout.
Native constraints operate in geomeTRIC's internal coordinate solver (enforce=0.1),
not as added energy terms. Its custom Engine forwards unmodified Hartree energies
and positive dE/dx gradients in Hartree/Bohr. Final constraints are independently
checked; a last iterate is never accepted as a converged optimization.
Native convergence energy thresholds are in Hartree, gradient thresholds in
Hartree/Bohr, and displacement thresholds (drms/dmax) in Angstrom.

Primary API sources: https://geometric.readthedocs.io/en/latest/engines.html and
https://geometric.readthedocs.io/en/latest/constraints.html, plus pinned 1.1.1
geometric.optimize.Optimizer / OPT_STATE and tests/test_customengine.py.
"""
from __future__ import annotations

import argparse
from contextlib import contextmanager
from copy import deepcopy
import hashlib
import json
import math
import os
from pathlib import Path
import signal
import subprocess
import sys
import threading
import time
import traceback

HERE = Path(__file__).resolve().parent
RUNTIME = Path('/Users/ashujo/.cache/dynamol-runtimes/covalent-geometry-v1')
SCOPE = 'small_closed_shell_numerical_qualification'
DEFAULTS = dict(convergence_energy=1e-6, convergence_grms=3e-4,
                convergence_gmax=4.5e-4, convergence_drms=0.0012,
                convergence_dmax=0.0018)
CONVERGENCE_UNITS = dict(convergence_energy='Hartree', convergence_grms='Hartree/Bohr',
                        convergence_gmax='Hartree/Bohr', convergence_drms='Angstrom',
                        convergence_dmax='Angstrom')
FROZEN_TOLERANCE_BOHR = 1e-5
ELEMENTS = {'H': 1, 'C': 6, 'N': 7, 'O': 8, 'F': 9, 'S': 16}
CALC_KEYS = {'schema_version', 'method', 'atom_ids', 'elements', 'charge', 'spin',
             'coords_bohr', 'accuracy', 'max_iter', 'temperature_hartree'}


def canonical(value):
    return json.dumps(value, sort_keys=True, separators=(',', ':'), allow_nan=False).encode()


def digest(value):
    return hashlib.sha256(canonical(value)).hexdigest()


def file_sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def write_json(path, value):
    path = Path(path)
    temporary = path.with_name(path.name + '.tmp')
    temporary.write_text(json.dumps(value, indent=2, sort_keys=True, allow_nan=False) + '\n')
    temporary.replace(path)


def finite_number(value):
    return type(value) in (int, float) and math.isfinite(value)


def dihedral_degrees(coordinates):
    """Same orientation/singularity convention as backend/qm_worker.py."""
    import numpy as np
    xyz = np.asarray(coordinates, dtype=float)
    if xyz.shape != (4, 3) or not np.isfinite(xyz).all():
        raise ValueError('Dihedral requires finite 4x3 coordinates')
    left, axis, right = xyz[0]-xyz[1], xyz[2]-xyz[1], xyz[3]-xyz[2]
    lengths = np.array([np.linalg.norm(left), np.linalg.norm(axis), np.linalg.norm(right)])
    if lengths.min() <= 1e-8:
        raise ValueError('Singular dihedral: coincident defining atoms')
    axis /= lengths[1]
    first, second = left-np.dot(left, axis)*axis, right-np.dot(right, axis)*axis
    if np.linalg.norm(first)/lengths[0] <= 1e-3 or np.linalg.norm(second)/lengths[2] <= 1e-3:
        raise ValueError('Singular dihedral: linear or nearly linear defining angle')
    return float(np.degrees(np.arctan2(np.dot(np.cross(axis, first), second), np.dot(first, second))))


def _coords(value, count):
    import numpy as np
    if (not isinstance(value, list) or len(value) != count or
        any(not isinstance(row, list) or len(row) != 3 or
            any(not finite_number(x) for x in row) for row in value)):
        raise ValueError('Finite explicit Nx3 coordinates required')
    return np.array(value, dtype=float)


def validate(request):
    """Strict admission; supplied graph is used, never inferred/replaced."""
    import numpy as np
    if not isinstance(request, dict) or set(request) != {
        'schema_version', 'scope', 'calculator', 'identity', 'optimization'}:
        raise ValueError('Unknown or missing adapter input fields')
    if type(request['schema_version']) is not int or request['schema_version'] != 1 or request['scope'] != SCOPE:
        raise ValueError('Only explicit small-molecule numerical qualification scope is admitted')
    calc = request['calculator']
    if not isinstance(calc, dict) or set(calc)-CALC_KEYS:
        raise ValueError('Unknown calculator setting')
    if type(calc.get('schema_version')) is not int or calc['schema_version'] != 1 or calc.get('method') != 'GFN2-xTB':
        raise ValueError('Explicit GFN2-xTB method and schema required')
    ids, elements = calc.get('atom_ids'), calc.get('elements')
    if not isinstance(ids, list) or not 2 <= len(ids) <= 32 or any(not isinstance(x, str) or not x for x in ids) or len(set(ids)) != len(ids):
        raise ValueError('Qualification requires 2-32 unique atom IDs; full-parent admission is not implemented')
    if not isinstance(elements, list) or len(elements) != len(ids) or any(not isinstance(x, str) or x not in ELEMENTS for x in elements):
        raise ValueError('Explicit H/C/N/O/F/S elements required')
    if type(calc.get('charge')) is not int or type(calc.get('spin')) is not int or calc['spin'] != 0:
        raise ValueError('Explicit integer charge and closed-shell spin 0 required')
    electrons = sum(ELEMENTS[x] for x in elements)-calc['charge']
    if electrons <= 0 or electrons % 2:
        raise ValueError('Closed-shell state requires positive even electron count')
    xyz = _coords(calc.get('coords_bohr'), len(ids))
    distances = np.linalg.norm(xyz[:, None]-xyz[None, :], axis=2)
    np.fill_diagonal(distances, np.inf)
    if distances.min() < 0.3:
        raise ValueError('Coincident or implausibly close atoms')
    accuracy, cycles, temperature = calc.get('accuracy', .001), calc.get('max_iter', 250), calc.get('temperature_hartree', .00095)
    if not finite_number(accuracy) or not 1e-5 <= accuracy <= 1:
        raise ValueError('Invalid SCC accuracy')
    if type(cycles) is not int or not 1 <= cycles <= 250:
        raise ValueError('Invalid SCC iteration limit')
    if not finite_number(temperature) or not 0 <= temperature <= .001:
        raise ValueError('Invalid electronic temperature')
    identity = request['identity']
    if not isinstance(identity, dict) or set(identity) != {'bonds', 'stereocenters'}:
        raise ValueError('Explicit bond graph and stereocenter declarations required')
    bonds = identity['bonds']
    if not isinstance(bonds, list) or not bonds:
        raise ValueError('A supplied nonempty bond graph is required')
    seen = set()
    neighbors = {i: set() for i in ids}
    for bond in bonds:
        if not isinstance(bond, dict) or set(bond) != {'atom_ids', 'order'}:
            raise ValueError('Each bond requires atom_ids and explicit order')
        pair = bond['atom_ids']
        if not isinstance(pair, list) or len(pair) != 2 or any(x not in ids for x in pair) or pair[0] == pair[1]:
            raise ValueError('Invalid bond atom IDs')
        key = tuple(sorted(pair))
        if key in seen or type(bond['order']) not in (int, float) or bond['order'] not in (1, 1.5, 2, 3):
            raise ValueError('Duplicate bond or unsupported explicit bond order')
        seen.add(key)
        neighbors[pair[0]].add(pair[1]); neighbors[pair[1]].add(pair[0])
    reached = {ids[0]}
    while True:
        expanded = reached.union(*(neighbors[i] for i in reached))
        if expanded == reached:
            break
        reached = expanded
    if reached != set(ids):
        raise ValueError('Only a single connected supplied molecular graph is admitted')
    stereo = identity['stereocenters']
    if not isinstance(stereo, list):
        raise ValueError('stereocenters must be an explicit list')
    centers = set()
    for row in stereo:
        if not isinstance(row, dict) or set(row) != {'center_atom_id', 'ordered_neighbor_ids', 'signed_volume_sign'}:
            raise ValueError('Invalid stereocenter declaration')
        center, selected = row['center_atom_id'], row['ordered_neighbor_ids']
        if center not in ids or center in centers or not isinstance(selected, list) or len(selected) != 4 or len(set(selected)) != 4 or set(selected) != neighbors[center]:
            raise ValueError('Stereocenter must declare four distinct bonded neighbors')
        if type(row['signed_volume_sign']) is not int or row['signed_volume_sign'] not in (-1, 1):
            raise ValueError('Stereocenter requires explicit signed-volume sign')
        centers.add(center)
    opt = request['optimization']
    allowed = {'max_steps', 'max_evaluations', 'max_wall_seconds', 'freeze_atom_ids',
               'dihedral_constraints', 'dihedral_tolerance_degrees'}
    if not isinstance(opt, dict) or set(opt)-allowed:
        raise ValueError('Unknown optimization setting; convergence thresholds are fixed')
    normalized = dict(max_steps=100, max_evaluations=300, max_wall_seconds=600,
                      freeze_atom_ids=[], dihedral_constraints=[], dihedral_tolerance_degrees=.1)
    normalized.update(opt)
    for key, maximum in [('max_steps', 200), ('max_evaluations', 600), ('max_wall_seconds', 600)]:
        if type(normalized[key]) is not int or not 1 <= normalized[key] <= maximum:
            raise ValueError('Invalid bounded ' + key)
    freezes = normalized['freeze_atom_ids']
    if not isinstance(freezes, list) or any(x not in ids for x in freezes) or len(set(freezes)) != len(freezes):
        raise ValueError('Freeze IDs must uniquely identify existing atoms')
    if len(freezes) == len(ids):
        raise ValueError('Cannot freeze every atom for optimization')
    tolerance = normalized['dihedral_tolerance_degrees']
    if not finite_number(tolerance) or not 0 < tolerance <= .5:
        raise ValueError('Periodic dihedral tolerance must be positive and at most 0.5 degrees')
    constraints = normalized['dihedral_constraints']
    if not isinstance(constraints, list) or len(constraints) > 8:
        raise ValueError('At most eight explicit dihedral constraints are admitted')
    seen = set()
    for row in constraints:
        if not isinstance(row, dict) or set(row) != {'atom_ids', 'target_degrees'}:
            raise ValueError('Each dihedral requires atom_ids and target_degrees')
        selected, target = row['atom_ids'], row['target_degrees']
        if not isinstance(selected, list) or len(selected) != 4 or any(x not in ids for x in selected) or len(set(selected)) != 4:
            raise ValueError('Dihedral requires four distinct existing atom IDs')
        if not finite_number(target) or abs(target) > 360:
            raise ValueError('Dihedral target must be finite and within -360 to 360 degrees')
        key = min(tuple(selected), tuple(selected[::-1]))
        if key in seen:
            raise ValueError('Duplicate or reversed dihedral constraint')
        if set(selected) <= set(freezes):
            raise ValueError('A fully frozen dihedral is redundant or conflicting')
        if any(b not in neighbors[a] for a, b in zip(selected[:-1], selected[1:])):
            raise ValueError('Dihedral definition must follow the supplied bonded graph')
        seen.add(key)
        dihedral_degrees(xyz[[ids.index(x) for x in selected]])
    result = deepcopy(request)
    result['optimization'] = normalized
    # Sign verification is a geometric check, not a CIP assignment/validation.
    _stereo_rows(xyz, ids, stereo)
    return result


def _stereo_rows(xyz, ids, definitions):
    import numpy as np
    rows = []
    for row in definitions:
        points = np.asarray(xyz)[[ids.index(x) for x in row['ordered_neighbor_ids']]]
        volume = float(np.linalg.det(points[:3]-points[3]))
        sign = 1 if volume > 0 else -1
        passed = abs(volume) > 1e-5 and sign == row['signed_volume_sign']
        rows.append({**row, 'signed_volume_bohr3': volume, 'satisfied': passed})
        if not passed:
            raise ValueError('Declared stereocenter sign is inverted or degenerate: ' + row['center_atom_id'])
    return rows


def constraint_text(request):
    ids = request['calculator']['atom_ids']
    opt = request['optimization']
    lines = []
    if opt['freeze_atom_ids']:
        lines += ['$freeze', 'xyz ' + ','.join(str(ids.index(x)+1) for x in opt['freeze_atom_ids'])]
    if opt['dihedral_constraints']:
        lines += ['$set']
        for row in opt['dihedral_constraints']:
            lines.append('dihedral ' + ' '.join(str(ids.index(x)+1) for x in row['atom_ids']) +
                         f" {math.remainder(row['target_degrees'], 360.0):.15g}")
    return '\n'.join(lines) + ('\n' if lines else '')


def verify_geometry(request, final):
    import numpy as np
    calc, opt = request['calculator'], request['optimization']
    ids = calc['atom_ids']
    xyz = _coords(final, len(ids))
    initial = np.asarray(calc['coords_bohr'], dtype=float)
    frozen = []
    for atom_id in opt['freeze_atom_ids']:
        displacement = float(np.linalg.norm(xyz[ids.index(atom_id)]-initial[ids.index(atom_id)]))
        frozen.append({'atom_id': atom_id, 'displacement_bohr': displacement,
                       'satisfied': displacement <= FROZEN_TOLERANCE_BOHR})
    dihedrals = []
    for row in opt['dihedral_constraints']:
        measured = dihedral_degrees(xyz[[ids.index(x) for x in row['atom_ids']]])
        error = math.remainder(measured-row['target_degrees'], 360.0)
        dihedrals.append({**row, 'measured_degrees': measured, 'signed_periodic_error_degrees': error,
                          'satisfied': abs(error) <= opt['dihedral_tolerance_degrees']})
    stereo = _stereo_rows(xyz, ids, request['identity']['stereocenters'])
    # A broad geometry consistency check catches bond loss; this is not a
    # chemical-state or force-field validity test and imposes no energy term.
    bonds = []
    for row in request['identity']['bonds']:
        i, j = [ids.index(x) for x in row['atom_ids']]
        ratio = float(np.linalg.norm(xyz[i]-xyz[j])/np.linalg.norm(initial[i]-initial[j]))
        bonds.append({**row, 'final_to_initial_distance_ratio': ratio, 'satisfied': .65 <= ratio <= 1.45})
    passed = all(row['satisfied'] for group in (frozen, dihedrals, stereo, bonds) for row in group)
    return {'passed': passed, 'frozen_tolerance_bohr': FROZEN_TOLERANCE_BOHR,
            'dihedral_tolerance_degrees': opt['dihedral_tolerance_degrees'],
            'frozen_atoms': frozen, 'dihedrals': dihedrals, 'stereocenters': stereo, 'bond_geometry': bonds}


class EvaluationSession:
    def __init__(self, request, output, evaluator):
        self.request, self.output, self.evaluator = request, Path(output), evaluator
        self.count = 0
        self.deadline = time.monotonic()+request['optimization']['max_wall_seconds']
        self.provider_binding = None
        (self.output/'evaluations').mkdir()

    def check_time(self):
        if time.monotonic() >= self.deadline:
            raise TimeoutError('Optimization wall limit exceeded')

    def evaluate(self, coords):
        import numpy as np
        self.check_time()
        if self.count >= self.request['optimization']['max_evaluations']:
            raise RuntimeError('Optimization evaluation budget exhausted')
        self.count += 1
        calc = deepcopy(self.request['calculator'])
        calc['coords_bohr'] = np.asarray(coords, dtype=float).reshape(-1, 3).tolist()
        _coords(calc['coords_bohr'], len(calc['atom_ids']))
        prefix = self.output/'evaluations'/f'{self.count:04d}'
        write_json(str(prefix)+'-request.json', calc)
        result = self.evaluator(deepcopy(calc))
        # Preserve malformed returns in their own file too, if JSON-serializable.
        try:
            write_json(str(prefix)+'-result.json', result)
        except (ValueError, TypeError):
            Path(str(prefix)+'-invalid-result.txt').write_text(repr(result))
        self.check_time()
        if not isinstance(result, dict) or result.get('accepted') is not True or result.get('physical_acceptance') is not False:
            raise ValueError('Calculator must report numerical acceptance without physical acceptance')
        for key in ('atom_ids', 'elements', 'charge', 'spin', 'method', 'coords_bohr'):
            if result.get(key) != calc[key]:
                raise ValueError('Calculator changed requested identity/state/geometry: ' + key)
        if type(result['charge']) is not int or type(result['spin']) is not int:
            raise ValueError('Calculator electronic state must retain explicit integer types')
        if result.get('units') != {'coordinates': 'Bohr', 'energy': 'Hartree', 'gradient': 'Hartree/Bohr; positive dE/dx'}:
            raise ValueError('Calculator energy/gradient units or gradient sign convention are not explicit')
        energy = result.get('energy_hartree')
        if not finite_number(energy):
            raise ValueError('Nonfinite calculator energy')
        gradient = _coords(result.get('gradient_hartree_per_bohr'), len(calc['atom_ids']))
        binding = {key: result.get(key) for key in ('calculator_sha256', 'runtime', 'native_version', 'settings')}
        if any(value is None for value in binding.values()):
            raise ValueError('Calculator provenance is incomplete')
        if self.provider_binding is None:
            self.provider_binding = deepcopy(binding)
        elif binding != self.provider_binding:
            raise ValueError('Calculator runtime/settings binding changed during optimization')
        return {'energy': energy, 'gradient': gradient.ravel()}, result


def _native_optimizer(request, output, session):
    import numpy as np
    import geometric
    from geometric.engine import Engine
    from geometric.internal import DelocalizedInternalCoordinates
    from geometric.molecule import Molecule
    from geometric.nifty import bohr2ang
    from geometric.optimize import Optimizer, OPT_STATE
    from geometric.params import OptParams
    from geometric.prepare import parse_constraints
    if geometric.__version__ != '1.1.1':
        raise ValueError('Only inspected geomeTRIC 1.1.1 API is admitted')
    ids = request['calculator']['atom_ids']
    xyz = np.array(request['calculator']['coords_bohr'], dtype=float)
    molecule = Molecule()
    molecule.elem = request['calculator']['elements']
    molecule.xyzs = [xyz*bohr2ang]
    molecule.comms = ['DynaMol small-molecule numerical qualification; physical_acceptance=false']
    molecule.bonds = sorted(tuple(sorted(ids.index(x) for x in row['atom_ids'])) for row in request['identity']['bonds'])
    molecule.top_settings['read_bonds'] = True
    molecule.build_topology(force_bonds=False)

    class NativeEngine(Engine):
        def calc_new(self, coords, dirname):
            return session.evaluate(coords)[0]

    engine = NativeEngine(molecule)
    text = constraint_text(request)
    (output/'constraints.txt').write_text(text)
    cons, values = parse_constraints(molecule, text) if text else (None, None)
    coordinates = DelocalizedInternalCoordinates(molecule, build=True, connect=True,
        addcart=False, constraints=cons, cvals=values[0] if values is not None else None)
    parameters = OptParams(maxiter=request['optimization']['max_steps'], enforce=.1,
        check=1, **DEFAULTS)
    optimizer = Optimizer(xyz.ravel(), molecule, coordinates, engine,
                          str(output/'native-work'), parameters, print_info=False)
    try:
        optimizer.optimizeGeometry()
    finally:
        write_json(output/'optimizer-state.json', {'state': optimizer.state,
            'native_converged': optimizer.state == OPT_STATE.CONVERGED,
            'iterations': optimizer.Iteration, 'coords_bohr': optimizer.X.reshape(-1, 3).tolist()})
    session.check_time()
    return {'native_converged': optimizer.state == OPT_STATE.CONVERGED,
            'native_state': optimizer.state, 'iterations': optimizer.Iteration,
            'coords_bohr': optimizer.X.reshape(-1, 3).tolist(),
            'optimizer': 'geomeTRIC', 'optimizer_version': geometric.__version__,
            'constraint_algorithm': 'DLC with enforce=0.1; no restraint energy'}


def _execute(request, output, evaluator, runner, evidence_kind):
    """Worker core; injected runners/evaluators are only synthetic contract tests."""
    output = Path(output)
    started = time.monotonic()
    session = None
    result = {'accepted': False, 'physical_acceptance': False, 'simulation_ready': False,
              'scope': SCOPE, 'evidence_kind': evidence_kind, 'adapter_sha256': file_sha(__file__)}
    try:
        normalized = validate(request)
        write_json(output/'request.json', request)
        write_json(output/'normalized-request.json', normalized)
        result.update(input_sha256=digest(request), identity_sha256=digest({
            key: normalized['calculator'][key] for key in ('atom_ids', 'elements', 'charge', 'spin')} | normalized['identity']),
            initial_geometry_sha256=digest(normalized['calculator']['coords_bohr']))
        session = EvaluationSession(normalized, output, evaluator)
        optimized = runner(normalized, output, session)
        write_json(output/'optimizer-result.json', optimized)
        if optimized.get('native_converged') is not True:
            raise ValueError('Native optimizer did not converge; last iterate rejected')
        if type(optimized.get('iterations')) is not int or not 0 <= optimized['iterations'] <= normalized['optimization']['max_steps']:
            raise ValueError('Native optimizer exceeded step budget')
        checks = verify_geometry(normalized, optimized.get('coords_bohr'))
        write_json(output/'final-geometry-checks.json', checks)
        if not checks['passed']:
            raise ValueError('Independent final constraint/bond geometry checks failed')
        # Deliberate uncached single point at the final accepted coordinates.
        _, final = session.evaluate(optimized['coords_bohr'])
        if file_sha(__file__) != result['adapter_sha256']:
            raise ValueError('Adapter source changed during execution')
        result.update(accepted=True, status='numerically_complete_research_geometry_only',
            native_optimization=optimized, geometry_checks=checks,
            atom_ids=final['atom_ids'], elements=final['elements'], charge=final['charge'], spin=final['spin'],
            coords_bohr=final['coords_bohr'], energy_hartree=final['energy_hartree'],
            gradient_hartree_per_bohr=final['gradient_hartree_per_bohr'],
            provider_binding=session.provider_binding, unbiased_electronic_energy=True,
            added_restraint_energy_hartree=0.0, convergence_thresholds=DEFAULTS,
            convergence_threshold_units=CONVERGENCE_UNITS,
            limitations='Small-molecule numerical qualification only. No full-parent admission, charge replacement, force-field validation, or app integration.')
    except Exception as error:
        result.update(status='rejected', error_type=type(error).__name__, error=str(error))
        (output/'failure-traceback.txt').write_text(traceback.format_exc())
    result.update(elapsed_seconds=time.monotonic()-started,
                  evaluations=session.count if session else 0)
    write_json(output/'result.json', result)
    return result


def run_synthetic_fixture(request, output, evaluator, runner):
    """No native launch. All injected tests carry an explicit synthetic label."""
    output = Path(output)
    output.mkdir(parents=True, exist_ok=False)
    return _execute(request, output, evaluator, runner, 'synthetic_engine_contract_fixture')


class Cancelled(Exception):
    """A supervisor cancellation, recorded only after its group is cleaned up."""


@contextmanager
def _cancellation_scope():
    """Latch cancellation without raising asynchronously at any boundary.

    Ordinary code checks the flag while polling. Signals therefore cannot abort
    Popen assignment, a finally block, or entry into the cleanup context itself.
    The child does not inherit a blocked signal mask.
    """
    if threading.current_thread() is not threading.main_thread():
        raise RuntimeError('The native process supervisor must run on the main thread')
    signals = {signal.SIGINT, signal.SIGTERM}
    state = {'requested': None}
    def cancel(signum, frame):
        if state['requested'] is None:
            state['requested'] = signum
    old_mask = signal.pthread_sigmask(signal.SIG_BLOCK, signals)
    handlers = {sig: signal.signal(sig, cancel) for sig in signals}
    try:
        signal.pthread_sigmask(signal.SIG_SETMASK, old_mask)
        yield state
    finally:
        signal.pthread_sigmask(signal.SIG_BLOCK, signals)
        for sig, handler in handlers.items():
            signal.signal(sig, handler)
        signal.pthread_sigmask(signal.SIG_SETMASK, old_mask)


@contextmanager
def _uninterrupted_cleanup():
    """Repeated cancellation cannot interrupt termination/reaping of our group."""
    signals = {signal.SIGINT, signal.SIGTERM}
    old_mask = signal.pthread_sigmask(signal.SIG_BLOCK, signals)
    handlers = {sig: signal.signal(sig, signal.SIG_IGN) for sig in signals}
    try:
        signal.pthread_sigmask(signal.SIG_SETMASK, old_mask)
        yield
    finally:
        signal.pthread_sigmask(signal.SIG_BLOCK, signals)
        for sig, handler in handlers.items():
            signal.signal(sig, handler)
        signal.pthread_sigmask(signal.SIG_SETMASK, old_mask)


def _group_is_alive(group_id):
    # A leader may already be reaped while its descendants remain. macOS can
    # return EPERM for killpg(..., 0) during teardown; inspect live group members.
    rows = subprocess.check_output(['ps', '-axo', 'pgid=,stat='], text=True,
                                   timeout=3).splitlines()
    return any(len(parts := row.split()) == 2 and int(parts[0]) == group_id and
               not parts[1].startswith('Z') for row in rows)


def _terminate_owned_group(child, grace_seconds=1.0):
    """Terminate descendants even if their leader exited; always reap the leader.

    macOS does not expose Linux subreapers: grandchildren are reaped by their
    parent or launchd. We verify there are no live executable group members.
    """
    try:
        try:
            os.killpg(child.pid, signal.SIGTERM)
        except ProcessLookupError:
            return
        except PermissionError:
            if not _group_is_alive(child.pid):
                return
            raise
        deadline = time.monotonic()+grace_seconds
        while time.monotonic() < deadline:
            child.poll()
            if not _group_is_alive(child.pid):
                return
            time.sleep(.025)
        try:
            os.killpg(child.pid, signal.SIGKILL)
        except ProcessLookupError:
            pass
        except PermissionError:
            if _group_is_alive(child.pid):
                raise
        deadline = time.monotonic()+3
        while _group_is_alive(child.pid):
            if time.monotonic() >= deadline:
                raise RuntimeError('Owned process group still has live members after SIGKILL')
            child.poll()
            time.sleep(.025)
    finally:
        child.wait(timeout=3)


def _run_owned_process(command, log_path, timeout, on_start, environment):
    """Own every post-spawn path, including artifact errors and cancellation."""
    child = None
    timed_out = False
    deadline = time.monotonic()+timeout
    with _cancellation_scope() as cancellation:
        with Path(log_path).open('w') as log:
            try:
                child = subprocess.Popen(command, stdout=log, stderr=subprocess.STDOUT,
                    start_new_session=True, env=environment)
                if cancellation['requested'] is not None:
                    raise Cancelled(f"Supervisor cancellation signal {cancellation['requested']}")
                on_start(child.pid)
                while True:
                    if cancellation['requested'] is not None:
                        raise Cancelled(f"Supervisor cancellation signal {cancellation['requested']}")
                    remaining = deadline-time.monotonic()
                    if remaining <= 0:
                        timed_out = True
                        raise TimeoutError('Native process group exceeded hard wall limit')
                    try:
                        code = child.wait(timeout=min(.1, remaining))
                        break
                    except subprocess.TimeoutExpired:
                        pass
            finally:
                if child is not None:
                    with _uninterrupted_cleanup():
                        _terminate_owned_group(child, grace_seconds=0.0 if timed_out else 1.0)
            # A first cancellation can arrive after normal completion but before
            # cleanup installs ignored handlers. Honor that request after cleanup.
            if cancellation['requested'] is not None:
                raise Cancelled(f"Supervisor cancellation signal {cancellation['requested']}")
            return code


def run_native(request, output, root):
    """Hard-bounded, process-isolated native entry; never overwrite an old run."""
    output = Path(output).resolve()
    output.mkdir(parents=True, exist_ok=False)
    try:
        normalized = validate(request)
        write_json(output/'launch-input.json', normalized)
        binding = {'input_sha256': file_sha(output/'launch-input.json'),
                   'adapter_sha256': file_sha(__file__), 'calculator_sha256': file_sha(HERE/'calculator.py'),
                   'runtime_manifest_sha256': file_sha(HERE/'runtime-manifest.json')}
        write_json(output/'launch-binding.json', binding)
        command = [str(RUNTIME/'bin/python'), str(Path(__file__).resolve()), '--worker',
                   '--root', str(Path(root).resolve()), '--output', str(output)]
        environment = dict(os.environ)
        for key in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS', 'VECLIB_MAXIMUM_THREADS', 'NUMEXPR_NUM_THREADS'):
            environment[key] = '1'
        def record_process(pid):
            write_json(output/'process.json', {'pid': pid, 'command': command,
                'max_wall_seconds': normalized['optimization']['max_wall_seconds']})
        code = _run_owned_process(command, output/'native.log',
            normalized['optimization']['max_wall_seconds'], record_process, environment)
        if code not in (0, 1) or not (output/'result.json').exists():
            raise RuntimeError(f'Native worker exited {code} without a complete result')
        result = json.loads((output/'result.json').read_text())
        if (not isinstance(result, dict) or type(result.get('accepted')) is not bool or
            result.get('physical_acceptance') is not False or
            result['accepted'] and code != 0):
            raise ValueError('Worker result and exit state are inconsistent')
        for key, path in [('input_sha256', output/'launch-input.json'), ('adapter_sha256', __file__),
                          ('calculator_sha256', HERE/'calculator.py'), ('runtime_manifest_sha256', HERE/'runtime-manifest.json')]:
            if file_sha(path) != binding[key]:
                raise ValueError('Launch input/source/runtime binding changed: '+key)
        return result
    except (Exception, KeyboardInterrupt) as error:
        result = {'accepted': False, 'physical_acceptance': False, 'simulation_ready': False,
                  'status': 'rejected', 'evidence_kind': 'native_small_molecule_numerical_qualification',
                  'error_type': type(error).__name__, 'error': str(error)}
        if (output/'result.json').exists():
            (output/'result.json').rename(output/'worker-result-before-supervisor-rejection.json')
        write_json(output/'result.json', result)
        (output/'supervisor-failure-traceback.txt').write_text(traceback.format_exc())
        return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--input', type=Path)
    parser.add_argument('--output', required=True, type=Path)
    parser.add_argument('--root', required=True, type=Path)
    parser.add_argument('--worker', action='store_true', help=argparse.SUPPRESS)
    args = parser.parse_args()
    if args.worker:
        # The supervisor owns the fresh directory and starts this process in a
        # new process group. Direct worker invocation is not a public API.
        import calculator
        binding = json.loads((args.output/'launch-binding.json').read_text())
        for key, path in [('input_sha256', args.output/'launch-input.json'), ('adapter_sha256', __file__),
                          ('calculator_sha256', HERE/'calculator.py'), ('runtime_manifest_sha256', HERE/'runtime-manifest.json')]:
            if file_sha(path) != binding[key]:
                raise ValueError('Worker launch binding mismatch: '+key)
        request = json.loads((args.output/'launch-input.json').read_text())
        with calculator.execution_guard(args.root):
            result = _execute(request, args.output, calculator.evaluate, _native_optimizer,
                              'native_small_molecule_numerical_qualification')
    else:
        if args.input is None:
            parser.error('--input is required')
        result = run_native(json.loads(args.input.read_text()), args.output, args.root)
    print(json.dumps({key: result.get(key) for key in ('accepted', 'status', 'error', 'evaluations')}))
    return 0 if result['accepted'] else 1


if __name__ == '__main__':
    sys.exit(main())
