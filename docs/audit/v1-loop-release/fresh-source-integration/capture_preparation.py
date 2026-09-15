#!/usr/bin/env python3
"""One research-only fresh preparation, intercepted before loop refinement.

No prepared archive is an input. Snapshot hashes establish identity/serialization,
not chemical or geometric acceptance. Native AM1-BCC/SQM may run as part of the
actual existing ligand preparation. The app's default and installed app are untouched.
"""
from __future__ import annotations
import argparse
import hashlib
import json
import os
from pathlib import Path, PurePosixPath
import posixpath
import shutil
import signal
import subprocess
import sys
import time
import traceback
from xml.etree import ElementTree

REPO = Path('/Users/ashujo/Documents/Science/DynaMol')
ROOT = Path('/Users/ashujo/.cache/dynamol-research/v1-loop-release/fresh-preparation-snapshot-v1-attempt-02')
PYTHON = Path('/Users/ashujo/.cache/dynamol-runtimes/covalent-v1-focused/bin/python')
AMBER = Path('/Users/ashujo/Library/Application Support/DynaMol/Runtimes/ambertools-b96b42253d7d46c5')
DATASET = '3a661cfc670048c5'
JOB = 'fresh-1ua2-capture'
LIMIT = 4 * 1024**3
THREADS = 2
DEADLINE = 600
SOURCE_NAMES = ['preparation_worker', 'loop_modeling', 'promod_loop_worker', 'loop_snapshot',
                'loop_process', 'loop_search', 'ligands', 'ligand_repair', 'modified_residues',
                'prepared_system', 'forcefield_identity', 'complex_topology', 'added_atom_stereo',
                'sequence_evidence', 'preparation', 'worker', 'storage', 'residue_identity', 'config']


def digest(path):
    with Path(path).open('rb') as stream:
        return hashlib.file_digest(stream, 'sha256').hexdigest()


def write(path, value):
    Path(path).parent.mkdir(parents=True, exist_ok=True)
    temporary = Path(str(path) + '.tmp')
    temporary.write_text(json.dumps(value, indent=2, allow_nan=False) + '\n')
    temporary.replace(path)


def verify(records, path_key='path'):
    for item in records:
        if digest(item[path_key]) != item['sha256']:
            raise ValueError('Frozen bytes changed: ' + item[path_key])


def setup_environment():
    os.environ.update(DYNAMOL_DATA_DIR=str(ROOT / 'workspace'), DYNAMOL_AMBERTOOLS=str(AMBER),
                      DYNAMOL_PROMOD3='/Users/ashujo/Library/Application Support/DynaMol/Runtimes/promod3-dfd2a559551cd2f2',
                      DYNAMOL_CPU_THREADS='2', OPENMM_CPU_THREADS='2', OPENBLAS_NUM_THREADS='2',
                      OMP_NUM_THREADS='2', MKL_NUM_THREADS='2', VECLIB_MAXIMUM_THREADS='2',
                      NUMEXPR_NUM_THREADS='2', PYTHONDONTWRITEBYTECODE='1')
    sys.path.insert(0, str(REPO))
    sys.path.insert(0, str(ROOT / 'overlay'))


def parameter_map(job, names):
    """Resolve explicit loader inputs once; enumerate every recursive XML include."""
    from openmm import app
    base = Path(app.__file__).parent / 'data'
    result = {}
    def add(name, path, ancestry=()):
        if name != PurePosixPath(name).as_posix() or name.startswith('/') or '..' in name.split('/'):
            raise ValueError('Unsafe logical parameter name')
        if name in ancestry:
            raise ValueError('Recursive XML include cycle')
        path = Path(path)
        if path.is_symlink() or not path.is_file():
            raise ValueError('Missing regular parameter file ' + str(path))
        if name in result:
            if result[name] != path:
                raise ValueError('Conflicting parameter source')
            return
        result[name] = path
        tree = ElementTree.fromstring(path.read_bytes())
        for node in tree.findall('Include'):
            value = node.get('file')
            if not value or Path(value).is_absolute() or '\\' in value:
                raise ValueError('Unsafe XML include')
            logical = posixpath.normpath(posixpath.join(posixpath.dirname(name), value))
            add(logical, path.parent / value, (*ancestry, name))
    for name in names:
        local = job / name
        add(name, local if local.is_file() else base / name)
    manifest = job / 'residue-parameters/manifest.json'
    if manifest.is_file():
        record = json.loads(manifest.read_text())
        for filename in ['manifest.json', *record['files']]:
            if filename != Path(filename).name:
                raise ValueError('Unsafe modified-residue sidecar name')
            path = manifest.parent / filename
            if path.is_symlink() or not path.is_file():
                raise ValueError('Missing modified-residue sidecar')
            if filename != 'manifest.json' and digest(path) != record['files'][filename]['sha256']:
                raise ValueError('Modified-residue sidecar changed')
            result['residue-parameters/' + filename] = path
    return result


def source_mapping(source, prepared, ligands):
    from backend.preparation_worker import atom_key
    from backend.storage import WATERS
    source_atoms, atoms = list(source.atoms()), list(prepared.atoms())
    targets = {atom_key(atom): atom.index for atom in atoms}
    if len(targets) != len(atoms):
        raise ValueError('Ambiguous prepared atom identities')
    renamed = {}
    for ligand in ligands:
        for item in ligand.atom_map:
            index = item['original_index']
            if index is None:
                continue  # explicitly modeled ligand atom; has no deposited source row
            if type(index) is not int or not 0 <= index < len(source_atoms):
                raise ValueError('Ligand original index is invalid')
            raw = source_atoms[index]
            if (raw.name != item['original_name'] or atom_key(raw)[:4] != tuple(ligand.original_residue_key)
                    or index in renamed):
                raise ValueError('Ligand atom map does not bind to the raw source')
            renamed[index] = (*ligand.original_residue_key, item['prepared_name'])
    mapping, events = [], []
    for atom in source_atoms:
        identity = atom_key(atom)
        if atom.element.atomic_number == 1:
            index, reason = None, 'Existing hydrogen explicitly removed before fixed-state reassignment'
        else:
            target = renamed.get(atom.index, identity)
            index = targets.get(target)
            if index is None:
                if atom.residue.name.upper() not in WATERS:
                    raise ValueError('Unexpected missing raw heavy atom: ' + repr(identity))
                reason = 'Water explicitly removed by remove_waters=true'
            else:
                if atoms[index].element != atom.element:
                    raise ValueError('Mapped element changed')
                reason = 'Explicit prepared_ligands.atom_map rename' if target != identity else 'Exact heavy atom identity'
        mapping.append(index)
        events.append({'source_index': atom.index, 'source_identity': identity, 'prepared_index': index,
                       'prepared_identity': atom_key(atoms[index]) if index is not None else None, 'reason': reason})
    if len(set(i for i in mapping if i is not None)) != sum(i is not None for i in mapping):
        raise ValueError('Source map is not injective')
    return mapping, events


def backbone_references(topology, mapping):
    from backend.preparation_worker import atom_key
    from backend.residue_identity import STANDARD_PROTEINS
    from backend.modified_residues import PARENT_RESIDUES
    residues = [r for r in topology.residues() if PARENT_RESIDUES.get(r.name, r.name) in STANDARD_PROTEINS]
    named = {r.index: {a.name: a for a in r.atoms()} for r in residues}
    previous, following = {}, {}
    for a, b in topology.bonds():
        if a.residue.index == b.residue.index or {a.name, b.name} != {'C', 'N'}:
            continue
        c, n = (a, b) if a.name == 'C' else (b, a)
        if c.residue.index in named and n.residue.index in named:
            if c.residue.index in following or n.residue.index in previous:
                raise ValueError('Ambiguous peptide adjacency')
            previous[n.residue.index] = c.residue.index
            following[c.residue.index] = n.residue.index
    observed = {i for i in mapping if i is not None}
    result = []
    for r in residues:
        own, prev, nxt = named[r.index], named.get(previous.get(r.index), {}), named.get(following.get(r.index), {})
        definitions = {'phi': [prev.get('C'), own.get('N'), own.get('CA'), own.get('C')],
                       'psi': [own.get('N'), own.get('CA'), own.get('C'), nxt.get('N')],
                       'preceding_omega': [prev.get('CA'), prev.get('C'), own.get('N'), own.get('CA')]}
        for label, atoms in definitions.items():
            if any(a is None for a in atoms):
                continue
            new = [atom_key(a) for a in atoms if a.index not in observed]
            if new:
                result.append({'residue_key': atom_key(next(r.atoms()))[:4], 'torsion': label,
                               'support': [atom_key(a) for a in atoms], 'new_support': new})
    return result


class CaptureComplete(Exception):
    pass


def capture(worker, modeller, forcefield, modeled_keys, context):
    import numpy as np
    import openmm as mm
    import pdbfixer
    from openmm import app, unit
    from backend.loop_snapshot import create_snapshot, load_snapshot
    from backend.preparation_worker import atom_key
    from backend.modified_residues import PARENT_RESIDUES, SUPPORTED_MODIFIED
    from backend.residue_identity import STANDARD_PROTEINS
    folder = ROOT / 'captured'
    folder.mkdir()
    source = app.PDBFile(str(worker.folder / 'input.pdb'))
    if source.topology.getNumAtoms() != 2322:
        raise ValueError('The raw input atom count changed')
    # Canonicalization in the current worker must be visible; no guessed equivalence.
    if [atom_key(a) for a in source.topology.atoms()] != [atom_key(a) for a in context['source_topology'].atoms()]:
        raise ValueError('Worker normalized a raw source identity without an explicit capture mapping')
    mapping, events = source_mapping(source.topology, modeller.topology, context['prepared_ligands'])
    references = backbone_references(modeller.topology, mapping)
    raw_xyz = np.asarray(source.positions.value_in_unit(unit.nanometer), dtype=float)
    xyz = np.asarray(modeller.positions.value_in_unit(unit.nanometer), dtype=float)
    movement = []
    for item in events:
        if item['prepared_index'] is not None:
            movement.append({**item, 'displacement_angstrom': float(np.linalg.norm(
                raw_xyz[item['source_index']] - xyz[item['prepared_index']]) * 10)})
    distances = np.asarray([r['displacement_angstrom'] for r in movement])
    movement_report = {'reference': 'Original raw input coordinates; common coordinates are not all experimentally observed',
                       'mapped_heavy_atoms': len(movement), 'rms_angstrom': float(np.sqrt(np.mean(distances**2))),
                       'maximum_angstrom': float(distances.max()), 'per_atom': movement}
    write(folder / 'source-movement.json', movement_report)
    write(folder / 'source-mapping.json', events)
    with (folder / 'prepared-topology.cif').open('x') as stream:
        app.PDBxFile.writeFile(modeller.topology, modeller.positions, stream, keepIds=True)
    # Complete physical model, without constraints or loop restraints. No Context or MD.
    system = forcefield.createSystem(modeller.topology, nonbondedMethod=app.NoCutoff,
                                    constraints=None, rigidWater=False)
    (folder / 'native-system.xml').write_text(mm.XmlSerializer.serialize(system))
    atoms = list(modeller.topology.atoms())
    parents = {a.index: [] for a in atoms if a.element.atomic_number == 1}
    for a, b in modeller.topology.bonds():
        if a.index in parents: parents[a.index].append(b.index)
        if b.index in parents: parents[b.index].append(a.index)
    hydrogens = [{'identity': atom_key(atoms[i]), 'index': i,
                  'parents': [{'index': j, 'identity': atom_key(atoms[j])} for j in p]}
                 for i, p in sorted(parents.items())]
    primary = json.loads((worker.folder / 'loop-model-input.json').read_text())
    requested = [(r['chain'], r['resid'], r['insertion_code'], r['residue'])
                 for chain in primary['chains'] for r in chain['residues'] if r['observed'] is False]
    if set(requested) != set(modeled_keys):
        raise ValueError('Common preparation only partially covers the primary sequence request')
    state = {'status': 'captured_before_refinement', 'capture_only': True, 'simulation_ready': False,
             'physical_model_validated': False, 'settings': worker.settings,
             **{k: context[k] for k in ['solvent', 'forcefield_files', 'parameter_state', 'protonation_states',
                                        'selected_variants', 'sidechain_adjustment', 'loop_modeler']},
             'protected_residues': sorted(context['protected_residues']),
             'newly_defined_backbone_residue_keys': sorted({tuple(r['residue_key']) for r in references}),
             'newly_defined_backbone_references': references,
             'hydrogen_inventory': hydrogens, 'source_mapping_file': 'source-mapping.json',
             'mapping_choices': 'Raw H removed; heavy identities exact except explicit prepared_ligands.atom_map renames. Removed waters recorded; no other heavy atom removal accepted.',
             'requested_modeled_keys': sorted(requested), 'modeled_keys': sorted(modeled_keys),
             'native_system': {'nonbonded_method': 'NoCutoff', 'constraints': None, 'rigid_water': False,
                               'particles': system.getNumParticles(), 'forces': [f.__class__.__name__ for f in system.getForces()],
                               'context_created_for_capture': False}}
    write(folder / 'parameter-state.json', state)
    raw = ROOT / 'workspace/datasets' / DATASET
    files = {'original-input.pdb': worker.folder / 'input.pdb',
             'prepared-topology.cif': folder / 'prepared-topology.cif',
             'native-system.xml': folder / 'native-system.xml', 'parameter-state.json': folder / 'parameter-state.json',
             'source-mapping.json': folder / 'source-mapping.json', 'source-movement.json': folder / 'source-movement.json',
             'loop-model-context.pdb': worker.folder / 'loop-model-context.pdb',
             'loop-model-input.json': worker.folder / 'loop-model-input.json',
             'sequence-source.cif': raw / 'sequence-source.cif',
             'original-source.cif': raw / 'original-source/aa3e69b79109.cif',
             'raw-metadata.json': raw / 'metadata.json', 'raw-provenance.json': raw / 'provenance.json',
             'raw-atom-map.json': raw / 'atom-map.json', 'capture-source.py': Path(__file__),
             'capture-plan.json': ROOT / 'plan.json'}
    template_root = Path(pdbfixer.__file__).parent / 'templates'
    for name in sorted({PARENT_RESIDUES.get(r.name, r.name) for r in modeller.topology.residues()
                        if r.name not in SUPPORTED_MODIFIED and PARENT_RESIDUES.get(r.name, r.name) in STANDARD_PROTEINS and r.name != 'GLY'}):
        files['templates/' + name + '.pdb'] = template_root / (name + '.pdb')
    for prefix in ['residue-parameters', 'ligands']:
        for p in sorted((worker.folder / prefix).rglob('*')):
            if p.is_file(): files['preparation-evidence/' + str(p.relative_to(worker.folder))] = p
    for p in sorted((ROOT / 'workspace/chemistry/ccd').glob('*.cif')):
        files['ccd/' + p.name] = p
    for name in ['loop-model-provenance.json', 'loop-model-output.json', 'loop-construction.json', 'added-atom-stereochemistry.json']:
        if (worker.folder / name).is_file(): files[name] = worker.folder / name
    snapshot = create_snapshot(ROOT / 'snapshot', source_topology=source.topology, source_positions=source.positions,
                               topology=modeller.topology, positions=modeller.positions,
                               modeled_keys=sorted(modeled_keys), requested_modeled_keys=sorted(requested),
                               source_to_prepared=mapping, preparation_provenance=state,
                               parameter_files=parameter_map(worker.folder, context['forcefield_files']), source_files=files)
    load_snapshot(snapshot.manifest_path.parent, expected_sha256=snapshot.sha256)
    write(ROOT / 'capture-result.json', {'status': 'captured_before_refinement', 'snapshot_sha256': snapshot.sha256,
          'snapshot_manifest': str(snapshot.manifest_path), 'source_atoms': len(snapshot.source.atoms),
          'prepared_atoms': len(snapshot.prepared.atoms), 'modeled_residues': len(snapshot.modeled_keys),
          'source_movement_rms_angstrom': movement_report['rms_angstrom'],
          'source_movement_maximum_angstrom': movement_report['maximum_angstrom'],
          'newly_defined_backbone_residue_keys': state['newly_defined_backbone_residue_keys'],
          'source_files': len(snapshot.source_files), 'parameter_files': len(snapshot.parameter_files),
          'simulation_ready': False, 'complete_quality_checks_performed': False, 'app_published': False})
    raise CaptureComplete('Immutable common preparation captured before refinement/publication')


def child():
    setup_environment()
    import psutil
    from backend.preparation_worker import PreparationWorker
    from backend.preparation import PREPARATION_DEFAULTS
    from backend.worker import Cancelled
    class CaptureWorker(PreparationWorker):
        def _refine_loops(self, modeller, forcefield, modeled_keys, *, context):
            return capture(self, modeller, forcefield, modeled_keys, context)
        def _before_publish_preparation(self, modeller, preparation):
            raise RuntimeError('Capture-only publication boundary must be unreachable')
    def cancel(signum, frame):
        raise Cancelled('Capture worker cancellation signal ' + str(signum))
    for sig in [signal.SIGTERM, signal.SIGINT]: signal.signal(sig, cancel)
    job = ROOT / 'workspace/jobs' / JOB
    job.mkdir(parents=True)
    shutil.copyfile(ROOT / 'workspace/datasets' / DATASET / 'topology.pdb', job / 'input.pdb')
    settings = {**PREPARATION_DEFAULTS, 'dataset_id': DATASET, 'name': 'Fresh 1UA2 common-preparation capture',
                'operation': 'prepare', 'build_missing_residues': True, 'ligand_actions': {}}
    write(job / 'status.json', {'id': JOB, 'config': settings, 'worker_pid': os.getpid(), 'status': 'starting',
                                'logs': [], 'total_steps': 7, 'completed_steps': 0, 'progress': 0})
    write(job / 'config.json', settings)
    original_popen = subprocess.Popen
    def recorded_popen(*args, **kwargs):
        process = original_popen(*args, **kwargs)
        record = {'pid': process.pid, 'parent_pid': os.getpid(), 'command': args[0] if args else kwargs.get('args'),
                  'started_wall': time.time(), 'start_new_session': kwargs.get('start_new_session', False)}
        try:
            p = psutil.Process(process.pid)
            record.update(create_time=p.create_time(), pgid=os.getpgid(p.pid))
        except psutil.NoSuchProcess:
            record['exited_before_identity_sample'] = True
        with (ROOT / 'child-process-events.jsonl').open('a') as stream: stream.write(json.dumps(record) + '\n')
        return process
    subprocess.Popen = recorded_popen
    worker = None
    try:
        verify(json.loads((ROOT / 'plan.json').read_text())['source_pins'])
        worker = CaptureWorker(JOB)
        worker.prepare()
        raise RuntimeError('Preparation returned without the required capture interception')
    except CaptureComplete as exc:
        worker.update(status='captured_before_refinement', stage='Common preparation captured', message=str(exc))
        write(job / 'provenance.json', worker.provenance)
        if sorted(p.name for p in (ROOT / 'workspace/datasets').iterdir()) != [DATASET]:
            raise RuntimeError('Capture unexpectedly published a dataset')
        if any((job / n).exists() for n in ['prepared.pdb', 'preparation.json', 'loop-refinement-final.npz']):
            raise RuntimeError('Capture crossed the refinement/publication boundary')
        return 0
    except BaseException as exc:
        write(ROOT / 'capture-failure.json', {'type': type(exc).__name__, 'error': str(exc), 'traceback': traceback.format_exc(),
                                             'simulation_ready': False, 'app_published': False})
        if worker is not None:
            worker.job.update(status='capture_failed', error=str(exc))
            write(job / 'status.json', worker.job)
        traceback.print_exc()
        return 1
    finally:
        subprocess.Popen = original_popen


class OwnedTree:
    """Track observed descendant identities, including legacy nested sessions."""
    def __init__(self):
        self.known = {}
        self.samples = []
        self.leader = None
    def record(self, process):
        import psutil
        try:
            with process.oneshot():
                row = dict(pid=process.pid, create_time=process.create_time(), parent_pid=process.ppid(),
                           pgid=os.getpgid(process.pid), command=process.cmdline(),
                           rss_bytes=process.memory_info().rss, status=process.status())
            self.known[(row['pid'], row['create_time'])] = row
            return row
        except (psutil.NoSuchProcess, ProcessLookupError): return None
    def matching(self, row):
        import psutil
        try:
            p = psutil.Process(row['pid'])
            return p if p.create_time() == row['create_time'] and p.status() != psutil.STATUS_ZOMBIE else None
        except psutil.NoSuchProcess: return None
    def sample(self, sample):
        import psutil
        self.leader = sample['group_id']
        try: descendants = [psutil.Process(self.leader), *psutil.Process(self.leader).children(recursive=True)]
        except psutil.NoSuchProcess: descendants = []
        for p in descendants: self.record(p)
        live = []
        for row in list(self.known.values()):
            p = self.matching(row)
            if p is not None:
                latest = self.record(p)
                if latest is not None: live.append(latest)
        rss = sum(r['rss_bytes'] for r in live)
        self.samples.append({'elapsed_seconds': sample['elapsed_seconds'], 'rss_bytes': rss, 'members': live})
        write(ROOT / 'tree-resource.json', {'samples': self.samples, 'peak_sampled_rss_bytes': max(r['rss_bytes'] for r in self.samples)})
        if rss > LIMIT: raise MemoryError('Complete descendant tree exceeded 4 GiB')
        if shutil.disk_usage(ROOT).free < 2 * 1024**3: raise OSError('Capture running disk reserve below 2 GiB')
    def cleanup(self):
        import psutil
        events = ROOT / 'child-process-events.jsonl'
        if events.exists():
            for line in events.read_text().splitlines():
                row = json.loads(line)
                if row.get('parent_pid') == self.leader and 'create_time' in row:
                    self.known[(row['pid'], row['create_time'])] = row
        actions = []
        for sig, seconds in [(signal.SIGTERM, 2), (signal.SIGKILL, 3)]:
            current = [p for row in list(self.known.values()) if (p := self.matching(row)) is not None]
            for p in current:
                for descendant in p.children(recursive=True): self.record(descendant)
            groups = {os.getpgid(p.pid) for p in current if p.is_running()}
            for group in groups:
                members = []
                for p in psutil.process_iter(['pid', 'create_time', 'status']):
                    try:
                        if os.getpgid(p.pid) == group and p.info['status'] != psutil.STATUS_ZOMBIE:
                            if (p.pid, p.info['create_time']) not in self.known:
                                raise RuntimeError('Refusing to signal a group containing an unowned/reused PID')
                            members.append(p.pid)
                    except (psutil.NoSuchProcess, ProcessLookupError): pass
                if members:
                    try: os.killpg(group, sig)
                    except ProcessLookupError: pass
                    actions.append({'pgid': group, 'signal': sig.name, 'verified_pids': members})
            deadline = time.monotonic() + seconds
            while time.monotonic() < deadline and any(self.matching(r) for r in self.known.values()): time.sleep(.05)
        remaining = [r for r in self.known.values() if self.matching(r) is not None]
        write(ROOT / 'descendant-cleanup.json', {'actions': actions, 'no_live_owned_descendants': not remaining,
                                               'remaining': remaining, 'identity_checks': 'PID + creation time + recorded ancestry/group membership'})
        if remaining: raise RuntimeError('Capture left live owned descendants')


def launch():
    setup_environment()
    from backend.loop_process import supervise_loop_child, _uninterrupted_cleanup
    frozen = ROOT / 'implementation/capture_preparation.py'
    frozen.parent.mkdir(exist_ok=True)
    if (ROOT / 'plan.json').exists(): raise ValueError('A native attempt plan already exists')
    shutil.copyfile(__file__, frozen)
    paths = [REPO / ('backend/' + n + '.py') for n in SOURCE_NAMES] + [frozen, REPO / 'pyproject.toml']
    pins = [{'path': str(p), 'sha256': digest(p)} for p in paths]
    previous = Path('/Users/ashujo/.cache/dynamol-research/1cll-ethanol-parameters-v1/installed-runtime-preflight.json')
    if digest(previous) != '80261708282c0162b99fdb75368004f3b2a92a0fafc3995f408e9f53fdfdcb54':
        raise ValueError('Qualified AmberTools preflight source changed')
    amber = json.loads(previous.read_text())
    verify(amber['pins'])
    write(ROOT / 'amber-runtime-recheck.json', {'prior_manifest': str(previous), 'prior_sha256': digest(previous),
          'pins': amber['pins'], 'status': 'same_qualified_parameterization_files'})
    raw = json.loads((ROOT / 'raw-source-provenance.json').read_text())
    verify(raw['files'], 'copy')
    overlay = json.loads((ROOT / 'overlay-provenance.json').read_text())
    verify(overlay['files'], 'copy')
    meta = json.loads((ROOT / 'workspace/datasets' / DATASET / 'metadata.json').read_text())
    if meta.get('preparation') is not None or meta['n_atoms'] != 2322 or meta['n_frames'] != 1:
        raise ValueError('Input is not the frozen unprepared single-frame source')
    free = shutil.disk_usage(ROOT).free
    if free < 4 * 1024**3: raise OSError('Less than 4 GiB free at bounded capture entry')
    plan = {'scope': 'One fresh actual 1UA2 preparation, capture before refinement and publication',
            'source_pins': pins, 'raw_manifest_sha256': digest(ROOT / 'raw-source-provenance.json'),
            'overlay_manifest_sha256': digest(ROOT / 'overlay-provenance.json'), 'imports_sha256': digest(ROOT / 'imports.json'),
            'settings': {'build_missing_residues': True, 'remove_waters': True, 'remove_heterogens': False,
                         'optimize_sidechains': True, 'ph': 7.0, 'seed': 2026},
            'native_attempts': 1, 'timeout_seconds': DEADLINE, 'max_tree_rss_bytes': LIMIT, 'threads': THREADS,
            'entry_disk_bytes': free, 'entry_reserve_bytes': 4 * 1024**3, 'running_reserve_bytes': 2 * 1024**3,
            'sampling_seconds': 1, 'short_lived_peak_limitation': 'Sampled RSS does not measure between-sample peaks',
            'SQM_AM1_BCC_allowed': True, 'MD_allowed': False, 'refinement_allowed': False,
            'chemistry_validation': 'Pending complete validator; hash integrity is not chemical acceptance',
            'app_default_changed': False, 'app_publication_allowed': False}
    write(ROOT / 'plan.json', plan)
    tree = OwnedTree()
    error, code = None, None
    started = time.monotonic()
    try:
        code = supervise_loop_child([str(PYTHON), '-B', str(frozen), '--child'], ROOT / 'supervisor',
                                    environment=dict(os.environ), timeout_seconds=DEADLINE,
                                    max_rss_bytes=LIMIT, on_sample=tree.sample)
    except BaseException as exc:
        error = {'type': type(exc).__name__, 'error': str(exc), 'traceback': traceback.format_exc()}
    finally:
        with _uninterrupted_cleanup():
            try: tree.cleanup()
            except BaseException as exc: error = {'type': type(exc).__name__, 'error': str(exc), 'previous': error}
    try:
        verify(pins); verify(raw['files'], 'copy'); verify(raw['files'], 'source'); verify(amber['pins']); verify(overlay['files'], 'copy')
    except BaseException as exc: error = {'type': type(exc).__name__, 'error': str(exc), 'previous': error}
    captured = ROOT / 'capture-result.json'
    success = code == 0 and error is None and captured.is_file()
    sqm_outputs = sorted(str(p.relative_to(ROOT)) for p in (ROOT / 'workspace/jobs').rglob('sqm.out'))
    write(ROOT / 'result.json', {'status': 'captured_before_refinement' if success else 'failed', 'returncode': code,
          'error': error, 'elapsed_seconds': time.monotonic()-started, 'SQM_output_files': sqm_outputs,
          'QM_claim': 'Native AM1-BCC/SQM was part of actual ligand preparation' if sqm_outputs else 'Inspect native logs; no SQM output observed',
          'snapshot': json.loads(captured.read_text()) if captured.is_file() else None,
          'simulation_ready': False, 'complete_quality_checks_performed': False, 'app_published': False})
    records = []
    for p in sorted(ROOT.rglob('*')):
        if p.is_file() and p.name != 'artifact-manifest.json' and '__pycache__' not in p.parts:
            records.append({'path': str(p.relative_to(ROOT)), 'sha256': digest(p), 'bytes': p.stat().st_size})
    write(ROOT / 'artifact-manifest.json', {'files': records, 'total_bytes': sum(r['bytes'] for r in records)})
    print(json.dumps({'status': 'captured_before_refinement' if success else 'failed', 'result': str(ROOT / 'result.json')}))
    return 0 if success else 1


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--child', action='store_true')
    args = parser.parse_args()
    raise SystemExit(child() if args.child else launch())
