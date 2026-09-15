#!/usr/bin/env python3
"""One configured private full preparation, with validated coordinate handoff.

No prepared archive is an input. Snapshot hashes establish identity/serialization,
not chemical or geometric acceptance. Native AM1-BCC/SQM may run as part of the
actual existing ligand preparation. The app's default and installed app are untouched.
"""
from __future__ import annotations
import argparse
import hashlib
import importlib
import re
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

SOURCE_REPO = Path('/Users/ashujo/Documents/Science/DynaMol')
PYTHON = Path('/Users/ashujo/.cache/dynamol-runtimes/covalent-v1-focused/bin/python')
AMBER = Path('/Users/ashujo/Library/Application Support/DynaMol/Runtimes/ambertools-b96b42253d7d46c5')
PROMOD = Path('/Users/ashujo/Library/Application Support/DynaMol/Runtimes/promod3-dfd2a559551cd2f2')
ROOT = REPO = DATASET = JOB = CASE = EXPECTED_RAW_ATOMS = SEED = CONFIGURATION = None
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


def configure(config_path=None):
    """Explicit JSON configuration, or the five DYNAMOL_CAPTURE_* variables."""
    global SOURCE_REPO, ROOT, REPO, DATASET, JOB, CASE, EXPECTED_RAW_ATOMS, SEED, CONFIGURATION
    config_path = config_path or os.environ.get('DYNAMOL_CAPTURE_CONFIG')
    if config_path:
        value = json.loads(Path(config_path).read_text())
    else:
        value = {key: os.environ.get('DYNAMOL_CAPTURE_' + env) for key, env in
                 [('case', 'CASE'), ('root', 'ROOT'), ('dataset_id', 'DATASET'),
                  ('expected_raw_atoms', 'RAW_ATOMS'), ('seed', 'SEED')]}
        for key in ['expected_raw_atoms', 'seed']:
            if value[key] is not None:
                value[key] = int(value[key])
    if not isinstance(value, dict):
        raise ValueError('Capture configuration must be a JSON object')
    required = {'case', 'root', 'dataset_id', 'expected_raw_atoms', 'seed'}
    if set(value) - required - {'source_repo'} or not required <= set(value):
        raise ValueError('Capture needs case, root, dataset_id, expected_raw_atoms and seed; source_repo is optional')
    for key in ['case', 'dataset_id']:
        if not isinstance(value[key], str) or not re.fullmatch(r'[A-Za-z0-9][A-Za-z0-9_-]{0,79}', value[key]):
            raise ValueError('Invalid single-component ' + key)
    for key in ['expected_raw_atoms', 'seed']:
        if type(value[key]) is not int or value[key] < (1 if key == 'expected_raw_atoms' else 0):
            raise ValueError(key + ' must be an explicit valid integer')
    if value['seed'] > 2**32 - 1:
        raise ValueError('seed must fit the existing NumPy preparation seed')
    for key in ['root', 'source_repo']:
        if key in value and (not isinstance(value[key], str) or not Path(value[key]).expanduser().is_absolute()):
            raise ValueError(key + ' must be an absolute path')
    ROOT = Path(value['root']).expanduser()
    if any(p.is_symlink() for p in (ROOT, *ROOT.parents)):
        raise ValueError('Capture root cannot contain symlinks')
    SOURCE_REPO = Path(value.get('source_repo', str(SOURCE_REPO))).expanduser()
    REPO = ROOT / 'implementation/source'
    CASE, DATASET = value['case'], value['dataset_id']
    EXPECTED_RAW_ATOMS, SEED = value['expected_raw_atoms'], value['seed']
    JOB = 'fresh-' + CASE.lower() + '-capture'
    CONFIGURATION = {**value, 'root': str(ROOT), 'source_repo': str(SOURCE_REPO)}


def freeze_sources():
    """Copy source and bundled data once; configuration keeps data/runtimes explicit."""
    if REPO.exists():
        raise ValueError('Frozen implementation already exists; use a new attempt root')
    paths = sorted((SOURCE_REPO / 'backend').glob('*.py'))
    paths += sorted(p for p in (SOURCE_REPO / 'backend/data').rglob('*') if p.is_file())
    paths += [SOURCE_REPO / 'pyproject.toml']
    records = []
    for path in paths:
        if path.is_symlink() or not path.is_file():
            raise ValueError('Missing regular implementation source: ' + str(path))
        target = REPO / path.relative_to(SOURCE_REPO)
        target.parent.mkdir(parents=True, exist_ok=True)
        shutil.copyfile(path, target)
        checksum = digest(path)
        if digest(target) != checksum:
            raise ValueError('Source changed while freezing: ' + str(path))
        target.chmod(0o444)
        records.append({'source': str(path), 'copy': str(target), 'sha256': checksum})
    write(ROOT / 'frozen-source-provenance.json', {'files': records})


def raw_snapshot_files(raw, sequence_source):
    """Bind every supplied raw file, including all original-source evidence."""
    records = json.loads((ROOT / 'raw-source-provenance.json').read_text())['files']
    files = {}
    for row in records:
        path = Path(row['copy'])
        if not path.is_relative_to(raw) or path.is_symlink() or not path.is_file():
            raise ValueError('Raw manifest copy must be a regular file inside the configured dataset')
        name = 'raw-dataset/' + path.relative_to(raw).as_posix()
        if name in files:
            raise ValueError('Repeated raw manifest file')
        files[name] = path
    actual = {p for p in raw.rglob('*') if p.is_file()}
    if actual != set(files.values()):
        raise ValueError('Raw manifest must enumerate every supplied dataset file')
    if sequence_source is None:
        raise ValueError('No actual sequence evidence was selected for this case')
    sequence_source = Path(sequence_source)
    if sequence_source not in actual:
        raise ValueError('Actual sequence evidence must be one of the pinned raw dataset files')
    files['sequence-evidence' + sequence_source.suffix.lower()] = sequence_source
    return files


def preflight_imports():
    """Imports and native startup occur inside the supervised child; no calculations."""
    records = {}
    for name in ['dimorphite_dl', 'loguru', 'rdkit', 'openmm', 'pdbfixer', 'numpy',
                 'mdtraj', 'parmed', 'gemmi', 'psutil', 'backend.preparation_worker']:
        module = importlib.import_module(name)
        records[name] = {'file': module.__file__, 'sha256': digest(module.__file__)}
        if hasattr(module, '__version__'):
            records[name]['version'] = module.__version__
    write(ROOT / 'imports.json', records)
    from backend.loop_modeling import loop_runtime_status
    runtime = loop_runtime_status()
    if not runtime['available']:
        raise ValueError('Installed ProMod3 preflight failed: ' + str(runtime))
    command = [str(PROMOD / 'bin/python'), '-I', '-B', '-c',
               "import sys,json,promod3,ost,openmm; from promod3 import modelling,loop; "
               "print(json.dumps({'python':sys.version,'promod3':promod3.__version__,"
               "'ost':ost.__version__,'openmm':openmm.__version__}))"]
    started = time.monotonic()
    process = subprocess.Popen(command, stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                               text=True, start_new_session=True, env=dict(os.environ))
    try:
        stdout, stderr = process.communicate(timeout=30)
    except BaseException:
        try:
            os.killpg(process.pid, signal.SIGKILL)
        except ProcessLookupError:
            pass
        process.communicate()
        raise
    result = {'status': 'passed' if process.returncode == 0 else 'failed',
              'returncode': process.returncode, 'stdout': stdout, 'stderr': stderr,
              'command': command, 'elapsed_seconds': time.monotonic() - started,
              'executable_sha256': digest(command[0]), 'native_calculations': False}
    write(ROOT / 'installed-promod-import-preflight.json', result)
    if process.returncode != 0:
        raise ValueError('Installed ProMod3 native imports failed; see preflight evidence')


def setup_environment():
    os.environ.update(DYNAMOL_DATA_DIR=str(ROOT / 'workspace'), DYNAMOL_AMBERTOOLS=str(AMBER),
                      DYNAMOL_PROMOD3=str(PROMOD),
                      DYNAMOL_CAPTURE_CONFIG=str(ROOT / 'capture-configuration.json'),
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


def topology_from_bundle(bundle):
    import importlib.util
    spec = importlib.util.spec_from_file_location('private_snapshot_validator', ROOT / 'implementation/validate_snapshot.py')
    module = importlib.util.module_from_spec(spec); spec.loader.exec_module(module)
    sources = {record.name: record for record in bundle.snapshot.source_files}
    return module.topology_from_snapshot(bundle.snapshot, sources['prepared-topology.cif'].path)[0]


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
    if source.topology.getNumAtoms() != EXPECTED_RAW_ATOMS:
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
             'capture-source.py': Path(__file__), 'capture-plan.json': ROOT / 'plan.json',
             'capture-configuration.json': ROOT / 'capture-configuration.json',
             'frozen-source-provenance.json': ROOT / 'frozen-source-provenance.json',
             'imports.json': ROOT / 'imports.json',
             'installed-promod-import-preflight.json': ROOT / 'installed-promod-import-preflight.json'}
    files.update(raw_snapshot_files(raw, context['source_sequence_path']))
    for p in sorted(REPO.rglob('*')):
        if p.is_file(): files['implementation/' + p.relative_to(REPO).as_posix()] = p
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
            import io
            import numpy as np
            from openmm import unit
            from backend.loop_process import supervise_loop_child
            from backend.loop_snapshot import verify_accepted_bundle, _inventory, _read
            try:
                capture(self, modeller, forcefield, modeled_keys, context)
            except CaptureComplete:
                pass
            recorded = json.loads((ROOT / 'capture-result.json').read_text())
            environment = dict(os.environ)
            environment['PYTHONPATH'] = str(REPO)
            self.update(stage='Validating loop candidates', message='Checking complete candidates before saving; original quality and source-movement limits retained.')
            code = supervise_loop_child([str(PYTHON), '-B', str(ROOT / 'implementation/run_snapshot_search.py'),
                str(ROOT / 'snapshot'), recorded['snapshot_sha256'], str(ROOT / 'candidate-search'),
                '--seed', str(SEED)], ROOT / 'search-supervisor', environment=environment,
                check_cancel=self.check_cancel, timeout_seconds=400)
            if code != 0:
                raise ValueError('No fully validated loop candidate was available; private preparation not published.')
            search = json.loads((ROOT / 'candidate-search/search/search.json').read_text())
            chosen = ROOT / 'candidate-search/search' / search['accepted_candidate']
            decision = json.loads((chosen / 'validation-result.json').read_text())
            bundle = verify_accepted_bundle(decision['accepted_bundle'], expected_sha256=decision['accepted_bundle_sha256'])
            current = _inventory(modeller.topology)
            expected = _inventory(topology_from_bundle(bundle))
            if current != expected or not np.array_equal(np.asarray(modeller.positions.value_in_unit(unit.nanometer)), bundle.snapshot.prepared_xyz_nm):
                raise ValueError('Current worker topology or common coordinates changed before final handoff')
            artifact = bundle.artifacts['coordinates'][0]
            with np.load(io.BytesIO(_read(artifact.path, artifact.sha256)), allow_pickle=False) as archive:
                final = archive['xyz_nm'].copy()
            report_files = {f.name: f for f in bundle.artifacts['provenance']}
            result = json.loads(_read(report_files['result.json'].path, report_files['result.json'].sha256))
            refinement = json.loads(_read(report_files['refinement.json'].path, report_files['refinement.json'].sha256))
            geometry = json.loads(_read(report_files['geometry.json'].path, report_files['geometry.json'].sha256))
            if not result['accepted'] or not all(v is True for v in result['checks'].values()) or not decision['original_source_flank_cap_passed']:
                raise ValueError('Bound final decision does not pass every check')
            refinement['flank_relaxation'] = {'residues': result['observed_context_residue_keys'],
                'geometry': geometry, 'within_displacement_limit': True,
                'maximum_displacement_nm': result['maximum_observed_heavy_displacement_angstrom']/10}
            refinement['accepted_loop_bundle'] = {'path': str(bundle.manifest_path), 'sha256': bundle.sha256,
                'checks': result['checks'], 'original_source_flank_cap_passed': True,
                'selected_candidate': search['accepted_candidate']}
            self.accepted_bundle = bundle
            self.accepted_xyz = final
            return final * unit.nanometer, refinement

        def _before_publish_preparation(self, modeller, preparation):
            import numpy as np
            from openmm import app, unit
            from backend.loop_snapshot import verify_accepted_bundle
            from backend.preparation_worker import atom_key
            verify_accepted_bundle(self.accepted_bundle.manifest_path, expected_sha256=self.accepted_bundle.sha256)
            if not np.array_equal(np.asarray(modeller.positions.value_in_unit(unit.nanometer)), self.accepted_xyz):
                raise ValueError('Worker changed validated coordinates after handoff')
            saved = app.PDBFile(str(self.folder / 'prepared.pdb'))
            if [atom_key(a) for a in saved.topology.atoms()] != [atom_key(a) for a in modeller.topology.atoms()]:
                raise ValueError('Final worker PDB changed the ordered atom identities')
            error = float(np.max(np.linalg.norm(np.asarray(saved.positions.value_in_unit(unit.nanometer)) - self.accepted_xyz, axis=1)))
            if error > .00009:
                raise ValueError('Final worker PDB exceeds the separately validated serialization allowance')
            write(ROOT / 'before-publication-verification.json', {'status': 'passed',
                'coordinates_exactly_unchanged_after_validation': True, 'ordered_pdb_identities_equal': True,
                'maximum_pdb_rounding_nm': error, 'bundle_sha256': self.accepted_bundle.sha256})
    def cancel(signum, frame):
        raise Cancelled('Capture worker cancellation signal ' + str(signum))
    for sig in [signal.SIGTERM, signal.SIGINT]: signal.signal(sig, cancel)
    job = ROOT / 'workspace/jobs' / JOB
    job.mkdir(parents=True)
    shutil.copyfile(ROOT / 'workspace/datasets' / DATASET / 'topology.pdb', job / 'input.pdb')
    settings = {**PREPARATION_DEFAULTS, 'dataset_id': DATASET, 'name': 'Fresh ' + CASE + ' common-preparation capture', 'seed': SEED,
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
        preflight_imports()
        worker = CaptureWorker(JOB)
        worker.prepare()
        dataset = ROOT / 'workspace/datasets' / worker.job['dataset_id']
        destination = dataset / 'loop-validation'
        shutil.copytree(worker.accepted_bundle.manifest_path.parent, destination)
        from backend.loop_snapshot import verify_accepted_bundle
        saved_bundle = verify_accepted_bundle(destination / 'bundle.json', expected_sha256=worker.accepted_bundle.sha256)
        if (dataset / 'prepared.pdb').read_bytes() != (worker.folder / 'prepared.pdb').read_bytes():
            raise ValueError('Dataset did not preserve the final prepared PDB bytes')
        metadata = json.loads((dataset / 'metadata.json').read_text())
        if not metadata['preparation']['simulation_ready']:
            raise ValueError('Worker did not complete prepared metadata')
        for parameter in saved_bundle.artifacts['parameters']:
            local = dataset / parameter.name
            if parameter.name.startswith(('ligands/', 'residue-parameters/')) and (not local.is_file() or digest(local) != parameter.sha256):
                raise ValueError('Dataset lost or changed a prepared parameter sidecar: ' + parameter.name)
        write(ROOT / 'full-preparation-result.json', {'status': 'completed_private_full_preparation',
            'dataset_id': metadata['id'], 'dataset_path': str(dataset), 'n_atoms': metadata['n_atoms'],
            'prepared_pdb_sha256': digest(dataset / 'prepared.pdb'),
            'accepted_bundle_sha256': saved_bundle.sha256, 'parameter_sidecars_preserved': True,
            'bundle_copied_and_verified_in_dataset': True, 'app_default_enabled': False,
            'installed_app_modified': False, 'MD_performed': False})
        return 0
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
    if (ROOT / 'plan.json').exists(): raise ValueError('A native attempt plan already exists')
    if not ROOT.is_dir(): raise ValueError('Input root must already exist')
    if sorted(p.name for p in (ROOT / 'workspace/datasets').iterdir()) != [DATASET]:
        raise ValueError('Capture requires exactly the configured raw dataset')
    write(ROOT / 'capture-configuration.json', CONFIGURATION)
    freeze_sources()
    helper_records = []
    for name in ['run_snapshot_search.py', 'validate_snapshot.py']:
        source = SOURCE_REPO / 'docs/audit/v1-loop-release/fresh-source-integration' / name
        target = ROOT / 'implementation' / name
        content = source.read_text()
        if name == 'run_snapshot_search.py':
            old_root = "REPO = Path('/Users/ashujo/Documents/Science/DynaMol')"
            if content.count(old_root) != 1:
                raise ValueError('Private runner root binding changed')
            content = content.replace(old_root, 'REPO = Path(' + repr(str(REPO)) + ')')
        target.write_text(content)
        helper_records.append({'source': str(source), 'source_sha256': digest(source),
                               'copy': str(target), 'sha256': digest(target),
                               'transformation': 'Bind runner REPO to exact private frozen implementation' if name.startswith('run_') else 'none'})
    write(ROOT / 'helper-source-provenance.json', {'files': helper_records})
    setup_environment()
    from backend.loop_process import supervise_loop_child, _uninterrupted_cleanup
    frozen = ROOT / 'implementation/capture_preparation.py'
    frozen.parent.mkdir(exist_ok=True)
    if (ROOT / 'plan.json').exists(): raise ValueError('A native attempt plan already exists')
    shutil.copyfile(__file__, frozen)
    source_freeze = json.loads((ROOT / 'frozen-source-provenance.json').read_text())
    verify(source_freeze['files'], 'copy')
    paths = [Path(row['copy']) for row in source_freeze['files']] + [frozen, ROOT / 'capture-configuration.json', *[Path(row['copy']) for row in helper_records]]
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
    if meta.get('preparation') is not None or meta['n_atoms'] != EXPECTED_RAW_ATOMS or meta['n_frames'] != 1:
        raise ValueError('Input is not the frozen unprepared single-frame source')
    free = shutil.disk_usage(ROOT).free
    if free < 4 * 1024**3: raise OSError('Less than 4 GiB free at bounded capture entry')
    plan = {'scope': 'One fresh actual ' + CASE + ' preparation with validated private dataset publication',
            'case': CASE, 'dataset_id': DATASET, 'expected_raw_atoms': EXPECTED_RAW_ATOMS,
            'source_pins': pins, 'raw_manifest_sha256': digest(ROOT / 'raw-source-provenance.json'),
            'overlay_manifest_sha256': digest(ROOT / 'overlay-provenance.json'),
            'frozen_source_manifest_sha256': digest(ROOT / 'frozen-source-provenance.json'),
            'import_preflight': 'Recorded inside supervised child before actual preparation',
            'settings': {'build_missing_residues': True, 'remove_waters': True, 'remove_heterogens': False,
                         'optimize_sidechains': True, 'ph': 7.0, 'seed': SEED},
            'native_attempts': 1, 'timeout_seconds': DEADLINE, 'max_tree_rss_bytes': LIMIT, 'threads': THREADS,
            'entry_disk_bytes': free, 'entry_reserve_bytes': 4 * 1024**3, 'running_reserve_bytes': 2 * 1024**3,
            'sampling_seconds': 1, 'short_lived_peak_limitation': 'Sampled RSS does not measure between-sample peaks',
            'SQM_AM1_BCC_allowed': True, 'MD_allowed': False, 'refinement_allowed': True,
            'chemistry_validation': 'Pending complete validator; hash integrity is not chemical acceptance',
            'app_default_changed': False, 'private_dataset_publication_allowed': True, 'installed_app_modified': False}
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
    success = code == 0 and error is None and (ROOT / 'full-preparation-result.json').is_file()
    sqm_outputs = sorted(str(p.relative_to(ROOT)) for p in (ROOT / 'workspace/jobs').rglob('sqm.out'))
    write(ROOT / 'result.json', {'status': 'completed_private_full_preparation' if success else 'failed', 'returncode': code,
          'error': error, 'elapsed_seconds': time.monotonic()-started, 'SQM_output_files': sqm_outputs,
          'QM_claim': 'Native AM1-BCC/SQM was part of actual ligand preparation' if sqm_outputs else 'Inspect native logs; no SQM output observed',
          'snapshot': json.loads(captured.read_text()) if captured.is_file() else None,
          'simulation_ready': success, 'complete_quality_checks_performed': success, 'private_dataset_published': success, 'installed_app_modified': False, 'app_published': False})
    records = []
    for p in sorted(ROOT.rglob('*')):
        if p.is_file() and p.name != 'artifact-manifest.json' and '__pycache__' not in p.parts:
            records.append({'path': str(p.relative_to(ROOT)), 'sha256': digest(p), 'bytes': p.stat().st_size})
    write(ROOT / 'artifact-manifest.json', {'files': records, 'total_bytes': sum(r['bytes'] for r in records)})
    print(json.dumps({'status': 'completed_private_full_preparation' if success else 'failed', 'result': str(ROOT / 'result.json')}))
    return 0 if success else 1


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--child', action='store_true')
    parser.add_argument('--config', help='Absolute JSON file; otherwise DYNAMOL_CAPTURE_CONFIG or explicit DYNAMOL_CAPTURE_* variables')
    args = parser.parse_args()
    try:
        configure(args.config)
        code = child() if args.child else launch()
    except BaseException as exc:
        if ROOT is not None and ROOT.is_dir() and not (ROOT / 'early-failure.json').exists():
            write(ROOT / 'early-failure.json', {'type': type(exc).__name__, 'error': str(exc),
                  'traceback': traceback.format_exc(), 'native_attempted': (ROOT / 'plan.json').exists(),
                  'simulation_ready': False, 'app_published': False})
        raise
    raise SystemExit(code)
