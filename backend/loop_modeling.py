"""Bounded private loop construction; refinement and quality gates are separate."""
from __future__ import annotations

import hashlib
import json
import os
from pathlib import Path
import re
import signal
import string
import subprocess
import time

import gemmi
import numpy as np
from openmm import app, unit

from . import config
from .modified_residues import PARENT_RESIDUES
from .promod_loop_worker import validate_input
from .residue_identity import STANDARD_PROTEINS, protein_residue_keys, residue_key

_PINS = {'promod3': '3.6.0', 'openstructure': '2.11.1', 'openmm': '8.5.1'}
_TIMEOUT_SECONDS = 180


def _digest(path):
    with Path(path).open('rb') as handle:
        return hashlib.file_digest(handle, 'sha256').hexdigest()


def loop_runtime_status():
    """Read-only installation check, without importing the private native modules."""
    prefix = Path(os.environ.get('DYNAMOL_PROMOD3', config.ROOT / '.tools/promod3')).expanduser().resolve()
    result = {'available': False, 'prefix': str(prefix), 'versions': {}, 'metadata': []}
    try:
        for name, version in _PINS.items():
            records = [(p, json.loads(p.read_text())) for p in (prefix / 'conda-meta').glob(name + '-*.json')]
            records = [(p, r) for p, r in records if r.get('name') == name]
            if len(records) != 1 or records[0][1].get('version') != version:
                raise ValueError(f'Loop building needs the private {name} {version} runtime.')
            path, record = records[0]
            result['versions'][name] = version
            result['metadata'].append({'path': str(path.relative_to(prefix)), 'sha256': _digest(path), 'build': record.get('build')})
        if not os.access(prefix / 'bin/python', os.X_OK) or not (prefix / 'lib/plugins').is_dir():
            raise ValueError('The private loop-building Python or OpenMM plugins are missing.')
        if not any((prefix / 'share/promod3').rglob('*.dat')):
            raise ValueError('The private ProMod3 fragment databases are missing.')
        result['available'] = True
    except (OSError, ValueError) as exc:
        result['error'] = str(exc)
    return result


def _key(atom):
    return (*residue_key(atom.residue), atom.name)


def _xyz(positions):
    return np.asarray(positions.value_in_unit(unit.nanometer) if hasattr(positions, 'value_in_unit') else positions, dtype=float)


def _context(topology, xyz, observed, environment):
    """Only observed residues enter the adapter; context renaming never escapes."""
    proteins = protein_residue_keys(None, topology)
    atoms = list(topology.atoms())
    protein_atoms = [a for a in atoms if residue_key(a.residue) in proteins]
    if len({_key(a) for a in protein_atoms}) != len(protein_atoms):
        raise ValueError('Protein atom identities are ambiguous for loop building.')
    modeled = {residue_key(r) for r in topology.residues() if residue_key(r) in proteins and not any(_key(a) in observed for a in r.atoms())}
    if any(key[3] not in STANDARD_PROTEINS for key in modeled):
        raise ValueError('Only wholly missing standard amino acids can be modeled.')
    target = {_key(a): a.index for a in protein_atoms if residue_key(a.residue) in modeled and a.element and a.element.symbol not in {'H', 'D'}}
    chains, used = [], set()
    for chain in topology.chains():
        residues = [r for r in chain.residues() if residue_key(r) in proteins]
        if not residues:
            continue
        if not re.fullmatch('[A-Za-z0-9]', chain.id) or chain.id in used:
            raise ValueError('Loop building needs unique canonical single-character protein chain IDs.')
        used.add(chain.id)
        rows = []
        residue_ids = set()
        for r in residues:
            key = residue_key(r)
            if key[1:3] in residue_ids:
                raise ValueError('Protein residue numbers and insertion codes are ambiguous for loop building.')
            residue_ids.add(key[1:3])
            if not re.fullmatch(r'-?\d+', key[1]) or str(int(key[1])) != key[1] or not -999 <= int(key[1]) <= 9999 or not re.fullmatch('[A-Za-z0-9]?', key[2]):
                raise ValueError('Protein residue identities cannot be represented exactly in the loop context PDB.')
            letter = gemmi.find_tabulated_residue(PARENT_RESIDUES.get(r.name, r.name)).one_letter_code.upper()
            if letter not in 'ACDEFGHIKLMNPQRSTVWY':
                raise ValueError(f'No unambiguous parent sequence letter is available for {r.name}.')
            rows.append(dict(chain=chain.id, resid=key[1], insertion_code=key[2], residue=r.name, one_letter=letter, observed=key not in modeled))
        chains.append(dict(chain_id=chain.id, residues=rows))
    validate_input(dict(input_pdb='context.pdb', chains=chains, max_res_extension=0))
    context, coordinates, mapping = app.Topology(), [], []
    free = iter(c for c in string.ascii_uppercase + string.ascii_lowercase + string.digits if c not in used)
    sources = [('scaffold', topology, xyz)]
    if environment is not None:
        sources.append(('environment', environment.topology, _xyz(environment.positions)))
    for label, source, positions in sources:
        if positions.shape != (source.getNumAtoms(), 3) or not np.isfinite(positions).all():
            raise ValueError('Loop context coordinates must be finite and match their topology.')
        atom_map = {}
        for chain in source.chains():
            groups = [True, False] if label == 'scaffold' else [False]
            for is_protein in groups:
                residues = [r for r in chain.residues() if (label == 'scaffold' and residue_key(r) in proteins) == is_protein and (not is_protein or residue_key(r) not in modeled)]
                if not residues:
                    continue
                name = chain.id if is_protein else next(free, None)
                if name is None or len(residues) > 9999:
                    raise ValueError('Too many context chains or residues for unambiguous loop building.')
                dest_chain = context.addChain(name)
                for i, r in enumerate(residues, 1):
                    resid, insertion = (r.id, residue_key(r)[2]) if is_protein else (str(i), '')
                    dest = context.addResidue(r.name, dest_chain, resid, insertion)
                    if not is_protein:
                        mapping.append({'source': label, 'chain_index': chain.index, 'source_residue': list(residue_key(r)), 'context_residue': [name, resid, insertion, r.name]})
                    for atom in r.atoms():
                        if not atom.element:
                            raise ValueError('Loop context atoms require known elements.')
                        if atom.element.symbol not in {'H', 'D'}:
                            atom_map[atom] = context.addAtom(atom.name, atom.element, dest)
                            coordinates.append(positions[atom.index])
        for a, b in source.bonds():
            if a in atom_map and b in atom_map:
                context.addBond(atom_map[a], atom_map[b])
    return context, np.asarray(coordinates), chains, target, mapping


def _stop(process):
    try:
        os.killpg(process.pid, signal.SIGTERM)
    except ProcessLookupError:
        pass
    try:
        process.wait(timeout=2)
    except subprocess.TimeoutExpired:
        try:
            os.killpg(process.pid, signal.SIGKILL)
        except ProcessLookupError:
            pass
        process.wait(timeout=5)


def _run(command, folder, runtime, check_cancel):
    env = {k: v for k, v in os.environ.items() if not k.startswith(('PYTHON', 'DYLD_', 'OPENMM', 'PM3_'))}
    env.update(OPENMM_PLUGIN_DIR=str(Path(runtime['prefix']) / 'lib/plugins'), DYNAMOL_CPU_THREADS='2',
               OPENMM_CPU_THREADS='2', PM3_OPENMM_CPU_THREADS='2', OMP_NUM_THREADS='2', OPENBLAS_NUM_THREADS='2', MKL_NUM_THREADS='2')
    check_cancel()
    with (folder / 'loop-model-worker.log').open('w') as log:
        process = subprocess.Popen(command, cwd=folder, env=env, stdout=log, stderr=subprocess.STDOUT, start_new_session=True)
        started = time.monotonic()
        try:
            while process.poll() is None:
                check_cancel()
                if time.monotonic() - started > _TIMEOUT_SECONDS:
                    raise TimeoutError('Loop construction exceeded its 180 second limit. No prepared model was accepted; review the recorded search stages or supply a repaired structure.')
                time.sleep(0.1)
            check_cancel()
            return process.returncode
        except BaseException:
            _stop(process)
            raise


def _transplant(result, target, xyz):
    if not isinstance(result, dict) or result.get('status') != 'candidate':
        raise ValueError('Loop fragment construction failed: ' + str(result.get('error', 'no candidate') if isinstance(result, dict) else 'malformed output'))
    if not isinstance(result.get('atoms'), list):
        raise ValueError('Loop candidate has no atom inventory.')
    updated, seen = xyz.copy(), set()
    for row in result['atoms']:
        if not isinstance(row, dict) or not isinstance(row.get('identity'), list) or len(row['identity']) != 5 or not all(isinstance(x, str) for x in row['identity']):
            raise ValueError('Loop candidate atom identity is malformed.')
        key = tuple(row['identity'])
        if key not in target or key in seen:
            raise ValueError('Loop candidate contains duplicate or unexpected atom identities.')
        coords = np.asarray(row.get('xyz_nm'), dtype=float)
        if coords.shape != (3,) or not np.isfinite(coords).all():
            raise ValueError('Loop candidate coordinates must be finite triples.')
        updated[target[key]] = coords
        seen.add(key)
    if seen != set(target):
        raise ValueError('Loop candidate is missing required heavy atoms.')
    return updated


def generate_loop_model(topology_complete, positions_scaffold, observed_keys, folder, seed, *, environment=None, check_cancel=lambda: None, on_progress=lambda msg: None):
    check_cancel()
    runtime = loop_runtime_status()
    if not runtime['available']:
        raise ValueError(runtime.get('error', 'The private ProMod3 loop runtime is unavailable.'))
    xyz = _xyz(positions_scaffold)
    if xyz.shape != (topology_complete.getNumAtoms(), 3) or not np.isfinite(xyz).all():
        raise ValueError('Loop scaffold coordinates must be finite and match the topology.')
    context, positions, chains, target, mapping = _context(topology_complete, xyz, set(observed_keys), environment)
    folder = Path(folder).resolve()
    folder.mkdir(parents=True, exist_ok=True)
    pdb, request, output = [folder / ('loop-model-' + name) for name in ('context.pdb', 'input.json', 'output.json')]
    if output.exists():
        raise ValueError('Loop model output already exists; preserve it and use a fresh attempt folder.')
    with pdb.open('w') as handle:
        app.PDBFile.writeFile(context, positions * unit.nanometer, handle, keepIds=True)
    request.write_text(json.dumps(dict(input_pdb=str(pdb), chains=chains, max_res_extension=0, requested_seed=seed), indent=2) + '\n')
    progress_file = folder / 'loop-model-progress.json'
    last_message = None
    def check_worker():
        nonlocal last_message
        check_cancel()
        try:
            if not progress_file.is_file() or progress_file.stat().st_size > 128 * 1024:
                return
            progress = json.loads(progress_file.read_text())
            message = progress.get('message')
            if isinstance(message, str) and len(message) <= 512 and message != last_message:
                last_message = message
                on_progress(message)
        except (OSError, ValueError):
            # A progress report is advisory; final bounded output is mandatory.
            pass
    worker = Path(__file__).with_name('promod_loop_worker.py')
    provenance = dict(status='pending', runtime=runtime, requested_seed=seed, seed_applied=False, modeled_heavy_atoms=len(target), environment_context_mapping=mapping,
                      context_limit='Only observed residues plus missing-atom repairs enter the quantized PDB context. Unattached ligands/waters may be omitted by ProMod3; complete-complex refinement and final environment checks are required.',
                      validation_scope='Loop candidate only; requested missing atom inventory never expands. Temporary context residues may move but only missing coordinates are transferred before fixed-heavy refinement and independent geometry/chirality gates.', files={})
    try:
        for path in sorted((Path(runtime['prefix']) / 'share/promod3').rglob('*.dat')):
            check_cancel()
            provenance.setdefault('databases', []).append(dict(path=str(path.relative_to(runtime['prefix'])), bytes=path.stat().st_size, sha256=_digest(path)))
        on_progress('Finding fragment candidates for missing loops…')
        code = _run([str(Path(runtime['prefix']) / 'bin/python'), '-I', '-B', str(worker), str(request), str(output)], folder, runtime, check_worker)
        if not output.exists() or output.stat().st_size > 16 * 1024 * 1024:
            raise ValueError('Loop worker did not produce a bounded candidate report.')
        result = json.loads(output.read_text())
        provenance['worker_report'] = result
        provenance['native_build_attempts'] = result.get('construction', {}).get('attempts', [])
        if code:
            raise ValueError('Loop fragment construction failed: ' + str(result.get('error', f'worker exit {code}')))
        if result.get('versions') != {'promod3': _PINS['promod3'], 'ost': _PINS['openstructure']}:
            raise ValueError('The loop worker loaded unexpected ProMod3/OpenStructure versions.')
        updated = _transplant(result, target, xyz)
        unchanged = np.ones(len(xyz), dtype=bool)
        unchanged[list(target.values())] = False
        if not np.array_equal(updated[unchanged], xyz[unchanged]):
            raise ValueError('Loop candidate changed an observed or non-loop atom.')
        provenance.update(status='candidate', observed_and_non_loop_coordinates_unchanged=True, worker_report=result)
        return updated * unit.nanometer, provenance
    except BaseException as exc:
        provenance.update(status='failed', error=str(exc), error_type=type(exc).__name__)
        raise
    finally:
        try:
            if progress_file.is_file() and progress_file.stat().st_size <= 128 * 1024:
                last_progress = json.loads(progress_file.read_text())
                provenance['last_worker_progress'] = last_progress
                if not provenance.get('native_build_attempts'):
                    provenance['native_build_attempts'] = last_progress.get('attempts', [])
        except (OSError, ValueError):
            pass
        for path in (pdb, request, output, progress_file, worker, folder / 'loop-model-worker.log'):
            if path.is_file():
                provenance['files'][path.name] = dict(sha256=_digest(path), bytes=path.stat().st_size)
        (folder / 'loop-model-provenance.json').write_text(json.dumps(provenance, indent=2, allow_nan=False) + '\n')
