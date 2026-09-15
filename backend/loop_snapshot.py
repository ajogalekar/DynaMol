"""Immutable loop-search inputs and verified handoff bundles.

These helpers establish serialization and identity integrity, never chemical
validity. Scientific decisions belong to the complete candidate validator.
File maps use explicit relative names; no force-field discovery or fitting is
performed. Raw and prepared inventories deliberately need not be identical.
"""
from __future__ import annotations

from dataclasses import dataclass
import hashlib
import io
import json
import os
from pathlib import Path, PurePosixPath
import posixpath
import stat
from types import MappingProxyType
from typing import Mapping
from xml.etree import ElementTree
import zipfile

import numpy as np

from .loop_search import REQUIRED_CHECKS


_MAX_FILE_BYTES = 512 * 1024**2
_MAX_MANIFEST_BYTES = 32 * 1024**2
_ROLES = frozenset({'coordinates', 'topology', 'parameters', 'provenance'})


@dataclass(frozen=True)
class AtomRecord:
    identity: tuple[str, str, str, str, str]
    element: str
    atomic_number: int


@dataclass(frozen=True)
class Inventory:
    atoms: tuple[AtomRecord, ...]
    # Canonical undirected endpoints, native bond type label and order.
    bonds: tuple[tuple[int, int, str | None, float | None], ...]
    # All source H neighbors are recorded literally; prepared H requires one
    # heavy parent. Raw PDB hydrogen connectivity can legitimately be absent.
    hydrogen_parents: tuple[tuple[int, tuple[int, ...]], ...]


@dataclass(frozen=True)
class FileRecord:
    name: str
    path: Path
    sha256: str
    size: int
    original_path: str


@dataclass(frozen=True)
class Snapshot:
    manifest_path: Path
    sha256: str
    source: Inventory
    prepared: Inventory
    source_xyz_nm: np.ndarray
    prepared_xyz_nm: np.ndarray
    modeled_keys: tuple[tuple[str, str, str, str], ...]
    requested_modeled_keys: tuple[tuple[str, str, str, str], ...]
    # Explicit source-index -> prepared-index or None; renames are never inferred.
    source_to_prepared: tuple[int | None, ...]
    removed_source_identities: tuple[tuple[str, str, str, str, str], ...]
    preparation_provenance: Mapping
    parameter_files: tuple[FileRecord, ...]
    source_files: tuple[FileRecord, ...]


@dataclass(frozen=True)
class VerifiedBundle:
    manifest_path: Path
    sha256: str
    snapshot: Snapshot
    artifacts: Mapping[str, tuple[FileRecord, ...]]
    validator_checks: Mapping[str, bool]


def _hash(data):
    return hashlib.sha256(data).hexdigest()


def _digest(value):
    if not isinstance(value, str) or len(value) != 64 or any(c not in '0123456789abcdef' for c in value):
        raise ValueError('An exact SHA256 is required')
    return value


def _name(value):
    if (not isinstance(value, str) or not value or '\\' in value or '\x00' in value
            or value != PurePosixPath(value).as_posix() or PurePosixPath(value).is_absolute()
            or any(p in {'', '.', '..'} for p in value.split('/'))):
        raise ValueError('Manifest names must be canonical relative paths without escapes')
    return value


def _plain_path(path):
    path = Path(path).absolute()
    if '..' in path.parts:
        raise ValueError('Path traversal is forbidden')
    for component in (path, *path.parents):
        if component.is_symlink():
            raise ValueError('Symlink paths are forbidden')
    return path


def _read(path, expected=None, limit=_MAX_FILE_BYTES):
    path = _plain_path(path)
    fd = os.open(path, os.O_RDONLY | getattr(os, 'O_NOFOLLOW', 0))
    with os.fdopen(fd, 'rb') as handle:
        before = os.fstat(handle.fileno())
        if not stat.S_ISREG(before.st_mode) or before.st_size > limit:
            raise ValueError('Expected a bounded regular file')
        data = handle.read(limit + 1)
        after = os.fstat(handle.fileno())
    current = _plain_path(path).stat()
    identity = lambda s: (s.st_dev, s.st_ino, s.st_size, s.st_mtime_ns)
    if len(data) > limit or identity(before) != identity(after) or identity(after) != identity(current):
        raise ValueError('A captured file changed during reading')
    if expected is not None and _hash(data) != _digest(expected):
        raise ValueError('File SHA256 does not match the frozen manifest')
    return data


def _json(data):
    def unique(pairs):
        result = {}
        for key, value in pairs:
            if key in result:
                raise ValueError('Duplicate JSON field')
            result[key] = value
        return result
    return json.loads(data, object_pairs_hook=unique,
                      parse_constant=lambda value: (_ for _ in ()).throw(ValueError('Nonfinite JSON value')))


def _encode(value):
    return (json.dumps(value, sort_keys=True, indent=2, allow_nan=False) + '\n').encode()


def _freeze(value):
    if isinstance(value, dict):
        return MappingProxyType({key: _freeze(item) for key, item in value.items()})
    if isinstance(value, list):
        return tuple(_freeze(item) for item in value)
    return value


def _key(value, length):
    if (not isinstance(value, (list, tuple)) or len(value) != length
            or any(not isinstance(x, str) for x in value)
            or not value[1] or not value[3] or (length == 5 and not value[4])):
        raise ValueError('Malformed exact atom/residue identity')
    return tuple(value)


def _inventory(topology):
    atoms = list(topology.atoms())
    if [a.index for a in atoms] != list(range(len(atoms))):
        raise ValueError('Topology atoms must have their full ordered indices')
    rows = []
    for a in atoms:
        if a.element is None:
            raise ValueError('Every snapshot atom needs a known element')
        r = a.residue
        rows.append({'identity': [r.chain.id, r.id, (r.insertionCode or '').strip(), r.name, a.name],
                     'element': a.element.symbol, 'atomic_number': a.element.atomic_number})
    bonds = []
    for bond in topology.bonds():
        a, b = bond
        bonds.append([*sorted((a.index, b.index)),
                      str(bond.type) if getattr(bond, 'type', None) is not None else None,
                      getattr(bond, 'order', None)])
    bonds.sort(key=lambda b: b[:2])
    return {'atoms': rows, 'bonds': bonds}


def _parse_inventory(record, *, prepared):
    from openmm.app import element
    if not isinstance(record, dict) or set(record) - {'atoms', 'bonds', 'hydrogen_parents'}:
        raise ValueError('Malformed topology inventory')
    atoms = []
    for row in record['atoms']:
        if set(row) != {'identity', 'element', 'atomic_number'} or type(row['atomic_number']) is not int:
            raise ValueError('Malformed atom element inventory')
        try:
            actual = element.Element.getBySymbol(row['element'])
        except (KeyError, ValueError, TypeError, AttributeError) as exc:
            raise ValueError('Unknown atom element') from exc
        if row['element'] != actual.symbol or row['atomic_number'] != actual.atomic_number:
            raise ValueError('Atom element and atomic number disagree')
        atoms.append(AtomRecord(_key(row['identity'], 5), actual.symbol, row['atomic_number']))
    if not atoms or len({a.identity for a in atoms}) != len(atoms):
        raise ValueError('Empty or duplicate atom identities')
    bonds, seen = [], set()
    neighbors = [[] for _ in atoms]
    for row in record['bonds']:
        if not isinstance(row, (list, tuple)) or len(row) != 4:
            raise ValueError('Malformed bond inventory')
        a, b, kind, order = row
        if type(a) is not int or type(b) is not int or not 0 <= a < b < len(atoms) or (a, b) in seen:
            raise ValueError('Duplicate, invalid or unordered bond identity')
        if kind is not None and (not isinstance(kind, str) or not kind):
            raise ValueError('Malformed bond type')
        if order is not None and (type(order) not in (int, float) or not np.isfinite(order) or order <= 0):
            raise ValueError('Malformed bond order')
        seen.add((a, b)); neighbors[a].append(b); neighbors[b].append(a)
        bonds.append((a, b, kind, order))
    if bonds != sorted(bonds, key=lambda b: b[:2]):
        raise ValueError('Bond inventory is not canonical')
    hydrogens = tuple((i, tuple(sorted(neighbors[i]))) for i, a in enumerate(atoms) if a.atomic_number == 1)
    if prepared and any(len(parents) != 1 or atoms[parents[0]].atomic_number == 1 for _, parents in hydrogens):
        raise ValueError('Each prepared hydrogen requires one exact heavy parent')
    encoded = [[i, list(parents)] for i, parents in hydrogens]
    if 'hydrogen_parents' in record and record['hydrogen_parents'] != encoded:
        raise ValueError('Hydrogen parent inventory differs from the full topology')
    record['hydrogen_parents'] = encoded
    return Inventory(tuple(atoms), tuple(bonds), hydrogens)


def _coordinates(value, count):
    if hasattr(value, 'value_in_unit'):
        from openmm import unit
        value = value.value_in_unit(unit.nanometer)
    xyz = np.asarray(value)
    if xyz.shape != (count, 3) or xyz.dtype.kind not in 'fiu' or not np.isfinite(xyz).all():
        raise ValueError('Coordinates must be finite real atom-count by 3 arrays in nm')
    # Immutable bytes backing prevents callers from re-enabling array writes.
    return np.frombuffer(np.asarray(xyz, dtype='<f8').tobytes(), dtype='<f8').reshape(count, 3)


def _array(data, count):
    """Check the NPY header before any allocation using its declared shape."""
    stream = io.BytesIO(data)
    version = np.lib.format.read_magic(stream)
    if version not in {(1, 0), (2, 0)}:
        raise ValueError('Unsupported coordinate NPY format')
    reader = np.lib.format.read_array_header_1_0 if version == (1, 0) else np.lib.format.read_array_header_2_0
    shape, fortran, dtype = reader(stream, max_header_size=4096)
    if shape != (count, 3) or dtype != np.dtype('<f8'):
        raise ValueError('Expected an exact atom-count by 3 real float64 coordinate array')
    if len(data) - stream.tell() != count * 3 * 8:
        raise ValueError('Coordinate payload size differs from its declared array')
    value = np.frombuffer(data, dtype=dtype, offset=stream.tell()).reshape(shape, order='F' if fortran else 'C')
    return _coordinates(value, count)


def _coverage(source, prepared, modeled, requested):
    modeled, requested = tuple(_key(k, 4) for k in modeled), tuple(_key(k, 4) for k in requested)
    if (not requested or len(set(modeled)) != len(modeled) or len(set(requested)) != len(requested)
            or set(modeled) != set(requested)):
        raise ValueError('Every requested missing residue requires exactly one modeled identity')
    source_residues = {a.identity[:4] for a in source.atoms}
    if set(modeled) & source_residues:
        raise ValueError('A modeled missing residue already has source atoms')
    for key in modeled:
        names = {a.identity[4] for a in prepared.atoms if a.identity[:4] == key and a.atomic_number != 1}
        if not {'N', 'CA', 'C', 'O'} <= names:
            raise ValueError('A requested loop has missing prepared backbone atoms')
    return tuple(sorted(modeled)), tuple(sorted(requested))


def _mapping(source, prepared, values):
    values = tuple(values)
    if len(values) != len(source.atoms):
        raise ValueError('Source mapping must explicitly cover every source atom')
    targets = [i for i in values if i is not None]
    if any(type(i) is not int or not 0 <= i < len(prepared.atoms) for i in targets) or len(set(targets)) != len(targets):
        raise ValueError('Invalid or ambiguous source-to-prepared mapping')
    if any((source.atoms[i].atomic_number, source.atoms[i].element) !=
           (prepared.atoms[j].atomic_number, prepared.atoms[j].element)
           for i, j in enumerate(values) if j is not None):
        raise ValueError('Explicit source mapping changes an atom element')
    return values


def _unmapped_loops(prepared, mapping, modeled):
    if any(prepared.atoms[i].identity[:4] in modeled for i in mapping if i is not None):
        raise ValueError('A wholly missing modeled residue cannot receive a mapped source atom')


def _includes(files):
    """Check recursive XML closure entirely within the supplied file map."""
    graph = {}
    for name, data in files.items():
        if not name.lower().endswith('.xml'):
            continue
        try:
            root = ElementTree.fromstring(data)
        except ElementTree.ParseError as exc:
            raise ValueError('Invalid parameter XML') from exc
        graph[name] = []
        for row in root.iter('Include'):
            target = row.get('file')
            if not target or '\\' in target or PurePosixPath(target).is_absolute():
                raise ValueError('Parameter Include needs an explicit relative target')
            resolved = _name(posixpath.normpath(posixpath.join(posixpath.dirname(name), target)))
            if resolved not in files:
                raise ValueError('Recursive parameter Include is absent from the explicit file list: ' + resolved)
            graph[name].append(resolved)
    visiting, done = set(), set()
    def visit(name):
        if name in visiting:
            raise ValueError('Cyclic parameter Include graph')
        if name in done:
            return
        visiting.add(name)
        for child in graph.get(name, ()):
            visit(child)
        visiting.remove(name); done.add(name)
    for name in graph:
        visit(name)


def _capture(files):
    if not isinstance(files, Mapping) or not files:
        raise ValueError('An explicit nonempty file map is required')
    return { _name(name): (_plain_path(path), _read(path)) for name, path in files.items() }


def _write(root, name, data):
    path = root / _name(name)
    path.parent.mkdir(parents=True, exist_ok=True)
    _plain_path(path)
    with path.open('xb') as handle:
        handle.write(data)


def _files(root, prefix, captured):
    rows = []
    for name, (original, data) in sorted(captured.items()):
        relative = prefix + '/' + name
        _write(root, relative, data)
        rows.append({'name': name, 'path': relative, 'sha256': _hash(data), 'size': len(data), 'original_path': str(original)})
    return rows


def _verify_files(root, rows, prefix):
    result, data = [], {}
    if not isinstance(rows, list) or not rows:
        raise ValueError('Missing declared files')
    for row in rows:
        if set(row) != {'name', 'path', 'sha256', 'size', 'original_path'}:
            raise ValueError('Malformed file record')
        name = _name(row['name'])
        if name in data or row['path'] != prefix + '/' + name or not isinstance(row['original_path'], str):
            raise ValueError('Duplicate or escaping manifest file path')
        content = _read(root / _name(row['path']), row['sha256'])
        if type(row['size']) is not int or row['size'] != len(content):
            raise ValueError('Manifest file size differs')
        data[name] = content
        result.append(FileRecord(name, root / row['path'], row['sha256'], len(content), row['original_path']))
    return tuple(result), data


def _header(kind):
    return {'schema_version': 1, 'kind': kind, 'integrity_only': True,
            'chemistry_validated': False, 'physical_model_validated': False}


def _manifest(path, digest, kind):
    path = _plain_path(path)
    data = _json(_read(path, digest, _MAX_MANIFEST_BYTES))
    if (type(data.get('schema_version')) is not int or data['schema_version'] != 1
            or data.get('kind') != kind or data.get('integrity_only') is not True
            or data.get('chemistry_validated') is not False or data.get('physical_model_validated') is not False):
        raise ValueError('Invalid integrity-only manifest header')
    return path, data


def _exact_tree(root, allowed):
    observed = set()
    for path in root.rglob('*'):
        if path.is_symlink():
            raise ValueError('Symlink artifact is forbidden')
        if path.is_file():
            observed.add(path.relative_to(root).as_posix())
        elif not path.is_dir():
            raise ValueError('Nonregular artifact is forbidden')
    if observed != set(allowed):
        raise ValueError('Missing or undeclared artifact in the frozen bundle')


def create_snapshot(folder, *, source_topology, source_positions, topology, positions,
                    modeled_keys, requested_modeled_keys, source_to_prepared,
                    preparation_provenance, parameter_files, source_files):
    """Freeze current preparation inputs, including explicit retention/rename map.

    ``source_to_prepared`` has one prepared index or None for each ordered raw
    source atom. ``preparation_provenance`` supplies the caller's state and
    retention decisions; this helper does not infer equivalent chemical names.
    File-map keys are relative names preserving explicitly listed XML Includes.
    The new folder must not exist; failed captures are never overwritten.
    """
    source_record, prepared_record = _inventory(source_topology), _inventory(topology)
    source = _parse_inventory(source_record, prepared=False)
    prepared = _parse_inventory(prepared_record, prepared=True)
    source_xyz = _coordinates(source_positions, len(source.atoms))
    prepared_xyz = _coordinates(positions, len(prepared.atoms))
    modeled, requested = _coverage(source, prepared, modeled_keys, requested_modeled_keys)
    mapping = _mapping(source, prepared, source_to_prepared)
    _unmapped_loops(prepared, mapping, modeled)
    if not isinstance(preparation_provenance, Mapping) or not preparation_provenance:
        raise ValueError('Explicit preparation state/retention provenance is required')
    provenance = _json(_encode(dict(preparation_provenance)))
    parameters, sources = _capture(parameter_files), _capture(source_files)
    _includes({name: data for name, (_, data) in parameters.items()})
    root = _plain_path(folder); root.mkdir(parents=True, exist_ok=False)
    record = _header('loop-preparation-snapshot')
    record.update(source=source_record, prepared=prepared_record, modeled_keys=modeled,
                  requested_modeled_keys=requested, source_to_prepared=mapping,
                  preparation_provenance=provenance,
                  parameter_files=_files(root, 'parameters', parameters), source_files=_files(root, 'sources', sources))
    record['coordinates'] = {}
    for name, xyz in [('source', source_xyz), ('prepared', prepared_xyz)]:
        stream = io.BytesIO(); np.save(stream, xyz, allow_pickle=False); data = stream.getvalue()
        path = 'coordinates/' + name + '.npy'; _write(root, path, data)
        record['coordinates'][name] = {'path': path, 'sha256': _hash(data), 'size': len(data), 'units': 'nanometer'}
    data = _encode(record); _write(root, 'snapshot.json', data)
    return load_snapshot(root, expected_sha256=_hash(data))


def load_snapshot(folder, *, expected_sha256):
    """Verify copied bytes and every inventory before returning immutable data."""
    path, record = _manifest(Path(folder) / 'snapshot.json', expected_sha256, 'loop-preparation-snapshot')
    root = path.parent
    if any('hydrogen_parents' not in record[name] for name in ('source', 'prepared')):
        raise ValueError('Both exact hydrogen inventories are required')
    source = _parse_inventory(record['source'], prepared=False)
    prepared = _parse_inventory(record['prepared'], prepared=True)
    modeled, requested = _coverage(source, prepared, record['modeled_keys'], record['requested_modeled_keys'])
    mapping = _mapping(source, prepared, record['source_to_prepared'])
    _unmapped_loops(prepared, mapping, modeled)
    if not isinstance(record['preparation_provenance'], dict) or not record['preparation_provenance']:
        raise ValueError('Missing state/retention provenance')
    parameters, parameter_data = _verify_files(root, record['parameter_files'], 'parameters')
    sources, _ = _verify_files(root, record['source_files'], 'sources')
    _includes(parameter_data)
    coordinates = {}
    if set(record['coordinates']) != {'source', 'prepared'}:
        raise ValueError('Both exact coordinate references are required')
    for name, inventory in [('source', source), ('prepared', prepared)]:
        row = record['coordinates'][name]
        if set(row) != {'path', 'sha256', 'size', 'units'} or row['path'] != f'coordinates/{name}.npy' or row['units'] != 'nanometer':
            raise ValueError('Invalid coordinate artifact declaration')
        data = _read(root / row['path'], row['sha256'])
        if type(row['size']) is not int or row['size'] != len(data):
            raise ValueError('Coordinate artifact size differs')
        coordinates[name] = _array(data, len(inventory.atoms))
    allowed = {'snapshot.json', *(r['path'] for r in record['coordinates'].values()),
               *(r['path'] for r in record['parameter_files']), *(r['path'] for r in record['source_files'])}
    _exact_tree(root, allowed)
    return Snapshot(path, expected_sha256, source, prepared, coordinates['source'], coordinates['prepared'],
                    modeled, requested, mapping, tuple(a.identity for a, target in zip(source.atoms, mapping) if target is None),
                    _freeze(record['preparation_provenance']), parameters, sources)


def _checks(checks):
    if not isinstance(checks, Mapping) or set(checks) != REQUIRED_CHECKS or any(v is not True for v in checks.values()):
        raise ValueError('Every scientific check must be explicitly reported true by the validator')
    return dict(checks)


def _handoff_coordinates(files, atom_count):
    """Final coordinates are explicit nm arrays, separate from topology files."""
    for name, data in files.items():
        if name.endswith('.npz'):
            with zipfile.ZipFile(io.BytesIO(data)) as archive:
                entries = archive.infolist()
                if len(entries) != 1 or entries[0].filename != 'xyz_nm.npy':
                    raise ValueError('Coordinate NPZ must contain only xyz_nm')
                if entries[0].file_size > atom_count * 3 * 8 + 8192:
                    raise ValueError('Coordinate archive exceeds its expected uncompressed array size')
                data = archive.read(entries[0])
        elif not name.endswith('.npy'):
            raise ValueError('Coordinates require NPY or single xyz_nm NPZ artifacts')
        _array(data, atom_count)


def create_accepted_bundle(folder, *, snapshot, artifacts, checks):
    """Copy a validator-approved bundle; a hash does not establish chemistry.

    ``artifacts`` maps coordinates/topology/parameters/provenance to nonempty
    relative-name -> path maps. Parameters must exactly match the snapshot's
    complete parameter file map. The frozen snapshot is copied into the bundle.
    """
    if not isinstance(snapshot, Snapshot):
        raise TypeError('A verified Snapshot is required')
    snapshot = load_snapshot(snapshot.manifest_path.parent, expected_sha256=snapshot.sha256)
    checks = _checks(checks)
    if not isinstance(artifacts, Mapping) or set(artifacts) != _ROLES:
        raise ValueError('Coordinate, topology, parameter and provenance roles are mandatory')
    captured = {role: _capture(files) for role, files in artifacts.items()}
    _handoff_coordinates({name: data for name, (_, data) in captured['coordinates'].items()}, len(snapshot.prepared.atoms))
    if {name: _hash(data) for name, (_, data) in captured['parameters'].items()} != {f.name: f.sha256 for f in snapshot.parameter_files}:
        raise ValueError('Accepted parameter bytes differ from the frozen preparation')
    root = _plain_path(folder); root.mkdir(parents=True, exist_ok=False)
    # Copy only verified snapshot members, then verify the captured copy again.
    for member in snapshot.manifest_path.parent.rglob('*'):
        if member.is_file():
            _write(root, 'snapshot/' + member.relative_to(snapshot.manifest_path.parent).as_posix(), _read(member))
    copied = load_snapshot(root / 'snapshot', expected_sha256=snapshot.sha256)
    record = _header('accepted-loop-bundle')
    record.update(snapshot={'path': 'snapshot/snapshot.json', 'sha256': copied.sha256}, validator_checks=checks,
                  artifacts={role: _files(root, 'artifacts/' + role, files) for role, files in sorted(captured.items())},
                  acceptance_scope='Caller-reported complete validator checks plus verified serialization/identity integrity; no physical-model claim')
    data = _encode(record); _write(root, 'bundle.json', data)
    return verify_accepted_bundle(root / 'bundle.json', expected_sha256=_hash(data))


def verify_accepted_bundle(manifest_path, *, expected_sha256):
    """Recheck every bound member before the caller promotes any prepared output."""
    path, record = _manifest(manifest_path, expected_sha256, 'accepted-loop-bundle')
    root = path.parent
    if set(record['snapshot']) != {'path', 'sha256'} or record['snapshot']['path'] != 'snapshot/snapshot.json':
        raise ValueError('Invalid bound snapshot path')
    snapshot = load_snapshot(root / 'snapshot', expected_sha256=record['snapshot']['sha256'])
    checks = _checks(record['validator_checks'])
    if set(record['artifacts']) != _ROLES:
        raise ValueError('Missing accepted artifact roles')
    verified = {role: _verify_files(root, rows, 'artifacts/' + role) for role, rows in record['artifacts'].items()}
    roles = {role: result[0] for role, result in verified.items()}
    _handoff_coordinates(verified['coordinates'][1], len(snapshot.prepared.atoms))
    if {f.name: f.sha256 for f in roles['parameters']} != {f.name: f.sha256 for f in snapshot.parameter_files}:
        raise ValueError('Handoff parameters differ from the snapshot')
    snapshot_paths = {p.relative_to(root).as_posix() for p in (root / 'snapshot').rglob('*') if p.is_file()}
    allowed = {'bundle.json', *snapshot_paths, *(r['path'] for rows in record['artifacts'].values() for r in rows)}
    _exact_tree(root, allowed)
    return VerifiedBundle(path, expected_sha256, snapshot, MappingProxyType(roles), MappingProxyType(checks))
