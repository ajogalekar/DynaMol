"""Research-only recovery journal; not integrated with backend/qm_worker.py.

Each record is one *completed SCF and gradient evaluation*, not an accepted
optimizer step or an optimized geometry. A line-search trial can be recorded.
Even the last record remains UNCONVERGED and unaccepted. Full optimization,
constraints, identity/stereo and final-property checks remain the worker's job.

Future integration contract: create a fresh journal from the worker's resolved
validated request and its original input/worker hashes. In the native callback,
pass coords in Bohr, positive dE/dx in Hartree/Bohr, energy in Hartree, and check
the observed scanner atom order/state against stable IDs before passing them.
The requested-method fingerprint binds declared settings; it is not runtime
attestation. Native SCF convergence must be true: the PySCF callback is invoked
BEFORE its own assert_convergence check. No SCF checkpoint is read or endorsed.

Only committed evaluation-XXXXXX.json files count. Temporary files are preserved
on failure and ignored on recovery. Publication uses a flushed/fsynced file and
an atomic, no-clobber hard link on a local filesystem supporting POSIX semantics.
An exception after publication may leave a valid committed record; verify it on
recovery. Power-loss durability still depends on the filesystem. There is no
mutable latest pointer, overwrite, rollback, checkpoint export or acceptance API.

This helper has no numerical dependencies and performs no physical calculations.
"""

import hashlib
import json
import math
import numbers
import os
from pathlib import Path
import re
import tempfile


UNITS = {"coordinates": "bohr", "gradient": "hartree/bohr", "energy": "hartree"}
STATUS = {
    "status": "UNCONVERGED",
    "accepted": False,
    "optimization_converged": False,
    "physical_acceptance": False,
    "scf_checkpoint_same_geometry_validity": "NOT_ESTABLISHED",
    "checkpoint_reuse_authorized": False,
}
IDENTITY_KEYS = {"atom_ids", "elements", "charge", "spin"}
RECORD_KEYS = {
    "schema_version", "kind", "cycle", "identity", "units", "coords_bohr",
    "gradient_hartree_per_bohr", "energy_hartree", "scf_converged",
    "manifest_sha256", "geometry_sha256", "gradient_sha256",
    "requested_method_fingerprint_sha256", *STATUS,
}


def _canonical(value):
    return json.dumps(value, sort_keys=True, separators=(",", ":"),
                      allow_nan=False).encode("utf-8")


def fingerprint(value):
    """SHA-256 of strict canonical JSON, not a file-byte hash."""
    return hashlib.sha256(_canonical(value)).hexdigest()


def _sha(value):
    if not isinstance(value, str) or not re.fullmatch("[a-f0-9]{64}", value):
        raise ValueError("Expected a lowercase SHA-256 digest")
    return value


def _unaccepted(value):
    # bools compare equal to 0/1 in Python; the on-disk contract is stricter.
    return all(type(value.get(k)) is type(v) and value[k] == v for k, v in STATUS.items())


def _identity(value):
    if not isinstance(value, dict) or set(value) != IDENTITY_KEYS:
        raise ValueError("Identity requires exact atom_ids/elements/charge/spin fields")
    ids, elements = value["atom_ids"], value["elements"]
    if (not isinstance(ids, list) or not ids or len(ids) > 1000
            or any(not isinstance(x, str) or not x for x in ids)
            or len(set(ids)) != len(ids)):
        raise ValueError("Atom IDs must be unique ordered nonempty strings")
    if (not isinstance(elements, list) or len(elements) != len(ids)
            or any(not isinstance(x, str) or not re.fullmatch("[A-Z][a-z]?", x)
                   for x in elements)):
        raise ValueError("Elements must be an ordered explicit symbol list")
    for name in ("charge", "spin"):
        if type(value[name]) is not int:
            raise ValueError("Charge and spin must be explicit integers")
    if value["spin"] < 0:
        raise ValueError("Spin must be nonnegative")
    return json.loads(_canonical(value))


def make_binding(validated_request, *, source_input_sha256, worker_source_sha256):
    """Bind a caller-validated resolved request, without selecting any method.

    The fingerprint conservatively includes EVERY non-coordinate/atom-array
    request field, including optimization thresholds and constraints. This is
    deliberately stricter than a chemistry-method label. Raw input-byte and
    resolved-request digests are separately named; supplied hashes are caller
    provenance and not independently verified against external source files.
    """
    request = json.loads(_canonical(validated_request))
    identity = _identity({k: request[k] for k in IDENTITY_KEYS})
    if not isinstance(request.get("method"), str) or not request["method"]:
        raise ValueError("An explicit validated requested method is required")
    if not request.get("basis"):
        raise ValueError("An explicit validated requested basis is required")
    if ("coords_bohr" in request) == ("coords_angstrom" in request):
        raise ValueError("Resolved request requires exactly one coordinate unit")
    key = "coords_bohr" if "coords_bohr" in request else "coords_angstrom"
    _matrix(request[key], len(identity["atom_ids"]), key, flat=False)
    settings = {k: v for k, v in request.items() if k not in
                {"atom_ids", "elements", "coords_bohr", "coords_angstrom"}}
    return {
        "schema_version": 1,
        "kind": "qm_evaluation_journal_binding",
        "identity": identity,
        "source_input_file_sha256": _sha(source_input_sha256),
        "worker_source_file_sha256": _sha(worker_source_sha256),
        "resolved_request_sha256": fingerprint(request),
        "requested_method_and_settings": settings,
        "requested_method_fingerprint_sha256": fingerprint(settings),
        "provenance_scope": "caller_declared_request_not_runtime_attestation",
        **STATUS,
    }


def _number(value):
    if isinstance(value, bool) or not isinstance(value, numbers.Real):
        raise ValueError("Expected a finite real number")
    value = float(value)
    if not math.isfinite(value):
        raise ValueError("Expected a finite real number")
    return value


def _matrix(value, n, label, *, flat):
    if hasattr(value, "tolist"):
        value = value.tolist()
    if not isinstance(value, (tuple, list)):
        raise ValueError(f"{label} must be an N by 3 array")
    if flat and len(value) == 3 * n and all(isinstance(x, numbers.Real) for x in value):
        value = [value[i:i + 3] for i in range(0, len(value), 3)]
    if len(value) != n:
        raise ValueError(f"{label} atom count mismatch")
    result = []
    for row in value:
        if hasattr(row, "tolist"):
            row = row.tolist()
        if not isinstance(row, (tuple, list)) or len(row) != 3:
            raise ValueError(f"{label} must be an N by 3 array")
        result.append([_number(x) for x in row])
    return result


def _write_all(fd, raw):
    remaining = memoryview(raw)
    while remaining:
        written = os.write(fd, remaining)
        if written <= 0:
            raise OSError("Journal write made no progress")
        remaining = remaining[written:]


def _sync_directory(path):
    fd = os.open(path, os.O_RDONLY)
    try:
        os.fsync(fd)
    finally:
        os.close(fd)


def _publish(path, payload):
    """No-clobber commit; preserve temporary evidence on every failure."""
    envelope = {"payload": payload, "payload_sha256": fingerprint(payload)}
    raw = _canonical(envelope) + b"\n"
    fd, temporary = tempfile.mkstemp(prefix=f".{path.name}.", suffix=".tmp", dir=path.parent)
    try:
        _write_all(fd, raw)
        os.fsync(fd)
    finally:
        os.close(fd)
    # link(), unlike replace(), atomically refuses an existing destination.
    os.link(temporary, path)
    _sync_directory(path.parent)
    os.unlink(temporary)
    return {"path": str(path), "file_sha256": hashlib.sha256(raw).hexdigest()}


def _read(path):
    if path.is_symlink() or not path.is_file():
        raise ValueError("Journal artifact must be a regular nonsymlink file")
    size = path.stat().st_size
    if not 0 < size <= 16 * 1024 * 1024:
        raise ValueError("Journal artifact is empty or oversized")
    raw = path.read_bytes()
    if len(raw) != size:
        raise ValueError("Journal artifact read was incomplete")
    envelope = json.loads(raw)
    if (not isinstance(envelope, dict) or set(envelope) != {"payload", "payload_sha256"}
            or fingerprint(envelope["payload"]) != envelope["payload_sha256"]):
        raise ValueError("Journal payload hash mismatch")
    return envelope["payload"]


class EvaluationJournal:
    """One immutable binding per fresh run, one record per native cycle number."""

    def __init__(self, directory, binding, *, create=False):
        self.directory = Path(directory)
        self.binding = json.loads(_canonical(binding))
        self.manifest_sha256 = fingerprint(self.binding)
        if create:
            self.directory.mkdir(parents=True, exist_ok=False)
            _publish(self.directory / "manifest.json", self.binding)
        self._check_manifest()

    def _check_manifest(self):
        actual = _read(self.directory / "manifest.json")
        if actual != self.binding or fingerprint(actual) != self.manifest_sha256:
            raise ValueError("Journal manifest binding mismatch")
        if (actual.get("kind") != "qm_evaluation_journal_binding"
                or type(actual.get("schema_version")) is not int or actual["schema_version"] != 1
                or not _unaccepted(actual)
                or fingerprint(actual["requested_method_and_settings"]) !=
                actual["requested_method_fingerprint_sha256"]):
            raise ValueError("Invalid journal binding")
        _identity(actual["identity"])

    def append(self, *, cycle, coords, gradient, energy, units, observed_identity,
               requested_method_fingerprint_sha256, scf_converged):
        self._check_manifest()
        if type(cycle) is not int or not 1 <= cycle <= 999999:
            raise ValueError("Cycle must be an integer from 1 to 999999")
        if scf_converged is not True:
            raise ValueError("Only completed, natively converged SCF evaluations may be journaled")
        if units != UNITS:
            raise ValueError("Explicit coordinate/gradient/energy units must be bohr/hartree/bohr/hartree")
        identity = _identity(observed_identity)
        if identity != self.binding["identity"]:
            raise ValueError("Observed atom order or molecular state identity mismatch")
        if requested_method_fingerprint_sha256 != self.binding["requested_method_fingerprint_sha256"]:
            raise ValueError("Requested method fingerprint mismatch")
        n = len(identity["atom_ids"])
        xyz = _matrix(coords, n, "coords", flat=True)
        grad = _matrix(gradient, n, "gradient", flat=False)
        payload = {
            "schema_version": 1, "kind": "completed_scf_gradient_evaluation",
            "cycle": cycle, "identity": identity, "units": dict(UNITS),
            "coords_bohr": xyz, "gradient_hartree_per_bohr": grad,
            "energy_hartree": _number(energy), "scf_converged": True,
            "manifest_sha256": self.manifest_sha256,
            "requested_method_fingerprint_sha256": requested_method_fingerprint_sha256,
            "geometry_sha256": fingerprint({"identity": identity, "unit": "bohr", "coords": xyz}),
            "gradient_sha256": fingerprint({"unit": "hartree/bohr", "gradient": grad}),
            **STATUS,
        }
        return _publish(self.directory / f"evaluation-{cycle:06d}.json", payload)

    def records(self):
        """Verify every committed record; refuse corruption, never select a fallback."""
        self._check_manifest()
        records = []
        for path in sorted(self.directory.glob("evaluation-*.json")):
            if not re.fullmatch(r"evaluation-[0-9]{6}\.json", path.name):
                raise ValueError("Malformed committed evaluation name")
            record = _read(path)
            if not isinstance(record, dict) or set(record) != RECORD_KEYS:
                raise ValueError("Invalid evaluation schema")
            _identity(record["identity"])
            if (type(record["schema_version"]) is not int or record["schema_version"] != 1
                    or record["kind"] != "completed_scf_gradient_evaluation"
                    or record["scf_converged"] is not True
                    or type(record["cycle"]) is not int or not 1 <= record["cycle"] <= 999999
                    or path.name != f"evaluation-{record['cycle']:06d}.json"
                    or record["manifest_sha256"] != self.manifest_sha256
                    or record["identity"] != self.binding["identity"]
                    or record["requested_method_fingerprint_sha256"] != self.binding["requested_method_fingerprint_sha256"]
                    or record["units"] != UNITS
                    or not _unaccepted(record)):
                raise ValueError("Evaluation binding or unaccepted-status mismatch")
            n = len(self.binding["identity"]["atom_ids"])
            xyz = _matrix(record["coords_bohr"], n, "coords", flat=False)
            grad = _matrix(record["gradient_hartree_per_bohr"], n, "gradient", flat=False)
            _number(record["energy_hartree"])
            if (record["geometry_sha256"] != fingerprint({"identity": record["identity"], "unit": "bohr", "coords": xyz})
                    or record["gradient_sha256"] != fingerprint({"unit": "hartree/bohr", "gradient": grad})):
                raise ValueError("Evaluation geometry or gradient hash mismatch")
            records.append(record)
        return records
