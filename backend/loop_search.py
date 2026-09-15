"""Bounded candidate orchestration, opt-in until native validation is complete.

A closed chain is only a proposal. Each generator starts again from the same
source; the next strategy runs when complete candidate validation rejects it.
Native work is supplied by supervised callbacks, not imported in this module.
"""
from __future__ import annotations

from dataclasses import dataclass
import hashlib
import json
from pathlib import Path
import re
import time
from types import MappingProxyType
from typing import Mapping


REQUIRED_CHECKS = frozenset({
    'identity', 'retained_environment', 'geometry', 'stereochemistry',
    'backbone_reference', 'observed_displacement', 'serialized_output',
    'parameter_integrity',
})


class CandidateUnavailable(Exception):
    """A supervised generator found no candidate; another strategy may run."""


class SearchCancelled(Exception):
    """Cancellation must never be treated as a reason to try another method."""


class SearchExhausted(Exception):
    """All scheduled attempts failed to supply an accepted complete model."""


@dataclass(frozen=True)
class Attempt:
    name: str
    generator: str


@dataclass(frozen=True)
class ValidatedArtifact:
    """A validator's final coordinate artifact and its SHA256, not a proposal."""
    path: Path | str
    sha256: str

    def __post_init__(self):
        if not isinstance(self.path, (str, Path)) or not str(self.path).strip():
            raise ValueError('A validated artifact requires a file path.')
        if not isinstance(self.sha256, str) or not re.fullmatch(r'[0-9a-f]{64}', self.sha256):
            raise ValueError('A validated artifact requires its SHA256.')
        object.__setattr__(self, 'path', Path(self.path))


@dataclass(frozen=True)
class Validation:
    checks: Mapping[str, bool | None]
    reasons: tuple[str, ...] = ()
    artifact: ValidatedArtifact | None = None

    def __post_init__(self):
        checks = dict(self.checks)
        if set(checks) != REQUIRED_CHECKS:
            raise ValueError('Complete loop validation must report every required check.')
        if any(value is not None and type(value) is not bool for value in checks.values()):
            raise ValueError('Validation checks must be boolean or explicitly untested.')
        reasons = tuple(self.reasons)
        if any(not isinstance(reason, str) or not reason.strip() for reason in reasons):
            raise ValueError('Validation reasons must be nonempty strings.')
        object.__setattr__(self, 'checks', MappingProxyType(checks))
        object.__setattr__(self, 'reasons', reasons)
        if self.artifact is not None and not isinstance(self.artifact, ValidatedArtifact):
            raise TypeError('Validation artifact must be an explicit ValidatedArtifact.')
        if self.accepted and self.artifact is None:
            raise ValueError('Successful validation must identify the final validated artifact.')

    @property
    def accepted(self):
        return all(self.checks.get(name) is True for name in REQUIRED_CHECKS)


def _verify_artifact(artifact, private, checkpoint):
    """Resolve within this attempt and verify bytes at the acceptance boundary."""
    if private.is_symlink():
        raise ValueError('The private attempt directory cannot be replaced by a symlink.')
    path = artifact.path if artifact.path.is_absolute() else private / artifact.path
    root = private.resolve(strict=True)
    resolved = path.resolve(strict=True)
    if not resolved.is_relative_to(root) or resolved == root:
        raise ValueError('Validated artifact must belong to its private attempt directory.')
    # Reject symlinks even when they currently resolve inside the attempt.
    lexical = Path(path.absolute())
    if not lexical.is_relative_to(root):
        raise ValueError('Validated artifact path escapes its private attempt directory.')
    while lexical != root:
        if lexical.is_symlink():
            raise ValueError('Validated artifact paths cannot contain symlinks.')
        lexical = lexical.parent
    if not resolved.is_file():
        raise ValueError('Validated artifact must be an existing regular file.')
    digest = hashlib.sha256()
    with resolved.open('rb') as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b''):
            checkpoint()
            digest.update(chunk)
    if digest.hexdigest() != artifact.sha256:
        raise ValueError('Validated artifact SHA256 does not match the accepted coordinates.')
    return ValidatedArtifact(resolved, artifact.sha256)


def run_candidate_search(attempts, folder, generate, validate, *, source_sha256,
                         check_cancel=lambda: None, on_progress=lambda event: None,
                         max_attempts=8, deadline_seconds=900):
    """Return the verified final ValidatedArtifact or raise with saved evidence.

Callbacks receive (Attempt/candidate, private_folder, checkpoint). They must
supervise native child processes and invoke checkpoint while waiting. Timeout
and cancellation stop the search. Only CandidateUnavailable and explicit
quality rejection advance it; identity/parameter/programming exceptions are
fatal. Successful validators must return the saved final artifact and its hash
in Validation.artifact. The file must be inside their private attempt folder;
relative paths are resolved there. The unrefined generator object is never
returned. This function is not used by the default preparation worker yet.
"""
    attempts = tuple(attempts)
    if not 1 <= max_attempts <= 16 or not 1 <= len(attempts) <= max_attempts:
        raise ValueError('Loop search requires a bounded, nonempty attempt plan.')
    if not 0 < deadline_seconds <= 1800:
        raise ValueError('Loop search deadline must be between zero and 1800 seconds.')
    if not re.fullmatch(r'[0-9a-f]{64}', source_sha256):
        raise ValueError('Loop search requires the immutable source SHA256.')
    if any(not isinstance(a, Attempt) or not re.fullmatch(r'[a-z0-9][a-z0-9_-]{0,63}', a.name)
           or not isinstance(a.generator, str) or not a.generator.strip() for a in attempts):
        raise ValueError('Invalid candidate attempt identity.')
    if len({a.name for a in attempts}) != len(attempts):
        raise ValueError('Candidate attempts require distinct names.')
    folder = Path(folder).resolve()
    folder.mkdir(parents=True, exist_ok=False)
    started = time.monotonic()
    report = {'schema_version': 2, 'source_sha256': source_sha256,
              'plan': [a.__dict__ for a in attempts], 'deadline_seconds': deadline_seconds,
              'status': 'running', 'accepted_candidate': None,
              'accepted_artifact': None, 'attempts': []}

    def save():
        temporary = folder / 'search.json.tmp'
        temporary.write_text(json.dumps(report, indent=2, allow_nan=False) + '\n')
        temporary.replace(folder / 'search.json')

    def checkpoint():
        check_cancel()
        if time.monotonic() - started >= deadline_seconds:
            raise TimeoutError('The overall loop candidate-search budget was exhausted.')

    def progress(row, stage):
        row['stage'] = stage
        save()
        on_progress({'attempt': row['name'], 'generator': row['generator'],
                     'stage': stage, 'number': len(report['attempts']), 'total': len(attempts)})

    save()
    current = None
    try:
        for attempt in attempts:
            checkpoint()
            private = folder / attempt.name
            private.mkdir()
            current = {'name': attempt.name, 'generator': attempt.generator, 'status': 'running'}
            report['attempts'].append(current)
            progress(current, 'generating')
            try:
                candidate = generate(attempt, private, checkpoint)
            except CandidateUnavailable as exc:
                checkpoint()
                current.update(status='unavailable', reason=str(exc))
                progress(current, 'complete')
                continue
            checkpoint()
            progress(current, 'validating')
            decision = validate(candidate, private, checkpoint)
            if not isinstance(decision, Validation):
                raise TypeError('Candidate validation did not return an explicit Validation decision.')
            checkpoint()
            current.update(checks=dict(decision.checks), reasons=list(decision.reasons))
            if not decision.accepted:
                current['status'] = 'rejected'
            progress(current, 'complete')
            if decision.accepted:
                checkpoint()
                artifact = _verify_artifact(decision.artifact, private, checkpoint)
                current.update(status='accepted', artifact={'path': str(artifact.path), 'sha256': artifact.sha256})
                report.update(status='accepted', accepted_candidate=attempt.name,
                              accepted_artifact=current['artifact'])
                save()
                return artifact
        report['status'] = 'exhausted'
        save()
        raise SearchExhausted('No candidate passed all complete loop-preparation checks.')
    except SearchExhausted:
        raise
    except BaseException as exc:
        status = 'cancelled' if isinstance(exc, (SearchCancelled, KeyboardInterrupt)) else 'failed'
        if isinstance(exc, TimeoutError):
            status = 'timed_out'
        report.update(status=status, error=str(exc), error_type=type(exc).__name__)
        if current is not None and current['status'] == 'running':
            current.update(status=status, error=str(exc), error_type=type(exc).__name__)
        save()
        raise
