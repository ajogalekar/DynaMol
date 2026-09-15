import hashlib
import json

import pytest

from backend.loop_search import (Attempt, CandidateUnavailable, REQUIRED_CHECKS,
                                 SearchCancelled, SearchExhausted, Validation, ValidatedArtifact,
                                 run_candidate_search)


SOURCE = 'a' * 64
PLAN = [Attempt('fragment', 'database'), Attempt('fallback', 'torsion-mc')]


def artifact(path, content='refined coordinates'):
    path.write_text(content)
    return ValidatedArtifact(path, hashlib.sha256(path.read_bytes()).hexdigest())


def accepted(folder):
    return Validation(dict.fromkeys(REQUIRED_CHECKS, True), artifact=artifact(folder / 'refined.pdb'))


def test_closed_but_refinement_rejected_candidate_reaches_fallback(tmp_path):
    built, validated, events = [], [], []
    def generate(attempt, folder, checkpoint):
        built.append(attempt.name)
        return {'name': attempt.name, 'closed': True, 'source': SOURCE}
    def validate(candidate, folder, checkpoint):
        validated.append(candidate['name'])
        checks = dict.fromkeys(REQUIRED_CHECKS, True)
        checks['backbone_reference'] = candidate['name'] == 'fallback'
        return Validation(checks, artifact=artifact(folder / 'refined.pdb') if checks['backbone_reference'] else None)
    folder = tmp_path / 'run'
    result = run_candidate_search(PLAN, folder, generate, validate, source_sha256=SOURCE,
                                  on_progress=events.append)
    report = json.loads((folder / 'search.json').read_text())
    assert result.path == folder / 'fallback/refined.pdb'
    assert result.path.read_text() == 'refined coordinates'
    assert report['accepted_artifact'] == {'path': str(result.path), 'sha256': result.sha256}
    assert built == validated == ['fragment', 'fallback']
    assert [a['status'] for a in report['attempts']] == ['rejected', 'accepted']
    assert [e['stage'] for e in events] == ['generating', 'validating', 'complete'] * 2


def test_untested_environment_cannot_be_accepted(tmp_path):
    checks = dict.fromkeys(REQUIRED_CHECKS, True)
    checks['retained_environment'] = None
    with pytest.raises(SearchExhausted):
        run_candidate_search(PLAN, tmp_path / 'run', lambda a, *_: a,
                             lambda *_: Validation(checks), source_sha256=SOURCE)
    report = json.loads((tmp_path / 'run/search.json').read_text())
    assert report['accepted_candidate'] is None
    assert len(report['attempts']) == 2


@pytest.mark.parametrize('failure', [SearchCancelled('cancelled'), TimeoutError('budget'), ValueError('atom map changed')])
def test_cancellation_timeout_and_integrity_failure_do_not_start_another_generator(tmp_path, failure):
    called = []
    def generate(attempt, *_):
        called.append(attempt.name)
        return attempt
    def validate(*_):
        raise failure
    with pytest.raises(type(failure)):
        run_candidate_search(PLAN, tmp_path / 'run', generate, validate, source_sha256=SOURCE)
    assert called == ['fragment']
    report = json.loads((tmp_path / 'run/search.json').read_text())
    assert report['status'] in {'cancelled', 'timed_out', 'failed'}
    assert report['accepted_candidate'] is None


def test_unavailable_generator_advances_but_success_stops(tmp_path):
    called = []
    def generate(attempt, *_):
        called.append(attempt.name)
        if attempt.name == 'fragment':
            raise CandidateUnavailable('No closed proposals')
        return attempt
    plan = PLAN + [Attempt('unused', 'other')]
    result = run_candidate_search(plan, tmp_path / 'run', generate, lambda _, folder, __: accepted(folder), source_sha256=SOURCE)
    assert result.path.parent.name == 'fallback'
    assert called == ['fragment', 'fallback']


@pytest.mark.parametrize('checks', [{}, dict.fromkeys(REQUIRED_CHECKS, 1), dict.fromkeys(REQUIRED_CHECKS, 'passed')])
def test_incomplete_or_truthy_nonboolean_decision_is_invalid(checks):
    with pytest.raises(ValueError):
        Validation(checks)


def test_artifacts_are_not_overwritten(tmp_path):
    folder = tmp_path / 'run'; folder.mkdir()
    (folder / 'search.json').write_text('prior failure')
    with pytest.raises(FileExistsError):
        run_candidate_search(PLAN, folder, lambda *_: None, lambda _, private, __: accepted(private), source_sha256=SOURCE)
    assert (folder / 'search.json').read_text() == 'prior failure'


def test_checks_are_copied_and_cannot_be_cleared_to_accept():
    checks = dict.fromkeys(REQUIRED_CHECKS, True)
    checks['geometry'] = False
    reasons = ['Geometry rejected']
    decision = Validation(checks, reasons)
    checks.clear()
    reasons.clear()
    assert set(decision.checks) == REQUIRED_CHECKS
    assert decision.checks['geometry'] is False
    assert decision.reasons == ('Geometry rejected',)
    assert not decision.accepted
    with pytest.raises(TypeError):
        decision.checks['geometry'] = True
    with pytest.raises(AttributeError):
        decision.checks.clear()


def test_success_requires_validated_artifact():
    with pytest.raises(ValueError, match='final validated artifact'):
        Validation(dict.fromkeys(REQUIRED_CHECKS, True))


def test_returns_refined_artifact_instead_of_generator_proposal(tmp_path):
    proposal = None
    def generate(attempt, private, checkpoint):
        nonlocal proposal
        proposal = private / 'proposal.pdb'
        proposal.write_text('unrefined coordinates')
        return proposal
    def validate(candidate, private, checkpoint):
        assert candidate == proposal
        return accepted(private)
    returned = run_candidate_search(PLAN, tmp_path / 'run', generate, validate, source_sha256=SOURCE)
    assert isinstance(returned, ValidatedArtifact)
    assert returned.path != proposal
    assert returned.path.read_text() == 'refined coordinates'
    assert returned.sha256 == hashlib.sha256(returned.path.read_bytes()).hexdigest()


@pytest.mark.parametrize('invalid', ['wrong_hash', 'different_file', 'outside_attempt', 'symlink', 'directory'])
def test_artifact_integrity_failure_is_fatal_and_never_advances(tmp_path, invalid):
    called = []
    def generate(attempt, *_):
        called.append(attempt.name)
        return attempt
    def validate(candidate, private, checkpoint):
        verified = artifact(private / 'refined.pdb')
        if invalid == 'wrong_hash':
            verified = ValidatedArtifact(verified.path, '0' * 64)
        elif invalid == 'different_file':
            other = private / 'proposal.pdb'; other.write_text('unrefined coordinates')
            verified = ValidatedArtifact(other, verified.sha256)
        elif invalid == 'outside_attempt':
            # A valid file/hash from another directory is not this attempt's output.
            verified = artifact(tmp_path / 'other-attempt.pdb')
        elif invalid == 'symlink':
            link = private / 'linked.pdb'; link.symlink_to(verified.path)
            verified = ValidatedArtifact(link, verified.sha256)
        else:
            verified = ValidatedArtifact(private, verified.sha256)
        return Validation(dict.fromkeys(REQUIRED_CHECKS, True), artifact=verified)
    folder = tmp_path / 'run'
    with pytest.raises(ValueError):
        run_candidate_search(PLAN, folder, generate, validate, source_sha256=SOURCE)
    report = json.loads((folder / 'search.json').read_text())
    assert called == ['fragment']
    assert report['status'] == report['attempts'][0]['status'] == 'failed'
    assert report['accepted_candidate'] is report['accepted_artifact'] is None


def test_file_is_rechecked_after_completion_progress_callback(tmp_path):
    saved = None
    def validate(candidate, private, checkpoint):
        nonlocal saved
        result = accepted(private); saved = result.artifact.path
        return result
    def progress(event):
        if event['stage'] == 'complete':
            saved.write_text('changed after validation')
    with pytest.raises(ValueError, match='SHA256'):
        run_candidate_search(PLAN, tmp_path / 'run', lambda a, *_: a, validate,
                             source_sha256=SOURCE, on_progress=progress)
    report = json.loads((tmp_path / 'run/search.json').read_text())
    assert report['accepted_artifact'] is None


def test_relative_artifact_path_is_bound_to_private_attempt(tmp_path):
    def validate(candidate, private, checkpoint):
        result = accepted(private)
        return Validation(result.checks, artifact=ValidatedArtifact('refined.pdb', result.artifact.sha256))
    result = run_candidate_search(PLAN, tmp_path / 'run', lambda a, *_: a, validate, source_sha256=SOURCE)
    assert result.path == tmp_path / 'run/fragment/refined.pdb'
    assert result.path.is_absolute()
