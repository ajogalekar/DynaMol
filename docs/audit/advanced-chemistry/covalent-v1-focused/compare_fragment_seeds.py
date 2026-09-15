"""Research-only alternative ProMod3 seeds; run with the separate native Python.

The preserved context worker provides identity mapping and full seed export.
This module replaces only its candidate selection with a bounded, recorded
choice. Neither a native score nor fragment closure accepts a prepared model.
"""
import argparse
import hashlib
import importlib.util
import json
import math
import os
from pathlib import Path
import time


EXPECTED_WORKER_SHA256 = '400f2e4cf5860773ebb838c64671958d0427d8e48159c5ba9b5007a7aff28975'


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def search(model, modelling, fdb, sdb, sampler, observed, maps, choices):
    if not modelling.IsBackboneScoringSetUp(model):
        modelling.SetupDefaultBackboneScoring(model)
    records = []
    used_choices = set()
    for original in [gap.Copy() for gap in model.gaps]:
        if original.IsTerminal():
            raise ValueError('This experiment only handles requested internal gaps')
        chain_index = original.GetChainIndex()
        name = model.model.chains[chain_index].name
        current = original.Copy()
        maximum = min(original.length + 2, fdb.MaxFragLength() - 2)
        extender = modelling.FullGapExtender(current, model.seqres[chain_index], maximum)
        insertions = modelling.CountEnclosedInsertions(model, original)
        alternatives = []
        first = True
        record = {'gap': str(original), 'contexts': [], 'maximum_candidates': 40}
        while current.length <= maximum and len(alternatives) < 40:
            if not first and not extender.Extend():
                break
            first = False
            if not current.before.IsValid() or not current.after.IsValid():
                break
            if modelling.CountEnclosedInsertions(model, current) != insertions:
                continue
            candidates = modelling.LoopCandidates.FillFromDatabase(
                current.before, current.after, current.full_seq, fdb, sdb, True)
            entry = {'gap': str(current), 'database_candidates': len(candidates)}
            record['contexts'].append(entry)
            if not len(candidates):
                continue
            candidates.ApplyCCD(current.before, current.after, sampler)
            entry['closed_candidates'] = len(candidates)
            if not len(candidates):
                continue
            modelling.FilterCandidates(candidates, model.model, current)
            entry['ring_filtered_candidates'] = len(candidates)
            if not len(candidates):
                continue
            start = current.before.number.num
            scores = modelling.ScoreContainer()
            candidates.CalculateBackboneScores(scores, model, start, chain_index)
            values = scores.LinearCombine(modelling.ScoringWeights.GetWeights())
            for candidate, native_score in zip(candidates, values):
                if len(alternatives) >= 40:
                    break
                distances, backbone = [], []
                for index in range(len(candidate)):
                    row = maps[name][start + index]
                    for atom in ('N', 'CA', 'C', 'O'):
                        point = getattr(candidate, 'Get' + atom)(index)
                        position = [point.x, point.y, point.z]
                        backbone.append({'identity': [name, row['resid'], row['insertion_code'], row['residue'], atom],
                                         'xyz_angstrom': position})
                        old = observed.get((name, start + index, atom))
                        if old is not None and atom != 'O':
                            distances.append(math.dist(old, position) ** 2)
                if not distances or not all(math.isfinite(v) for v in distances + [native_score]):
                    continue
                if not all(math.isfinite(v) for row in backbone for v in row['xyz_angstrom']):
                    raise ValueError('Nonfinite native candidate coordinates')
                identifier = hashlib.sha256(json.dumps([str(current), backbone], sort_keys=True).encode()).hexdigest()
                alternatives.append({'native': candidate.Copy(), 'gap': current.Copy(),
                    'record': {'id': identifier, 'context_gap': str(current),
                               'native_backbone_score': float(native_score),
                               'observed_backbone_squared_displacement': sum(distances),
                               'observed_backbone_rmsd_angstrom': math.sqrt(sum(distances) / len(distances)),
                               'backbone_atoms': backbone}})
        if not alternatives:
            raise ValueError('No native closed alternative found for ' + str(original))
        if len({a['record']['id'] for a in alternatives}) != len(alternatives):
            raise ValueError('Duplicate native candidate IDs; cannot establish unique choices')
        anchor_order = sorted(alternatives, key=lambda a: (
            a['record']['observed_backbone_squared_displacement'], a['record']['native_backbone_score'], a['record']['id']))
        native_order = sorted(alternatives, key=lambda a: (
            a['record']['native_backbone_score'], a['record']['observed_backbone_squared_displacement'], a['record']['id']))
        choice = choices.get(str(original), choices.get('default', {'ranking': 'native', 'rank': 0}))
        ranking, rank = choice['ranking'], choice['rank']
        if ranking not in {'native', 'anchor'} or type(rank) is not int or not 0 <= rank < len(alternatives):
            raise ValueError('Choice is outside the recorded native candidate inventory')
        if str(original) in choices:
            used_choices.add(str(original))
        selected = (native_order if ranking == 'native' else anchor_order)[rank]
        record.update(candidate_count=len(alternatives), choice=choice,
                      selected_id=selected['record']['id'], selected_context_gap=str(selected['gap']),
                      native_order=[a['record']['id'] for a in native_order],
                      anchor_order=[a['record']['id'] for a in anchor_order],
                      candidates=[a['record'] for a in alternatives])
        modelling.InsertLoopClearGaps(model, selected['native'], selected['gap'])
        records.append(record)
    if used_choices != set(choices) - {'default'}:
        raise ValueError('Requested choice did not map to an actual gap')
    return records


def run(worker_path, source, output, choices):
    if digest(worker_path) != EXPECTED_WORKER_SHA256:
        raise ValueError('Preserved research worker changed')
    original_config = json.loads(source.read_text())
    input_pdb = Path(original_config['input_pdb'])
    hashes = {str(path): digest(path) for path in (worker_path, source, input_pdb, Path(__file__))}
    output.mkdir(exist_ok=False)
    config = dict(original_config)
    config.pop('_progress_path', None)
    config['_progress_path'] = str(output / 'loop-model-progress.json')
    (output / 'input.json').write_text(json.dumps(config, indent=2) + '\n')
    (output / 'choices.json').write_text(json.dumps(choices, indent=2) + '\n')
    (output / 'implementation.py').write_text(Path(__file__).read_text())
    spec = importlib.util.spec_from_file_location('preserved_context_worker', worker_path)
    worker = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(worker)
    maps = {chain['chain_id']: {i + 1: r for i, r in enumerate(chain['residues'])} for chain in config['chains']}

    def alternatives_search(model, modelling, fdb, sdb, sampler, extension, observed):
        return search(model, modelling, fdb, sdb, sampler, observed, maps, choices)

    worker.context_fragment_search = alternatives_search
    def bounded_construction(model, modelling, fdb, sdb, sampler, on_progress=lambda report: None,
                             context_search=None):
        report = {'attempts': [], 'max_context_extension_residues': 2,
                  'observed_context_changes_exported': True,
                  'scope': 'Research seed with recorded native context; no app admission.'}
        attempt = {'method': 'bounded_alternative_fragment_database',
                   'context_extension_residues': 2, 'gaps_before': [str(g) for g in model.gaps],
                   'status': 'running'}
        report['attempts'].append(attempt)
        on_progress(report)
        attempt['gap_searches'] = context_search(2)
        attempt.update(status='complete', gaps_after=[str(g) for g in model.gaps])
        report['gaps_after'] = attempt['gaps_after']
        on_progress(report)
        if model.gaps:
            raise ValueError('The bounded alternative search left unresolved gaps')
        return report

    worker.construct_loops = bounded_construction
    started = time.monotonic()
    try:
        result = worker.run(config)
        result['method'] = 'Bounded native alternative selection, retaining up to two observed context residues and 40 closed candidates per unresolved gap. Research coordinates only.'
        result['construction']['alternative_selection_experiment'] = True
        result['construction']['acceptance'] = 'Unaccepted research seed; final original-coordinate, environment, stereochemistry, geometry and native Ramachandran checks are required.'
        (output / 'output.json').write_text(json.dumps(result, indent=2, allow_nan=False) + '\n')
        if hashes != {name: digest(Path(name)) for name in hashes}:
            raise ValueError('Source changed during the native alternative experiment')
        (output / 'provenance.json').write_text(json.dumps({'source_sha256': hashes, 'choices': choices,
            'elapsed_seconds': time.monotonic() - started, 'app_ready': False,
            'physical_model_validated': False}, indent=2) + '\n')
    except Exception as exc:
        (output / 'failure.json').write_text(json.dumps({'type': type(exc).__name__, 'error': str(exc)}, indent=2) + '\n')
        raise


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    for name in ['worker', 'source', 'output', 'choices']:
        parser.add_argument('--' + name, type=Path, required=True)
    args = parser.parse_args()
    os.environ['PM3_OPENMM_CPU_THREADS'] = '2'
    os.environ['OMP_NUM_THREADS'] = '2'
    run(args.worker, args.source, args.output, json.loads(args.choices.read_text()))
