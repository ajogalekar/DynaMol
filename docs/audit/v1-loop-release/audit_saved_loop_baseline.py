"""Read-only native reference audit of frozen, previously accepted loop outputs.

Run in the isolated OpenMM/CCTBX research runtime. This does not prepare a
structure, execute dynamics, modify an old artifact, or grant app admission.
"""
from __future__ import annotations

import argparse
from collections import Counter
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path
import shutil
import sys
import traceback

import numpy as np
from openmm import app, unit

ROOT = Path(__file__).resolve().parents[3]
OLD = ROOT / 'docs/audit/loop-fallback'
REFERENCE = ROOT / 'docs/audit/advanced-chemistry/covalent-v1-focused'
sys.path.insert(0, str(REFERENCE))
from rama_reference import rama_report, residue_key


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def write(path, value):
    path.write_text(json.dumps(value, indent=2) + '\n')


def atom_key(atom):
    return (*residue_key(atom.residue), atom.name,
            atom.element.symbol if atom.element else None)


def read_json(path, hashes):
    hashes[str(path)] = digest(path)
    return json.loads(path.read_text())


def target_rows(report, selection):
    return [row for row in report['rows'] if tuple(row['residue']) in selection]


def assess(spec, destination):
    destination.mkdir(exist_ok=False)
    hashes = {}
    if spec['kind'] == 'panel':
        record = read_json(OLD / spec['case']['final_report'], hashes)
        if not record['passed']:
            raise ValueError('Frozen previously accepted case is no longer marked passed')
        target = Path(record['workspace']) / 'datasets' / record['prepared_dataset_id']
        source_pdb = Path(record['workspace']) / 'datasets' / record['selected_dataset_id'] / 'topology.pdb'
        loop = record['preparation']['loop_construction']
    else:
        target = Path(spec['trial']['isolated_job_directory'])
        loop = read_json(target / 'loop-construction.json', hashes)
        # The copied repeat job has no source dataset link. Its exact refinement
        # entry still permits attribution to changed versus fixed observed atoms.
        source_pdb = None
    if not loop['accepted']:
        raise ValueError('Original loop output was not accepted')
    artifacts = loop['refinement_artifacts']
    for name in artifacts.values():
        hashes[str(target / name)] = digest(target / name)
        if hashes[str(target / name)] != loop['refinement_artifacts_sha256'][name]:
            raise ValueError(f'Original refinement artifact hash changed: {name}')
    topology_file = target / artifacts['topology']
    topology = app.PDBxFile(str(topology_file)).topology
    before = np.load(target / artifacts['input'])['xyz_nm']
    after = np.load(target / artifacts['final'])['xyz_nm']
    parameters = read_json(target / artifacts['parameters'], hashes)
    atoms = list(topology.atoms())
    keys = [atom_key(a) for a in atoms]
    if len(set(keys)) != len(keys):
        raise ValueError('Duplicate exact atom identities prevent a reliable mapping')
    if before.shape != after.shape or after.shape != (len(atoms), 3):
        raise ValueError('Refinement arrays do not match saved atom inventory')
    if not np.isfinite(before).all() or not np.isfinite(after).all():
        raise ValueError('Nonfinite refinement coordinates')
    modeled = {tuple(row['identity'][:4]) for row in loop['modeler']['worker_report']['atoms']}
    if modeled != set(map(tuple, parameters['modeled_residue_keys'])):
        raise ValueError('Worker/refinement modeled residue identities disagree')
    flanks = set(map(tuple, loop.get('observed_flanks_remodeled', [])))
    relaxation = loop['refinement'].get('flank_relaxation')
    if flanks != set(map(tuple, relaxation['residues'] if relaxation else [])):
        raise ValueError('Recorded movable context identities disagree')
    heavy = [a.index for a in atoms if a.element != app.element.hydrogen]
    changed = {i for i in heavy if not np.array_equal(before[i], after[i])}
    moved_observed = {keys[i][:4] for i in changed if keys[i][:4] not in modeled}
    if not moved_observed <= flanks or modeled & flanks:
        raise ValueError('Observed atom movement is outside the declared context')
    selected = modeled | moved_observed
    available = {residue_key(r) for r in topology.residues()}
    if not selected <= available:
        raise ValueError('Selected residue identity absent from saved topology')
    reports = {}
    for label, xyz in [('refinement_entry', before), ('refinement_final', after)]:
        reports[label] = rama_report(topology, xyz * unit.nanometer)
    final_by_key = {tuple(r['residue']): r for r in reports['refinement_final']['rows']}
    entry_by_key = {tuple(r['residue']): r for r in reports['refinement_entry']['rows']}
    selected_missing = [list(key) for key in sorted(selected - set(final_by_key))]
    pdb_path = target / 'prepared.pdb'
    roundtrip = None
    if pdb_path.exists():
        hashes[str(pdb_path)] = digest(pdb_path)
        pdb = app.PDBFile(str(pdb_path))
        if keys != [atom_key(a) for a in pdb.topology.atoms()]:
            raise ValueError('Serialized PDB atom identities/order differ')
        pdb_xyz = np.asarray(pdb.positions.value_in_unit(unit.nanometer))
        delta = float(np.max(np.abs(after - pdb_xyz)))
        if delta > .00005 + 1e-12:
            raise ValueError('Serialized PDB differs by more than PDB rounding')
        reports['serialized_pdb'] = rama_report(pdb.topology, pdb.positions)
        saved_by_key = {tuple(r['residue']): r for r in reports['serialized_pdb']['rows']}
        mismatches = [{'residue': list(k), 'npz': final_by_key[k]['classification'],
                       'pdb': saved_by_key.get(k, {}).get('classification')}
                      for k in selected if k in final_by_key and
                      final_by_key[k]['classification'] != saved_by_key.get(k, {}).get('classification')]
        roundtrip = {'maximum_component_rounding_nm': delta,
                     'selected_classification_changes': mismatches}
    source = None
    if source_pdb is not None:
        hashes[str(source_pdb)] = digest(source_pdb)
        pdb = app.PDBFile(str(source_pdb))
        source_atoms = list(pdb.topology.atoms())
        source_keys = [atom_key(a) for a in source_atoms]
        if len(source_keys) != len(set(source_keys)):
            raise ValueError('Source atom identities ambiguous')
        source = (dict(zip(source_keys, range(len(source_keys)))),
                  np.asarray(pdb.positions.value_in_unit(unit.nanometer)))
    classified = []
    affected = set(selected)
    for row in reports['refinement_final']['rows']:
        key = tuple(row['residue'])
        # phi/psi atoms cover the adjoining peptide carbon and nitrogen; source
        # attribution never treats a newly defined gap-boundary torsion as an
        # unchanged experimental torsion merely because its own residue is fixed.
        defining = set(row['phi_atoms'] + row['psi_atoms'])
        touches_modeled = any(keys[i][:4] in modeled for i in defining)
        touches_moved = bool(defining & changed)
        if touches_modeled or touches_moved:
            affected.add(key)
        if key in modeled:
            category = 'modeled'
        elif key in moved_observed:
            category = 'moved_observed_context'
        elif touches_modeled or touches_moved:
            category = 'fixed_observed_boundary_with_new_or_changed_torsion'
        else:
            category = 'unchanged_observed_during_loop_refinement'
        original_match = None
        if source is not None and key not in modeled:
            lookup, source_xyz = source
            present = all(keys[i] in lookup for i in defining)
            original_match = present and all(
                np.max(np.abs(after[i] - source_xyz[lookup[keys[i]]])) <= 1e-12
                for i in defining)
        entry = entry_by_key[key]
        classified.append({**row, 'scope_category': category,
                           'refinement_entry_classification': entry['classification'],
                           'refinement_entry_score': entry['score'],
                           'all_phi_psi_atoms_present_and_unchanged_from_imported_source': original_match})
    reports['refinement_final']['rows'] = classified
    reports['refinement_final']['outliers'] = [r for r in classified if r['classification'] == 'outlier']
    # Retain full reference output (including unrelated observed outliers) to
    # allow independent review without counting those as modeled failures.
    for label, report in reports.items():
        write(destination / (label + '.json'), report)
    final_selected = target_rows(reports['refinement_final'], selected)
    final_affected = target_rows(reports['refinement_final'], affected)
    result = {
        'id': spec['id'], 'kind': spec['kind'], 'status': 'audited',
        'original_target': str(target), 'modeled_residues': sorted(modeled),
        'declared_mobile_context_residues': sorted(flanks),
        'actually_moved_observed_context_residues': sorted(moved_observed),
        'modeled_count': len(modeled), 'selected_count': len(selected),
        'selected_unavailable': selected_missing,
        'selected_outliers': [r for r in final_selected if r['classification'] == 'outlier'],
        'affected_boundary_outliers': [r for r in final_affected if r['classification'] == 'outlier'
                                       and tuple(r['residue']) not in selected],
        'unchanged_observed_outliers': [r for r in classified if r['classification'] == 'outlier'
                                        and tuple(r['residue']) not in affected],
        'selected_native_reference_pass': not selected_missing and bool(final_selected)
                                          and all(r['classification'] != 'outlier' for r in final_selected),
        'selected_and_affected_boundary_native_reference_pass': not selected_missing and bool(final_affected)
                                          and all(r['classification'] != 'outlier' for r in final_affected),
        'refinement_entry_selected_outlier_count': sum(r['classification'] == 'outlier'
                                    for r in target_rows(reports['refinement_entry'], selected)),
        'serialized_pdb_comparison': roundtrip,
        'source_comparison_limit': None if source is not None else
             'Repeat job lacks a saved selected-source link; unchanged refers to refinement entry only.',
        'native_reference_version': reports['refinement_final']['cctbx_version'],
        'native_reference_library_sha256': reports['refinement_final']['native_library_sha256'],
        'artifact_sha256': hashes, 'original_artifacts_unchanged': True,
        'app_admission_changed': False, 'physical_model_validated': False,
    }
    if hashes != {name: digest(name) for name in hashes}:
        raise ValueError('Original source changed during read-only audit')
    write(destination / 'result.json', result)
    return result


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=False)
    summary_path = OLD / 'summary.json'
    summary = json.loads(summary_path.read_text())
    manifest = OLD / 'benchmark-manifest.json'
    if digest(manifest) != summary['frozen_manifest_sha256']:
        raise ValueError('Frozen panel manifest changed')
    specs = [{'id': c['id'], 'kind': 'panel', 'case': c} for c in summary['cases']
             if c['final_status'] == 'prepared']
    repeat_path = OLD / '8k5r-full-prep-final/summary.json'
    repeats = json.loads(repeat_path.read_text())
    specs += [{'id': f"8K5R-repeat-{t['seed']}", 'kind': 'repeat', 'trial': t}
              for t in repeats['trials'] if t['status'] == 'completed']
    sources = {str(p): digest(p) for p in [Path(__file__), summary_path, manifest,
                                         repeat_path, REFERENCE / 'rama_reference.py']}
    shutil.copyfile(__file__, args.output / 'implementation.py')
    shutil.copyfile(REFERENCE / 'rama_reference.py', args.output / 'reference-implementation.py')
    write(args.output / 'plan.json', {'created_at': datetime.now(timezone.utc).isoformat(),
          'source_sha256': sources, 'specs': specs,
          'excluded_originally_blocked_cases': [c['id'] for c in summary['cases'] if c['final_status'] != 'prepared']})
    results = []
    for spec in specs:
        print('Auditing', spec['id'], flush=True)
        try:
            result = assess(spec, args.output / spec['id'])
        except Exception as error:
            result = {'id': spec['id'], 'kind': spec['kind'], 'status': 'audit_error',
                      'error': str(error), 'traceback': traceback.format_exc()}
            write(args.output / spec['id'] / 'audit-error.json', result)
        results.append(result)
        print(json.dumps({k: result.get(k) for k in ['id', 'status', 'modeled_count',
              'selected_native_reference_pass', 'selected_and_affected_boundary_native_reference_pass', 'error']}), flush=True)
        write(args.output / 'progress.json', results)
    if sources != {name: digest(name) for name in sources}:
        raise ValueError('Audit implementation or frozen manifest changed during run')
    result = {'completed_at': datetime.now(timezone.utc).isoformat(), 'source_sha256': sources,
              'original_panel_size': len(summary['cases']), 'audited_panel_case_count': len(specs)-2,
              'audit_error_count': sum(r['status'] != 'audited' for r in results),
              'results': results, 'scope': 'Saved static structure reference audit only; no native-loop accuracy, MD stability, or app admission claim.'}
    write(args.output / 'result.json', result)


if __name__ == '__main__':
    main()
