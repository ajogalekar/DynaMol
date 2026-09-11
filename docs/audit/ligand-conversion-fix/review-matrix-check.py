"""Rerun retained ligand conversion checks without altering original evidence."""
from __future__ import annotations

import datetime
import hashlib
import importlib.metadata
import json
import shutil
import sys
import tempfile
from pathlib import Path

from rdkit import Chem

ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT))
from backend.ligands import _validate_parameter_conversion


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main():
    parent = ROOT / 'data/integration-checks/ligand-conversion-fix'
    parent.mkdir(parents=True, exist_ok=True)
    output = Path(tempfile.mkdtemp(prefix='retained-matrix-', dir=parent))
    source_files = ['backend/ligands.py', 'backend/forcefield_identity.py', 'backend/prepared_system.py']
    hashes_before = {name: digest(ROOT / name) for name in source_files}
    fixtures = [(path.parent.parent.name, path.parent) for path in sorted(
        (ROOT / 'data/integration-checks/release-chemistry/chemistry-matrix').glob('*/ligand-*/ligand.prmtop'))]
    fixtures.extend(('9AX6-' + json.loads((path.parent / 'provenance.json').read_text())['component_id'], path.parent)
                    for path in sorted((ROOT / 'data/integration-checks/ligand-smoke/9ax6').glob('ligand-*/ligand.prmtop')))
    results = []
    for label, original in fixtures:
        target = output / label
        target.mkdir()
        inputs = {name: digest(original / name) for name in ['ligand.prmtop', 'ligand.inpcrd', 'ligand.xml', 'input.sdf']}
        for name in inputs:
            shutil.copy2(original / name, target / name)
        item = {'fixture': label, 'original': str(original.relative_to(ROOT)),
                'copied_folder': str(target.relative_to(ROOT)), 'input_sha256': inputs}
        try:
            result = _validate_parameter_conversion(target, target / 'ligand.xml')
            mol = next(iter(Chem.SDMolSupplier(str(target / 'input.sdf'), removeHs=False)))
            item.update(passed=result['passed'], atoms=mol.GetNumAtoms(), formal_charge_e=Chem.GetFormalCharge(mol),
                        native_mapping_atoms=len(result['atom_mapping']['atoms']), conformations=result['conformations'],
                        energy_max_abs_kj_mol=max(abs(x['energy_delta_kj_mol']) for x in result['conformations']),
                        force_max_abs_kj_mol_nm=max(x['maximum_force_delta_kj_mol_nm'] for x in result['conformations']),
                        validation_sha256=digest(target / 'conversion-validation.json'))
        except Exception as exc:
            item.update(passed=False, error=str(exc), error_type=type(exc).__name__)
        item['original_inputs_unchanged'] = all(digest(original / name) == hash_ for name, hash_ in inputs.items())
        results.append(item)
        print(label + ': ' + ('PASS' if item['passed'] else 'FAIL ' + item['error']), flush=True)
    report = {'created_at': datetime.datetime.now(datetime.timezone.utc).isoformat(),
              'scope': 'Twelve retained native Amber ligand bundles evaluated with the current production conversion validator and explicit identity adapter. No parameter regeneration, new chemistry, or trajectory convergence claim.',
              'versions': {name: importlib.metadata.version(name) for name in ['openmm', 'parmed', 'numpy', 'rdkit']},
              'source_sha256': hashes_before,
              'source_unchanged_during_check': hashes_before == {name: digest(ROOT / name) for name in source_files},
              'script_sha256': digest(Path(__file__)),
              'platform': 'OpenMM Reference', 'nonbonded_method': 'NoCutoff', 'constraints': None,
              'pose_perturbation_seeds': [2027, 2028], 'pose_perturbation_sigma_angstrom': 0.003,
              'energy_tolerance_kj_mol': 0.0001, 'force_tolerance_kj_mol_nm': 0.001,
              'fixture_count': len(results), 'passed_count': sum(item['passed'] for item in results),
              'original_inputs_unchanged': all(item['original_inputs_unchanged'] for item in results),
              'results': results}
    successes = [x for x in results if x['passed']]
    report['maximum_energy_error_kj_mol'] = max((x['energy_max_abs_kj_mol'] for x in successes), default=None)
    report['maximum_force_error_kj_mol_nm'] = max((x['force_max_abs_kj_mol_nm'] for x in successes), default=None)
    destination = ROOT / 'docs/audit/ligand-conversion-fix/matrix-results.json'
    destination.write_text(json.dumps(report, indent=2) + '\n')
    print(json.dumps({k: report[k] for k in ['fixture_count', 'passed_count', 'maximum_energy_error_kj_mol', 'maximum_force_error_kj_mol_nm', 'original_inputs_unchanged', 'source_unchanged_during_check']}, indent=2))
    return 0 if all(item['passed'] for item in results) and report['source_unchanged_during_check'] else 1


if __name__ == '__main__':
    raise SystemExit(main())
