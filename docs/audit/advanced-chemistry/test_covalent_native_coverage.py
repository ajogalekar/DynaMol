import importlib.util
from pathlib import Path

import pytest

BASE = Path(__file__).resolve().parent
spec = importlib.util.spec_from_file_location('coverage_audit', BASE/'covalent_native_coverage.py')
audit = importlib.util.module_from_spec(spec)
spec.loader.exec_module(audit)


def test_unmarked_zero_periodicity_rejected_but_zero_barrier_retained(tmp_path):
    path = tmp_path/'test.frcmod'
    path.write_text('remark\nDIHE\nc3-c3-n -c    0    0.000         0.000           0.000\nc3-c3-n -h1   6    0.000         0.000           2.000\n')
    parsed = audit.parse_frcmod(path)
    assert len(parsed['invalid_terms']) == 1
    assert parsed['entries'][1]['invalid'] is None


def test_all_analogy_scores_retained_without_interpreting_as_accuracy(tmp_path):
    path = tmp_path/'test.frcmod'
    path.write_text('remark\nDIHE\nca-ca-cc-nd   1    0.820       180.000           1.000      same as ca-ce-ce-n2, penalty score=774.0\nc3-c3-n -h1   6    0.000         0.000           2.000      same as X -c3-n -X , penalty score=  0.0\nIMPROPER\nca-ca-ca-h4         1.1          180.0         2.0          Using the default value\n')
    parsed = audit.parse_frcmod(path)
    assert parsed['max_analogy_penalty'] == 774.0
    assert [x['analogy_penalty'] for x in parsed['analogy_terms']] == [774.0, 0.0]
    assert len(parsed['explicit_default_terms']) == 1
    assert not parsed['invalid_terms']


def test_actual_typing_roundtrip_and_coordinate_mutation_rejected(tmp_path):
    folder = audit.ROOT/'build/advanced-chemistry/covalent/native-type-coverage-probe'
    if not folder.exists():
        pytest.skip('Native audit output not present')
    result = audit.verify_typing(folder/'input.sdf', folder/'typed.mol2')
    assert result['atom_count'] == 94
    assert result['charge_columns_are_unusable_zero_placeholders']
    changed = tmp_path/'changed.mol2'
    changed.write_text((folder/'typed.mol2').read_text().replace('-1.8620', '-2.8620', 1))
    with pytest.raises(ValueError, match='coordinates'):
        audit.verify_typing(folder/'input.sdf', changed)


def test_actual_typing_rejects_bond_order_reassignment(tmp_path):
    folder = audit.ROOT/'build/advanced-chemistry/covalent/native-type-coverage-probe'
    if not folder.exists():
        pytest.skip('Native audit output not present')
    lines = (folder/'typed.mol2').read_text().splitlines()
    i = lines.index('@<TRIPOS>BOND') + 1
    values = lines[i].split()
    values[3] = '2' if values[3] == '1' else '1'
    lines[i] = ' '.join(values)
    changed = tmp_path/'changed.mol2'
    changed.write_text('\n'.join(lines)+'\n')
    with pytest.raises(ValueError, match='bond orders'):
        audit.verify_typing(folder/'input.sdf', changed)
