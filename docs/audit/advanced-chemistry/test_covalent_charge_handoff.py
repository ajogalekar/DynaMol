"""Synthetic contract fixtures only: none of these vectors are QM charges."""
import copy
import importlib.util
import json
from pathlib import Path

import pytest

BASE = Path(__file__).resolve().parent
spec = importlib.util.spec_from_file_location('charge_handoff', BASE/'covalent_charge_handoff.py')
handoff = importlib.util.module_from_spec(spec)
spec.loader.exec_module(handoff)
CAP = handoff.ROOT/'build/advanced-chemistry/covalent'


@pytest.fixture
def charge_fixture():
    source = handoff.read(CAP/'capped-adduct.json')
    constraints = handoff.read(CAP/'resp/constraints.json')
    fixed = {r['index']: r['charge_e'] for r in constraints['frozen_atoms']}
    # Labeled synthetic bookkeeping vector; never exported to a biological model.
    remainder = (source['formal_charge']-sum(fixed.values()))/(len(source['atom_map'])-len(fixed))
    charges = [fixed.get(i, remainder) for i in range(len(source['atom_map']))]
    fit = {'constraint_checks_passed': True, 'charges': charges}
    return source, constraints, copy.deepcopy(constraints), fit, list(charges)


def test_joint_vector_is_retained_without_repartition_or_renormalization(charge_fixture):
    value = handoff.validate_charges(*charge_fixture)
    assert value['charges_e'] == charge_fixture[3]['charges']
    assert value['canonical_fixed_count'] == 18
    assert len(value['removed_cap_indices']) == 12
    assert len(value['atom_ids']) == 94
    assert abs(value['retained_charge_sum_e']) < 1e-10


@pytest.mark.parametrize('mutation,match', [
    ('missing', 'Missing or unfinished'), ('nonfinite', 'Invalid fitted'),
    ('fixed', 'Canonical cap/backbone'), ('reordered', 'native RESP stage2'),
    ('sum', 'Joint adduct charge'), ('constraints', 'Canonical charge constraints'),
    ('indices', 'Source atom indices'), ('equivalence', 'methyl equivalence'),
])
def test_invalid_charge_evidence_is_rejected(charge_fixture, mutation, match):
    source, constraints, expected, fit, raw = charge_fixture
    if mutation == 'missing': fit['constraint_checks_passed'] = False
    elif mutation == 'nonfinite': fit['charges'][20] = raw[20] = float('nan')
    elif mutation == 'fixed':
        fit['charges'][0] += .01; fit['charges'][20] -= .01; raw[:] = fit['charges']
    elif mutation == 'reordered': raw.reverse()
    elif mutation == 'sum': fit['charges'][20] += .01; raw[:] = fit['charges']
    elif mutation == 'constraints': constraints['frozen_atoms'][0]['charge_e'] += .01
    elif mutation == 'indices': source['atom_map'][1]['index'] = 0
    elif mutation == 'equivalence':
        fit['charges'][63] += .01; fit['charges'][64] -= .01; raw[:] = fit['charges']
    with pytest.raises(ValueError, match=match):
        handoff.validate_charges(source, constraints, expected, fit, raw)


@pytest.fixture
def terminal_fixture(tmp_path):
    plan = {'cap_graph_sha256': 'graph', 'cap_sdf_sha256': 'sdf',
            'optimizer_input_sha256': 'input', 'parent': str(tmp_path/'parent'),
            'allowed_continuation_sha256': ['code']}
    state = {**{k: plan[k] for k in ('cap_graph_sha256', 'cap_sdf_sha256', 'optimizer_input_sha256')},
             'status': 'research_candidate', 'stage': 'awaiting_full_adduct_torsion_evidence',
             'full_preparation_ready': False, 'script_sha256': 'code', 'optimizer': plan['parent'],
             'stages': {'resp': {'fit': {'constraint_checks_passed': True}}}}
    return state, plan


@pytest.mark.parametrize('key,value', [('status', 'running'), ('status', 'failed'),
    ('stage', 'constrained_native_resp'), ('full_preparation_ready', True),
    ('cap_graph_sha256', 'other'), ('cap_sdf_sha256', 'other'),
    ('optimizer_input_sha256', 'other'), ('script_sha256', 'unknown'), ('optimizer', '/other/parent')])
def test_terminal_parent_binding_is_required(terminal_fixture, key, value):
    state, plan = terminal_fixture
    state[key] = value
    with pytest.raises(ValueError): handoff.validate_terminal_state(state, plan)


def test_terminal_state_is_not_physical_acceptance(terminal_fixture):
    state, plan = terminal_fixture
    handoff.validate_terminal_state(state, plan)
    assert state['full_preparation_ready'] is False


def test_incomplete_source_read_is_rejected(tmp_path):
    path = tmp_path/'empty.json'
    path.write_text('')
    with pytest.raises(ValueError, match='incompletely downloaded'): handoff.read_bytes(path)


def test_pending_run_preserves_failure_and_does_not_assemble(tmp_path):
    plan = tmp_path/'plan.json'
    plan.write_text(json.dumps({'continuation': str(tmp_path/'pending')}))
    output = tmp_path/'handoff'
    result = handoff.run(plan, output, assemble=True)
    assert result['status'] == 'failed'
    assert result['simulation_ready'] is False
    assert result['assembled'] is False
    assert (output/'result.json').exists()
    assert not (output/'hybrid.mol2').exists()
    with pytest.raises(FileExistsError): handoff.run(plan, output, assemble=True)


def test_native_rewrite_changes_only_coordinates_charge_and_metadata(charge_fixture):
    source, _, _, fit, _ = charge_fixture
    original = handoff.read_bytes(CAP/'interface-native-scope-v4/hybrid.mol2').decode()
    coords = [[i+.01, i+.02, i+.03] for i in range(94)]
    rewritten = handoff.rewrite_mol2(original, source, fit['charges'], coords)
    left, right = original.splitlines(), rewritten.splitlines()
    section = None
    for a, b in zip(left, right):
        if a.startswith('@<TRIPOS>'): section = a
        if section == '@<TRIPOS>ATOM' and not a.startswith('@') and a.strip():
            aa, bb = a.split(), b.split()
            assert aa[:2]+aa[5:8]+aa[9:] == bb[:2]+bb[5:8]+bb[9:]
            index = int(bb[0])-1
            assert [float(x) for x in bb[2:5]] == coords[index]
            assert abs(float(bb[8])-fit['charges'][index]) < 1e-10
        elif a != 'bcc': assert a == b
    assert 'USER_CHARGES' in rewritten


def test_native_rewrite_refuses_atom_alias_changes(charge_fixture):
    source, _, _, fit, _ = charge_fixture
    original = handoff.read_bytes(CAP/'interface-native-scope-v4/hybrid.mol2').decode().replace('1 C1 ', '1 OTHER ', 1)
    with pytest.raises(ValueError, match='identity/order'):
        handoff.rewrite_mol2(original, source, fit['charges'], [[0, 0, 0]]*94)


@pytest.mark.parametrize('parameter', ['bond', 'torsion', 'one_four', 'lj', 'mass'])
def test_charge_only_handoff_detects_other_mechanical_changes(parameter):
    import parmed
    native = parmed.load_file(str(CAP/'interface-native-scope-v4/hybrid.prmtop'))
    original = handoff.mechanical_inventory(native)
    for atom in native.atoms: atom.charge += .001
    assert handoff.mechanical_inventory(native) == original
    if parameter == 'bond': native.bonds[0].type.req += .01
    elif parameter == 'torsion': native.dihedrals[0].type.phi_k += .01
    elif parameter == 'one_four': native.dihedrals[0].type.scee += .01
    elif parameter == 'lj': native.atoms[0].epsilon += .01
    elif parameter == 'mass': native.atoms[0].mass += .01
    assert handoff.mechanical_inventory(native) != original
