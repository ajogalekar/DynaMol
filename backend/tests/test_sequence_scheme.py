"""Dense declared sequence, engineered author numbering and exact loop identity."""
from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pytest
from openmm import app, unit

from backend import preparation
from backend.sequence_evidence import load_scheme, restoration_plan, restore_residue_identities, recover_observed_insertion_codes

NAMES = ['GLY', 'ILE', 'ASP', 'ARG', 'ALA', 'SER', 'LYS', 'GLY']
IDS = ['240', '241', '242', '243', '272', '273', '274', '275']


def scheme_file(tmp_path, *, names=None, ids=None, insertions=None):
    names, ids = names or NAMES, ids or IDS
    insertions = insertions or ['.'] * len(ids)
    columns = ['asym_id', 'seq_id', 'mon_id', 'pdb_seq_num', 'auth_seq_num', 'pdb_strand_id', 'pdb_ins_code']
    text = 'data_sequence\nloop_\n' + '\n'.join('_pdbx_poly_seq_scheme.' + c for c in columns) + '\n'
    text += '\n'.join(f'L {i + 1} {name} {rid} {rid if i in (0, len(ids) - 1) else "?"} AUTH {code}'
                      for i, (name, rid, code) in enumerate(zip(names, ids, insertions))) + '\n'
    path = tmp_path / 'scheme.cif'
    path.write_text(text)
    return path


def read_scheme(path, names=None, chain='A'):
    return load_scheme(path, {'L': 'AUTH'}, {'AUTH': chain}, {chain}, [SimpleNamespace(chainId=chain, residues=names or NAMES)])


def fixer(rows, observed=(0, 7), *, chain='A'):
    topology = app.Topology()
    new_chain = topology.addChain(chain)
    positions = []
    for offset, index in enumerate(observed):
        row = rows[index]
        residue = topology.addResidue(row['residue'], new_chain, row['resid'], row['insertion_code'])
        for name, element, point in [('N', app.element.nitrogen, [0, 0, 0]), ('CA', app.element.carbon, [.13, .05, 0]),
                                     ('C', app.element.carbon, [.25, 0, .02]), ('O', app.element.oxygen, [.3, -.1, 0])]:
            topology.addAtom(name, element, residue)
            positions.append(np.asarray(point) + [offset * .38, 0, 0])
    return SimpleNamespace(topology=topology, positions=np.asarray(positions) * unit.nanometer,
                           sequences=[SimpleNamespace(chainId=chain, residues=[r['residue'] for r in rows])], sequence_scheme={chain: rows})


def test_dense_scheme_finds_six_residues_without_inventing_number_gap(tmp_path):
    rows = read_scheme(scheme_file(tmp_path))['A']
    model = fixer(rows)
    identity = [(r.id, r.insertionCode, r.name) for r in model.topology.residues()]
    xyz = model.positions.value_in_unit(unit.nanometer).copy()
    preparation.find_missing_residues_preserving_identity(model)
    assert model.missingResidues == {(0, 1): NAMES[1:7]}
    assert [r['resid'] for r in model.missingResidueIdentities[(0, 1)]] == IDS[1:7]
    assert identity == [(r.id, r.insertionCode, r.name) for r in model.topology.residues()]
    np.testing.assert_array_equal(xyz, model.positions.value_in_unit(unit.nanometer))


def test_inserted_author_ids_restore_and_reinspection_builds_nothing_twice(tmp_path):
    rows = read_scheme(scheme_file(tmp_path))['A']
    model = fixer(rows)
    preparation.find_missing_residues_preserving_identity(model)
    plan = restoration_plan(model, model.missingResidues)
    reconstructed = fixer(rows, observed=range(8))
    # Emulate PDBFixer's consecutive IDs immediately before the following275.
    for residue, rid in zip(list(reconstructed.topology.residues())[1:7], range(269, 275)):
        residue.id = str(rid)
    restored = restore_residue_identities(reconstructed.topology, plan)
    assert [r.id for r in reconstructed.topology.residues()] == IDS
    assert [row['resid'] for row in restored] == IDS[1:7]
    preparation.find_missing_residues_preserving_identity(reconstructed)
    assert reconstructed.missingResidues == {}


def test_source_label_author_and_monomer_chain_mapping_is_explicit(tmp_path):
    scheme = read_scheme(scheme_file(tmp_path), chain='B')
    assert set(scheme) == {'B'}
    model = fixer(scheme['B'], chain='B')
    preparation.find_missing_residues_preserving_identity(model)
    assert model.missingResidues == {(0, 1): NAMES[1:7]}
    assert all(row['label_asym_id'] == 'L' for row in model.missingResidueIdentities[(0, 1)])


def test_insertion_codes_remain_part_of_the_exact_identity(tmp_path):
    ids = ['240', '240', *IDS[2:]]
    rows = read_scheme(scheme_file(tmp_path, ids=ids, insertions=['.', 'A', '.', '.', '.', '.', '.', '.']))['A']
    model = fixer(rows, observed=(0, 1, 7))
    preparation.find_missing_residues_preserving_identity(model)
    assert model.missingResidues == {(0, 2): NAMES[2:7]}
    list(model.topology.residues())[1].insertionCode = ''
    with pytest.raises(ValueError, match='exact polymer-scheme identity'):
        preparation.find_missing_residues_preserving_identity(model)


@pytest.mark.parametrize('case', ['duplicate', 'sequence_mismatch', 'observed_mismatch'])
def test_ambiguous_or_mismatched_evidence_fails_closed(tmp_path, case):
    if case == 'duplicate':
        with pytest.raises(ValueError, match='duplicate'):
            read_scheme(scheme_file(tmp_path, ids=['240', '240', *IDS[2:]]))
    elif case == 'sequence_mismatch':
        with pytest.raises(ValueError, match='disagree'):
            read_scheme(scheme_file(tmp_path, names=['ALA', *NAMES[1:]]))
    else:
        model = fixer(read_scheme(scheme_file(tmp_path))['A'])
        list(model.topology.residues())[0].name = 'ALA'
        with pytest.raises(ValueError, match='exact polymer-scheme identity'):
            preparation.find_missing_residues_preserving_identity(model)


def test_missing_modified_residue_is_never_replaced_by_its_parent(tmp_path):
    names = ['GLY', 'TPO', *NAMES[2:]]
    rows = read_scheme(scheme_file(tmp_path, names=names), names=names)['A']
    model = fixer(rows)
    preparation.find_missing_residues_preserving_identity(model)
    assert model.missingResidues[(0, 1)][0] == 'TPO'
    assert 'no parent-residue substitution' in preparation.unsupported_missing_residue_message('A', model.missingResidues[(0, 1)])


def test_retained_4bvn_scheme_has_six_exact_internal_omissions(monkeypatch):
    from backend import config
    root = Path(__file__).resolve().parents[2] / 'data/integration-checks/five-system-openmm'
    if not (root / 'datasets/2d151914298241d9/source.cif').is_file():
        pytest.skip('Optional retained4BVN audit fixture unavailable; synthetic controls remain active')
    monkeypatch.setattr(config, 'DATA_ROOT', root)
    monkeypatch.setattr(config, 'DATASETS_DIR', root / 'datasets')
    model, _ = preparation.current_fixer('2d151914298241d9')
    before = model.positions.value_in_unit(unit.nanometer).copy()
    preparation.find_missing_residues_preserving_identity(model)
    chains = list(model.topology.chains())
    internal = [rows for (chain, position), rows in model.missingResidueIdentities.items()
                if position not in {0, len(list(chains[chain].residues()))}]
    assert len(internal) == 1 and [row['resid'] for row in internal[0]] == IDS[1:7]
    assert [row['residue'] for row in internal[0]] == NAMES[1:7]
    np.testing.assert_array_equal(before, model.positions.value_in_unit(unit.nanometer))


def test_distinct_number_name_pairs_restore_three_codes_and_pdb_roundtrip(tmp_path):
    names = ['GLY', 'TYR', 'GLY', 'LYS', 'ALA', 'GLN']
    ids = ['184', '184', '188', '188', '221', '221']
    rows = read_scheme(scheme_file(tmp_path, names=names, ids=ids, insertions=['.', 'A', '.', 'A', '.', 'A']), names=names)['A']
    model = fixer(rows, observed=range(6))
    for residue in model.topology.residues():
        residue.insertionCode = ''
    before = model.positions.value_in_unit(unit.nanometer).copy()
    recovered = recover_observed_insertion_codes(model)
    assert [(r['resid'], r['residue'], r['insertion_code']) for r in recovered] == [('184', 'TYR', 'A'), ('188', 'LYS', 'A'), ('221', 'GLN', 'A')]
    preparation.find_missing_residues_preserving_identity(model)
    assert model.missingResidues == {}
    np.testing.assert_array_equal(before, model.positions.value_in_unit(unit.nanometer))
    assert recover_observed_insertion_codes(model) == []
    path = tmp_path / 'prepared.pdb'
    with path.open('w') as handle:
        app.PDBFile.writeFile(model.topology, model.positions, handle, keepIds=True)
    loaded = app.PDBFile(str(path))
    assert [(r.id, r.name, r.insertionCode.strip()) for r in loaded.topology.residues()] == [(rid, name, code) for rid, name, code in zip(ids, names, ['', 'A', '', 'A', '', 'A'])]


@pytest.mark.parametrize('case', ['same_name', 'wrong_order'])
def test_insertion_recovery_refuses_ambiguity_or_wrong_order_atomically(tmp_path, case):
    names = ['TYR', 'GLY', 'GLY'] if case == 'same_name' else ['TYR', 'GLY', 'ALA']
    rows = read_scheme(scheme_file(tmp_path, names=names, ids=['100', '184', '184'], insertions=['A', '.', 'A']), names=names)['A']
    model = fixer(rows, observed=range(3) if case == 'same_name' else [0, 2, 1])
    for residue in model.topology.residues():
        residue.insertionCode = ''
    with pytest.raises(ValueError, match='no unique|do not follow'):
        recover_observed_insertion_codes(model)
    assert all(r.insertionCode == '' for r in model.topology.residues())


def test_retained_3ptb_recovers_only_unique_source_codes(monkeypatch):
    from backend import config
    root = Path(__file__).resolve().parents[2] / 'data/integration-checks/five-system-openmm'
    if not (root / 'datasets/59ff9ca69b454df1/source.cif').is_file():
        pytest.skip('Optional retained3PTB audit fixture unavailable; synthetic controls remain active')
    monkeypatch.setattr(config, 'DATA_ROOT', root)
    monkeypatch.setattr(config, 'DATASETS_DIR', root / 'datasets')
    model, _ = preparation.current_fixer('59ff9ca69b454df1')
    assert not getattr(model, 'sequence_scheme_error', None)
    assert [(r['resid'], r['residue'], r['insertion_code']) for r in model.recovered_insertion_codes] == [('184', 'GLY', 'A'), ('188', 'GLY', 'A'), ('221', 'ALA', 'A')]
    preparation.find_missing_residues_preserving_identity(model)
    assert model.missingResidues == {}
