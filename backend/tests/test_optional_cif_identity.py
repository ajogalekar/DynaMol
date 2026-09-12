"""Valid label identifiers remain usable when optional author fields are absent."""
import pytest
from openmm import app

from backend import ligands
from backend.residue_identity import _source_polymer_keys


@pytest.mark.parametrize('author_columns', [True, False])
def test_component_and_polymer_identity_with_optional_author_columns(tmp_path, monkeypatch, author_columns):
    top = app.Topology()
    residue = top.addResidue('ABC', top.addChain('A'), '23')
    top.addAtom('CA', app.element.carbon, residue)
    columns = ['id', 'type_symbol', 'label_atom_id', 'label_comp_id', 'label_asym_id', 'label_seq_id', 'label_entity_id']
    row = ['1', 'C', 'CA', 'ABCDE', 'A', '23', '1']
    if author_columns:
        columns += ['auth_atom_id', 'auth_asym_id', 'auth_seq_id']
        row += ['CA', 'A', '23']
    source = tmp_path / 'source.cif'
    source.write_text("data_optional\nloop_\n_chem_comp.id\n_chem_comp.type\nABCDE .\n"
                      "loop_\n_entity_poly.entity_id\n_entity_poly.type\n1 'polypeptide(L)'\nloop_\n"
                      + '\n'.join('_atom_site.' + c for c in columns) + '\n' + ' '.join(row) + '\n')
    monkeypatch.setattr(ligands, '_source_cif', lambda _: source)
    mapping, path = ligands._component_ids('fixture', top)
    assert path == source
    assert mapping['A:23::ABC'] == 'ABCDE'
    assert _source_polymer_keys('fixture', top) == {('A', '23', '', 'ABC')}
