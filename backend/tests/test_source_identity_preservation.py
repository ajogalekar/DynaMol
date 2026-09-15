"""Native author/label chains and insertion codes survive import and selection."""
import hashlib
import io
from pathlib import Path

import gemmi
import numpy as np
import pytest
from openmm import app, unit

from backend import config, sources, storage, preparation, monomers


@pytest.fixture(autouse=True)
def isolated(tmp_path, monkeypatch):
    monkeypatch.setattr(config, 'DATA_ROOT', tmp_path)
    monkeypatch.setattr(config, 'DATASETS_DIR', tmp_path / 'datasets')
    monkeypatch.setattr(config, 'JOBS_DIR', tmp_path / 'jobs')
    config.DATASETS_DIR.mkdir(); config.JOBS_DIR.mkdir()


def protein(topology, chain, identities):
    coordinates = []
    for i, (number, insertion) in enumerate(identities):
        residue = topology.addResidue('ALA', chain, str(number), insertion)
        for name, symbol, point in [('N','N',[0,0,0]),('CA','C',[.14,.05,0]),('C','C',[.26,0,.02]),('O','O',[.31,-.10,0]),('CB','C',[.14,.12,.13])]:
            topology.addAtom(name, app.element.Element.getBySymbol(symbol), residue)
            coordinates.append(np.asarray(point) + [i*.4,0,0])
    return coordinates


@pytest.mark.parametrize('suffix', ['pdb', 'cif'])
def test_same_author_number_same_residue_names_preserve_insertion_codes(tmp_path, suffix):
    top = app.Topology(); chain = top.addChain('A')
    xyz = protein(top, chain, [(10,''),(10,'A'),(10,'B'),(11,'')])
    path = tmp_path / ('insertions.' + suffix)
    with path.open('w') as handle:
        (app.PDBFile if suffix=='pdb' else app.PDBxFile).writeFile(top, np.asarray(xyz)*unit.nanometer, handle, keepIds=True)
    before = hashlib.sha256(path.read_bytes()).hexdigest()
    meta = sources.import_structure(path)
    assert meta['n_atoms'] == 20
    assert len(meta['atoms']) == 20
    assert storage.load_physical(meta['id']).n_atoms == 20
    native = app.PDBFile(str(preparation.exact_input_path(meta['id'])))
    expected = [('10',''),('10','A'),('10','B'),('11','')]
    assert [(r.id,r.insertionCode.strip()) for r in native.topology.residues()] == expected
    selected = monomers.create_monomer(meta['id'], monomers.MonomerRequest(chain_index=0))
    again = app.PDBFile(str(preparation.exact_input_path(selected['id'])))
    assert [(r.id,r.insertionCode.strip()) for r in again.topology.residues()] == expected
    assert selected['n_atoms'] == 20
    np.testing.assert_allclose(storage.load_physical(selected['id']).xyz[0], xyz, atol=1e-7)
    assert hashlib.sha256(path.read_bytes()).hexdigest() == before


def cif_namespaces(tmp_path, *, symmetry=None):
    top = app.Topology(); chain = top.addChain('A')
    xyz = protein(top, chain, [(101,''),(103,'')])
    residue = top.addResidue('LIG', top.addChain('B'), '501')
    top.addAtom('C1', app.element.carbon, residue); xyz.append([.2,.3,0])
    out = io.StringIO(); app.PDBxFile.writeFile(top,np.asarray(xyz)*unit.nanometer,out,keepIds=True)
    document = gemmi.cif.read_string(out.getvalue()); block = document.sole_block()
    atoms = block.get_mmcif_category('_atom_site.')
    atoms['auth_asym_id'] = ['X'] * len(atoms['id'])
    atoms['label_seq_id'] = [('1' if n=='101' else '3') if c=='A' else False for n,c in zip(atoms['auth_seq_id'],atoms['label_asym_id'])]
    atoms['label_entity_id'] = ['1' if c=='A' else '2' for c in atoms['label_asym_id']]
    block.set_mmcif_category('_atom_site.',atoms)
    block.set_mmcif_category('_struct_asym.',{'id':['A','B'],'entity_id':['1','2']})
    block.set_mmcif_category('_entity_poly.',{'entity_id':['1'],'type':['polypeptide(L)'],'pdbx_strand_id':['X']})
    block.set_mmcif_category('_entity_poly_seq.',{'entity_id':['1']*3,'num':['1','2','3'],'mon_id':['ALA']*3})
    block.set_mmcif_category('_pdbx_poly_seq_scheme.',{'asym_id':['A']*3,'entity_id':['1']*3,'seq_id':['1','2','3'],'mon_id':['ALA']*3,'pdb_seq_num':['101','102','103'],'auth_seq_num':['101',None,'103'],'pdb_strand_id':['X']*3,'pdb_ins_code':[False]*3})
    if symmetry:
        block.set_mmcif_category('_struct_conn.',{
            'id':['link1'],'conn_type_id':['covale'],
            'ptnr1_label_asym_id':['A'],'ptnr1_label_seq_id':['1'],'ptnr1_label_comp_id':['ALA'],'ptnr1_label_atom_id':['CA'],
            'ptnr1_auth_asym_id':['X'],'ptnr1_auth_seq_id':['101'],'pdbx_ptnr1_PDB_ins_code':[None],'ptnr1_symmetry':['1_555'],
            'ptnr2_label_asym_id':['B'],'ptnr2_label_seq_id':[False],'ptnr2_label_comp_id':['LIG'],'ptnr2_label_atom_id':['C1'],
            'ptnr2_auth_asym_id':['X'],'ptnr2_auth_seq_id':['501'],'pdbx_ptnr2_PDB_ins_code':[None],'ptnr2_symmetry':[symmetry],
        })
    path = tmp_path/'namespaces.cif'; document.write_file(str(path)); return path


def test_sequence_mapping_follows_native_label_namespace_not_shared_author(tmp_path):
    path = cif_namespaces(tmp_path)
    meta = sources.import_structure(path)
    selected = monomers.create_monomer(meta['id'],monomers.MonomerRequest(chain_index=0))
    fixer,_ = preparation.current_fixer(selected['id'])
    preparation.find_missing_residues_preserving_identity(fixer)
    assert fixer.missingResidues == {(0,1):['ALA']}
    assert fixer.missingResidueIdentities[(0,1)][0]['resid'] == '102'
    assert 'A' in fixer.sequence_scheme
    assert any(r['name']=='LIG' for r in selected['monomer_selection']['retained_associated_molecules'])


@pytest.mark.parametrize('symmetry,expected_bond', [('1_555',True),('2_555',False),('?',True)])
def test_cif_symmetry_filter_preserves_source_and_same_or_unknown_image_links(tmp_path,symmetry,expected_bond):
    path = cif_namespaces(tmp_path,symmetry=symmetry); original=path.read_bytes()
    meta=sources.import_structure(path)
    native=app.PDBFile(str(preparation.exact_input_path(meta['id'])))
    cross = [(a,b) for a,b in native.topology.bonds() if a.residue.chain != b.residue.chain]
    assert bool(cross) == expected_bond
    assert path.read_bytes()==original
    assert (storage.dataset_dir(meta['id'])/'source.cif').read_bytes()==original
    selected=monomers.create_monomer(meta['id'],monomers.MonomerRequest(chain_index=0))
    assert any(r['name']=='LIG' for r in selected['monomer_selection']['retained_associated_molecules'])
    if expected_bond:
        assert selected['monomer_selection']['retained_covalent_connections']
    else:
        assert selected['monomer_selection']['excluded_symmetry_connections']
        assert any('symmetry mates' in text for text in meta['warnings'])
