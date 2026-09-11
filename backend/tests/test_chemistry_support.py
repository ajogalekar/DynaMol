"""Failures found by the representative release chemistry matrix."""
from pathlib import Path
import json
import numpy as np
import pytest
from openmm import app, unit
from rdkit import Chem
from rdkit.Chem import AllChem

from backend import config, preparation, sources, storage
from backend.ions import inspect_ions, ION_STATES
from backend.ligands import (_graph, _select_state, _restore_protonation_maps,
                             _parameterize, _mapped_variant, _graph_from_override,
                             inspect_ligands, ligand_runtime_status)
from backend.preparation_worker import protonation_inventory
from backend.prepared_system import load_prepared_forcefield


@pytest.fixture
def isolated(tmp_path, monkeypatch):
    monkeypatch.setattr(config, 'DATA_ROOT', tmp_path)
    monkeypatch.setattr(config, 'DATASETS_DIR', tmp_path/'datasets')
    monkeypatch.setattr(config, 'JOBS_DIR', tmp_path/'jobs')
    config.DATASETS_DIR.mkdir(); config.JOBS_DIR.mkdir()
    return tmp_path


@pytest.mark.parametrize('name', sorted(ION_STATES))
def test_declared_ions_match_actual_explicit_forcefield_charge(name):
    import openmm as mm
    symbol, charge = ION_STATES[name]
    top=app.Topology(); r=top.addResidue(name,top.addChain('I'),'1')
    top.addAtom(name,app.element.Element.getBySymbol(symbol),r)
    row=inspect_ions(top)[0]
    assert row['supported'] and row['formal_charge']==charge
    system=app.ForceField('amber14/tip3p.xml').createSystem(top)
    force=next(f for f in system.getForces() if isinstance(f,mm.NonbondedForce))
    assert force.getParticleParameters(0)[0].value_in_unit(unit.elementary_charge)==charge


@pytest.mark.parametrize('name,elements', [('NA',['Mg']),('MG',['C','O']),('FE',['Fe']),('ZN',['Zn'])])
def test_ion_identity_or_unregistered_oxidation_state_is_rejected(name,elements):
    top=app.Topology(); r=top.addResidue(name,top.addChain('I'),'1')
    for i,e in enumerate(elements): top.addAtom(e+str(i),app.element.Element.getBySymbol(e),r)
    row=inspect_ions(top)[0]
    assert not row['supported'] and row['error'] and row['formal_charge'] is None


def test_ion_provenance_cannot_be_loaded_with_implicit_solvent(tmp_path):
    with pytest.raises(ValueError,match='ions require explicit'):
        load_prepared_forcefield(tmp_path,{'ions':[{'element':'Mg','formal_charge':2}]},solvent='implicit')


def test_explicit_deposited_ion_charge_is_not_overwritten():
    top=app.Topology();r=top.addResidue('MG',top.addChain('I'),'1')
    top.addAtom('MG',app.element.magnesium,r,formalCharge=1)
    row=inspect_ions(top)[0]
    assert not row['supported'] and row['deposited_formal_charge']==1
    assert 'no charge was changed' in row['error']


def test_thiolate_and_disulfide_are_different_inventory_states():
    top=app.Topology(); chain=top.addChain('A')
    sulfur=[]
    for i in range(4):
        r=top.addResidue('CYS',chain,str(i)); sulfur.append(top.addAtom('SG',app.element.sulfur,r))
        if i==0:
            hydrogen=top.addAtom('HG',app.element.hydrogen,r); top.addBond(sulfur[0],hydrogen)
    top.addBond(sulfur[2],sulfur[3])
    assert [r['state'] for r in protonation_inventory(top,[])]==['CYS','CYM','CYX','CYX']


def test_standard_protonation_aliases_only_normalize_same_heavy_graph():
    top=app.Topology(); chain=top.addChain('A')
    for i,name in enumerate(['LYN','CYM','MSE','ALY','TPO']):
        r=top.addResidue(name,chain,str(i)); top.addAtom('CA',app.element.carbon,r)
    records=preparation.canonicalize_protonation_aliases(top)
    assert [r.name for r in top.residues()]==['LYS','CYS','MSE','ALY','TPO']
    assert [r['input_alias'] for r in records]==['LYN','CYM']


def test_readiness_validation_uses_submission_rules_without_job_or_dataset_writes(isolated,monkeypatch):
    from backend import config as cfg
    source=cfg.ROOT/'docs/audit/preparation-fixtures/six_residues_intact.pdb'
    meta=sources.import_structure(source)
    def forbidden(*args,**kwargs): raise AssertionError('Readiness submitted a job')
    monkeypatch.setattr(preparation,'_submit',forbidden)
    before={str(p.relative_to(isolated)):p.read_bytes() for p in isolated.rglob('*') if p.is_file()}
    settings,inspection=preparation.validate_preparation({'dataset_id':meta['id']})
    after={str(p.relative_to(isolated)):p.read_bytes() for p in isolated.rglob('*') if p.is_file()}
    assert settings['dataset_id']==meta['id'] and inspection['can_prepare']
    assert before==after and list(cfg.JOBS_DIR.iterdir())==[]


def test_unregistered_protein_crosslink_is_an_upfront_blocker(isolated,monkeypatch):
    source=config.ROOT/'docs/audit/preparation-fixtures/six_residues_intact.pdb'
    meta=sources.import_structure(source)
    original=preparation.current_fixer
    def with_link(*args,**kwargs):
        fixer,path=original(*args,**kwargs)
        atoms=list(fixer.topology.atoms())
        sulfur=next(a for a in atoms if a.name=='SD')
        nitrogen=next(a for a in atoms if a.name=='NZ')
        fixer.topology.addBond(sulfur,nitrogen)
        return fixer,path
    monkeypatch.setattr(preparation,'current_fixer',with_link)
    inspection=preparation.inspect_preparation(meta['id'])
    assert any('Unsupported covalent protein crosslink' in b for b in inspection['blockers'])
    with pytest.raises(ValueError,match='specialized template'):
        preparation.validate_preparation({'dataset_id':meta['id']})
    assert list(config.JOBS_DIR.iterdir())==[]


def test_built_in_caps_are_not_misclassified_as_unsupported_uaa(isolated):
    source=config.ROOT/'docs/audit/modified-residues/native/TPO.pdb'
    meta=sources.import_structure(source)
    settings,inspection=preparation.validate_preparation({'dataset_id':meta['id']})
    assert inspection['can_prepare'] and not inspection['ligands']
    assert {r['residue'] for r in inspection['modified_residues']}=={'TPO'}


def test_missing_protonation_map_recovered_only_with_unique_chemical_identity():
    original=Chem.MolFromSmiles('[CH3:1][O:2][P:3](=[O:4])([O-:5])[OH:6]')
    proposed=Chem.MolFromSmiles('[CH3:1][O:2][P:3](=[O:4])([O-])[O-:6]')
    fixed=_restore_protonation_maps(original,proposed)
    assert {a.GetAtomMapNum() for a in fixed.GetAtoms()}==set(range(1,7))
    ambiguous=Chem.MolFromSmiles('[CH3:1][O:2][P:3](=[O:4])([O-])[O-]')
    with pytest.raises(ValueError,match='uniquely'):
        _restore_protonation_maps(original,ambiguous)
    conflict=Chem.MolFromSmiles('[CH3:1][O:2][P:3](=[O:4])([O-:5])[O-:5]')
    with pytest.raises(ValueError,match='conflicting'):
        _restore_protonation_maps(original,conflict)


def test_radical_input_and_override_have_explicit_unsupported_state_messages(isolated):
    meta=sources.create_smiles('[CH3]')
    assert 'closed-shell' in inspect_ligands(meta['id'])[0]['error']
    with pytest.raises(ValueError,match='closed-shell'):
        _mapped_variant(Chem.MolFromSmiles('C'),'[CH3]')
    top=app.Topology();r=top.addResidue('LIG',top.addChain('A'),'1');top.addAtom('C1',app.element.carbon,r)
    with pytest.raises(ValueError,match='closed-shell'):
        _graph_from_override(r,np.zeros((1,3)),'[CH3]')


def test_nad_phosphate_map_and_stereochemistry_survive_actual_ph_heuristic():
    import gemmi
    path=config.ROOT/'docs/audit/release-readiness/chemistry-inputs/NAD.cif'
    block=gemmi.cif.read_file(str(path)).sole_block(); table=block.get_mmcif_category('_chem_comp_atom.')
    top=app.Topology(); residue=top.addResidue('NAD',top.addChain('A'),'1'); xyz=[]
    for i,name in enumerate(table['atom_id']):
        if table['type_symbol'][i]=='H': continue
        top.addAtom(name,app.element.Element.getBySymbol(table['type_symbol'][i].title()),residue)
        xyz.append([float(table[f'pdbx_model_Cartn_{axis}_ideal'][i]) for axis in 'xyz'])
    mol,atoms,names,stereo=_graph(residue,np.array(xyz),block=block)
    selected,report=_select_state(mol,names,'NAD',7,None)
    assert Chem.GetFormalCharge(selected)==-1 and len(stereo)==8
    assert selected.GetNumAtoms()==44
    assert [a.GetChiralTag() for a in mol.GetAtoms()]==[a.GetChiralTag() for a in selected.GetAtoms()]
    np.testing.assert_array_equal(mol.GetConformer().GetPositions(),selected.GetConformer().GetPositions())
    assert 'unique element/connectivity' in report['warnings'][0]


@pytest.mark.skipif(not ligand_runtime_status()['available'],reason='Native AmberTools unavailable')
def test_halogen_aromatic_native_reversed_improper_converts_exactly(tmp_path):
    mol=Chem.AddHs(Chem.MolFromSmiles('Fc1cc(Cl)cc(Br)c1'))
    assert AllChem.EmbedMolecule(mol,randomSeed=2026)==0
    xml,names,charges=_parameterize(mol,tmp_path/'native','DML_halogen',None,None)
    report=json.loads((xml.parent/'conversion-validation.json').read_text())
    assert report['passed'] and len(report['conformations'])==3
    assert max(abs(r['energy_delta_kj_mol']) for r in report['conformations'])<1e-4
    assert max(r['maximum_force_delta_kj_mol_nm'] for r in report['conformations'])<1e-3
