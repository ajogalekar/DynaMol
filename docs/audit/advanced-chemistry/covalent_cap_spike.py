"""Isolated 6OIM capped-adduct feasibility probe; not an application prep path.

Use only the declared human KRAS CYS12–MOV303 bond. The CCD already describes
the saturated bound ligand. Preserve its heavy-atom graph and every observed
adduct heavy coordinate, and derive hydrogens from the connected valence graph.
The native GAFF2/AM1-BCC capped molecule is a preliminary parameterization
experiment, not yet an ff14SB-compatible residue or an accuracy validation.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import sys
import time
from pathlib import Path

import gemmi
import numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem

ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT))


def build_cap(folder):
    structure = gemmi.cif.read_file(str(folder / '6OIM.cif')).sole_block()
    ccd = gemmi.cif.read_file(str(folder / 'MOV.cif')).sole_block()
    rows = structure.get_mmcif_category('_atom_site.')
    observed = {}
    for i in range(len(rows['id'])):
        if rows['pdbx_PDB_model_num'][i] != '1' or rows['auth_asym_id'][i] != 'A':
            continue
        if rows['label_alt_id'][i] not in (False, None, '.', 'A'):
            continue
        key = (rows['auth_seq_id'][i], rows['label_comp_id'][i], rows['label_atom_id'][i])
        observed[key] = np.array([float(rows[axis][i]) for axis in ('Cartn_x', 'Cartn_y', 'Cartn_z')])
    connection = structure.get_mmcif_category('_struct_conn.')
    declared = [i for i, kind in enumerate(connection['conn_type_id']) if kind == 'covale']
    if len(declared) != 1:
        raise ValueError('This isolated probe requires exactly the documented 6OIM bond.')
    i = declared[0]
    ends = [(connection[f'ptnr{side}_auth_seq_id'][i], connection[f'ptnr{side}_label_comp_id'][i],
             connection[f'ptnr{side}_label_atom_id'][i]) for side in (1, 2)]
    if ends != [('12', 'CYS', 'SG'), ('303', 'MOV', 'C25')]:
        raise ValueError('The deposited covalent link differs from the reviewed probe input.')

    mol = Chem.RWMol()
    names, coords, origins = [], [], []
    def atom(name, element, point, origin):
        a = Chem.Atom(element)
        a.SetProp('source_name', name)
        idx = mol.AddAtom(a)
        names.append(name); coords.append(point); origins.append(origin)
        return idx
    def bond(a, b, kind=Chem.BondType.SINGLE):
        mol.AddBond(names.index(a), names.index(b), kind)
    # ACE and NME caps use observed neighboring peptide backbone coordinates as
    # geometric proxies. Cap atoms are explicitly identified and never exported
    # back to the protein. The central residue and ligand heavy atoms are exact.
    for name, symbol, resid, resname, original in [
        ('ACE_CH3','C','11','ALA','CA'), ('ACE_C','C','11','ALA','C'), ('ACE_O','O','11','ALA','O'),
        ('N','N','12','CYS','N'), ('CA','C','12','CYS','CA'), ('C','C','12','CYS','C'),
        ('O','O','12','CYS','O'), ('CB','C','12','CYS','CB'), ('SG','S','12','CYS','SG'),
        ('NME_N','N','13','GLY','N'), ('NME_CH3','C','13','GLY','CA')]:
        key=(resid,resname,original)
        if key not in observed:
            # Obtain the known neighboring standard identity from the source;
            # never substitute a coordinate from another residue or atom.
            matches=[k for k in observed if k[0]==resid and k[2]==original]
            if len(matches)!=1: raise ValueError(f'Ambiguous cap source atom {key}')
            key=matches[0]
        atom(name,symbol,observed[key], {'chain':'A','resid':key[0],'resname':key[1], 'atom':key[2],
                                        'role':'observed-adduct' if resid=='12' else 'temporary-cap'})
    for left,right in [('ACE_CH3','ACE_C'),('ACE_C','N'),('N','CA'),('CA','C'),('CA','CB'),('CB','SG'),('C','NME_N'),('NME_N','NME_CH3')]:
        bond(left,right)
    for left,right in [('ACE_C','ACE_O'),('C','O')]:bond(left,right,Chem.BondType.DOUBLE)
    atoms=ccd.get_mmcif_category('_chem_comp_atom.')
    for i,name in enumerate(atoms['atom_id']):
        symbol=atoms['type_symbol'][i]
        if symbol=='H': continue
        index=atom('MOV_'+name,symbol,observed[('303','MOV',name)],
                   {'chain':'A','resid':'303','resname':'MOV','atom':name,'role':'observed-adduct'})
        mol.GetAtomWithIdx(index).SetFormalCharge(int(atoms['charge'][i]))
    bonds=ccd.get_mmcif_category('_chem_comp_bond.')
    order={'SING':Chem.BondType.SINGLE,'DOUB':Chem.BondType.DOUBLE,'TRIP':Chem.BondType.TRIPLE,'AROM':Chem.BondType.AROMATIC}
    for a,b,o in zip(bonds['atom_id_1'],bonds['atom_id_2'],bonds['value_order']):
        if 'MOV_'+a in names and 'MOV_'+b in names:bond('MOV_'+a,'MOV_'+b,order[o])
    bond('SG','MOV_C25')
    result=mol.GetMol();Chem.SanitizeMol(result)
    conformer=Chem.Conformer(len(coords))
    for i,xyz in enumerate(coords):conformer.SetAtomPosition(i,xyz)
    result.AddConformer(conformer)
    Chem.AssignStereochemistryFrom3D(result)
    Chem.AssignStereochemistry(result, cleanIt=True, force=True)
    expected={'CA':'R','MOV_C20':'S'}
    for name,cip in expected.items():
        a=result.GetAtomWithIdx(names.index(name))
        if not a.HasProp('_CIPCode') or a.GetProp('_CIPCode')!=cip:
            raise ValueError(f'Observed stereochemistry disagrees at {name}.')
    result=Chem.AddHs(result,addCoords=True)
    for i in range(len(names),result.GetNumAtoms()):
        parent=result.GetAtomWithIdx(i).GetNeighbors()[0].GetIdx()
        names.append('H'+str(i+1));origins.append({'role':'generated-hydrogen','parent_index':parent})
    # Relax only generated H coordinates for charge calculation input.
    props=AllChem.MMFFGetMoleculeProperties(result)
    ff=AllChem.MMFFGetMoleculeForceField(result,props)
    for i in range(len(coords)):ff.AddFixedPoint(i)
    ff.Minimize(maxIts=500)
    if not np.array_equal(np.asarray(result.GetConformer().GetPositions())[:len(coords)], np.asarray(coords)):
        raise ValueError('Hydrogen placement moved an observed heavy atom.')
    hydrogens={name:sum(n.GetAtomicNum()==1 for n in result.GetAtomWithIdx(names.index(name)).GetNeighbors())
               for name in ('SG','MOV_C24','MOV_C25')}
    if hydrogens!={'SG':0,'MOV_C24':2,'MOV_C25':2}:
        raise ValueError(f'Unexpected thioether product H inventory: {hydrogens}')
    report={'source':'6OIM','source_sha256':hashlib.sha256((folder/'6OIM.cif').read_bytes()).hexdigest(),
            'ccd':'MOV','ccd_sha256':hashlib.sha256((folder/'MOV.cif').read_bytes()).hexdigest(),
            'declared_link':ends,'link_distance_angstrom':float(np.linalg.norm(observed[ends[0]]-observed[ends[1]])),
            'reaction_state':'Deposited saturated bound MOV graph; no heavy bond-order change applied.',
            'formal_charge':Chem.GetFormalCharge(result),'atoms':result.GetNumAtoms(),'heavy_atoms':len(coords),
            'bound_heavy_coordinates_preserved_exactly':True,'product_hydrogens':hydrogens,
            'stereochemistry':expected,'atom_map':[dict(index=i,name=name,**origin) for i,(name,origin) in enumerate(zip(names,origins))],
            'scope':'Capped adduct feasibility only; not a prepared full protein or validated ff14SB interface.'}
    return result,report


def main():
    parser=argparse.ArgumentParser();parser.add_argument('--folder',type=Path,default=ROOT/'build/advanced-chemistry/covalent');parser.add_argument('--parameterize',action='store_true')
    args=parser.parse_args();folder=args.folder.resolve();started=time.monotonic()
    mol,report=build_cap(folder)
    (folder/'capped-adduct.sdf').write_text(Chem.MolToMolBlock(mol)+'\n$$$$\n')
    (folder/'capped-adduct.json').write_text(json.dumps(report,indent=2)+'\n')
    if args.parameterize:
        from backend.ligands import _parameterize,_amber_versions
        param=folder/'capped-native'
        xml,names,charges=_parameterize(mol,param,'DML_COV_SPIKE',lambda text:print(text,flush=True),lambda:None)
        report.update(native_parameters={'path':str(param),'xml':str(xml),'charge_sum_e':sum(charges),
                                        'ambertools':_amber_versions(),'validation':json.loads((param/'conversion-validation.json').read_text())})
    report['elapsed_seconds']=time.monotonic()-started
    (folder/'capped-adduct.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps({k:v for k,v in report.items() if k not in ('atom_map','native_parameters')},indent=2))


if __name__=='__main__':main()
