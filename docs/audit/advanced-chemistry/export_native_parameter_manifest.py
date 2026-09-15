"""Export native Amber parameter evidence independently of an OpenMM System.

Only the identity binding uses AmberPrmtopFile.topology; every parameter value
comes from the independently generated native prmtop read through ParmEd.
Supports ordinary Amber 12-6 terms, not Chamber, polarizable, or CMAP models.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import math
from pathlib import Path

import parmed
from openmm import app


def reject_unsupported_native_terms(raw):
    """Do not let an ordinary-term manifest silently omit a native CMAP.

    Some writers include an explicitly empty count. That does not add an
    energy term; inconsistent counts or populated orphan map arrays do.
    Both Amber and CHARMM spellings must be inspected before any export.
    """
    if any(key in raw for key in ['CHARMM_UREY_BRADLEY_COUNT', 'AMOEBA_FORCEFIELD']):
        raise ValueError('This independent exporter supports ordinary Amber 12-6 models only.')
    for prefix in ('', 'CHARMM_'):
        key = prefix+'CMAP_COUNT'
        count = raw.get(key)
        related = [name for name in raw if name.startswith(prefix+'CMAP_') and name != key]
        populated = any(len(raw[name]) for name in related)
        if count is None:
            if populated:
                raise ValueError(f'Unsupported or malformed native CMAP arrays without {key}.')
            continue
        if (not isinstance(count, (list, tuple)) or len(count) != 2 or
                any(isinstance(value, bool) or not isinstance(value, int) or value < 0 for value in count)):
            raise ValueError(f'Malformed native {key}; cannot establish complete parameter coverage.')
        if count != [0, 0] and count != (0, 0):
            raise ValueError(f'Native {key} contains unsupported CMAP terms; no partial manifest was exported.')
        if populated:
            raise ValueError(f'Native {key} is zero but populated CMAP arrays remain; no partial manifest was exported.')


def export(prmtop,output,source_map=None):
    prmtop=Path(prmtop).resolve();output=Path(output).resolve()
    native=parmed.load_file(str(prmtop));raw=native.parm_data
    reject_unsupported_native_terms(raw)
    topology=app.AmberPrmtopFile(str(prmtop)).topology
    topology_atoms=list(topology.atoms())
    if len(topology_atoms)!=len(native.atoms):raise ValueError('Native and imported topology inventories differ.')
    original=json.loads(Path(source_map).read_text()) if source_map else None
    if original and len(original['atom_map'])!=len(native.atoms):raise ValueError('Source mapping length differs.')
    atom_ids=[f'native:{i+1}' for i in range(len(native.atoms))]
    report={'schema_version':1,'scope':'full_system','model_id':'isolated-connected-capped-adduct-'+hashlib.sha256(prmtop.read_bytes()).hexdigest()[:16],
        'native_reference':{'path':str(prmtop),'sha256':hashlib.sha256(prmtop.read_bytes()).hexdigest(),'reader':'ParmEd '+parmed.__version__},
        'chemical_state':{'status':'parameter-interface research specimen; not full prepared protein or validated chemistry',
                          'source_mapping':str(Path(source_map).resolve()) if source_map else None},
        'atoms':[],'bonds':[],'angles':[],'torsions':[],'exceptions':[],
        'units_and_conventions':{'bond_k':'native kcal/mol/angstrom^2 multiplied by 2*4.184*100; harmonic half-k convention',
          'angle_k':'native kcal/mol/radian^2 multiplied by 2*4.184',
          'torsion_k':'native kcal/mol amplitude multiplied by 4.184; preserve ordered quartet and improper marker',
          'nonbonded':'native charge_e, native stored A/B 12-6 coefficients and per-dihedral SCEE/SCNB; all zero exclusion pairs retained'}}
    for atom,target in zip(native.atoms,topology_atoms):
        # AmberPrmtopFile explicitly normalizes its PDB naming aliases (e.g.
        # native N-terminal H1 -> H). Apply that exact table, preserving the
        # original native labels below; never permit an arbitrary atom reorder.
        imported_residue_name=app.PDBFile._residueNameReplacements.get(atom.residue.name,atom.residue.name)
        imported_atom_name=app.PDBFile._atomNameReplacements.get(imported_residue_name,{}).get(atom.name,atom.name)
        if imported_residue_name!=target.residue.name or imported_atom_name!=target.name or atom.atomic_number!=target.element.atomic_number:raise ValueError('Native identity order differs.')
        row={'id':atom_ids[atom.idx],'identity':{'chain':str(target.residue.chain.id),'resid':str(target.residue.id),
              'insertion':str(target.residue.insertionCode or '').strip(),'resname':str(target.residue.name),'atomname':str(target.name)},
              'element':target.element.symbol,'mass_da':float(atom.mass),'charge_e':float(atom.charge),
              'sigma_nm':float(atom.sigma)*.1,'epsilon_kj_mol':float(atom.epsilon)*4.184,
              'native_atom_type':str(atom.type),'native_index':atom.idx,
              'native_atom_name':str(atom.name),'native_residue_name':str(atom.residue.name)}
        if original:row['source_record']=original['atom_map'][atom.idx]
        report['atoms'].append(row)
    for bond in native.bonds:
        term={'atoms':[atom_ids[x.idx] for x in [bond.atom1,bond.atom2]],
            'length_nm':bond.type.req*.1,'k_kj_mol_nm2':bond.type.k*836.8}
        if (bond.atom1.atomic_number==bond.atom2.atomic_number==1 and
            bond.atom1.residue is bond.atom2.residue and bond.atom1.residue.name=='WAT'):
            # Native Amber stores this physical geometry term, but OpenMM
            # deliberately omits an H-H chemical edge from water topology.
            term['topology_role']='native_water_hh_geometry'
        report['bonds'].append(term)
    for angle in native.angles:
        report['angles'].append({'atoms':[atom_ids[x.idx] for x in [angle.atom1,angle.atom2,angle.atom3]],
            'angle_radian':math.radians(angle.type.theteq),'k_kj_mol_rad2':angle.type.k*8.368})
    for torsion in native.dihedrals:
        report['torsions'].append({'atoms':[atom_ids[x.idx] for x in [torsion.atom1,torsion.atom2,torsion.atom3,torsion.atom4]],
            'periodicity':int(torsion.type.per),'phase_radian':math.radians(torsion.type.phase),
            'k_kj_mol':torsion.type.phi_k*4.184,'improper':bool(torsion.improper),'ignore_end':bool(torsion.ignore_end)})
    # Amber excludes the full bonded-neighbor list from general nonbonded pairs,
    # then reinstates explicitly flagged 1-4 pairs using stored A/B and scaling.
    exceptions={};offset=0
    for i,count in enumerate(raw['NUMBER_EXCLUDED_ATOMS']):
        for one_based_j in raw['EXCLUDED_ATOMS_LIST'][offset:offset+count]:
            if one_based_j>0:
                pair=tuple(sorted((i,one_based_j-1)))
                exceptions[pair]={'atoms':[atom_ids[x] for x in pair],'chargeprod_e2':0.,'sigma_nm':1.,'epsilon_kj_mol':0.,'native_kind':'zero_exclusion'}
        offset+=count
    ntypes=native.ptr('ntypes');scaled={}
    for torsion in native.dihedrals:
        if torsion.improper or torsion.ignore_end:continue
        i,j=torsion.atom1.idx,torsion.atom4.idx;pair=tuple(sorted((i,j)))
        ti,tj=raw['ATOM_TYPE_INDEX'][i]-1,raw['ATOM_TYPE_INDEX'][j]-1
        term=raw['NONBONDED_PARM_INDEX'][ntypes*ti+tj]-1
        if term<0:raise ValueError('Native 10-12 hydrogen-bond pairs require a separate supported exporter.')
        a=raw['LENNARD_JONES_ACOEF'][term];b=raw['LENNARD_JONES_BCOEF'][term]
        sigma=(a/b)**(1/6)*.1 if a and b else 1.
        epsilon=b*b/(4*a)*4.184 if a and b else 0.
        if torsion.type.scee<=0 or torsion.type.scnb<=0:raise ValueError('Invalid native scaling for a physical 1-4 pair.')
        values={'atoms':[atom_ids[x] for x in pair],'chargeprod_e2':native.atoms[i].charge*native.atoms[j].charge/torsion.type.scee,
                'sigma_nm':sigma,'epsilon_kj_mol':epsilon/torsion.type.scnb,'native_kind':'scaled_1_4',
                'native_scee':torsion.type.scee,'native_scnb':torsion.type.scnb}
        if pair in scaled and scaled[pair]!=values:raise ValueError('Conflicting duplicate native 1-4 definitions.')
        scaled[pair]=values
    exceptions.update(scaled);report['exceptions']=[exceptions[pair] for pair in sorted(exceptions)]
    report['counts']={key:len(report[key]) for key in ['atoms','bonds','angles','torsions','exceptions']}
    report['exporter_sha256']=hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
    output.write_text(json.dumps(report,indent=2)+'\n')
    return report


if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('prmtop',type=Path);parser.add_argument('output',type=Path);parser.add_argument('--source-map',type=Path)
    args=parser.parse_args();print(json.dumps(export(args.prmtop,args.output,args.source_map)['counts'],indent=2))
