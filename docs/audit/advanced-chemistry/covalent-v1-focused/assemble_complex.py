"""Assemble native candidate complexes without discarding retained environment.

Research-only preparation, with exact adduct inventory/charge checks. A valid
native topology does not establish physical flexibility or release readiness.
"""
from __future__ import annotations

import argparse
from collections import OrderedDict
import json
from pathlib import Path

import numpy as np
import parmed
from scipy.spatial import cKDTree

from prepare_adduct import digest, mol2_sections, native


def pdb_residues(path):
    residues=OrderedDict()
    for line in path.read_text().splitlines():
        if line.startswith(('ATOM  ','HETATM')):
            key=(line[21:22],line[22:26].strip(),line[26:27].strip(),line[17:20].strip())
            residues.setdefault(key,[]).append(line)
    return residues


def source_identity(record):
    if 'source' in record:
        source=record['source']; identity=source['identity']
        return (identity['author_chain'],identity['author_residue_id'],
                identity.get('insertion_code') or '',record['resname'],source['atom'])
    return (record['chain'],record['resid'],'',record['resname'],record['atom'])


def assemble(root,code,amber,output,exclude_nonadduct_drug_copies=False,protein_path=None):
    output.mkdir(parents=True,exist_ok=False)
    prefix=code.lower()
    baseline=root/('baseline-6oim' if prefix=='6oim' else prefix+'-baseline-v1')
    prep=root/(prefix+'-prepgen-v2')
    protein=protein_path if protein_path is not None else root/(prefix+'-protein-v1')
    if not (prep/'result.json').exists() and (root/(prefix+'-prepgen-v1')/'result.json').exists():
        prep=root/(prefix+'-prepgen-v1')
    if not (protein/'result.json').is_file():
        raise ValueError('Protein intermediate has not passed geometry/stereochemistry checks')
    protein_report=json.loads((protein/'result.json').read_text())
    if protein_report.get('stage') != 'protein_intermediate_complete':
        raise ValueError('Protein result is not an accepted coordinate intermediate')
    if protein_report.get('loop_geometry') is not None and not protein_report['loop_geometry']['accepted']:
        raise ValueError('Protein intermediate failed its loop geometry screen')
    environment_check=protein_report.get('retained_environment_check')
    if environment_check is not None and not environment_check['passed']:
        raise ValueError('Protein intermediate failed its retained-environment screen')
    for name,expected in protein_report.get('outputs_sha256',{}).items():
        if digest(protein/name)!=expected:
            raise ValueError('Screened protein intermediate changed: '+name)
    report=json.loads((prep/'result.json').read_text())
    changes=report['charge_changes']
    target=source_identity(next(x['source'] for x in changes if x['name']=='N'))[:4]
    ligand_records=[x for x in changes if x['source'].get('role')=='observed-adduct'
                    and x['source'].get('resname')!='CYS']
    components={x['source']['resname'] for x in ligand_records}
    if len(components)!=1:
        raise ValueError('Expected exactly one fully mapped drug component')
    component=components.pop()
    gdp_source=root/'gdp-scoped-v6'
    gdp_loaded=False
    raw=mol2_sections(baseline/'interface-native-scope-v4/hybrid.mol2')['ATOM']
    protein_residues=pdb_residues(protein/'protein.pdb')
    if target not in protein_residues:
        raise ValueError('Mapped cysteine is missing from the repaired protein')
    rows=[]; identities=[]; residue_inventory=[]; serial=0; resnum=0
    def add_atom(name,element,xyz,source=None):
        nonlocal serial
        serial+=1
        if len(name)>4: raise ValueError('Atom name exceeds PDB capacity')
        label=f' {name:<3}' if len(element)==1 and len(name)<4 else f'{name:<4}'
        x,y,z=xyz
        rows.append(f'ATOM  {serial:5d} {label} {resname:>3} A{resnum:4d}    '
                    f'{x:8.3f}{y:8.3f}{z:8.3f}  1.00  0.00          {element:>2}  ')
        identities.append({'input_serial':serial,'assembled_residue_index':resnum-1,
                           'name':name,'element':element,'source':source})
    max_cys_shift=0.; protein_chain=None
    for key,lines in protein_residues.items():
        if protein_chain is not None and key[0]!=protein_chain: rows.append('TER')
        protein_chain=key[0];resnum+=1
        resname='COV' if key==target else key[3]
        residue_inventory.append({'index':resnum-1,'name':resname,'source':key})
        if key==target:
            heavy_by_name={line[12:16].strip():np.array([float(line[a:b]) for a,b in [(30,38),(38,46),(46,54)]])
                           for line in lines if line[76:78].strip()!='H'}
            for item in changes:
                i=item['capped_index'];atom=raw[i];xyz=np.array([float(v) for v in atom[2:5]])
                source=item['source'];element=''.join(c for c in atom[1] if c.isalpha())
                if source.get('resname')=='CYS':
                    name=source['atom'];shift=float(np.linalg.norm(heavy_by_name[name]-xyz))
                    max_cys_shift=max(max_cys_shift,shift)
                    if shift>.002:
                        raise ValueError('Repair moved the covalent residue; a recorded pose reconciliation is required')
                add_atom(item['name'],element,xyz,source)
        else:
            for line in lines:
                if line[76:78].strip()=='H':continue
                add_atom(line[12:16].strip(),line[76:78].strip(),
                         [float(line[a:b]) for a,b in [(30,38),(38,46),(46,54)]],
                         {'residue':key,'atom':line[12:16].strip()})
    rows.append('TER')
    environment=pdb_residues(protein/'retained-environment.pdb')
    declared_drug=source_identity(ligand_records[0]['source'])
    ligand_residues=[(k,v) for k,v in environment.items() if k[1:4]==declared_drug[1:4]]
    if len(ligand_residues)!=1:
        raise ValueError('Expected one retained drug residue for mapped replacement')
    expected_drug={x['source']['atom'] for x in ligand_records}
    actual_drug={x[12:16].strip() for x in ligand_residues[0][1] if x[76:78].strip()!='H'}
    if expected_drug!=actual_drug:
        raise ValueError('Not every original drug heavy atom maps to the adduct')
    excluded=[]
    for key,lines in environment.items():
        if key==ligand_residues[0][0]:continue
        if key[3]==component and exclude_nonadduct_drug_copies:
            excluded.append({'identity':key,'atoms':len(lines),
                'reason':'Explicit focused-control selection retains the deposited covalent adduct and excludes other noncovalently associated copies of the same drug. Original input/environment files remain intact.'})
            continue
        resnum+=1
        aliases={'HOH':'WAT','CL':'Cl-','MG':'MG','GDP':'GDP'}
        if key[3] not in aliases:
            raise ValueError('Retained cofactor needs explicit native parameters: '+str(key))
        resname=aliases[key[3]]
        residue_inventory.append({'index':resnum-1,'name':resname,'source':key})
        if resname=='GDP':
            if prefix!='6oim' or not (gdp_source/'result.json').exists():
                raise ValueError('Retained GDP needs its qualified source parameter model')
            gdp_atoms=json.loads((gdp_source/'bound-atoms.json').read_text())
            actual={line[12:16].strip() for line in lines if line[76:78].strip()!='H'}
            if actual!={atom['source_name'] for atom in gdp_atoms if atom['source_heavy']}:
                raise ValueError('Not every observed GDP heavy atom maps to the source model')
            for atom in gdp_atoms:
                add_atom(atom['name'],atom['element'],atom['xyz_angstrom'],
                         {'residue':key,'atom':atom['source_name'],'source_heavy':atom['source_heavy']})
            rows.append('TER');gdp_loaded=True
            continue
        for line in lines:
            if line[76:78].strip()=='H':continue
            name='Cl-' if resname=='Cl-' else line[12:16].strip()
            add_atom(name,line[76:78].strip(),[float(line[a:b]) for a,b in [(30,38),(38,46),(46,54)]],
                     {'residue':key,'atom':line[12:16].strip()})
        rows.append('TER')
    rows.append('END')
    (output/'complex-input.pdb').write_text('\n'.join(rows)+'\n')
    for file in ['covalent.prepi','gaff-adduct.frcmod']:
        (output/file).write_bytes((prep/file).read_bytes())
    extra_setup=''
    if gdp_loaded:
        for name in ['gdp-scoped.prep','gdp-scoped.frcmod','gdp-elements.leap']:
            (output/name).write_bytes((gdp_source/name).read_bytes())
        extra_setup='loadamberparams gdp-scoped.frcmod\nloadamberprep gdp-scoped.prep\nsource gdp-elements.leap\n'
    setup=('source leaprc.protein.ff14SB\nsource leaprc.gaff2\nsource leaprc.water.tip3p\n'
           'loadamberparams gaff-adduct.frcmod\nloadamberprep covalent.prepi\n'+extra_setup+
           'model=loadpdb complex-input.pdb\ncheck model\n')
    (output/'dry.leap').write_text(setup+'saveamberparm model dry.prmtop dry.inpcrd\n'
                                 'savepdb model dry.pdb\nquit\n')
    commands=[native(amber,output,['tleap','-f','dry.leap'])]
    (output/'tleap-dry.log').write_bytes((output/'tleap.log').read_bytes())
    dry=parmed.load_file(str(output/'dry.prmtop'),xyz=str(output/'dry.inpcrd'))
    charge=float(sum(a.charge for a in dry.atoms))
    if abs(charge-round(charge))>1e-4:
        raise ValueError('Full-complex charge is not consistent with an integer state')
    if gdp_loaded:
        from scope_gdp import evaluate
        source_gdp=parmed.load_file(str(gdp_source/'gdp.prmtop'))
        extracted=dry[':GDP'];extracted.save(str(output/'gdp-extracted.prmtop'))
        if [a.name for a in extracted.atoms]!=[a.name for a in source_gdp.atoms]:
            raise ValueError('Full assembly changed GDP atom order')
        positions=np.load(gdp_source/'bound-positions-nm.npy')
        e1,f1=evaluate(gdp_source/'gdp.prmtop',positions)
        e2,f2=evaluate(output/'gdp-extracted.prmtop',positions)
        if abs(e2-e1)>1e-4 or np.max(abs(f1-f2))>1e-3:
            raise ValueError('Protein force-field loading changed GDP mechanics')
        gdp_types={a.type for a in extracted.atoms}
        if any(a.type in gdp_types for a in dry.atoms if a.residue.name!='GDP'):
            raise ValueError('GDP types overlap other molecule parameter namespaces')
        mg=[a for a in dry.atoms if a.residue.name=='MG']
        if len(mg)!=1 or mg[0].atomic_number!=12 or abs(mg[0].charge-2)>1e-6:
            raise ValueError('The declared magnesium ion was not retained as Mg(II)')
    modified=[r for r in dry.residues if r.name=='COV']
    if len(modified)!=1:raise ValueError('Native assembly did not retain one adduct')
    cov=modified[0];by_name={a.name:a for a in cov.atoms}
    if set(by_name)!={x['name'] for x in changes}:raise ValueError('Native assembly changed adduct atoms')
    for item in changes:
        if abs(by_name[item['name']].charge-item['native_prepgen_charge'])>1e-8:
            raise ValueError('Full protein assembly changed an adduct charge')
    peptide=parmed.load_file(str(prep/'peptide.prmtop'))
    ref=next(r for r in peptide.residues if r.name=='COV')
    def bond_inventory(model,residue):
        return {tuple(sorted((b.atom1.name,b.atom2.name))):(b.type.k,b.type.req)
                for b in model.bonds if b.atom1.residue is residue and b.atom2.residue is residue}
    if bond_inventory(dry,cov)!=bond_inventory(peptide,ref):
        raise ValueError('Full assembly changed adduct bond parameters/connectivity')
    boundary=[b for b in dry.bonds if (b.atom1.residue is cov)!=(b.atom2.residue is cov)]
    if len(boundary)!=2 or any(b.type is None for b in boundary):
        raise ValueError('Adduct must retain both peptide connections')
    # Check that every input heavy atom survives native assembly at the same position.
    mapped=[];input_xyz=np.array([[float(line[a:b]) for a,b in [(30,38),(38,46),(46,54)]]
                for line in rows if line.startswith('ATOM')])
    for entry,xyz in zip(identities,input_xyz):
        residue=dry.residues[entry['assembled_residue_index']]
        if residue.name!=residue_inventory[entry['assembled_residue_index']]['name']:
            raise ValueError('Native residue order changed')
        atoms=[a for a in residue.atoms if a.name==entry['name']]
        if len(atoms)!=1:raise ValueError('Native identity mapping is ambiguous')
        atom=atoms[0]
        if np.linalg.norm(dry.coordinates[atom.idx]-xyz)>.0001:
            raise ValueError('Native assembly moved an input atom')
        mapped.append(dict(entry,native_index=atom.idx))
    solvation=setup+'solvatebox model TIP3PBOX 10.0\n'
    replacements=[]
    ions=''
    if round(charge):
        # Native addions can replace a deposited water. Select exclusively
        # newly added bulk solvent here, so the retained source stays intact.
        (output/'water.leap').write_text(solvation+
            'saveamberparm model water.prmtop water.inpcrd\nquit\n')
        commands.append(native(amber,output,['tleap','-f','water.leap'],seconds=120))
        water=parmed.load_file(str(output/'water.prmtop'),xyz=str(output/'water.inpcrd'))
        box=water.box[:3]
        if not np.allclose(water.box[3:],90):raise ValueError('Bulk ion placement currently requires an orthorhombic box')
        retained_heavy=[a.idx for a in water.atoms[:len(dry.atoms)] if a.atomic_number>1]
        tree=cKDTree(water.coordinates[retained_heavy]%box,boxsize=box)
        candidates=[r for r in water.residues[len(dry.residues):] if r.name=='WAT']
        rng=np.random.default_rng(2026)
        for index in rng.permutation(len(candidates)):
            residue=candidates[index];oxygen=next(a for a in residue.atoms if a.atomic_number==8)
            position=water.coordinates[oxygen.idx]
            if tree.query(position%box)[0]<5:continue
            if any(np.linalg.norm((position-np.array(x['position_angstrom'])+box/2)%box-box/2)<5 for x in replacements):continue
            replacements.append({'water_residue_index':residue.idx,'position_angstrom':position.tolist()})
            if len(replacements)==abs(round(charge)):break
        if len(replacements)!=abs(round(charge)):raise ValueError('Insufficient newly added bulk water for retained-environment-safe ion placement')
        for item in sorted(replacements,key=lambda x:x['water_residue_index'],reverse=True):
            ions+=f"remove model model.{item['water_residue_index']+1}\n"
        for i,item in enumerate(replacements):
            position=' '.join(f'{x:.9f}' for x in item['position_angstrom'])
            ions+=f'ion{i}=copy {"Cl-" if charge>0 else "Na+"}\n'
            ions+=f'set ion{i}.1.1 position {{ {position} }}\nmodel=combine {{ model ion{i} }}\n'
        ions+='set model box { '+' '.join(f'{x:.9f}' for x in box)+' }\n'
    (output/'solvate.leap').write_text(solvation+ions+
        'check model\nsaveamberparm model solvated.prmtop solvated.inpcrd\n'
        'savepdb model solvated.pdb\nquit\n')
    commands.append(native(amber,output,['tleap','-f','solvate.leap'],seconds=120))
    solvated=parmed.load_file(str(output/'solvated.prmtop'),xyz=str(output/'solvated.inpcrd'))
    from native_elements import inventory
    inventory(output/'dry.prmtop')
    inventory(output/'solvated.prmtop')
    if abs(sum(a.charge for a in solvated.atoms))>1e-4:
        raise ValueError('Solvated model was not neutralized')
    # LEaP centers the box and moves counterions before water residues. Match
    # the retained coordinates and names instead of assuming stable indices.
    translation=solvated.coordinates[0]-dry.coordinates[0]
    distances,indices=cKDTree(solvated.coordinates).query(dry.coordinates+translation)
    if np.max(distances)>.0001 or len(set(indices))!=len(dry.atoms):
        raise ValueError('Solvation removed or moved a retained atom beyond a common translation')
    for old,new in zip(dry.atoms,indices):
        atom=solvated.atoms[int(new)]
        if (old.name,old.residue.name)!=(atom.name,atom.residue.name):
            raise ValueError('Solvation changed a retained atom identity')
    for item in mapped:
        item['dry_native_index']=item['native_index']
        item['native_index']=int(indices[item['native_index']])
    (output/'source-mapping.json').write_text(json.dumps(mapped,indent=2)+'\n')
    report={'stage':'full_complex_native_assembly_complete','case':code,
        'physical_model_validated':False,'app_ready':False,'dynamics_completed':False,
        'native_commands':commands,'protein_atoms_in_intermediate':protein_report['atoms'],
        'protein_intermediate':str(protein),
        'dry_atoms':len(dry.atoms),'solvated_atoms':len(solvated.atoms),
        'dry_charge_e':charge,'covalent_residue_index':cov.idx,'covalent_atom_indices':[a.idx for a in cov.atoms],
        'modeled_residues':protein_report['modeled_residues'],
        'all_input_atoms_mapped':True,'native_elements_explicitly_verified':True,
        'retained_environment':[(k,len(v)) for k,v in environment.items() if not any(tuple(x['identity'])==k for x in excluded)],
        'deposited_environment':[(k,len(v)) for k,v in environment.items()],
        'input_scope':'Explicitly selected complex; any excluded deposited molecules are listed separately',
        'excluded_environment':excluded,'all_deposited_environment_retained':not excluded,
        'solvation_translation_angstrom':translation.tolist(),
        'retained_native_atoms_after_solvation':len(indices),
        'counterion_placement':{'method':'Replace newly added bulk waters only, at least 5 angstrom from retained heavy atoms and each other under periodic distances',
                               'seed':2026,'replacements':replacements},
        'maximum_cysteine_heavy_shift_angstrom':max_cys_shift,
        'solvation':'TIP3P rectangular box, 10 angstrom clearance; counterions only, no additional salt',
        'source_sha256':{str(p):digest(p) for p in [protein/'protein.pdb',protein/'retained-environment.pdb',
                    protein/'result.json',prep/'covalent.prepi',prep/'gaff-adduct.frcmod']},
        'outputs_sha256':{n:digest(output/n) for n in ['dry.prmtop','dry.inpcrd','solvated.prmtop','solvated.inpcrd','source-mapping.json']}}
    (output/'result.json').write_text(json.dumps(report,indent=2)+'\n')
    return {k:v for k,v in report.items() if k not in ['native_commands','retained_environment','deposited_environment','covalent_atom_indices','source_sha256','outputs_sha256']}


if __name__=='__main__':
    parser=argparse.ArgumentParser()
    for n in ['root','amber','output']:parser.add_argument('--'+n,type=Path,required=True)
    parser.add_argument('--case',required=True)
    parser.add_argument('--exclude-nonadduct-drug-copies',action='store_true')
    parser.add_argument('--protein',type=Path)
    args=parser.parse_args()
    print(json.dumps(assemble(args.root,args.case,args.amber,args.output,args.exclude_nonadduct_drug_copies,args.protein),indent=2))
