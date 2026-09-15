"""Isolated complete 1HCK assembly with the published Hu2024 ATP/Mg model.

No MD, no QM, no application integration. The final SI files are immutable and
have a CC BY-NC license. CMAP and pair-specific LJ are required model terms.
"""
from __future__ import annotations
import argparse
import hashlib
import json
import os
from pathlib import Path
import subprocess
import sys

import gemmi
import numpy as np
import parmed

BASE = Path(__file__).resolve().parent
ROOT = BASE.parents[2]
PARAM = BASE/'prototype-1hck-mg-atp/published/Parameters'
PRIOR = ROOT/'docs/audit/loop-fallback/runs/review-final/1HCK/workspace/jobs/6677c28a341e45f6'
SOURCE = BASE/'inputs/1HCK.cif'
SOURCE_SHA = '62807748b242e7e8026668edbc4242d1e46de8770e784bcfdf97b73219976d5a'


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def dump(path, obj):
    path.write_text(json.dumps(obj, indent=2)+'\n')


def pdb_record(serial, name, resname, residue, xyz, element, chain='A'):
    atom_name = f' {name:<3}' if len(element)==1 and len(name)<4 else f'{name:<4}'
    return (f'ATOM  {serial:5d} {atom_name} {resname:>3} {chain}{residue:4d}    '
            f'{xyz[0]:8.3f}{xyz[1]:8.3f}{xyz[2]:8.3f}  1.00  0.00          {element:>2}  ')


def main():
    parser=argparse.ArgumentParser();parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args();out=args.output.resolve();out.mkdir(exist_ok=False)
    if sha(SOURCE)!=SOURCE_SHA or sha(PRIOR/'sequence-source.cif')!=SOURCE_SHA:
        raise ValueError('Frozen source or prior loop-repair source changed')
    block=gemmi.cif.read_file(str(SOURCE)).sole_block()
    raw=block.get_mmcif_category('_atom_site.')
    source={};excluded=[]
    for i in range(len(raw['id'])):
        record={k:v[i] for k,v in raw.items()}
        if record['pdbx_PDB_model_num']!='1':raise ValueError('Unexpected multiple models')
        key=(record['label_asym_id'],record['auth_seq_id'],record['label_comp_id'],record['label_atom_id'])
        if record['label_alt_id'] not in (False,None,'A'):
            if key[:3]!=('A','131','GLN') or record['label_alt_id']!='B':
                raise ValueError('Unreviewed alternate conformer')
            excluded.append(record);continue
        if key in source:raise ValueError('Duplicate selected source identity')
        source[key]=record
    if len(source)!=2506 or len(excluded)!=4:raise ValueError('Frozen source inventory changed')
    # Reuse a previously accepted short-loop model, with exact observed-heavy
    # comparisons. This is a provisional conformation, not native-loop truth.
    residues={};histidines=[]
    for line in (PRIOR/'prepared.pdb').read_text().splitlines():
        if not line.startswith(('ATOM  ','HETATM')) or line[21]!='A':continue
        rid=line[22:26].strip();name=line[12:16].strip();resname=line[17:20].strip()
        xyz=np.array([float(line[30:38]),float(line[38:46]),float(line[46:54])])
        element=line[76:78].strip().upper()
        residues.setdefault(rid,[]).append(dict(name=name,resname=resname,xyz=xyz,element=element))
    if list(residues)!=[str(x) for x in range(1,299)]:raise ValueError('Protein residue map differs')
    lines=[];map_input={};serial=0;max_delta=0.;matched=set()
    def emit(native_resid, name, resname, xyz, element, record):
        nonlocal serial
        key=(native_resid,name)
        if key in map_input:raise ValueError('Ambiguous native input identity')
        serial+=1;lines.append(pdb_record(serial,name,resname,native_resid,xyz,element))
        map_input[key]=dict(record,input_name=name,input_residue_name=resname,
                           input_xyz_angstrom=np.asarray(xyz).tolist())
    for rid,atoms in residues.items():
        names={a['name'] for a in atoms};resname=atoms[0]['resname']
        native_name=resname
        if resname=='HIS':
            native_name='HIP' if {'HD1','HE2'}<=names else 'HID' if 'HD1' in names else 'HIE' if 'HE2' in names else None
            if native_name is None:raise ValueError('Unassigned histidine state')
            histidines.append({'author_chain':'A','author_resid':rid,'variant':native_name,
                              'selection':'reused pH7 geometric hydrogen template; not a pKa calculation'})
        for a in atoms:
            key=('A',rid,resname,a['name']);record=source.get(key)
            xyz=a['xyz'];name=a['name']
            if record is not None:
                original=np.array([float(record[f'Cartn_{c}']) for c in 'xyz'])
                delta=float(np.max(np.abs(xyz-original)));max_delta=max(max_delta,delta)
                if delta>.0005001:raise ValueError(f'Observed source heavy atom moved: {key} {delta}')
                xyz=original;matched.add(key)
                role='observed source heavy atom'
            else:
                role='added hydrogen' if a['element']=='H' else 'modeled loop heavy atom' if 37<=int(rid)<=40 else 'added terminal heavy atom'
                if a['element']!='H' and not (37<=int(rid)<=40 or (rid=='298' and name=='OXT')):
                    raise ValueError(f'Unreviewed new protein heavy atom: {key}')
            if rid=='1' and name=='H' and {'H2','H3'}<=names:name='H1'
            emit(int(rid),name,native_name,xyz,a['element'],{
                'source_or_added_id':':'.join(key),'role':role,'source_atom_site':record})
    lines.append('TER');protein_lines=list(lines);protein_atoms=serial
    # Restore every deposited heterogen/water from original source, including
    # all 108 water oxygens. Native LEaP will add template hydrogens only.
    nonprotein={}
    for key,record in source.items():
        if key[0]!='A':nonprotein.setdefault(key[:3],[]).append((key,record))
    if sorted((k[0],k[2],len(v)) for k,v in nonprotein.items()).count(('B','MG',1))!=1:
        raise ValueError('Expected exactly one Mg')
    heterogen_residues={}
    for native_resid,(reskey,atoms) in enumerate(nonprotein.items(),start=299):
        heterogen_residues[native_resid]=reskey
        for key,record in atoms:
            xyz=np.array([float(record[f'Cartn_{c}']) for c in 'xyz'])
            emit(native_resid,key[3],key[2],xyz,record['type_symbol'],{
                'source_or_added_id':':'.join(key),'role':'observed source heavy atom','source_atom_site':record})
            matched.add(key)
        lines.append('TER')
    if matched!=set(source):raise ValueError('Not all selected original heavy atoms retained')
    (out/'prepared-input.pdb').write_text('\n'.join(lines+['END'])+'\n')
    (out/'protein-input.pdb').write_text('\n'.join(protein_lines+['END'])+'\n')
    dump(out/'input-selection.json',{'source_sha256':SOURCE_SHA,'source_atom_rows':2510,
         'selected_source_heavy_atoms':2506,'alternate_selection':'GLN131 A occupancy0.80; all four B occupancy0.20 atom records retained below',
         'excluded_alternate_atom_records':excluded,'all_crystal_water_oxygens_retained':108,
         'Mg_occupancy':0.60,'occupancy_interpretation':'explicit occupied microstate; full Mg(II) +2 model, no occupancy-scaled charge',
         'loop_repair':{'residues':'A37–40 LDTE','source':'previous independently accepted preparation',
                        'prepared_pdb_sha256':sha(PRIOR/'prepared.pdb'),'provenance_sha256':sha(PRIOR/'loop-model-provenance.json'),
                        'maximum_observed_heavy_difference_angstrom':max_delta},
         'histidines':histidines,'source_struct_conn':block.get_mmcif_category('_struct_conn.'),
         'interpretation':'Explicit template state and candidate loop geometry; not experimental protonation or loop accuracy.'})
    amber=ROOT/'.tools/ambertools'
    leap=['source leaprc.protein.ff14SB','source leaprc.water.tip3p',
          'protein = loadpdb protein-input.pdb','saveamberparm protein protein-reference.prmtop protein-reference.inpcrd',
          'addAtomTypes { { "OY" "O" "sp3" } { "O3" "O" "sp2" } }',
          f'loadAmberPrep {PARAM/"ATP-B3.prepi"}',f'loadAmberParams {PARAM/"ATP-B3.frcmod"}',
          'mol = loadpdb prepared-input.pdb','check mol','savepdb mol baseline-native.pdb',
          'saveamberparm mol baseline.prmtop baseline.inpcrd','quit']
    (out/'tleap.in').write_text('\n'.join(leap)+'\n')
    env=dict(os.environ,AMBERHOME=str(amber),OMP_NUM_THREADS='2',OPENBLAS_NUM_THREADS='2')
    with (out/'tleap.log').open('w') as log:
        run=subprocess.run([str(amber/'bin/tleap'),'-s','-f','tleap.in'],cwd=out,env=env,stdout=log,stderr=subprocess.STDOUT,timeout=180)
    if run.returncode or not (out/'baseline.prmtop').exists():raise ValueError('Native LEaP failed; retained log')
    native=parmed.load_file(str(out/'baseline.prmtop'),xyz=str(out/'baseline.inpcrd'))
    atom_map=[];matched_native=set();max_native_delta=0.
    for atom in native.atoms:
        key=(atom.residue.idx+1,atom.name);record=map_input.get(key)
        if record is None:
            if atom.atomic_number!=1 or atom.residue.idx+1 not in heterogen_residues:
                raise ValueError(f'Unexpected newly generated native atom {key}')
            reskey=heterogen_residues[atom.residue.idx+1]
            record={'role':'native template added hydrogen','source_or_added_id':':'.join(reskey+(atom.name,)),
                    'source_atom_site':None,'input_xyz_angstrom':None}
        else:
            delta=float(np.max(np.abs(native.coordinates[atom.idx]-record['input_xyz_angstrom'])))
            max_native_delta=max(max_native_delta,delta)
            if delta>.0005001:raise ValueError('Native LEaP moved input atoms')
            if record['role']=='observed source heavy atom':matched_native.add(record['source_or_added_id'])
        atom_map.append(dict(record,native_index=atom.idx,native_atom_name=atom.name,
                             native_residue_name=atom.residue.name,native_residue_index=atom.residue.idx))
    if len(matched_native)!=2506:raise ValueError('Native source retention failure')
    dump(out/'source-atom-map.json',{'atom_map':atom_map})
    with (out/'published-modification.log').open('w') as log:
        run=subprocess.run([sys.executable,str(PARAM/'mod.py'),'-top','baseline.prmtop','-out','published.prmtop','-method','B3LYP'],
                           cwd=out,env=env,stdout=log,stderr=subprocess.STDOUT,timeout=120)
    if run.returncode or not (out/'published.prmtop').exists():raise ValueError('Published mandatory postprocessor failed')
    (out/'published.inpcrd').write_bytes((out/'baseline.inpcrd').read_bytes())
    # Record independent native type table changes. Never omit CMAP/pair terms
    # merely because a standard-term validator does not support them.
    modified=parmed.load_file(str(out/'published.prmtop'))
    expected_pairs={'O3':(56958.78,69.26),'P':(133097.78,154.94),'OS':(26407.56,66.18),'O2':(56958.78,69.26),'OY':(26407.56,66.18)}
    pair_changes=[];pd=native.parm_data;md=modified.parm_data;ntypes=native.ptr('ntypes')
    for i in range(ntypes):
        for j in range(i,ntypes):
            idx=pd['NONBONDED_PARM_INDEX'][i*ntypes+j]-1
            before=[pd['LENNARD_JONES_ACOEF'][idx],pd['LENNARD_JONES_BCOEF'][idx]]
            after=[md['LENNARD_JONES_ACOEF'][idx],md['LENNARD_JONES_BCOEF'][idx]]
            if before!=after:
                types_i=sorted({a.type for a in native.atoms if a.nb_idx==i+1})
                types_j=sorted({a.type for a in native.atoms if a.nb_idx==j+1})
                pair_changes.append({'native_lj_type_indices':[i+1,j+1],'atom_types':[types_i,types_j],
                                     'before_A_B_amber':before,'after_A_B_amber':after,
                                     'affected_atom_indices':[[a.idx for a in native.atoms if a.nb_idx==i+1],[a.idx for a in native.atoms if a.nb_idx==j+1]]})
    expected_indices={}
    mg_type=next(a.nb_idx for a in native.atoms if a.type=='Mg2+')
    for kind,coefficients in expected_pairs.items():
        kinds={a.nb_idx for a in native.atoms if a.type==kind}
        if len(kinds)!=1:raise ValueError('Published script assumes one native LJ type per atom-type name')
        pair=tuple(sorted((mg_type,next(iter(kinds)))))
        if pair in expected_indices and expected_indices[pair]!=coefficients:
            raise ValueError('Published atom-type overrides conflict through a shared LJ table index')
        expected_indices[pair]=coefficients
    if {tuple(x['native_lj_type_indices']) for x in pair_changes}!=set(expected_indices):
        raise ValueError('Published native pair table changed an unrequested index')
    for change in pair_changes:
        if not np.allclose(change['after_A_B_amber'],expected_indices[tuple(change['native_lj_type_indices'])],atol=1e-7,rtol=0):
            raise ValueError('Published pair coefficients differ')
    if md['CMAP_COUNT']!=[1,1] or md['CMAP_RESOLUTION']!=[24] or len(md['CMAP_PARAMETER_01'])!=576:
        raise ValueError('Missing complete published CMAP')
    atom_types={a.idx:a for a in native.atoms};cmap_atoms=[atom_types[i-1] for i in md['CMAP_INDEX'][:5]]
    if [a.name for a in cmap_atoms]!=['O3B','PB','O3A','PA',"O5'"] or len({a.residue.idx for a in cmap_atoms})!=1:
        raise ValueError('CMAP atom mapping differs')
    changed_flags=[key for key in pd if pd[key]!=md.get(key)]
    if set(changed_flags)-{'LENNARD_JONES_ACOEF','LENNARD_JONES_BCOEF'}:
        raise ValueError(f'Published postprocessor changed unrelated native arrays: {changed_flags}')
    protein=parmed.load_file(str(out/'protein-reference.prmtop'))
    # Exact ff14SB identity and parameter arrays for the protein-only unit must
    # remain equal after nucleotide parameters were loaded into the same LEaP.
    if len(protein.atoms)!=protein_atoms:raise ValueError('Protein-only native atom count changed')
    for a,b in zip(protein.atoms,native.atoms[:protein_atoms],strict=True):
        if (a.name,a.type,a.charge,a.mass,a.epsilon,a.rmin)!=(b.name,b.type,b.charge,b.mass,b.epsilon,b.rmin):
            raise ValueError('Loading ATP files changed protein nonbonded parameters')
    def terms(top, attr, atom_count):
        result=[]
        for item in getattr(top,attr):
            atoms=[getattr(item,k) for k in ['atom1','atom2','atom3','atom4'] if hasattr(item,k)]
            if not all(a.idx<atom_count for a in atoms):continue
            typ=item.type
            names={'bonds':['k','req'],'angles':['k','theteq'],'dihedrals':['phi_k','per','phase','scee','scnb']}[attr]
            result.append((tuple(a.idx for a in atoms),tuple(float(getattr(typ,k)) for k in names),getattr(item,'improper',False)))
        return sorted(result)
    for attr in ['bonds','angles','dihedrals']:
        if terms(protein,attr,protein_atoms)!=terms(native,attr,protein_atoms):raise ValueError(f'ATP file changed protein {attr}')
    mg=[a for a in modified.atoms if a.atomic_number==12]
    if len(mg)!=1 or abs(mg[0].charge-2)>1e-7 or any(mg[0] in [b.atom1,b.atom2] for b in modified.bonds):
        raise ValueError('The chosen published Mg model must be nonbonded +2')
    atp=[a for a in native.atoms if a.residue.name=='ATP']
    if len(atp)!=43 or abs(sum(a.charge for a in atp)+4)>2e-6:raise ValueError('ATP4− template inventory changed')
    report={'status':'complete native research assembly; independent CMAP/pair transport and physical validation pending',
            'model':'Hu2024 final ACS SI ATP-B3 + mandatory B3LYP CMAP/Mg pair overrides; ff14SB/TIP3P/Li-Merz compromise Mg(II)',
            'native_atoms':len(native.atoms),'native_residues':len(native.residues),'protein_atoms':protein_atoms,
            'observed_heavy_atoms_retained':len(matched_native),'water_residues':sum(r.name=='WAT' for r in native.residues),
            'Mg_charge_e':mg[0].charge,'ATP_charge_e':sum(a.charge for a in atp),'total_charge_e':sum(a.charge for a in native.atoms),
            'maximum_native_input_rounding_angstrom':max_native_delta,'protein_internal_parameters_unchanged':True,
            'cmap':{key:md[key] for key in ['CMAP_COUNT','CMAP_RESOLUTION','CMAP_INDEX','CMAP_PARAMETER_01']},
            'pair_specific_changes':pair_changes,'pair_scope':'The literal published script changes three native table pairs: O/O2/O3 share one LJ index; OS/OY share another. Thus all protein carbonyl O and carboxylate O2 are affected, not only ATP. This scope is preserved and is a scientific transferability question, not evidence of protein-site accuracy.',
            'license':'Final source SI CC BY-NC 4.0; not unrestricted open-source redistribution',
            'source_files':{str(p.relative_to(ROOT)):sha(p) for p in [SOURCE,PARAM/'ATP-B3.prepi',PARAM/'ATP-B3.frcmod',PARAM/'mod.py',Path(__file__)]},
            'output_hashes':{name:sha(out/name) for name in ['baseline.prmtop','baseline.inpcrd','published.prmtop','source-atom-map.json']},
            'dynamics':'not run','quantum_calculations':'not run','application_integration':False}
    dump(out/'native-assembly-report.json',report)
    print(json.dumps({k:v for k,v in report.items() if k not in ['cmap','pair_specific_changes']},indent=2))


if __name__=='__main__':main()
