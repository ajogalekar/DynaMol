"""Prepare/run classic Amber RESP with frozen canonical caps and backbone.

Inputs must be the explicit capped product from covalent_cap_spike.py. An ESP
file must contain the same atoms in the same order, in native RESP format.
This module does not generate QM evidence or claim a fit is chemically valid.
"""
from __future__ import annotations
import argparse
import hashlib
import json
import os
import subprocess
import shutil
from pathlib import Path
import numpy as np
from rdkit import Chem
from parmed.amber import AmberOFFLibrary

ROOT=Path(__file__).resolve().parents[3]


def prepare(folder):
    folder=folder.resolve();data=json.loads((folder/'capped-adduct.json').read_text())
    mol=Chem.SDMolSupplier(str(folder/'capped-adduct.sdf'),removeHs=False)[0]
    if mol is None or mol.GetNumAtoms()!=len(data['atom_map']):raise ValueError('Capped atom inventory is inconsistent.')
    records=data['atom_map'];names={r['name']:r['index'] for r in records}
    files={key:ROOT/'.tools/ambertools/dat/leap/lib'/name for key,name in
           [('CYS','amino12.lib'),('ACE','aminont12.lib'),('NME','aminoct12.lib')]}
    templates={key:{a.name:a for a in AmberOFFLibrary.parse(str(path))[key].atoms} for key,path in files.items()}
    heavy={'ACE_CH3':('ACE','CH3'),'ACE_C':('ACE','C'),'ACE_O':('ACE','O'),
           'NME_N':('NME','N'),'NME_CH3':('NME','C'),
           'N':('CYS','N'),'CA':('CYS','CA'),'C':('CYS','C'),'O':('CYS','O')}
    hsource={'ACE_CH3':('ACE','H1'),'NME_N':('NME','H'),'NME_CH3':('NME','H1'),
             'N':('CYS','H'),'CA':('CYS','HA')}
    fixed={};atom_sources={}
    for record in records:
        source=heavy.get(record['name'])
        if record['role']=='generated-hydrogen':
            source=hsource.get(records[record['parent_index']]['name'])
        if source:
            entry=templates[source[0]][source[1]]
            fixed[record['index']]=float(entry.charge)
            atom_sources[record['index']]={'template':source[0],'atom':source[1],'type':entry.type,'charge_e':float(entry.charge)}
    cap_indices=[r['index'] for r in records if r['role']=='temporary-cap' or
                 (r['role']=='generated-hydrogen' and records[r['parent_index']]['role']=='temporary-cap')]
    if not all(index in fixed for index in cap_indices):raise ValueError('Some cap atoms lack a canonical fixed charge.')
    if abs(sum(fixed[index] for index in cap_indices))>1e-8:raise ValueError('Canonical cap charges are not neutral.')
    qin=np.zeros(mol.GetNumAtoms())
    for i,q in fixed.items():qin[i]=q
    ivary1=[-1 if i in fixed else 0 for i in range(mol.GetNumAtoms())]
    # Conservative stage2: refine only rotatable methyl groups and equivalence
    # their H charges; do not equate potentially diastereotopic methylene H.
    ivary2=[-1]*mol.GetNumAtoms();methyl=[]
    for atom in mol.GetAtoms():
        hs=[a.GetIdx() for a in atom.GetNeighbors() if a.GetAtomicNum()==1]
        if atom.GetAtomicNum()==6 and len(hs)==3 and atom.GetIdx() not in fixed and not(set(hs)&fixed.keys()):
            ivary2[atom.GetIdx()]=0;ivary2[hs[0]]=0
            for index in hs[1:]:ivary2[index]=hs[0]+1
            methyl.append([atom.GetIdx(),*hs])
    out=folder/'resp';out.mkdir(exist_ok=True)
    for stage,ivary,weight in [(1,ivary1,.0005),(2,ivary2,.001)]:
        lines=[f'DynaMol isolated capped adduct RESP stage{stage}',
               f' &cntrl nmol=1, iqopt=2, ihfree=1, irstrnt=1, qwt={weight}, ioutopt=1, &end',
               '1.0','6OIM connected capped product',f'{data["formal_charge"]:5d}{mol.GetNumAtoms():5d}']
        lines += [f'{a.GetAtomicNum():5d}{ivary[a.GetIdx()]:5d}' for a in mol.GetAtoms()]
        (out/f'stage{stage}.in').write_text('\n'.join(lines)+'\n\n\n')
    (out/'canonical.qin').write_text('\n'.join(''.join(f'{x:10.6f}' for x in qin[i:i+8]) for i in range(0,len(qin),8))+'\n')
    report={'source':'AmberTools ff14SB residue libraries, exact atoms matched by cap/backbone role',
            'source_sha256':{str(path.relative_to(ROOT)):hashlib.sha256(path.read_bytes()).hexdigest() for path in files.values()},
            'atom_count':mol.GetNumAtoms(),'formal_charge':data['formal_charge'],
            'frozen_atoms':[{'index':i,'capped_name':records[i]['name'],**atom_sources[i]} for i in sorted(fixed)],
            'cap_indices':cap_indices,'cap_charge_e':sum(fixed[i] for i in cap_indices),
            'stage2_methyl_equivalence_groups':methyl,
            'stage2_policy':'Frozen stage1 charges except noncanonical methyl C/H; methyl H equivalenced, diastereotopic methylene H not forced equal.',
            'scope':'Inputs and exact boundary constraints only. ESP provenance, conformation treatment, residuals and chemistry validation are separate requirements.'}
    (out/'constraints.json').write_text(json.dumps(report,indent=2)+'\n')
    return out,report


def run(folder,esp,executable=None):
    out,report=prepare(folder);amber=ROOT/'.tools/ambertools'
    selection=ROOT/'build/advanced-chemistry/resp-runtime/selected-runtime.json'
    if executable is None and selection.exists():
        selected=json.loads(selection.read_text())
        executable=Path(selected['executable'])
        if not selected.get('validated') or hashlib.sha256(executable.read_bytes()).hexdigest()!=selected['sha256']:
            raise ValueError('Selected audit RESP runtime lacks validation or its executable hash changed.')
    executable=Path(executable or amber/'bin/resp').resolve()
    executable_hash=hashlib.sha256(executable.read_bytes()).hexdigest()
    env=dict(os.environ,AMBERHOME=str(amber),DYLD_FALLBACK_LIBRARY_PATH=str(amber/'lib'),OMP_NUM_THREADS='2',OPENBLAS_NUM_THREADS='2')
    # Classic RESP has short Fortran filename fields. Keep every file argument
    # local, even when the user's application workspace has a long path.
    shutil.copyfile(esp,out/'esp.dat')
    for stage in (1,2):
        qin='canonical.qin' if stage==1 else 'stage1.qout'
        command=[str(executable),'-O','-i',f'stage{stage}.in','-o',f'stage{stage}.out',
                 '-p',f'stage{stage}.punch','-q',qin,'-t',f'stage{stage}.qout','-e','esp.dat','-s',f'stage{stage}.esout']
        process=subprocess.run(command,cwd=out,env=env,capture_output=True,text=True,timeout=120)
        (out/f'stage{stage}.log').write_text(process.stdout+'\n'+process.stderr)
        if process.returncode or not (out/f'stage{stage}.qout').exists():raise RuntimeError(f'Native RESP stage{stage} failed; inspect its log/output.')
    values=np.array([float(x) for x in (out/'stage2.qout').read_text().split()])
    fixed_error=max(abs(values[x['index']]-x['charge_e']) for x in report['frozen_atoms'])
    checks={'esp_sha256':hashlib.sha256(esp.read_bytes()).hexdigest(),'charges':values.tolist(),
            'resp_executable':str(executable),'resp_executable_sha256':executable_hash,
            'fitter_source_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
            'charge_sum_e':float(values.sum()),'maximum_fixed_charge_error_e':float(fixed_error),
            'cap_charge_e':float(values[report['cap_indices']].sum()),
            'constraint_checks_passed':bool(len(values)==report['atom_count'] and np.isfinite(values).all() and fixed_error<1e-6 and abs(values.sum()-report['formal_charge'])<1e-4),
            'scope':'Constraint and native-fitting execution checks only, not an accuracy verdict.'}
    (out/'fit-constraints.json').write_text(json.dumps(checks,indent=2)+'\n')
    if not checks['constraint_checks_passed']:raise ValueError('RESP charge constraints failed.')
    return checks


if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('folder',type=Path);parser.add_argument('--esp',type=Path);parser.add_argument('--executable',type=Path)
    args=parser.parse_args()
    result=run(args.folder,args.esp,args.executable) if args.esp else prepare(args.folder)[1]
    print(json.dumps({k:v for k,v in result.items() if k not in ['charges','frozen_atoms']},indent=2))
