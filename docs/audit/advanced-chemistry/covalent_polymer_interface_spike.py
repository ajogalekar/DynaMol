"""Explicit native peptide-interface assembly for an isolated covalent adduct.

This audit helper requires an explicit charge-evidence kind. A synthetic ESP
fixture can exercise topology mechanics, but can never yield usable parameters.
It does not prepare a full structure, integrate a metal site, or enable the app.
"""
from __future__ import annotations
import argparse
import hashlib
import json
import math
import os
import shutil
import subprocess
from pathlib import Path

import numpy as np
import parmed

ROOT=Path(__file__).resolve().parents[3]


def run(base,interface,fitfolder,output,evidence_kind):
    base,interface,fitfolder,output=[path.resolve() for path in [base,interface,fitfolder,output]]
    output.mkdir(parents=True,exist_ok=False)
    if evidence_kind not in ('synthetic_execution_fixture','quantum_parameter_research'):
        raise ValueError('Charge evidence must be explicitly labeled.')
    mapping=json.loads((base/'capped-adduct.json').read_text())['atom_map']
    constraints=json.loads((fitfolder/'resp/constraints.json').read_text())
    fit=json.loads((fitfolder/'resp/fit-constraints.json').read_text())
    charges=np.array(fit['charges']);remove=set(constraints['cap_indices'])
    if not fit['constraint_checks_passed'] or len(charges)!=len(mapping) or not np.isfinite(charges).all():raise ValueError('Invalid native charge fit inventory.')
    if abs(float(charges[list(remove)].sum()))>1e-6:raise ValueError('Cap removal changes the model charge.')
    for row in constraints['frozen_atoms']:
        if abs(charges[row['index']]-row['charge_e'])>1e-6:raise ValueError('Canonical boundary charge changed.')
    atoms=[];bonds=[];section=None
    for line in (interface/'hybrid.mol2').read_text().splitlines():
        if line.startswith('@<TRIPOS>'):section=line;continue
        if not line.strip():continue
        if section=='@<TRIPOS>ATOM':atoms.append(line.split())
        if section=='@<TRIPOS>BOND':bonds.append(line.split())
    if len(atoms)!=len(mapping):raise ValueError('Native atom map length differs.')
    retained=[i for i in range(len(atoms)) if i not in remove];new={old:i+1 for i,old in enumerate(retained)}
    names={};canonical={row['index']:row for row in constraints['frozen_atoms']}
    for old in retained:
        record=mapping[old]
        if old in canonical:names[old]=canonical[old]['atom']
        elif record.get('resname')=='CYS':names[old]=record['atom']
        elif record['role']=='generated-hydrogen':names[old]=f'H{old+1}'
        else:names[old]=f'L{old+1:03d}'
    if len(set(names.values()))!=len(names):raise ValueError('Template atom names are not unique.')
    selected_bonds=[(new[int(row[1])-1],new[int(row[2])-1],row[3]) for row in bonds if int(row[1])-1 in new and int(row[2])-1 in new]
    lines=['@<TRIPOS>MOLECULE','COV',f'{len(retained)} {len(selected_bonds)} 1 0 0','SMALL','USER_CHARGES','','@<TRIPOS>ATOM']
    for old in retained:
        row=atoms[old];lines.append(f'{new[old]} {names[old]} {row[2]} {row[3]} {row[4]} {row[5]} 1 COV {charges[old]:.10f}')
    lines+=['@<TRIPOS>BOND']+[f'{i+1} {a} {b} {order}' for i,(a,b,order) in enumerate(selected_bonds)]
    lines+=['@<TRIPOS>SUBSTRUCTURE','1 COV 1 RESIDUE 1 A COV 1']
    (output/'covalent-residue.mol2').write_text('\n'.join(lines)+'\n')
    shutil.copyfile(interface/'gaff-adduct.frcmod',output/'gaff-adduct.frcmod')
    (output/'leap.in').write_text('source leaprc.protein.ff14SB\nsource leaprc.gaff2\nloadamberparams gaff-adduct.frcmod\nCOV=loadmol2 covalent-residue.mol2\nset COV head COV.1.N\nset COV tail COV.1.C\nset COV.1 restype protein\nmodel=sequence {ACE ALA COV GLY NME}\ncheck model\nsaveamberparm model peptide.prmtop peptide.inpcrd\nquit\n')
    amber=ROOT/'.tools/ambertools';env=dict(os.environ,AMBERHOME=str(amber),DYLD_FALLBACK_LIBRARY_PATH=str(amber/'lib'))
    process=subprocess.run([str(amber/'bin/tleap'),'-f','leap.in'],cwd=output,env=env,capture_output=True,text=True,timeout=60)
    (output/'tleap.log').write_text(process.stdout+'\n'+process.stderr)
    if process.returncode or not (output/'peptide.prmtop').exists():raise RuntimeError('Native peptide assembly failed; inspect tleap.log.')
    model=parmed.load_file(str(output/'peptide.prmtop'))
    residues=[r for r in model.residues if r.name=='COV']
    if len(residues)!=1 or len(residues[0].atoms)!=len(retained):raise ValueError('Native covalent residue identity/inventory changed.')
    residue=residues[0];central={a.name:a for a in residue.atoms}
    if set(central)!=set(names.values()):raise ValueError('Native covalent atom names changed.')
    native_bonds={tuple(sorted((bond.atom1.name,bond.atom2.name))) for bond in model.bonds if bond.atom1.residue is residue and bond.atom2.residue is residue}
    expected_bonds={tuple(sorted((names[retained[a-1]],names[retained[b-1]]))) for a,b,_ in selected_bonds}
    if native_bonds!=expected_bonds:raise ValueError('Native assembly changed connected adduct bonds.')
    crosslinks=[];counts={}
    for group,count in [('bonds',2),('angles',3),('dihedrals',4)]:
        terms=[]
        for term in getattr(model,group):
            members=[getattr(term,f'atom{i+1}') for i in range(count)]
            if any(a.residue is residue for a in members) and any(a.residue is not residue for a in members):
                if term.type is None:raise ValueError('A cross-boundary parameter is missing.')
                values={name:float(getattr(term.type,name)) for name in (('req','k') if group=='bonds' else ('theteq','k') if group=='angles' else ('per','phase','phi_k'))}
                if not all(math.isfinite(x) for x in values.values()):raise ValueError('A cross-boundary parameter is nonfinite.')
                terms.append({'atoms':[[a.residue.idx,a.residue.name,a.name] for a in members],'native_parameters':values})
        counts[group]=len(terms);crosslinks.append({'kind':group,'terms':terms})
    if counts['bonds']!=2:raise ValueError('The central residue must have exactly two explicit peptide connections.')
    peptide_charge=float(sum(a.charge for a in model.atoms))
    if abs(peptide_charge-round(peptide_charge))>1e-4:raise ValueError('Native peptide total charge is nonintegral after neutral cap removal.')
    report={'assembled':True,'usable_parameters':False,'charge_evidence_kind':evidence_kind,
        'scope':'Native topology assembly and complete cross-boundary mechanics only; synthetic charges are not physical parameters. Full quantum model quality, force validation and original complex preparation remain separate.',
        'charge_fit_sha256':hashlib.sha256((fitfolder/'resp/fit-constraints.json').read_bytes()).hexdigest(),
        'native_prmtop_sha256':hashlib.sha256((output/'peptide.prmtop').read_bytes()).hexdigest(),
        'removed_cap_atom_indices':sorted(remove),'remaining_adduct_atoms':len(retained),'total_peptide_atoms':len(model.atoms),
        'peptide_charge_e':peptide_charge,'cross_boundary_term_counts':counts,'cross_boundary_terms':crosslinks,
        'adduct_mapping':[{'capped_index':i,'template_name':names[i],'peptide_index':central[names[i]].idx,'source':mapping[i]} for i in retained]}
    (output/'result.json').write_text(json.dumps(report,indent=2)+'\n')
    return report


if __name__=='__main__':
    parser=argparse.ArgumentParser()
    for name in ('base','interface','fitfolder','output'):parser.add_argument(name,type=Path)
    parser.add_argument('--evidence-kind',required=True,choices=['synthetic_execution_fixture','quantum_parameter_research'])
    args=parser.parse_args();print(json.dumps({k:v for k,v in run(args.base,args.interface,args.fitfolder,args.output,args.evidence_kind).items() if k not in ['cross_boundary_terms','adduct_mapping']},indent=2))
