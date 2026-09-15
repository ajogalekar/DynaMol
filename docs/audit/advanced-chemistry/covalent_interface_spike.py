"""Native ff14SB-backbone/GAFF2-adduct typing feasibility, not app support.

Keeps joint capped AM1-BCC charges for this mechanics experiment. Actual
polymer integration must instead use reviewed constrained RESP charges,
explicit head/tail topology, cross-boundary validation and QM term assessment.
"""
from __future__ import annotations
import argparse
from collections import Counter
import hashlib
import json
import os
import re
import subprocess
from pathlib import Path

import parmed
import numpy as np

ROOT=Path(__file__).resolve().parents[3]


def compare_canonical_backbone(out,constraints):
    """Compare every wholly canonical-subset term with native ACE-CYS-NME."""
    hybrid=parmed.load_file(str(out/'hybrid.prmtop'))
    reference=parmed.load_file(str(out/'reference.prmtop'))
    reference_ids={(atom.residue.name,atom.name):atom.idx for atom in reference.atoms}
    mapping={};methyl_counts={}
    for row in constraints['frozen_atoms']:
        template,name=row['template'],row['atom']
        if template in ('ACE','NME') and name=='H1':
            methyl_counts[template]=methyl_counts.get(template,0)+1
            name=f'H{methyl_counts[template]}'
        mapping[row['index']]=reference_ids[(template,name)]
    if len(mapping)!=18 or len(set(mapping.values()))!=18:raise ValueError('Canonical reference atom mapping is incomplete or ambiguous.')
    def inventory(model,selected,indexmap):
        result={}
        for group,count,fields in [('bonds',2,('req','k')),('angles',3,('theteq','k')),('dihedrals',4,('per','phase','phi_k'))]:
            terms=[]
            for term in getattr(model,group):
                indices=tuple(getattr(term,f'atom{i+1}').idx for i in range(count))
                if not set(indices)<=selected:continue
                mapped=tuple(indexmap[i] for i in indices)
                values=tuple(round(float(getattr(term.type,name)),8) for name in fields)
                if group=='dihedrals':values+=bool(term.improper),bool(term.ignore_end)
                terms.append((min(mapped,mapped[::-1]),values))
            result[group]=Counter(terms)
        return result
    left=inventory(hybrid,set(mapping),mapping)
    right=inventory(reference,set(mapping.values()),{i:i for i in mapping.values()})
    checks={group:{'hybrid_terms':sum(left[group].values()),'reference_terms':sum(right[group].values()),
                    'exact_native_terms_match':left[group]==right[group],
                    'missing':list((right[group]-left[group]).elements()),'extra':list((left[group]-right[group]).elements())}
            for group in left}
    report={'accepted':all(row['exact_native_terms_match'] for row in checks.values()),'checks':checks,
            'scope':'Every mechanical term wholly within the 18 canonical cap/backbone atoms matches native ff14SB ACE-CYS-NME. Charges, modified sidechain/interface terms and polymer integration are separate requirements.',
            'parameter_rounding_digits':8,'canonical_index_to_reference_index':mapping,
            'hybrid_sha256':hashlib.sha256((out/'hybrid.prmtop').read_bytes()).hexdigest(),
            'reference_sha256':hashlib.sha256((out/'reference.prmtop').read_bytes()).hexdigest()}
    (out/'backbone-reference-validation.json').write_text(json.dumps(report,indent=2)+'\n')
    if not report['accepted']:raise ValueError('Hybrid mechanics changed a canonical backbone term.')
    return report


def run(base,out):
    base=base.resolve();out=out.resolve();out.mkdir(parents=True,exist_ok=False)
    constraints=json.loads((base/'resp/constraints.json').read_text())
    overrides={row['index']:row['type'] for row in constraints['frozen_atoms']}
    lines=(base/'capped-native/charged.mol2').read_text().splitlines();inside=False;changed=[]
    for i,line in enumerate(lines):
        if line.startswith('@<TRIPOS>'):inside=line=='@<TRIPOS>ATOM';continue
        if inside and line.strip():
            fields=line.split();index=int(fields[0])-1
            if index in overrides:
                changed.append({'index':index,'gaff2_type':fields[5],'protein_type':overrides[index]})
                fields[5]=overrides[index];lines[i]=' '.join(fields)
    (out/'hybrid.mol2').write_text('\n'.join(lines)+'\n')
    amber=ROOT/'.tools/ambertools';environment=dict(os.environ,AMBERHOME=str(amber),DYLD_FALLBACK_LIBRARY_PATH=str(amber/'lib'))
    environment['PATH']=str(amber/'bin')+os.pathsep+environment.get('PATH','')
    commands=[]
    def native(args):
        process=subprocess.run([str(amber/'bin'/args[0]),*args[1:]],cwd=out,env=environment,capture_output=True,text=True,timeout=60)
        log=f'{len(commands)+1:02d}-{args[0]}.log';(out/log).write_text(process.stdout+'\n'+process.stderr)
        commands.append({'arguments':args,'returncode':process.returncode,'log':log})
        if process.returncode:raise RuntimeError(f'Native command failed; inspect {log}')
    native(['parmchk2','-i','hybrid.mol2','-f','mol2','-o','protein-all.frcmod','-s','parm10','-frc','ff14SB','-a','Y'])
    native(['parmchk2','-i','hybrid.mol2','-f','mol2','-o','gaff-all.frcmod','-s','gaff2','-a','Y'])
    # Keep canonical protein terms solely from native ff14SB. Only GAFF/mixed
    # terms may enter the adduct supplement. Protein parmchk2 output is retained
    # for audit only: its unmarked placeholders/default impropers for GAFF atom
    # types are not legitimate overrides of the GAFF parameters.
    supplement=[];omitted=[];section=None;canonical=set(overrides.values())|{'X'}
    for line in (out/'gaff-all.frcmod').read_text().splitlines():
        if line.strip() in ('MASS','BOND','ANGLE','DIHE','IMPROPER','NONBON'):
            section=line.strip();supplement.append(line);continue
        if not line.strip() or section is None:
            supplement.append(line);continue
        width={'MASS':2,'BOND':5,'ANGLE':8,'DIHE':11,'IMPROPER':11}.get(section)
        types=[x.strip() for x in line[:width].split('-')] if width else [line.split()[0]]
        if all(kind in canonical for kind in types):
            omitted.append(line);continue
        invalid_torsion=False
        if section in ('DIHE','IMPROPER') and '-' in line[:11]:
            fields=line[11:].split()
            try:
                if section=='DIHE':
                    divisor,amplitude,phase,periodicity=map(float,fields[:4])
                    invalid_torsion=divisor<=0 or periodicity==0
                else:
                    amplitude,phase,periodicity=map(float,fields[:3])
                    invalid_torsion=periodicity<=0
            except (ValueError,IndexError):invalid_torsion=True
        if 'ATTN' in line or invalid_torsion:raise ValueError('An adduct/mixed parameter is missing or invalid; no replacement was accepted: '+line)
        supplement.append(line)
    (out/'gaff-adduct.frcmod').write_text('\n'.join(supplement)+'\n')
    (out/'leap.in').write_text('source leaprc.protein.ff14SB\nsource leaprc.gaff2\nloadamberparams gaff-adduct.frcmod\nadduct=loadmol2 hybrid.mol2\ncheck adduct\nsaveamberparm adduct hybrid.prmtop hybrid.inpcrd\nquit\n')
    native(['tleap','-f','leap.in'])
    (out/'reference-leap.in').write_text('source leaprc.protein.ff14SB\nref=sequence {ACE CYS NME}\ncheck ref\nsaveamberparm ref reference.prmtop reference.inpcrd\nquit\n')
    native(['tleap','-f','reference-leap.in'])
    comparison=compare_canonical_backbone(out,constraints)
    import openmm as mm
    from openmm import app,unit
    topology=app.AmberPrmtopFile(str(out/'hybrid.prmtop'))
    system=topology.createSystem(nonbondedMethod=app.NoCutoff,constraints=None)
    integrator=mm.VerletIntegrator(.001)
    context=mm.Context(system,integrator,mm.Platform.getPlatformByName('Reference'))
    context.setPositions(app.AmberInpcrdFile(str(out/'hybrid.inpcrd')).positions)
    state=context.getState(getEnergy=True,getForces=True)
    energy=float(state.getPotentialEnergy().value_in_unit(unit.kilojoule_per_mole))
    forces=np.asarray(state.getForces(asNumpy=True).value_in_unit(unit.kilojoule_per_mole/unit.nanometer))
    if not np.isfinite(energy) or not np.isfinite(forces).all():raise ValueError('Native model cannot create a finite OpenMM Context.')
    del context,integrator
    report={'scope':'Mixed ff14SB-backbone/GAFF2 connected capped-adduct typing and complete-term feasibility only; charges remain joint capped AM1BCC. Not constrainedRESP, polymer integration, or force-field accuracy validation.',
            'changed_atom_types':changed,'commands':commands,'canonical_gaff_overrides_omitted':omitted,
            'parameter_precedence':'Wholly canonical protein-type terms come exclusively from native ff14SB; only terms touching a GAFF atom type enter gaff-adduct.frcmod. Protein parmchk2 output is audit-only.',
            'gaff_missing_attention_terms':[line for line in (out/'gaff-all.frcmod').read_text().splitlines() if 'ATTN' in line],
            'native_topology_created':(out/'hybrid.prmtop').exists(),
            'native_openmm_context':{'created':True,'finite_energy_and_forces':True,'energy_kj_mol':energy},
            'canonical_backbone_reference':{'accepted':comparison['accepted'],'report':'backbone-reference-validation.json'}}
    (out/'result.json').write_text(json.dumps(report,indent=2)+'\n')
    if not report['native_topology_created'] or report['gaff_missing_attention_terms']:raise RuntimeError('No complete native interface model was created.')
    return report


if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('folder',type=Path);parser.add_argument('output',type=Path)
    args=parser.parse_args();result=run(args.folder,args.output)
    print(json.dumps({k:v for k,v in result.items() if k not in ['changed_atom_types','protein_override_omitted']},indent=2))
