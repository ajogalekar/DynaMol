"""Independent saved-coordinate screen. No generator calls or refinement."""
import hashlib
import json
from pathlib import Path
import sys

import numpy as np
from openmm import app, unit

ROOT=Path(__file__).resolve().parents[4]
sys.path.insert(0,str(ROOT/'docs/audit/advanced-chemistry/covalent-v1-focused'))
sys.path.insert(0,str(ROOT/'docs/audit/v1-loop-release/native-sampling'))
from rama_reference import rama_report, residue_key, dihedral64
from screen_kic_backbones import geometry

RUN=Path.home()/'.cache/dynamol-research/v1-loop-release/hybrid-native-mc-v1'


def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def save(p,d):p.write_text(json.dumps(d,indent=2,allow_nan=False)+'\n')


def source_atoms(path):
    result={}
    for l in path.read_text().splitlines():
        if l[:6].strip() not in ('ATOM','HETATM'):continue
        key=(l[21].strip(),l[22:26].strip(),l[26].strip(),l[17:20].strip(),l[12:16].strip())
        if key in result:raise ValueError('Ambiguous original atom')
        result[key]=np.array([float(l[30:38]),float(l[38:46]),float(l[46:54])])/10
    return result


def topology_and_positions(config,source,backbone):
    topology=app.Topology();coords=[];keys=[]
    for chainrow in config['chains']:
        chain=topology.addChain(chainrow['chain_id']);previous=None
        for index,row in enumerate(chainrow['residues']):
            rk=(chain.id,row['resid'],row['insertion_code'],row['residue'])
            atomkeys=[rk+(a,) for a in ('N','CA','C','O')]
            if not any(k in source or k in backbone for k in atomkeys):
                previous=None;continue
            if not all(k in source or k in backbone for k in atomkeys):raise ValueError('Incomplete defining backbone')
            residue=topology.addResidue(row['residue'],chain,row['resid'],row['insertion_code']);atoms={}
            for k in atomkeys:
                element=app.element.nitrogen if k[-1]=='N' else app.element.oxygen if k[-1]=='O' else app.element.carbon
                atoms[k[-1]]=topology.addAtom(k[-1],element,residue);keys.append(k)
                coords.append(backbone[k] if k in backbone else source[k])
            for a,b in [('N','CA'),('CA','C'),('C','O')]:topology.addBond(atoms[a],atoms[b])
            if previous is not None:topology.addBond(previous,atoms['N'])
            previous=atoms['C']
    return topology,np.array(coords),keys


def screen(path,case):
    candidate=json.loads(path.read_text());config=json.loads(Path(case['input']).read_text())
    source=source_atoms(Path(case['source']));bb={tuple(a['identity']):np.array(a['xyz_nm']) for a in candidate['backbone_atoms']}
    if len(bb)!=len(candidate['backbone_atoms']):raise ValueError('Duplicate candidate identity')
    for a in candidate['backbone_atoms']:
        if np.max(np.abs(np.asarray(a['xyz_angstrom'])/10-bb[tuple(a['identity'])]))>1e-12:raise ValueError('Å/nm mismatch')
    top,xyz,keys=topology_and_positions(config,source,bb)
    modeled=set(map(tuple,candidate['modeled_residue_keys']))
    changed={i for i,k in enumerate(keys) if k[:4] in modeled or (k in source and np.linalg.norm(xyz[i]-source[k])*10>1e-4)}
    all_rama=rama_report(top,xyz*unit.nanometer)
    affected=modeled|{tuple(r['residue']) for r in all_rama['rows'] if changed.intersection(r['phi_atoms']+r['psi_atoms'])}
    rama=rama_report(top,xyz*unit.nanometer,affected)
    first,last=candidate['settings']['span']
    local=[{'identity':list(k),'xyz_angstrom':(xyz[i]*10).tolist()} for i,k in enumerate(keys)
           if k[0]=='A' and first-1<=int(k[1])<=last+1]
    gross=geometry(local)
    displacements=[np.linalg.norm(bb[k]-source[k])*10 for k in bb if k[:4] not in modeled]
    cap=max(displacements)<=1
    local_lookup={(int(r['identity'][1]),r['identity'][4]):r for r in local}
    observed_links=[]
    for link in gross['omega_checks']:
        n,m=link['first'],link['second']
        defining=[local_lookup[n,'CA'],local_lookup[n,'C'],local_lookup[m,'N'],local_lookup[m,'CA']]
        if not all(tuple(r['identity']) in source for r in defining):continue
        original_omega=dihedral64(np.array([source[tuple(r['identity'])] for r in defining]))
        actual_omega=link['degrees']
        observed_links.append({'first':n,'second':m,'source_degrees':original_omega,'candidate_degrees':actual_omega,
                               'basin_preserved':(-90<original_omega<90)==(-90<actual_omega<90)})
    omega=all(r['basin_preserved'] for r in observed_links)
    native_keys={tuple(a['identity']):np.array(a['xyz_nm']) for a in candidate.get('modeled_heavy_atoms',[])+candidate.get('observed_context_atoms',[])}
    overlaps={k for k in bb if k in native_keys}
    coherence=max([float(np.linalg.norm(bb[k]-native_keys[k])*10) for k in overlaps] or [0])
    context=set(map(tuple,candidate['context_residue_keys']))
    outside=sorted({tuple(a['identity'][:4]) for a in candidate.get('observed_context_atoms',[]) if tuple(a['identity'][:4]) not in context})
    report={'candidate':str(path),'candidate_sha256':sha(path),'id':candidate['id'],'case':candidate['case'],'generator':candidate['generator'],
            'seed':candidate['settings']['seed'],'span':candidate['settings']['span'],
            'observed_backbone_cap_pass':bool(cap),'maximum_observed_backbone_displacement_A':float(max(displacements)),
            'observed_omega_basins_preserved':omega,'gross_backbone_including_joins':gross,'rama':rama,
            'observed_omega_checks_including_outer_joins':observed_links,
            'affected_residue_keys':sorted(affected),'raw_vs_sidechain_backbone_maximum_difference_A':coherence,
            'observed_heavy_changes_outside_declared_context':outside,
            'plausible_backbone_screen_pass':bool(cap and omega and gross['accepted'] and rama['all_selected_scored_without_outliers'] and coherence<1e-4 and not outside),
            'sidechain_full_complex_validated':False,'app_ready':False}
    if 'heldout_control_source' in case:
        hidden=source_atoms(Path(case['heldout_control_source']))
        deviations=[np.sum((bb[k]-hidden[k])**2) for k in bb if k[:4] in modeled]
        report['heldout_modeled_backbone_RMSD_A']=float(np.sqrt(np.mean(deviations))*10)
        report['heldout_coordinates_used_for_generation_or_ranking']=False
    save(path.parent/'independent-backbone-screen.json',report)
    return report


def main():
    plan=json.loads((RUN/'plan.json').read_text());cases={c['case']:c for c in plan['cases']}
    paths=sorted(RUN.glob('*/*/native/*/candidate.json'));hashes={str(p):sha(p) for p in paths}
    reports=[screen(p,cases[p.relative_to(RUN).parts[0]]) for p in paths]
    assert hashes=={str(p):sha(p) for p in paths}
    compact=[]
    for r in reports:
        compact.append({k:v for k,v in r.items() if k not in ('gross_backbone_including_joins','rama')}|
                       {'gross_pass':r['gross_backbone_including_joins']['accepted'],
                        'rama_pass':r['rama']['all_selected_scored_without_outliers'],
                        'rama_outliers':[x['residue'] for x in r['rama']['outliers']]})
    result={'plan_sha256':sha(RUN/'plan.json'),'source_sha256':hashes,'candidate_count':len(reports),
            'plausible_backbone_count':sum(r['plausible_backbone_screen_pass'] for r in reports),
            'candidates':compact,'all_physical_models_unvalidated':True,
            'scope':'Local static backbone screen including joins/boundaries; no full-atom or app admission.'}
    save(RUN/'independent-screen-summary.json',result)
    save(Path(__file__).with_name('native-comparison-summary.json'),result)
    print(json.dumps({'candidate_count':len(reports),'plausible_backbone_count':result['plausible_backbone_count']}))


if __name__=='__main__':main()
