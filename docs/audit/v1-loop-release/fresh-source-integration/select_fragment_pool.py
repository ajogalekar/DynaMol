"""Select coherent seeds from a bounded native pool against a frozen snapshot.

Private integration helper, not an application admission path. Ranking is
case-independent and never waives final complete-complex validation.
"""
import hashlib
import importlib.util
import json
from pathlib import Path
import sys

import numpy as np
from openmm import unit
from backend.loop_backbone_reference import rama_report
from backend.loop_context import outer_torsion_anchor_indices, outside_observed_torsions
from backend.loop_snapshot import load_snapshot, _read, _json
from backend.residue_identity import residue_key

spec=importlib.util.spec_from_file_location('snapshot_validator',Path(__file__).with_name('validate_snapshot.py'))
validator=importlib.util.module_from_spec(spec);spec.loader.exec_module(validator)
BB=('N','CA','C','O')


def fingerprint(value):
    return hashlib.sha256(json.dumps(value,sort_keys=True,separators=(',',':'),allow_nan=False).encode()).hexdigest()


def proper_fit(candidate_anchors,reference_anchors):
    first=np.asarray(candidate_anchors,dtype=np.float64)
    second=np.asarray(reference_anchors,dtype=np.float64)
    if first.shape!=second.shape or first.ndim!=2 or first.shape[1]!=3 or len(first)<3 or not np.isfinite([first,second]).all():
        raise ValueError('Malformed rigid placement anchors')
    u,s,vt=np.linalg.svd((first-first.mean(0)).T@(second-second.mean(0)))
    if s[1]<1e-12:
        raise ValueError('Degenerate rigid placement anchors')
    correction=np.diag([1.,1.,np.linalg.det(u@vt)])
    rotation=u@correction@vt
    translation=second.mean(0)-first.mean(0)@rotation
    if not np.allclose(rotation.T@rotation,np.eye(3),atol=1e-12) or not np.isclose(np.linalg.det(rotation),1.,atol=1e-12):
        raise ValueError('Placement must preserve chirality')
    return rotation,translation


def select(request):
    if set(request)!={'snapshot','snapshot_sha256','pool','pool_sha256','native_input','native_input_sha256'}:
        raise ValueError('Exact snapshot, pool and native-input hash bindings are required')
    snapshot=load_snapshot(request['snapshot'],expected_sha256=request['snapshot_sha256'])
    pool=_json(_read(Path(request['pool']),request['pool_sha256']))
    config=_json(_read(Path(request['native_input']),request['native_input_sha256']))
    sources={f.name:f for f in snapshot.source_files}
    original_config=_json(_read(sources['loop-model-input.json'].path,sources['loop-model-input.json'].sha256))
    source_sha=sources['loop-model-context.pdb'].sha256
    if config['chains']!=original_config['chains'] or config.get('source_sha256')!=source_sha or pool.get('source_sha256')!=source_sha:
        raise ValueError('Pool sequence/context differs from the frozen complete snapshot')
    if pool.get('status')!='candidate' or pool.get('schema_version')!=1:
        raise ValueError('Expected a bounded native candidate pool')
    top,_=validator.topology_from_snapshot(snapshot,sources['prepared-topology.cif'].path)
    atoms=list(top.atoms());index={(*residue_key(a.residue),a.name):a.index for a in atoms}
    rows={(c['chain_id'],i+1):r for c in config['chains']for i,r in enumerate(c['residues'])}
    def key(chain,pos):
        r=rows[chain,pos]
        return (chain,r['resid'],r['insertion_code'],r['residue'])
    def atom_id(chain,pos,name):return (*key(chain,pos),name)
    modeled=set(map(tuple,snapshot.modeled_keys))
    covered=set();reports=[];selections={}
    for attempt in pool['attempts']:
        catalog=attempt.get('candidate_backbones')
        if not isinstance(catalog,list) or not 0<len(catalog)<=40:
            raise ValueError('Each gap needs a nonempty pool bounded at 40 native candidates')
        example=catalog[0];chain=example['chain'];left=example['original_start'];right=example['original_end']
        if not isinstance(chain,str) or type(left)is not int or type(right)is not int or not 1<=right-left-1<=12:
            raise ValueError('Malformed original sequence gap')
        segment={key(chain,p)for p in range(left+1,right)}
        flanks={key(chain,left),key(chain,right)}
        if not segment<=modeled or covered&segment or any(rows[chain,p]['observed']for p in range(left+1,right)) or not all(rows[chain,p]['observed']for p in (left,right)):
            raise ValueError('Pool gap does not match exact original missing identities')
        covered|=segment
        selected=segment|flanks
        selected_ids=[index[atom_id(chain,p,name)]for p in range(left,right+1)for name in BB]
        fixed=sorted(outer_torsion_anchor_indices(top,segment,flanks))
        if len(fixed)!=4 or not set(fixed)<=set(selected_ids):
            raise ValueError('Unexpected connected outer torsion support')
        observations=[index[atom_id(chain,p,name)]for p in (left,right)for name in BB]
        candidates=[]
        for candidate in catalog:
            if (candidate['chain'],candidate['original_start'],candidate['original_end'])!=(chain,left,right):
                raise ValueError('Mixed original gap identities in native catalog')
            start=candidate['start'];sequence=candidate['sequence']
            if type(start)is not int or sequence!=''.join(rows[chain,p]['one_letter']for p in range(start,start+len(sequence))):
                raise ValueError('Native candidate sequence mismatch')
            expected={(p,name)for p in range(start,start+len(sequence))for name in BB}
            coordinates={}
            for row in candidate['atoms']:
                pair=(row['position'],row['name']);values=np.asarray(row['xyz_nm'],dtype=np.float64)
                if pair not in expected or pair in coordinates or values.shape!=(3,) or not np.isfinite(values).all():
                    raise ValueError('Malformed or ambiguous candidate backbone atoms')
                coordinates[pair]=values
            if set(coordinates)!=expected or not set((p,name)for p in range(left,right+1)for name in BB)<=expected:
                raise ValueError('Candidate backbone inventory is incomplete')
            xyz=snapshot.prepared_xyz_nm.copy()
            for p in range(left,right+1):
                for name in BB:xyz[index[atom_id(chain,p,name)]]=coordinates[p,name]
            rotation,translation=proper_fit(xyz[fixed],snapshot.prepared_xyz_nm[fixed])
            residual=np.linalg.norm(xyz[fixed]@rotation+translation-snapshot.prepared_xyz_nm[fixed],axis=1)
            xyz[selected_ids]=xyz[selected_ids]@rotation+translation
            xyz[fixed]=snapshot.prepared_xyz_nm[fixed]
            outside=outside_observed_torsions(top,snapshot.prepared_xyz_nm,xyz,segment,flanks)
            if not outside['all_defining_coordinates_exactly_preserved']:
                raise ValueError('Candidate placement changed an outside torsion definition')
            reference=rama_report(top,xyz*unit.nanometer,selected)
            new_peptides=[r for r in reference['rows']if tuple(r['residue'])!=key(chain,left)]
            peptide_ok=len(new_peptides)==len(segment)+1 and all(
                (r['residue'][3]=='PRO' or abs(r['preceding_omega_degrees'])>90+1e-6)
                and min(abs(r['preceding_omega_degrees']),180-abs(r['preceding_omega_degrees']))<=35
                for r in new_peptides)
            displacement=np.linalg.norm(xyz[observations]-snapshot.prepared_xyz_nm[observations],axis=1)*10
            seed_ok=reference['all_selected_scored_without_outliers'] and peptide_ok
            candidates.append({'rank':candidate['rank'],'candidate_sha256':fingerprint(candidate),
                'seed_backbone_qualified':bool(seed_ok),'backbone_reference':reference,
                'new_peptide_states_and_planarity_passed':bool(peptide_ok),
                'maximum_observed_backbone_seed_displacement_angstrom':float(max(displacement)),
                'rms_observed_backbone_seed_displacement_angstrom':float(np.sqrt(np.mean(displacement**2))),
                'fit':{'rotation':rotation.tolist(),'translation_nm':translation.tolist(),
                       'anchor_rms_residual_angstrom':float(np.sqrt(np.mean(residual**2))*10)},
                'complete_candidate_accepted':False})
        eligible=[r for r in candidates if r['seed_backbone_qualified']]
        ranking=lambda r:(r['maximum_observed_backbone_seed_displacement_angstrom'],r['rms_observed_backbone_seed_displacement_angstrom'],r['candidate_sha256'])
        eligible.sort(key=ranking)
        if not eligible:
            reports.append({'gap':attempt['gap'],'candidates':candidates,'selected':None})
            continue
        chosen=eligible[0]
        catalog_record=next(c for c in catalog if fingerprint(c)==chosen['candidate_sha256'])
        selection_key=chain+':'+str(left)
        if selection_key in selections:raise ValueError('Duplicate original gap selection')
        selections[selection_key]={'rank':chosen['rank'],'candidate_sha256':chosen['candidate_sha256'],
                                    'candidate':catalog_record,'placement':chosen['fit']}
        reports.append({'gap':attempt['gap'],'candidates':candidates,'selected':chosen['rank']})
    if covered!=modeled:raise ValueError('Catalog omitted originally missing residues')
    return {'status':'selected' if len(selections)==len(reports)else'unavailable',
        'inputs':request,'source_sha256':source_sha,'selections':selections,'gaps':reports,
        'method':'Proper component anchor fit, exact outside support, qualified starting backbone/peptide state, then minimum observed backbone seed displacement. No case identifiers or rank constants.',
        'source_movement_gate':'1 A is enforced after full refinement; this seed ranking does not waive it.',
        'scope':'Backbone candidate selection only. Retained chemistry, sidechains, geometry, chirality, energy and serialized complete output still require validation.',
        'app_default_enabled':False,'complete_candidate_accepted':False}


if __name__=='__main__':
    request=_json(_read(Path(sys.argv[1])))
    result=select(request)
    with Path(sys.argv[2]).open('x')as handle:
        json.dump(result,handle,indent=2,allow_nan=False);handle.write('\n')
    print(json.dumps({'status':result['status'],'selected_ranks':{key:value['rank']for key,value in result['selections'].items()}}))
    raise SystemExit(0 if result['status']=='selected'else 2)
