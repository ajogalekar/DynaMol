"""Independent all-output static screen; failed native closures cannot pass."""
import importlib.util,json,hashlib,math
from pathlib import Path
import numpy as np
ROOT=Path('/Users/ashujo/Documents/Science/DynaMol')
RUN=Path('/Users/ashujo/.cache/dynamol-research/v1-loop-release/shifted-anchor-native-v1')
MODULE=ROOT/'docs/audit/v1-loop-release/hybrid-fallback/screen_native_comparison.py'
spec=importlib.util.spec_from_file_location('existing_screen',MODULE);sc=importlib.util.module_from_spec(spec);spec.loader.exec_module(sc)
def get(p):return json.loads(Path(p).read_text())
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def save(p,d):Path(p).write_text(json.dumps(d,indent=2,allow_nan=False)+'\n')
plan=get(RUN/'plan.json');cases={c['case']:c for c in plan['cases']}
assert all(sha(p)==h for p,h in plan['frozen_sha256'].items())
paths=sorted(RUN.glob('*/native/*/candidate.json'));hashes={str(p):sha(p) for p in paths}
trials=sorted(RUN.glob('*/native/*/trial-coordinates.json'));trialhashes={str(p):sha(p) for p in trials}
rows=[]
for attempt in plan['attempts']:
    case=cases[attempt['case']];folder=RUN/case['case']/'native'/attempt['name'];outcome=get(folder/'attempt-result.json');diagnostics=get(folder/'trial-coordinates.json')
    settings=get(folder/'attempt-plan.json');old=get(attempt['input_candidate']);source=sc.source_atoms(Path(case['source']))
    assert diagnostics['native_closure_calls']==1
    assert max(diagnostics['prefix_displacements_during_closure_A'].values())==0
    assert settings['input_candidate_sha256']==sha(attempt['input_candidate']) and settings['carbonyl_rotation_sequential']
    assert settings['observed_coordinate_restoration'] is False and settings['resampling'] is False
    assert old['source_sha256']==sha(case['source'])
    candidate=folder/'candidate.json';success=diagnostics['native_closure_success']
    sidechain_covered=False;coherence=None
    if success:
        c=get(candidate);assert c['source_sha256']==sha(case['source']) and c['native_closure_calls']==1
        assert c['modeled_residue_keys']==case['modeled_residue_keys'] and c['backbone_atoms']==diagnostics['recombined']
        assert c['sidechain_reconstruction_status']=='complete'
        assert c['id']==hashlib.sha256(json.dumps([case['case'],settings,c['backbone_atoms']],sort_keys=True).encode()).hexdigest()
        bb={tuple(a['identity']):np.array(a['xyz_nm']) for a in c['backbone_atoms']};modeled=set(map(tuple,c['modeled_residue_keys']))
        sa=c['modeled_heavy_atoms']+c['observed_context_atoms'];scatoms={tuple(a['identity']):np.array(a['xyz_nm']) for a in sa};assert len(sa)==len(scatoms)
        sidechain_covered=all(k in scatoms for k in bb if k[:4] in modeled)
        assert sidechain_covered
        coherence=float(max(np.linalg.norm(bb[k]-(scatoms[k] if k in scatoms else source[k]))*10 for k in bb))
    else:
        assert not candidate.exists()
        candidate=folder/'failed-backbone-for-screen.json'
        c={'id':attempt['name']+'-native-closure-failed','case':case['case'],'generator':'failed_shifted_anchor_ccd_trial',
           'settings':settings,'backbone_atoms':diagnostics['recombined'],'modeled_residue_keys':case['modeled_residue_keys'],
           'context_residue_keys':old['context_residue_keys'],'observed_internal_residue_keys':[],
           'native_closure_success':False,'app_ready':False,'physical_model_validated':False,
           'scope':'Failed closure coordinates for independent diagnosis only; not a generated candidate.'}
        save(candidate,c)
    r=sc.screen(candidate,case)
    r['backbone_only_diagnostic_pass']=r['plausible_backbone_screen_pass']
    r['native_closure_success']=success
    r['native_and_backbone_screen_pass']=bool(success and r['plausible_backbone_screen_pass'] and sidechain_covered and coherence<=1e-4)
    if not success:r['plausible_backbone_screen_pass']=False
    save(folder/'independent-backbone-screen.json',r)
    modeled=set(map(tuple,case['modeled_residue_keys']))
    compact={k:v for k,v in r.items() if k not in ['rama','gross_backbone_including_joins']}
    compact.update(name=attempt['name'],native_closure_success=success,native_closure_calls=1,
       gross_pass=r['gross_backbone_including_joins']['accepted'],gross_failures=[dict(check_type=t,**v) for t in ['bond_checks','angle_checks','omega_checks'] for v in r['gross_backbone_including_joins'][t] if not v['accepted']],
       rama_pass=r['rama']['all_selected_scored_without_outliers'],modeled_outliers=[x['residue'] for x in r['rama']['outliers'] if tuple(x['residue']) in modeled],
       boundary_outliers=[x['residue'] for x in r['rama']['outliers'] if tuple(x['residue']) not in modeled],
       boundary_reference_rows=[x for x in r['rama']['rows'] if tuple(x['residue']) not in modeled],
       native_sidechain_coverage=sidechain_covered,all_sidechain_backbone_max_difference_A=coherence,
       native_and_backbone_screen_pass=bool(success and r['plausible_backbone_screen_pass'] and sidechain_covered and coherence<=1e-4),
       locked_prefix_source_displacements_A=diagnostics['source_n_stem_N_CA_C_O_displacement_A']['locked'],
       locked_prefix_displacements_during_closure_A=diagnostics['prefix_displacements_during_closure_A'],
       virtual_anchor_N_CA_C_displacements_A=diagnostics['virtual_anchor_N_CA_C_displacements_A'],
       locked_carbonyl_orientation_error_radians=diagnostics['locked_orientation_error_radians'],
       first_join_locked=diagnostics['joins_locked'][0],first_join_recombined=diagnostics['joins_recombined'][0],
       source_candidate=attempt['input_candidate'],source_candidate_sha256=attempt['input_candidate_sha256'])
    blookup={tuple(v['identity']):np.array(v['xyz_nm']) for v in diagnostics['recombined']}
    last=tuple(diagnostics['recombined'][-1]['identity'][:4])
    compact['c_stem_N_CA_C_RMSD_A']=float(np.sqrt(np.mean([np.sum((blookup[last+(a,)]-source[last+(a,)])**2) for a in ['N','CA','C']]))*10)
    # A failed native closure has no sidechain acceptance; do not let a geometric diagnostic imply otherwise.
    if not success:compact['plausible_backbone_screen_pass']=False
    rows.append(compact)
assert len(rows)==9 and sum(r['native_closure_calls'] for r in rows)==9
assert hashes=={str(p):sha(p) for p in paths} and trialhashes=={str(p):sha(p) for p in trials}
assert all(sha(p)==h for p,h in plan['frozen_sha256'].items())
summary={'plan_sha256':sha(RUN/'plan.json'),'implementation_sha256':sha(__file__),'existing_screen_sha256':sha(MODULE),
 'candidate_sha256':hashes,'trial_sha256':trialhashes,'frozen_source_implementation_hashes_verified':True,
 'attempts':9,'native_closure_calls':9,'native_successes':sum(r['native_closure_success'] for r in rows),
 'native_and_backbone_screen_passes':sum(r['native_and_backbone_screen_pass'] for r in rows),'rows':rows,
 'scope':'One shifted-anchor native closure per saved candidate, then independent local static screen of every output including failed closures. Full complex remains unvalidated.',
 'app_ready':False,'physical_model_validated':False}
save(RUN/'independent-screen-summary.json',summary);save(Path(__file__).with_name('shifted-summary.json'),summary)
for r in rows:print(json.dumps({k:r[k] for k in ['case','name','native_closure_success','gross_pass','maximum_observed_backbone_displacement_A','observed_backbone_cap_pass','modeled_outliers','boundary_outliers','native_and_backbone_screen_pass']}))
