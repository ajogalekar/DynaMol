"""Private automatic bounded fragment selection, placement and full validation.

CLI matches the full-worker qualification seam. No PDB-specific decisions,
manual rank inputs, MD, parameter fitting or application enablement here.
"""
import argparse,hashlib,json,os,shutil,sys
from pathlib import Path
from backend.loop_process import supervise_loop_child
from backend.loop_snapshot import load_snapshot

REPO = Path('/Users/ashujo/Documents/Science/DynaMol')
PYTHON=Path('/Users/ashujo/.cache/dynamol-runtimes/covalent-v1-focused/bin/python')
PROMOD=Path('/Users/ashujo/Library/Application Support/DynaMol/Runtimes/promod3-dfd2a559551cd2f2')
REFERENCE=Path('/Users/ashujo/.cache/dynamol-runtimes/cctbx-rama-reference-v1')


def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def write(p,value):
    with Path(p).open('x')as h:json.dump(value,h,indent=2,allow_nan=False);h.write('\n')


def native_catalog_sources(destination):
    """Explicitly archived research adapter around the bounded native worker."""
    legacy=(REPO/'backend/promod_loop_worker.py').read_text()
    marker="        record['candidate_count']=len(alternatives)"
    catalog="""        record['candidate_backbones']=[{'rank':rank,'context':str(item[3]),'start':item[3].before.number.num,
            'original_start':original.before.number.num,'original_end':original.after.number.num,
            'chain':name,'sequence':item[3].full_seq,'native_score':item[1],
            'anchor_squared_displacement_sum':item[0],
            'atoms':[{'position':item[3].before.number.num+j,'name':atom,
                'xyz_nm':[v/10 for v in (getattr(item[2],'Get'+atom)(j).x,getattr(item[2],'Get'+atom)(j).y,getattr(item[2],'Get'+atom)(j).z)]}
                for j in range(len(item[2]))for atom in ('N','CA','C','O')]}
            for rank,item in enumerate(sorted(alternatives,key=lambda candidate:candidate[:2]))]
"""
    if legacy.count(marker)!=1:raise ValueError('Native catalog adapter insertion point changed')
    legacy=legacy.replace(marker,marker+'\n'+catalog)
    marker='            best=min(alternatives,key=lambda candidate:candidate[:2])'
    choose="""            selection_key=name+':'+str(original.before.number.num)
            expected=_RESEARCH_EXPECTED_SELECTED.get(selection_key)
            if expected is None:
                best=min(alternatives,key=lambda candidate:candidate[:2])
            else:
                matches=[r['rank']for r in record['candidate_backbones']
                    if all(r[k]==expected[k]for k in ('atoms','sequence','start'))]
                if not matches:
                    raise ValueError('Selected backbone is absent from the same bounded source pool')
                best=sorted(alternatives,key=lambda candidate:candidate[:2])[matches[0]]
                record['automatically_selected_rank']=matches[0]
"""
    if legacy.count(marker)!=1:raise ValueError('Native selection adapter insertion point changed')
    legacy=legacy.replace(marker,choose)
    worker=(REPO/'backend/promod_candidate_worker.py').read_text()
    marker='            attempts = legacy.context_fragment_search(model, native.modelling, native.loop.LoadFragDB(),'
    if worker.count(marker)!=1:raise ValueError('Native proposal adapter insertion point changed')
    worker=worker.replace(marker,"            legacy._RESEARCH_EXPECTED_SELECTED=config.get('research_expected_selected',{})\n"+marker)
    (destination/'promod_loop_worker.py').write_text(legacy)
    (destination/'promod_candidate_worker.py').write_text(worker)


def run(snapshot_path,snapshot_hash,output,seed):
    snapshot=load_snapshot(snapshot_path,expected_sha256=snapshot_hash)
    out=Path(output).resolve();out.mkdir(parents=True,exist_ok=False)
    implementation=out/'implementation';implementation.mkdir()
    native_catalog_sources(implementation)
    helper_names=['validate_snapshot.py','select_fragment_pool.py','place_fragment_proposal.py']
    for name in helper_names:shutil.copyfile(Path(__file__).with_name(name),implementation/name)
    shutil.copyfile(__file__,implementation/'run_qualified_fragment_search.py')
    pins={str(p):sha(p)for p in [*sorted(implementation.glob('*.py')),*sorted((REPO/'backend').glob('*.py'))]}
    write(out/'implementation-pins.json',pins)
    def verify():
        for p,digest in pins.items():
            if sha(p)!=digest:raise ValueError('Pinned search implementation changed: '+p)
    env=dict(os.environ)
    env.update(PYTHONPATH=os.pathsep.join([str(REPO),str(REFERENCE)]),OPENMM_CPU_THREADS='2',DYNAMOL_CPU_THREADS='2',
        PM3_OPENMM_CPU_THREADS='2',OPENBLAS_NUM_THREADS='2',OMP_NUM_THREADS='2',MKL_NUM_THREADS='2',VECLIB_MAXIMUM_THREADS='2',PYTHONDONTWRITEBYTECODE='1')
    def child(command,folder,timeout=180):
        verify();code=supervise_loop_child(command,folder,environment=env,timeout_seconds=timeout,check_cancel=verify)
        verify();return code
    sources={f.name:f for f in snapshot.source_files}
    native_input=json.loads(sources['loop-model-input.json'].path.read_text())
    native_input.update(input_pdb=str(sources['loop-model-context.pdb'].path),source_sha256=sources['loop-model-context.pdb'].sha256,
                        strategy='fragment',seed=seed,max_res_extension=2)
    pool=out/'pool';pool.mkdir();write(pool/'input.json',native_input)
    write(out/'plan.json',{'snapshot_sha256':snapshot_hash,'source_sha256':native_input['source_sha256'],
        'max_examined_native_candidates_per_gap':40,'native_context_extension':2,'seed':seed,
        'native_generations_maximum':2,'full_refinements_maximum':1,'source_files_sha256':pins,
        'selection':'Qualified placed backbone then minimum observed seed movement; final cap unchanged.',
        'manual_rank_selection':False,'app_default_enabled':False,'MD_or_parameter_fitting':False})
    code=child([str(PROMOD/'bin/python'),'-I',str(implementation/'promod_candidate_worker.py'),str(pool/'input.json'),str(pool/'proposal.json')],pool/'process')
    if code:
        detail=json.loads((pool/'proposal.json').read_text())
        status='unavailable'if code==2 and detail.get('status')=='unavailable'else'failed'
        return {'status':status,'stage':'bounded_native_pool','returncode':code,'detail':detail}
    selection={'snapshot':str(snapshot.manifest_path.parent),'snapshot_sha256':snapshot_hash,'pool':str(pool/'proposal.json'),
        'pool_sha256':sha(pool/'proposal.json'),'native_input':str(pool/'input.json'),'native_input_sha256':sha(pool/'input.json')}
    write(out/'selection-input.json',selection)
    code=child([str(PYTHON),'-B',str(implementation/'select_fragment_pool.py'),str(out/'selection-input.json'),str(out/'selection.json')],out/'selection-process')
    if code:return {'status':'unavailable'if code==2 else'failed','stage':'coherent_seed_selection','returncode':code}
    selected=json.loads((out/'selection.json').read_text())
    search=out/'search';search.mkdir();attempt=search/'qualified-fragment';attempt.mkdir()
    native_input['research_expected_selected']={k:v['candidate']for k,v in selected['selections'].items()}
    write(attempt/'native-input.json',native_input)
    code=child([str(PROMOD/'bin/python'),'-I',str(implementation/'promod_candidate_worker.py'),str(attempt/'native-input.json'),str(attempt/'native-proposal.json')],attempt/'generation-process')
    if code:
        detail=json.loads((attempt/'native-proposal.json').read_text())
        status='unavailable'if code==2 and detail.get('status')=='unavailable'else'failed'
        return {'status':status,'stage':'selected_native_reconstruction','returncode':code,'detail':detail}
    placement={'snapshot':str(snapshot.manifest_path.parent),'snapshot_sha256':snapshot_hash,
        'proposal':str(attempt/'native-proposal.json'),'proposal_sha256':sha(attempt/'native-proposal.json')}
    write(attempt/'placement-input.json',placement)
    code=child([str(PYTHON),'-B',str(implementation/'place_fragment_proposal.py'),str(attempt/'placement-input.json'),str(attempt/'proposal.json')],attempt/'placement-process')
    if code:return {'status':'failed','stage':'coherent_placement','returncode':code}
    validation={'snapshot':str(snapshot.manifest_path.parent),'snapshot_sha256':snapshot_hash,'proposal':str(attempt/'proposal.json'),
        'proposal_sha256':sha(attempt/'proposal.json'),'output':str(attempt/'complete'),'preserve_seed_torsions':True}
    write(attempt/'validation-input.json',validation)
    code=child([str(PYTHON),'-B',str(implementation/'validate_snapshot.py'),str(attempt/'validation-input.json'),str(attempt/'validation-result.json')],attempt/'validation-process',300)
    decision=json.loads((attempt/'validation-result.json').read_text())
    if code not in (0,2) or decision.get('status')not in {'static_checks_passed','rejected'}:
        return {'status':'failed','stage':'complete_validation','returncode':code,'decision':decision}
    accepted=code==0 and decision['status']=='static_checks_passed'
    write(search/'search.json',{'status':'accepted'if accepted else'rejected','accepted_candidate':'qualified-fragment'if accepted else None,
        'method':'Automatic bounded native catalog qualification, coherent placement and complete validation.',
        'full_refinements':1,'decision':decision})
    return {'status':'static_checks_passed'if accepted else'rejected','decision':decision,
            'manual_rank_selection':False,'app_default_enabled':False}


if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('snapshot');parser.add_argument('snapshot_sha256');parser.add_argument('output')
    parser.add_argument('--seed',type=int,default=42);args=parser.parse_args()
    try:result=run(args.snapshot,args.snapshot_sha256,args.output,args.seed)
    except Exception as exc:
        result={'status':'failed','error_type':type(exc).__name__,'error':str(exc),'app_default_enabled':False}
        if Path(args.output).is_dir():write(Path(args.output)/'result.json',result)
        raise
    write(Path(args.output)/'result.json',result)
    print(json.dumps({'status':result['status'],'result':str(Path(args.output)/'result.json')}))
    raise SystemExit(0 if result['status']=='static_checks_passed'else 2)
