"""Immutable cold conventional-HF ESP / native RESP charge-evidence branch.

All public workflow entry points are read-only and never spawn processes. A
validated numerical charge vector remains unaccepted as a physical force field.
The old checkpoint-only continuation and native handoff are unchanged.
"""
import ast
import hashlib
import json
import math
import re
from pathlib import Path
import types

CORE_SHA='e160aee512211df052261b1ef11b4277e81fd192cf71369dbd5c5f466674e187'
CORE_PATH=Path(__file__).with_name('recovery_cold_evidence_core.py')
source=CORE_PATH.read_bytes()
if hashlib.sha256(source).hexdigest()!=CORE_SHA:raise ValueError('Unreviewed cold-evidence core')
core=types.ModuleType('pinned_cold_evidence_core');core.__file__=str(CORE_PATH)
exec(compile(source,str(CORE_PATH),'exec'),core.__dict__)
require=core.require;data=core.data;sha=core.sha;load=core.load;digest=core.digest


def worker_module(context):
    module=context['materializer'];path=context['source']/'qm_worker.py'
    namespace=types.ModuleType('pinned_qm_validation');namespace.__file__=str(path)
    exec(compile(data(path,module.WORKER_SHA),str(path),'exec'),namespace.__dict__)
    return namespace


def read_native_esp(folder,materialized,context):
    """Immutable native-format reader. Caller must first admit launch receipt."""
    import numpy as np
    module=context['materializer'];folder=Path(folder)
    files=core.snapshot(folder,core.ESP_NATIVE_FILES);pins=core.hashes(files)
    result=json.loads(files['result.json']);request=json.loads(files['input.json'])
    require(files['input.json']==materialized['files']['exact-esp-input.json'],
            'Native ESP captured a different input file')
    require(digest(request)==digest(materialized['request']),'Native ESP request differs from its materialized parent')
    worker=worker_module(context);resolved=worker.validate_input(request)
    require(digest(json.loads(files['resolved-input.json']))==digest(resolved),'Resolved ESP input changed')
    n=len(context['request']['atom_ids']);points=materialized['request']['esp_points_bohr'];count=len(points)
    arrays=core.real_arrays(files['arrays.npz'],{'coords_bohr':(n,3),'coords_angstrom':(n,3),
        'atomic_numbers':(n,),'effective_nuclear_charges':(n,),'esp_points_bohr':(count,3),'esp_hartree_per_e':(count,)})
    binding={'input_sha256':pins['input.json'],'worker_sha256':module.WORKER_SHA,'runtime_sha256':module.RUNTIME_SHA,
        'arrays_sha256':pins['arrays.npz'],'resolved_method_file_sha256':pins['resolved-method.json']}
    module.validate_cold_esp(request,result,context['request'],context['method'],context['arrays']['coords_bohr'].tolist(),
        points,{k:v.tolist() for k,v in arrays.items()},module.COLD,binding)
    require(result.get('phase')=='completed' and digest(json.loads(files['progress.json']))==
            digest({'phase':'completed','status':'completed','accepted':True}),'Native ESP is not terminal')
    require(pins['runtime-manifest.json']==module.RUNTIME_SHA and result.get('checkpoint_file')=='scf.chk' and
            result.get('checkpoint_sha256')==pins['scf.chk'],'Native ESP runtime/checkpoint artifact mismatch')
    # scf.chk is checked only as immutable output bytes. No density/geometry
    # consistency is asserted and this reader never imports or reuses it.
    method=json.loads(files['resolved-method.json'])
    expected_method={k:v for k,v in context['method'].items() if k not in ('auxbasis_resolved','auxbasis_sha256')}
    expected_method['density_fit']=False
    require(digest(method)==digest(expected_method),'Complete resolved conventional ESP method differs from its parent basis/state')
    require(method.get('density_fit') is False and not method.get('ecp_resolved') and
            'auxbasis_resolved' not in method and 'auxbasis_sha256' not in method,
            'Native ESP is not conventional all-electron HF')
    for field in ('basis','ecp'):
        require(method.get(field+'_sha256')==digest(method.get(field+'_resolved')) and
                method[field+'_sha256']==context['method'][field+'_sha256'],'Expanded native ESP method differs')
    require(digest(result['method'])==digest({k:v for k,v in method.items() if not k.endswith('_resolved')}),
            'Expanded/native ESP method reports disagree')
    require(digest(result.get('optimization'))==digest({'requested':False,'converged':None,'settings':resolved['optimization'],'steps':[]}),
            'ESP unexpectedly performed optimization')
    require(digest(result.get('arrays'))==digest({k:{'shape':list(v.shape),'dtype':str(v.dtype)} for k,v in arrays.items()}),
            'Native array metadata differs from actual stored arrays')
    for key in ('explicit_electron_count','nao'):
        require(type(result.get(key)) is int and result[key]==context['result'][key],
                'Native ESP electron/AO inventory differs: '+key)
    electrons=result.get('density_electron_count')
    require(type(electrons) in (int,float) and math.isfinite(electrons) and
            abs(electrons-result['explicit_electron_count'])<=1e-6,'Native ESP density has the wrong electron count')
    for key in ('atomic_numbers','effective_nuclear_charges'):
        require(np.array_equal(arrays[key],context['arrays'][key]),'Native ESP nuclear identity differs: '+key)
    require(result['units'].get('coords_angstrom')=='angstrom' and result['units'].get('effective_nuclear_charges')=='e' and
            np.allclose(arrays['coords_angstrom'],arrays['coords_bohr']*.52917721092,rtol=0,atol=1e-10),
            'Native ESP coordinate/charge unit metadata differs')
    for report_ in (result,method):
        for key in ('requested_threads','actual_pyscf_threads'):
            require(type(report_.get(key)) is int and report_[key]==request['threads'],'Native ESP thread report changed')
    require(digest(result.get('resources'))==digest({'threads':request['threads'],'max_memory_mb':request['max_memory_mb'],
        'max_wall_seconds':request['max_wall_seconds'],'requested_threads':request['threads'],'actual_pyscf_threads':request['threads'],
        'memory_limit_kind':'PySCF allocation setting, not an operating-system hard limit'}),
        'Native ESP resource metadata changed')
    check=result.get('esp',{})
    error=check.get('max_absolute_error_hartree_per_e')
    expected_indices=np.unique(np.linspace(0,count-1,min(resolved['esp_crosscheck_points'],count),dtype=int)).tolist()
    batch=min(resolved['esp_batch_size'],max(1,64_000_000//(8*result['nao']**2)))
    require(check.get('passed') is True and check.get('tolerance_hartree_per_e')==1e-8 and
            type(error) in (int,float) and math.isfinite(error) and 0<=error<=1e-8 and
            digest(check.get('checked_point_indices'))==digest(expected_indices) and
            type(check.get('batch_size')) is int and check['batch_size']==batch and
            check.get('method')=='Batched int3c2e/fakemol_for_charges; independent int1e_rinv',
            'Native ESP independent integral crosscheck changed or failed')
    require(core.hashes(core.snapshot(folder,core.ESP_NATIVE_FILES))==pins,'Native ESP evidence changed during read')
    return {'folder':folder,'hashes':pins,'request':request,'result':result,'resolved_method':method,'arrays':arrays,
            'checkpoint_reused':False,'physical_acceptance':False}


def native_esp_text(coords_bohr,points_bohr,potentials):
    """Original native RESP formatting, never an alternate electrostatic model."""
    lines=[f'{len(coords_bohr):5d}{len(points_bohr):5d}{0:5d}']
    lines+=[' '*17+''.join(f'{float(x):16.7E}' for x in row) for row in coords_bohr]
    lines+=[' '+''.join(f'{float(x):16.7E}' for x in [potential,*point])
            for point,potential in zip(points_bohr,potentials,strict=True)]
    return ('\n'.join(lines)+'\n').encode()


def check_native_esp_roundtrip(blob,esp):
    """Require the exact original writer representation and its recorded precision."""
    import numpy as np
    arrays=esp['arrays'];xyz=arrays['coords_bohr'];points=arrays['esp_points_bohr'];potentials=arrays['esp_hartree_per_e']
    require(blob==native_esp_text(xyz,points,potentials),'Native RESP ESP file differs from exact formatting of admitted ESP')
    lines=blob.decode().splitlines();n=len(xyz)
    require([int(x) for x in lines[0].split()]==[n,len(points),0],'Native RESP ESP inventory differs')
    nuclei=np.asarray([[float(x) for x in row.split()] for row in lines[1:1+n]])
    samples=np.asarray([[float(x) for x in row.split()] for row in lines[1+n:]])
    require(nuclei.shape==xyz.shape and samples.shape==(len(points),4) and
            np.allclose(nuclei,xyz,rtol=5.1e-8,atol=5.1e-8) and
            np.allclose(samples[:,0],potentials,rtol=5.1e-8,atol=5.1e-12) and
            np.allclose(samples[:,1:],points,rtol=5.1e-8,atol=5.1e-8),'Native RESP ESP precision roundtrip failed')


def original_charge_validator(context):
    path=context['source']/'covalent_charge_handoff.py'
    return context['materializer'].extract_function(path,'validate_charges',{'require':require,'math':math},core.HANDOFF_SHA)


LAUNCHER_SHA='66127b114c0f795a892d0d0ad0b4ff4bbb1c8961f934279363a2d235d67f5980'
LAUNCHER=Path('/Users/ashujo/.cache/dynamol-research/6oim-resp-recovery-v1/prepared/cold-esp-launch-v1/launch_cold_esp.py')
LAUNCH_AUDIT=Path('/Users/ashujo/.cache/dynamol-research/reviews/6oim-cold-esp-launch-v3')
LAUNCH_MANIFEST_SHA='2ab394c6176fbbec757142d74e168eb41dfd99a6415e13556630510c6fa9308c'
LAUNCH_REVIEW=Path('/Users/ashujo/.cache/dynamol-research/reviews/6oim-cold-esp-launch-independent-final-v1/review.json')
LAUNCH_REVIEW_SHA='6e7f4c72b8c1d3fe285133675d223897beb3e8d691c908e4ff2cd13588a8f2f5'
SUPERVISOR_SHA='e3f1ac0b86d2d9d13e2796ac564439eea76a93a39e554033f93b5a36cb1eec18'
QM_INVENTORY_SHA='683b958788605f35ce2fc9a5b7ac94edb03f5316bd0b166cfcb415ffaea1f416'
LAUNCH_LIMITS=dict(native_memory_mb=2000,threads=2,max_process_rss_bytes=8_000_000_000,
    max_group_rss_bytes=8_000_000_000,native_wall_seconds=14400,outer_wall_seconds=14520,
    entry_disk_bytes=35*1024**3,running_disk_bytes=8*1024**3)
LAUNCH_FILES={'controller-source.py','launch-spec.json','input.json','qm_worker.py',
              'controller-native.log','controller-progress.json','controller-result.json'}


def launcher_admission():
    """Bind separate review and enablement; never import or execute the launcher."""
    data(LAUNCHER,LAUNCHER_SHA);review=load(LAUNCH_REVIEW,LAUNCH_REVIEW_SHA)
    manifest=load(LAUNCH_AUDIT/'artifact-manifest.json',LAUNCH_MANIFEST_SHA)
    for name,pin in manifest.items():
        require(Path(name).name==name,'Invalid launch-audit artifact name')
        blob=data(LAUNCH_AUDIT/name,pin['sha256']);require(len(blob)==pin['bytes'],'Launch-audit size changed')
    enablement=load(LAUNCH_AUDIT/'enablement-review.json')
    require(enablement['source_sha256']==LAUNCHER_SHA and enablement['all_functions_AST_exact'] is True and
            enablement['native_quantum_launched'] is False and review['core_sha256']==CORE_SHA,
            'Launcher enablement/review binding changed')
    return {'controller_sha256':LAUNCHER_SHA,'review_sha256':LAUNCH_REVIEW_SHA,
            'enablement_sha256':manifest['enablement-review.json']['sha256'],'manifest_sha256':LAUNCH_MANIFEST_SHA}


def expected_launch_spec(materialized,context):
    module=context['materializer']
    return dict(schema_version=1,parent=str(module.PARENT),plan_sha256=core.MATERIALIZER_PLAN_SHA,
        core_sha256=CORE_SHA,materializer_sha256=core.MATERIALIZER_SHA,materialized_hashes=materialized['hashes'],
        input_sha256=materialized['hashes']['exact-esp-input.json'],worker_sha256=module.WORKER_SHA,
        controller_sha256=LAUNCHER_SHA,runtime_manifest_sha256=module.RUNTIME_SHA,
        runtime_inventory_sha256=QM_INVENTORY_SHA,supervisor_sha256=SUPERVISOR_SHA,limits=dict(LAUNCH_LIMITS))


def validate_launch_receipt(receipt,progress,spec,spec_hash,native_hashes):
    """Pure terminal-receipt check. Hashes are observed from immutable files."""
    require(digest(receipt)==digest(progress),'Final launch progress and result differ')
    require(type(receipt.get('schema_version')) is int and receipt['schema_version']==1 and
        receipt.get('status')=='numerically_complete_pending_esp_review' and
        'child_pid' in receipt and receipt['child_pid'] is None and
        receipt.get('physical_acceptance') is False and receipt.get('simulation_ready') is False and
        type(receipt.get('exit_code')) is int and receipt['exit_code']==0 and
        receipt.get('cleanup_state')=='verified_no_live_group','ESP launch is not cleanly terminal')
    expected={k:spec[k] for k in ('core_sha256','materializer_sha256','controller_sha256','worker_sha256',
        'runtime_manifest_sha256','runtime_inventory_sha256','limits','materialized_hashes','input_sha256')}
    expected.update(spec_sha256=spec_hash,native_hashes=native_hashes,result_sha256=native_hashes['result.json'])
    require(digest({k:receipt.get(k) for k in expected})==digest(expected),'ESP terminal source/input/output binding differs')
    require(type(receipt.get('controller_pid')) is int and receipt['controller_pid']>0,'Invalid controller process identity')
    times=[receipt.get(k) for k in ('started_unix','qm_started_unix','completed_unix','updated_unix')]
    require(all(type(t) in (int,float) and math.isfinite(t) and t>0 for t in times) and
            times==sorted(times),'Invalid launch time ordering')
    for field,bound in [('peak_sampled_process_rss_bytes','max_process_rss_bytes'),
                        ('peak_sampled_group_rss_bytes','max_group_rss_bytes')]:
        value=receipt.get(field)
        require(type(value) is int and 0<=value<=LAUNCH_LIMITS[bound],'Invalid/over-limit sampled launch resources')
    sample=receipt.get('resource_sample')
    if sample is not None:
        for key,peak in [('maximum_process_rss_bytes','peak_sampled_process_rss_bytes'),
                         ('process_group_rss_bytes','peak_sampled_group_rss_bytes')]:
            require(type(sample.get(key)) is int and 0<=sample[key]<=receipt[peak],'Resource sample contradicts sampled peak')
        require(type(sample.get('free_disk_bytes')) is int and sample['free_disk_bytes']>=LAUNCH_LIMITS['running_disk_bytes'] and
            type(sample.get('sample_unix')) in (int,float) and math.isfinite(sample['sample_unix']) and
            times[1]<=sample['sample_unix']<=times[2],'Invalid terminal resource sample')
    allowed=set(expected)|{'schema_version','status','child_pid','controller_pid','physical_acceptance','simulation_ready',
        'exit_code','cleanup_state','started_unix','qm_started_unix','completed_unix','updated_unix',
        'peak_sampled_process_rss_bytes','peak_sampled_group_rss_bytes','resource_sample'}
    require(set(receipt)<=allowed,'Unexpected or contradictory terminal receipt fields')


def read_launch_receipt(job,materialized,context):
    job=Path(job);require(job.is_dir() and not job.is_symlink() and job.resolve().is_relative_to(context['materializer'].CACHE),
        'ESP launch evidence must be in stable recovery storage')
    require({p.name for p in job.iterdir()}==LAUNCH_FILES|{'materialized','native','scratch'},'ESP launch inventory differs')
    for name in ('materialized','native','scratch'):
        require((job/name).is_dir() and not (job/name).is_symlink(),'Invalid ESP launch directory')
    files={name:data(job/name) for name in LAUNCH_FILES};pins=core.hashes(files)
    require(pins['controller-source.py']==LAUNCHER_SHA and pins['qm_worker.py']==context['materializer'].WORKER_SHA and
            files['input.json']==materialized['files']['exact-esp-input.json'],'Captured launch source or input changed')
    spec=json.loads(files['launch-spec.json']);expected=expected_launch_spec(materialized,context)
    require(digest(spec)==digest(expected),'Complete launch specification differs')
    native_hashes=core.hashes(core.snapshot(job/'native',core.ESP_NATIVE_FILES))
    receipt=json.loads(files['controller-result.json'])
    validate_launch_receipt(receipt,json.loads(files['controller-progress.json']),spec,pins['launch-spec.json'],native_hashes)
    require(core.hashes({name:data(job/name) for name in LAUNCH_FILES})==pins,'Launch files changed during read')
    return {'folder':job,'hashes':pins,'native_hashes':native_hashes,'receipt':receipt,'spec':spec}


def read_completed_esp(plan_path,job):
    """Actual-parent-only entry point; no output writes, launches or checkpoint reuse."""
    context=core.revalidate_actual_parent(plan_path)  # Refuse unfinished parent before downstream reads.
    admission=launcher_admission();job=Path(job)
    materialized=core.verify_materialized_request(job/'materialized',context)
    launch=read_launch_receipt(job,materialized,context)
    esp=read_native_esp(job/'native',materialized,context)
    require(esp['hashes']==launch['native_hashes'],'Native ESP changed after terminal receipt admission')
    final_context=core.revalidate_actual_parent(plan_path)
    final_materialized=core.verify_materialized_request(job/'materialized',final_context)
    final_launch=read_launch_receipt(job,final_materialized,final_context)
    require(final_launch['hashes']==launch['hashes'] and final_launch['native_hashes']==esp['hashes'],
            'Launch/parent evidence changed during full read')
    return {'status':'bound_numerical_esp_evidence','accepted':False,'physical_acceptance':False,'simulation_ready':False,
        'checkpoint_reused':False,'context':context,'materialized':materialized,'launch':launch,'esp':esp,
        'launcher_admission':admission,'scope':'Converged exact-parent conventional-HF ESP evidence; no fitted charges or force field.'}


def parse_charge_text(blob,count):
    values=[float(x) for x in blob.decode('ascii').split()]
    require(len(values)==count and all(math.isfinite(x) for x in values),'Missing or nonfinite native RESP charges')
    return values


def parse_native_resp_stage(files,stage_input,initial_charge_bytes,atomic_numbers,esp_bytes):
    """Read completed native-format bytes only, without claiming their launch origin.

    This lower-level parser cannot admit a fit: the stable-cache executor and its
    reviewed terminal receipts are a separate remaining obligation.
    """
    esp_lines=esp_bytes.decode('ascii').splitlines();esp_header=[int(x) for x in esp_lines[0].split()]
    count=len(atomic_numbers)
    require(len(esp_header)==3 and esp_header[0]==count and esp_header[1]>count and esp_header[2]==0 and
            len(esp_lines)==1+count+esp_header[1],'RESP stage ESP inventory differs')
    for i,line in enumerate(esp_lines[1:]):
        row=line.split();require(len(row)==(3 if i<count else 4) and all(math.isfinite(float(x)) for x in row),
                               'RESP stage ESP coordinate/potential is invalid')
    require(set(files)=={'out','punch','qout','esout','log'},'RESP stage output inventory differs')
    require(all(type(v) is bytes and (len(v)>0 or k=='log') and len(v)<=core.MAX_ARTIFACT_BYTES for k,v in files.items()),'Missing native stage output')
    values=parse_charge_text(files['qout'],count);initial=parse_charge_text(initial_charge_bytes,count)
    lines=stage_input.decode('ascii').splitlines()
    require([int(x) for x in lines[4].split()][1:]==[count] and len(lines)>=5+count,'RESP stage input inventory differs')
    ivary=[]
    for i,line in enumerate(lines[5:5+count]):
        pair=[int(x) for x in line.split()];require(pair[0]==int(atomic_numbers[i]) and len(pair)==2,'RESP input nuclear order differs')
        ivary.append(pair[1])
    text=files['punch'].decode('ascii')
    controls={}
    for key in ('iqopt','irstrnt','ihfree','qwt'):
        matches=re.findall(r'\b'+key+r'\s*=\s*([-+0-9.Ee]+)',lines[1])
        require(len(matches)==1,'Native RESP input control missing or duplicated')
        controls[key]=float(matches[0])
        require(math.isfinite(controls[key]),'Native RESP input control is nonfinite')
    header_matches=re.findall(r'iqopt\s+irstrnt\s+ihfree\s+qwt\s*\n([^\n]+)',text)
    require(len(header_matches)==1,'Native RESP reported fit controls missing or duplicated')
    reported=[float(x) for x in header_matches[0].split()]
    require(reported==[controls[k] for k in ('iqopt','irstrnt','ihfree','qwt')],
            'Native RESP reported hyperbolic-restraint controls differ')
    marker='Point charges before & after optimization'
    require(text.count(marker)==1,'Native RESP point-charge table missing or duplicated')
    table=text.split(marker)[1].split('Statistics of the fitting:')[0].splitlines()[2:]
    rows=[line.split() for line in table if line.strip()]
    require(len(rows)==count,'Native RESP point-charge table is incomplete')
    for i,row in enumerate(rows):
        require(len(row)==6 and int(row[0])==i+1 and int(row[1])==int(atomic_numbers[i]) and int(row[4])==ivary[i],
                'Native RESP point-charge identity/order/constraint differs')
        numbers=[float(row[j]) for j in (2,3,5)]
        require(all(math.isfinite(x) for x in numbers) and abs(numbers[0]-initial[i])<=5.1e-7 and
                abs(numbers[1]-values[i])<=5.1e-7,'Native RESP punch and native charge vector differ')
        if ivary[i]==-1:require(abs(values[i]-initial[i])<=1e-7,'Native RESP frozen stage charge changed')
        elif ivary[i]>0:
            require(ivary[i]<=count and abs(values[i]-values[ivary[i]-1])<=1e-7,'Native RESP stage equivalence changed')
    out=files['out'].decode('ascii');log=files['log'].decode('ascii')
    require(not re.search(r'Error on OPEN|STOP|NaN|Infinity|convergence failure',out+'\n'+log,re.I),
            'Native RESP reports an error despite output files')
    metrics={}
    for label,key in [('The initial sum of squares (ssvpot)','initial_sum_squares'),
        ('The residual sum of squares (chipot)','residual_sum_squares'),
        ('The std err of estimate (sqrt(chipot/N))','standard_error'),
        ('ESP relative RMS (SQRT(chipot/ssvpot))','relative_rms')]:
        matches=re.findall(re.escape(label)+r'\s+([-+0-9.Ee]+)',out)
        require(len(matches)==1,'Native RESP completed fit statistics missing or duplicated')
        value=float(matches[0]);require(math.isfinite(value) and value>=0,'Native RESP fit statistic is invalid')
        metrics[key]=value
    eslines=files['esout'].decode('ascii').splitlines();header=eslines[0].split()
    require(len(header)==4 and int(header[0])>0 and int(header[1])==6 and len(eslines)==int(header[0])+1 and int(header[0])==esp_header[1],
            'Native RESP ESP diagnostic inventory is incomplete')
    require(all(math.isfinite(float(x)) for x in header[2:]),'Native RESP ESP diagnostic header is nonfinite')
    for line in eslines[1:]:
        row=line.split();require(len(row)==6 and all(math.isfinite(float(x)) for x in row),
            'Native RESP ESP diagnostic row is invalid')
    # Diagnostic format/finite checks only; no unreviewed unit conversion or fit-quality threshold.

    return {'charges':values,'ivary':ivary,'printed_fit_statistics':metrics,
            'output_sha256':{k:hashlib.sha256(v).hexdigest() for k,v in files.items()},
            'input_sha256':hashlib.sha256(stage_input).hexdigest(),'initial_charge_sha256':hashlib.sha256(initial_charge_bytes).hexdigest(),
            'esp_sha256':hashlib.sha256(esp_bytes).hexdigest(),
            'execution_provenance_admitted':False,'physical_acceptance':False}


def validate_joint_charge_candidate(context,esp,stage1,stage2,esp_bytes,source, constraints):
    """Retain exact jointly fitted charges; this is not the execution admission gate."""
    module=context['materializer'];root=context['source'];pins=context['plan']['source_pins']
    for name in ('capped-adduct.json','resp/constraints.json','resp/stage1.in','resp/stage2.in','resp/canonical.qin'):
        data(root/'cap'/name,pins['cap/'+name])
    expected_source=load(root/'cap/capped-adduct.json');expected_constraints=load(root/'cap/resp/constraints.json')
    require(digest(source)==digest(expected_source) and digest(constraints)==digest(expected_constraints),
            'Complete joint charge graph or canonical constraints differ')
    require(stage1['input_sha256']==pins['cap/resp/stage1.in'] and stage2['input_sha256']==pins['cap/resp/stage2.in'] and
            stage1['initial_charge_sha256']==pins['cap/resp/canonical.qin'] and
            stage2['initial_charge_sha256']==stage1['output_sha256']['qout'],'Native RESP two-stage input chain differs')
    require(stage1['esp_sha256']==stage2['esp_sha256']==hashlib.sha256(esp_bytes).hexdigest(),
            'Native RESP stages do not share the admitted exact ESP file')
    check_native_esp_roundtrip(esp_bytes,esp)
    charge_ids=[f'{r["index"]:03d}:{r["name"]}' for r in source['atom_map']]
    require(digest(charge_ids)==digest(context['request']['atom_ids']),'Charge source atom identity differs from parent')
    charges=stage2['charges']
    # Exactly the old pure charge gate, with no zero-fill, fragment refit, or renormalization.
    fit={'charges':charges,'constraint_checks_passed':True}
    values=original_charge_validator(context)(source,constraints,expected_constraints,fit,charges)
    for relative,expected in constraints['source_sha256'].items():
        data(Path('/Users/ashujo/Documents/Science/DynaMol')/relative,expected)
    retained=[{'atom_id':atom_id,'charge_e':charge} for i,(atom_id,charge) in enumerate(zip(values['atom_ids'],values['charges_e'],strict=True))
              if i not in values['removed_cap_indices']]
    return {'status':'parsed_joint_charge_candidate_pending_execution_receipts','accepted':False,
        'physical_acceptance':False,'simulation_ready':False,'execution_provenance_admitted':False,
        'charge_validation':values,'retained_atoms':retained,'source_graph_sha256':digest(source),
        'scope':'Exact native output parsing and unchanged canonical charge constraints only; no launch provenance or mechanical model acceptance.'}
