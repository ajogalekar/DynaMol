"""No-launch materializer for one reviewed 6OIM recovery and explicit cold HF ESP.

Old continuation, checkpoint, charge-handoff and force-field artifacts are never
modified. A materialized request is not an accepted result. Future ESP/RESP and
native assembly require separate launch/admission; this module never runs them.
"""
import argparse
import ast
import copy
import hashlib
import json
import math
import os
from pathlib import Path
import sys
import tempfile

CACHE=Path('/Users/ashujo/.cache/dynamol-research/6oim-resp-recovery-v1')
PREP_PREFIX=Path('/Users/ashujo/.cache/dynamol-runtimes/covalent-prep-local-v1')
PARENT=Path('/Users/ashujo/.cache/dynamol-research/6oim-local-recovery-v2-memory2000/run')
SPEC_SHA='3a2077e216b79b636afa6eac354ebf75e754d0b9eb3da7582c0dece1aaf6918a'
INPUT_SHA='bf48d681b57232bf0eb416ca5f4bde948902971561850716249a1b77db28154e'
WORKER_SHA='ddb67a1f8e9be6ed05c3c95e9864e20b79f0233604e27926e30e5ffc09295bc1'
CONTROLLER_SHA='e3f1ac0b86d2d9d13e2796ac564439eea76a93a39e554033f93b5a36cb1eec18'
RUNTIME_SHA='91a060defde8a3506861dffbb5277a4672bbfd041208c2cbbfeedb66dfe0c1b9'
CHECKER_SHA='e3810f144cb85035f391262bd6a0520f2450b235828dc2d171a9ddead36058c3'
SOURCE_MANIFEST_SHA='0a5954d5d23cc8878fb8a1912433e8c8f9199e0cb7497e68cabeae3d0c6b73b5'
PREP_QUALIFICATION_SHA='a4bf3f2f4e296f7333596673c23deae7755733460c95489c2ccbce0d1368617e'
PREP_ARTIFACT_MANIFEST_SHA='cc5de227f9526b77c3107b38301942473bd287d0a208a23d4ccb9707b11afaf3'
PREP_INDEPENDENT_REVIEW_SHA='e4913028d6c8692ec069fe79e33d82655d12b8816e1d8be545e7ea17cfc5e259'
RESP_PREFIX=Path('/Users/ashujo/.cache/dynamol-runtimes/resp-heap-local-v1')
RESP_PINS={'qualification_sha256':'dfe338698f80029eb15668cb44480b695e875ac0b85f37794b0e422448889851',
           'inventory_sha256':'80a80b1e6716aade0a42150957ee9056e56a49a16fa398e93d93212f167fb0c7',
           'artifact_manifest_sha256':'f519f2266c1de508994292303192b14c612b9610e079f72d894879fcae7cf189',
           'independent_review_sha256':'92d59bc4bcdb50586aa3bcc39ca9f4cae6d8aa53f3f104be2f1514098de012a5',
           'executable_sha256':'2d759392c83363f4294c01bc7c8602134e51b3a58f6b350044e1becce61048a3'}
METHOD_PINS={'basis_sha256':'da5c390b26164095e90eb5ae8693bfe2e4b5fd88ec6682a8d3079d81da3fb5ed',
             'ecp_sha256':'44136fa355b3678a1146ad16f7e8649e94fb4fc21fe77e8310c060f61caaff8a',
             'auxbasis_sha256':'8f3cb962eb094db5d454011cdf8a71fbc9649c972222cda7e4ce8094bc344b1c'}
COLD='cold_scf_no_checkpoint'
NATIVE_FILES={'input.json','resolved-input.json','native.log','runtime-manifest.json','progress.json',
              'constraints.txt','optimization-report.json','resolved-method.json','arrays.npz','scf.chk','result.json'}


def require(condition,message):
    if not condition:raise ValueError(message)


def raw(path,expected=None):
    path=Path(path)
    require(path.is_file() and not path.is_symlink(),'Missing/nonregular source: '+str(path))
    data=path.read_bytes()
    require(len(data)==path.stat().st_size and (data or expected==hashlib.sha256(b'').hexdigest()),
            'Incomplete source read: '+str(path))
    if expected is not None:require(hashlib.sha256(data).hexdigest()==expected,'Source hash mismatch: '+str(path))
    return data


def load(path,expected=None):return json.loads(raw(path,expected))
def sha(path):return hashlib.sha256(raw(path)).hexdigest()
def digest(value):return hashlib.sha256(json.dumps(value,sort_keys=True,separators=(',',':'),allow_nan=False).encode()).hexdigest()


def publish(path,value):
    """Atomic no-clobber JSON publication. Failed temporary files are retained."""
    data=(json.dumps(value,indent=2,allow_nan=False)+'\n').encode()
    fd,temporary=tempfile.mkstemp(prefix='.'+path.name+'.',suffix='.tmp',dir=path.parent)
    with os.fdopen(fd,'wb') as stream:
        stream.write(data);stream.flush();os.fsync(stream.fileno())
    os.link(temporary,path)
    directory=os.open(path.parent,os.O_RDONLY)
    try:os.fsync(directory)
    finally:os.close(directory)
    os.unlink(temporary)


def validate_controller(controller):
    require(type(controller.get('schema_version')) is int and controller['schema_version']==1,
            'Invalid terminal controller schema')
    require(controller.get('status')=='numerically_complete_pending_geometry_review' and
            'child_pid' in controller and controller['child_pid'] is None,'Parent is running, failed or nonterminal')
    for key in ('physical_acceptance','simulation_ready'):
        require(controller.get(key) is False,'Unexpected parent physical acceptance flag')
    require(controller.get('spec_sha256')==SPEC_SHA and controller.get('controller_sha256')==CONTROLLER_SHA,
            'Parent controller/source specification mismatch')
    for key,value in {'native_memory_mb':2000,'max_rss_bytes':8_000_000_000,'threads':2,'outer_wall_seconds':21720}.items():
        observed=controller.get('limits',{}).get(key)
        require(type(observed) is int and observed==value,'Unreviewed resource setting: '+key)
    peak=controller.get('peak_sampled_rss_bytes')
    require(type(peak) is int and 0<peak<=8_000_000_000,'Invalid parent observed memory peak')


def validate_identity(observed,expected):
    for key in ('atom_ids','elements'):
        require(observed.get(key)==expected[key],'Atom identity/order mismatch: '+key)
    for key in ('charge','spin'):
        require(type(observed.get(key)) is int and observed[key]==expected[key],'Electronic state mismatch: '+key)


def validate_preparation_runtime(admission):
    require(admission.get('root_reviewed') is True,'Preparation runtime independent admission is pending')
    require(admission.get('sha256') and admission.get('path'),'Preparation runtime qualification is pending')
    require(admission['sha256']==PREP_QUALIFICATION_SHA and
            admission.get('artifact_manifest_sha256')==PREP_ARTIFACT_MANIFEST_SHA and
            admission.get('independent_review_sha256')==PREP_INDEPENDENT_REVIEW_SHA,
            'Preparation admission differs from the independently reviewed artifacts')
    raw(admission['artifact_manifest_path'],PREP_ARTIFACT_MANIFEST_SHA)
    raw(admission['independent_review_path'],PREP_INDEPENDENT_REVIEW_SHA)
    qualification=load(admission['path'],admission['sha256'])
    require(qualification.get('status')=='qualified_for_scoped_preparation' and
            qualification.get('required_sys_prefix')==str(PREP_PREFIX) and
            Path(sys.prefix).resolve()==PREP_PREFIX,'Use the separately qualified stable preparation runtime')
    inventory=load(qualification['inventory_path'],qualification['inventory_sha256'])
    for item in inventory['files']:
        path=Path(item['path'])
        require(path.resolve().is_relative_to(PREP_PREFIX),'Runtime payload escaped the reviewed prefix')
        raw(path,item['sha256'])
    for item in inventory['symlinks']:
        path=Path(item['path'])
        require(path.is_symlink() and os.readlink(path)==item['link_target'] and str(path.resolve())==item['resolved_path'],
                'Runtime interpreter link changed')
        raw(path.resolve(),item['resolved_sha256'])
    for key in ('package_record_comparison','checks','import_observations'):
        raw(qualification[key+'_path'],qualification[key+'_sha256'])
    require(qualification.get('check_geometry_original_source_sha256')==CHECKER_SHA,'Runtime qualified a different geometry checker')
    return qualification


def validate_resp_runtime(admission):
    """Verify the separately reviewed local provider; never execute or qualify a fit."""
    require(admission.get('root_reviewed') is True and admission.get('prefix')==str(RESP_PREFIX),
            'Stable RESP runtime independent admission is pending')
    require(all(admission.get(k)==v for k,v in RESP_PINS.items()),'Unreviewed RESP provider evidence')
    require(RESP_PREFIX.resolve()==RESP_PREFIX and not RESP_PREFIX.is_symlink(),'RESP runtime path changed')
    qualification=load(admission['qualification_path'],RESP_PINS['qualification_sha256'])
    review=load(admission['independent_review_path'],RESP_PINS['independent_review_sha256'])
    require(qualification.get('qualified') is True and qualification.get('runtime_prefix')==str(RESP_PREFIX) and
            qualification.get('executable')==str(RESP_PREFIX/'bin/resp') and
            qualification.get('executable_sha256')==RESP_PINS['executable_sha256'],
            'RESP qualification does not bind this provider')
    require(review.get('status')=='no_blocking_finding_for_scoped_RESP_relocation' and
            review.get('qualification_sha256')==RESP_PINS['qualification_sha256'],
            'RESP independent review does not bind the qualification')
    # The original qualification's pending-review field is immutable. Admission
    # comes from the separately pinned review above, never by rewriting it.
    inventory=load(qualification['runtime_inventory'],RESP_PINS['inventory_sha256'])
    require({x['path'] for x in inventory['files']}==
            {'bin/resp','lib/libgcc_s.1.1.dylib','lib/libgfortran.5.dylib','lib/libquadmath.0.dylib'} and
            len(inventory['files'])==4,'Unexpected RESP runtime payload inventory')
    for item in inventory['files']:
        path=RESP_PREFIX/item['path']
        require(path.resolve().is_relative_to(RESP_PREFIX),'RESP dependency escaped its local runtime')
        require(len(raw(path,item['sha256']))==item['bytes'],'RESP payload size changed')
    artifacts=load(admission['artifact_manifest_path'],RESP_PINS['artifact_manifest_sha256'])
    require(len(artifacts['artifacts'])==91,'RESP qualification artifact inventory changed')
    artifact_root=Path(admission['artifact_manifest_path']).parent
    for item in artifacts['artifacts']:
        path=artifact_root/item['path']
        require(path.resolve().is_relative_to(artifact_root.resolve()),'RESP evidence path escaped its audit')
        require(len(raw(path,item['sha256']))==item['bytes'],'RESP evidence size changed')
    return {'provider_admitted':True,'prefix':str(RESP_PREFIX),'pins':dict(RESP_PINS),
            'fit_launched':False,'charges_accepted':False,'physical_acceptance':False,
            'required_executor':'Reviewed generic supervisor with fresh private cwd, short relative filenames, '
                '3 GB sampled RSS stop and parsed native outputs/charge checks; not the parity harness.'}


def validate_cold_request(request,parent_request,coords_bohr,points_bohr,policy):
    require(policy==COLD,'Checkpoint policy must be explicitly cold SCF')
    require('initial_checkpoint' not in request,'Hidden or explicit checkpoint reuse is forbidden in cold branch')
    expected=build_cold_request(parent_request,coords_bohr,points_bohr)
    require(request==expected,'Cold ESP request changed method, state, geometry, grid or reviewed settings')
    require(type(request.get('schema_version')) is int and request['schema_version']==1 and
            request.get('density_fit') is False,'Malformed cold request schema/method flag')
    for key in ('threads','max_memory_mb','max_wall_seconds','esp_batch_size'):
        require(type(request.get(key)) is int,'Malformed cold request resource: '+key)
    for key in ('max_scf_cycles','max_cphf_cycles','esp_crosscheck_points','grid_level'):
        if key in request:require(type(request[key]) is int,'Malformed cold request integer setting: '+key)
    validate_identity(request,parent_request)
    for name,matrix,rows in [('coordinates',coords_bohr,len(parent_request['atom_ids'])),
                             ('grid',points_bohr,len(points_bohr)),
                             ('request coordinates',request.get('coords_bohr'),len(parent_request['atom_ids'])),
                             ('request grid',request.get('esp_points_bohr'),len(points_bohr))]:
        require(type(matrix) is list and len(matrix)==rows and rows>0 and
                all(type(row) is list and len(row)==3 and all(type(x) in (int,float) and math.isfinite(x) for x in row)
                    for row in matrix),'Invalid cold ESP '+name)


def build_cold_request(parent_request,coords_bohr,points_bohr):
    """Preserve original conventional-HF ESP chemistry; resource allocation only differs."""
    request={k:copy.deepcopy(v) for k,v in parent_request.items()
             if k not in ('coords_angstrom','coords_bohr','optimization','initial_checkpoint','auxbasis')}
    request.update(coords_bohr=copy.deepcopy(coords_bohr),density_fit=False,operations=['esp'],
        esp_points_bohr=copy.deepcopy(points_bohr),esp_batch_size=64,
        max_wall_seconds=14400,max_memory_mb=2000,threads=2)
    return request


def validate_cold_esp(request,result,parent_request,parent_method,coords_bohr,points_bohr,arrays,policy,binding):
    """New cold branch only; hash binding of files is also mandatory at integration.

    binding is supplied by a separately reviewed immutable-file reader, and
    contains the observed input, worker, runtime, arrays and resolved-method
    SHA256 values. This pure validator is not that file reader or a launch path.
    It returns no charges or acceptance certificate.
    """
    validate_cold_request(request,parent_request,coords_bohr,points_bohr,policy)
    require(type(result.get('schema_version')) is int and result['schema_version']==1 and
            result.get('accepted') is True and result.get('status')=='completed' and
            result.get('scf',{}).get('converged') is True,'Conventional HF ESP did not converge')
    require('initial_checkpoint' not in result,'Cold ESP result reports checkpoint reuse')
    validate_identity({**result,'spin':result.get('spin_2S')},parent_request)
    require(set(binding)=={'input_sha256','worker_sha256','runtime_sha256','arrays_sha256','resolved_method_file_sha256'},
            'Missing or unknown cold ESP evidence binding')
    for key,value in binding.items():
        require(type(value) is str and len(value)==64 and all(c in '0123456789abcdef' for c in value),'Invalid evidence hash: '+key)
        observed=result.get('runtime_manifest',{}).get('sha256') if key=='runtime_sha256' else result.get(key)
        require(observed==value,'Cold ESP evidence mismatch: '+key)
    require(binding['worker_sha256']==WORKER_SHA and binding['runtime_sha256']==RUNTIME_SHA,'Unreviewed ESP worker/runtime')
    require(result.get('arrays_file')=='arrays.npz' and result.get('resolved_method_file')=='resolved-method.json',
            'Unexpected ESP artifact paths')
    scf=result['scf']
    require(scf.get('conv_tol_hartree')==parent_request.get('scf_conv_tol',1e-10) and
            scf.get('conv_tol_grad')==parent_request.get('scf_conv_tol_grad',1e-7) and
            type(scf.get('max_cycles')) is int and scf['max_cycles']==parent_request.get('max_scf_cycles',100) and
            type(scf.get('cycles')) is int and 0<scf['cycles']<=scf['max_cycles'],'ESP SCF settings changed')
    require(type(result.get('energy_hartree')) in (int,float) and math.isfinite(result['energy_hartree']),
            'Nonfinite ESP energy')
    units=result.get('units',{})
    require(all(units.get(k)==v for k,v in {'energy':'hartree','coords_bohr':'bohr','esp_points_bohr':'bohr',
            'esp_hartree_per_e':'hartree/e','bohr_to_angstrom':.52917721092}.items()),'ESP units changed')
    method=result.get('method',{})
    require(method.get('method')=='RHF' and method.get('density_fit') is False and
            method.get('basis_requested','').lower()=='6-31g*' and
            method.get('basis_sha256')==parent_method['basis_sha256'] and
            method.get('ecp_sha256')==parent_method['ecp_sha256'] and
            method.get('spherical_basis') is True and method.get('symmetry') is False and
            method.get('functional_requested') is None and method.get('grid_level') is None and
            method.get('ecp_requested') is None and 'auxbasis_sha256' not in method,'Conventional ESP resolved method changed')
    require(arrays.get('coords_bohr')==coords_bohr and arrays.get('esp_points_bohr')==points_bohr,
            'ESP result uses different geometry or grid')
    for key in ('coords_bohr','esp_points_bohr'):
        require(type(arrays[key]) is list and all(type(row) is list and len(row)==3 and
            all(type(x) in (int,float) and math.isfinite(x) for x in row) for row in arrays[key]),
            'Malformed ESP coordinate array: '+key)
    values=arrays.get('esp_hartree_per_e',[])
    require(len(values)==len(points_bohr) and len(values)>len(parent_request['atom_ids']) and
            all(type(x) in (int,float) and math.isfinite(x) for x in values),'Nonfinite/malformed ESP values')


def extract_function(path,name,globals_,expected):
    data=raw(path,expected);tree=ast.parse(data)
    selected=[x for x in tree.body if isinstance(x,ast.FunctionDef) and x.name==name]
    require(len(selected)==1,'Missing original reusable function: '+name)
    exec(compile(ast.Module(body=selected,type_ignores=[]),str(path),'exec'),globals_)
    return globals_[name]


def materialize(plan_path,output):
    output=Path(output).resolve()
    require(output.is_relative_to(CACHE),'All new output/scratch must remain in the stable local cache')
    output.mkdir(parents=True,exist_ok=False)
    report={'status':'checking','accepted':False,'physical_acceptance':False,'request_emitted':False,
            'QM_launched':False,'RESP_launched':False,'checkpoint_reused':False,'script_sha256':sha(__file__)}
    try:
        plan=load(plan_path);report['plan_sha256']=sha(plan_path)
        require(plan.get('checkpoint_policy')==COLD,'Missing/unsupported checkpoint branch')
        require(plan.get('parent')==str(PARENT),'Unreviewed parent path')
        # Reject an active or failed parent before native arrays/checkpoint or imports.
        require((PARENT/'controller-result.json').is_file(),'Accepted terminal parent absent; no request emitted')
        controller=load(PARENT/'controller-result.json');validate_controller(controller)
        admission=plan.get('preparation_runtime_admission',{})
        qualification=validate_preparation_runtime(admission)
        resp_provider=validate_resp_runtime(plan.get('resp_runtime_admission',{}))
        source=Path(plan['sources']).resolve()
        require(source.is_relative_to(CACHE),'Preparation sources must be in the stable cache')
        manifest=load(source/'source-manifest.json',SOURCE_MANIFEST_SHA)
        expected_pins={k:v['sha256'] for k,v in manifest.items()}
        require(plan.get('source_pins')==expected_pins,'Missing or changed frozen chemistry sources')
        for relative,value in expected_pins.items():raw(source/relative,value)
        import numpy as np
        import rdkit
        from rdkit import Chem
        require(all(Path(x.__file__).resolve().is_relative_to(PREP_PREFIX) for x in (np,rdkit,Chem)),
                'Preparation imports escaped the qualified runtime')
        native=PARENT/'optimization'
        require(not native.is_symlink() and {p.name for p in native.iterdir()}==NATIVE_FILES,
                'Missing/extra native parent artifacts')
        before={name:sha(native/name) for name in NATIVE_FILES}
        require(before['result.json']==controller.get('result_sha256'),'Terminal controller/result mismatch')
        raw(PARENT/'spec.json',SPEC_SHA);raw(PARENT/'controller-source.py',CONTROLLER_SHA)
        raw(PARENT/'qm_worker.py',WORKER_SHA);raw(PARENT/'input.json',INPUT_SHA)
        request=load(native/'input.json',INPUT_SHA);result=load(native/'result.json')
        require(type(result.get('schema_version')) is int and result['schema_version']==1 and result.get('accepted') is True and result.get('status')=='completed'
                and result.get('scf',{}).get('converged') is True,'Parent numerical result is unaccepted')
        require(result.get('input_sha256')==INPUT_SHA and result.get('worker_sha256')==WORKER_SHA,'Parent input/worker binding mismatch')
        require('initial_checkpoint' not in request and 'initial_checkpoint' not in result,'Recovery was not cold SCF')
        validate_identity({**result,'spin':result.get('spin_2S')},request)
        require(len(request['atom_ids'])==94 and request['charge']==0 and request['spin']==0,'Wrong full adduct inventory/state')
        worker_globals={'__name__':'frozen_validation','__file__':str(source/'qm_worker.py')}
        exec(compile(raw(source/'qm_worker.py',WORKER_SHA),str(source/'qm_worker.py'),'exec'),worker_globals)
        resolved=worker_globals['validate_input'](request)
        require(load(native/'resolved-input.json')==resolved,'Resolved parent input changed')
        require(result['scf'].get('max_cycles')==resolved['max_scf_cycles'] and
                result['scf'].get('conv_tol_hartree')==resolved['scf_conv_tol'] and
                result['scf'].get('conv_tol_grad')==resolved['scf_conv_tol_grad'],'Parent SCF thresholds changed')
        opt=result.get('optimization',{})
        require(opt.get('requested') is True and opt.get('converged') is True and opt.get('settings')==resolved['optimization'],
                'Optimization did not pass original requested settings')
        require(load(native/'optimization-report.json')==opt and
                load(native/'progress.json')=={'phase':'completed','status':'completed','accepted':True},'Nonterminal/inconsistent native reports')
        constraints=('$freeze\nxyz '+','.join(str(request['atom_ids'].index(x)+1) for x in request['optimization']['freeze_atom_ids'])+'\n').encode()
        require(raw(native/'constraints.txt')==constraints and opt.get('constraint_file_sha256')==hashlib.sha256(constraints).hexdigest(),
                'Actual cap constraints changed')
        require(opt.get('dihedral_constraints_satisfied') is True and opt.get('final_dihedral_constraints')==[],
                'Unexpected/unverified dihedral constraints')
        for name,key in [('arrays.npz','arrays_sha256'),('resolved-method.json','resolved_method_file_sha256'),('scf.chk','checkpoint_sha256')]:
            require(before[name]==result.get(key),'Parent artifact hash mismatch: '+name)
        require(before['runtime-manifest.json']==RUNTIME_SHA and result.get('runtime_manifest',{}).get('sha256')==RUNTIME_SHA,
                'Parent runtime mismatch')
        method=load(native/'resolved-method.json')
        require(method.get('method')=='RHF' and method.get('basis_requested','').lower()=='6-31g*' and
                method.get('density_fit') is True and not method.get('ecp_resolved') and
                request.get('auxbasis')=='def2-universal-jkfit' and
                all(method.get(k)==v for k,v in METHOD_PINS.items()),'Parent resolved method changed')
        for field in ('basis','ecp','auxbasis'):
            require(method.get(field+'_sha256')==digest(method.get(field+'_resolved')),'Parent expanded method payload hash changed')
        require(result.get('method')=={k:v for k,v in method.items() if not k.endswith('_resolved')},'Parent method reports disagree')
        with np.load(native/'arrays.npz',allow_pickle=False) as archive:
            require(set(archive.files)=={'coords_bohr','coords_angstrom','gradient_hartree_per_bohr','atomic_numbers','effective_nuclear_charges'},
                    'Unexpected parent arrays')
            arrays={k:archive[k].copy() for k in archive.files}
        for key,value in arrays.items():
            shape=(94,) if key in ('atomic_numbers','effective_nuclear_charges') else (94,3)
            require(value.shape==shape and value.dtype.kind in 'if' and np.isfinite(value).all(),'Malformed parent array: '+key)
        atomic_numbers=[Chem.GetPeriodicTable().GetAtomicNumber(e) for e in request['elements']]
        require(arrays['atomic_numbers'].tolist()==atomic_numbers and arrays['effective_nuclear_charges'].tolist()==atomic_numbers,
                'Parent nuclei or all-electron chemistry changed')
        require(type(result.get('explicit_electron_count')) is int and result['explicit_electron_count']==388 and
                type(result.get('nao')) is int and result['nao']==816,'Parent electron/AO inventory changed')
        require(type(result.get('energy_hartree')) in (int,float) and math.isfinite(result['energy_hartree']),
                'Nonfinite parent energy')
        require(result.get('units',{}).get('bohr_to_angstrom')==.52917721092 and
                result['units'].get('gradient_hartree_per_bohr')=='hartree/bohr','Parent units changed')
        check_geometry=extract_function(source/'continue_covalent_qm.py','check_geometry',
            {'json':json,'np':np,'Chem':Chem},CHECKER_SHA)
        optimized,geometry=check_geometry(source/'cap',request,result,arrays)
        reference=load(source/'initial-reference-comparison.json')
        require(reference.get('accepted') is True and reference.get('scf_converged') is True and
                reference.get('atom_ids')==request['atom_ids'] and reference.get('elements')==request['elements'],
                'Original conventional reference is not qualified')
        surface_globals={'__name__':'frozen_surface'}
        exec(compile(raw(source/'vdw_surface.py'),str(source/'vdw_surface.py'),'exec'),surface_globals)
        shells=[surface_globals['vdw_surface'](arrays['coords_angstrom'],[x.upper() for x in request['elements']],scale,1.0,{})[0]
                for scale in (1.4,1.6,1.8,2.0)]
        points=(np.concatenate(shells)/result['units']['bohr_to_angstrom']).tolist()
        coords=arrays['coords_bohr'].tolist()
        esp_request=build_cold_request(request,coords,points)
        validate_cold_request(esp_request,request,coords,points,plan['checkpoint_policy'])
        worker_globals['validate_input'](esp_request)
        require({name:sha(native/name) for name in NATIVE_FILES}==before and
                load(PARENT/'controller-result.json')==controller,'Parent evidence changed during admission')
        publish(output/'geometry.json',geometry)
        publish(output/'parent-binding.json',{'controller':controller,'native_hashes':before,
            'geometry_sha256':digest(coords),'method_sha256':digest(method),'checkpoint_geometry_validated':False,
            'checkpoint_policy':COLD,'accepted':False,'preparation_runtime_admission':admission,
            'preparation_runtime_inventory_sha256':qualification['inventory_sha256'],'source_manifest_sha256':SOURCE_MANIFEST_SHA})
        publish(output/'downstream-contract.json',{'accepted':False,'fit_inputs':{k:v for k,v in expected_pins.items() if '/resp/' in k},
            'unchanged_native_assembly_policy':load(source/'original-handoff-plan.json')['native_assembly_policy'],
            'resp_provider':resp_provider,
            'cold_branch_handoff_implemented':False,
            'requirement':'A new reviewed file-binding reader must call validate_cold_esp, then retain the original native RESP and charge/signature validation gates. The old checkpoint-only handoff is intentionally unchanged.'})
        publish(output/'exact-esp-input.json',esp_request)
        report.update(status='cold_esp_request_materialized_not_launched',request_emitted=True,
            input_sha256=sha(output/'exact-esp-input.json'),parent_result_sha256=before['result.json'],
            remaining=['Independent ESP launch/resource admission','Converged exact-geometry conventional ESP',
                'Pinned native94atom RESP with18fixedcharges','New cold-branch full native handoff review and torsion validation'])
    except Exception as error:
        report.update(status='refused',error_type=type(error).__name__,error=str(error))
    publish(output/'result.json',report)
    return report


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--plan',type=Path,required=True);parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args();outcome=materialize(args.plan,args.output)
    print(json.dumps(outcome,indent=2));raise SystemExit(0 if outcome['request_emitted'] else 1)
