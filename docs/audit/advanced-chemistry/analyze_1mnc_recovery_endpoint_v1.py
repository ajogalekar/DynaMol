"""Strict admission of the one reviewed 1MNC cycle-8 cold recovery.

No QM execution, fabricated endpoint, Hessian, parameter extraction or model
acceptance. Source-frame chemistry is delegated to the unmodified pinned original
analyzer. Only numerical endpoint input/runtime/memory pins differ in its copied
reader; all original method, state, constraint, derivative and geometry gates stay.
"""
from __future__ import annotations
import argparse
import copy
import hashlib
import json
import math
from pathlib import Path
import time
import types
import numpy as np

ROOT = Path(__file__).resolve().parents[3]
AUDIT = ROOT/'docs/audit/advanced-chemistry'
RECOVERY = Path('/Users/ashujo/.cache/dynamol-research/1mnc-local-recovery-v2-memory2000')
PRIOR = RECOVERY.parent/'1mnc-local-recovery-v1/prepared'
RUNTIME = Path('/Users/ashujo/.cache/dynamol-runtimes/qm-parallel-local-v1')
RECOVERED = Path('/Users/ashujo/Library/Application Support/DynaMol/ResearchRuns/1mnc-resource-stop-v1/cycle-008-recovered-geometry.json')
ORIGINAL_ANALYZER_HASH = '291abdc7971a8676912fcf1bda655edd5439c75fe37535633406072e128085d0'
RECOVERY_INPUT_HASH = '51f3fc98f826eb00706219e5dad76716b8d8f36b411eb4d8b84ccfcaa47b1c0c'
RECOVERY_RUNTIME_HASH = '91a060defde8a3506861dffbb5277a4672bbfd041208c2cbbfeedb66dfe0c1b9'
SPEC_HASH = 'ec77d2cf20f5e7b3f77408fb2ad9c03239c4e167f22d2c2c9f31030d1f43ee06'
CONTROLLER_HASH = 'e3f1ac0b86d2d9d13e2796ac564439eea76a93a39e554033f93b5a36cb1eec18'
PINS = {
    RECOVERY/'spec.json': SPEC_HASH,
    RECOVERY/'prepared/input.json': RECOVERY_INPUT_HASH,
    RECOVERY/'prepared/input-difference-review.json': '70b850e8a52ba3bec23028db87e9135d19071eb1df464b2afc041d27f4167734',
    PRIOR/'input.json': '60d28261a2e8d7f2cd2662d4419b2c01da04afb62b690a9936417cb85227a68c',
    PRIOR/'preparation-review.json': '6aa7298e4a622b17c817d4302811fd035f0d1da249242a5ebfd50004be5abc27',
    RECOVERED: 'eb8445aca56e40b66eb2e2d82d95dabd44ed78069762711b7128c2f3bedc1315',
    RUNTIME/'dynamol-qm-runtime.json': RECOVERY_RUNTIME_HASH,
    RUNTIME/'local-recovery-evidence/runtime-files.json': '683b958788605f35ce2fc9a5b7ac94edb03f5316bd0b166cfcb415ffaea1f416',
    RUNTIME/'local-recovery-evidence/water-qualification-v1/comparison.json': '4d801678f2bdecc067b8df4ac6ea17503c508bb0d59c97c116eeb26bf2fcc05a',
}


def original_analyzer():
    path = AUDIT/'analyze_1mnc_anchored_endpoint.py'
    data = path.read_bytes()
    if hashlib.sha256(data).hexdigest() != ORIGINAL_ANALYZER_HASH:
        raise ValueError('Original analyzer changed; a separate review is required')
    module = types.ModuleType('pinned_original_1mnc_endpoint')
    module.__file__ = str(path)
    exec(compile(data, str(path), 'exec'), module.__dict__)
    return module


core = original_analyzer()
load, raw, sha, check = core.load, core.raw, core.sha, core.check
EvidenceError, BOHR, Z, WORKER_HASH = core.EvidenceError, core.BOHR, core.Z, core.WORKER_HASH


def validate_lineage(original, recovered, request):
    """Pure comparison; successful admission describes an input, never an endpoint."""
    check(recovered.get('status') == 'completed_gradient_cycle_recovered_not_optimized_endpoint'
          and recovered.get('completed_cycle') == 8 and recovered.get('accepted') is False
          and recovered.get('optimization_converged') is False, 'Wrong recovery cycle or acceptance lineage')
    check(recovered.get('source_input_sha256') == core.INPUT_HASH, 'Recovery source input differs')
    check(recovered.get('source_native_log_sha256') == '5412ea436e880004ed15122ac2cd8d502f2a3c720bc55d4358dc274d99aa43cf',
          'Recovery native log lineage differs')
    for key in ['atom_ids','elements','charge','spin']:
        check(recovered.get(key) == original[key], 'Recovery source state differs: '+key)
    xyz = np.asarray(recovered['coords_bohr'],dtype=float)*BOHR
    check(xyz.shape == (97,3) and np.isfinite(xyz).all(), 'Invalid cycle-8 geometry')
    original_xyz = np.asarray(original['coords_angstrom'])
    anchors = original['optimization']['freeze_atom_ids']
    check(anchors == ['cap:small:15:CH3','cap:small:54:CH3','cap:small:92:CH3'], 'Original cap inventory differs')
    indices = [original['atom_ids'].index(a) for a in anchors]
    corrections = np.linalg.norm(xyz[indices]-original_xyz[indices],axis=1)
    check(float(corrections.max()) <= 6e-11, 'Anchor restoration exceeds reviewed printed-coordinate correction')
    xyz[indices] = original_xyz[indices]
    expected = copy.deepcopy(original)
    expected.pop('initial_checkpoint')
    expected.update(coords_angstrom=xyz.tolist(), max_memory_mb=2000)
    check(request == expected, 'Recovery request changed beyond cycle-8 coordinates, exact anchors, cold SCF and memory')
    return {'initial_geometry_source':'completed cycle 8, not a converged endpoint',
            'free_atom_count':94,'fixed_atom_count':3,'maximum_anchor_restore_angstrom':float(corrections.max()),
            'comparison_frame':'Original crystal cluster and omitted source protein; never restart-centered',
            'allowed_differences':['coords_angstrom: reviewed cycle 8 plus exact three original cap anchors',
                                   'max_memory_mb: 8000 to 2000','initial_checkpoint: removed for cold SCF']}


def load_lineage(ctx):
    evidence = {str(p):load(p,digest) for p,digest in PINS.items()}
    request = evidence[str(RECOVERY/'prepared/input.json')]
    lineage = validate_lineage(ctx['request'], evidence[str(RECOVERED)], request)
    prior_expected = copy.deepcopy(request); prior_expected['max_memory_mb'] = 5000
    check(evidence[str(PRIOR/'input.json')] == prior_expected, '5000 to 2000 MB lineage differs')
    spec = evidence[str(RECOVERY/'spec.json')]
    for key,value in [('input',str(RECOVERY/'prepared/input.json')),('worker',str(RECOVERY/'prepared/qm_worker.py')),
                      ('output',str(RECOVERY/'run')),('runtime',str(RUNTIME/'bin/python')),
                      ('native_memory_mb',2000),('max_rss_bytes',8_000_000_000),('threads',2),
                      ('outer_wall_seconds',21720),('entry_disk_bytes',35*1024**3),('running_disk_bytes',8*1024**3)]:
        check(spec.get(key) == value, 'Recovery execution specification differs: '+key)
    check(spec['pins'][spec['input']] == RECOVERY_INPUT_HASH and spec['pins'][spec['worker']] == WORKER_HASH,
          'Recovery spec input/worker binding differs')
    runtime = evidence[str(RUNTIME/'dynamol-qm-runtime.json')]
    old_runtime = load(ctx['base']/'source-manifest.json',core.SOURCE_MANIFEST_HASH)[str(ROOT/'.tools/qm-parallel/dynamol-qm-runtime.json')]['sha256']
    check(runtime['source_runtime']['manifest_sha256'] == old_runtime, 'Relocation source runtime differs')
    inventory = evidence[str(RUNTIME/'local-recovery-evidence/runtime-files.json')]
    check(inventory['runtime_manifest_sha256'] == RECOVERY_RUNTIME_HASH and inventory['file_count'] == 9051
          and not inventory['dataless_files'], 'Runtime inventory lineage differs')
    for row in inventory['files']:
        check(spec['pins'].get(row['path']) == row['sha256'], 'Runtime file not bound by controller: '+row['path'])
    regression = evidence[str(RUNTIME/'local-recovery-evidence/water-qualification-v1/comparison.json')]
    check(regression['status'] == 'passed' and regression['same_coordinates_exact'] is True
          and regression['same_resolved_physics_exact'] is True
          and regression['runtime_manifest_sha256'] == RECOVERY_RUNTIME_HASH
          and regression['worker_sha256'] == WORKER_HASH
          and regression['absolute_energy_difference_hartree'] <= 1e-10
          and regression['max_absolute_gradient_difference_hartree_per_bohr'] <= 1e-10,
          'Relocated-runtime small-water implementation regression not qualified')
    lineage.update(input_sha256=RECOVERY_INPUT_HASH, runtime_manifest_sha256=RECOVERY_RUNTIME_HASH,
                   original_input_sha256=core.INPUT_HASH, original_policy_sha256=core.POLICY_HASH,
                   original_analyzer_sha256=ORIGINAL_ANALYZER_HASH, specification_sha256=SPEC_HASH,
                   runtime_qualification_scope=regression['scope'], source_pins={str(p):h for p,h in PINS.items()})
    return spec, lineage


NATIVE_FILES = {'input.json','resolved-input.json','native.log','runtime-manifest.json','progress.json',
                'constraints.txt','optimization-report.json','resolved-method.json','arrays.npz','scf.chk','result.json'}


def validate_native_inventory(folder):
    folder = Path(folder)
    check(not folder.is_symlink(), 'Endpoint directory must not be a symlink')
    check({p.name for p in folder.iterdir()} == NATIVE_FILES, 'Missing or extra native endpoint artifacts')
    check(all(p.is_file() and not p.is_symlink() for p in folder.iterdir()), 'Native artifacts must be regular files')


def validate_controller(controller, native_hash):
    check(type(controller.get('schema_version')) is int and controller['schema_version'] == 1
          and controller.get('status') == 'numerically_complete_pending_geometry_review'
          and 'child_pid' in controller and controller['child_pid'] is None,
          'No terminal completed controller record')
    for key in ['physical_acceptance','simulation_ready']:
        check(controller.get(key) is False, 'Terminal controller flag is not explicit false: '+key)
    for key,expected in [('spec_sha256',SPEC_HASH),('controller_sha256',CONTROLLER_HASH),('result_sha256',native_hash)]:
        check(key in controller and controller[key] == expected, 'Terminal controller metadata differs: '+key)
    limits=controller.get('limits',{})
    for key,expected in [('native_memory_mb',2000),('max_rss_bytes',8_000_000_000),('threads',2),('outer_wall_seconds',21720)]:
        check(type(limits.get(key)) is int and limits[key] == expected, 'Terminal resource limits differ: '+key)
    peak=controller.get('peak_sampled_rss_bytes')
    check(type(peak) is int and 0 < peak <= 8_000_000_000, 'Missing or excessive observed resource peak')


def validate_real_arrays(arrays):
    # The pinned worker emits real arrays. Do not let the original generic
    # numeric reader's complex dtype acceptance reach a float coercion.
    for name,value in arrays.items():
        dtype = np.asarray(value).dtype
        check(np.issubdtype(dtype,np.integer) or np.issubdtype(dtype,np.floating),
              'Endpoint array is not real integer/floating data: '+name)


def read_native_endpoint(ctx, folder):
    folder = Path(folder)
    if not (folder/'result.json').is_file():
        return None
    result = load(folder/'result.json')
    check(result.get('schema_version') == 1 and result.get('status') == 'completed' and result.get('accepted') is True,
          'Endpoint worker did not report a completed numerical result')
    req = ctx['request']
    for k,expected in [('atom_ids',req['atom_ids']),('elements',req['elements']),('charge',1),('spin_2S',0),
                       ('explicit_electron_count',372),('nao',778),('input_sha256',RECOVERY_INPUT_HASH),('worker_sha256',WORKER_HASH)]:
        check(result.get(k) == expected, f'Endpoint metadata differs: {k}')
    check(sha(folder/'input.json') == RECOVERY_INPUT_HASH, 'Actual worker input differs')
    check(result.get('scf',{}).get('converged') is True, 'Final SCF not converged')
    check(result['scf'].get('conv_tol_hartree') == req['scf_conv_tol'] and result['scf'].get('conv_tol_grad') == req['scf_conv_tol_grad'],
          'Final SCF thresholds differ')
    opt = result.get('optimization',{})
    check(opt.get('requested') is True and opt.get('converged') is True, 'Geometry optimization not converged')
    expected_opt = load(ctx['base']/'validated-input.json')['optimization']
    check(opt.get('settings') == expected_opt, 'Actual optimizer settings differ')
    check(opt.get('constraint_atom_indices_zero_based') == [ctx['index'][a] for a in expected_opt['freeze_atom_ids']],
          'Actual Cartesian cap constraint indices differ')
    constraint_bytes=('$freeze\nxyz '+','.join(str(ctx['index'][a]+1) for a in expected_opt['freeze_atom_ids'])+'\n').encode()
    check(raw(folder/'constraints.txt') == constraint_bytes and opt.get('constraint_file_sha256') == hashlib.sha256(constraint_bytes).hexdigest(),
          'Actual constraint file differs from the three-carbon-only request')
    check(opt.get('dihedral_constraints_satisfied') is True and opt.get('final_dihedral_constraints') == [],
          'Unexpected or unverified dihedral constraints')
    check(opt.get('steps') and all(s.get('scf_converged') is True and math.isfinite(s.get('energy_hartree',float('nan'))) for s in opt['steps']),
          'Unconverged/nonfinite optimization step evidence')
    check(math.isfinite(result.get('energy_hartree',float('nan'))), 'Nonfinite endpoint energy')
    density_count=result.get('density_electron_count',float('nan'))
    check(math.isfinite(density_count) and abs(density_count-372)<1e-6, 'Endpoint density electron count differs')
    check(result.get('units',{}).get('bohr_to_angstrom') == BOHR and
          result['units'].get('gradient_hartree_per_bohr') == 'hartree/bohr', 'Endpoint derivative/coordinate units differ')
    for name,key in [('arrays.npz','arrays_sha256'),('scf.chk','checkpoint_sha256'),('resolved-method.json','resolved_method_file_sha256')]:
        check(result.get({'arrays.npz':'arrays_file','scf.chk':'checkpoint_file','resolved-method.json':'resolved_method_file'}[name]) == name,
              'Unexpected result artifact layout')
        check(sha(folder/name) == result.get(key), f'Endpoint artifact checksum failed: {name}')
    reference_path = Path(req['initial_checkpoint']['result_path']).parent/'resolved-method.json'
    reference_result = load(Path(req['initial_checkpoint']['result_path']),req['initial_checkpoint']['result_sha256'])
    reference = load(reference_path,reference_result['resolved_method_file_sha256'])
    method = load(folder/'resolved-method.json')
    for k in ['method','basis_sha256','ecp_sha256','spherical_basis','symmetry','density_fit',
              'functional_resolved_sha256','grid_level','auxbasis_sha256']:
        check(method.get(k) == reference.get(k), f'Resolved endpoint method changed: {k}')
    runtime = RECOVERY_RUNTIME_HASH
    check(result.get('runtime_manifest',{}).get('sha256') == runtime, 'Declared runtime manifest differs')
    check(sha(folder/'runtime-manifest.json') == runtime, 'Copied runtime manifest differs')
    check(result.get('actual_pyscf_threads') == 2 and result.get('resources',{}).get('max_wall_seconds') == 21600,
          'Actual resource settings differ')
    check(result['resources'].get('max_memory_mb') == 2000 and result.get('versions',{}).get('pyscf') == '2.14.0',
          'Requested memory/provider version differs')
    arrays = {}
    with np.load(folder/'arrays.npz', allow_pickle=False) as payload:
        expected_arrays=['coords_bohr','coords_angstrom','gradient_hartree_per_bohr','atomic_numbers','effective_nuclear_charges']
        check(set(payload.files) == set(expected_arrays), 'Missing or unrequested endpoint arrays')
        for key in expected_arrays:
            check(key in payload, f'Missing endpoint array: {key}')
            value = np.asarray(payload[key])
            shape = (97,3) if key in ['coords_bohr','coords_angstrom','gradient_hartree_per_bohr'] else (97,)
            check(value.shape == shape and np.issubdtype(value.dtype,np.number) and np.isfinite(value).all(), f'Malformed/nonfinite endpoint array: {key}')
            arrays[key] = value.copy()
    numbers = np.asarray([Z[e] for e in ctx['elements']])
    check(np.array_equal(arrays['atomic_numbers'],numbers) and np.array_equal(arrays['effective_nuclear_charges'],numbers),
          'Atomic number/effective nuclear charge map differs')
    check(np.max(np.abs(arrays['coords_bohr']*BOHR-arrays['coords_angstrom'])) <= 1e-12, 'Bohr/Angstrom coordinate arrays disagree')
    return {'folder':folder,'result':result,'arrays':arrays,'hashes':{
        name:sha(folder/name) for name in ['result.json','input.json','arrays.npz','scf.chk','resolved-method.json','runtime-manifest.json']}}


def read_recovery_endpoint(ctx, spec):
    run = Path(spec['output'])
    if not (run/'controller-result.json').is_file():
        return None  # Never inspect live result, gradient arrays or checkpoint.
    controller = load(run/'controller-result.json')
    check(controller.get('status') == 'numerically_complete_pending_geometry_review',
          'Recovery terminal controller did not complete numerically')
    folder = run/'optimization'
    validate_native_inventory(folder)
    before = {name:sha(folder/name) for name in NATIVE_FILES}
    validate_controller(controller,before['result.json'])
    check(load(run/'spec.json') == spec, 'Captured controller specification differs')
    for path,digest in [(run/'controller-source.py',CONTROLLER_HASH),(run/'qm_worker.py',WORKER_HASH),
                        (run/'input.json',RECOVERY_INPUT_HASH)]:
        check(sha(path) == digest, 'Captured execution source changed: '+str(path))
    endpoint=read_native_endpoint(ctx,folder)
    check(endpoint is not None, 'Terminal native result missing')
    validate_real_arrays(endpoint['arrays'])
    check(endpoint['result'].get('initial_checkpoint') is None, 'Unexpected checkpoint reuse')
    expected_resolved = load(ctx['base']/'validated-input.json')
    expected_resolved.pop('initial_checkpoint')
    expected_resolved.update(coords_angstrom=load(run/'input.json')['coords_angstrom'],max_memory_mb=2000)
    check(load(folder/'resolved-input.json') == expected_resolved, 'Resolved recovery input changed')
    check(load(folder/'optimization-report.json') == endpoint['result']['optimization'], 'Optimizer report differs')
    check(load(folder/'progress.json') == {'phase':'completed','status':'completed','accepted':True}, 'Worker progress not terminal')
    check({name:sha(folder/name) for name in NATIVE_FILES} == before, 'Native evidence changed during read')
    check(load(run/'controller-result.json') == controller, 'Controller evidence changed during read')
    endpoint['hashes'].update({'controller-result.json':sha(run/'controller-result.json')})
    return endpoint


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args()
    args.output.mkdir(parents=True,exist_ok=False)
    report={'checked_unix':time.time(),'QM_launched':False,'model_accuracy_validated':False,
            'full_preparation_ready':False,'analyzer_sha256':sha(Path(__file__))}
    try:
        ctx=core.load_context()  # Original request must remain original crystal coordinates.
        spec,lineage=load_lineage(ctx)
        report['lineage']=lineage
        endpoint=read_recovery_endpoint(ctx,spec)
        if endpoint is None:
            report['status']='lineage_ready_terminal_endpoint_not_present_no_native_result_read'
        else:
            analysis=core.chemistry_report(ctx,endpoint['arrays']['coords_angstrom'],endpoint['arrays']['gradient_hartree_per_bohr'])
            report.update(status=analysis['status'],analysis=analysis,endpoint_hashes=endpoint['hashes'],
                          endpoint_energy_hartree=endpoint['result']['energy_hartree'])
        report['downstream_request_generated']=False
        report['next_gate']='Only after all original endpoint checks pass: separately reviewed exact-endpoint conventional gradient request; no automatic launch or Hessian.'
    except Exception as error:
        report.update(status='evidence_or_analysis_failed',error_type=type(error).__name__,error=str(error))
    core.write_new(args.output/'result.json',report)
    print(json.dumps({k:v for k,v in report.items() if k not in ['analysis','lineage']},indent=2))
    if report['status']=='evidence_or_analysis_failed':
        raise SystemExit(1)


if __name__ == '__main__':
    main()
