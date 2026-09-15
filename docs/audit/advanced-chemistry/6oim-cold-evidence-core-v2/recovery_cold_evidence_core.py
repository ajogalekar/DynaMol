"""Read-only recovery-aware conventional-HF ESP and native RESP evidence branch.

No process is spawned. No checkpoint is reused. Returning numerically bound
charges does not validate a force field or authorize a simulation.
"""
import ast
import copy
import hashlib
import importlib.util
import io
import json
import math
from pathlib import Path
import zipfile

MATERIALIZER_SHA='698c6d6727e365860e3d00debbf060c216425e31f6f7e1a6ac7d55eb0237fab2'
MATERIALIZER_PLAN_SHA='8753021db5fe96d71fffcb29e7198b4ecf7b8668343389c9f110b81d76e43a20'
HANDOFF_SHA='b9548748219357d76478158a91af955fe9794885cb9f12d84a96bb6182a3d1af'
MAX_ARTIFACT_BYTES=64*1024*1024
MATERIALIZED_FILES={'result.json','exact-esp-input.json','parent-binding.json','geometry.json','downstream-contract.json'}
ESP_NATIVE_FILES={'input.json','resolved-input.json','native.log','runtime-manifest.json','progress.json',
                  'resolved-method.json','arrays.npz','scf.chk','result.json'}


def require(condition,message):
    if not condition:raise ValueError(message)


def data(path,expected=None):
    path=Path(path)
    require(path.is_file() and not path.is_symlink() and path.stat().st_size<=MAX_ARTIFACT_BYTES,
            'Missing, nonregular or oversized evidence: '+str(path))
    value=path.read_bytes();require(len(value)==path.stat().st_size,'Incomplete evidence read')
    if expected is not None:require(hashlib.sha256(value).hexdigest()==expected,'Evidence digest mismatch: '+str(path))
    return value


def sha(path):return hashlib.sha256(data(path)).hexdigest()
def load(path,expected=None):return json.loads(data(path,expected))
def digest(value):return hashlib.sha256(json.dumps(value,sort_keys=True,separators=(',',':'),allow_nan=False).encode()).hexdigest()


def materializer_module():
    path=Path(__file__).with_name('materialize_6oim_recovery_esp.py')
    source=data(path,MATERIALIZER_SHA)
    import types
    module=types.ModuleType('pinned_no_launch_materializer');module.__file__=str(path)
    exec(compile(source,str(path),'exec'),module.__dict__)
    return module


def readonly_parent_function(module):
    """Reuse the exact pinned admission statements before publication, unmodified.

    The only new AST nodes are the function signature and a return of existing
    validated locals. This function neither copies nor weakens scientific gates.
    """
    source=data(module.__file__,MATERIALIZER_SHA)
    function=next(x for x in ast.parse(source).body if isinstance(x,ast.FunctionDef) and x.name=='materialize')
    block=next(x for x in function.body if isinstance(x,ast.Try)).body
    require(isinstance(block[0],ast.Assign) and isinstance(block[1],ast.Assign),'Materializer preamble changed')
    require(ast.unparse(block[0])=='plan = load(plan_path)' and
            ast.unparse(block[1])=="report['plan_sha256'] = sha(plan_path)",'Materializer read-only extraction preamble changed')
    end=next(i for i,x in enumerate(block) if isinstance(x,ast.Expr) and isinstance(x.value,ast.Call)
             and isinstance(x.value.func,ast.Name) and x.value.func.id=='publish')
    selected=copy.deepcopy(block[2:end])
    require(not any(isinstance(x,ast.Name) and x.id in ('report','output','publish') for s in selected for x in ast.walk(s)),
            'Read-only admission unexpectedly references publication state')
    names=['plan','source','request','result','arrays','method','geometry','optimized','esp_request','before','controller',
           'qualification','resp_provider']
    returned=ast.Return(value=ast.Dict(keys=[ast.Constant(x) for x in names],values=[ast.Name(id=x,ctx=ast.Load()) for x in names]))
    template=ast.parse('def read_parent(plan):\n    pass\n').body[0]
    template.body=selected+[returned]
    namespace=dict(module.__dict__)
    exec(compile(ast.fix_missing_locations(ast.Module(body=[template],type_ignores=[])),module.__file__,'exec'),namespace)
    return namespace['read_parent']


def revalidate_actual_parent(plan_path):
    """Recheck actual terminal parent and runtime/source/method/geometry evidence.

    Requires the admitted preparation interpreter. The returned dictionary
    contains data for admission, never physical acceptance or fitted charges.
    An unfinished parent raises before native arrays/checkpoints or imports.
    """
    module=materializer_module()
    plan=load(plan_path,MATERIALIZER_PLAN_SHA)
    context=readonly_parent_function(module)(plan)
    context['materializer']=module
    context['plan_path']=Path(plan_path)
    return context


def snapshot(folder,names):
    folder=Path(folder)
    require(folder.is_dir() and not folder.is_symlink() and {p.name for p in folder.iterdir()}==set(names),
            'Missing or extra evidence inventory: '+str(folder))
    return {name:data(folder/name) for name in names}


def hashes(snapshot_):return {k:hashlib.sha256(v).hexdigest() for k,v in snapshot_.items()}


def expected_materialized_payloads(context):
    """Evaluate the exact pinned publication expressions against rechecked locals."""
    module=context['materializer'];source=data(module.__file__,MATERIALIZER_SHA)
    function=next(x for x in ast.parse(source).body if isinstance(x,ast.FunctionDef) and x.name=='materialize')
    block=next(x for x in function.body if isinstance(x,ast.Try)).body
    namespace=dict(module.__dict__);namespace.update(context)
    namespace.update(coords=context['arrays']['coords_bohr'].tolist(),
                     admission=context['plan']['preparation_runtime_admission'],expected_pins=context['plan']['source_pins'])
    outputs={}
    for statement in block:
        if not (isinstance(statement,ast.Expr) and isinstance(statement.value,ast.Call) and
                isinstance(statement.value.func,ast.Name) and statement.value.func.id=='publish'):continue
        call=statement.value
        require(len(call.args)==2 and isinstance(call.args[0],ast.BinOp) and
                isinstance(call.args[0].right,ast.Constant),'Pinned materializer publication changed')
        name=call.args[0].right.value
        outputs[name]=eval(compile(ast.Expression(body=copy.deepcopy(call.args[1])),module.__file__,'eval'),namespace)
    require(set(outputs)==MATERIALIZED_FILES-{'result.json'},'Pinned materialized payload set changed')
    return outputs


def expected_materialized_report(folder,context):
    module=context['materializer'];source=data(module.__file__,MATERIALIZER_SHA)
    function=next(x for x in ast.parse(source).body if isinstance(x,ast.FunctionDef) and x.name=='materialize')
    initial=next(x.value for x in function.body if isinstance(x,ast.Assign) and
                 any(isinstance(t,ast.Name) and t.id=='report' for t in x.targets))
    block=next(x for x in function.body if isinstance(x,ast.Try)).body
    final=next(x.value for x in block if isinstance(x,ast.Expr) and isinstance(x.value,ast.Call) and
               isinstance(x.value.func,ast.Attribute) and isinstance(x.value.func.value,ast.Name) and
               x.value.func.value.id=='report' and x.value.func.attr=='update')
    namespace=dict(module.__dict__);namespace.update(context);namespace['output']=Path(folder)
    report=eval(compile(ast.Expression(body=copy.deepcopy(initial)),module.__file__,'eval'),namespace)
    report['plan_sha256']=MATERIALIZER_PLAN_SHA
    require(not final.args and all(k.arg is not None for k in final.keywords),'Pinned successful report changed')
    report.update({k.arg:eval(compile(ast.Expression(body=copy.deepcopy(k.value)),module.__file__,'eval'),namespace)
                   for k in final.keywords})
    return report


def verify_materialized_request(folder,context):
    """Bind the complete unlaunched request set back to independently rechecked parent."""
    folder=Path(folder);module=context['materializer']
    require(folder.resolve().is_relative_to(module.CACHE),'Materialized request is outside stable local cache')
    files=snapshot(folder,MATERIALIZED_FILES);pins=hashes(files)
    parsed={k:json.loads(v) for k,v in files.items()}
    result=parsed['result.json'];request=parsed['exact-esp-input.json'];binding=parsed['parent-binding.json']
    require(result.get('status')=='cold_esp_request_materialized_not_launched' and result.get('request_emitted') is True,
            'Materializer did not finish creating an ESP request')
    for key in ('accepted','physical_acceptance','QM_launched','RESP_launched','checkpoint_reused'):
        require(result.get(key) is False,'Unexpected materialization acceptance/execution flag: '+key)
    require(result.get('script_sha256')==MATERIALIZER_SHA and result.get('plan_sha256')==MATERIALIZER_PLAN_SHA and
            result.get('input_sha256')==pins['exact-esp-input.json'] and
            result.get('parent_result_sha256')==context['before']['result.json'],'Materializer source/input/parent binding changed')
    expected_payloads=expected_materialized_payloads(context)
    for name,expected in expected_payloads.items():
        require(digest(parsed[name])==digest(expected),'Complete materialized payload changed: '+name)
    require(digest(result)==digest(expected_materialized_report(folder,context)),
            'Complete materializer terminal report changed')
    require(binding.get('native_hashes')==context['before'] and binding.get('controller')==context['controller'] and
            binding.get('geometry_sha256')==digest(context['arrays']['coords_bohr'].tolist()) and
            binding.get('method_sha256')==digest(context['method']),'Materializer parent evidence changed')
    require(binding.get('checkpoint_policy')==module.COLD and binding.get('checkpoint_geometry_validated') is False and
            binding.get('accepted') is False,'Materialized checkpoint branch changed')
    require(parsed['geometry.json']==context['geometry'],'Materialized geometry check changed')
    require(request==context['esp_request'],'Materialized request changed after parent admission')
    module.validate_cold_request(request,context['request'],context['arrays']['coords_bohr'].tolist(),
                                 context['esp_request']['esp_points_bohr'],module.COLD)
    require(hashes(snapshot(folder,MATERIALIZED_FILES))==pins,'Materialized evidence changed during read')
    return {'folder':folder,'files':files,'hashes':pins,'request':request,'result':result,'parent_binding':binding}


def real_arrays(blob,expected_shapes):
    """Read only finite real native arrays of the exact requested small inventory."""
    import numpy as np
    with zipfile.ZipFile(io.BytesIO(blob)) as archive:
        require(sum(x.file_size for x in archive.infolist())<=MAX_ARTIFACT_BYTES,'Oversized uncompressed native arrays')
    with np.load(io.BytesIO(blob),allow_pickle=False) as archive:
        require(set(archive.files)==set(expected_shapes),'Native array inventory changed')
        result={}
        for name,shape in expected_shapes.items():
            array=archive[name]
            require(array.shape==shape and array.dtype.kind in 'if' and np.isfinite(array).all(),
                    'Malformed/nonfinite real native array: '+name)
            result[name]=array.copy()
    return result
