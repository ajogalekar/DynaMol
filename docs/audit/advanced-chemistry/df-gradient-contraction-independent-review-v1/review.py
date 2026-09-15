"""Independent tiny-array review. NumPy only; no PySCF or molecular evaluation."""
import dataclasses
import hashlib
import importlib.util
import json
import math
from pathlib import Path
import sys
from unittest.mock import patch
import numpy as np

HERE=Path(__file__).resolve().parent
SOURCE=HERE.parent/'df-gradient-contraction-candidate-v1/tiled_contract.py'
EXPECTED='8d9e6bbb64e9e525ac2d41cb2783f699050e32c0516e89cb821a1229b1012060'
raw=SOURCE.read_bytes();assert hashlib.sha256(raw).hexdigest()==EXPECTED
spec=importlib.util.spec_from_file_location('reviewed_tiled_contract',SOURCE)
candidate=importlib.util.module_from_spec(spec);sys.modules[spec.name]=candidate;spec.loader.exec_module(candidate)
results=[]


def reference(a,b):
    values=np.zeros((a.shape[0],b.shape[0]));bounds=values.copy()
    for p in range(a.shape[0]):
        for q in range(b.shape[0]):
            terms=[float(a[p,i,j])*float(b[q,j,i]) for i in range(a.shape[1]) for j in range(a.shape[2])]
            values[p,q]=math.fsum(terms)
            bounds[p,q]=8*np.finfo(float).eps*(len(terms)+1)*math.fsum(abs(t) for t in terms)+8*np.finfo(float).eps
    return values,bounds


def check(label,a,b,block):
    expected,bound=reference(a,b)
    left_before,right_before=a.copy(),b.copy()
    got=candidate.auxiliary_exchange_contract(a,b,aux_block_size=block)
    assert np.all(np.abs(got-expected)<=bound)
    assert np.array_equal(a,left_before) and np.array_equal(b,right_before)
    assert got.flags.c_contiguous and got.dtype==np.dtype('float64')
    results.append({'label':label,'a_shape':list(a.shape),'b_shape':list(b.shape),'block':block,
        'max_absolute_error':float(np.max(np.abs(got-expected))),
        'max_roundoff_envelope_ratio':float(np.max(np.abs(got-expected)/bound))})
    return got


def main():
    rng=np.random.default_rng(195501)
    for p,q,i,j in [(1,1,1,1),(1,5,3,1),(6,1,1,4),(5,7,2,5),
                    (7,5,5,2),(3,4,3,3),(9,10,1,7),(10,9,7,1)]:
        a=rng.normal(size=(p,i,j));b=rng.normal(size=(q,j,i))
        for block in (1,2,3,64):check('shape-singleton-tail',a,b,block)
    base_a=rng.normal(size=(10,4,10));base_b=rng.normal(size=(14,10,4))
    layouts=[('C/F',base_a[::2,::2,::2].copy(),np.asfortranarray(base_b[::2,::2,::2])),
        ('negative-strides',base_a[::-2,::-2,::-2],base_b[::-2,::-2,::-2]),
        ('broadcast-zero-strides',np.broadcast_to(rng.normal(size=(1,2,5)),(5,2,5)),
         np.broadcast_to(rng.normal(size=(1,5,2)),(7,5,2)))]
    for label,a,b in layouts:
        for block in (1,2,3,64):check(label,a,b,block)
    a=np.arange(50.,dtype=np.float64).reshape(5,2,5)+1
    b=(np.arange(70.,dtype=np.float64).reshape(7,5,2)+3)*2
    got=check('rectangular-nonsymmetric-index-discriminator',a,b,3)
    wrong=a.reshape(5,10)@b.reshape(7,10).T
    wrong_gap=float(np.max(np.abs(got-wrong)))
    assert wrong_gap>1 and np.array_equal(got,reference(a,b)[0])

    plan=candidate.plan_workspace(a.shape,b.shape,aux_block_size=3)
    allocations=[];finite_tiles=[];matmul_tiles=[];copy_shapes=[]
    original_empty,original_finite,original_matmul,original_copy=np.empty,np.isfinite,np.matmul,np.copyto
    def empty(*args,**kwargs):
        value=original_empty(*args,**kwargs);allocations.append(value);return value
    def finite(array,*,out):
        assert np.shares_memory(out,allocations[4]);assert out.nbytes<=plan.finite_mask_bytes
        finite_tiles.append([list(array.shape),out.nbytes]);return original_finite(array,out=out)
    def matmul(left,right,*,out):
        assert np.shares_memory(left,allocations[1]) and np.shares_memory(right,allocations[2])
        assert np.shares_memory(out,allocations[3])
        assert left.flags.c_contiguous and right.flags.f_contiguous and out.flags.c_contiguous
        matmul_tiles.append(list(out.shape));return original_matmul(left,right,out=out)
    def copyto(dest,src,**kwargs):
        assert any(np.shares_memory(dest,v) for v in allocations)
        copy_shapes.append(list(dest.shape));return original_copy(dest,src,**kwargs)
    with patch.object(candidate.np,'empty',empty),patch.object(candidate.np,'isfinite',finite),\
         patch.object(candidate.np,'matmul',matmul),patch.object(candidate.np,'copyto',copyto):
        candidate.auxiliary_exchange_contract(a,b,aux_block_size=3,
            max_workspace_bytes=plan.workspace_bytes,max_output_bytes=plan.output_bytes)
    assert len(allocations)==5 and sum(x.nbytes for x in allocations)==plan.managed_array_bytes
    assert [2,1] in matmul_tiles
    allocation_report={'plan':dataclasses.asdict(plan),'allocated_bytes':[x.nbytes for x in allocations],
        'matmul_output_shapes':matmul_tiles,'finite_tiles':finite_tiles,'copy_shapes':copy_shapes,
        'tracked_explicit_arrays_only':True,'blas_internal_memory_or_rss_measured':False}

    rejection_cases=[]
    def rejected(label,left,right,**kw):
        try:candidate.auxiliary_exchange_contract(left,right,aux_block_size=3,**kw)
        except ValueError:rejection_cases.append(label)
        else:raise AssertionError(label+' unexpectedly admitted')
    for key,value in [('max_workspace_bytes',plan.workspace_bytes-1),('max_output_bytes',plan.output_bytes-1)]:
        with patch.object(candidate.np,'empty',side_effect=AssertionError('budget failure allocated')):
            rejected(key,a,b,**{key:value})
    for label,bad in [('float32',a.astype('float32')),('complex',a.astype('complex128')),
                     ('non-native-endian',a.astype('>f8' if sys.byteorder=='little' else '<f8')),
                     ('subclass',a.view(type('ArraySubclass',(np.ndarray,),{}))),
                     ('empty',a[:0])]:rejected(label,bad,b)
    for operand in ('left','right'):
        for value in (np.nan,np.inf,-np.inf):
            aa,bb=a.copy(),b.copy();(aa if operand=='left' else bb)[-1,-1,-1]=value
            rejected(operand+' late '+str(value),aa,bb)
    rejected('finite operands overflow',np.full(a.shape,1e308),np.full(b.shape,1e308))
    # Dimension planning must not allocate any real-size tensor.
    with patch.object(candidate.np,'empty',side_effect=AssertionError('large tensor forbidden')):
        real_plan=candidate.plan_workspace((4739,194,194),(4739,194,194))
    assert real_plan.workspace_bytes==40980736 and real_plan.output_bytes==179664968
    report={'schema_version':1,'verdict':'no_blocker_for_standalone_small_array_candidate',
        'candidate_sha256':EXPECTED,'numpy_version':np.__version__,'numeric_cases':results,
        'numeric_case_count':len(results),'max_absolute_error':max(x['max_absolute_error'] for x in results),
        'rectangular_wrong_index_formula_max_difference':wrong_gap,
        'rejection_cases':rejection_cases,'allocation_review':allocation_report,
        'real_shape_plan_only':dataclasses.asdict(real_plan),'full_size_arrays_allocated':False,
        'pyscf_imported':any(k=='pyscf' or k.startswith('pyscf.') for k in sys.modules),
        'native_qm_evaluations':0,'accepted':False,'physical_acceptance':False,
        'limits':['Explicit-array budget excludes input tensors, Python/allocator overhead and BLAS/native workspaces.',
          'Float64 reduction ordering need not be bitwise identical to the original einsum.',
          'Full DF-gradient dispatch, finite-difference parity, actual peak RSS and full-scale speed remain unqualified.']}
    assert report['pyscf_imported'] is False
    (HERE/'result.json').write_text(json.dumps(report,indent=2)+'\n')
    (HERE/'reviewed-source.py').write_bytes(raw)
    print(json.dumps({k:report[k] for k in ['verdict','candidate_sha256','numpy_version','numeric_case_count',
        'max_absolute_error','rectangular_wrong_index_formula_max_difference','native_qm_evaluations']},indent=2))


if __name__=='__main__':main()
