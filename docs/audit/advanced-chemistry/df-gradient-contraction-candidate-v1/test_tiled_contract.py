"""Small deterministic array tests only; no quantum evaluations or full tensors."""
import dataclasses
import hashlib
import json
import math
import os
from pathlib import Path
import time
import unittest
from unittest.mock import patch
import numpy as np
from pyscf import lib
import tiled_contract as candidate

lib.num_threads(1)
CASE_RESULTS = []


def scalar_reference(a, b):
    result = np.empty((a.shape[0], b.shape[0]))
    magnitude = np.empty_like(result)
    for p in range(a.shape[0]):
        for q in range(b.shape[0]):
            terms = [float(a[p,i,j])*float(b[q,j,i])
                     for i in range(a.shape[1]) for j in range(a.shape[2])]
            result[p,q] = math.fsum(terms)
            magnitude[p,q] = math.fsum(abs(x) for x in terms)
    return result, magnitude


class ContractionTests(unittest.TestCase):
    def test_layouts_blocks_native_and_scalar_reference(self):
        rng = np.random.default_rng(94173)
        a, b = rng.normal(size=(11,9,8)), rng.normal(size=(13,8,9))
        strided_a = rng.normal(size=(22,18,16))[::2,::2,::2]
        strided_b = rng.normal(size=(26,16,18))[::2,::2,::2]
        readonly_a, readonly_b = a.copy(), b.copy()
        readonly_a.setflags(write=False); readonly_b.setflags(write=False)
        shared = rng.normal(size=(9,8,8))
        cases = [('CC',a,b), ('CF',a,np.asfortranarray(b)),
                 ('FC',np.asfortranarray(a),b), ('FF',np.asfortranarray(a),np.asfortranarray(b)),
                 ('strided',strided_a,strided_b), ('reversed',a[::-1,::-1,::-1],b[::-1,::-1,::-1]),
                 ('transposed_inner',a.swapaxes(1,2),b.swapaxes(1,2)),
                 ('readonly',readonly_a,readonly_b), ('aliased_nonsymmetric',shared,shared),
                 ('scaled_signs',a*1e5,b*1e-5)]
        for name, left, right in cases:
            self.assertGreater(left.size, lib.numpy_helper.EINSUM_MAX_SIZE)
            self.assertGreater(right.size, lib.numpy_helper.EINSUM_MAX_SIZE)
            native = lib.einsum('pij,qji->pq', left, right)
            scalar, magnitude = scalar_reference(left,right)
            # Full reduction, no scientific tolerance: conservative forward-error
            # envelope for K float64 products/additions, compared to math.fsum.
            k = left.shape[1]*left.shape[2]
            bound = 8*np.finfo(float).eps*(k+1)*magnitude + 8*np.finfo(float).eps
            self.assertTrue(np.all(np.abs(native-scalar) <= bound))
            before_a, before_b = left.copy(), right.copy()
            for block in (1,2,4,64):
                with self.subTest(case=name,block=block):
                    got = candidate.auxiliary_exchange_contract(left,right,aux_block_size=block)
                    self.assertEqual(got.dtype,np.dtype(np.float64))
                    self.assertTrue(got.flags.c_contiguous)
                    self.assertTrue(np.all(np.abs(got-scalar) <= bound))
                    self.assertTrue(np.all(np.abs(got-native) <= 2*bound))
                    np.testing.assert_array_equal(left,before_a)
                    np.testing.assert_array_equal(right,before_b)
                    CASE_RESULTS.append({'case':name,'block':block,'a_shape':left.shape,'b_shape':right.shape,
                        'max_abs_vs_native':float(np.max(np.abs(got-native))),
                        'max_abs_vs_scalar':float(np.max(np.abs(got-scalar))),
                        'max_scalar_error_bound_ratio':float(np.max(np.abs(got-scalar)/bound))})

    def test_index_order_no_symmetry_or_conjugation_assumption(self):
        a=np.array([[[1.,2.],[3.,7.]],[[5.,11.],[13.,17.]],[[19.,23.],[29.,31.]]])
        b=np.array([[[2.,5.],[7.,11.]],[[13.,17.],[19.,23.]],[[29.,31.],[37.,41.]],[[43.,47.],[53.,59.]]])
        expected,_=scalar_reference(a,b)
        got=candidate.auxiliary_exchange_contract(a,b,aux_block_size=2)
        np.testing.assert_array_equal(got,expected)
        wrong=np.einsum('pij,qij->pq',a,b)
        self.assertGreater(float(np.max(np.abs(got-wrong))),1.)

    def test_managed_allocations_and_contiguous_blas_edges(self):
        a=np.arange(5*3*2,dtype=float).reshape(5,3,2)
        b=np.arange(7*2*3,dtype=float).reshape(7,2,3)
        plan=candidate.plan_workspace(a.shape,b.shape,aux_block_size=3)
        allocations=[];finite=[];products=[]
        original_empty,original_isfinite,original_matmul=np.empty,np.isfinite,np.matmul
        def empty(*args,**kwargs):
            value=original_empty(*args,**kwargs);allocations.append((value.shape,value.nbytes));return value
        def isfinite(array,*args,**kwargs):
            self.assertIn('out',kwargs);finite.append(kwargs['out'].nbytes)
            return original_isfinite(array,*args,**kwargs)
        def matmul(left,right,*,out):
            products.append((left.flags.c_contiguous,right.flags.f_contiguous,out.flags.c_contiguous,out.nbytes))
            return original_matmul(left,right,out=out)
        with patch.object(candidate.np,'empty',empty),patch.object(candidate.np,'isfinite',isfinite),patch.object(candidate.np,'matmul',matmul):
            got=candidate.auxiliary_exchange_contract(a,b,aux_block_size=3)
        self.assertEqual(len(allocations),5)
        self.assertEqual(sum(n for _,n in allocations),plan.managed_array_bytes)
        self.assertLessEqual(max(finite),plan.finite_mask_bytes)
        self.assertTrue(all(l and r and o and n<=plan.product_bytes for l,r,o,n in products))
        np.testing.assert_array_equal(got,scalar_reference(a,b)[0])

    def test_invalid_arrays_rank_dtype_shape_and_empty_refused(self):
        a=np.ones((2,3,4));b=np.ones((5,4,3))
        bad=[a.tolist(),np.ma.array(a),a.astype(np.float32),a.astype(np.int64),a.astype(np.complex128),
             a.astype('>f8'),a[0],a[...,None],np.empty((0,3,4)),np.empty((2,0,4))]
        for value in bad:
            with self.subTest(kind=str(type(value))),self.assertRaises(ValueError):
                candidate.auxiliary_exchange_contract(value,b)
        with self.assertRaises(ValueError):candidate.auxiliary_exchange_contract(a,np.ones((5,3,4)))
        with self.assertRaises(ValueError):candidate.auxiliary_exchange_contract(a,b.astype(np.float32))

    def test_budget_rejection_precedes_any_allocation(self):
        a=np.ones((2,3,4));b=np.ones((5,4,3));p=candidate.plan_workspace(a.shape,b.shape)
        for kwargs in ({'max_workspace_bytes':p.workspace_bytes-1},{'max_output_bytes':p.output_bytes-1},
                       {'aux_block_size':0},{'aux_block_size':True},{'aux_block_size':1.0},
                       {'max_workspace_bytes':False},{'max_output_bytes':-1}):
            with self.subTest(kwargs=kwargs),patch.object(candidate.np,'empty',side_effect=AssertionError('unexpected allocation')):
                with self.assertRaises(ValueError):candidate.auxiliary_exchange_contract(a,b,**kwargs)

    def test_nonfinite_inputs_in_late_tiles_refused(self):
        for operand in (0,1):
            for value in (float('nan'),float('inf'),-float('inf')):
                a=np.ones((5,3,2));b=np.ones((7,2,3));(a,b)[operand][-1,-1,-1]=value
                with self.subTest(operand=operand,value=value),self.assertRaisesRegex(ValueError,'Nonfinite'):
                    candidate.auxiliary_exchange_contract(a,b,aux_block_size=2)

    def test_overflow_refused_without_returning_partial_result(self):
        with self.assertRaisesRegex(ValueError,'Nonfinite'):
            candidate.auxiliary_exchange_contract(np.full((2,2,2),1e308),np.full((3,2,2),1e308))

    def test_real_case_planning_uses_no_tensor_allocation(self):
        with patch.object(candidate.np,'empty',side_effect=AssertionError('full-size allocation forbidden')):
            p=candidate.plan_workspace((4739,194,194),(4739,194,194))
        self.assertEqual(p.left_bytes,19269632)
        self.assertEqual(p.right_bytes,19269632)
        self.assertEqual(p.output_bytes,179664968)
        self.assertLess(p.workspace_bytes,64*1024**2)


if __name__=='__main__':
    start=time.monotonic();program=unittest.main(exit=False,verbosity=2)
    report={'status':'passed' if program.result.wasSuccessful() else 'failed',
        'test_methods':program.result.testsRun,'failures':len(program.result.failures),'errors':len(program.result.errors),
        'numerical_cases':CASE_RESULTS,'case_count':len(CASE_RESULTS),'elapsed_seconds':time.monotonic()-start,
        'numpy_version':np.__version__,'pyscf_threads':lib.num_threads(),'einsum_backend':lib.numpy_helper.EINSUM_BACKEND,
        'candidate_sha256':hashlib.sha256(Path(candidate.__file__).read_bytes()).hexdigest(),
        'real_case_plan':dataclasses.asdict(candidate.plan_workspace((4739,194,194),(4739,194,194))),
        'scope':'Small synthetic arrays and scalar shape arithmetic only; no full-size tensors, quantum evaluations or runtime integration.'}
    if path:=os.environ.get('DYNAMOL_CONTRACTION_TEST_REPORT'):Path(path).write_text(json.dumps(report,indent=2)+'\n')
    raise SystemExit(0 if program.result.wasSuccessful() else 1)
