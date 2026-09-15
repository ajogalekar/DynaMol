"""Standalone research candidate for real float64 C[p,q] = sum_ij A[p,i,j] B[q,j,i].

No symmetry, conjugation, approximation, scientific unit conversion, or global
transpose/reshape materialization is performed. Only positive, compatible 3-D
plain NumPy arrays are admitted; this matches the nonempty real rhok_oo tensors
at the reviewed PySCF DF-gradient call. This is NOT integrated into PySCF.

The plan bounds explicitly allocated NumPy arrays, excluding pre-existing inputs,
Python/NumPy objects, BLAS/native-library workspaces, and process RSS. Input/output
finite checks reuse one bounded Boolean tile, never a full-operand Boolean copy.
Concurrent mutation of input arrays is unsupported. Inputs are not modified.
"""
from dataclasses import dataclass
import numpy as np


@dataclass(frozen=True)
class WorkspacePlan:
    a_shape: tuple
    b_shape: tuple
    block_p: int
    block_q: int
    inner_size: int
    left_bytes: int
    right_bytes: int
    product_bytes: int
    finite_mask_bytes: int
    workspace_bytes: int
    output_bytes: int
    managed_array_bytes: int


def _positive_integer(value, label):
    if type(value) is not int or value < 1:
        raise ValueError(f'{label} must be a positive integer')
    return value


def plan_workspace(a_shape, b_shape, *, aux_block_size=64,
                   max_workspace_bytes=64 * 1024**2, max_output_bytes=256 * 1024**2):
    """Validate dimensions/budgets and calculate bytes without allocating tensors."""
    for shape in (a_shape, b_shape):
        if not isinstance(shape, tuple) or len(shape) != 3:
            raise ValueError('Each input shape must be a three-dimensional tuple')
        for dimension in shape:
            _positive_integer(dimension, 'Input dimension')
    if a_shape[1:] != b_shape[1:][::-1]:
        raise ValueError('Expected A[P,I,J] and B[Q,J,I] with matching inner dimensions')
    block = _positive_integer(aux_block_size, 'aux_block_size')
    workspace_limit = _positive_integer(max_workspace_bytes, 'max_workspace_bytes')
    output_limit = _positive_integer(max_output_bytes, 'max_output_bytes')
    p, i, j = a_shape
    q = b_shape[0]
    bp, bq, k = min(block, p), min(block, q), i*j
    left, right, product = 8*bp*k, 8*bq*k, 8*bp*bq
    mask = max(bp*k, bq*k, bp*bq) * np.dtype(np.bool_).itemsize
    workspace, output = left+right+product+mask, 8*p*q
    if workspace > workspace_limit:
        raise ValueError(f'Explicit tile workspace requires {workspace} bytes, over limit {workspace_limit}')
    if output > output_limit:
        raise ValueError(f'Contraction output requires {output} bytes, over limit {output_limit}')
    return WorkspacePlan(a_shape, b_shape, bp, bq, k, left, right, product,
                         mask, workspace, output, workspace+output)


def _array(value, name):
    if type(value) is not np.ndarray:
        raise ValueError(f'{name} must be a plain NumPy ndarray; implicit conversion is forbidden')
    if value.dtype != np.dtype(np.float64) or not value.dtype.isnative:
        raise ValueError(f'{name} must have native real float64 dtype')
    if value.ndim != 3:
        raise ValueError(f'{name} must have rank three')
    return value


def _finite(array, mask, label):
    scratch = mask[:array.size].reshape(array.shape)
    np.isfinite(array, out=scratch)
    if not bool(scratch.all()):
        raise ValueError(f'Nonfinite {label}; no result accepted')


def auxiliary_exchange_contract(a, b, *, aux_block_size=64,
                                max_workspace_bytes=64 * 1024**2,
                                max_output_bytes=256 * 1024**2):
    """Return the exact contraction algebra, with ordinary float64 rounding.

    All P,Q entries are calculated independently. A and B may be C-/F-contiguous,
    read-only, sliced, reversed or aliased. There is no symmetry assumption.
    The I,J reduction is not truncated or split, and no physical term is omitted.
    A failed finite/overflow check raises; a partial result is never returned.
    """
    a, b = _array(a, 'A'), _array(b, 'B')
    plan = plan_workspace(a.shape, b.shape, aux_block_size=aux_block_size,
                          max_workspace_bytes=max_workspace_bytes,
                          max_output_bytes=max_output_bytes)
    p, i, j = a.shape
    q, bp, bq, k = b.shape[0], plan.block_p, plan.block_q, plan.inner_size
    result = np.empty((p, q), dtype=np.float64)
    left = np.empty((bp, k), dtype=np.float64)
    right = np.empty((bq, k), dtype=np.float64)
    product = np.empty(bp*bq, dtype=np.float64)
    finite_mask = np.empty(plan.finite_mask_bytes, dtype=np.bool_)
    for p0 in range(0, p, bp):
        p1 = min(p0+bp, p)
        pw = p1-p0
        left_view = left[:pw]
        np.copyto(left_view.reshape(pw, i, j), a[p0:p1], casting='no')
        _finite(left_view, finite_mask, 'left operand')
        for q0 in range(0, q, bq):
            q1 = min(q0+bq, q)
            qw = q1-q0
            right_view = right[:qw]
            # transpose() is a view; copy only this auxiliary-axis tile directly
            # into its final packed layout. Never reshape the whole operand.
            np.copyto(right_view.reshape(qw, i, j),
                      b[q0:q1].transpose(0, 2, 1), casting='no')
            _finite(right_view, finite_mask, 'right operand')
            # The shortened edge tile remains C-contiguous; do not use a strided
            # rectangular view of a larger output matrix as the BLAS destination.
            tile = product[:pw*qw].reshape(pw, qw)
            try:
                with np.errstate(over='raise', invalid='raise'):
                    np.matmul(left_view, right_view.T, out=tile)
            except FloatingPointError as error:
                raise ValueError('Nonfinite contraction result; no result accepted') from error
            _finite(tile, finite_mask, 'contraction result')
            np.copyto(result[p0:p1, q0:q1], tile, casting='no')
    return result
