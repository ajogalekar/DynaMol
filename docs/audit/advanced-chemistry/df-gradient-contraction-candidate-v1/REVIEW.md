# Standalone tiled DF-gradient contraction candidate

**Small-array qualification passed; not integrated.** The candidate replaces only the algebra of PySCF's `rhok_oo` contraction `pij,qji->pq` with auxiliary-axis tiles. It does not modify the active runtime, worker, gradient, chemistry or any simulation output.

Helper: `tiled_contract.py`, SHA-256 `8d9e6bbb64e9e525ac2d41cb2783f699050e32c0516e89cb821a1229b1012060`. Tests: `test_tiled_contract.py`, SHA-256 `798b3083997a012be4fbfe87bc6e56686957730e290769d65a5988130791a528`.

## Exact algebra and admitted inputs

For `A[P,I,J]` and `B[Q,J,I]`, every output entry is `C[p,q] = sum(i,j) A[p,i,j] B[q,j,i]`. There is **no symmetry assumption**, conjugation, averaging, coefficient adjustment, truncation, physical unit conversion or omitted term. In the original source, `rhok_oo` is explicitly nonsymmetric in its occupied indices because only one index carries the occupancy factor. Both inputs are native real float64 plain NumPy arrays with positive compatible dimensions, as produced by the reviewed call. Complex, float32, foreign-endian, object-like and empty inputs fail explicitly; there is no implicit conversion or approximate fallback.

Each left tile is copied into its final contiguous `(P_tile,I*J)` buffer. A transpose **view of one right tile** is copied directly into `(Q_tile,I*J)` in the matching index order. Their matrix product writes into a reusable contiguous product buffer, then into the corresponding output slice. Every operand copy is bounded by the configured tile, including when a small operand fits in one tile. The entire I,J reduction is retained; ordinary BLAS rounding can differ from the original order, so bitwise identity is not promised.

## Managed memory

With tile size64, input shapes `(4739,194,194)` and `(4739,194,194)`, shape-only planning gives:

| Explicit array | Bytes |
|---|---:|
| Left operand tile | 19,269,632 |
| Right operand tile | 19,269,632 |
| Product tile | 32,768 |
| Reusable finite mask | 2,408,704 |
| **Temporary workspace** | **40,980,736** |
| Required full output | 179,664,968 |
| **Total managed arrays** | **220,645,704** |

The original layout can temporarily create two additional full occupied/auxiliary copies totaling 2,853,712,064 bytes. The candidate uses no whole-input or whole-output Boolean finite mask: it checks packed tiles and result tiles into the reusable bounded mask. Budget checks happen before all array allocation; default limits are64MiB explicit workspace and256MiB output. Larger shapes must fit explicitly configured limits or raise.

These are bounds on the helper's explicit NumPy arrays, **not a total RSS guarantee**. Existing rhok_oo/metric tensors, prior live gradient buffers, allocator overhead and native-library workspaces remain outside the accounting. Matrix inputs and all edge-tile outputs are contiguous to avoid unnecessary adapter copies, but a BLAS implementation may still use private workspace. No real-size tensor was allocated or benchmarked.

## Small-array evidence

Eight tests passed in 0.026340 seconds. Forty deterministic comparisons cover C/C, C/F, F/C, F/F, strided, reversed, inner-transposed, read-only, aliased nonsymmetric and scaled signed arrays, at block sizes1,2,4,64. Native comparisons use the original PySCF contraction above its small-array fallback threshold; an independent `math.fsum` scalar reference checks the explicit i,j index mapping. Maximum absolute discrepancy was **1.42108547152e-14** versus native and **1.42108547152e-14** versus the scalar reference.

The scalar comparison uses a documented float64 forward-error envelope `8*eps*(K+1)*sum(abs(products)) + 8*eps`; native versus candidate is checked within twice that envelope. The maximum candidate/scalar error used 0.00253653 of the envelope. This is arithmetic validation, not a changed physical acceptance tolerance. An exact integer-valued nonsymmetric fixture rejects the incorrect unswapped `qij` mapping.

Allocation observers confirm precisely five explicit arrays, bounded finite masks and contiguous BLAS edge outputs. Negative tests cover rank/dtype/shape/budget failures before allocation, late-tile NaN/infinity, and overflowing results. No partial result is returned after a failure. The runtime reported one PySCF thread; thread environment variables were also set to one. Both inspected PySCF source files match the frozen runtime inventory.

Detailed retained results and source copies: `/Users/ashujo/.cache/dynamol-research/df-gradient-contraction-candidate-v1`. See `INTEGRATION-PLAN.md` for the still-required checks before any isolated future integration.
