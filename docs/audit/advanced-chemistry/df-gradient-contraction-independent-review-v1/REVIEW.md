# Independent tiled-contraction review

**Verdict: no blocker for the standalone small-array candidate.** This does not qualify full density-fitted gradients, an active runtime integration, peak RSS or full-scale performance.

Reviewed source SHA-256: `8d9e6bbb64e9e525ac2d41cb2783f699050e32c0516e89cb821a1229b1012060`.

The packing implements `C[p,q] = sum(i,j) A[p,i,j] B[q,j,i]`: each left tile is packed in `(i,j)` order, and the corresponding right tile is transposed to that order before multiplication. It preserves independent P/Q entries and makes no occupied-index symmetry assumption. The primary-source target call uses `rhok_oo[k]` and `rhok_oo[l]`, so admitting distinct nonsymmetric operands is appropriate.

I independently ran **45 tiny-array cases** with NumPy 2.4.6 and a scalar `math.fsum` reference. Cases cover rectangular occupied dimensions, singleton axes, C/F layouts, negative and zero strides, and simultaneous auxiliary tile tails. Maximum absolute difference was `1.7763568394002505e-15`, within the explicit floating-point envelope. The nonsymmetric rectangular discriminator matched its scalar result exactly and differed from the wrong untransposed flattening by 60. Inputs remained unchanged.

Independent allocation instrumentation found exactly the five planned arrays: output, two packed operand tiles, product tile and Boolean finite mask. Actual bytes matched the plan at its exact budget boundary. Finite checks reused that mask; left/right/product views shared their preallocated storage, and the final 2-by-1 product tile was contiguous. Both one-byte-under-budget cases rejected before allocation. Native-endian float64 restrictions, array-subclass/empty rejection, late NaN/Inf detection and finite-operand overflow refusal passed.

The real-shape **arithmetic-only** plan gives 40,980,736 bytes of explicit tile workspace plus 179,664,968 output bytes. Existing input tensors, allocator/Python overhead and BLAS/internal workspaces are excluded. No full-size tensor was allocated and no actual RSS saving is established by that arithmetic. The right tile is repacked for each left tile, so large-scale timing also remains open.

I verified all ten frozen author artifacts and reviewed the preserved eight-test/40-case evidence (maximum native/scalar difference `1.4210854715202004e-14`); I did not rerun that PySCF-importing suite. My checks imported no PySCF and ran no quantum calculation. No candidate, backend or runtime source was edited. Full-gradient dispatch/parity, finite differences and a bounded resource probe remain separate future qualifications.
