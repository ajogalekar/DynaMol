# Independent review of the journal worker candidate

**Verdict: no blocking finding for the isolated candidate and its stated observer contract. Native numerical qualification is still required before use.** This review did not modify the candidate, application worker, running job or installed runtime.

The reviewed worker SHA-256 is `780b5bb968fa9d45987597be9412331d071e17ea38b2299749dc5e92be9f12b5`; the pinned journal helper is `b3a0543b5b2394373e7a787ad01ecd77865cc507f85bec8ec99efe6efed26528`. All 14 artifacts in the author's mock-audit manifest verified. Eleven actual PySCF 2.14.0 Python sources were copied and checked against the previously frozen local-runtime inventory `683b958788605f35ce2fc9a5b7ac94edb03f5316bd0b166cfcb415ffaea1f416`.

## Scanner family and method binding

`scf/hf.py:1612–1613` constructs the scanner with bases `(SCF_Scanner, mf.__class__)`. `lib/misc.py:1374–1376` obtains that scanner through `g.base.as_scanner()`. Thus the candidate's `isinstance(scanner.base, native_family)` check admits the actual dynamic scanner subclass while retaining the original method-family identity; it does not require an incorrect exact-type equality. The gradient scanner likewise inherits the concrete gradient class through `grad/rhf.py:300–302`.

Six independent cases constructed actual RHF and RKS objects, each without DF, with default DF and with an explicit auxiliary basis. They exercised the actual `SCF_Scanner.__call__`, `SCF_GradScanner.__call__` and `PySCFEngine.calc_new` dispatch using explicitly synthetic SCF and gradient kernels. Every case passed the family check and unchanged method fingerprint, including resolved orbital basis, ECP, auxiliary-basis declaration, functional/grid level and the configured SCF/resource settings. RKS/DF gradient factories retained their expected concrete gradient classes.

Native computational entry points were guarded to raise if reached. No SCF, integrals, DF tensor build, gradient computation or geometry optimization ran. These cases validate object and callback wiring, not any numerical energy or gradient.

## Coordinates, derivatives and convergence

The frozen source establishes the evaluation chain:

1. `geomopt/geometric_solver.py:70–96` reshapes incoming coordinates, sets the engine molecule in Bohr, obtains `(energy, gradients)` from the gradient scanner and calls the callback with those locals. The worker constructs nonsymmetric molecules, so this path does not symmetrize the coordinates. PySCF's convergence assertion follows the callback.
2. `grad/rhf.py:308–318` resets the gradient scanner to that molecule, evaluates its SCF scanner and returns the gradient kernel result. `scf/hf.py:1620–1650` resets the SCF scanner to that same molecule before its kernel.
3. `grad/rhf.py:458–465` assigns the full positive nuclear gradient to `self.de` and returns `self.de` after any native additions. `grad/rks.py` and both DF gradient subclasses inherit this kernel. The returned callback gradient and `scanner.de` consequently refer to the same result. The observer checks exact values; it neither negates nor rescales them.
4. `gto/mole.py:3125–3135` stores nonsymmetric ndarray geometry using the Bohr unit factor of one; `3258–3264` returns the stored Bohr coordinates. The six dispatcher cases confirmed exact coordinate equality across the callback, engine molecule, gradient scanner and SCF scanner, and confirmed that all three callback molecule references were the same object. They also confirmed `env['gradients'] is scanner.de` and preservation of the original mean-field geometry.

The default SCF kernel initializes and assigns Python `False`/`True` constants (`scf/hf.py:137,196,230,232`). The only other assignments invoke the optional custom convergence hook; that hook is `None` in all six configured object cases and is not installed by the worker. There is no evidence for a NumPy-boolean incompatibility in this path. A false convergence flag creates no journal evaluation; the original native assertion remains enabled. The helper's explicit-true requirement and disagreement rejection are appropriate for this contract.

## Persistence and unchanged calculation path

The 20 existing mock integration tests independently passed in 0.095 seconds. They cover stale coordinates/energy/gradient, changed method/basis/state, incomplete or reordered gradient indices, nonfinite values, convergence disagreement, source/native coordinate bindings, and failures before and after atomic publication. A persistence error reaches the worker's failure path, without accepted result arrays or a progress reference claiming a successful write.

An independent AST comparison reproduced the author's invariant: after removing only `_create_evaluation_journal` and the three exact initialization/callback/reference hook statements, the candidate is structurally identical to original worker `ddb67a1f8e9be6ed05c3c95e9864e20b79f0233604e27926e30e5ffc09295bc1`. Numerical methods, optimizer arguments, constraints, convergence thresholds, final checks and acceptance code are unchanged. The observer reads arrays and writes its journal; it does not modify native coordinates or derivatives and makes no evaluation calls.

All records remain `UNCONVERGED`, unaccepted and unauthorized for checkpoint reuse, including after a later successful optimization. A completed gradient evaluation may be an optimizer trial that was subsequently rejected. Neither a committed journal record nor these tests establishes an accepted optimization step, a usable checkpoint, an optimized minimum or chemical-model validity.

## Remaining limits and next gate

Stable atom IDs are a validated positional map; PySCF provides elements and state, not an independently tracked stable-ID array. A coordinated same-element permutation cannot be independently detected by this observer. Declared/native settings fingerprints do not attest loaded binaries or arbitrary unrecorded monkey patches. The helper's previously reviewed local-filesystem atomic publication and durability limits still apply.

This is a source and noncomputing contract review. Small, separately resource-approved native RHF and RKS qualifications should demonstrate the journal against actual converged gradients and exercise an induced journal failure before any adoption. No full-parent run, full-gradient memory qualification, active-job replacement or release acceptance follows from this verdict.

Evidence, synthetic journals, frozen sources and reproduction scripts are retained at `/Users/ashujo/.cache/dynamol-research/qm-journal-worker-independent-review-v1`. `evidence-pointer.json` identifies the final manifest and digest. Synthetic energy `-123.456` and synthetic derivative arrays are explicitly test values and have no molecular interpretation.
