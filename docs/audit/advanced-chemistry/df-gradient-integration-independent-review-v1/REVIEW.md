# Independent review of the small gradient integration

No blocker was found for the isolated numerical candidate at the tested settings.
This review ran no quantum calculations and made no active-runtime changes.

All 175 frozen artifacts matched their manifest. Reversing the one import and
one contraction substitution reproduces the original PySCF source byte for byte
and as an AST. Source inspection confirms separate gradient objects and copied
outputs: the matching arrays are not aliases of one reused result. The RKS path
retains native grid and auxiliary-response handling, with explicit private J/K
dispatch. Both paths use the same converged SCF state to isolate the substitution;
their central energy equality is therefore not an independent energy test.

Recalculation from retained arrays confirms seven exactly matching full-gradient
pairs. All 24 central-difference estimates, their errors and unchanged tolerances
were independently recomputed from the saved displaced energies. Orbital-array
hashes still match the initial shared state. Nonzero auxiliary contributions and
the expected single patched contraction are present in every case. The retained
postprocessing equality failure at 2.71e-20 was not hidden by rerunning a case or
changing the predeclared full-gradient tolerance.

The parent-basis supplement checks 6-31G*/def2-universal-jkfit with native RKS
grid-response defaults. It does not provide a finite-difference claim for that
default approximation. The main finite-difference cases use grid response.

These checks establish small-case numerical consistency, not physical force-field
accuracy, large-molecule peak RSS or speed. A separate bounded large-case probe
and its source/runtime admission are still required. The live optimizer must
remain unchanged.

Reproduction: `review.py` in this directory, using a NumPy-equipped Python.
Evidence: `/Users/ashujo/.cache/dynamol-research/reviews/df-gradient-integration-independent-v1/result.json`.
