# Journal worker candidate

This is an isolated copy of `backend/qm_worker.py`, with a pinned copy of the reviewed evaluation-journal helper. It is **not integrated into the application or an active runtime**, and has no native QM qualification yet.

The only original worker changes initialize an evaluation journal before optimization, invoke its observer first in the existing callback, and attach a committed record reference to progress. The observer makes no energy, gradient or optimization calls. The calculation, optimizer arguments, thresholds, constraints, final SCF/property checks and acceptance path remain unchanged; the AST audit verifies this after removing those three statements.

The observer requires explicit native SCF convergence. A false flag writes no evaluation record and leaves PySCF's following convergence assertion intact. For a completed evaluation it checks current engine, gradient-scanner and SCF-scanner coordinates in Bohr, element order/state, full gradient atom order, native gradient units, returned energy/gradient consistency, and current method/basis/settings. It then records the full positive gradient and coordinates without conversion. The input unit, PySCF Bohr conversion constant, normalized element mapping, original resolved input, native initial geometry and method fingerprints are bound in the immutable manifest.

Each record stays UNCONVERGED and unaccepted, including after the whole worker eventually succeeds. A callback may represent a rejected optimizer trial; it is not an accepted step or optimized minimum. An I/O failure propagates through the existing worker failure path, with no accepted result arrays. The journal never reads or validates an SCF checkpoint.

Stable atom IDs are inherited from the validated positional map. PySCF exposes element order, not a separate stable-ID map; a same-element permutation cannot be independently identified by this observer. Method fingerprints are observations of declared/native settings, not loaded-binary attestation. Atomic publication and durability retain the reviewed helper's local-filesystem limitations.

Run the mock integration tests with `python3 test_integration.py -v`. They import no NumPy or PySCF and perform no native calculations. `run_mock_audit.py` records a fresh audit with source hashes, an exact AST comparison, test output and an explicitly synthetic evaluation specimen. The helper remains byte-identical to the independently reviewed 35-test version.

Before any later use, this candidate still needs independent review and resource-approved native qualification on small closed-shell RHF and RKS examples, including convergence, coordinate/state/method binding and an induced journal failure. No full-parent calculation or active-run replacement is authorized by these mock results.
