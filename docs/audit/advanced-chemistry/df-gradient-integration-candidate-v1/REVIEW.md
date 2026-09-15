# Small native qualification of the tiled DF-gradient contraction

**Seven small full-gradient comparisons passed. This qualifies the isolated numerical integration at the tested settings, not large-molecule memory use or any covalent/metal force field.** No installed runtime, application source or active calculation was changed.

The new provider and complete evidence live at `/Users/ashujo/.cache/dynamol-research/df-gradient-integration-candidate-v1`. The first six-case plan was frozen before evaluation with SHA-256 `ffd703ebd1c91907abff62eb03f805b9eeb2bebc4528ef1c86fedba979983bd5`. A separately predeclared one-case supplement, SHA-256 `0b7cfe61b980d17f305bdd08705adbe339dd78f347b515dbcc62632e7b033742`, exercised the requested parent basis settings and native RKS grid-response default within the original 600-second allowance.

## Exact source and dispatch scope

The private `baseline_rhf.py` is byte-identical to pinned PySCF 2.14.0 `df/grad/rhf.py`, SHA-256 `a52525f269974884b28846b55d021c169d430048e0755b55557480f733128735`. The private candidate, SHA-256 `4480bdc1df639cadd7b0e89f86caac7943fc11765291c9b6e191b92431421df0`, adds one helper-bridge import and replaces only:

```python
lib.einsum('pij,qji->pq', rhok_oo[k], rhok_oo[l])
```

An AST restoration check proves that removing the added import and restoring that single call reproduces the original source AST exactly. The Apache source header is preserved. No global `einsum` replacement or installed-module mutation occurs.

The bridge verifies and loads the previously reviewed helper, SHA-256 `8d9e6bbb64e9e525ac2d41cb2783f699050e32c0516e89cb821a1229b1012060`, with auxiliary block size 64, 64 MiB explicit workspace and 256 MiB output limits. Its observations contain shapes, strides, dtype, arithmetic workspace plan and elapsed time; they do not copy full input tensors.

RHF gradients instantiate the private module's native gradient class. For hybrid RKS, a private subclass retains the installed native RKS `get_veff`, grid processing and `extra_force` methods, binding only `get_jk/get_j/get_k` to the private RHF module. The same construction is used for baseline and candidate. Runtime assertions check those exact function/global identities. A Python call observer records the existing auxiliary `extra_force` return path without replacing it.

Each candidate full gradient called the patched contraction exactly once. Actual operands were plain real float64 arrays of shapes `71×5×5`, `82×5×5`, `120×8×8` and, in the supplement, `188×8×8`. All cases exercised multiple auxiliary tiles and a tail. No full-adduct tensor was allocated. Nonzero auxiliary-basis response contributions were observed for every atom in both gradient paths.

## Numerical plan and results

The main six cases use three deliberately nonsymmetric, neutral, closed-shell structures: water, ammonia and formaldehyde. Each was evaluated with RHF and hybrid RKS/B3LYPG, orbital basis 6-31G, Weigend auxiliary basis, and auxiliary-basis response enabled. RKS uses grid level 3 with grid response enabled for finite-difference consistency. B3LYPG's observed range/hybrid coefficients are `(0, 0.2, 0.2)`.

The seventh case is formaldehyde with B3LYPG/6-31G*, `def2-universal-jkfit`, grid level 3 and the native **unchanged `grid_response=False` default**, verified by assertion rather than reassignment. Its auxiliary response remains enabled. It is a full-gradient parity check; it does not claim finite-difference qualification of the default grid-response approximation.

All SCFs use a 512 MB allocation hint, one PySCF thread, energy convergence `1e-12` Hartree, gradient convergence `1e-8` and at most 80 SCF iterations. Atom coordinates, element order, charge, spin, basis and settings are frozen in the plans. Positive gradients are in Hartree/Bohr; no sign or unit conversion occurs.

Predeclared full-gradient tolerances were maximum absolute error `1e-10` and L2 error `3e-10` Hartree/Bohr. **Both errors were zero for all seven comparisons.** Each pair shares one converged SCF state to isolate the changed contraction; central energy is therefore unchanged by construction. Orbital arrays and their hashes are retained, together with both complete gradients.

The six main cases also checked two representative atom-coordinate components at central-difference steps `0.001` and `0.0005` Bohr. The 48 displaced energies came from independent same-method SCFs. The 24 derivative estimates were compared against both gradient paths. Their predeclared maximum errors were `2e-6` Hartree/Bohr for RHF and `1e-5` for RKS; no tolerance was changed.

| Method | Largest error at 0.001 Bohr | Largest error at 0.0005 Bohr |
| --- | ---: | ---: |
| RHF | 5.1423741176e-7 | 1.2856192022e-7 |
| B3LYPG | 4.9974492633e-7 | 1.2391755710e-7 |

The roughly fourfold decrease on halving the step is consistent with central-difference truncation error in these fixtures. It is evidence for these calculated gradients, not an independent validation of a molecular force field.

A postprocessing-only exact-equality assertion on individual auxiliary contributions failed for the two water cases: the maximum auxiliary difference was `2.7105054312e-20` Hartree/Bohr while their final rounded full gradients were exactly equal. The initial script and failure are preserved in `summarize-initial.py` and `initial-summary-failure.json`. The summary now reports that measured difference explicitly. No native calculation, scientific setting, original result or predeclared acceptance threshold was changed or rerun.

## Resources and retained evidence

The controllers held the shared launch lock and allowed at most one of these native workers at a time. The sole other observed quantum worker was the existing PID 37089. All seven audit children exited successfully and were subsequently checked absent.

- Total work: **55 SCFs, 14 complete gradients, seven patched contractions**; no optimization or full-adduct calculation.
- Combined controller execution: **11.8511 seconds**. A conservative bound including the separately prepared supplement is **140.692 seconds**, below the original 600-second allowance. Each native worker completed in 0.516–2.132 seconds, below its 120-second limit.
- Largest process-tree RSS sampled every 50 ms: **367,640,576 bytes**. Largest self-reported process RSS high-water mark: **371,163,136 bytes**. Both are below the enforced 1,000,000,000-byte limit.
- PySCF and the actually loaded OpenMP library reported **one thread**. BLAS-related environment limits were also set to one, including `VECLIB_MAXIMUM_THREADS`; the loaded-library query exposed no separate BLAS thread-count getter. A stronger native BLAS thread observation was not obtained.
- Starting available memory was over 10 GB, with over 29 GB free disk. The guard required 2 GiB available memory, 25 GiB starting disk and an 8 GiB running disk reserve.
- Before launches, all **9,051 runtime files and six interpreter links** matched the pinned inventory. Installed Python/native code pins and all local provider pins matched again after completion.

PySCF's RKS VXC code retains its native `max(2000, ...)` block-budget floor despite the 512 MB user allocation hint. That hint is not an allocator or RSS bound. The sampled limit, process-group timeout and retained process high-water mark are the resource evidence for these small runs; none proves a bound for the 94-atom calculation. Timing comparisons between the baseline and candidate are also not a large-scale benchmark and include execution-order/cache effects.

All plans, source deltas, complete source copies, provider/runtime pin references, launch PIDs, native logs, finite-difference energies, full gradients, orbital arrays, memory samples, resource results and the postprocessing failure are retained. `summary.json` contains per-case details. `evidence-pointer.json` identifies the final immutable-file manifest and digest.

## Next gate

This small native integration should be independently reviewed before any new resource probe. A later separately bounded probe would need the actual parent input/state, source/runtime hashes and measured RSS under the unchanged 8e9-byte external cap. No such probe was started here. The installed worker, its journal behavior, model acceptance and release preparation remain outside this candidate.
