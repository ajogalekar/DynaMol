# 6OIM density-fitted gradient memory diagnosis

**Verdict:** A 2,000 MB allocation setting is a reasonable single resource-only trial, but it cannot bound this gradient below 8 GB. The pinned implementation has large allocations and temporary copies that ignore `max_memory`. Keep the external 8e9-byte RSS cap. No additional quantum calculation or runtime/source edit was made for this review.

## What the failed run actually established

The 5,000 MB run converged its first density-fitted RHF SCF to -2794.06939286907709 Ha. At gradient entry the native log reported 2,831 MB current memory. The supervisor terminated the owned process at **8,035,172,352 RSS bytes**, 267.003 seconds after launch. There is no completed gradient or optimizer callback. The result's `peak_sampled_rss_bytes` field (6,344,540,160) is stale; the actual terminating observation is retained inside its error message.

The case has 816 AOs, 4,739 auxiliary functions and 194 occupied orbitals. These are read from the saved checkpoint, retained DF tensor metadata and frozen auxiliary-basis data. The SCF DF tensor is already disk-backed at 12,637,440,576 file bytes; lowering the memory setting does not eliminate a previously in-memory full SCF tensor. The second 5.6 GiB scratch file has an unreadable HDF5 object header after termination. Its failure is retained; it was not repaired or treated as a complete gradient artifact.

## Dispatch and auxiliary response are correct

The log first says `Create scanner for <class 'pyscf.df.grad.rhf.Gradients'>` (line765). The [DF factory](/Users/ashujo/.cache/dynamol-runtimes/qm-parallel-local-v1/lib/python3.12/site-packages/pyscf/df/df_jk.py:194) selects the DF RHF gradient. [Scanner construction](/Users/ashujo/.cache/dynamol-runtimes/qm-parallel-local-v1/lib/python3.12/site-packages/pyscf/grad/rhf.py:301) retains the original gradient class as a base; [dynamic class metadata](/Users/ashujo/.cache/dynamol-runtimes/qm-parallel-local-v1/lib/python3.12/site-packages/pyscf/lib/misc.py:960) explains why the later label mentions `pyscf.grad.rhf`. The DF `get_jk`, `get_veff` and `extra_force` remain in dispatch. [Auxiliary-basis response defaults to true](/Users/ashujo/.cache/dynamol-runtimes/qm-parallel-local-v1/lib/python3.12/site-packages/pyscf/df/grad/rhf.py:542) and is included in both J/K derivatives and final per-atom contributions (lines565–578). Nothing in this worker disables it.

## Large allocations that the setting cannot shrink

| Allocation | Decimal MB | Source |
|---|---:|---|
| One rhok_oo occupied/auxiliary tensor | 1426.856 | `df/grad/rhf.py:175; full allocation, independent of max_memory` |
| Three coexisting original/reshape/Fortran tensors | 4280.568 | `df/grad/rhf.py:229 and lib/numpy_helper.py:270-273; 12-double layout discriminator confirms both copies` |
| Full auxiliary metric | 179.665 | `df/grad/rhf.py:487 and :426; factorization may coexist with original metric` |
| Full auxiliary metric derivative | 538.995 | `df/grad/rhf.py:224; independent of max_memory` |
| Full auxiliary contraction output | 179.665 | `df/grad/rhf.py:229; output allocated after operand conversion` |
| Two J/K derivative matrices | 31.961 | `df/grad/rhf.py:122-123` |
| ip1 AO derivative integral block at floor20 | 319.611 | `df/grad/rhf.py:126,134; lower floor remains even when available memory is negative` |
| ip1 AO-occupied transform at floor20 | 75.986 | `df/grad/rhf.py:143` |
| get_rhok block at floor20 | 25.329 | `df/grad/rhf.py:521` |
| Entire disk-backed rhok tensor | 6001.621 | `df/grad/rhf.py:492,511; stored in HDF5, retrieved in blocks; not all held in RAM by this function` |
| Entire disk-backed SCF DF tensor | 12637.434 | `df/df.py:167-197; this case already exceeds 5000 MB and 2000 MB allocation settings` |

The most consequential identified transient is [the `pij,qji->pq` auxiliary-response contraction](/Users/ashujo/.cache/dynamol-runtimes/qm-parallel-local-v1/lib/python3.12/site-packages/pyscf/df/grad/rhf.py:229). Its input `rhok_oo` occupies 1,426,856,032 bytes. The pinned runtime uses the default PySCF einsum backend. [The contraction implementation](/Users/ashujo/.cache/dynamol-runtimes/qm-parallel-local-v1/lib/python3.12/site-packages/pyscf/lib/numpy_helper.py:270) takes the reversed F-contiguous view, reshapes in NumPy's default C order, then converts back to Fortran order. A **12-double layout discriminator** shows that both steps allocate copies. The original, reshape copy and final Fortran copy can coexist: **4,280,568,096 bytes**, before adding the 538,994,904-byte metric derivative and other state. This is a source-supported possible transient, **not proof of the exact failed allocation**: the native log uses INFO level, so it lacks intermediate DF-gradient timers.

## What changing 5,000 to 2,000 MB can and cannot do

The worker assigns the requested setting to both `Mole` and the mean-field object. DF J/K and construction code, plus [gradient initialization](/Users/ashujo/.cache/dynamol-runtimes/qm-parallel-local-v1/lib/python3.12/site-packages/pyscf/grad/rhf.py:333), inherit it. The following block formulas use remaining memory `requested_MB - current_RSS_MB`: SCF tagged-density block (`df/df_jk.py:360`), metric solve AO partition (`df/grad/rhf.py:494–496`), AO derivative block (`:125–127`), and occupied-transformed retrieval (`:171–172`). At a hypothetical common 1,000 MB current RSS, lowering 5,000→2,000 reduces the respective targets from 225→56 auxiliary functions, 64.649→16.162 AOs, 125→31 auxiliary functions, and 1,579→394 auxiliary functions. Exact shell-balanced ranges can differ. Detailed scalar scenarios are in `analysis.json`.

The reduction is not universal. The derivative/retrieval formulas retain minimum blocks of20, and a nonpositive AO target partitions at individual shells (largest five AOs here). The metric/factorization, full occupied/auxiliary tensor, full metric derivative and contraction copies ignore this setting. In addition, the retrieval budget is chosen **before** allocating the full occupied/auxiliary tensor (`:171–175`), and the earlier AO-derivative partition is reused for the auxiliary derivative **after** that allocation (`:192`). Native workspaces, HDF5 buffers, earlier live arrays and allocator retention are not captured by a single scalar memory hint.

A lower setting can reduce the surrounding baseline enough to make the unchanged calculation fit, but only the bounded trial can establish that. If the single 2,000 MB trial also crosses8GB, the next step should be an equivalent bounded contraction/layout or disk-spilling review with numerical parity tests, rather than repeated hint reductions. Auxiliary derivatives, scientific method, atom set, source parameters and acceptance tolerances must remain intact.

## Evidence and limits

All nine inspected PySCF source/basis files match the frozen runtime inventory. Exact source copies and SHA-256 values are in `/Users/ashujo/.cache/dynamol-research/6oim-local-recovery-v1-memory-review-v1/source-manifest.json`; failed-run snapshots, HDF5 metadata/errors, basis counts, byte arithmetic and layout results accompany `/Users/ashujo/.cache/dynamol-research/6oim-local-recovery-v1-memory-review-v1/analysis.json`. No large arrays were allocated during this diagnosis. No quantum job was launched, and neither the failed scratch nor the running replacement was touched.
