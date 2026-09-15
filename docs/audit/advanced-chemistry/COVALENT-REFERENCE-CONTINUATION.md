# Prepared parent DFT continuation

**Live update, 2026-09-13 21:17 UTC:** root reviewed the version 2 corrections
and started detached waiting controller **18420** at
`build/advanced-chemistry/covalent/parent-dft-continuation-v1`. Its original
eight-hour wait deadline is 1789363076.926211. The initial progress record
confirms `waiting_for_parent`, with no child calculation. The preparation
and review history below is preserved; do not repeat its launch command.
The resident source SHA is
`8124f88d0ed844a96d3768b4077500c624b17780d5385ddbff9c38f096e57768`.

The separate `continue_covalent_reference.py` controller is implemented and tested. **It has not been launched.** The qualified reference provider, original parent plan, pinned materializer, resident RHF optimizer, and existing RESP controller are unchanged.

The prepared configuration is `covalent-reference-continuation-plan.json`. A read-only preflight verified 48 pinned files and correctly returned `waiting_for_parent`, with only the existing two quantum workers present. It creates no claim, request, queue, or quantum output.

## Dependency and scientific contract

The controller waits for the actual parent optimization to be accepted and converged. It also requires the existing `resp-continuation-v1` to finish with `research_candidate`, binding its input, parent result, coordinate arrays, cap graph, ESP checkpoint/result/arrays, and fitting provenance to that same parent. A failed dependency or changed hash stops the continuation and preserves evidence. A `research_candidate` status does not imply physical validation.

Only the pinned `materialize_covalent_reference.py` creates the coordinate request; its existing identity, electronic-state, stereochemistry, frozen-cap, and bonded-geometry checks remain authoritative. The controller never guesses coordinates, bypasses failed checks, or substitutes the DFT calculation for HF ESP/RESP charges.

The one planned calculation is **B3LYP-D3(BJ)/DZVP at the accepted DF-RHF/6-31G* parent geometry**: a mixed-level single point, not a DFT-optimized structure. One geometry cannot establish torsion-profile or force-field accuracy.

## Lifetime and resource contract

- Combined dependency/resource wait: at most eight hours from controller start.
- Actual native calculation: 14,400 seconds; outer process-group timeout: 14,520 seconds.
- Two threads and an 8,000 MB provider memory hint, with common native thread environment variables also limited to two. The memory setting is not an operating-system memory cap.
- At least 25 GiB free before launch, at least 8 GiB while running, and at most three recognized quantum workers.
- The shared QM launch lock is tested and released before provider startup. The provider itself holds that lock during execution; the controller does not deadlock by holding it across launch. A resource race causes a retained failure rather than an automatic duplicate retry.
- Materialization and native execution each use their own process group. Timeout, low disk, and SIGINT/SIGTERM cancellation terminate that group, including descendants, with bounded escalation.
- A per-parent/method lifetime lock and a persisted exclusive claim prevent duplicates across output directories. Claims remain after success, failure, or interruption; restarting requires explicit review of the retained record. Existing output directories are never overwritten.
- Progress and final records include exact child PIDs/process-group IDs/commands, timestamps, source/config hashes, dependency hashes, and bounded resource observations. Native output is preserved.

The static preflight checks provider/runtime manifests and their referenced native-library/source files, the materializer and executable, qualification fixture results/arrays, and parent inputs. These checks repeat while waiting and before launch; completed output is checked against the unchanged dependency and input identities. The materialized request is hashed from its first completed byte snapshot. That original digest is required immediately before provider launch and again at completion, including the provider's actual input copy. A self-consistent replacement request cannot redefine the accepted geometry.

Process-group cleanup temporarily ignores repeated SIGINT/SIGTERM signals while termination, escalation and reaping complete; both handlers are changed/restored under an atomic signal mask. The wait deadline is checked before and after dependency/resource verification, before readiness can trigger materialization. Delayed readiness therefore cannot bypass an expired wait.

## Tests and artifacts

Thirteen tests use small local Python stubs, not quantum calculations. They cover a complete successful dependency → materializer → provider chain, native failure, bounded dependency wait, wrong input/result/ESP bindings, changed pending parent files, changed materializer, persisted/live duplicate claims, process-group timeout with a descendant, running disk reserve, and actual controller SIGTERM cancellation. A macOS teardown race discovered by the tests is preserved in the initial failure log and fixed by checking remaining executable process-group members instead of signal-0 during group teardown.

The original evidence remains in `build/advanced-chemistry/covalent/reference-continuation-review-v1/`. Independent review then identified the mutable-request digest, repeat-cancellation cleanup and late-readiness deadline cases. The expanded suite passed **18 tests**, adding prelaunch request mutation, a self-consistently changed request during provider startup, an altered provider input after the request is restored, repeated SIGTERM/SIGINT with termination-ignoring descendants, and readiness arriving after the wait expired. The current source/config/test hashes, passing log and fresh read-only preflight are recorded in `reference-continuation-review-v2/review.json`. The provider, materializer and configuration remain unchanged. No large reference calculation or continuation queue was started.

## Concrete commands for root review

Read-only preflight:

```sh
.venv/bin/python docs/audit/advanced-chemistry/continue_covalent_reference.py docs/audit/advanced-chemistry/covalent-reference-continuation-plan.json --preflight
```

After the root agent reviews the controller, start the durable continuation with a fresh output path:

```sh
.venv/bin/python docs/audit/advanced-chemistry/continue_covalent_reference.py docs/audit/advanced-chemistry/covalent-reference-continuation-plan.json --output build/advanced-chemistry/covalent/parent-dft-continuation-v1
```

These commands are relative to `/Users/ashujo/Documents/Science/DynaMol`. The second command waits for the actual parent and RESP prerequisites; it does not start a calculation merely because it was invoked. It should run under a detached local parent process if it must survive a terminal closing.
