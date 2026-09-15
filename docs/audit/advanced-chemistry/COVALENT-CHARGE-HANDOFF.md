# Covalent charge and torsion handoff

The new `covalent_charge_handoff.py` adapter is a research gate between the completed sotorasib QM/RESP calculation and native Amber topology assembly. It does not enable advanced chemistry in the app and does not label numerical completion as a physically validated force field.

## Current evidence

- **28 contract tests passed** in `build/advanced-chemistry/covalent/charge-handoff-contract-tests-v3.log`. These cover pending/failed calculations, changed parent identity, missing/reordered/nonfinite charges, altered canonical charges and methyl equivalences, MOL2 atom aliases, and unexpected changes to bond, torsion, 1–4, LJ, or mass parameters. Charge vectors in these tests are explicitly synthetic bookkeeping fixtures.
- A real attempt against the still-running parent was refused before native assembly. Its unchanged failure record is `build/advanced-chemistry/covalent/charge-handoff-pending-refusal-v1/result.json`. No charge substitution or biological-model output occurred.
- The complete accepted-QM-to-native branch has **not run with actual sotorasib RESP charges yet**. Its checks and native charge transport still need that first real execution. The previous v4 native mechanics specimens use AM1-BCC or synthetic charges and remain separately labeled.

The initial test log v1 contains one failed case caused by capitalization in the expected error text. The actual rejection worked. Logs v1/v2/v3 remain available.

## Admission and native transport

`covalent-charge-handoff-plan.json` pins the capped molecular graph and SDF, original parent input, canonical library constraints, v4 native adduct supplement/reference topology, code, native TLEAP executable, and Amber parameter/library sources **before** final charges are available. This is an audit manifest for accidental changes and provenance; it is not a digital signature or independent proof that the QM model is physically adequate.

The adapter requires the exact terminal continuation, a converged and accepted DF-RHF parent, the parent-bound conventional RHF ESP calculation, matching expanded basis and electronic state, full atom identities, coordinate-array/checkpoint hashes, and matching fitted RESP artifacts. It recomputes the existing frozen-cap, molecular stereochemistry, and gross bond-geometry checks. The rounded native ESP text must reproduce the accepted QM ESP arrays. A final `accepted` flag by itself cannot supply this evidence.

All 94 joint charges are carried by their original source atom identities. The 18 fixed ff14SB cap/backbone charges are verified against the pinned canonical constraint inventory. The 12 temporary cap atoms are neutral together, so removing them retains the original integral charge without renormalization. Methyl H equivalence is checked; methylene H charges are not forced equal. The 82-atom residue keeps the complete reacted ligand and cysteine atoms. There is no independent ligand-only recharging, AM1-BCC fallback, zero-charge fallback, or second reaction of the deposited saturated graph.

The native MOL2 receives the accepted optimized coordinates and RESP charges. Its old AM1-BCC metadata becomes `USER_CHARGES`; its original topology and type records stay intact. Optional TLEAP assembly then verifies every charge, atom name/order, and optimized coordinate, and compares all non-charge native mechanics against the pinned v4 model. Canonical mechanics also undergo the independent ff14SB ACE–CYS–NME comparison. The v4 precedence remains native ff14SB plus the adduct-only GAFF supplement; protein `parmchk2` placeholders never override it.

The optional ACE–ALA–COV–GLY–NME construction verifies native polymer connections and every retained adduct charge. That scaffold is a topology/interface test, not the original 6OIM protein geometry. Full-complex assembly, all-force transport, and physical motion validation remain separate requirements. Every generated result keeps `simulation_ready: false` and `physical_parameter_acceptance: false`.

Every output directory must be new. An unfinished input produces a retained failure report, and later attempts use a new output directory. Pinned input hashes are checked again after assembly to detect concurrent changes. Empty or incomplete file-provider reads are rejected.

## Limited torsion fitting

`covalent-torsion-fitting-handoff.json` records two explicit candidate proper torsions:

1. Cys CB–SG–ligand C25–C24, about the S–C attachment.
2. SG–ligand C25–C24–C23, about the adjacent reacted linker.

Their current native type patterns each have one occurrence in the capped adduct. The complete inventory also shows the hydrogen-ended torsions about each axis, which are preserved. The starting plan permits at most six Fourier amplitudes across the two selected heavy-atom tuples, while charges, all other proper/improper terms, bond/angle terms, LJ parameters, masses, and 1–4 conventions remain fixed. A full-protein export must recheck the scope of each type pattern to prevent unintended parameter changes elsewhere.

Training and held-out angular offsets are declared separately before fitting. They are a **dataset contract, not a launched scan campaign**. The previously consulted six coarse diagnostics cannot be relabeled as unseen validation. The full-parent dispersion-corrected method, geometry quality, resource budget, convergence, and actual charge handoff must pass their gates first. Failed geometries remain in the denominator.

Each fitting point must carry graph/parent/charge hashes, stable atom IDs, exact functional/basis/dispersion provenance, coordinates and units, electronic/dispersion/total energies and gradients, and the unchanged native-MM evaluation at those same coordinates. Native Paramfit has only synthetic numerical qualification so far. The selected-term inventory and fit rank/conditioning must be verified before using it on this chemistry.

A common QM/MM energy-origin nuisance offset can be recorded during parameter fitting; it must never be used as an engine-conversion correction. Acceptance will require held-out energy and gradient/torque behavior, periodic continuity, plausible minima/barriers, separate conformational context, and full-complex unrestrained checks. Improved training error or a bond that remains connected is insufficient.
