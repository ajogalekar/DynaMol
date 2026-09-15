# One-call shifted-anchor construction probe

**Five of six 1UA2 inputs closed with the shifted anchor, and all five now pass the original observed-backbone displacement cap and gross geometry checks. None passes the independent Ramachandran screen. All three masked-control inputs failed native closure and remain rejected.** No complete-complex acceptance is claimed.

The plan included all six original-span 1UA2 candidates and all three masked controls from `boundary-conditioned-native-v1`, with no selection by quality outcomes. Each received exactly one native CCD closure call using its saved seed and the native maximum of 1,000 iterations. No Monte Carlo draws, additional retries, geometry fitting after closure, new forces, refinement, MD or QM occurred.

## Outcomes

The named protocol identifies the saved input to the shifted-anchor step; every new closure uses native torsion-aware CCD. Displacements include all observed N/CA/C/O atoms in the original span.

| Case / saved input | Seed | Native closure | Largest observed displacement (Å) | Gross geometry | Modeled reference outliers | Affected boundary outliers |
|---|---:|---|---:|---|---|---|
| 1UA2 torsion-source-dirty | 2026 | Converged | 0.191 | Pass | LYS44, LEU45, SER49, GLU50, ALA51, LYS52 | ASN56, ARG57 |
| 1UA2 torsion-source-dirty | 2027 | Converged | 0.235 | Pass | GLY54 | ARG57 |
| 1UA2 torsion-source-ccd | 2026 | Failed | 0.537 | Fail | ALA51, LYS52 | ARG57 |
| 1UA2 torsion-source-ccd | 2027 | Converged | 0.160 | Pass | LYS52 | ARG57 |
| 1UA2 fragment-source-ccd | 2026 | Converged | 0.121 | Pass | LYS52 | ARG57 |
| 1UA2 fragment-source-ccd | 2027 | Converged | 0.170 | Pass | ASP53 | ARG57 |
| Masked control torsion-source-dirty | 2026 | Failed | 0.567 | Pass | None | None |
| Masked control torsion-source-ccd | 2026 | Failed | 0.376 | Pass | None | None |
| Masked control fragment-source-ccd | 2026 | Failed | 0.476 | Pass | None | None |

The converged 1UA2 inputs previously moved ILE43 O by 1.022–2.162 Å. After coherent orientation locking, ILE43 O is within 0.094–0.100 Å of the source. That small remaining difference reflects incoming versus source frame/carbonyl geometry; this procedure does not copy the oxygen coordinate independently.

The five successful native calls have c-stem N/CA/C RMSD 0.0984–0.0999 Å, satisfying native CCD's 0.1 Å criterion. The failed 1UA2 call remains at 0.2394 Å and creates an ASN56 C–ARG57 N distance of 1.145 Å, which also fails the unchanged gross join check. The three control calls remain at 0.4715, 0.1268 and 0.4137 Å respectively. Although those control outputs pass coarse geometry and the local reference, they do **not** satisfy native closure and are not candidates for admission.

ARG57 was already an outlier in the original source; its defining preceding carbonyl changes here, so it remains an affected boundary rather than a wholly unchanged source outlier. Its interpretation follows the existing source-boundary audit. No threshold or exemption was introduced.

## Construction and integrity checks

The saved full backbone is loaded by exact source atom identity. Its first N/CA/C frame is fitted to the observed stem, then first psi is rotated with `sequential=True`, moving the carbonyl and downstream chain coherently. A two-residue temporary entity supplies the first modeled residue and its preceding context to the native closer. The suffix closes against the original far-end stem and its original next residue; recombination copies the closed suffix without a new transform.

For all nine attempts:

- The locked prefix N/CA/C/O and virtual first-modeled N/CA/C coordinates remain exactly unchanged during suffix closure.
- The first join's C–N length, CA–C–N/O–C–N/C–N–CA angles and omega remain exactly unchanged during closure. The largest carbonyl-orientation residual is 0.00000572 radians.
- Original observed-only omega basins remain preserved in the independently evaluated region.
- Before/aligned/locked/recombined coordinates, atom identities, temporary-anchor maps and native success/failure are retained, including failed calls.

All five converged outputs completed native sidechain reconstruction. Their exported backbone matches the raw closure output exactly, and no exported observed heavy-atom change falls outside the declared context. TPO/ligand chemistry remains bound to the unchanged original full source by reference; these local checks do not validate full-complex parameterization, contacts, rotamers or energetics.

The frozen source, candidate and implementation hashes verified before and after generation/screening. All prior artifacts remained unchanged. No native API/type correction or rerun was required. The two supervised native batches took 0.944 seconds combined, with peak sampled RSS 232,980,480 bytes, below the fixed 120-second/4-GiB/two-thread limits.

## Candidate handoff

The five converged candidates, with modeled heavy atoms and declared context moves, are:

- [1UA2 torsion-source-dirty, seed 2026](/Users/ashujo/.cache/dynamol-research/v1-loop-release/shifted-anchor-native-v1/1UA2/native/torsion-source-dirty-seed-2026/candidate.json)
- [1UA2 torsion-source-dirty, seed 2027](/Users/ashujo/.cache/dynamol-research/v1-loop-release/shifted-anchor-native-v1/1UA2/native/torsion-source-dirty-seed-2027/candidate.json)
- [1UA2 torsion-source-ccd, seed 2027](/Users/ashujo/.cache/dynamol-research/v1-loop-release/shifted-anchor-native-v1/1UA2/native/torsion-source-ccd-seed-2027/candidate.json)
- [1UA2 fragment-source-ccd, seed 2026](/Users/ashujo/.cache/dynamol-research/v1-loop-release/shifted-anchor-native-v1/1UA2/native/fragment-source-ccd-seed-2026/candidate.json)
- [1UA2 fragment-source-ccd, seed 2027](/Users/ashujo/.cache/dynamol-research/v1-loop-release/shifted-anchor-native-v1/1UA2/native/fragment-source-ccd-seed-2027/candidate.json)

Machine-readable evidence: [shifted-summary.json](/Users/ashujo/Documents/Science/DynaMol/docs/audit/v1-loop-release/hybrid-fallback/shifted-anchor/shifted-summary.json). New implementation files are `native_probe.py`, `run_probe.py` and `screen_probe.py` in this directory. The cache `/Users/ashujo/.cache/dynamol-research/v1-loop-release/shifted-anchor-native-v1` holds the frozen plan, implementation copies, input hashes, supervisor records and all outputs.

Plan SHA256: `bc8d4c1cccfcba2350c046d37d174654327037224c31eb810bc56e51ee81f4b9`. Independent summary SHA256: `2589f4c73f5ec0748bba0ce34928eaae64ef83baa60f446a0eb7d55e810a96b8`.

Every output was independently screened with the existing CCTBX and geometry implementation. Failed closures have separately named diagnostic coordinate files and explicit native-failure flags; a diagnostic geometry pass is not promoted to success. The next complete-complex validation is parent-owned. This probe establishes a local geometric improvement on these five outputs, not a validated loop-repair workflow or release readiness.
