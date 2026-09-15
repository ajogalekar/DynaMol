# 1MNC: first anchored catalytic-zinc optimization request

Prepared 2026-09-13. **No quantum calculation launched by this preparation. No fitted model, Hessian, full-complex preparation or physical validation is claimed.**

The selected research hypothesis holds exactly three native cap carbons at the source His218/222/228 alpha-carbon positions. The other **94 of 97 atoms remain free**, including every cap hydrogen, the zinc, all five donors and the entire **51-atom PLH inhibitor**. No metal–donor distance, angle or dihedral is restrained. The native unconstrained MCPB Gaussian/GAMESS requests are preserved separately as an alternative for a future model-sensitivity comparison.

## Why this boundary choice

The official MCPB workflow optimizes a small cluster before obtaining its same-level Hessian and explicitly requires inspecting the optimized coordination geometry. Our generated native request has no cap anchors. Its unconstrained gas-phase model could let the three cut histidine attachments move without the missing protein. [Official Amber MCPB tutorial](https://ambermd.org/tutorials/advanced/tutorial20/mcpbpy.php).

Anchoring truncated alpha-carbon positions has primary research precedent in metal-site cluster models. Moubarak and colleagues used such coordinate anchoring for a different Zn/Fe enzyme; their study does **not** validate a three-carbon-only 1MNC model or our method. Choosing only the three independently mapped carbons, while allowing the native cap-hydrogen guesses to relax, is our explicitly documented inference. [Moubarak et al., 2022](https://doi.org/10.3389/fmolb.2022.945415).

This remains an isolated gas-phase first-shell hypothesis. It omits nonlocal protein electrostatics, pocket sterics, solvent and potentially relevant hydrogen bonds. Free distal PLH motion must be inspected against the omitted protein in the original coordinate frame. Large distortions require a better environment/model; numerical convergence does not validate them. Structural Zn282 and Ca283 remain separate unfinished obligations for complete 1MNC preparation.

## Exact request

Artifacts are under `build/advanced-chemistry/metal/qm-1mnc-anchored-optimization-v1/`. The reproducible preparation script is `docs/audit/advanced-chemistry/prepare_1mnc_anchored_optimization.py`; it refuses an existing output directory and never calls a QM executable.

| Item | Frozen value |
|---|---|
| Electronic method | DF-RKS/B3LYPG; PySCF spherical 6-31G*; grid level 4; def2-universal-jkfit auxiliary basis |
| Chemical state | Zn(II), neutral HID/NE2 donors, O1-deprotonated PLH hydroxamate; cluster +1, closed-shell singlet |
| Initial coordinates | Exactly the original 97-atom native small model; no deletion, truncation of PLH or coordinate relaxation during preparation |
| Initial density | Verified completed same-geometry DF checkpoint, re-converged by the worker; not a transferred energy result |
| Frozen IDs | `cap:small:15:CH3`, `cap:small:54:CH3`, `cap:small:92:CH3` |
| Numerical SCF | Energy 1e-10 Hartree, orbital-gradient criterion 1e-7, at most 100 cycles |
| Optimization | At most 100 steps; energy 1e-6 Hartree; RMS/max gradient 3e-4/4.5e-4 Hartree/Bohr; RMS/max displacement 1.2e-3/1.8e-3 Angstrom |
| Resources | Two threads, 8,000 MB PySCF hint; six-hour native deadline and 21,720-second outer deadline |
| Derivative request | Final gradient only; no Hessian, ESP or force-field extraction |

The optimizer criteria and Cartesian constraint syntax are documented by [geomeTRIC](https://geometric.readthedocs.io/en/latest/how-it-works.html) and its [constraint documentation](https://geometric.readthedocs.io/en/latest/constraints.html). Its pinned local implementation calculates RMS/max atomic gradient norms after constraint projection. The audit additionally requires explicit free-atom and cap-atom Cartesian gradient reporting so a full gradient norm cannot be misread as constrained nonconvergence.

## Checks declared before results

The machine-readable `optimization-policy.json` is authoritative. Its geometry numbers are prospective engineering alarms, not experimentally justified accuracy bounds or fitting targets.

- Preserve all identities, the full named ligand bond graph, and the proposed proton arrangement. Both source stereocenters remain C3(R)/C9(S), checked by fixed-neighbor determinant signs and an independent final graph-based assignment.
- Report each of the five Zn distances and all ten Zn-centered angles. Distances outside 1.7–2.8 Angstrom, a new N/O atom within 2.8 Angstrom, or an angle change over 30 degrees requires model review. No generic minimum angle is imposed on the bidentate hydroxamate.
- Record all cap displacements and gradients, source-frame atomic displacements, protein-attachment geometry, whole and distal PLH changes, proper ligand torsions, and contacts against omitted protein atoms. The freeze tolerance is 1e-5 Bohr. Free atomic gradient RMS/max must separately satisfy 3e-4/4.5e-4 Hartree/Bohr; cap reaction forces remain visible.
- After a checked endpoint, obtain a conventional-integral gradient at **exactly the same geometry and electronic method**, using the DF density only as an initial guess. Preset DF/conventional gradient-difference RMS/max atomic-vector limits are 3e-5/1e-4 Hartree/Bohr. Preserve every raw difference; do not adjust energies to pass. Initial-geometry agreement does not establish endpoint or conformational accuracy.
- An alarm blocks automatic parameter extraction. Passing these checks still does not establish chemical accuracy. A later Hessian must respect that this is a constrained stationary point, including boundary reaction forces and the free subspace; ordinary full-space minimum/normal-mode claims would be invalid. No artificial anchor force constants may enter a fitted physical model.

Future launch requires root review, current PID/resource checks and the shared `QM-LAUNCH.lock`. At most three QM workers and two large DF jobs are permitted. Entry requires 22 GiB free with no other large DF job, or 35 GiB with one; the controller must stop its own child below an 8 GiB running reserve. Memory pressure, swap, hydration and time bounds must also be checked. These are allocation/resource controls, not an operating-system memory cap.

## Pinned handoff

- Input SHA-256: `3b6fc0c51200556bd78bf884760b77c7feeeb079d614ceee253480e9887a5cc5`
- Policy SHA-256: `61c06fa5679fbaafd7db06afb225b57b021985b7976fa3e7aca8a94da3ebd93f`
- Worker SHA-256: `ddb67a1f8e9be6ed05c3c95e9864e20b79f0233604e27926e30e5ffc09295bc1`

`source-manifest.json` binds the originating structure, CCD, atom/cap mapping, original unconstrained requests, completed DF artifacts, runtime declaration and optimizer implementation. `artifact-manifest.json` binds the prepared files. The unchanged initial donor distances are 1.980, 2.039, 2.051, 2.293 and 2.107 Angstrom for His218, His222, His228, PLH O1 and PLH O2 respectively.
