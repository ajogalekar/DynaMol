# Independent endpoint analyzer

`analyze_1mnc_anchored_endpoint.py` implements the frozen `optimization-policy.json` without editing the request, worker or policy. It does not run quantum chemistry. The actual optimization endpoint was not present when its readiness check ran; no partial trajectory or raw-crystal geometry was substituted.

The analyzer first verifies pinned request/policy/manifests, atom/state identity, actual constraint file, numerical convergence, resolved method, runtime declaration, units and checksums. It rejects mismatched, missing, nonfinite or unrequested arrays. Source hashes are evidence identity checks; they do not authenticate a loaded binary or prove a calculation's chemistry is accurate.

Its chemical-geometry diagnostics cover:

- All 97 gradient vectors, with full/free/frozen-cap statistics and distinct atomic-vector versus Cartesian-component RMS values. Cap reaction forces remain visible and are never zeroed to force a pass.
- Exactly five named Zn distances and ten Zn-centered angles, plus every other N/O atom near zinc.
- All 96 intramolecular covalent edges and 50 named hydrogen parents. HID connectivity comes from the pinned native OpenMM Amber ff14SB residue template; cap bonds come from the preserved MCPB mapping. No template charges or numerical bonded parameters are transplanted by this analysis.
- Both PLH stereocenters using signed neighbor determinants and RDKit reassignment from actual coordinates. Inherited chiral tags/CIP annotations are cleared first. The fixed ligand graph, bond orders, formal charges and atom order are independently checked before reassignment.
- Every source-frame atomic displacement, cap-pair distance, Zn/donor/core displacement, whole/distal ligand displacement and all 36 ligand heavy-atom proper torsions. No ligand alignment hides translation or rotation relative to the pocket.
- Contacts against 1,212 omitted original protein heavy atoms, using label-chain IDs, author residue numbers, explicit HIS/HID identity mapping and unchanged coordinates. Every retained heavy atom has a nearest omitted-protein contact; the detailed list includes contacts within 4 Angstrom in either geometry. This 4 Angstrom reporting window is not an acceptance threshold. The frozen 1 Angstrom collision alarm is unchanged.

Hydrogen-parent tests are geometric checks of the declared state. They do not establish protonation equilibria or infer electronic bond orders. A zero-alarm result remains a constrained-cluster numerical specimen, not an accepted physical force field or full preparation.

## Verification

The synthetic test suite passed **37 tests**. Its positive fixture uses explicitly artificial zero gradients and a deliberately non-native test checkpoint. It is never passed to a QM executable and never counted as a real endpoint. Tests cover cap forces versus free forces, maximum versus RMS criteria, anchor movement, donor loss, new coordination, preserved-distance angle distortion, proton transfer, stereo inversion despite inherited source tags, intracluster/pocket collisions, ligand drift, wrong units/identities/methods, incomplete/nonfinite arrays, hidden constraints and downstream-request blocking.

The first 30-test log remains intact; the expanded 37-test log is `build/advanced-chemistry/metal/1mnc-endpoint-analyzer-tests-v2.log`. Both are software verification only. `endpoint-analyzer-audit.json` records their hashes and the actual readiness result.

## Handoff after the worker finishes

Run using the application analysis runtime, against the completed `optimization` directory, with a new audit output directory. Output creation refuses to overwrite existing evidence.

```sh
.venv/bin/python docs/audit/advanced-chemistry/analyze_1mnc_anchored_endpoint.py \
  --endpoint build/advanced-chemistry/metal/qm-1mnc-anchored-optimization-v1/optimization \
  --output build/advanced-chemistry/metal/1mnc-endpoint-review-v1 \
  --write-conventional-input
```

The optional conventional-gradient input is generated only if every declared endpoint check passes. It uses the actual final Bohr coordinate array exactly, preserves the method/state, disables density fitting, and binds the completed endpoint density/checkpoint as an initial guess. Any alarm or evidence failure prevents that input from being written. The analyzer never launches it; the same resource/controller review remains necessary. Convergence and the later conventional check still do not authorize a Hessian or physical parameter extraction automatically.
