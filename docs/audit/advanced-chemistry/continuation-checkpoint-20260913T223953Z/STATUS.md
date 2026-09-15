# Advanced chemistry development status

User objective: carefully support metalloproteins (including metalloproteases)
and covalent protein–ligand adducts (including sotorasib), then test a fixed set
of 15 complexes in each group. Reliability and defensible parameterization take
priority over an immediate release. The standalone Mac window and release
packaging remain deferred. PatAgent is not used, following the user's request.

## Completed foundations

- Two original-input panels are frozen in `panel-freeze.json`, `covalent-panel.json`
  and `metal-panel.json`. See [the panel](PANEL.md). The original CIFs/CCDs,
  assemblies, endpoints, omissions and exclusions are retained.
- Primary-source method reviews: [covalent adducts](covalent-method-review.md)
  and [metal sites](metal-method-review.md).
- Actual 6OIM combined capped Cys–sotorasib product graph and a native Amber
  GAFF2 mechanics prototype, with source atom identities and observed geometry
  retained. This is not yet a validated parameter set for the app. High-penalty
  assigned torsions and ff14SB interface charges still need further work.
- Native MCPB array/ESP bridge checked against preserved 1OKL reference files.
  Matching native parser/parameter outputs tests the bridge; it does not validate
  a newly calculated metal site. Historical RESP-charge differences are retained
  as an unresolved comparison.
- Private PySCF 2.14.0 / geomeTRIC 1.1.1 runtime installed under `.tools/qm`.
  Initial computations use at most two threads and 8 GB per job on this 32 GB Mac.
- The neutral QM worker passed 41 small-molecule/protocol tests, including
  independent gradient/Hessian finite differences and electrostatic-potential
  integrals, verified same-geometry checkpoint initialization and explicit
  fixed-dihedral optimization with final periodic-angle checks. See
  [worker details](QM-WORKER.md). A larger calculation has not
  yet established an accepted new physical parameter model.
- The macOS RESP storage-only rebuild passed ten two-stage 113-atom reference
  fits (20 successful launches), matching every printed charge exactly. A second,
  constrained 94-atom synthetic fixture also matched both stages and every frozen
  charge. Source/build/validation records are retained under
  `build/advanced-chemistry/resp-runtime`; the installed Amber binary is unchanged.
- An immutable research snapshot loader preserves native topology, coordinates,
  chemical-state declarations, original structure and exact parameter manifests.
  Forty integrity/identity/round-trip tests pass after independent review. It
  explicitly cannot authorize app MD. See [production gaps](ATOM-MAPPING-REVIEW.md).
  The corrected 111-atom synthetic covalent polymer specimen also survived copy,
  reload and finite-force evaluation with identical System XML and coordinates.
  Its synthetic charges remain unsuitable for physical simulation.
- The full 4,534-atom 1CA2 snapshot survives copying/reopening with identical
  System XML and coordinates. Its source display mapping restores all 2,207
  original heavy-atom identities/coordinates and every native bond. Saved
  measurements now require their original mapping fingerprint; bare in-range
  indices cannot silently attach to another model.
- A separate `.tools/qm-parallel` runtime recompiles the same PySCF 2.14.0
  wheel C sources with OpenMP. The original active runtime is unchanged. All 31
  then-current native/protocol tests passed, and the full 94-atom fixed-geometry
  DF energy/gradient agrees with the serial runtime to 2.18e-11 Hartree and
  1.90e-12 Hartree/bohr. Two effective threads are now honored; concurrent-run
  timings do not establish a material speedup. Exact build/library hashes and
  numerical evidence are in `build/advanced-chemistry/qm-runtime-parallel`.
- Latest full backend checkpoint: **739 passed, 12 optional native QM tests
  skipped**, 9 expected fixture warnings. Native QM tests were run separately.
  Final affected identity/mechanics suites passed **63 tests** after review.
- A separate B3LYP-D3(BJ)/DZVP reference provider passed small-molecule energy,
  gradient and basis-identity qualification. The exact LibXC B3LYP convention,
  Psi4 basis source and two-body D3(BJ) definition are pinned. This qualifies the
  implementation, not the accuracy of an adduct model. Evidence is under
  `build/advanced-chemistry/covalent/reference-provider/qualification-v3`;
  [the provider checkpoint](COVALENT-REFERENCE-PROVIDER.md) includes launch-lock,
  resource-accounting and unfinished-parent rejection checks. No large DFT
  parent calculation or scan has launched.
- The isolated double-precision GROMACS comparison is complete. At unchanged
  coordinates and TPR, the one-versus-two-thread energy difference shrank from
  41.71875 kJ/mol to 2.56e-7 kJ/mol. Raw double-precision GROMACS minus OpenMM
  Reference energy is -0.408748 kJ/mol; a -0.457012 kJ/mol dispersion-tail
  component leaves +0.048263 kJ/mol unassigned. No offset was subtracted from
  the acceptance comparison. Force RMS/max differences are 0.00124345 /
  0.00706216 kJ/mol/nm. Exact scalar energy parity and NPT pressure equivalence
  are not claimed. See [the general contract](PERIODIC-VALIDATION.md) and
  `prototype-1ca2-periodic-diagnostic/double-precision-discriminator.md`.

## In progress

1. Actual sotorasib quantum charge/interface work. The full 94-atom DF-RHF
   optimization is running in `.tools/qm-parallel` with the five original frozen
   cap heavy atoms. The separate conventional reference calculation initialized
   from the accepted same-geometry DF density completed in 3,808 seconds.
   DF minus conventional energy is +0.0013747982 Hartree; maximum/RMS gradient
   component differences are 5.53992e-5 / 1.45855e-5 Hartree/bohr. These measure
   one geometry's approximation error, not relative conformer energies or
   physical force-field accuracy.
   The candidate remains conditional on exact-reference comparison and later
   geometry, charge and attachment-torsion checks. The original slow optimization
   was deliberately stopped with its evidence retained.
2. Final covalent parameter precedence and larger-peptide assembly. Independent
   Context/engine checks caught invalid zero-periodicity torsion overrides and
   protein-database defaults overwriting ligand terms. Earlier matching manifests
   did not catch these. Failed/superseded specimens are retained; canonical
   ff14SB and ligand/interface scopes are now separated by actual atom types,
   with no protein parmchk2 defaults overriding ligand terms. The final v4
   research specimens passed native backbone, Context and independent-engine
   comparisons; those are still not physical flexibility validation.
3. The 1MNC catalytic Zn prototype has 97-atom small and 139-atom large models,
   preserving its complete hydroxamate ligand and five declared coordination
   edges. The small-model conventional B3LYPG/6-31G* fixed-geometry gradient
   completed in 6,506 seconds with 22 SCF cycles and the declared +1 singlet state.
   Maximum gradient component is 0.120615 Hartree/bohr, so this experimental
   cluster geometry is not an optimized minimum. A same-geometry DF comparison
   completed with the exact accepted density used only as an initial guess.
   Native runtime was 735 seconds; DF minus conventional energy is
   -0.000487164787 Hartree, with maximum/RMS gradient component differences
   of 4.59547e-5 / 9.78536e-6 Hartree/bohr at identical ordered coordinates.
   This is an approximation measurement at one geometry, not an accuracy or
   conformational-energy validation.
   The original waiting controller was intentionally stopped before any child
   launched; its replacement preserves the scientific request and original
   deadline, permitting two large DF jobs only with at least 35 GiB free space.
   No Hessian or new physical parameters are claimed.
   The chosen deprotonated hydroxamate is an explicit model hypothesis.
   Structural Zn and Ca remain separate full-complex obligations.
   The next 97-atom optimization is now running with exactly three verified
   source-CA cap-carbon anchors; all remaining atoms, including Zn, its donors,
   cap hydrogens and complete PLH, remain free. This is an explicitly isolated
   constrained-cluster hypothesis. Its geometry, proton/stereo, pocket-contact,
   free-gradient and cap-force checks were frozen before launch; no Hessian is
   requested. See `prototype-1mnc-anchored-optimization-v1/README.md`.
   Its independent endpoint analyzer passed 37 synthetic checks, including
   free versus anchored gradients, both ligand stereocenters, all 96 covalent
   edges and contacts with the omitted protein. No actual endpoint exists yet.
   The analyzer will only write a conventional-gradient request after every
   frozen policy check passes; it never launches a Hessian or accepts a model.
   See `prototype-1mnc-anchored-optimization-v1/ENDPOINT-ANALYZER.md`.
4. The complete 1CA2 published ZAFF water-state model passed full native mechanics
   and static native/OpenMM comparisons. Amber electrostatic and near-pi torsion
   conventions were diagnosed against pinned source code; original raw failures
   and diagnostic copies are retained, with production parameters unchanged.
   Two 5-ps OpenMM seeds and one 5-ps GROMACS seed passed initial numerical/site
   smoke checks on 34,414 atoms. These were not equilibrated trajectories: the
   initial water-box density was only 0.8141 g/mL. A 75-ps GROMACS NPT warm-up
   completed, followed by 20 ps each of unrestrained OpenMM and GROMACS NVT from
   identical physical coordinates and box (1.02455 g/mL). Mean temperatures
   were 300.02/300.67 K, finite-force and broad coordination checks passed,
   and maximum aligned C-alpha RMSDs were 0.707/0.728 Angstrom. This is a short
   numerical/site check, not evidence of thermodynamic convergence.
   The periodic numerical discrepancy is substantially explained by precision
   and accumulation behavior in direct-space Coulomb bookkeeping. The completed
   double-precision discriminator and independent review are recorded above;
   original parameters and all raw failed results remain unchanged. The isolated
   numerical analyzer has 19 adversarial tests and no default tolerances or
   authority to approve a chemical model or release.
   Neutral bound water is an explicit model hypothesis; the hydroxide alternative
   has different published parameters.
5. Graph-only screening has constructed product candidates for **14/15 covalent
   cases**, including complete peptide-ligand expansion for 2H5I. Boronate remains
   a separate state/model obligation. These are not preparation or force-field
   successes. All 14 unique accepted graph chemistries now pass native GAFF2
   atom-type/lookup coverage screening, with no ATTN records or invalid torsion
   placeholders. This used zero placeholder charges and does not establish a
   physical model; analogy penalties reach 780. Eighteen graph/coverage contract
   tests pass. Evidence: `frozen15-graph-v4` and `frozen15-native-coverage-v2`
   under `build/advanced-chemistry/covalent`.
6. [The 15-metal method map](metal-panel-method-map.md) distinguishes published
   reusable sites from new QM/state obligations. Matched open-shell/oxyheme/SOD
   parameter sources are being investigated; no unmatched P450, Cu/Cu, bare-ion
   or different-redox model is substituted for the requested native site.
   Actual published His-ligated oxy-myoglobin parameter files were retrieved
   from Dryad with repository checksums and CC0 metadata. Their native baseline
   and updated-charge System definitions were assembled and reviewed separately;
   they are not ff14SB-compatible snippets. Complete 1A6M preparation remains
   blocked by a compatible additive sulfate model; original sulfates were retained.
   See `prototype-1a6m-charmm-source-audit/ASSEMBLY-PLAN.md`.
   A broader sulfate source audit found complete sulfate in the separate Drude
   package, but no verified compatible additive sulfate bundle. A 2026 primary
   study provides a relevant Cannon/CUFIX lead; its public repository lacks
   the required topology and pair corrections. Legacy CHARMM19 or Drude terms
   were not transplanted. See `prototype-1a6m-sulfate-followup-v1/VERDICT.md`.
   A separate 1HCK Mg/ATP transport prototype retains all original waters and
   examines published CMAP and pair corrections. Model transferability and its
   CC BY-NC source license remain explicit applicability/bundling limits.
   The full 5,216-atom native static discrepancy is now explained across three
   poses. Fresh Sander processes avoid repeated CMAP allocation failure;
   diagnostic copies matching documented Amber Coulomb and near-pi phase
   conventions agree to 9.22e-9 kJ/mol in energy and 6.45e-9 kJ/mol/nm in maximum
   force component. Raw failures and original models remain unchanged. This
   validates numerical transport of CMAP and special LJ terms, not their protein
   transferability. See `prototype-1hck-mg-atp/native-parity-v3/REVIEW.md`.
7. Full-system mechanics validation now rejects unvalidated custom/CMAP forces
   and requires constraints to match an independent protocol inventory. A valid
   physical bond plus an undeclared restraint or constraint cannot pass. The
   combined identity/snapshot/mechanics suite passed **80 tests** after these
   changes and an independent review finding: an undeclared virtual site could
   otherwise alter coordinates and force redistribution while passing a native
   term comparison. Unsupported virtual sites within the validated scope now
   reject; outside partial-scope sites remain visible as unvalidated context.
   Hashes and command are in `full-force-constraint-coverage-checkpoint-v2.json`;
   the preceding 77-test checkpoint is preserved.
   A prior partial test run was intentionally interrupted during File Provider
   import stalls and is preserved without a pass/fail claim. The native exporter
   also rejects unsupported Amber CMAP maps, with 13 focused guard tests passing.
8. Native Paramfit executed the two unmodified upstream NMA fitting controls.
   These are synthetic classical-energy fixtures, not new QM-fitted parameters.
   Their original post-analysis failed on a four-column output format and is
   retained. The [independent corrected review](paramfit-native-corrected-review/README.md)
   confirms the current upstream assertions for all 500/1,000 rows, with actual
   fitted-minus-target RMS errors of 0.000242033 / 0.000058798 kcal/mol and no
   recentering. It preserves sampling warnings, incompatible historical saved
   references and a source-confirmed misleading fourth-column header. No native
   calculation was rerun; this does not qualify a new molecular parameter model.
9. A separate [charge handoff adapter](COVALENT-CHARGE-HANDOFF.md) now validates
   the entire accepted parent/HF-ESP/RESP chain before native adduct assembly.
   It retains joint charges, all 18 canonical fixed charges and v4 parameter
   precedence. Twenty-eight contract tests passed, and a real attempt against
   pending RESP was correctly refused without a model. Actual charge-to-native
   execution awaits completion. The two candidate proper-torsion axes and
   separate training/held-out targets are declared before fitting; no new scan
   campaign is launched by that plan.
10. The separate DFT continuation is now waiting as controller **18420**.
    Eighteen tests passed after independent review fixed request-hash binding,
    repeated-cancellation cleanup and late-readiness deadline handling.
    It requires accepted parent geometry and the same-geometry RESP result
    before materializing one B3LYP-D3(BJ)/DZVP single point. The DFT calculation
    itself has not started; this is a mixed-level reference, not a torsion fit.
    See [continuation details](COVALENT-REFERENCE-CONTINUATION.md).
11. A separate GFN2-xTB energy/gradient adapter passed numerical qualification
    on small neutral H/C/N/O/F/S molecules. All 63 derivative components at two
    finite-difference steps passed, with maximum difference 1.40e-8 Hartree/bohr.
    Upstream regression, error, identity/transformation and one-thread controls
    passed independent review. The initial loose-SCC regression-protocol failure
    remains retained. A guard-only correction stops counting a test supervisor
    as another quantum worker; evaluated calculator functions are unchanged.
    The current interface admits at most 32 atoms; no full-adduct admission or
    physical-accuracy claim follows. See `covalent-geometry-provider/QUALIFICATION.md`.
    The separate constrained optimizer passed two native four-atom H2O2 checks
    (60 degrees and periodic -181/179 degrees) in 8.18 seconds total. Its
    62 contract tests include eight actual child-process cases, with cancellation,
    timeout, artifact failure and descendants all cleaned up. Independent review
    confirmed the numerical results and final lifecycle fixes; all failures and
    old sources remain preserved. Numerical optimization code was unchanged by
    the supervisor fixes. Current source, results and next admission gates are
    recorded in `LIVE-JOBS.json`; no full-parent geometry or scan was launched.
12. Auditing retained GDP exposed a generic validator bug: published signed
    periodic Fourier amplitudes were rejected by an inappropriate nonnegative
    bound. The correction preserves exact signs and phases, with all other
    mechanical checks unchanged. The focused suite passed 47 tests; independent
    review checked all eight negative GDP terms, twelve analytic energy/derivative
    controls and seven negative controls. The unmodified published GDP model's
    28-heavy-atom graph and all native parameter terms match. Its template
    coordinates have a severe close contact and are not a simulation starting
    structure; GDP/Mg model applicability and full 6OIM assembly remain pending.
    See `prototype-6oim-gdp-source-v1/REVIEW.md`.

## Still required before claiming app support

- Independent molecular-model quality checks, distinct from engine conversion
  agreement. Native bond survival alone does not establish physical accuracy.
- Complete topology and parameter preservation through prep, solvation, saving,
  reloading and simulation, including external bonds, exclusions and all interface
  terms. Metal formal oxidation state must remain separate from fitted partial charge.
- Clear UI review of the detected chemical site and selected model, with bounded
  background progress, cancellation and retained diagnostics.
- Native Amber/OpenMM comparisons and any supported GROMACS conversion checked
  without discarding terms or substituting approximate parameters silently.
- All 30 fixed cases tested and every failure retained; successful models then
  receive unrestrained motion checks with explicit solvent and repeated seeds.
- Updated documentation and package verification after these checks. The current
  downloaded app does not include this work.

This work has not yet established broad covalent or metal preparation support.
Fixed-topology MD models an intact chemical state; reaction kinetics, bond
formation/breaking, redox changes and coordination exchange require other methods.

## Continuing work and live calculations

A 30-minute thread follow-up named **DynaMol chemistry validation** is active
(`dynamol-chemistry-validation`). It continues useful work and reports meaningful
changes; it should not repeatedly report unchanged calculations. Pause it when
the agreed work is complete or the user asks to stop.

- [LIVE-JOBS.json](LIVE-JOBS.json) records the covalent optimizer and lightweight
  ESP/RESP continuation, output paths, bounds and next steps. Inspect actual
  PIDs/results before launching work. Conventional reference calculation is done.
- Root-owned 1MNC DF comparison v2 **completed**: former controller 96205 and
  worker 96452; `build/advanced-chemistry/metal/qm-1mnc-df-comparison-v2/result.json`.
  It started after disk space recovered above 60 GiB. Old controller 71211 was
  intentionally stopped while waiting; its frozen source and stop record remain
  in the v1 folder. Inspect actual PIDs/results before assuming either is active.
  Source script: [compare_1mnc_density_fit.py](compare_1mnc_density_fit.py).
  Six-hour maximum resource wait, two-hour native calculation plus independent
  7,320-second timeout, at most two compute threads and an 8 GB memory hint.
  Uses the shared QM launch lock, at most two large DF jobs and three QM jobs
  total. Entry reserve is 22 GiB without another DF job or 35 GiB with one; an
  8 GiB remaining-disk stop condition applies to this controller's own child.
- [METAL-STATUS.json](METAL-STATUS.json) and
  [warmed comparison](prototype-1ca2-warmed-comparison.json) record completed
  1CA2 simulations. No further metal dynamics are active. The independent
  periodic-energy review is in
  `prototype-1ca2-periodic-diagnostic/review-periodic-localization.md`.
- Active 1MNC anchored optimization: controller **12299**, worker **12304**;
  `build/advanced-chemistry/metal/qm-1mnc-anchored-optimization-v1/controller-progress.json`.
  Native six-hour limit and outer deadline 1789354347.475484; two threads,
  8,000 MB hint, at least 35 GiB entry reserve alongside the existing DF job,
  8 GiB running reserve. Twelve prepared-artifact and nineteen source hashes
  passed at launch. The controller/source copy and original unconstrained inputs
  are preserved. Inspect actual PIDs/results before starting another large job.
- DFT waiting controller **18420**: original eight-hour wait deadline
  1789363076.926211, at least 25 GiB before its one planned calculation,
  four-hour native limit plus 120-second outer margin and 8 GiB running reserve.
  Its persisted claim prevents duplicate jobs even after failure or cancellation.
  Current state is in `parent-dft-continuation-v1/progress.json` under the
  covalent build directory; verify actual state rather than assuming launch.
- Current large DF optimizer uses an approximately 12.6 GB temporary integral
  file; free disk space recovered from 11–22 GiB to 62 GiB, then was about 51 GiB
  with both jobs running. Do not delete active
  scratch or launch multiple large DF jobs simply because CPU slots are free.
  File Provider is also evicting small source/runtime files: observed dataless
  placeholders block read/import calls for minutes. Targeted Foundation download
  requests succeeded; exact storage provider is unconfirmed. This is separate
  from scientific computation failure. The active optimizer remains untouched.
- Practical integration finding: existing app workers explicitly run fixed-volume
  NVT without pressure equilibration. A validated density-equilibration stage is
  needed in the new production workflow; a short finite smoke run is insufficient.
- Packaging, publication, source-chemistry acceptance, full app lifecycle/UI
  integration and the final 30-case report remain unfinished. No advanced model
  has been promoted into app readiness merely because it serialized or ran.
