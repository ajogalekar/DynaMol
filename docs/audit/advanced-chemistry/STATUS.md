# Advanced chemistry development status

User objective: carefully support metalloproteins (including metalloproteases)
and covalent protein–ligand adducts (including sotorasib), then test a fixed set
of 15 complexes in each group. On 14 September the user split development into
two parallel tracks: fix internal loops up to 12 residues per loop and prepare
v1 for release; continue advanced chemistry separately without making it a
release blocker. See [the current two-track plan](../v1-loop-release/PLAN.md).
Complete-complex support includes all required ligands, ions and cofactors;
isolated covalent-bond success does not qualify a metal-containing complex.
No advanced-chemistry prototype is promoted by this scope change. Packaging
preparation follows the loop release gates; publication is not yet authorized.
PatAgent is not used, following the user's request.

### Latest release-track update

The [fresh-source integration](../v1-loop-release/fresh-source-integration/REPORT.md)
now reaches complete common preparation from original 1UA2 input, including
all 12 missing residues, retained ATP/TPO, actual ATP AM1-BCC/SQM and hydrogen
assignment. Pinned physical parameters reload exactly after accounting for
the undirected representation of harmonic bonds. This capture intentionally
stopped before loop refinement/publication. It exposed preliminary sidechain
sampling moving an observed loop neighbor beyond compatible source/reference
caps. The narrow fix preserves those immediate neighbors until complete loop
refinement; other sidechain and metal handling stays intact. Fresh validation
and the frozen panel still gate app enablement, packaging and the video.
No loop model was newly accepted, no MD ran, and all capture children exited.

The [outer-torsion anchor comparison](../v1-loop-release/hybrid-fallback/OUTER-TORSION-ANCHORS.md)
now yields two private complete archived 1UA2 static passes from eleven fixed
refinements. Both model the same 12-residue gap and retain ATP/TPO. The general
anchor rule preserves every fully observed outside torsion definition exactly;
all twenty earlier reference selections are unchanged by the selection rule.
No quality or movement limit was relaxed. A real orchestration regression then
rejects four candidates, accepts the fifth and returns its final validated
coordinates with progress events. These are archived-complex checks, not
fresh-source application preparation, broad panel success or release readiness.
Independent audits confirm both passing saved structures and separately the
handoff endpoint, including unchanged ATP/TPO chemistry and exact physical
parameter equivalence after PDB reload. Repeated native refinement is not
bitwise identical; each endpoint is assessed separately.
Next integrate complete candidate validation and state-preserving coherent
handoff in an opt-in fresh preparation, then run the frozen panel. The app
default, packaging and queued video remain gated on that work.

On the separate chemistry track, [public DMP archive recovery](complete-complex-track/fit-artifact-recovery/red-db-followup/REPORT.md)
finds four candidate RESP charge graphs linked to the source cited by Panteva.
Neither those records nor the bounded [manuscript follow-up](complete-complex-track/fit-artifact-recovery/dupradeau-manuscript-followup/REPORT.md)
identifies the exact complete fitting input; another orientation variant has
the same citation. No parameter substitution or new simulation was performed.
The retrieval branch is complete pending a genuinely new primary artifact;
the complete-complex model remains unqualified and is not a v1 release blocker.

The subsequent [source-conditioned construction/replay](../v1-loop-release/hybrid-fallback/CONDITIONED-COMPLETE.md)
corrects omitted real boundary inputs in the research sampler and separately
tests native torsion-aware closure. Fifteen native attempts and 12 complete
refinements are finished with no admitted candidate. One 1UA2 output now passes
the modeled and newly defined boundary reference, but still alters a defining
coordinate of an existing outside outlier. Raw cap failures are localized to
the n-stem carbonyl. An intact control supports a concrete shifted-anchor
closure test next; it is not a target-protein validation. The release gates,
complete retained chemistry and queued video remain unchanged.

The approved [hybrid loop fallback comparison](../v1-loop-release/hybrid-fallback/COMPLETE-COMPARISON.md)
is complete. The bounded runner now retries after full candidate rejection,
preserves cancellation/integrity failures and returns only a verified final
artifact; 69 relevant software checks pass. Eight archived complete-complex
refinements finish with exact retained chemistry. One 12-residue native 1UA2
proposal passes the modeled backbone reference, but its affected joins still
fail, and no proposal is admitted. DiSGro remains a private evaluation with
unresolved inherited redistribution provenance. No app default, advanced
chemistry admission or publication changed. Next improve native candidate
construction at the observed joins, then repeat complete preparation gates.
No new MD/QM was run; all bounded native work is complete. The video remains
queued after loop validation.

### First results after the track split

The release-track [saved loop audit](../v1-loop-release/SAVED-BASELINE.md)
completed all 12 coordinate checks with intact source identities/hashes.
Only two of ten previously accepted main-panel outputs pass the added native
backbone reference when newly defined boundary torsions are included. The
earlier gross-geometry passes remain historical evidence, not a new release
qualification. Carbonyl-aware ranking passes software checks but is not a
sufficient construction fix; it remains an explicit research opt-in. Native
closure/refinement work and the final release audit are still required.

Two bounded 8K5R closure comparisons now isolate the next construction issue.
The transplanted refinement seed has malformed join geometry that torsion
sampling preserves. Starting from a coherent native fragment removes that
defect, but none of 36 extended-context closure solutions is within the
unchanged 1 Å observed-backbone cap; the closest maximum displacement is
4.010 Å. Three pass the local backbone reference, which is insufficient without
the other checks. Intact controls pass. See the
[coherent-seed report](../v1-loop-release/native-sampling/COHERENT-SEED.md).
The subsequent [bounded fragment pool](../v1-loop-release/native-sampling/FRAGMENT-POOL.md)
retains 40 coherent alternatives across five contexts and closes a frozen
ten-fragment subset. None of the 40 returned closures passes the observed
1 Å cap or the full local reference. The closest maximum displacement is
1.892 Å at an observed residue inside an extended span. A separate
[paired 1UA2 refinement](../v1-loop-release/CONFORMATION-REFINEMENT.md) shows
that preserving allowed seed torsions alone leaves a boundary outlier and
a 62.508° O–C–N join angle. Neither protocol is accepted. The next construction
work must account for observed internal context during generation, not only
at the final cap check. No failed candidate was accepted and no MD/QM was run.

The separate [complete KRAS coordination audit](complete-complex-track/README.md)
searches all 8,746 O/N atoms across 121 existing saved frames. The same six Mg
donors remain nearest, including GDP and retained waters; no additional donor
is within 3.2 Å. The GDP contact is still shortened relative to the deposited
pose, with a changed Ser17–Mg–GDP angle. This narrows the diagnostic question
without establishing an accurate solution-phase model. The next bounded task
checks every relevant solute/ion pair of a published Mg alternative against
the actual complete-complex topology. No new MD/QM or app admission occurred.

That [microMg applicability inventory](complete-complex-track/MICROMG-APPLICABILITY.md)
is now complete. GDP's coordinating O2B maps from scoped QR to original O3;
the pinned microMg overrides have no O3 entry. Matching Lennard–Jones values
does not make the differently charged RNA O2 atom interchangeable. Protein
and covalent-adduct type coverage is incomplete, and exact TIP3P oxygen
parameters also differ. All 29 published pair rows agree across the pinned
files. This is not an accepted replacement for the complete KRAS system.
That inventory measures published fitted-override coverage; an absent
override alone does not establish an absent mathematical model when a
documented default combining rule applies. No default or override was assigned
to the current model by that audit.

The [Panteva modified 12-6-4 assessment](complete-complex-track/PANTEVA-APPLICABILITY.md)
is now complete. It establishes documented source-rule coverage for 71 of
72 native types after recorded GDP alias resolution. GDP O3 and Ser OH have
defaults; four drug atoms of GAFF2 type c6 lack a source polarizability. The
remaining issues include transfer of fitted DMP/guanosine corrections to the
GDP/protein site, a consistent TIP4P-Ew water/counterion model, splitting shared
LJ type indices for selective corrections, and an exact GROMACS C4 path.
OpenMM's parser contains a C4 implementation; current ParmEd/GROMACS standard
conversion does not carry that term. These are source-level findings, not
executed candidate validation. Amber parameter data have explicit
public-domain provenance, separate from the article/SI and code licenses.
All 29 inspected inputs remain unchanged. Next resolve c6 provenance and GDP
transfer before building a separate static-test reference topology. Do not
tune a radius to one deposited pose or remove the nucleotide/metal.

The [focused c6/GDP source audit](complete-complex-track/c6-gdp-transfer/REPORT.md)
checks the current official Amber source. c6 is still absent from the C4 table
and no implemented/documented fallback was found. All 15 guanine-base charges
of retained GDP match the Amber GN reference within rounding precision;
this positive local evidence does not qualify fitted transfer to the charged
complete complex. Default O2/O3 phosphate coverage remains defined, while
the exact DMP fit-charge artifact remains unrecovered. All 26 inputs are
unchanged. No topology, parameter, MD/QM job or app admission resulted.

## Advanced-chemistry track: focused covalent validation

The user previously approved a practical candidate scope of known human KRAS, BTK and EGFR
cysteine drug adducts. Current cases are **6OIM–sotorasib, 5P9J–ibrutinib,
and 6JXT–osimertinib**. Unknown-chemistry QM controls are deferred to v2.
The original 15+15 panels remain preserved; no full-panel success is claimed.
This work now continues separately from the next v1 release milestone. The
historical results below retain their original artifact names, including `v1`.

### Latest verified results (14 September; supersede the earlier KRAS protocol qualification)

The corrected KRAS model completed its 100 ps OpenMM test, all 100 saved-frame
checks, and a verified portable restart. It then completed a 20 ps native
GROMACS continuation with 21 finite-force frames, preserved mapped stereochemistry
and no LINCS warnings. The exported and native topologies now both have 26,041
constraints; Reference-reader energy/force differences are 2.66e-5 kJ/mol and
1.42e-5 kJ/mol/nm. Evidence: `6oim-stability-v2` and
`6oim-gromacs-continuation-v2`.

The Mg–GDP concern persists: the corrected OpenMM run samples a mean Mg–O2B
distance of 1.887 Å; GROMACS samples 1.899 Å, versus 2.157 Å initially. Numerical
agreement has not resolved the physical-model question.

The original 6JXT loop-water collision is resolved in a new research candidate,
`6jxt-coherent-context-v3`. All 87 retained nonprotein heavy atoms now serve as
fixed steric obstacles during construction only. The previous GLY A873 O–water
D1226 O contact increases from 1.022 Å to 2.899 Å with the water unchanged.
Independent PDB-roundtrip/environment checks pass with no external soft overlaps,
unchanged unselected observed atoms, retained stereochemistry, and maximum
observed-context displacement 0.78935 Å. The rejected v2 candidate and its
failed environment screen remain preserved.

The screened coordinate intermediate `6jxt-protein-v5` assembled successfully
as `6jxt-complex-v1`: 42,900 solvated atoms, all four rebuilt residues, the full
covalent osimertinib adduct, original chloride and all 49 deposited waters.
Source identity mapping and complete native-loop geometry/stereochemistry pass;
the original covalent cysteine heavy atoms remain unchanged during assembly.
Temporary construction forces are absent from the native MD system. They are
geometric aids, not newly parameterized ligand/ion interactions.

The first monitored full-complex run, `6jxt-stability-v1`, is **rejected at
1 ps**: ALA750–THR751 peptide omega departed 36.934° from planarity, beyond
the unchanged 35° gross limit. Full-complex minimization passed but left that
peptide 32.299° from planarity, versus 7.280° in the assembled pose. Stereo,
clash and bond-angle checks passed; this is residual peptide strain, not a
return of the water collision. All failure evidence is preserved.

The separate `6jxt-stability-v2` test, which frees only the mapped remodeled
region during warm-up, is also rejected: the same peptide reaches a 38.383°
deviation at 7 ps. That change alone did not resolve the loop concern.

New read-only control checks found isolated deviations beyond 35° in all three
completed controls: 7/28,800 peptide-frame samples in 6JX0, 25/26,500 in BTK and
8/16,900 in corrected KRAS. A static preparation cutoff is not automatically an
appropriate per-frame MD acceptance criterion. However, the repaired EGFR peptide
is persistently more distorted than the corresponding control peptide during
the same first seven ps (mean 27.96° versus 2.65°). No failure was reclassified.

The **diagnostic-only** `6jxt-stability-diagnostic-v3` completed all 100 ps.
Independent saved-file checks verified 100 frames, final coordinates/box and a
portable restart with zero Reference energy/force difference. The restart and
audit preserve `diagnostic_only=true` and `qualified_for_handoff=false`.
Twenty saved frames have recorded loop peptide-geometry failures; completion
does not clear them. All 18 monitored peptide angles at all 100 saves agree
with an independent DCD replay within 0.00055 degrees.

At ALA750–THR751, mean departure from planarity is 30.02 degrees during heating,
35.03 degrees during restrained NPT and 11.20 degrees during the final 70 ps of
unrestrained NPT. The decline coincides with restraint release, but elapsed
relaxation time is a confounder. The matched diagnostic comparison
`6jxt-stability-diagnostic-v4` has now also completed 100 ps, releasing restraints
before heating with the same original inputs, native parameters, seed and
temperature/pressure schedule. All saved frames and its portable restart are
verified, retaining unqualified diagnostic status. There are 13 flagged
loop/context frames versus 20 in v3. At ALA750–THR751 the mean deviation improves
from 35.03 to 17.00 degrees during 11–30 ps, but worsens from 11.20 to 23.14 degrees
during 31–100 ps. Removing warm-up restraints is not accepted as a fix.

The native parameter comparison against 6JX0 matches exactly for all 24 atoms,
25 bonds, 47 angles, 137 proper torsion terms and six improper terms touching
ALA750/THR751, including neighboring endpoints. This rules out a difference in
those local terms between the two assembled models; it does not validate the
repaired coordinates or their surrounding environment. No MD or QM process is
currently active. Next assess the remodeled local conformation/environment and
the outstanding covalent charge/torsion and Mg–GDP model questions; do not simply
repeat the same loop diagnostic. Neither endpoint qualifies for a GROMACS
handoff or app admission.

A separate wholly unrestrained native minimization, `6jxt-native-relaxation-v1`,
did not improve the target peptide: its deviation is 33.41 degrees, and observed
context moved up to 4.39 Å from the assembled input. Its output is preserved,
not promoted to a preparation replacement or used for the comparison run.
Six native-frame admission checks and two GROMACS endpoint-refusal checks pass.
See `covalent-v1-focused/PEPTIDE-PLANARITY-DIAGNOSTIC.md` for evidence and limits.
No new QM, packaging or publication was introduced.

The next structural audit identifies a concrete refinement defect: THR751 is
favored in the coherent fragment seed but becomes a native CCTBX Ramachandran
outlier after refinement. The previously screened candidate has four selected
backbone outliers (THR751, SER752, ALA866 and GLY874). The 15-residue seed itself
has one PRO753 outlier. All 45 residue classifications and double-precision
scores across three saved structures match the complete independent CCTBX
PDB/ramalyze path; no statistical cutoff was fitted or relaxed.

A bounded conformation-protection experiment keeps THR751 favored but fails
peptide omega and two other backbone scores, and remains rejected. The research
exporter now checks the saved PDB with native CCTBX before creating an
intermediate: the old four-outlier candidate is correctly refused. Original
intermediates and diagnostics remain preserved and unqualified. No new MD/QM
run was started. Next compare bounded native fragment alternatives with
conformation quality assessed before and after refinement. Evidence and scope:
`covalent-v1-focused/RAMACHANDRAN-REFINEMENT.md`. The application and release
bundle have not received this research change.

A bounded seven-trial native fragment comparison has now completed. Four seeds
have no selected backbone-reference outliers, but no refined candidate passes
all geometry/reference checks. All rejected outputs and the frozen trial plan
are retained under `6jxt-fragment-alternatives-v1`; no MD or QM was launched.
The comparison exposes a narrower construction issue: temporary peptide-plane
restraints omit links solely between observed context residues that refinement
allows to move. SER752–PRO753 is omitted and fails in five trials. The native
physical force-field terms remain present. Next test construction coverage over
the explicitly mobile context while preserving original observed peptide
cis/trans basins, strengths and all acceptance limits. This is a testable
correction, not yet an accepted fix. See the updated Ramachandran report.

The context-peptide correction has now been implemented as an opt-in helper
option (off by default) and passes 27 analytical/software tests. A paired replay
of all seven saved seeds, with identical initial atoms and hydrogens, completed
14 refinements. None passes all final checks. Three corrected candidates pass
gross geometry but fail stereochemistry/reference checks; others retain bond
or peptide strain. No app workflow, native assembly or MD run accepts them.

A native KIC control reproduces an intact backbone, but a probe on one EGFR
seed exposes another boundary issue: matching the stem N/CA/C atoms leaves its
observed carbonyl oxygen displaced 2.0705 Å. Aligning that carbonyl orientation
before closure gives no solution for this seed's fixed internal geometry.
Next use native torsion sampling/closure with complete observed stem geometry
and full-atom candidate screening, rather than further restraint escalation.
Evidence is in `6jxt-context-peptide-replay-v1`, `6jxt-kic-feasibility-v2` and
`6jxt-kic-carbonyl-v1`; earlier failures and the frozen panels are preserved.
No MD or QM process is active, and release packaging remains deferred.

### Earlier handoff findings and implementation background

BTK and the 6JX0 EGFR control completed actual full-complex 20 ps unrestrained
GROMACS continuations, with finite saved forces and retained stereochemistry.
Native-value text precision fixed an export-rounding failure without fitting
parameters. See `covalent-v1-focused/HANDOFF-AND-LOOPS.md` for evidence.

The KRAS handoff caught missing element metadata after GDP atom-type scoping.
OpenMM had omitted 12 GDP hydrogen constraints. The old run remains preserved
but is superseded as a qualification of the intended constraint protocol.
Explicit native element declarations now preserve all 40 GDP elements in
`gdp-scoped-v6`; complete dry-model unconstrained energies/forces are unchanged.
The corrected `6oim-complex-v3` contains 27,477 atoms after fresh solvation and
has the corrected completed 100 ps run in `6oim-stability-v2`, audited above.
Do not perform a corrected KRAS handoff using the old endpoint or input hashes.

A separate coherent-context EGFR loop candidate now passes the existing local
geometry/stereochemistry checks within the unchanged 1 Å displacement cap.
It uses temporary construction forces at native backbone-angle targets;
`6jxt-coherent-context-v2` is not yet a full prepared complex or an app fix.
Its completed environment check found the collision recorded above; the new
v3 candidate and v5 intermediate supersede it for the active assembly/run.
Earlier failures,
physical attachment charge/torsion assessment and Mg–GDP model assessment remain
retained or unfinished. No app integration or publication was performed.

The native AM1-BCC/PREPGEN modified-residue route has now assembled each
adduct between two peptide neighbors: all three preserve their complete
adduct atom/bond inventory and have native parameters for both peptide
connections. Final charges agree with the declared neutral/neutral/+1 states
within 7e-6 e. These are actual native candidates, not synthetic RESP models.
They have not yet established physical flexibility accuracy or app readiness.

The BTK/EGFR charging check initially stopped on rounded SQM atomic charges.
Saved SQM molecular totals were correct, and native BCC conserved the rounded
atomic sums. Admission now checks that evidence explicitly, without modifying
charges; native PREPGEN's final charge check remains strict. Both candidates
subsequently assembled successfully without repeating the charge calculations.

The existing DynaMol ProMod3 workflow repaired KRAS's three-residue internal
loop and BTK's three-residue internal loop. Both complete complexes are now
assembled in native Amber topologies: 27,474 atoms for KRAS, including GDP,
Mg(II) and all original waters; 36,093 atoms for BTK. Source identity mapping
passes after solvation. GDP's source mechanics are isolated from protein
parameter names and preserved in the full KRAS assembly.

The KRAS and BTK 100-ps OpenMM tests completed under `6oim-stability-v1` and
`5p9j-stability-v1`. Each includes 70 ps of unrestrained NPT. Mean unrestrained
temperatures were 300.30 K (KRAS) and 300.22 K (BTK), with final densities
1.024 and 1.029 g/mL. Finite-force and mapped-stereochemistry checks passed.
Independent reads of all 100 saved frames per run agree with the final engine
coordinates/box within DCD precision. The drug attachment is unconstrained;
its heavy-atom torsion sampled spans of 51.3 degrees (KRAS) and 37.6 degrees
(BTK). These are numerical/site checks, not force-field accuracy validation.

KRAS retained its six initial magnesium donors, but the GDP oxygen contact
shortened from 2.157 angstrom to a sampled mean of 1.874 angstrom. This needs
physical assessment of the Mg/phosphate model; keeping the ion nearby is not
enough. BTK's stale final-PDB box metadata was corrected in the separately
verified `final-with-current-box.pdb`; the original export is preserved.
Portable restarts for all three completed runs are verified under each run's
`portable-restart-v2/`. They include Context parameters and a zero serialized
warm-up-restraint default. The original runner omitted parameters from State
XML and retained the initial nonzero System default, which could re-enable
restraints on reload. Reference-platform reloads of the corrected exports
preserve energy and forces exactly. Original artifacts and failed CPU-comparison
attempts remain preserved. The EGFR runner predates this export fix; its completed
output has now been corrected and independently reloaded in `portable-restart-v2`.

The original 6JXT EGFR adduct template is ready for research assembly, but its complete protein
still fails loop geometry. Four attempts, including experimental two-neighbor
relaxation and a restart from the original candidate, are preserved and rejected.
Next resolve consistent fragment/context
geometry without weakening acceptance checks. See the
[focused execution checkpoint](covalent-v1-focused/README.md) for results,
remaining work, and the corrected BTK export.

An additional EGFR control, 6JX0 from the same deposited study, has a verified
Cys797–osimertinib bond and no internal sequence gaps. Its protein intermediate
passed. Native charge and peptide assembly completed under `6jx0-baseline-v1`
and `6jx0-prepgen-v1`; the complete selected model has 53,791 atoms under
`6jx0-complex-v1`. Its 100-ps test completed in `6jx0-stability-v1`, with 70 ps
unrestrained NPT, mean unrestrained temperature 300.51 K, final density
1.016 g/mL and maximum aligned C-alpha RMSD 1.124 Å. The attachment ranged
from 1.763 to 1.935 Å without an attachment constraint. All 100 saved frames,
final coordinates/box and mapped stereochemistry passed the independent audit;
the attachment torsion sampled a 37.4-degree span. These are short numerical
and geometry checks, not physical force-field validation. This control does not
replace the rejected 6JXT loop test. The focused control retained the linked
drug, ions and waters and explicitly excluded the two other noncovalently
associated drug copies; the complete original CIF/environment remain preserved.

Current implementation: `covalent-v1-focused/`; actual native outputs:
`/Users/ashujo/.cache/dynamol-research/covalent-v1-focused/`.
Use the dedicated `covalent-v1-focused` cache runtime. Do not restart the
large 94-atom HF optimization or held 1MNC job merely because a slot is free.

Native GROMACS capped-peptide energies/forces now pass the fixed tolerances for
BTK, KRAS and the 6JX0 EGFR control on three coordinate probes each. The original
6JXT export exposed decimal-rounding sensitivity: preserving native LJ/charge
precision fixes the Reference-reader mismatch. Its mixed-precision native
comparison remains outside the fixed absolute tolerance at the high-force
source pose. Failed evidence remains retained; see
`covalent-v1-focused/CROSS-ENGINE.md`. Full-complex cross-engine checks remain.
The KRAS metal-site follow-up confirmed the actual loaded Li/Merz CM/TIP3P model
and pinned a published microMg implementation for compatibility assessment;
no ion or GDP parameters changed. See `covalent-v1-focused/MG-GDP-ASSESSMENT.md`.

## Earlier reassessment and retained diagnostics

The sotorasib 2,000 MB recovery **failed after 23 completed energy/gradient
observations** when its process-monitor `ps` command timed out. Both processes
exited; the optimization did not converge and no new charges were accepted.
The original run and diagnostics are retained. The prepared 1MNC job is held
because it uses the same supervisor.

Following the user's request for fresh thinking, the primary-source review is
complete: [covalent MD reassessment](covalent-literature-reassessment-v1/REVIEW.md).
The next milestone is an explicit protein charge/interface protocol using the
existing combined-adduct native baseline, followed by targeted assessment of
uncertain terms and complete-complex tests. A whole-adduct ab initio optimization
is no longer the automatic prerequisite for every case. No chemistry acceptance
criteria or app-readiness claims are relaxed by this research decision.

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
  Initial computations requested two threads and an 8,000 MB allocation hint;
  this was not a resident-memory cap. See the local-storage recovery below for
  the observed overrun, missing memory-accounting dependency and correction.
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
   optimization failed in its 2000 MB recovery process monitor after a previous
   checkpoint read timeout. It retained five frozen cap heavy atoms; no
   converged endpoint exists. See the current research direction above. The separate
   conventional reference calculation initialized
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
   The 97-atom optimization was stopped after an observed memory overrun and
   awaits a lower-allocation local recovery. It retains exactly three verified
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
10. The separate DFT continuation **18420 stopped after its parent failed**.
    Eighteen tests passed after independent review fixed request-hash binding,
    repeated-cancellation cleanup and late-readiness deadline handling.
    It requires accepted parent geometry and the same-geometry RESP result
    before materializing one B3LYP-D3(BJ)/DZVP single point. The DFT calculation
    itself has not started; this is a mixed-level reference, not a torsion fit.
    Recovery must explicitly bind the new accepted parent and local runtime;
    the old failed continuation and its persistent claim must remain preserved.
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
- The original 6OIM optimizer **failed** after 57 completed SCF/gradient
  evaluations: its HDF5 `scf.chk` read timed out with errno 60 while macOS marked
  the file `dataless`. Last proposed cycle 58 was not scored. Original results
  remain failed; the RESP and DFT continuations both stopped without launching
  downstream calculations. Do not restart those old controllers.
- The original 1MNC optimizer was **stopped through its controller** after its
  observed resident memory exceeded 10 GB, despite an 8,000 MB PySCF hint.
  Controller 12299 and worker 12304 both exited; no group members remain.
  Cycle 8 geometry/gradient are recoverable. Its checkpoint mixes uncompleted
  cycle 9 coordinates with cycle 8 SCF data, so it must not initialize a
  same-geometry density calculation. Original files and stop evidence remain in
  `build/advanced-chemistry/metal/1mnc-resource-stop-v1`.
- **Lower-memory local recovery trials were attempted and failed.** Jobs, inputs, worker copies,
  checkpoints and scratch use `~/.cache/dynamol-research`; the isolated runtime
  uses `~/.cache/dynamol-runtimes/qm-parallel-local-v1`. Original runtimes are
  unchanged. The copied native code payloads match; obsolete Documents RPATHs
  were removed only from copies. The added `psutil` dependency fixes macOS
  PySCF memory reporting, which previously silently returned zero. Water RHF
  energy/gradient comparison passed against the original qualified calculation.
- Recovery reuses independently checked last-completed geometries, with a cold
  SCF and unchanged method, state, atom identities and convergence criteria.
  Five sotorasib cap coordinates are restored to the exact original anchors
  (maximum correction 1.54e-8 A); all 89 free atoms retain cycle 57 positions.
  Native allocation is reduced to 5,000 MB. The new supervisor runs recoveries
  serially, samples RSS each second and terminates on an observed 8e9-byte
  crossing; this is explicitly a sampled monitor, not an OS allocation cap.
  Two threads, six hours plus a two-minute outer margin, 35 GiB entry reserve
  and 8 GiB running reserve remain enforced. Actual launch status is recorded
  in `LIVE-JOBS.json`; no recovered geometry is labeled converged.
  The independently reviewed supervisor passed **54 tests**. Recovery
  controller **34843** stopped worker **34848** after the first gradient
  crossed the monitored bound (8,035,172,352 RSS bytes). Both exited. The
  initial SCF converged; no gradient or geometry result was accepted. The
  allocation hint was still insufficient despite correct memory reporting;
  its 17 GiB scratch and all diagnostics remain retained. Source review found
  working blocks controlled by the hint plus a fixed tensor contraction whose
  original and copies can briefly total about 4.28 GB. DF and auxiliary-basis
  response dispatch are correct. This is a plausible contributor, not a sampled
  allocation-stack diagnosis.
  Nine pinned source files and two small layout/arithmetic checks support this
  review; including the metric derivative, the possible transient is 4.82 GB
  before surrounding live state. See `runtime-local-recovery-v1/6OIM-DF-GRADIENT-MEMORY.md`.
  Exactly one unchanged-method **2,000 MB allocation trial** ran as
  controller **37084**, worker **37089**, outside Documents. Both have exited
  after a process-monitor timeout. The trial retained the
  original 8e9-byte monitored stop, two threads, six-hour native limit and disk
  reserves; current bounds and results are in `LIVE-JOBS.json`. The lower hint
  can reduce surrounding blocks, not the fixed tensor spike. If it still fails,
  investigate the contraction rather than blindly repeating smaller hints.
  The first complete gradient of the 2,000 MB trial passed its predeclared
  near-geometry comparison against original cycle 57: energy difference
  4.10e-10 Hartree and maximum gradient-component difference 3.22e-8
  Hartree/bohr. Free coordinates match; restored anchors differ by at most
  1.54e-8 A. Sampled peak at that check was 6.97 GB. This establishes initial
  numerical consistency within the resource limit, not optimization convergence
  or model accuracy. Evidence: `~/.cache/dynamol-research/reviews/6oim-memory2000-first-gradient-v1`.
  A corresponding 2,000 MB 1MNC spec is prepared but unlaunched. Both former
  5,000 MB specs remain retained; see `METAL-STATUS.json` for the current queue.
- Practical integration finding: existing app workers explicitly run fixed-volume
  NVT without pressure equilibration. A validated density-equilibration stage is
  needed in the new production workflow; a short finite smoke run is insufficient.
- Packaging, publication, source-chemistry acceptance, full app lifecycle/UI
  integration and the final 30-case report remain unfinished. No advanced model
  has been promoted into app readiness merely because it serialized or ran.

## Recovery journal candidate

A separate `qm_evaluation_journal.py` helper passed 35 synthetic tests for exact
units, identities, state and provenance binding, completed-SCF checks, atomic
no-clobber snapshots and interrupted/corrupt writes. It records coordinates,
gradients and energy while explicitly retaining an unconverged, unaccepted
state. It does not infer that an HDF5 checkpoint matches those coordinates.
Independent review found and fixed recovery-time boolean/float identity and
schema confusion. Four new regressions fail against the retained old helper;
all malformed specimens now reject correctly. The v1 and v2 audits remain in
`qm-evaluation-journal-audit-v1` and `qm-evaluation-journal-audit-v2`. No remaining
blocker was found for this standalone bookkeeping candidate; it is not yet
integrated and no production or active worker was changed.

## Latest bounded candidate checks

The 2000 MB sotorasib recovery failed after 23 completed energy/gradient
evaluations, without optimization convergence. The terminal error was a
three-second `ps` timeout, not an observed RSS-limit crossing. Peak sampled RSS
was 7,991,640,064 bytes. Both processes exited. The zinc recovery is held for
monitor reliability; no additional large QM job launched.

A separate tiled density-fit contraction passed seven small native full-gradient
comparisons (RHF and B3LYPG), with exactly matching gradient arrays, and 24
finite-difference estimates. Independent review verified all 175 evidence files
and the exact source change. It retains the nonsymmetric contraction with no
omitted or approximated terms. Shape-only accounting gives 40,980,736 bytes of
explicit workspace plus 179,664,968 output bytes, excluding surrounding state
and BLAS workspace. Full94-atom RSS savings and speed remain untested.
See `df-gradient-integration-candidate-v1/REVIEW.md` and
`df-gradient-integration-independent-review-v1/REVIEW.md`. The candidate did not change the worker or runtime used by the failed run.

The journal-worker candidate passed two genuine three-atom constrained optimizer
pairs, ten cold single-point record replays and a first-record ENOSPC injection.
Original/candidate optimizer arrays match exactly; replay differences are at
most 4.27e-14 Hartree and 1.17e-9 Hartree/Bohr. Independent review verified all
200 artifacts and the earlier exact AST comparison preserves numerical code.
Every journal record remains explicitly unaccepted and does not establish
checkpoint reuse. Children were shorter than the one-second RSS sampling
interval, so these runs do not establish actual peak memory. See
`qm-journal-native-qualification-v1/RESULT.md` and
`qm-journal-native-independent-review-v1/REVIEW.md`. No active worker changed.

The queued 1MNC recovery endpoint wrapper now passes 75 software tests and
independent review. Review fixed strict terminal metadata and real-array
validation; the original crystal/pocket frame, chemistry functions and all
scientific thresholds remain unchanged. No terminal endpoint exists and no
downstream request was generated.
See `1MNC-RECOVERY-ENDPOINT-ADMISSION-V1.md`.

Preparation dependencies and the qualified RESP executable have separate
verified local runtimes outside Documents. Independent review checked the
preparation package inventory and all four relocated native RESP stage fits;
every reference charge file is byte-identical, with all 18 fixed charges exact.
RESP needs short relative filenames and complete parsed outputs: long input
names can truncate and return zero without results. This was a qualification
harness finding; the existing fitter already uses short names. Failed tests
remain retained. The new RESP fits peaked at 2.076 GB sampled RSS; future fits
must use reviewed lifecycle controls, not the one-off qualification harness.
These qualify local research runtimes, not portable packaging or new fitted
molecular parameters. See `covalent-prep-local-runtime-v1/REVIEW.md` and
`resp-local-runtime-v1/REVIEW.md`.

The recovery-aware ESP materializer now passes 37 software tests and independent
review. Its read-only core passes 18 tests and a separate review, with exact
reuse of the parent admission statements, complete binding of all five emitted
objects and array-header checks before allocation. Original type/parser failures
remain retained. The then-running, unfinished parent was correctly refused before native
reads. The cold conventional-HF branch never reuses the old checkpoint and
preserves the original charge constraints. See
`6oim-resp-recovery-materializer-v1/README.md` and
`6oim-cold-evidence-independent-review-v1/REVIEW.md`.

The separate cold-ESP launch wrapper passes 17 focused checks and independent
review. It uses the unchanged reviewed process supervisor, exact runtime and
artifact binding, isolated interpreter flags and a fresh bytecode-cache prefix.
Two tiny process tests check flags and timeout cleanup; these are not quantum
calculations. The enabled launcher correctly refused the actual unfinished
parent without a child. Its genuine full-adduct native path remains unexercised.
The separate cold native ESP/charge reader now passes 25 tests and independent
scoped review; actual fitting execution and handoff remain unfinished. See
`6oim-cold-esp-launch-v1/REVIEW.md` and `6oim-cold-native-reader-v2/REVIEW.md`.
Actual native fitting and charge handoff remain unfinished. No new ESP, RESP
parameter model or DFT calculation has launched.

A primary-source boronate follow-up found useful parameterization leads but no
complete admitted 5LF3 N-terminal Thr–bortezomib model. The original chemistry
remains retained and unresolved; see `boronate-source-followup-v1/REVIEW.md`.
The frozen 1CLL Ca/calmodulin case exposed an ethanol molecule missing its oxygen.
A separate narrow two-anchor repair candidate passed 27 checks and independent
review. Both observed carbons stay exact across 24 explicitly modeled alternatives;
the observed H21 record is retained as hydrogen. The original backend's three-anchor
gate remains unchanged. A hydrogen-admission revision passed 41 scoped checks.
The native named nine-atom AM1-BCC/GAFF2 ethanol component is now complete;
independent static comparison and complete Ca/protein/water preparation remain. Available
Ca/TIP3P source parameters do not establish EF-hand model accuracy. See
`1cll-two-anchor-repair-candidate-v1/REVIEW.md` and
`1cll-two-anchor-independent-review-v1/REVIEW.md`.

These are implementation and recovery milestones. Broad molecular validation,
the full 30-case preparation/motion panel and app integration remain unfinished.
