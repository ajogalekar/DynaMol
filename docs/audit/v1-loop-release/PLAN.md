# DynaMol: two development tracks

Direction authorized on 14 September 2026. This supersedes the earlier plan
that made advanced covalent chemistry part of the next release milestone.

**Status, September 2026 (v1.0.0):** the user chose to publish v1.0.0 with
the existing default loop path and an explicit statement of its
conformational limitation, rather than hold the release for Track A's
integration. The Track A candidate-search workflow below stays research code
under this directory; enabling it in the preparation worker is the planned
loop change after v1.0.0. The [released-source recheck](../loop-fallback/REPORT.md#v100-recheck-on-the-released-source)
and [LOOP_REPAIR.md](../../LOOP_REPAIR.md) record what shipped. Track B is
unchanged by the release.

## Track A: loop repair, then v1 readiness

The release scope is the existing viewer, analysis, engine and supported
metal-free standard/noncovalent preparation functionality, plus a validated loop repair
fix. The size target is **up to 12 missing residues per internal loop**, as
specified earlier by the user, not a limit of twelve separate loops. The
existing 96-residue total construction budget remains unchanged. Missing
terminal extensions and longer or unsupported gaps remain explicit limits.

1. Audit the saved successful preparations in the frozen 15-structure loop
   panel, including 1UA2 and 8K5R, with an independent backbone reference.
   Distinguish reconstructed residues and actually moved context from
   unchanged pre-existing structural outliers. Preserve the original verdicts
   and artifacts as historical evidence.
2. Correct the construction/refinement defects exposed by that audit.
   Recent research identifies a concrete boundary issue: matching observed
   N/CA/C coordinates can still lose the observed carbonyl orientation.
   Candidate selection must respect complete join geometry; increasing
   construction restraints or loosening acceptance limits is not a fix.
3. Re-run eligible complete-complex preparation cases, including both
   12-residue examples and multiple-loop cases. Check exact sequence and atom
   identities, retained molecules, connectivity, chirality, backbone
   conformation, environment collisions, observed-coordinate displacement,
   cancellation, progress reporting and saved-file round trips. Unsupported
   advanced chemistry stays clearly distinguished from loop-search failures.
4. After those gates pass, run the focused UI/engine regressions and build a
   reviewable v1 candidate with accurate support documentation. Verify the
   exact candidate's bundled runtimes and installation before recommending
   publication. No new GitHub publication is authorized by this plan.
5. After the 12-residue loop work is validated, update the existing inaugural
   product video to show the finished v1 workflow. Preserve its brief, polished
   style and prior requested highlights: bundled OpenMM/GROMACS; file input
   and PDB fetch; protein preparation including loop building; monomer
   selection; display modes with atom-type coloring for sticks/ball-and-stick;
   the visible water box; actual protein motion; and interactive distance or
   dihedral picking with the trajectory moving during measurement analysis.
   Use crisp, legible captions and capture actual working behavior. Preserve
   the earlier video and produce an updated local preview/export for review.
   Advanced chemistry claims must reflect the released scope. This video
   update is queued after loop validation, not an immediate detour.

Modeled loops remain provisional coordinate models. Software tests and
geometry checks do not establish the native loop conformation or an accurate
conformational ensemble. Do not advertise universal recovery of every gap
within the size limit.

## Track B: complete covalent and metal-containing complexes

Continue implementation and validation independently of Track A. The frozen
15 covalent and 15 metal-containing panels remain unchanged, as do all failed
and superseded artifacts. Current focused human-target candidates are
KRAS–sotorasib, BTK–ibrutinib and EGFR–osimertinib, with independent controls.

Support applies to **the complete complex**: required ligands, ions and
cofactors must remain present with a defensible model. For KRAS this includes
GDP and magnesium. A working attachment bond with an unsupported metal site
does not qualify the complete system. Numerical engine agreement and short
stability runs are useful checks, not physical-model validation.

This track is not a v1 release blocker. Experimental parameters and admission
paths stay out of the release until separately qualified. Unknown-chemistry
QM controls remain deferred. Do not restart the failed 94-atom HF or held
1MNC calculation; no paid services or PatAgent. Existing resource bounds
remain in force.

## Current work

- Active priority: [Abl–imatinib short-gap follow-up](fresh-source-integration/ABL-IMATINIB.md).
  The first targeted complete static pass repairs both 2HYY gaps with imatinib
  retained and maximum observed-heavy movement 0.767 Å. It combines better
  placed fragments with temporary preservation of initially accepted torsions.
  General automatic selection, independent result review and fresh full-worker
  qualification remain pending; the app fallback is still disabled. Earlier
  automated panel counts and failed candidates are preserved below.
- The [fresh representative panel](fresh-source-integration/PANEL.md) is complete:
  five metal-free static passes (both 12-residue examples), four unresolved
  multiple-gap cases, five expected unsupported-input checks and one deferred
  ATP–magnesium case. Metals are explicitly outside v1. The first sweep took
  5 minutes 17 seconds, excluding later debugging. Two native handling bugs
  were corrected without changing final quality or movement limits. Fresh
  full-worker 1UA2 and 1M17 publication handoffs both passed privately; the app fallback
  remains disabled and all failed candidates remain failed. Target public v1
  is tomorrow morning, contingent on honest scope and remaining app/package
  checks. The video follows a working, validated loop workflow.
- Latest Track A milestone: the [outer-torsion anchor comparison](hybrid-fallback/OUTER-TORSION-ANCHORS.md)
  produces two complete archived 1UA2 static passes from eleven fixed
  refinements. Both model the same 12-residue gap, retain ATP/TPO and preserve
  all fully observed outside torsion definitions exactly. This is a general
  connectivity-based anchor rule, with unchanged quality and movement gates.
  The new rule reproduces all 20 earlier reference selections; failed cases
  remain failed. A separate real orchestration regression rejects four
  candidates, accepts the fifth and returns the final validated coordinate
  artifact with progress events. This is not a fresh preparation or broad
  panel qualification, and the app default remains unchanged.
  Independent saved-artifact audits confirm both original passes and the
  separately refined accepted handoff endpoint. The latter has small seed-H
  variation and different final loop-sidechain coordinates, so no bitwise
  native reproducibility is claimed; its own checks pass.
  Next follow the [integration review](hybrid-fallback/integration-review/REVIEW.md):
  extract complete attempt preparation/validation after the frozen scaffold,
  reuse pinned chemistry/state bundles, supervise native work, bind all output
  sidecars, then exercise fresh-source 1UA2 through the actual preparation
  worker. Do not substitute more sampling of the same successful 1UA2 proposals
  for this integration step. After that, run the frozen multi-length panel,
  including all 8K5R gaps and the second 12-residue case, before enabling the
  automatic UI path. The following entries preserve the earlier comparisons.
- Track A: the [saved-output baseline](SAVED-BASELINE.md) is complete. Five of
  ten previously accepted panel outputs pass the modeled-residue reference;
  only two of ten also pass newly defined boundary torsions. Both additional
  8K5R repeats fail. All 12 audits pass identity/hash and serialization checks.
  A [carbonyl-ranking comparison](CARBONYL-RANKING.md) and 42 software checks
  are complete; ranking alone is insufficient and remains research-only.
  Two bounded native-closure comparisons on the smallest failing 8K5R gap are
  complete. The transplanted refinement seed contains malformed joins that
  torsion-only closure cannot repair. A coherent native fragment removes those
  inherited defects, but its closed alternatives still exceed the unchanged
  1 Å observed-backbone cap. No acceptable protein candidate resulted.
  The [bounded multi-fragment comparison](native-sampling/FRAGMENT-POOL.md)
  is now complete: 40 coherent raw fragments and 40 closure outputs from a
  frozen ten-fragment subset produce no combined pass. The closest closed
  result moves an observed internal context atom 1.892 Å and has reference
  outliers. A separate [paired 1UA2 refinement](CONFORMATION-REFINEMENT.md)
  confirms that protecting allowed seed torsions alone trades modeled
  outliers for a distorted join; neither protocol passes.
  The [saved endpoint audit](native-sampling/KIC-ENDPOINT-AUDIT.md) verifies
  158 returned solutions, including repeated controls: all actual closure
  endpoints stay within the cap. Excess motion occurs in observed residues
  made internal by extending the modeled span. No endpoint targeting,
  source-identity or unit-conversion error was found. All bounded audits are
  complete; no native sampling, refinement or MD/QM process remains active.
  Next make observed internal context geometry part of generation, rather
  than only an endpoint target and final rejection check, then apply full-atom
  environment/refinement checks and the reference gate. Do not repeat
  sampling of the malformed transplanted seed, simply increase torsion
  protection, or increase the displacement cap.
  Root owns construction changes and release gates; no loop fix is yet
  declared passed.
  The user's proposed short/long-loop split prompted an
  [alternative-method review](ALTERNATIVE-METHODS.md). Prefer automatic
  fragment-first repair with fallback after complete candidate rejection,
  using length to guide budgets rather than as the only routing rule. The
  existing native Monte Carlo stage currently runs only for unclosed gaps;
  late quality failures do not reach it. Prototype that routing correction
  and compare unused native fragment proposals. DiSGro/pyDisgro is the first
  external benchmark candidate, subject to implementation and inherited-data
  provenance review; no new dependency is selected for the release.
  The approved [hybrid fallback prototype and complete-complex comparison](hybrid-fallback/COMPLETE-COMPARISON.md)
  are now complete: 69 software checks pass, and eight frozen full-complex
  refinements demonstrate actual retry after a closed candidate is rejected.
  One native 1UA2 fragment proposal has zero modeled-residue reference
  outliers after refinement, but affected joins still fail; no candidate
  passes all gates. All alternatives meet the original observed-heavy cap
  after refinement. The DiSGro alternative remains private and unbundled:
  its inherited distribution rights are unresolved and its sampler omits
  nonprotein chemistry. Do not enable the runner in the app yet. Next make
  boundary torsions and complete observed join geometry part of native
  candidate construction, retaining the existing source-outlier distinction
  documented in [the boundary audit](hybrid-fallback/SOURCE-BOUNDARIES.md).
  All new native comparisons are complete with no remaining processes.
  The next [source-conditioned comparison](hybrid-fallback/CONDITIONED-COMPLETE.md)
  supplies the real flanking torsions and residue identities, which the
  installed default MC pipeline omits. Fifteen fixed native attempts yield
  nine candidates; all three masked controls pass locally. Twelve complete
  1UA2 refinements still yield no accepted model. One now has no modeled or
  newly defined boundary outliers, but changes the defining coordinate of an
  existing ARG57 outlier. All six raw 1UA2 movement failures are localized to
  ILE43 O. An [intact-control design check](hybrid-fallback/CARBONYL-LOCK-DESIGN.md)
  supports moving the numerical closure anchor to the first modeled residue
  after coherently locking the observed n-stem carbonyl. Next test this on a
  frozen small set of saved coherent proposals, then full-complex refinement;
  no target-protein result or new acceptance is claimed for that design yet.
  Any partial-anchor refinement must preserve complete outside torsion
  definitions, with actual coordinate checks and no exemption for affected
  source outliers. All completed comparisons retain the same gates.
- Track B: the [complete KRAS coordination audit](../advanced-chemistry/complete-complex-track/README.md)
  is complete. All required components and the six saved-frame Mg donors are
  retained; the GDP-contact contraction and angular differences remain a
  physical-model question. The [microMg compatibility inventory](../advanced-chemistry/complete-complex-track/MICROMG-APPLICABILITY.md)
  is also complete. Absence of a fitted pair override is not by itself absence
  of a model when a documented default rule exists. The next published
  [modified 12-6-4 assessment](../advanced-chemistry/complete-complex-track/PANTEVA-APPLICABILITY.md)
  finds default-rule coverage for 71 of 72 native types after documented GDP
  alias resolution; four drug atoms of type c6 lack a source polarizability.
  GDP O3 and Ser OH defaults exist, but fitted-correction transfer, changing
  all water to TIP4P-Ew, selective atom-type indexing and exact GROMACS C4
  transport remain unresolved. Source-code support in OpenMM is not an
  executed complete-complex validation. Next resolve parameter provenance
  and transfer before constructing a separate static-test topology.
  No fresh MD/QM run was required for these diagnostics.
  The [c6/GDP transfer audit](../advanced-chemistry/complete-complex-track/c6-gdp-transfer/REPORT.md)
  now confirms that the current official C4 source table still lacks c6 and
  has no documented fallback. All 15 retained GDP guanine-base charges match
  the Amber guanosine reference within rounding precision, supporting local
  similarity but not whole-GDP fitted transfer. Default phosphate O2/O3
  coverage remains defined; the exact DMP fit-charge artifact is unrecovered.
  No new parameter was assigned. Next recover the explicit carbon crosswalk
  and fit topology/charge evidence before constructing an alternative model.
  The bounded [public DMP artifact recovery](../advanced-chemistry/complete-complex-track/fit-artifact-recovery/red-db-followup/REPORT.md)
  found four candidate charge graphs tied to Dupradeau 2010, but no complete
  uniquely identified Panteva input. A [primary-source follow-up](../advanced-chemistry/complete-complex-track/fit-artifact-recovery/dupradeau-manuscript-followup/REPORT.md)
  found another orientation variant under the same citation; access failures
  prevent concluding what the unavailable manuscript/SI contains. End that
  retrieval branch until a new primary artifact is available. No charge/type
  substitution, complete-complex parameter admission or new MD/QM resulted.
- Shared release/chemistry follow-up now uses this plan. Root owns shared
  status/checkpoint updates; each worker writes new evidence in its own
  directory.
