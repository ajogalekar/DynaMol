# Missing-loop repair

**Release status (v1.0.0, September 2026):** the loop builder described
below is enabled in the app. On the released source, the frozen
[15-structure panel](audit/loop-fallback/REPORT.md#v100-recheck-on-the-released-source)
prepared **7 of the 9 metal-free repair candidates** through the real
preparation worker, including both 12-residue examples (1UA2, 1M17) and the
four-gap 8K5R; each accepted model passed the panel's 20 software, identity,
geometry and static force-field checks. 2CG9 and 3HEG were rejected by the
acceptance gates, and the remaining panel entries stay blocked by covalent
chemistry, a ligand reference conflict, the 12-residue limit or the deferred
metal case.

**Known limitation.** An independent Ramachandran check of the seven accepted
models with native CCTBX reference tables finds backbone outliers among
modeled residues in **4 of 7** (1UA2: 3 of 12 residues, 8K5R: 5 of 17, 1M17:
1 of 12, 2HYY: 1 of 4); 1FPU, 4HHY and 4HJO have none. The geometry, chirality
and collision gates below contain no backbone-conformation reference, so an
accepted loop can be sterically clean yet conformationally implausible. Treat
rebuilt loops as rough starting coordinates: inspect them, minimize and
equilibrate before analysis, and prefer a supplied complete model for regions
that matter to the question being asked. The earlier
[saved-output audit](audit/v1-loop-release/SAVED-BASELINE.md) reached the same
conclusion on the previous panel outputs.

A stricter workflow (bounded candidate search, preservation of qualified seed
torsions during refinement, and a CCTBX backbone-reference gate) passes 6 of
the 9 metal-free panel cases on all eight complete checks in private research
runs; see the [fresh-source panel](audit/v1-loop-release/fresh-source-integration/PANEL.md).
That code is in this repository under `docs/audit/v1-loop-release/` as
research scripts and is **not enabled in the app**: it depends on a separate
CCTBX runtime that the app does not bundle and on hard-coded local paths.
Integrating it into the preparation worker is the planned loop change after
v1.0.0.

Enable **Build supported missing loops / residues** before Prep. DynaMol uses
the original PDB/mmCIF sequence to model internal gaps of up to **12 standard
amino acids per gap**, with a local work budget of **96 residues across all
chains**. Thus repeated copies of a 12-residue gap can be prepared together.
Longer gaps and missing modified residues need a supplied complete model.
Unresolved terminal extensions are listed but remain omitted.

Rebuilt regions are listed in the preparation summary as provisional models.
Missing coordinates are not recoverable with experimental certainty from the
sequence alone. A passing geometry check does not establish the native pose.

## Method

[ProMod3](https://openstructure.org/promod3/) searches its local structural
fragment database and builds sidechains. The pipeline first searches with
the original stems, then checks sparse fragment matches. If a gap remains,
it retries with one and then two neighboring residues in its **temporary
construction context**. Remaining gaps receive a bounded native Monte Carlo
attempt (at most six candidates per context, 5,000 steps each). The complete
native construction process has a 180-second deadline and can be cancelled.

Short gaps can fail the original rigid-stem search even when their endpoints
are close enough: stem orientation and geometry matter as well as distance,
and the fragment database does not cover every fixed pair of stems. Temporary
context extension follows ProMod3's supported search mechanism. It does not
itself export changes to observed coordinates.
Context fragments still pass ProMod3's closure and ring checks. DynaMol
prefers candidates that displace the observed N/CA/C backbone atoms least, using the
native backbone score to break ties. This geometric ranking is designed
for restoring the original anchors; it is not a loop-accuracy score.
A research-only option also considers carbonyl oxygen because agreement on
N/CA/C alone can hide an incompatible oxygen orientation. Its first native
comparison did not resolve the backbone-reference failures, so the application
does not enable that ranking. The candidate budget is now enforced inside each
native collection, retaining at most 40 candidates per context search. Complete
construction/refinement changes still need the new release audit.

Only the requested missing atoms are transferred to DynaMol; changes that
ProMod3 makes to its temporary context are discarded. Its parent-residue
scoring context can omit ligands and modification atoms. The original
modified residues and retained molecules are preserved in the prepared system.

OpenMM then refines loop atoms and regenerated hydrogens against the complete
prepared force field, with the remaining heavy atoms fixed. This separate
construction step is limited to 1,000 minimizer iterations split between
collision relief and a geometry stage. If geometry still fails halfway
through, the second stage strengthens temporary peptide-planarity terms
fivefold. The prepared force field's equilibrium bond lengths, angles and
force constants stay unchanged. Full geometry checks still decide acceptance.
Temporary peptide
torsion restraints initialize new non-proline links as trans and retain the
candidate cis/trans basin for proline-like links. Those restraints are **not
exported to either MD engine**. Their energy is labeled separately; it is not
an unbiased force-field energy or a measure of model accuracy.

If the fixed-atom result still fails geometry, DynaMol can make one additional
local attempt with the immediate standard protein residue on each side of a
failed loop allowed to relax. These flanks receive strong positional
restraints and temporary chirality barriers. Their heavy atoms must remain
within **1 Å** of the coordinates entering refinement (after the separately
recorded sidechain sampling). A temporary radial penalty begins at 0.75 Å to
guide the minimizer away from this cap; the independent 1 Å acceptance limit
still applies. Ligands, ions, modified residues, observed
metal-contact residues and non-peptide crosslinked flanks stay fixed.
The fallback has its own 1,000-iteration limit, for **2,000 total at most**.
The original fixed-atom attempt is retained, and moved flank atoms and their
displacements are listed in the preparation diagnostics and summary. Flanks
that move become part of the provisional local model. Their geometry and
connections to the fixed protein must pass the same final checks; positions
are never clipped to make a failed displacement check appear to pass.

Every modeled loop receives connectivity, bond-length, backbone-angle,
peptide-planarity, chirality and retained-environment collision checks. Failed
candidates are retained in job diagnostics and are not accepted as prepared
datasets. Minimize the final simulation system before dynamics; this remains
required after loop construction.

The original structure remains available. Job and prepared-dataset files retain
the source identities, fragment-runtime versions and hashes, candidate atoms,
search stages and timings, discarded context displacements, refinement
settings, checks and warnings. The preparation monitor displays the active
search stage. Exhausted searches report the unresolved gaps without claiming
that they are impossible to model. The high-level native search APIs do not
expose the requested preparation seed; the Monte Carlo fallback uses native
default seed 0, recorded when it runs. The requested seed still controls other
preparation operations; cross-platform bitwise reproducibility is not claimed.

## Installation and evidence

The [15-structure preparation audit](audit/loop-fallback/REPORT.md) uses real
experimental missing regions, with its input set frozen before the final run.
Ten structures completed preparation, rebuilding 16 gaps and 90 residues;
five remained blocked by unsupported covalent chemistry, a ligand reference
conflict, or the 12-residue size limit. Both 12-residue examples and all four
8K5R gaps passed. Separate complete 8K5R repeats exercised the restrained-flank
fallback as well as the fixed-heavy-atom path. The report includes original
failures, final source hashes and independent saved-coordinate checks. No MD
was rerun for this preparation audit, and this small panel is not an estimate
of success across the whole PDB.

The rebuilt Mac bundle includes a separate ProMod3 runtime; no additional
account, server, GPU or installation is needed. Source-checkout users run
`bash scripts/install_loop_tools.sh`, or set `DYNAMOL_PROMOD3` to the matching
private runtime. It pins ProMod3 3.6.0, OpenStructure 2.11.1 and OpenMM 8.5.1.
The latter is private to loop modeling and does not replace DynaMol's MD engine.

The development checks rebuilt 1UA2's 12-residue gap while preserving its
observed protein, TPO and ATP coordinates during the fixed-context loop test.
Full preparation with ATP/TPO retained also passed. A separate 12-residue
deletion in ubiquitin passed the same geometry and chirality checks. These
are software feasibility checks, not a held-out loop-accuracy benchmark or
evidence of MD convergence; the fragment database may contain related structures.

Reference: Studer et al., [ProMod3—A versatile homology modelling toolbox](https://doi.org/10.1371/journal.pcbi.1008667),
*PLOS Computational Biology* (2021). ProMod3 is Apache-2.0 software; its
dependencies retain their respective licenses and bundled notices.
