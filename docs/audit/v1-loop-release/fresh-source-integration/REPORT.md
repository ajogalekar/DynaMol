# Fresh-source loop integration

## Fresh protected-flank replay and representative panel — 14 September

The corrected fresh 1UA2 preparation took 70.5 seconds and retained ATP/TPO,
all 12 requested residues and the fixed hydrogen/parameter state. All 16 raw
heavy atoms of the two immediate observed flanks stayed exactly in place
during preliminary sidechain sampling. The complete search then evaluated
three candidates in 22.0 seconds: fragment and the first torsion candidate
failed reference checks, while torsion sampling seed 2027 passed all eight
checks, including the original-source 1 Å flank limit and saved-file checks.

Evidence: `~/.cache/dynamol-research/v1-loop-release/fresh-protected-loop-search-v1/`.
This is a fresh common-preparation plus complete static validation result;
returning the accepted coordinates through the real worker's final publication
path and enabling the app fallback remain outstanding. No MD run or proof of
the native loop conformation is claimed.

The frozen 15-entry panel is now running two cases at a time. It has ten
repair candidates and five expected unsupported-input checks, retaining each
selected complex's ligands/ions. 1M17 and 8K5R are first; the fresh 1UA2 result
is included with its seed explicitly recorded. Results are at
`~/.cache/dynamol-research/v1-loop-release/frozen-panel-fresh-validation-v1/results.json`.
The target is a focused v1 release tomorrow morning, contingent on the
remaining preparation, app integration and exact installation checks.

The earlier failed captures and the superseded pre-fix account below are
retained as historical evidence.


The backend now has a supervised native proposal worker, strict context
admission, immutable complete preparation snapshots, and a complete static loop
validator with a bound output bundle. Automatic fallback is still private;
this is not an enabled app workflow or release qualification.

The focused backend regression run passed 228 tests. One additional existing
solvation test cannot read its offloaded local demo PDB and has not passed.
The private snapshot wrapper has 11 passing integrity tests with refinement
stubbed out. Independent review confirms exact numerical-core and native
output parity of the extracted backbone scorer against the previously audited
implementation. These software checks do not establish structural accuracy.

One fresh native proposal for the 12-residue 1UA2 gap completed in 1.994 seconds,
with 91 modeled heavy atoms and 16 observed-context atoms. Schema, exact atom
coverage and source checks passed. This used an earlier observed-source file;
its hash differs from the new preparation context, so that proposal must not
be relabeled or reused as a fresh-source-bound result.

## Actual preparation from the original structure

The private capture began from unprepared dataset `3a661cfc670048c5`, containing
2,322 atoms and the original sequence evidence. No prepared model archive was
an input. The first attempt failed because it selected the development
ProMod3 executable, which could not launch. That failed attempt remains intact.
One corrected retry explicitly selected the installed, tested ProMod3 runtime.

The retry completed common preparation in 51.61 seconds: 4,825 atoms, all 12
requested modeled residues, retained ATP/TPO, and 2,411 assigned hydrogens.
ATP was freshly parameterized with the existing AM1-BCC workflow, including a
native SQM calculation. No MD ran. Peak sampled aggregate memory was 206.6 MB;
all owned child processes exited. The capture contains 75 source files and
12 parameter files, including exact atom/H-parent mappings, state decisions,
recursive base force-field inputs, modified-residue sidecars and reference
templates. Its SHA256 is
`b6e228bf42119ad738e86a814cf33611dce23214add4177517f7f825e575b432`.

The capture deliberately stopped before loop refinement and publication.
It does not claim a completed preparation or accepted loop model.

## Boundary staging issue found

All original protein backbone coordinates remain fixed, and retained ATP
coordinates are unchanged within floating-point precision. However, ordinary
preliminary chi sampling changes observed sidechains, including ILE A43 next
to the missing loop. Its CD1 moves 3.720 Å before the complete loop validator
runs. The existing common-reference cap and the separate original-source cap
are both 1 Å. No final coordinate can satisfy both for an atom whose two
reference positions differ by more than 2 Å. See
[the exact movement evidence](fresh-flank-reference-conflict.json).

Complete refinement was therefore not launched against that incompatible
reference. The narrow correction is to protect immediate observed peptide
neighbors of modeled residues during preliminary chi sampling, preserving
existing metal-site protection. The complete candidate refinement will handle
their limited movement later. Other sidechain optimization remains available;
no acceptance threshold is relaxed.

That change is implemented in the preparation worker and has passed an
independent static review. Six focused tests pass, including actual chi movement
remaining blocked only for the intended flank/metal residues, while other and
modeled sidechains can still be optimized. Its fresh native preparation replay
is pending.

The exact pinned physical System also passes the
[static parameter reload check](fresh-static-parameter-reload-v2.json).
The initial comparison rejected 904 harmonic bonds solely because their
particle endpoints were reversed after canonical topology reconstruction.
The comparator now recognizes those undirected terms while preserving every
parameter value and term multiplicity. All other System fields and terms
match exactly. The reloaded System's SHA256 is unchanged between the failed
and corrected comparisons; no physical parameter was changed. Both initial
and corrected evidence are retained. The comparator and wrapper checks pass
20 focused tests; no Context or dynamics was needed for this diagnosis.

Next verify that correction, repeat fresh capture with preserved boundaries,
generate proposals from the exact new source context, and run the full loop
validator and saved-bundle handoff. Then exercise the frozen protein panel,
including every 8K5R gap and the second 12-residue case. App enablement,
packaging and the queued inaugural-video update remain behind those gates.

Full capture evidence:
`~/.cache/dynamol-research/v1-loop-release/fresh-preparation-snapshot-v1-attempt-02/`.
Independent capture review:
`~/.cache/dynamol-research/v1-loop-release/fresh-preparation-capture-review-v1/review.json`.
