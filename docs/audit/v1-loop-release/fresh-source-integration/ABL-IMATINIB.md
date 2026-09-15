# Abl–imatinib: short-gap diagnosis and first complete static pass

Updated 14 September 2026 (local time). The selected 2HYY chain C has missing
ASP276 and THR389–GLY390–ASP391. The retained complex includes imatinib.

## Result

A targeted private comparison passes all eight complete static checks, including
the original-source and common-reference 1 Å observed-heavy displacement caps.
Maximum loop-flank observed-heavy movement is **0.767282 Å**. All eight modeled/boundary
Ramachandran definitions are scored without outliers; two are allowed and six
favored. The saved structure, original physical parameters, retained ligand,
outside torsion support, stereochemistry and geometry pass their existing checks.
The 1 Å limit concerns loop flanks. Earlier common-preparation sidechain repair
moved ARG239 NH2 by 7.432 Å from the deposition; this loop step preserves that
already prepared atom exactly. It is not a 1 Å cap on all protein preparation.

This is the first complete 2HYY static pass, **not yet an automatic app workflow,
fresh full-worker publication, MD stability result or evidence of the native loop
conformation**. The earlier failed comparisons remain failed. The first automated
panel still has five metal-free passes and four unresolved cases; this is a
separate targeted follow-up awaiting general integration and independent audit.

Evidence directory:
`/Users/ashujo/.cache/dynamol-research/v1-loop-release/abl-seed-preservation-v1`

Accepted bundle SHA256:
`061724fd2c121f4677269a9feeef57119ae45b8b6cc08e7416ea53029dab4471`

Frozen fresh snapshot SHA256:
`67a13780c6460d688af2a4a21b63053aa72419261b20660a73fe38c33329dddb`

## What was wrong

1. Source identity was correct: all 2,120 source coordinates reproduce the
   deposited structure, with exact missing sequence identities. The deposited
   LYS274 outlier is pre-existing and its outside torsion support remains fixed.
   There is no established gross endpoint-distance impossibility.
2. Individually restoring partial observed anchors after inserting a fragment can
   distort the join. Proper rigid placement of each complete fragment and its
   immediate flanks before exact anchor restoration improves some seeds, but
   did not by itself produce an accepted result.
3. The nearest one-context fragment for the one-residue gap has a cis peptide
   before ASP276 (about 5°). That conflicts with the existing construction policy
   for new non-Pro links. The default minimizer targets trans; the conversion
   damages the backbone conformation. This does **not** prove the missing native
   peptide cannot be cis: its state is unobserved. No same-entry chain supplies a
   complete observed GLU275–ASP276–THR277 segment.
4. A read-only audit of the already bounded fragment pools found alternatives
   with acceptable starting torsions after proper placement. The chosen example
   uses pool rank 3 for each gap, with exact native backbone bytes checked before
   sidechain reconstruction. Its ordinary refinement still creates two outliers.
5. A paired private refinement protects the initially allowed/favored phi/psi
   values with temporary circular torsion restraints (200 kJ/mol). It passes the
   unchanged final checks. These construction forces are not included in the
   saved native physical System or MD parameters. It uses no new QM or MD.

The temporary restraints bias the starting model toward the selected fragment.
Validation of that model is not independent evidence that the fragment is the
native conformation. Using such restraints to preserve a good starting backbone
is consistent with the distinction made in the
[PHENIX refinement guidance](https://www.phenix-online.org/documentation/tutorials/lowres_restraints.html);
they should not be used merely to conceal a poor reconstruction.

## General correction already implemented

The opt-in fragment worker now rejects undefined, ambiguous or cis-basin **new
non-Pro** peptide candidates after CCD and before anchor ranking. It uses exact
sequence mapping, preserves observed–observed context states and exempts X–Pro
candidate states. Rejected candidates consume the existing 40-candidate examined
budget. Unavailable outcomes retain their rejection evidence. The legacy app and
MC policies remain unchanged. The focused software suite passes **71 tests**;
a native replay confirms rejection of both incompatible short-gap candidates.

## Next gate

Replace the private rank choices with a general, bounded selection rule based on
the placed candidate's complete joins and source movement. Preserve initially
qualified conformations during refinement without relaxing acceptance limits.
Then repeat from fresh source through the real preparation worker and test other
eligible cases. Do not enable the app or claim universal recovery from this one
successful comparison. Metal-containing cases remain outside v1.
