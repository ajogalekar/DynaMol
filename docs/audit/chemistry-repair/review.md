# Chemistry repair review

The scoped implementation checks pass. Two concrete review findings were fixed: removed molecules and old hydrogens no longer obstruct repair of retained ligands, and display categories now use the same named ion registry as preparation, including the CAL alias. No chemistry software finding remains open in this review.

This is a language-model synthesis of source records, code and observed test results, not professional certification or an endorsement by Pat Walters. The PatAgent preflight passed all six gates. Postflight has **seven passes and one failed uncertainty gate**: the saved checks do not establish the experimental positions of missing atoms or validate metal-site coordination. The automatic scientific-claim gate remains blocked; no confidence intervals, native-pose accuracy or MD convergence claim is made.

## Evidence

| Check | Observed result |
|---|---|
| Ligand compatibility regressions | 61 passed, four existing expensive cases deselected; the new repair module's 11 cases also passed independently |
| Ion parameters and final-System validation | 73 passed, six unrelated cases deselected |
| Shared ion display and roundtrip checks | 30 passed independently, covering named ions, CAL precedence, wrong elements and multiatom compounds |
| Native ligand parameter conversion | Repaired F:3004 retained 13 observed heavy atoms and gained three; actual antechamber/parmchk2/tleap and Amber-to-OpenMM validation passed three poses, with maximum energy/force differences of 1.124×10⁻⁸ kJ/mol and 1.922×10⁻⁷ kJ/mol/nm |
| Complete isolated 6A93 preparation | Completed in 66.56 seconds with 6,142 atoms, five organic residues and one verified Zn²⁺ ion; no dynamics integration |
| Independent saved-coordinate replay | All 98 original ligand/ion heavy atoms had exactly zero displacement in the stored physical frames; both original dataset file-hash sets were unchanged |
| Browser interaction checks | Five Chrome cases passed without skips or failures: missing-loop controls, repair/remove intent, resize persistence/reset, canvas bounds and fullscreen header actions |

Counts overlap and should not be added as a unique-test total. Exact scope, source/evidence hashes and commands are in [the postflight manifest](review-manifest.postflight.json). Native files and parameter logs remain in the isolated integration directories. [Independent replay](review-independent-native.json), [browser report](review-browser-results.json), [resolved findings](review-findings.json) and [PatAgent result](review-postflight.json) preserve the evidence. Browser tests use real dry inspection/readiness and inspect the outgoing Prep payload; actual preparation is established separately by the recorded native worker result.

## Identity and coordinate handling

[RCSB identifies 1PE as neutral pentaethylene glycol](https://www.rcsb.org/ligand/1PE), with 16 heavy atoms. A freshly retrieved official CCD file matched the local cache byte for byte. In the original structure and selected chain, F:3004 lacked C12/C22/OH2, G:3005 was complete, and H:3006 lacked C15/C16/C25/C26/OH6/OH7. Final output has 16 heavy atoms in each copy and neutral assembled charge within native numerical tolerance. Only F and H contain modeled-atom records; the complete G copy has no modeled additions. Their explicit per-residue repair choices and low-confidence warnings survive preparation. Source sequence records identify no missing protein residues in this 6A93 input, and no protein segments were invented.

Repair requires explicit opt-in, a compatible full CCD graph and stereo reference, at least three noncollinear observed heavy atoms and one connected observed core. It is limited to eight missing heavy atoms and at most half the ligand; missing ring atoms, ambiguous identities and unsupported covalent chemistry remain actionable errors. Eight seeded candidate attempts use fixed observed atoms, bounded local optimization and a gross-clash screen. Explicit removal applies only to the chosen validated noncovalent residue. Readiness inspection does not generate conformers or silently accept a default repair.

The retained observed coordinates are checked after completion and assembly. Modeled atom names, positions, seed, candidate outcomes and parameter mappings are recorded separately. The generated-coordinate record differs from canonical stored physical coordinates by at most 5.4×10⁻⁶ Å in this test, consistent with finite-precision serialization. This does not alter the exact observed-heavy-atom preservation result.

## Scientific limits

[RCSB's ZN component](https://www.rcsb.org/ligand/ZN) is monatomic Zn²⁺. The installed Amber14 TIP3P file provides the verified +2e template, σ=0.22646645415127425 nm and ε=0.01381916624 kJ/mol. Preparation checks actual final-System parameters and retained ion identity rather than accepting a residue name alone. The [OpenMM documentation](https://docs.openmm.org/latest/userguide/application/02_running_sims.html#amber14) pairs Amber14 water files with compatible ions; the associated [Li et al. parameter study](https://pubmed.ncbi.nlm.nih.gov/23914143/) concerns water-dependent nonbonded models. These checks support using that installed starting model, not arbitrary oxidation-state assignment, directional coordination, polarization or catalytic-site accuracy.

[RDKit documents constrained embedding](https://www.rdkit.org/docs/source/rdkit.Chem.AllChem.html#rdkit.Chem.AllChem.ConstrainedEmbed) as approximate core placement; this implementation additionally fixes and tests observed coordinates explicitly. Successful local optimization and native force-field translation do not recover experimental evidence for unseen atoms. Candidate local energies and the 1.6 Å gross-overlap cutoff are limited construction checks, not affinity, equilibrium probability or pose-quality measures. No production MD was run for this review. The source and claim analysis is detailed in [reference notes](review-reference-notes.md).
