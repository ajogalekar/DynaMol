# Covalent MD: literature-guided reassessment

Reviewed 2026-09-13 (America/Los_Angeles). This changes the research plan, not the application's accepted chemistry. No parameter model is approved by this document.

## Decision

Stop making a successful whole-adduct ab initio geometry optimization the prerequisite for every useful covalent benchmark. Use established preparation methods to build explicit reacted-residue candidates, then spend QM effort on identified parameter uncertainties. Preserve the protein force-field interface, source chemistry and complete-complex validation.

The existing 94-atom sotorasib optimization failed in its process monitor after 23 completed evaluations; it did not converge. Its files are retained. Neither restarting that calculation nor repairing more calculation infrastructure is the immediate priority.

## Primary evidence and its limits

| Source | What was actually demonstrated | Lesson for DynaMol |
| --- | --- | --- |
| Kong, Park and Im, 2024, [Covalent Ligand Modeling and Simulation](https://doi.org/10.1016/j.jmb.2024.168554) | Capped, combined amino-acid–ligand parameterization with CGenFF; 1,644 successful preparations. The exact 6OIM system was among the MD cases. Preparation success did not require low parameter penalties. | Borrow the combined-residue architecture and broad chemistry coverage. Preparation counts do not establish force-field accuracy. |
| 2023, [Covalent modifiers of the androgen receptor](https://doi.org/10.3390/biology12111442) | Modified cysteine building blocks using Antechamber AM1-BCC, PREPGEN and parameter checks; explicit inspection of backbone chemistry. The study used short peptides and assumed one linkage stereoisomer. | A fast, published route exists. Neither its stereochemical assumption nor its narrow validation should become a general app default. |
| 2019, [Sos1–Grb2 covalent conjugate](https://doi.org/10.1038/s41598-019-56078-7) | Selected small-molecule models, DFT geometry optimization, HF ESP/RESP charges, retained protein charge regions and GAFF2 supplements. Some related residues reused the calculated model data. | Custom QM can focus on a chemically justified model rather than every complete conjugate. Its separate artificial reaction-transition procedure is outside our fixed-adduct scope. |
| 2018, [Comparison of charge derivation methods](https://doi.org/10.1021/acsomega.8b00438) | Compared RESP and AM1-BCC on representative ordinary amino acids. Results often overlapped within uncertainty; unconstrained AM1-BCC did not preserve the specified backbone charges. Authors explicitly limited extrapolation to harder chemistry. | Full optimization and RESP are not universal proofs of better motion. Keep charge-interface validation; do not declare arbitrary charge redistribution equivalent to constrained fitting. |
| Pantsar, 2020, [KRAS–AMG 510 dynamics](https://doi.org/10.1038/s41598-020-68950-y) | OPLS3e/Desmond simulations with multiple independent trajectories; mobile protein switches despite stable ligand binding. The starting model was adapted from 6P8Y. | Protein motion away from the crystal pose is not automatically a failure. This is an external behavior reference, not an Amber parameter source or an unchanged-6OIM reproduction. |
| Li et al., 2022, [AMG 510 and KRAS dynamics](https://doi.org/10.1016/j.csbj.2022.02.018) | Describes 6OIM-derived ff14SB/GAFF/RESP simulations, with sequence changes and modeled gaps. The inspected methods do not fully specify every junction parameter and charge-boundary choice. | Relevant precedent, but insufficient detail to reproduce or approve our junction directly. Preserve our original panel rather than adopting their sequence alterations. |

These were read from publisher pages and PMC full-text XML, not inferred from abstracts alone. PubMed/PMC metadata identified the articles; metadata or failed fetches were not used as methodological evidence. No raw paper archive is reproduced here.

The [official CGenFF guidance](https://mackerell.umaryland.edu/cgenff_faq.php) specifically recommends smaller model compounds for difficult parameterization: distant intramolecular interactions can contaminate torsion fitting, and larger calculations can be less transferable as well as expensive. It also warns that CHARMM charges and other force-field families are not interchangeable.

The [official Amber modified-residue tutorial](https://ambermd.org/tutorials/basic/tutorial5/index.php) supplies a practical Antechamber/PREPGEN/LEaP workflow, including boundary atoms and parameters. It calls inferred parameters a starting point for validation. Its example charge redistribution is a method choice, not permission to silently change our current constrained-charge model.

## Revised implementation sequence

1. **Recover useful existing work first.** Inspect and reuse the preserved native combined-product AM1-BCC/GAFF2 calculation and actual native term inventory. The earlier local run completed in 176.5 seconds according to `../covalent-method-review.md`. That was an isolated candidate and conversion check, not a finished ff14SB residue. Re-running it or calling it validated would add no evidence.
2. **Make the residue boundary the next concrete milestone.** Produce one explicit, reviewable charge/interface protocol and exercise it in the existing neighboring-peptide assembly. AM1-BCC is a candidate baseline, not an automatic replacement for the unfinished constrained RESP model. Report the retained adduct charge, changes to backbone charges and local electrostatic effects. Do not silently mix charge protocols to force an integer total.
3. **Prioritize the uncertain terms.** Separate reaction-junction parameters from remote scaffold torsions. A high penalty identifies an analogy to examine; it is not a numerical error estimate. For needed QM refinement, use smaller capped models that retain the relevant electronic environment, then check transfer to the complete adduct. Fragment transfer is a hypothesis until checked.
4. **Keep geometry and charge jobs separable.** If constrained RESP remains necessary, choose and validate an affordable geometry-generation protocol before the expensive single-point ESP work. A cheaper geometry followed by HF ESP is a proposed comparison, not a result established here. Do not represent the failed optimizer checkpoint as an accepted geometry.
5. **Reach complete complexes sooner.** Retain 6OIM GDP, Mg and coordinating waters and qualify their model independently. Minimize and equilibrate the actual complex; inspect the junction, torsion distributions, stereochemistry, pocket contacts, temperature and solvent density. Use independent short trajectories for a numerical stability screen, with selected longer comparisons for motion. A 100-ps run cannot establish equilibrium behavior or parameter accuracy.
6. **Scale by chemistry, not by PDB ID.** Reusable templates must match reacted graph, stereochemistry, protonation, linkage and force-field version. Group the frozen 15 cases by these features, while retaining every original case and exception. One cysteine-adduct success cannot approve boronates, other residue chemistries or metals.

Candidate generation, numerical engine consistency, short-run stability and physical-model validation remain distinct results. Existing OpenMM/GROMACS transport checks are useful implementation evidence; they cannot select the more accurate force field.

## What is and is not reusable now

We have a useful published reference for the exact sotorasib case, but have not identified and verified a redistributable ready-made parameter bundle for it. Do not splice CGenFF terms into the Amber model. Treat the CHARMM route as a coherent external reference, and the installed AmberTools route as our local candidate path. New libraries or services are not required to carry out the next comparison.

The literature changes the order and scope of the calculations. It does not establish that DynaMol already handles all 15 covalent complexes, solve the separate metal-site work, or authorize an advanced-chemistry release.
