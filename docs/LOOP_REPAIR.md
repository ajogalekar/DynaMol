# Missing-loop repair

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
fragment database and builds sidechains. The pipeline uses one bounded
fragment search, without Monte Carlo or extension into observed residues.
Only the requested missing atoms are transferred to DynaMol; changes that
ProMod3 makes to its temporary context are discarded. Its parent-residue
scoring context can omit ligands and modification atoms. The original
modified residues and retained molecules are preserved in the prepared system.

OpenMM then refines loop atoms and regenerated hydrogens against the complete
prepared force field, with the remaining heavy atoms fixed. This separate
construction step is limited to 1,000 minimizer iterations. Temporary peptide
torsion restraints initialize new non-proline links as trans and retain the
candidate cis/trans basin for proline-like links. Those restraints are **not
exported to either MD engine**. Their energy is labeled separately; it is not
an unbiased force-field energy or a measure of model accuracy.

Every modeled loop receives connectivity, bond-length, backbone-angle,
peptide-planarity, chirality and retained-environment collision checks. Failed
candidates are retained in job diagnostics and are not accepted as prepared
datasets. Minimize the final simulation system before dynamics; this remains
required after loop construction.

The original structure remains available. Job and prepared-dataset files retain
the source identities, fragment-runtime versions and hashes, candidate atoms,
refinement settings, checks and warnings. The fragment path does not apply a
random seed; the requested seed still controls other preparation operations.

## Installation and evidence

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
