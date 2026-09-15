# 1A6M sulfate source follow-up

**The complete 1A6M assembly remains unresolved. No compatible, complete additive sulfate source was verified, and no System, Context, QM calculation or dynamics was launched. Both source sulfates remain required.** This is a bounded source finding, not proof that an appropriate model does not exist elsewhere.

## What the broader source audit found

The [official CHARMM download page](https://mackerell.umaryland.edu/charmm_ff.shtml) identifies the February 2026 distribution and its separate Drude content. The maintainer repository was pinned at `ed3bbdaa4eb0a0152e577480bb0a9d22cd9743c3`. An independent composition/alias search covered 911 residue/patch blocks in the legacy archive and 8,406 blocks in ordinary current topology/stream files. It searches sulfur/oxygen atom inventories using MASS definitions, rather than relying on SO4 residue names. It also inspects nested Drude archives, covering another 476 blocks. This is a candidate search, not a general CHARMM topology parser.

| Candidate | Exact chemistry/model finding | Decision |
|---|---|---|
| Legacy `toph19.inp` SO4 | Five atoms, net −2; S charge 0 and four OC charges −0.5. `param19.inp` provides S–O and O–S–O terms, wildcard `S*` LJ parameters, and explicitly labels the sulfate angle force constant as a guess. Its global electrostatic 1–4 factor is 0.4, with legacy water NBFIX rules. | This is not the published oxyheme CHARMM22/27 Hamiltonian, whose factor is 1.0. No unsupported transplant. |
| Current CGenFF and lipid `MSO4` | Nine-atom methyl sulfate, net −1, with a covalent O–CH3 group. Sulfated detergent and carbohydrate definitions are also derivatives or patches. | Removing methyl atoms and changing charges would create a new model; it is not a residue alias for SO4²⁻. |
| Nested `drude_toppar_2023` SO4 | A real five-nucleus, net−2 sulfate topology with explicit ALPHA/THOLE terms. The corresponding SD1A/OD2C2B bond, angle and LJ definitions are present in the Drude main stream. | Complete model family is polarizable. Omitting its polarization particles or mixing it into the additive oxyheme model is not justified. |

The nested Drude finding clarifies the previous ordinary-file search: the current download **does contain sulfate**, inside a separate polarizable package. It does not close the additive-model requirement. Source lines, units/comments and hashes are retained in `nested-drude-sulfate.json`, `candidate-parameter-lines.json` and `legacy-parameter-context-v2.json`; no parameter values were converted or installed. The legacy context correction retains section alias `THETAS`, wildcard `S*`, and the original 1–4/water rules, which simple text-hit labels did not parse.

## New primary-literature lead, still missing complete files

[Kuhn et al., Nature Communications 2026](https://www.nature.com/articles/s41467-026-75749-4) explicitly used a Cannon sulfate model with CHARMM36m/CHARMM TIP3P and a CUFIX modification for sulfate–protein interactions. This is relevant additive-model application evidence. The article makes raw data available on request and links an analysis repository.

That [author repository](https://jugit.fz-juelich.de/computational-neurophysiology/slc26a11-ion-binding), pinned at `33c0063f78ddf67f822d699fcf03d3c0d25fd501`, has 40 tree entries: two directories, 36 PDB structures, README, and the analysis script. Its complete pinned tree has no topology/parameter files or sulfate-specific nonbonded override table. Structures alone cannot recover those parameters.

The [official CUFIX page](https://bionano.physics.illinois.edu/CUFIX) supports CHARMM36/27/22 use with its specified ion and CHARMM TIP3P conventions. Its downloaded stream contains general pair corrections but no sulfate residue definition or explicit mapping of the 2026 study's sulfate types to those corrections. Neither its name nor the 2016 amine–phosphate citation establishes the complete sulfate model. No numerical sulfate–protein correction was guessed, and the study's newer protein force field was not substituted for the published oxyheme model.

## Remaining concrete requirement

Obtain the actual complete additive sulfate topology, internal bonded/constraint policy, atom charges, LJ and pair-specific tables with clear source provenance, then verify compatibility with the retained oxyheme/water Hamiltonian. A separately developed sulfate model would require its own water/ion/protein interaction validation. The promising 2026 author data is a retrieval lead; it is not a ready parameter bundle. No author message was sent.

The article declares CC BY 4.0. The public force-field and code repositories inspected here do not independently establish an unrestricted redistribution license for their complete contents; the official GitHub metadata reports no license. No new source material entered a DynaMol release. `retrieval.json`, the two pinned repository-tree records, and `artifact-manifest.json` preserve this audit's provenance.
