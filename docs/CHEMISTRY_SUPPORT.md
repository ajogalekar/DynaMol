# DynaMol chemistry support and release checks

DynaMol distinguishes **viewing**, **chemical identity**, **preparation** and
**numerical simulation checks**. Loading a molecule does not establish suitable
MD parameters. A successful short run establishes software continuity, not the
accuracy of a protonation state, force field or biological conclusion.

| Molecular type | Current preparation behavior | Evidence and limits |
| --- | --- | --- |
| Standard amino acids | ff14SB hydrogen/template assignment, missing-heavy-atom repair, recorded pH-dependent states | All 20 amino-acid types represented in capped native fixtures; 26 residue/state cases plus free termini and a cysteine disulfide checked |
| HIS protonation/tautomers | HID/HIE/HIP templates supported; OpenMM chooses protonation/neutral tautomer at preparation | Actual bonded-H inventory is recorded. Imported H is removed when Prep is requested; it is not a constant-pH or pKa calculation |
| ASP/GLU/LYS/CYS charge variants | ASH/ASP, GLH/GLU, LYN/LYS and CYM/CYS represented | LYN/CYM input naming is normalized to the same LYS/CYS heavy graph for explicit pH reassignment, with the original alias recorded. CYM thiolate is distinguished from CYX disulfide by actual S–S connectivity |
| Termini and caps | Ordinary charged N/C termini and built-in ACE/NME caps | Checked as complete fragments. This is not proof of every residue/state/terminal combination; missing or unsupported terminal templates remain errors |
| SEP, TPO, PTR, HYP | Exact curated covalent amino-acid templates, explicit TIP3P/OpenMM | Existing four-residue, three-pose native Amber parity suite; phosphate residues use fixed −2 states. HYP internal/CT states and supported stereochemistry only. [Registry details](MODIFIED_RESIDUES.md) |
| Other modified/unnatural amino acids | Recognized as polymer and retained; specific missing-template blockers | Negative identity/preservation checks include MSE, ALY, CSO, FME, DLY and unknown ZZZ. No conversion to parent residues or disconnected GAFF ligands. MSE/ALY's separate published `ff14SB_modAA` integration is still future work |
| Noncovalent organic ligands | Retained graph/CCD/explicit SMILES → recorded state → native GAFF2/AM1-BCC | Neutral alcohol, anionic carboxylate/phosphate, cationic amine, halogenated aromatic and chiral amide representatives. Every generated XML must reproduce native Amber energies and forces |
| ATP, ADP, NAD, FAD | Supported as complete, noncovalent organic cofactor graphs through the ligand route | At nominal pH 7, these fixtures selected net charges −4, −3, −1 and −2 respectively. These are whole-molecule net charges; NAD's pyridinium center is +1. They do not establish the dominant bound state or redox behavior |
| Na+, Cl−, K+, Mg2+, Ca2+ | Explicit TIP3P-compatible nonbonded ion models | Exact element/name/atom-count and template charges checked; SOD/CLA/POT/CAL aliases also checked. Five ions survived actual capped-protein preparation and solvation |
| Heme/organometallics, other metal states | Require specialized models; no guessed oxidation state or metal parameters | Iron-containing HEM explicitly rejected. Retaining nonbonded Mg/Ca and observed donor/water coordinates does not validate coordination energetics |
| Covalent inhibitors/cofactors and nonstandard crosslinks | Explicit blockers unless a matching supported covalent template exists | Original covalent FAD linkage and nonstandard protein crosslink checks; no bond is deleted to make a ligand appear independent |
| Radicals, incomplete/ambiguous ligand graphs, other elements | Explicit actionable rejection | Boron-containing ligand, radical methyl and ATP missing a heavy atom tested; no closed-shell or missing-atom replacement is invented |

Preparation and the readiness pane share `validate_preparation`, so known
structural gaps, missing templates, disabled required repairs and ligand chemistry
errors use the same eligibility rules. Inspection can retrieve a missing CCD
reference into its cache, but it does not create a job or alter a dataset.

Retained ions now carry their declared states in preparation metadata and require
explicit water even when no organic ligand is present. The native five-ion
preparation/solvation fixture contains 3,785 atoms. Supported prepared complexes,
PTMs and ion states are carried into OpenMM; GROMACS prepared-state transfer is
still not implemented.

## Findings corrected in this pass

The halogenated aromatic fixture exposed a legitimate native Amber improper
quartet whose central atom is second. Reversing all four atom indices restores
the central-third convention without changing the torsion, energy or forces.
The converter now handles this case and still requires per-ligand native parity.

The NAD fixture exposed an omitted atom-map number in Dimorphite's phosphate
normalization. DynaMol now recovers missing maps only when the complete
heavy-element/connectivity graph and every surviving map permit exactly one
solution. Conflicting maps, ambiguous mappings, changed connectivity and bounded
search exhaustion remain failures. Original coordinates and stereotags are
preserved; the recovery is disclosed in the inspection warning.

The matrix also found misleading CYM/CYX state reporting and an omitted
explicit-solvent requirement for ion-only protein complexes. Both were corrected,
and unnamed or element-mismatched ion states now receive a specific message.

## Reproduce and interpret the evidence

```bash
DYNAMOL_DATA_DIR="$PWD/data/integration-checks/release-chemistry" \
  .venv/bin/python scripts/check_chemistry_matrix.py
DYNAMOL_DATA_DIR="$PWD/data/integration-checks/release-chemistry" \
  .venv/bin/python scripts/check_chemistry_boundaries.py
.venv/bin/python -m pytest backend/tests/test_chemistry_support.py \
  backend/tests/test_residue_identity.py backend/tests/test_modified_residues.py
```

The [matrix](audit/release-readiness/chemistry-matrix.json),
[expected boundaries and actual ion pipeline](audit/release-readiness/chemistry-boundaries.json),
and [focused regression output](audit/release-readiness/chemistry-tests.txt)
record the checked outcomes. An initial matrix record is retained separately;
it includes the bugs above and test-fixture issues corrected during the audit.
In particular, the original heavy-coordinate assertion compared PDB-rounded
preparation input with the higher-precision display storage. The corrected test
compares against the exact input consumed by preparation; it does not loosen the
heavy-atom preservation tolerance. Two protonation aliases and an incorrectly
constructed native terminal fixture were corrected explicitly.

Native ligand logs, parameters, input coordinates and per-atom charges are retained
under the isolated data directory listed in each record. CCD ideal coordinates
are reproducible test geometry, not deposited bound poses. The ten ligand/cofactor
fixtures use five 0.2-fs nonperiodic OpenMM steps after bounded minimization;
these are numerical checks, not aqueous dynamics or equilibration. Existing
complex and modified-residue end-to-end checks provide separate explicit-solvent
continuity evidence. No statistical accuracy or convergence interval is inferred
from these deterministic tests.

Primary implementation references: [OpenMM hydrogen and model building](https://docs.openmm.org/latest/userguide/application/03_model_building_editing.html),
[OpenMM force-field resources](https://github.com/openmm/openmmforcefields),
[RCSB CCD reference downloads](https://www.rcsb.org/downloads/ligands), and
[Dimorphite-DL](https://github.com/durrantlab/dimorphite_dl).
The [evidence retrieval](audit/release-readiness/blog-search.json) searched the
complete cached Practical Cheminformatics corpus. The caution that a chosen
molecular representation need not be the populated experimental state is relevant
here; applying it to these MD preparation choices is this language model's
inference, not expert certification or an endorsement by Pat Walters.
See [The Trouble With Tautomers](https://patwalters.github.io/The-Trouble-With-Tautomers/).

The chemistry [automated postflight](audit/release-readiness/chemistry-postflight.json)
reports **7 passes and 1 failure**, with decision `block`, because scientifically
meaningful uncertainty intervals are explicitly absent. That gate is not claimed
to have passed: these fixtures cannot support accuracy, equilibrium or predictive
model promotion. Passing the recorded functional/numerical tests supports only
the narrower software behavior described above. The review is a model-generated
evidence synthesis with an automated checklist, not independent expert validation.
