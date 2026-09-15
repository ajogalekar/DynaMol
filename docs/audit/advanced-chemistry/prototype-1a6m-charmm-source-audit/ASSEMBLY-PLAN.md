# 1A6M: complete CHARMM assembly plan and static source audit

**The published His–oxyheme model is mechanically buildable as a separate CHARMM model, but the complete frozen 1A6M case is still blocked by sulfate coverage.** Both sulfates and all crystallographic waters remain in the plan. No 1A6M System, simulation, QM calculation, or application readiness change was made.

The reference is the [2015 oxy-myoglobin study](https://doi.org/10.1371/journal.pone.0128496) and its [deposited topology, parameter file, construction script, PSF and PDB](https://doi.org/10.5061/dryad.5g1k6). This is a legacy CHARMM22 protein/CHARMM27 lipid parameter lane with additive CHARMM TIP3P water, unscaled electrostatic 1–4 interactions, and the published atom-specific Lennard–Jones rules. It is not an Amber ff14SB add-on. The source files remain unchanged in `../non-zinc-models/dryad-oxy-myoglobin/`.

## What the new static audit establishes

- All **27,386** atoms in the deposited horse PSF match the supplied baseline residue atom types and partial charges after applying its declared GLYP/CTER terminal patches. There are no unexplained template mismatches. This includes 8,284 TIP3P waters, 25 sodium ions and 23 chloride ions in the published PSF.
- The supplied horse PDB has **2,486** solute atoms. Every atom maps to the PSF's solute prefix after precisely the two aliases declared in the original construction script: HIS→HSD and ILE CD1→CD. The 196 alias records are retained. The PDB does **not** supply the full PSF's solvent/ion coordinates or periodic box; it was never assigned to that complete System by index.
- Two distinct mechanical definitions were built using OpenMM **8.6.0.dev-c6173db**: baseline and the separately deposited updated charges. Both contain 27,386 particles, 38,552 nonbonded exceptions and zero constraints. Particle charges and masses reproduce their respective PSFs exactly. These are imported System definitions; no Context, energy evaluation or motion was run.
- Updating the topology's declared HEM/OXY charges changes exactly **27 particle charges and 189 exception charge products**. All other force definitions, masses, Lennard–Jones terms, exception pairs and exclusions remain identical. Every nonzero exception charge product follows the declared unscaled CHARMM 1–4 rule. This is an internal charge-variant check, not independent native-engine energy/force validation.
- All six duplicate-improper reader warnings were traced to the original parameter lines. Each duplicate has the **same numeric force constant and equilibrium angle** as the prior definition. Last-definition precedence therefore changes no numeric term for these six keys. The warnings and original lines are retained.

| Explicit version | HEM total partial charge | OXY total partial charge | Fe partial charge | O1 / O2 partial charges | Matches deposited PSF? |
|---|---:|---:|---:|---:|---|
| `dryad-2015-baseline` | −2.00 | 0.00 | +0.24 | +0.021 / −0.021 | Yes |
| `dryad-2015-updated-charges` | −1.50 | −0.50 | +1.42 | −0.18 / −0.32 | No; separate derived PSF |

Partial charges do not specify formal Fe oxidation state or quantum spin. The deposited construction script loads the baseline topology, so the deposited trajectory cannot be claimed as validation of the updated-charge variant. The published zero O2–O1–Fe angular force constant and zero O2–O1–Fe–porphyrin torsion remain exactly unchanged. Other bonded and nonbonded terms still contribute to the whole-site angular potential; that potential needs independent validation.

## Exact frozen 1A6M input and mapping

The [1A6M source structure](https://www.rcsb.org/structure/1A6M) is SHA-256 `25ae24d2b518a7cb7ef24aeb1eb3fc5c35707ce4e696956d11b8f9fdf9fa7218`. Assembly 1 is an identity-operation monomer containing all label chains A–F. No component or symmetry copy was silently removed.

The source contains 1,584 heavy-atom records, including 139 alternate-location pairs. A documented whole-residue highest-occupancy selection, with lexical tie-breaking, retains **1,445 distinct heavy atoms**; all 139 alternate records remain in the original CIF and alternate-selection audit. This covers 151 observed protein residues, 43 heme heavy atoms, 2 oxygen-ligand atoms, 10 sulfate atoms and 186 distinct crystal-water oxygens. No internal sequence gap is present.

| Source identity | Proposed CHARMM identity | Required treatment |
|---|---|---|
| Protein label A, author A1–A151 | Segment `PROA`, original residue numbers | Rebuild from the 1A6M sequence, preserving all selected observed heavy atoms. ILE CD1→CD is an explicit alias. |
| VAL A1 / TYR A151 | `PROA:1` NTER / `PROA:151` CTER | Do not copy horse GLYP: horse residue 1 is GLY. Terminal source O/OXT map explicitly to CHARMM OT1/OT2. |
| HIS A93 NE2 | `PROA:93:HSD:NE2` | Published proximal model is neutral ND1-protonated HSD; NE2 binds Fe. Other His protonation states still need assessment. |
| HEM author A154, label D | `HEME:154:HEM` | Preserve all 43 heavy-atom names and coordinates; add only declared template hydrogens. |
| OXY author A157, label E | `OXY:157:OXY` | Preserve source O1/O2 identities. Published model bonds only Fe–O1; the additional deposited Fe–O2 contact remains coordination evidence. |
| SO4 author A155, label B | `SUL1:155:SO4` | Retain S/O1/O2/O3/O4, source dianion state; compatible complete additive parameters unresolved. |
| SO4 author A156, label C | `SUL2:156:SO4` | Same obligation; no removal or conversion to sulfuric acid. |
| HOH label F, all source author IDs | Segment `WAT`, `TIP3`; O→OH2 | Retain all 186 selected oxygen positions and author identities; account explicitly for each added H. |

The exact source row, author and label identities, selected alternate, proposed CHARMM identity and coordinate are stored for every selected atom in `1a6m-identity-map.json`. Aliases never overwrite the source identity. The proposed all-HIS→HSD policy mirrors the original script and is **provisional** for nonproximal histidines; it is not a pKa calculation.

## Complete assembly sequence and remaining gates

1. Select **one named charge version** and its immutable topology hash. Use the supplied full `par_all27.inp`; preserve the CHARMM water and 1–4 conventions throughout. Do not load ff14SB, GAFF or a different protein force field around the published heme patch.
2. Resolve a complete additive all-atom sulfate topology and parameter source compatible with this Hamiltonian, including its charge state, bonded terms, Lennard–Jones interactions and redistribution terms. Both sulfate identities above are required before a runnable 1A6M assembly is attempted.
3. Construct protein segment PROA from the selected 1A6M coordinates with explicit NTER/CTER and histidine assignments. Build separate HEME, OXY, SUL1, SUL2 and water segments with no automatic terminal patches. Preserve source-heavy coordinates during hydrogen placement; record every addition and alias.
4. Follow the published patch ordering and term-generation policy: protein/heme auto angles and dihedrals, oxygen auto none; apply `PHEM PROA:93 HEME:154`, then `PLO2 OXY:157 HEME:154 PROA:93`. Check the specific patch-added terms and patch-deleted trans-heme angles against `published-site-terms.json`. Do not indiscriminately regenerate extra metal-site angles or torsions after patching.
5. Verify complete source retention and all prepared bonds, angles, proper/improper torsions, nonbonded exclusions, 1–4 parameters, charges and masses. Native CHARMM/NAMD versus imported-engine static energy/force parity remains to be performed for the complete assembled model. A successful OpenMM parser alone does not establish that parity.
6. Before describing the site model as validated, assess the **whole-site** oxygen bending and torsion energy behavior against independent reference evidence, along with Fe/His/porphyrin geometry and distal interactions. Keep the published zero terms unless a separately justified new model is declared. Only then consider a short dynamics workflow check; intact bonds alone are insufficient evidence of realistic motion.

## Sulfate coverage was checked in actual official files

The [official MacKerell CHARMM force-field page](https://mackerell.umaryland.edu/charmm_ff.shtml) supplied both a legacy `toppar_c31b1.tar.gz` archive and the current `toppar_c36_feb26.tgz`. Downloaded bytes and hashes are retained. We inspected all ordinary five-atom, net−2 residue definitions, as well as named sulfate entries: the only such sulfate template in either archive is **`toppar/toph19.inp` SO4**. This is the CHARMM19 model already excluded as an unvalidated substitute for the published CHARMM22/27 lane. No compatible additive sulfate model was found in this bounded source audit. This is not a claim that none exists elsewhere. A Drude sulfate model, sulfated detergent, or sulfate-ester fragment is a different model and is not silently substituted.

The Dryad dataset declares CC0. The additional official CHARMM archive downloads were inspected as research source material; this audit does not establish independent redistribution permission for those complete archives. No new files were integrated into DynaMol's release bundle.

`audit_assembly.py` reproduces the core source, alternate selection and two-model build checks. `check_variant_systems.py` reproduces the charge-only mechanical comparison. JSON evidence includes `published-template-mapping.json`, `published-coordinate-mapping.json`, `published-native-system-audit.json`, `variant-mechanics-comparison.json`, `improper-duplicate-audit.json`, `1a6m-assembly-plan.json`, `additional-retrieval.json` and `sulfate-archive-audit.json`. Original inputs and derived artifacts are separately hashed in `artifact-manifest.json`.
