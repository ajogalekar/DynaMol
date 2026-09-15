# Final chemistry audit — September 12, 2026

**Chemistry preparation passed on the representative cases.** The 6A93 water preview was refused at the default 100,000-atom limit; this is retained as a resource boundary, not counted as a completed preview.

Fresh evidence: **44/44 native matrix cases**, **6/6 boundary checks**, and **3/3 real complex preparations**. Explicit water previews: **2 completed; 1 rejected by the atom cap**. No chemistry implementation or validation tolerance was changed to achieve these outcomes.

| Input | Fresh preparation coverage | Prepared atoms | 1 nm water preview |
| --- | --- | ---: | --- |
| 6DBK | Protein, PTR, charged G5D ligand; sidechain and hydrogen preparation; native ligand parameter parity | 4,713 | Completed: 50,033 atoms |
| 6A93, observed chain A | Charged 8NU, cholesterol CLR, three 1PE molecules including two explicit repairs, and Zn2+; native ligand parameter parity | 6,142 | Refused: default atom cap exceeded |
| 1UA2, observed chain A | TPO, ATP, and a 12-residue sequence-supported internal loop; native ligand parameter parity | 4,825 | Completed: 58,392 atoms |

All three preparations used pH 7, seed 2026, the normal missing-heavy-atom and sidechain controls, retained supported partners, and preserved original source files. The combined workers checked 1,013 protein/modified-residue stereocenters with no violations. The 1UA2 loop and nine missing 1PE heavy atoms are labeled provisional, low-confidence modeled coordinates.

The matrix includes 28 standard amino-acid/protonation/cap/disulfide cases, six named ion states (Na, Cl, K, Mg, Ca, Zn), and ten organic ligand/cofactor cases: ethanol, acetate, methylammonium, a halogenated aromatic, phosphate, a chiral amide, ATP, ADP, NAD and FAD. Each ligand retained heavy coordinates and matched its selected formal charge; native Amber and converted XML energies/forces matched at three poses under the unchanged production thresholds.

The six boundary checks exercised boron chemistry outside the validated element set, an open-shell radical, iron-containing heme, covalently linked FAD, missing ATP heavy atoms without an explicit repair choice, and real capped-peptide preparation/solvation with five retained ions. Expected rejections remained explicit; no implicit removal or replacement was accepted.

## Limits

- GROMACS prepared-complex transfer remains unsupported; native GROMACS standard-protein stability is audited separately.
- No metal coordination or oxidation-state inference: Zn is the named divalent Amber TIP3P nonbonded model.
- SEP/TPO/PTR use fixed -2 states, not protein pKa calculations. HYP has only its documented template forms; MSE, ALY, and other unregistered amino acids are blocked, not mutated.
- Heme/iron, covalent cofactors, unsupported elements/radicals and ambiguous graphs remain outside automatic preparation.
- Only explicitly selected eligible incomplete ligands and sequence-supported standard internal loops are built; modeled atoms are low confidence. Terminal sequence extensions are not built.
- 6A93 at 1 nm padding exceeds the 100,000-atom default profile. Its failed water job is retained; no partner was trimmed and no atom cap or padding bound was relaxed.
- Organic matrix fixtures use fixed-seed SMILES or CCD ideal coordinates, not experimental complexes. Their 5 Reference-platform steps (0.001 ps, vacuum) establish finite numerical continuity only.
- Complex CIFs are retained previously fetched RCSB inputs copied into the audit, not fresh web retrievals. Fresh network fetches and uploads are audited separately.

## Evidence

- [Machine-readable summary](summary.json), [matrix](chemistry-matrix.json), [boundary outcomes](chemistry-boundaries.json), [raw complex outcomes including failed water attempt](complex-preparation.json).
- [Exact source inspection and selections](complex-inspection.json), [environment and source hashes](environment.json), [artifact SHA256 index](artifacts-sha256.json).
- Drivers: [matrix](run_matrix.py), [complex workers](run_complexes.py), [this summary](summarize.py). Native logs and parameters remain under the isolated data root recorded in the artifact index.
- Complete backend regressions, UI functions, uploads/fetches, and native engine stability are separate parts of the parent final audit; their totals are not duplicated here.
