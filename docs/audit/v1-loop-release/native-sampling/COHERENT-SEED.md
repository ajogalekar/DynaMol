# Coherent native fragment closure: bounded 8K5R comparison

**No candidate passed all unchanged local backbone checks.** Starting from a coherent native fragment fixes the inherited peptide-join defects seen with the transplanted seed, but this single-fragment KIC comparison did not find a backbone close enough to the original observed atoms.

| Case | Raw native solutions | Within 1 Å observed N/CA/C/O cap | Gross local backbone geometry pass | Local CCTBX pass | Combined pass |
|---|---:|---:|---:|---:|---:|
| Intact control | 12 | 4 | 12 | 4 | 4 |
| Original stems 176–180 | 0 | 0 | 0 | 0 | 0 |
| Extended context 175–181 | 36 | 0 | 36 | 3 | 0 |

All 36 extended-context solutions preserve the original observed-only peptide cis/trans basins. The closest candidate still moves SER180 O **4.010 Å**, N 2.141 Å, C 1.967 Å and CA 1.581 Å; LEU176 O moves 1.340 Å. The 1 Å cap was not relaxed. The original-stem seed required a 2.233-radian native psi rotation to align its observed n-stem carbonyl before closure, then yielded no solutions.

The coherent seed was generated with the preserved generic context worker (`400f2e4cf5860773ebb838c64671958d0427d8e48159c5ba9b5007a7aff28975`) on private copies of the frozen 8K5R input and context. Missing atoms and `research_context_atoms` were exported together before any original-coordinate restoration. Both proposed local spans passed bond, angle and omega checks before closure. An independent export replay reproduced the seed PDB byte for byte.

The native-selected fragment spans source residues **176–181**; its native internal numbering is 172–177. We kept the previously frozen original 176–180 and extended 175–181 comparison. The latter contains only two observed internal residues, LEU176 and SER180. Its four internal pivot sets use 24 total native coil `Draw` calls (12 conformers × two non-pivot residues), seed 2026. All input and solution coordinates, unsuccessful trials, native versions and hashes are retained. The native construction worker took 1.36 s and KIC took 0.27 s, each under a 150 s supervisor limit.

The independent CCTBX screen uses original observed coordinates outside the candidate span and identifies affected boundary torsions from the original observed baseline. The 0.0001 Å movement tolerance only separates native floating-point coordinate roundoff when selecting affected residues; it does not alter the 1 Å movement cap or any geometry/reference threshold. Observed cis/trans state is checked against the original source, not inferred from the changed seed.

These results support changing how candidates are generated and matched to the complete observed stem geometry. They do **not** show that the gap is impossible to model, nor justify increasing the movement cap. This was one native fragment and one small torsion-sampling comparison. Sidechain placement, chirality, complete retained-environment contacts, force-field refinement and serialized full-complex acceptance were not validated. No app or complete protein candidate was admitted.

[Compact results](coherent-seed-summary.json), [exporter](export_coherent_backbone.py), [native KIC sampler](sample_8k5r_kic.py), and [independent screen](screen_kic_backbones.py). Complete immutable inputs, implementation snapshots and results: `/Users/ashujo/.cache/dynamol-research/v1-loop-release/8k5r-coherent-native-sampling-v1`.
