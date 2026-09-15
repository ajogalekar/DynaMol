# Torsion-only closure from the transplanted seed

This single bounded comparison produced **no acceptable 8K5R backbone**. The intact control reproduced a valid backbone in all four internal pivot configurations. Some 8K5R candidates pass native CCTBX and the unchanged 1 Å observed N/CA/C/O cap, but every candidate retains invalid backbone joins already present in the refinement-entry seed.

| Case | Raw KIC solutions | Within observed cap and omega basins | Also passes local CCTBX | Also passes gross backbone geometry |
|---|---:|---:|---:|---:|
| Intact control | 12 | 4 | 4 | 4 |
| Original stems 176–180 | 2 | 2 | 1 | 0 |
| Context 175–181 | 44 | 7 | 6 | 0 |

Counts are raw solutions, including repeated conformations reached with different pivots; they are not counts of independent validated models. The context comparison uses four internal pivot sets and 24 total native coil sampler `Draw` calls (12 conformers × two non-pivot residues), seed 2026. All draws, failed checks and coordinate solutions are retained. Carbonyl orientation was represented before closure, and observed-only peptide cis/trans basins were independently checked against the original observed source. No bond-length, angle or force-constant fitting was introduced.

The source refinement seed already contains an ASN179 C–SER180 N distance of 2.06665 Å, an omega of −72.18°, an O176–C176–N177 angle of 42.99°, and a C179–N180–CA180 angle of 66.73°. Torsion-only closure does not repair those inherited internal geometries. This motivates a separate, explicitly bounded comparison starting from a coherent native fragment before original observed atoms are restored; increasing sampling of this malformed seed is not justified by these results.

The native subprocess completed in 0.39 s under a 150 s supervisor limit. This is backbone feasibility evidence only; no sidechain rebuilding, full-atom environment check, chirality assessment, force-field refinement or app admission is implied. The screen uses native CCTBX classifications and broad backbone bond/angle/omega limits; it does not validate native loop conformation.

[Compact results](transplanted-seed-summary.json), [native implementation](sample_8k5r_kic.py), and [independent screen](screen_kic_backbones.py). Complete artifacts and implementation snapshots: `/Users/ashujo/.cache/dynamol-research/v1-loop-release/8k5r-native-sampling-v1`.
