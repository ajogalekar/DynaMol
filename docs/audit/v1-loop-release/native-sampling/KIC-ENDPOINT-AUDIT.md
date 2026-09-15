# Saved native KIC endpoint and identity audit

The large observed-atom displacements in the bounded 8K5R experiments occur **inside extended closure spans**, not at the actual KIC endpoints. Independent coordinate recomputation agrees with the installed API and the matching ProMod3 3.6.0 source. This audit found no endpoint targeting, source-identity remapping, or Å/nm error in those saved experiments. It does not establish that any failed loop can be repaired under the original 1 Å observed-backbone cap.

## Scope and evidence

This was read-only: 158 saved returned solutions were examined, comprising 122 8K5R solutions and three repeated 12-solution intact-control replays. The repeated controls are not independent validation samples. No new closure, sampling, refinement, or MD was run. Original coordinates and experimental artifacts were preserved.

The missing residues are A:177 ALA, A:178 LYS, and A:179 ASN. In the original A:176–180 span, A:176 LEU and A:180 SER are the two KIC targets. For the A:175–181 span, the targets instead are A:175 and A:181; LEU176 and SER180 are now observed **internal** residues, and the closure solver does not independently pin them. Other bounded pool contexts follow the same rule. Exact chain, residue number, insertion code, residue name, and atom-name identities are retained in the result.

Independent checks read raw PDB coordinates in Å, compare them with direct saved native coordinates, and reconstruct both endpoint rigid fits with NumPy Kabsch superposition. All identity checks passed. Maximum disagreement with the saved displacement records was 0.0000041 Å. Maximum disagreement between predicted and returned endpoint N/CA/C coordinates was 0.0000052 Å. Predicted oxygen positions agreed within 0.0000052 Å.

| Saved comparison | Solutions | Maximum endpoint N/CA/C displacement (Å) | Maximum endpoint O displacement (Å) | Maximum observed internal displacement (Å) |
|---|---:|---:|---:|---:|
| Intact control, each replay | 12 | 0.000004 | 0.000003 | 3.529426 |
| Transplanted seed, original stems | 2 | 0.000003 | 0.051171 | No observed internal residues |
| Transplanted seed, extended span | 44 | 0.000007 | 0.082068 | 4.643971 |
| Coherent seed, extended span | 36 | 0.000378 | 0.082114 | 4.301718 |
| Multiple coherent fragments, all returned solutions | 40 | 0.067402 | 0.160328 | 6.322489 |

The coherent original-stem case and four selected pool fragments returned no solutions; they are not silently counted as passing endpoints. Every returned endpoint remains below the original 1 Å cap. The intact control demonstrates that alternate valid closure solutions can preserve endpoints while moving observed interior residues substantially.

## Native contract and implementation

Installed `loop/backbone.hh:116–139` constructs `GetTransform(residue)` from N, CA, and C and calls `MinRMSDSuperposition`. `FillKICParameters` independently fits the first and last backbone residue to the supplied n-stem and c-stem. `ApplySolution` applies those frames to the outer fragments. Consequently n-stem C and c-stem N/CA, together with the other endpoint N/CA/C atoms, follow their endpoint's minimum-RMSD rigid fit; they are not individually forced to exact source coordinates when fragment and target internal geometry differ. Pivots are strictly internal. [ProMod3 3.6.0 KIC source](https://git.scicore.unibas.ch/schwede/ProMod3/-/raw/3.6.0/modelling/src/kic.cc)

Oxygen is absent from that fit. The n-stem O follows the native frame after the experiment's explicit carbonyl-orientation preparation. The c-stem O is then explicitly reconstructed by `ReconstructCStemOxygen` when a valid following residue exists: native code places it 1.230 Å from C along the normalized bisector opposite CA and the following observed N. The resulting small c-stem O displacement is therefore expected native behavior, not a later atom remapping. The independent audit reproduces this operation from saved coordinates; it does not modify coordinates. [ProMod3 3.6.0 backbone source](https://git.scicore.unibas.ch/schwede/ProMod3/-/raw/3.6.0/loop/src/backbone.cc)

Code inspection confirms that `sample_8k5r_kic.py` records `GetN/GetCA/GetC/GetO` coordinates directly after `Close`. The screening adapter overlays those coordinates by exact atom identity and divides by 10 once for OpenMM's nm coordinates. No downstream rigid transform changes the saved closure endpoints. The matching release source and installed headers are hashed; this is source/API and coordinate-behavior evidence, not a binary reproducible-build attestation.

## Consequence for the next step

Keep the observed N/CA/C/O cap unchanged. Further work should distinguish endpoint fitting from movement of observed residues made internal by context extension, and retain coherent peptide/carbonyl geometry across joins. Endpoint correctness does not rescue the existing failed candidates: neither Ramachandran acceptance nor a backbone-only pass constitutes complete all-atom acceptance. No loop-fix, full-complex acceptance, or release-readiness claim follows from this audit.

## Reproduction and artifacts

- Implementation: [audit_kic_endpoints.py](audit_kic_endpoints.py)
- Compact source hashes, exact target identities, and per-case maxima: [kic-endpoint-audit-summary.json](kic-endpoint-audit-summary.json)
- Full per-solution evidence: `/Users/ashujo/.cache/dynamol-research/v1-loop-release/kic-endpoint-audit-v1/result.json`
- Full result SHA256: `81c4cbd38b84e25a41cf105d162cb45c98e7586de368c9337e477aa0a8150a3b`
- Archived upstream source: `source-release-3.6.0.cc` (SHA256 `6f5d88c5b8c52508d56c60f4db8c57980c13d129fe717680fc7cd3b465d7095a`) and `source-release-backbone-3.6.0.cc` (SHA256 `cd516c263fa9b674de62e7fa1949725a40cb31d09340ab2133ca2c1b6a58f937`) in that same cache directory. Fetch metadata, including unsuccessful earlier attempts, is preserved.

Run the audit from the repository with the isolated research Python and `OPENBLAS_NUM_THREADS=2 OMP_NUM_THREADS=2`. It uses only NumPy and saved evidence, and does not import the native sampler or modify previous results.
