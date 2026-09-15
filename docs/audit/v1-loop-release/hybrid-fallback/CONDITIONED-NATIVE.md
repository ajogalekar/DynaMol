# Independent screen of source-conditioned native construction

**Source-conditioned native CCD produces three 1UA2 candidates without modeled Ramachandran outliers, but none of the six 1UA2 candidates passes the existing local screen. Every one moves observed ILE43 O beyond 1 Å. All six 8K5R attempts fail native initialization. The three masked-control candidates pass the local backbone screen.** These are construction results before full-complex refinement, not accepted preparation outputs.

The frozen comparison contains 15 attempts: original stems, seeds 2026/2027 for 1UA2 and 8K5R, and seed 2026 for the masked control. It tests source-conditioned torsion sampling with DirtyCCD (`torsion-source-dirty`), source-conditioned torsion sampling with torsion-aware CCD (`torsion-source-ccd`), and fragment insertion with torsion-aware CCD (`fragment-source-ccd`). FragmentSampler itself does not accept the boundary-condition fields; in that branch the intervention is the closer. Each attempt resets the model to the original missing inventory. Native MC steps, score subset/weight retrieval, cooling, candidate count, geometry thresholds and context limits are unchanged.

## Matched 1UA2 original-stem comparison

The loop is source A:44–55, with ILE43 and ASN56 as the original stems. Each default below is the already saved candidate from `hybrid-native-mc-v1` with the same source, sampler mode, seed and original span. The two torsion protocols share their matched default. No baseline candidate was regenerated or overwritten.

| Protocol | Seed | Modeled outliers, default → new | Boundary outliers, default → new | Largest observed displacement, default → new (Å) | New local result |
|---|---:|---:|---|---:|---|
| Torsion + source fields + DirtyCCD | 2026 | 3 → 6 | ARG57 → ASN56, ARG57 | 1.848 → 1.022 | Fail |
| Torsion + source fields + DirtyCCD | 2027 | 5 → 1 | ARG57 → ASN56, ARG57 | 1.159 → 1.672 | Fail |
| Torsion + source fields + CCD | 2026 | 3 → 0 | ARG57 → ARG57 | 1.848 → 1.941 | Fail |
| Torsion + source fields + CCD | 2027 | 5 → 0 | ARG57 → ARG57 | 1.159 → 2.162 | Fail |
| Fragment insertion + CCD | 2026 | 1 → 0 | ARG57 → ARG57 | 1.838 → 1.885 | Fail |
| Fragment insertion + CCD | 2027 | 1 → 1 | ASN56, ARG57 → ARG57 | 2.199 → 2.066 | Fail |

All six new candidates pass gross backbone geometry including outer joins and preserve the original cis/trans basin of every evaluated observed-only peptide link. All four source-CCD candidates have allowed or favored ASN56. ASN56 was already allowed/favored in three of the four distinct matched defaults; only default fragment/2027 had an ASN56 outlier. Supplying source fields while retaining DirtyCCD creates an ASN56 outlier in both tested outputs, so conditioning the sampler alone is not a consistent improvement.

Every new 1UA2 displacement maximum is **ILE43 O**. The largest N/CA/C displacement per candidate ranges from 0.094 to 0.138 Å. This is therefore a carbonyl-orientation problem in these raw candidates, not a large error in the N/CA/C frame fit. The unchanged full N/CA/C/O cap correctly catches it.

ARG57 was an outlier in the original source, but its defining preceding carbonyl changes in these candidates. It must be distinguished from a wholly unchanged source outlier; its presence is not proof that the repair created its original outlier status. The source-versus-repaired distinction is recorded in [SOURCE-BOUNDARIES.md](/Users/ashujo/Documents/Science/DynaMol/docs/audit/v1-loop-release/hybrid-fallback/SOURCE-BOUNDARIES.md). These results do not justify changing its classification policy or any acceptance threshold.

## Other attempts

All six 8K5R source A:177–179 attempts report `Failed to initialize monte carlo sampling protocol!`, including both seeds for each of the three protocols. Their matched original-stem default torsion and fragment attempts also failed initialization. There is no candidate to assess, and no new retries or context expansion were performed.

| Masked control, seed 2026 | Default local result | New local result | New observed displacement maximum (Å) | Modeled backbone RMSD to held-out coordinates (Å) |
|---|---|---|---:|---:|
| Torsion + source fields + DirtyCCD | Pass | Pass | 0.501 | 0.763 |
| Torsion + source fields + CCD | Pass | Pass | 0.313 | 0.387 |
| Fragment insertion + CCD | Fail: VAL5 boundary outlier | Pass | 0.538 | 0.806 |

The control's masked coordinates were excluded from generation and scoring. Its native fragment database may contain homologous or identical fragments; this is not a held-out-database accuracy benchmark. Three local control passes do not establish complete physical acceptance or generalization to experimental gaps.

## Independent identity and coherence checks

All 15 frozen source/input/implementation/parent-plan hashes verified before screening and remained identical afterward. All saved baseline candidate hashes also verified unchanged. Every new candidate binds the expected source path/hash, exact modeled identities, original outer stems and two observed context residues; none declares an observed internal context residue. Candidate IDs were recomputed from the recorded settings and backbone coordinates.

Boundary values were independently recalculated from original PDB atom identities and compared with the recorded native values, including radians/degrees consistency. The maximum native-versus-double-precision difference is 0.00018834°. Neighbor identities are K/R around 1UA2 ILE43/ASN56 and S/Q around 8K5R LEU176/SER180. The control explicitly retains the native n-terminal default because its first stem has no preceding residue; its observed c-stem boundary uses LYS6. There is no silent terminal default for either real target.

All nine candidates report completed native sidechain reconstruction. The independently reconstructed backbone after sidechain export matches the raw candidate backbone exactly (maximum difference **0 Å**), including modeled-atom coverage and unchanged observed atoms supplied by their source coordinates. No exported observed heavy-atom change lies outside the declared context. This establishes local coordinate coherence; it does not validate sidechain energetics, ligands, modified residues or the full-complex force field. Their complete source files remain hash-bound, and full-complex evaluation is separate.

## Evidence and limits

Machine-readable result: [conditioned-native-summary.json](/Users/ashujo/Documents/Science/DynaMol/docs/audit/v1-loop-release/hybrid-fallback/conditioned-native-summary.json). Cache root: `/Users/ashujo/.cache/dynamol-research/v1-loop-release/boundary-conditioned-native-v1`.

- `screen_conditioned.py` calls the existing `screen_native_comparison.screen(path, case)` directly for all nine produced candidates. It does not call that module's hardcoded main. Each new candidate directory contains its independent screen; old reports remain unchanged.
- `independent-screen-summary.json` has SHA256 `fe457dc7e8231b51e9806ceb99431fada4769cfd744e883cc8d4c58522d952ed` and matches the repository summary. It records attempt failures, exact hashes, source-boundary checks, matched defaults and all affected reference rows.
- The first audit wrapper stopped before screening because it expected the sidechain status spelling `succeeded`; the native schema says `complete`. The corrected audit assertion, original wrapper and integration note are preserved. This was an audit-script mismatch, not a generation failure.
- The completed independent screen used the isolated research Python and native CCTBX reference: 0.161 seconds, 98,713,600 bytes peak RSS. It performed no search, refinement, MD or QM. The parent-owned native generation batches totaled 16.753 seconds with no supervisor stop; sampled peak RSS was 822,919,168 bytes.

No app code, shared checkpoint, source coordinate, preserved result or acceptance threshold was changed by this screen. The parent is evaluating all six 1UA2 candidates with the full-complex validator regardless of these local outcomes. Raw Rama improvements must not be presented as a successful loop repair before that validation.
