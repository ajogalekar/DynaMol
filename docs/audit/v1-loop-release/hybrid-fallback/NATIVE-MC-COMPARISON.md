# Forced native torsion versus fragment-insertion Monte Carlo

**Neither native mode produced an acceptable research backbone for 8K5R or 1UA2 in this bounded comparison.** One masked intact-loop control passed the local backbone checks. All outputs remain research candidates; complete-complex refinement and all-atom acceptance are separate.

## Frozen comparison

Each attempt rebuilt the original ProMod3 raw model and targeted the exact missing source identities. No database candidate was inserted first, and prior gap clearance could not suppress either fallback. The comparison used the installed ProMod3 3.6.0/OpenStructure 2.11.1 runtime, its native MC scoring weights, simulated-annealing schedule, and `DirtyCCDCloser`.

- 8K5R: missing A:177–179 (AKN), source spans 176–180 and 175–181.
- 1UA2: missing A:44–55 (KLGHRSEAKDGI), source spans 43–56 and 42–57.
- Masked six-residue intact control: A:2–4 (QIF) removed from generation/scoring input; span 1–5. Hidden coordinates were used only for the later RMSD measurement.
- Two seeds, 2026 and 2027; one returned conformation per seed/span/mode; 5,000 native MC steps per attempt. Maximum 20 attempts, without outcome-driven retries.
- Torsion mode used `PhiPsiSampler` with native default boundary conditioning. Fragment mode used `FragmentSampler`, five initialization insertions and 100 fragments per position. Proposal fragments were length 3 for 8K5R/control and the native default length 9 for 1UA2. No profiles or hidden control coordinates were used in fragment search.
- Original stems or at most two observed internal context residues; the original 1 Å observed-backbone cap was unchanged.

These modes differ in their proposal mechanism but share native closure and backbone scoring. This is not an independent algorithm/library replication. `DirtyCCDCloser` does not impose its own allowed-Ramachandran check, making the independent reference gate necessary. The native defaults were retained for this comparison; it is not an optimized protocol.

## Results

| Case and mode | Frozen attempts | Returned candidates | Gross backbone including outer joins | Complete local backbone pass |
|---|---:|---:|---:|---:|
| 8K5R torsion MC | 4 | 2 | 2 | 0 |
| 8K5R fragment MC | 4 | 2 | 2 | 0 |
| 1UA2 torsion MC | 4 | 4 | 4 | 0 |
| 1UA2 fragment MC | 4 | 4 | 4 | 0 |
| Masked control torsion MC | 2 | 2 | 2 | 1 |
| Masked control fragment MC | 2 | 2 | 2 | 0 |

Both modes failed to initialize with the original 8K5R stems, for both seeds. All four failures are preserved. The extended 8K5R candidates exceed the observed-backbone cap and have local reference outliers. The closest of these moves an observed backbone atom 3.596 Å.

For 1UA2, the closest original-stem torsion candidate moves observed backbone 1.159 Å and has reference outliers. The extended fragment candidate with seed 2027 has no selected CCTBX outliers but moves observed backbone 4.214 Å. Neither qualifies as a passing backbone. These are distinct failure modes; native closure or reference classification alone is insufficient.

The passing masked torsion-control candidate (seed 2026) has maximum observed-backbone displacement 0.485 Å and modeled N/CA/C/O RMSD 0.740 Å from the withheld coordinates, without superposition. It is still not a validated complete all-atom candidate. One favorable control does not establish an accuracy distribution; the native database was not purged of homologous or identical fragments.

All 16 returned candidates completed native sidechain reconstruction. Reconstructed backbone coordinates exactly match the raw candidate coordinates. No moved observed heavy-atom residue was found outside the explicitly declared context. Original observed-only peptide cis/trans basins were checked through the outer joins and preserved in all returned candidates.

The six native batches completed in 19.95 seconds total; the longest took 9.29 seconds. Peak sampled RSS was 1,093,943,296 bytes (1.02 GiB). Each batch had a 180-second deadline, a sampled 4 GiB RSS stop and two-thread environment limits. No resource stop occurred. CPU sampling and these small cases do not establish desktop latency for a whole protein panel.

## Preservation and validation boundaries

The original complete selected structures are copied and hashed. The 1UA2 source retains A:170 TPO and E:381 ATP; its temporary scoring input aliases ATP to B:1. The 8K5R source retains A:186 TPO and C:401 VQE; its temporary scoring input aliases VQE to B:1. Native standardization of TPO to THR is recorded and never exported as a replacement. Native backbone scoring can omit nonprotein chemistry; full source preservation is not a claim of complete-complex native scoring.

Independent screening uses native CCTBX classifications, the existing broad bond/angle/peptide-planarity bounds, the original displacement cap, and exact source-identity/Å-to-nm checks. It covers modeled residues and every fully defined phi/psi affected by candidate backbone movement, including outer observed boundaries. Observed-only omega basins are compared with original coordinates. Unbuilt distant gaps are not connected artificially for scoring. Sidechain chirality, rotamers, retained-environment contacts, force-field refinement and final admission remain the root validator's responsibility.

## Adapter contract and artifacts

`native_mc_adapter.py` emits `candidate.json` before sidechain reconstruction and updates it with the optional native heavy-atom result. It exposes:

- `backbone_atoms`: exact `[chain, resid, insertion_code, residue, atom]` identities and both `xyz_nm` and raw `xyz_angstrom`.
- Separate `modeled_residue_keys`, `context_residue_keys`, `observed_internal_residue_keys`, and `outer_stem_residue_keys`.
- `modeled_heavy_atoms` and `observed_context_atoms` after native sidechain construction. Context contains explicit moved observed heavy atoms, including outer stems. All other observed/source chemistry remains in the original source file. Raw backbone and all context changes must be considered coherently; silently restoring only the observed atoms would repeat the earlier join defect.
- Original source hashes, sequence, generator/settings/seed, candidate ID, per-atom displacement, observed omega checks, and native parent-residue omissions. `app_ready` and `physical_model_validated` remain false.

Implementations: [generator](native_mc_adapter.py), [frozen plan and supervisor](run_native_comparison.py), [independent local screen](screen_native_comparison.py). Compact results: [native-comparison-summary.json](native-comparison-summary.json).

Complete cache: `/Users/ashujo/.cache/dynamol-research/v1-loop-release/hybrid-native-mc-v1`. It contains private source/input snapshots, frozen `plan.json`, each failed or successful attempt, coherent candidate coordinates, per-candidate independent reports, logs/resource records, native database/library hashes, and `environment-preservation-map.json`.

First passing local control: `masked-control/torsion_mc/native/1-5-seed-2026/candidate.json` under that cache. Root can feed these candidate files into the shared full-complex validator; no native candidate has been admitted automatically.

Frozen plan SHA256: `1a1d36c14a38439a93ecbc293c98a1539b2ced0a798b08d8311c26f8b1620054`.
Independent summary SHA256: `3d1105bde085d1ad5e3c2dc3690665440fb7806993985593bc7d4c14779a8188`.

The bounded result supports acceptance-driven fallback as a routing correction, but does not establish that these two native proposal modes solve the geometry/context mismatch. No further sampling was performed after screening, no app backend was changed, and no MD/QM was run.
