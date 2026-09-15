# Native artifact and atom-mapping review

The research transport helpers preserve native parameters and coordinates, with explicit display permutations for fused covalent residues. This review fixed several concrete integrity gaps. It does not certify the source chemistry or make the models ready for application simulations.

## Fixes and evidence

- A saved measurement could previously reuse in-range indices after the map changed. `measurement_to_native` now requires `expected_mapping_sha256`, persisted when the atoms were selected. A newly generated current hash must never be substituted for a saved selection's missing hash.
- The mapping hash now includes complete native atom labels as well as source identities, element inventory, bonds and display order. Renaming native atoms with unchanged counts/elements therefore invalidates the old binding. Optional explicit `native_identity` records are checked against the stated particle index.
- Contradictory structured and legacy source IDs fail instead of silently choosing one. Empty, duplicate, malformed and out-of-range measurement selections are rejected. Complex-valued arrays are not accepted as physical coordinates.
- Native snapshot loading and copying require a pinned SHA-256 manifest digest. `None` can no longer bypass the pin. Boolean or floating-point native indices, malformed manifests, modified files, escaping paths and symlinked artifacts are rejected. Read-only inventory inspection can still explicitly omit a pin.
- Native/display coordinate inversion and cross-residue covalent bonds remain exact in the interleaved synthetic adduct control. The independent [4,534-atom 1CA2 check](identity-review-1ca2-roundtrip.json) verifies bitwise coordinate roundtrip and identical selected distances using the preserved native snapshot and source-map file.

The [identity review suite](identity-review-tests.log) passes 40 tests. The final [combined identity/mechanics suite](identity-and-mechanics-review-tests.log) passes 63 tests, including the explicit GROMACS zero-amplitude omission case. Current null-term policy is documented in [bonded-site validation](BONDED-SITE-VALIDATION.md): exact-zero torsions may be omitted, but positive periodicity, all nonzero terms and the complete independent 1–4 exception graph remain required.

## Before production lifecycle integration

1. Bind the source atom map, native parameter snapshot, chemical-state specification and each engine's particle inventory together as required artifacts. The current snapshot hashes optional map files if supplied, but does not require a source map or automatically join one to a `ResearchModel`.
2. Independently verify observed source identities/elements/coordinates and generated-atom provenance against the original structure and declared reaction/metal model. A complete, hashed map proves what was recorded; it does not prove that the recorder assigned every atom correctly. `chemical_state` and `source_structure` are preserved data, without automatic scientific consistency validation.
3. Persist the map/snapshot binding alongside frames and saved measurements through save, reopen, archive, engine export and preparation changes. Coordinate shape alone cannot distinguish two same-sized models. Explicitly migrate stable source-ID selections after preparation changes or report missing atoms; do not reuse bare indices. Only display coordinates may be permuted: the native System and restart must retain native order.
4. Keep display topology/System objects from being mutated behind an accepted report. The research helpers provide hashes, not operating-system immutability or automatic invalidation of every in-memory change. Revalidate modified native parameters and regenerate the associated bindings.
5. Integrate the existing runtime and scientific validation gates: finite forces, explicit site checks, energy/force conversion parity, protonation/adduct/metal-state provenance and held-out torsion-profile evidence. A serialized roundtrip or intact bond does not validate physical motion. The current `research_prototype` stage and `simulation_ready: false` remain unchanged.

This review intentionally does not add production storage/UI integration or silently infer source identities from geometry.
