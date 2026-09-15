# Independent saved-complex and accepted-handoff audit

The two saved 1UA2 results (`04-conditioned`, `08-shifted`) and the separate accepted-handoff endpoint pass the recomputed private static checks. No false acceptance or saved-file parameter change was found in this bounded audit. This is a milestone for the archived complete complex and its 12-residue A44–55 repair; it does not qualify the application, establish a general loop success rate, prove the predicted conformation, or validate dynamics.

No preparation, refinement, sampling, MD, or QM was run here. Original coordinates, candidates, parameters, and reports were read-only. Sources and results were hash-checked before and after the audit. Independent selection and exact identity/coordinate comparisons were combined with fresh invocations of the existing geometry, stereochemistry, and native CCTBX reference implementations; these are recomputations, not independent implementations of those scientific reference methods.

| Saved artifact | Atoms / indexed bonds | Selected residues | Largest observed-heavy displacement | Largest PDB rounding difference |
| --- | --- | --- | --- | --- |
| `04-conditioned` | 4,825 / 4,886 | A43–56 (14) | 0.754203 Å | 0.000841 Å |
| `08-shifted` | 4,825 / 4,886 | A43–56 (14) | 0.602455 Å | 0.000842 Å |
| Accepted handoff, torsion-source-ccd 2027 | 4,825 / 4,886 | A43–56 (14) | 0.754564 Å | 0.000860 Å |

All three retain the complete atom identities, elements, and indexed bonds after PDB reload. Their selected geometry and stereochemistry checks pass, no gross environment collisions are reported, and all 14 selected residues are scored without native CCTBX outliers before and after serialization. No acceptance threshold was changed.

## Exact boundary support and retained components

ARG57 remains a pre-existing native CCTBX outlier, with score 0.000296071995799. Its defining atom coordinates are **bitwise identical across the original PDB, archived initial coordinates, refined NPZ, and reloaded PDB**, in all three artifacts:

| Torsion | Exact defining atoms | Angle |
| --- | --- | --- |
| φ57 | ASN56 C; ARG57 N, CA, C | 67.7190902151° |
| ψ57 | ARG57 N, CA, C; THR58 N | −83.3662937518° |
| Preceding ω57 | ASN56 CA, C; ARG57 N, CA | 171.1305956392° |

The four fixed flank atoms (ILE43 N/CA and ASN56 CA/C) are also exact. Independent selection using actual φ/ψ/preceding-ω defining-coordinate changes and modeled definers selects A43–56, excluding ARG57 because its support is unchanged. This exclusion does not overlook a newly changed outlier. All observed heavy atoms outside the modeled residues and declared flanks are exact between archive and refined NPZ.

TPO A170 retains 11 heavy atoms, six hydrogens, and 18 incident bonds; ATP E381 retains 31 heavy atoms, 12 hydrogens, and 45 incident bonds. Original heavy-atom identities, archived hydrogen names/parent bonds, and all saved indexed bonds are retained. Archive→refined heavy coordinates are bitwise exact; original PDB→saved PDB heavy coordinates are also bitwise exact. PDB-versus-CIF parsing has only a maximum 8.9×10⁻¹⁵ Å numerical difference for ATP. Hydrogen **positions** are not fixed: existing refinement moves them, including ATP hydrogens by up to approximately 1.90 Å. Preserved hydrogen identity/connectivity must not be described as unchanged hydrogen coordinates.

## Fresh physical System and serialization

A freshly built physical System from the archived topology exactly matches the recorded native-System XML hash in all three results. All 4,825 particle masses are positive. Its forces are HarmonicBondForce, PeriodicTorsionForce, NonbondedForce, CMMotionRemover, and HarmonicAngleForce, with no construction-force classes or construction zero masses; this particular replay has zero constraints.

Fresh creation from the saved PDB produces a different XML byte hash. The difference is **only harmonic-bond term ordering**: exact multisets of all 4,886 indexed bond terms, including every parameter and term multiplicity, match. After removing that force, all remaining System fields and force terms are canonically identical. No coordinate-energy replay was needed to distinguish ordering from changed parameters.

- Archived fresh/recorded System SHA256: `b9f070712bb5d0307089e243237d9f1db0aa863b64f454c54c8f66a21ffa57d2`.
- Saved-PDB fresh System SHA256: `5269f48ebdedca5530d2a84ce939f7ccc0b83278beca2910c20418c54cb465ab`.

The replay uses the frozen current ff14SB/TIP3P files and the archived phosphoresidue/ATP parameter files. It establishes preservation and consistent instantiation of these parameters, not their general scientific accuracy or equivalence to an unrecorded historical protein/water force-field distribution.

## Separate accepted-coordinate handoff

The saved runner evidence confirms four complete-check rejections followed by acceptance of the fifth candidate, with 15 progress events. All eight fixed check names are present and true for acceptance. The returned path is the validated `refined.npz` inside the accepted attempt directory, and its SHA256 matches both the runner and validator manifests: `99e31cc3175ff7d0a8d6b606e4e4b3fb2456786c15576504bc9945ae3d8ecb23`.

The repeated endpoint was audited independently rather than inheriting the first endpoint's results. Its seed heavy coordinates match `04-conditioned` bitwise, while seed hydrogen coordinates differ by up to 0.001083 Å. Final coordinates differ by up to 0.578652 Å at modeled ASP53 OD1, with all-atom RMS difference 0.033587 Å. The first detectable difference therefore precedes refinement in hydrogen initialization; this comparison does not isolate every cause of final variation. Neither NPZ timestamp differences nor bitwise-reproducible minimization should be claimed. Both distinct final endpoints pass their own checks. This known-outcome handoff regression is not a new success-rate estimate.

## Evidence and audit corrections

- [Main detailed audit](audit-summary.json), [separate handoff audit](handoff-audit-summary.json), and [runner/variation evidence](handoff-variation.json).
- Reproducible read-only scripts: [audit_saved.py](audit_saved.py), [audit_handoff.py](audit_handoff.py), [audit_handoff_variation.py](audit_handoff_variation.py).
- Main evidence cache: `/Users/ashujo/.cache/dynamol-research/v1-loop-release/outer-anchors-independent-v1`; handoff cache: `/Users/ashujo/.cache/dynamol-research/v1-loop-release/handoff-independent-v1`.
- Input artifact caches: `outer-torsion-anchors-v1` and `accepted-coordinate-handoff-v1` under the same `v1-loop-release` cache root. Exact source, implementation, result, and artifact hashes are in the JSON reports.

The audit initially stopped on three overly strict reporting assertions: XML byte-order equality, JSON list versus in-memory tuple identity comparison, and PDB/CIF floating-point parser equality for ATP. These were corrected in the audit only after inspecting the underlying indexed terms, normalized identities, and numerical differences. The first script snapshot, both System XML files, and their ordering diff were preserved in the main evidence cache. No candidate coordinates, application validator, or acceptance policy was altered. Final native audits took 2.31 seconds (both initial artifacts) and 1.31 seconds (handoff), with peak sampled process RSS below 350 MB.
