# Modified amino-acid preparation

DynaMol retains supported modified amino acids in the protein chain. It does not
replace their sidechains with a standard parent residue or parameterize them as
free, disconnected ligands.

| Residue | Chemistry | Template | Internal net charge |
| --- | --- | --- | --- |
| SEP | Phosphoserine | phosaa14SB SEP | −2 e |
| TPO | Phosphothreonine | phosaa14SB TPO | −2 e |
| PTR | Phosphotyrosine | phosaa14SB PTR | −2 e |
| HYP | 4-hydroxyproline | ff14SB HYP | 0 e |

Phosphoresidues use fixed dianionic templates. Selecting a preparation pH does
not calculate phosphate pKa values or automatically change these states. The
inspection card and preparation record disclose the choice. C-terminal HYP uses
CHYP (−1 e); unsupported terminal phosphoresidues, N-terminal HYP and additional
covalent crosslinks require matching templates. All supported modifications
currently use OpenMM with explicit water.

Unknown modifications and unnatural amino acids remain visible and receive a
specific missing-template message. Sequence-supported gaps also retain their original residue names; missing modified or unknown amino acids require external modeling and are never rebuilt as standard parent residues. In particular, the presence of MSE or ALY in
an installed Amber library does not by itself establish that DynaMol can prepare
it: those require the separate ff14SB_modAA parameter set and an independently
checked integration. No generic force-field parameters or parent substitution
are silently applied.

## Sources and coordinate repair

The phosphorylated residue force field is vendored from the official
[OpenMM forcefields repository](https://github.com/openmm/openmmforcefields/blob/3f73b6ff9730e4e356fc6a0dbde6039517d59412/openmmforcefields/ffxml/amber/phosaa14SB.xml),
commit `3f73b6ff9730e4e356fc6a0dbde6039517d59412`, converted upstream from
AmberTools 24.8. It cites Raguette et al., *phosaa14SB and phosaa19SB: Updated Amber
Force Field Parameters for Phosphorylated Amino Acids*, JCTC 20, 7199–7209 (2024).
The upstream MIT notice accompanies the raw XML. HYP uses OpenMM's bundled
`amber14/protein.ff14SB.xml`.

Raw sources, the normalized XML, CCD records and hashes live in
`backend/data/modified_residues/`. Normalization changes atom-type identifiers to
match the OpenMM protein namespace and imports the exact phosphorus type and
Lennard-Jones entry omitted from protein-only ff14SB. It preserves the numeric
parameters and supported residue graphs. Reproduce it with:

```bash
.venv/bin/python scripts/generate_modified_residue_xml.py
```

PDBFixer receives explicitly registered heavy-atom templates from vendored RCSB
CCD ideal coordinates for [SEP](https://www.rcsb.org/ligand/SEP),
[TPO](https://www.rcsb.org/ligand/TPO), [PTR](https://www.rcsb.org/ligand/PTR), and
[HYP](https://www.rcsb.org/ligand/HYP). These coordinates support missing-atom
placement; they are modeled geometry, not experimental reconstruction evidence.
Hydrogen names/bonds come from the exact force-field templates. Stereo checks
include phosphate-residue CA, TPO CB and HYP CG, including insertion codes in
residue identity.

## Native Amber comparison and portability

Native TLEAP capped peptides (`ACE–residue–NME`) are retained in
`docs/audit/modified-residues/native/`. Parameter/source hashes and the original
TLEAP input/logs accompany them. Twelve coordinate comparisons (four residues,
three poses each) pass against native Amber energies and forces. The regression
suite also reverses atom and bond order, repairs a deliberately removed modified
heavy atom, adds hydrogens, checks charges and stereo, and verifies System XML
serialization.

The test found that generic OpenMM improper matching can permute two aromatic
neighbors in PTR relative to native TLEAP. DynaMol's
`register_modified_forcefield(forcefield)` restores the four native PTR improper
quartets by residue and atom name after matching. It leaves their numeric
parameters unchanged. This correction is required when rebuilding from the XML
files; the XML files alone do not encode the Python generator. The fully assembled
`modified-residue-system.xml` in a prepared dataset stores all corrected terms for
that unsolvated, nonperiodic preparation topology. It is not the final periodic MD
system; completed OpenMM jobs save their assembled MD model as `system.xml`.
The registry's source-code checksum is recorded in preparation provenance.

Reproduce the checks with:

```bash
.venv/bin/python docs/audit/modified-residues/verify_native.py
.venv/bin/pytest -q backend/tests/test_modified_residues.py
```

The [machine-readable native comparison](audit/modified-residues/native-parity.json)
reports the actual deviations. The thresholds (0.002 kJ/mol energy and
0.02 kJ/mol/nm per force component) accommodate native Amber's finite-precision
angular serialization. Passing these checks establishes implementation agreement
for the tested cases. It does not establish a correct bound protonation state,
rotamer, missing loop, coordination model, or converged simulation.

## 1UA2 and worker continuity

Inspection of the original 1UA2 dataset recognizes all four TPO170 residues and
their peptide links, separately from its four ATP ligands. Each protein chain also
has a genuine 12-residue internal sequence gap between residues 43 and 56. These
gaps are within the current optional builder's limit of 12 standard residues per
internal gap and 96 total. Subsequent full monomer and four-chain preparation
checks rebuilt 12 and 48 residues while retaining ATP/TPO; see the
[loop validation record](audit/loop-repair/REVIEW.md). These are provisional
starting models, not verified native loops. No phosphate removal or artificial
bond across an unbuilt gap bypasses the preparation requirements.

A documented 13-residue fragment, 1UA2 chain A residues 164–176, passed the actual
background preparation, explicit-solvation and OpenMM worker pipeline. All 11 TPO
heavy atoms, both peptide links, phosphate bonds and stereocenters survived
preparation, including when free heterogen removal was selected. TPO charge was
−2 e; its heavy atoms moved 0.0 Å during preparation. The parameter-file hashes and
solute charges survived into the 5,687-atom solvated system and three finite
frames over 0.004 ps. This bounded software test is not full-1UA2 preparation or a
convergence study. See [the exact worker record](audit/modified-residue-native-run.json)
and reproduce it with `scripts/check_modified_residue.py`.
