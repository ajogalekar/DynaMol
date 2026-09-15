# Paired 1UA2 loop-refinement replay

**Neither protocol passed the combined checks.** Protecting an initially
allowed backbone conformation by itself is not a sufficient loop fix: in this
12-residue example it removes modeled-residue reference outliers while leaving
a bad boundary and a distorted peptide join.

The comparison used the exact archived 4,825-atom refinement input, including
its saved hydrogens, ATP and phosphorylated threonine. Both protocols loaded
the same archived, provenance-checked native parameter bundle. All 2,323
nonmodeled heavy atoms remained bitwise fixed. The modeled sequence is
A44–55, KLGHRSEAKDGI; the reference screen also includes its immediate peptide
neighbors ILE43 and ASN56. The seed already has GLU50 and ASN56 outliers.

| Protocol | Modeled reference outliers | Boundary outliers | Gross geometry | Stereochemistry |
|---|---|---|---|---|
| Existing refinement | LYS44, GLY46, ARG48 | ASN56 | Pass | Pass |
| Temporary protection of allowed seed phi/psi | None | ASN56 | Fail | Pass |

The second protocol uses the previously defined research-only periodic
phi/psi tethers at 1,000 kJ/mol for initially allowed/favored selected torsions;
seed outliers remain unprotected. This was a paired comparison, not a force
constant search or physical parameter fit. The protected output has an
ILE43 O–C–LYS44 N angle of **62.508°**, outside the unchanged gross screening
range of 80–150°, and a nearby soft overlap. It is rejected despite its zero
modeled-residue Ramachandran outliers. Endpoint energies include temporary
construction forces and are not interpreted as native conformational energies.

Both saved PDB files preserve exact atom identity/order and the affected
reference classifications after round-trip loading. Both fresh native System
serializations are identical and contain no added construction torsion force.
The original-protocol replay reproduces the archived selected outlier
identities, but its final coordinates are not bitwise identical to the older
run (maximum atom displacement 0.208 Å); no exact archived-coordinate replay
claim is made. The two new protocols do have identical starting coordinates,
atom identities and native System hashes.

Native refinement took 1.69 and 2.47 seconds. Each child was supervised with
120 seconds, two CPU threads and a 4 GiB sampled-RSS limit; peak sampled RSS
was below 219 MB. Source hashing incurred a separate local filesystem read
delay before these supervised child runs. Source hashes were unchanged
throughout, and AST/targeted whitespace checks passed. No full preparation,
dynamics, force-field fitting, package build or app admission occurred.

The result reinforces the need to construct a coherent join against the
observed carbonyl and neighboring backbone geometry, then refine and check
the full retained environment. Keeping selected torsions near a seed cannot
compensate for an incompatible join. No threshold was relaxed and no failed
artifact was replaced.

[Compact evidence](conformation-refinement-summary.json) and
[reproducible comparison](replay_conformation_refinement.py). Full snapshots,
source hashes, coordinates, native Systems, native-reference reports and
supervisor logs are retained under
`/Users/ashujo/.cache/dynamol-research/v1-loop-release/1ua2-conformation-replay-v1`.
