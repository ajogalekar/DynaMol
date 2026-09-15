# Mg–GDP model assessment, 14 September 2026

The completed 6OIM baseline retained its six initial magnesium donors, but the
GDP O2B contact shortened from 2.157 Å in the input to a 1.874 Å sampled mean.
This is a reason to investigate the interaction model, not evidence of an
accurate metal site. No replacement parameters have been accepted or inserted
into the app.

The actual full-complex topology contains the bundled Li/Merz TIP3P 12–6 CM
magnesium parameters: charge +2, Rmin/2 1.360 Å and epsilon 0.01020237 kcal/mol.
GDP O2B has charge −0.9552 e and scoped type QR, preserving native GDP type O3.
Its Rmin/2 is 1.6612 Å and epsilon 0.2100 kcal/mol. The topology has no C4 terms.
These are verified loaded values, not an inferred choice based on a filename.
The saved evidence is `mg-gdp-model-assessment-v1/loaded-model.json` under the
focused research cache, including hashes of the topology and bundled source.

Two primary-source alternatives merit comparison, with explicit applicability
limits:

- Panteva, Giambașu and York's modified 12–6–4 model adjusts particular ion–solute
  interactions because the original model overbinds phosphate sites. Their
  reported calibration uses TIP4P-Ew and defined nucleic-acid models. It cannot
  be copied directly into the present TIP3P/Meagher GDP combination. Adding a C4
  term alone is not the published correction. [Original paper](https://theory.rutgers.edu/resources/pdfs/Panteva_JChemPhysB_2015_v119_p15460.pdf).
- Grotz, Cruz-Leon and Schwierz's microMg provides a TIP3P-compatible 12–6
  candidate with explicit pair overrides, and an author GROMACS implementation.
  Their phosphate calibration changes DMP oxygen charges, types and angles to
  match the RNA force field. Consequently its success for RNA does not establish
  transferability to the retained GDP model or protein donors. The separate
  nanoMg variant deliberately accelerates water exchange; that is not the
  default choice for assessing physical motion here. [Original paper](https://doi.org/10.1021/acs.jctc.0c01281),
  [author implementation](https://github.com/bio-phys/Magnesium-FFs).

The author implementation was retrieved at commit
`e71379d1a88b54221018af53ca54d3c18066ba05` into the cache, with exact file hashes
in `retrieval.json`. It includes chloride and RNA pair overrides, not merely a
new magnesium radius. The retrieved root has no explicit license file; this
snapshot is local research evidence, not a bundled redistribution decision.

Next establish a coherent nucleotide/ion/water combination and assess protein
donors, then test the candidate against the retained baseline with explicit
pairwise transport checks in both engines. Do not tune a radius to reproduce
one crystal distance, change GDP charges without a charge model, silently drop
pair overrides, or describe a short stable run as a binding-affinity validation.
