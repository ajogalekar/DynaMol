# Pinned microMg applicability inventory — 14 September 2026

**The published microMg files are not a direct replacement for the complete
current KRAS–sotorasib–GDP–Mg model.** The main gap is the retained nucleotide
interaction model; merely changing Mg's Lennard-Jones parameters would omit
published pair corrections. No app parameter or source topology was changed.

This audit inventories one candidate at author repository commit
`e71379d1a88b54221018af53ca54d3c18066ba05`. The complete 80-file tree was retrieved
without truncation. The three previously retained files match their recorded
SHA256 hashes and the pinned Git blob identities. Nine additional example
files were verified against that tree. All 29 microMg pair rows, including
self and chloride, agree numerically between the root file and both examples.

The complete corrected 6OIM topology contains 72 atom types: 2,693 protein
atoms, an 82-atom covalent drug residue, 40 GDP atoms, one Mg, seven sodium
counterions and 24,654 water atoms. The inventory records every component/type,
its native LJ values and charge range, source aliases, published pair lookup
and any ordinary mixing-rule fallback. A fallback is mathematically defined;
that does not establish its physical accuracy.

| Component | Exact finding | Applicability |
| --- | --- | --- |
| GDP | 21 scoped types match none of the raw override names. Undoing the documented aliases identifies 19 named counterparts; **CK and O3 remain absent**. | No automatic alias assignment made. |
| Coordinating GDP O2B | Native `QR → O3`, charge **−0.9552 e**. Published RNA phosphate oxygen is `O2`, **−0.7760 e**. Their LJ values match at source precision, but chemistry/charge do not. | Renaming O3 to O2 would conceal the unresolved model change. |
| Ser17 OG | Native `OH`, −0.6546 e, matching the source protein template. A published `mMg–OH` override exists; RNA ribose OH has −0.6139 e. | Numeric coverage exists; this does not validate the complete protein coordination site. |
| Remaining protein | 15 of 28 native type names occur in the RNA override list; 13 do not. | A raw include would alter only part of the protein's Mg interactions. |
| Covalent drug residue | Four of 25 types have matching override names; 21, principally GAFF types, do not. | No demonstrated complete-adduct transfer. |
| Sodium | Native `Na+` has no calibrated pair in this source. | The published chloride correction cannot be substituted for it. |
| TIP3P | Charges agree. Current OW σ is **0.3150752407 nm**, ε **0.6359679988 kJ/mol**; the pinned example uses **0.315061 nm**, **0.636386 kJ/mol**. | Same water family, small but real numerical variant; no silent replacement. |

The source's Ser17-compatible named pair has σ = 0.23238044375 nm and
ε = 4.61046323173 kJ/mol. Those values are copied into the inventory as source
evidence, not assigned to DynaMol. The source atom-type table's zero charge
fields for RNA atoms are placeholders; residue-specific charges were read
from the actual RTP templates.

The [original paper](https://pmc.ncbi.nlm.nih.gov/articles/PMC8047801/) specifies
TIP3P and separate corrections for chloride and RNA. Its DMP calibration changes
the phosphate donor's charge, type and angle parameters to the RNA model;
Table 1 limits the scaling factors to those named chloride/RNA combinations.
That evidence does not establish transfer to the retained Meagher GDP(−3)
diphosphate model or this mixed protein/cofactor/covalent system. Although the
paper discusses broader biomolecular uses, the complete KRAS site was not a
reported validation case.

## Redistribution evidence

No license-named file exists in the complete pinned repository tree, and no
explicit repository-wide grant was found in the retrieved Mg files/README.
Not every third-party file header was audited. Separately, the full article's
license explicitly links to **CC BY 4.0**. Its verified XML is preserved from
[Europe PMC](https://www.ebi.ac.uk/europepmc/webservices/rest/PMC8047801/fullTextXML).
This supplies a possible route to an independently written implementation of
the licensed article's equations/parameters with attribution and recorded
changes. It does not automatically license copied repository files. Published
table precision and the exact implementation values would need reconciliation.

## Outcome and next alternative

No complete-system substitution is admitted. The next bounded source assessment
can consider the [Panteva–Giambașu–York modified 12–6–4 model](https://theory.rutgers.edu/resources/pdfs/Panteva_JChemPhysB_2015_v119_p15460.pdf),
which calibrates specific phosphate/purine interactions against binding data.
It uses TIP4P-Ew and defined nucleic-acid models, so it also requires a coherent
water/cofactor choice and explicit C4 transport; it is not an immediate fix.
Neither that model nor new simulations were started here. The advanced chemistry
track remains independent of the v1 loop-repair release.

## Evidence

- Repository machine-readable inventory: `micromg-applicability.json`.
- Reproducible reader: `inventory_micromg.py`.
- Cache: `/Users/ashujo/.cache/dynamol-research/complete-complex-track/micromg-applicability-v1`.
- `compatibility-inventory.json` contains complete input SHA256 provenance;
  `source-files.json` contains the pinned URLs and Git blob identities.
- The prior `6oim-coordination-v1` through `v4` artifacts and frozen 15+15 panels
  remain unchanged. The unsuccessful institutional PDF retrieval was recorded;
  the same article was subsequently retrieved as full XML from Europe PMC.
