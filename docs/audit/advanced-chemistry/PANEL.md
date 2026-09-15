# DynaMol advanced-chemistry source panels

**Thirty actual deposited structures are frozen: 15 covalent-adduct complexes and 15 metal-containing proteins. No full-case preparation outcomes were used to select them.** These are deliberate chemistry challenges, not a random sample or an estimated PDB-wide success rate.

The JSON manifests preserve source SHA-256 hashes, deposited biological assembly 1 and its operators, author and label atom identifiers, chemical-component dictionaries, original covalent links, metal donors and coordinating waters, loop inventories, and every retained nonstandard residue. Model charges, redox/spin states and parameters remain unassigned at this curation stage.

- [Covalent manifest](covalent-panel.json) · [Metal manifest](metal-panel.json) · [Frozen hashes](panel-freeze.json)
- [Original-source verification and outside-assembly links](panel-source-verification.json) · [Excluded candidates](excluded-candidates.json)

## Covalent panel

Each entry links to its primary PDB record; the manifests also record the associated structure-publication DOI. “No gaps” refers to internal sequence-supported coordinate gaps in the selected deposited chains, not missing terminal tails or side-chain atoms.

| PDB | Complex | Resolution (Å) | Deposited product class | Internal gaps |
| --- | --- | ---: | --- | --- |
| [6OIM](https://www.rcsb.org/structure/6OIM) | KRAS–sotorasib | 1.65 | Cys thioether; acrylamide Michael product | 3 |
| [6UT0](https://www.rcsb.org/structure/6UT0) | KRAS–adagrasib | 1.94 | Cys thioether; fluoroacrylamide Michael product | None |
| [5P9J](https://www.rcsb.org/structure/5P9J) | BTK–ibrutinib | 1.08 | Cys thioether; acrylamide Michael product | 3 |
| [6J6M](https://www.rcsb.org/structure/6J6M) | BTK–zanubrutinib | 1.25 | Cys thioether; acrylamide Michael product | None |
| [8A4V](https://www.rcsb.org/structure/8A4V) | Cathepsin L–E64 | 1.65 | Cys thioether; epoxide-opening product | None |
| [5TDI](https://www.rcsb.org/structure/5TDI) | Cathepsin K–odanacatib | 1.4 | Cys thioimidate; nitrile-addition product | None |
| [6QBS](https://www.rcsb.org/structure/6QBS) | Cathepsin K–alkyne inhibitor | 1.7 | Cys vinyl thioether; alkyne-addition product | None |
| [3LJ7](https://www.rcsb.org/structure/3LJ7) | Humanized-rat FAAH–URB597 | 2.3 | Ser carbamoyl adduct | None |
| [3LJ6](https://www.rcsb.org/structure/3LJ6) | Humanized-rat FAAH–PF-3845 | 2.42 | Ser carbamoyl adduct from urea inhibitor | None |
| [6AX1](https://www.rcsb.org/structure/6AX1) | Monoacylglycerol lipase–carbamate | 2.26 | Ser carbamoyl adduct | None |
| [2H5I](https://www.rcsb.org/structure/2H5I) | Caspase-3–Ac-DEVD-CHO | 1.69 | Cys thiohemiacetal; peptide-aldehyde product | None |
| [3BJM](https://www.rcsb.org/structure/3BJM) | DPP4–saxagliptin | 2.35 | Ser imidate; nitrile-addition product | None |
| [5O4A](https://www.rcsb.org/structure/5O4A) | FGFR1–m-FSBA | 2.01 | Lys sulfonamide; sulfonyl-fluoride substitution | 11 |
| [5LF3](https://www.rcsb.org/structure/5LF3) | Human 20S proteasome–bortezomib | 2.1 | N-terminal Thr boronate adduct | 6 |
| [5F19](https://www.rcsb.org/structure/5F19) | COX2–aspirin acetylation | 2.04 | Ser O-acetyl ester encoded as one modified residue | None |

The product graph is essential: 6AX1 omits the CCD-declared carbamate leaving group; 5O4A omits the CCD-declared fluoride leaving atom; 8A4V uses an already opened E64 product; and 5TDI uses the imine-containing bound ligand. Reconstructing the free inhibitor in any of these cases would change the deposited chemistry. 2H5I represents the inhibitor as a peptide chain, and 5F19 encodes the acetylated serine as OAS with an internal OG–C1A bond.

5LF3 and 5F19 are combined-chemistry challenges. The full proteasome retains its other modified residues, salts and PEG; the COX2 dimer retains glycans and cobalt protoporphyrin. A successful isolated covalent-fragment prototype does not establish that either complete assembly is prepared or validated.

## Metal panel

Counts below are metal atoms in the selected **deposited coordinates**, before biological-assembly operator expansion. 1LUV has two deposited Mn sites and four in its declared tetramer. Contacts come from `struct_conn` and explicit intra-component CCD metal bonds; no donor is assigned solely by distance.

| PDB | Protein/site | Resolution (Å) | Deposited metals | Source-declared metal–water contacts |
| --- | --- | ---: | --- | ---: |
| [1MNC](https://www.rcsb.org/structure/1MNC) | MMP8–hydroxamate | 2.1 | CA × 1, ZN × 2 | 0 |
| [1HFC](https://www.rcsb.org/structure/1HFC) | MMP1–hydroxamate | 1.5 | CA × 1, ZN × 2 | 0 |
| [1GKC](https://www.rcsb.org/structure/1GKC) | MMP9–hydroxamate | 2.3 | CA × 5, ZN × 2 | 6 |
| [1XUC](https://www.rcsb.org/structure/1XUC) | MMP13–non-zinc-binding inhibitor | 1.7 | CA × 2, ZN × 2 | 3 |
| [1BKC](https://www.rcsb.org/structure/1BKC) | ADAM17–hydroxamate | 2 | ZN × 1 | 0 |
| [1CA2](https://www.rcsb.org/structure/1CA2) | Carbonic anhydrase II | 2 | ZN × 1 | 1 |
| [1T64](https://www.rcsb.org/structure/1T64) | HDAC8–trichostatin A | 1.9 | CA × 2, NA × 2, ZN × 1 | 13 |
| [1CLL](https://www.rcsb.org/structure/1CLL) | Calmodulin | 1.7 | CA × 4 | 4 |
| [1HET](https://www.rcsb.org/structure/1HET) | Horse liver alcohol dehydrogenase | 1.15 | ZN × 4 | 4 |
| [1A6Q](https://www.rcsb.org/structure/1A6Q) | PPM1A phosphatase | 2 | MN × 2 | 7 |
| [1HL5](https://www.rcsb.org/structure/1HL5) | Cu/Zn superoxide dismutase | 1.8 | CU × 2, ZN × 2 | 2 |
| [1HCK](https://www.rcsb.org/structure/1HCK) | CDK2–ATP/Mg | 1.9 | MG × 1 | 1 |
| [1A6M](https://www.rcsb.org/structure/1A6M) | Oxy-myoglobin | 1 | FE × 1 | 0 |
| [1LUV](https://www.rcsb.org/structure/1LUV) | Human MnSOD H30V | 1.85 | MN × 2 | 2 |
| [1OKL](https://www.rcsb.org/structure/1OKL) | Carbonic anhydrase II–dansylamide | 2.1 | HG × 1, ZN × 1 | 1 |

The formal ionic charges in CCD are recorded independently of the model choice. In particular, a Zn-bound solvent oxygen does not identify water versus hydroxide; neutral CCD hydroxamic acids do not settle their bound protonation state; Cu/Mn CCD labels do not establish a redox or spin model; and the Fe in HEM and Co in COH have no declared CCD atomic formal charge. 1A6M therefore requires a specific oxyheme model. 1HET also needs review of the publication’s NADH/hydroxide-adduct state. 1OKL retains its Hg site as a separate obligation alongside Zn.

## Scope and evidence limits

- All selected protein sources are human or mammalian. Expression-host names are recorded separately in the original source and are not the protein organism.
- All original CIF files remain unchanged. Assembly choices exclude other independent copies explicitly; molecular components within the chosen assembly are retained.
- Twenty-four of the thirty selected deposited assemblies have no internal sequence gaps. The others have gaps of 3, 3, 11, 6, 4 and 4 residues. Repair is still subject to the existing geometric safeguards.
- Nonstandard heavy-atom inventories match CCD after accounting for explicitly flagged leaving atoms, except for one incomplete crystallization ethanol in 1CLL (missing O). That omission remains visible and must be resolved explicitly rather than deleting the ethanol to obtain a pass.
- Source connection annotations can be incomplete or involve a crystal image. The supplementary verification preserves links outside the selected assembly; no claim of complete coordination is made from the listed contact count.
- No native protein preparation, simulation, coordinate mutation, charge assignment, or force-field validation was performed by the curation scripts. Separate implementation/prototype work is not a full-panel result.
- Source quality, chemical graph validity, complete-complex parameterization, energy/force parity and stability are distinct validation stages. A deposited structure and a successful parameter-generation command are not sufficient to claim a physically reliable model.

## Reproducibility

`inventory_sources.py` downloads only official RCSB CIF/CCD files and generates the raw inventories. `freeze_panels.py` records the predeclared choices and refuses to overwrite a frozen manifest. `panel-freeze.json` locks both manifests plus the source and CCD inputs. `panel-source-verification.json` verifies every locked file and every explicit target-link endpoint against the original atom table. Future changes to the panel must be versioned and explained, preserving this initial selection.
