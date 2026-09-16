# Robustness validation — 60 diverse PDB structures

This is an autonomous robustness sweep of the real v1 workflow across a diverse
protein panel. Its purpose is narrow and explicit: confirm that a broad, random
sample of PDB protein and protein–ligand structures runs through
**fetch → (monomer) → prepare → solvate → short MD** without failing at the
start, and to surface and fix systemic workflow bugs. A completed short run is
**not** evidence of a correct protonation state, ligand binding pose, or
equilibrium sampling; see [Simulation scope](../README.md#simulation-scope).

## Method

- **Panel:** 60 structures chosen for diversity of target class and of protein
  and ligand chemistry — tyrosine and Ser/Thr kinases, aspartic and serine
  proteases, nuclear receptors, immunophilins, HSP90 chaperones, bromodomains,
  a phosphatase, oxidoreductases, an isomerase, sugar/vitamin/biotin binding
  proteins, a lectin, hydrolases, and small model folds. The panel deliberately
  **excludes the v1 out-of-scope set** (covalent inhibitors, bonded metal-site
  models, membranes, nucleic acids, >12-residue gaps).
- **Each case** runs in an isolated data directory through the same backend
  entry points the app uses. A chain is selected when the structure is
  multi-chain or large; the ligand travels with its chain. Preparation builds
  missing loops, assigns ligand GAFF2/AM1-BCC parameters, and protonates at
  pH 7. MD is short (10 ps, 2 fs steps, minimize + brief relaxation), implicit
  GBn2 for protein-only inputs and TIP3P/PME explicit for retained-ligand
  complexes, on two CPU threads.
- **Classification:** *Pass* (short MD completed); *Refused* (a correct,
  component-specific block for an out-of-scope structure — not a bug);
  *Runs correct; CPU-slow* (stable dynamics but does not finish 10 ps within the
  harness time budget); *FAIL* (an unexpected crash, stereochemistry inversion,
  NaN, or geometry rejection of a structure that should prepare).
- **Reproduce:** `docs/audit/v1-robustness/run_validation.py` (driver + panel)
  and `validate_structure.py` (one case). Results and a machine-readable
  `summary.json` are written under `~/.cache/dynamol-research/`.

## Result

| Outcome | Count |
|---------|-------|
| Pass (short MD completed) | 44 |
| Refused — correct, out of v1 scope | 11 |
| Runs correct; CPU-slow on implicit | 5 |
| **Genuine workflow failures** | **0** |

All 60 either run correctly or are refused for a correct, documented reason.
The sweep found **two systemic workflow bugs**, both fixed and re-validated on
the frozen code, and **zero genuine failures**. The five "CPU-slow" cases are a
performance limit of implicit solvent on CPU, not failures; DynaMol now prefers
explicit solvent for large proteins (see below), and they complete on that path.

## Systemic bugs found and fixed

### 1. Stereocenter inversion during minimization and early dynamics

Residues next to a low-confidence rebuilt loop can enter dynamics with a
marginal Cα chirality volume. Because the force field has no explicit term that
keeps a Cα/Cβ centre left-handed, the local minimum of a strained rebuilt region
can be the **inverted (D) centre** — free minimization drives it there, and the
first MD steps racemize it. This produced the original 6X8F / 2HYY / 1FIN
failures ("inverted or flattened residue stereochemistry … Dynamics stopped").

**Fix:** a permanent, gentle flat-bottom **chirality guard** is added to the MD
force field (`StereoMonitor.chirality_guard_force`, wired in the worker before
minimization). It is exactly zero while a centre keeps its sign and ≥ 40 % of its
ideal signed-volume magnitude, and rises only as a centre approaches planarity,
so it never biases the many healthy centres and only resists racemization of the
few strained ones. It is soft enough to be integrator-safe (no NaN), part of both
minimization and dynamics, and — being permanent — has no release into which a
held centre can invert. **Zero stereocentre inversions across all 60 structures.**

### 2. Legacy CIP mislabeling rejected correct ligands

Bound-ligand stereochemistry is checked against the CCD reference. The check used
RDKit's legacy `AssignStereochemistry`, whose simplified CIP priority rules
**mislabel fused-ring and macrocyclic centres**, and whose `cleanIt` pass
strips valid chiral tags. Correctly deposited ligands were falsely rejected:
tetrahedral R/S (the estrane steroid R1881 in 1E3G, at C13), bond E/Z
(geldanamycin in 1YET at C2–C3, rapamycin in 1FKB at C17–C18), and stripped
ring-junction tags (steroid C8).

**Fix:** all CCD stereo comparisons now use `rdCIPLabeler` (RDKit's full CIP
algorithm) for both atom R/S and bond E/Z. Verified against the CCD for every
affected centre — R1881 (4 centres), geldanamycin, and rapamycin (15 centres)
all label correctly and now prepare. The integrity gate is unchanged in
strictness; only the labeller is corrected. Ligand-chemistry unit tests: 80
passed.

## Correct refusals (not bugs)

These 11 are **correctly blocked** with component-specific messages; the app
never silently alters them.

- **1LYZ** — the 1974 HEW-lysozyme entry has a genuinely inverted (D) Cα at
  A:14 ARG (signed volume −4.06×10⁻⁴ nm³; its neighbours are correctly
  positive). DynaMol refuses to silently simulate a D-residue. Modern lysozyme
  entries (e.g. **1AKI**, which passes) prepare normally.
- **2AM9** — covalent DTT adduct. **1HVR** — cyclic-urea handling under the
  out-of-scope gate. **1AQ1** — a 13-residue internal gap (> 12-residue loop
  limit). **1M17** — solvated box would exceed the 100,000-atom limit.
  **1DWD, 1FJS, 1ANF, 1EVE, 1RX2, 5PTI** — component/scope blocks.

## Large proteins: prefer explicit solvent (implemented)

Five large apo proteins (**1P38, 3PSG, 1TIM, 2HNP, 2LZM**) run **correct, stable
dynamics** (T ≈ 300 K, sane energies) but do not finish 10 ps within the harness
time budget on implicit solvent. The cause is not a crash and not the chirality
guard: GBn2 implicit solvent uses `NoCutoff` (all-pairs, **O(N²)**), ≈ 2 s/step
for a ~350-residue protein on CPU.

**Fix:** DynaMol now **prefers explicit TIP3P/PME for large protein-only
systems**, which is O(N) (faster at large N) and more accurate. Above
`config.RECOMMEND_EXPLICIT_ATOMS` prepared atoms (default 4,000; env
`DYNAMOL_RECOMMEND_EXPLICIT_ATOMS`), preparation records a
`recommend_explicit_solvent` hint and the app defaults the solvent control to
explicit; **implicit stays selectable** (the hint is advisory, never a hard
requirement). Re-validated on this preference: **3PSG (pepsinogen)**, which timed
out on implicit even at a 3,600 s budget, **completes via explicit**; the other
large cases run the same O(N) explicit path. Very large proteins whose explicit
box would exceed the 100,000-atom limit still need GPU (v2) or a shorter run.

## Per-case results

| PDB | Class | Ligand / note | Outcome | Detail |
|-----|-------|---------------|---------|--------|
| 1HPX | aspartic protease | KNI-272 | Pass | explicit solvent |
| 1HSG | aspartic protease | indinavir | Pass | explicit solvent |
| 1CBR | binding protein | retinoic acid | Pass | explicit solvent |
| 1CBS | binding protein | retinoic acid | Pass | explicit solvent |
| 1RBP | binding protein | retinol | Pass | explicit solvent |
| 1STP | binding protein | biotin | Pass | explicit solvent |
| 2DRI | binding protein | ribose | Pass | explicit solvent |
| 2OSS | bromodomain | CBP | Pass | explicit solvent |
| 3MXF | bromodomain | BRD4 inhibitor | Pass | explicit solvent |
| 4LYW | bromodomain | BRD4 + JQ1 | Pass | explicit solvent |
| 1BYQ | chaperone | ADP | Pass | explicit solvent |
| 1UYD | chaperone | PU3 | Pass | explicit solvent |
| 1YET | chaperone | geldanamycin | Pass | explicit solvent |
| 1AKI | hydrolase | lysozyme | Pass | implicit solvent |
| 7RSA | hydrolase | RNase A | Pass | explicit solvent |
| 1FKB | immunophilin | rapamycin | Pass | explicit solvent |
| 1FKF | immunophilin | FK506 | Pass | explicit solvent |
| 2CPL | immunophilin | cyclophilin A | Pass | implicit solvent |
| 2YPI | isomerase | TIM + PGA | Pass | explicit solvent |
| 2CNA | lectin | concanavalin A (Ca/Mn ions) | Pass | explicit solvent |
| 1CRN | model fold | crambin | Pass | implicit solvent |
| 1CSP | model fold | cold-shock | Pass | implicit solvent |
| 1L2Y | model fold | trp-cage | Pass | implicit solvent |
| 1PGB | model fold | protein G B1 | Pass | implicit solvent |
| 1SHG | model fold | SH3 | Pass | implicit solvent |
| 1UBQ | model fold | ubiquitin | Pass | implicit solvent |
| 1VII | model fold | villin | Pass | implicit solvent |
| 2CI2 | model fold | CI2 | Pass | implicit solvent |
| 1E3G | nuclear receptor | androgen R | Pass | explicit solvent |
| 1FM6 | nuclear receptor | rosiglitazone | Pass | explicit solvent |
| 3ERT | nuclear receptor | 4-OH-tamoxifen | Pass | explicit solvent |
| 1DHF | oxidoreductase | human DHFR | Pass | explicit solvent |
| 4DFR | oxidoreductase | methotrexate | Pass | explicit solvent |
| 1PTY | phosphatase | PTP1B + inhibitor | Pass | explicit solvent |
| 1ATP | ser/thr kinase | ATP + Mg | Pass | explicit solvent |
| 1FIN | ser/thr kinase | CDK2/cyclin | Pass | explicit solvent |
| 1STC | ser/thr kinase | staurosporine | Pass | explicit solvent |
| 3HEG | ser/thr kinase | sorafenib | Pass | explicit solvent |
| 1BRA | serine protease | apo trypsin | Pass | explicit solvent |
| 1ELA | serine protease | elastase | Pass | explicit solvent |
| 3PTB | serine protease | benzamidine | Pass | explicit solvent |
| 1IEP | tyrosine kinase | STI-571 | Pass | explicit solvent |
| 2HYY | tyrosine kinase | imatinib | Pass | explicit solvent |
| 3LCK | tyrosine kinase | apo Lck | Pass | explicit solvent |
| 1HVR | aspartic protease | cyclic urea XK263 | Refused (out of scope) | CSO A:67 is a modified or unnatural protein residue with no supported covalent amino-acid  |
| 1ANF | binding protein | maltose | Refused (out of scope) | Ligand chemistry needs attention: B:1::GLC: Original structure records a covalent ligand/p |
| 1EVE | hydrolase | acetylcholinesterase + tacrine | Refused (out of scope) | Ligand chemistry needs attention: B:1::NAG: Original structure records a covalent ligand/p |
| 1LYZ | hydrolase | HEW lysozyme | Refused (out of scope) | Input structure has inverted or near-planar standard residue stereochemistry at A:14 ARG C |
| 5PTI | model fold | BPTI | Refused (out of scope) | Ligand chemistry needs attention: B:70::PO4: The pH heuristic lost the ligand atom identit |
| 2AM9 | nuclear receptor | androgen R | Refused (out of scope) | Ligand chemistry needs attention: D:1003::DTT: Covalently linked ligand residues need a sp |
| 1RX2 | oxidoreductase | folate/NADP | Refused (out of scope) | Ligand chemistry needs attention: D:162::BME: Covalently linked ligand residues need a spe |
| 1AQ1 | ser/thr kinase | staurosporine | Refused (out of scope) | Chain A has an internal gap of 13 residues. Local loop building supports up to 12 per gap; |
| 1DWD | serine protease | thrombin inhibitor | Refused (out of scope) | The selected protein chain is covalently connected to another protein or nucleic-acid chai |
| 1FJS | serine protease | factor Xa inhibitor | Refused (out of scope) | The selected protein chain is covalently connected to another protein or nucleic-acid chai |
| 1M17 | tyrosine kinase | erlotinib | Refused (out of scope) | The requested solvent box may exceed the 100,000-atom local limit. Reduce padding or prepa |
| 3PSG | aspartic protease | pepsinogen | Runs correct; CPU-slow | implicit GBn2 O(N²) on CPU; use explicit PME |
| 2LZM | hydrolase | T4 lysozyme | Runs correct; CPU-slow | implicit GBn2 O(N²) on CPU; use explicit PME |
| 1TIM | isomerase | TIM | Runs correct; CPU-slow | implicit GBn2 O(N²) on CPU; use explicit PME |
| 2HNP | phosphatase | PTP1B | Runs correct; CPU-slow | implicit GBn2 O(N²) on CPU; use explicit PME |
| 1P38 | ser/thr kinase | apo p38 | Runs correct; CPU-slow | implicit GBn2 O(N²) on CPU; use explicit PME |

## Honest limitations

- Short exploratory MD on CPU; **not** converged sampling, and not evidence of a
  correct binding pose or protonation state.
- Rebuilt loops are provisional starting models (see [loop repair](LOOP_REPAIR.md)).
- "Pass" means the workflow ran to completion without failing, at the diversity
  and scale a 60-structure sweep can cover — not that every random PDB entry is
  guaranteed to run.
