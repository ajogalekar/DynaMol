# Covalent protein–ligand preparation: method review and native spike

Status: isolated prototype, 2026-09-13. This is not yet enabled in DynaMol.
The scope is fixed, explicitly specified covalent products; classical MD does
not model reaction formation, cleavage, reversibility, or reaction rates.

Research update: [the literature reassessment](covalent-literature-reassessment-v1/REVIEW.md)
supersedes the order of future calculations below. Existing results and scientific
limitations remain valid; the failed full-adduct optimization is not a required
prerequisite for all candidate comparisons.

## Evidence and recommended direction

The practical open-source route is to parameterize the **connected reacted
amino acid–ligand adduct**, using a capped molecular model for charge fitting
and local bonded parameters, then insert that residue into the protein with
an explicit, checked peptide interface. Parameterizing the free ligand and
adding a bond afterwards leaves charges, atom types, angles, torsions,
impropers, and nonbonded exclusions inconsistent with the product.

Amber's modified-residue tutorial supplies the relevant tools: antechamber,
prepgen, parmchk2, and LEaP. It explicitly distinguishes polymer insertion
from ordinary small-molecule setup, identifies the atoms removed at the
boundary and the neighboring atom types, and recommends validating inferred
parameters. Its AM1-BCC example is an approximation, not a general accuracy
guarantee. [Official Amber tutorial](https://ambermd.org/tutorials/basic/tutorial5/index.php)

For the user's accuracy-first extension, prefer **constrained RESP fitted to
HF/6-31G\* electrostatic potentials of the capped connected adduct**, keeping
the canonical ff14SB cap and backbone charges fixed. Use multiple conformers
where feasible and document the geometry-generation method. Published
nonstandard-residue work demonstrates cap/backbone constraints and the need
to distinguish charge-model choices; AM1-BCC does not automatically reproduce
ff14SB residue charges. [Charge-method comparison](https://pmc.ncbi.nlm.nih.gov/articles/PMC6641686/),
[Covalent conjugate parameterization](https://pmc.ncbi.nlm.nih.gov/articles/PMC6934455/),
[Constrained capped-residue RESP example](https://www.rsc.org/suppdata/d3/sc/d3sc05450k/d3sc05450k1.pdf)

Constrained RESP improves electrostatic consistency; it does **not** validate
new bonded terms. High-penalty parmchk2 assignments and the reaction center
still require QM geometry/force or torsional checks. Retain a clearly labeled
provisional AM1-BCC route only if the product state and the corresponding
parameter limitations are explicit.

## Actual human KRAS–sotorasib prototype: 6OIM

The deposited structure is human KRAS G12C covalently bound to AMG 510 at
1.65 Å resolution. [RCSB 6OIM](https://www.rcsb.org/structure/6OIM)

Independent inspection of the downloaded original mmCIF and CCD established:

- `struct_conn/covale1` records author `A:CYS12:SG` to author
  `A:MOV303:C25`, same-image symmetry `1_555`, distance 1.805 Å. Label
  identities differ: protein `A/26`, ligand `D`.
- The MOV CCD is explicitly the **bound form**, with a saturated C24–C25
  single bond. Applying another generic alkene-reduction rule would be wrong.
- The connected product has no S–H at the cysteine sulfur and has two H each
  on MOV C24 and C25. Hydrogens are derived after adding the explicit bond.
- The prototype retains CYS alpha-carbon R and MOV C20 S stereochemistry.
- The complete structure also contains GDP and Mg with recorded coordinating
  waters. These remain required for a later full-complex test; the cap probe
  makes no claim to prepare or validate those components.

[Original mmCIF](https://files.rcsb.org/download/6OIM.cif),
[MOV CCD](https://files.rcsb.org/ligands/download/MOV.cif)

The isolated script `covalent_cap_spike.py` builds a 94-atom, 52-heavy-atom,
neutral ACE–CYS(MOV)–NME model. The affected residue and ligand heavy coordinates
remain exact. Cap coordinates use the adjacent observed peptide backbone as
temporary geometric proxies, and only generated H atoms undergo initial
MMFF relaxation. These caps are never exported to a full protein. This is a
charge/parameter model, not a reconstruction of a new binding pose.

### Results actually obtained

Native AmberTools 24.8, GAFF2 2.2.20, AM1-BCC completed in **176.5 s** on the
local two-thread limit. The output contains the S–C bond and all linked
terms together: **one bond, four angles, nine proper torsion entries** involve
both linked atoms. Native types at the link are `ss` and `c3`.

The resulting Amber topology and the converted XML agree in OpenMM Reference
at the bound pose and two deterministic 0.003 Å perturbations:
maximum energy difference **4.3e-9 kJ/mol**, maximum force-component difference
**3.7e-7 kJ/mol/nm**. This is a successful atom-map and force-field conversion
test, not evidence that the force field predicts the chemistry accurately.

An additional independent **Sander** comparison was run. Raw Sander versus
OpenMM defaults differ by approximately 0.0446 kJ/mol here because Amber's
legacy stored-charge convention uses 18.2223, while OpenMM uses its current
electrostatic constant. A diagnostic copy matched these conventions and
agreed within 8.1e-6 kJ/mol and 0.00325 kJ/mol/nm. The diagnostic reports both
raw and adjusted results; production charges and files were not changed.
The scaling was derived independently from the documented prmtop reader and
a separate two-charge OpenMM probe, not fitted to this adduct.
[OpenMM Amber reader](https://github.com/openmm/openmm/blob/master/wrappers/python/openmm/app/internal/amber_file_parser.py),
[OpenMM physical constants](https://github.com/openmm/openmm/blob/master/platforms/reference/include/SimTKOpenMMRealType.h)

**Unresolved parameter quality:** parmchk2 emitted no missing `ATTN` entries,
but several aromatic torsional analogies have penalty scores **223–774**.
Thus absence of `ATTN` is insufficient to claim accurate parameters. These
assignments need an explicit uncertainty report and targeted QM assessment.
The existing application currently checks missing terms but does not promote
analogy penalties to a preparation-quality assessment.

Artifacts (relative to the repository):

- `build/advanced-chemistry/covalent/capped-adduct.json`: exact input hashes,
  product graph decisions, atom/cap map, hydrogen inventory and native results.
- `build/advanced-chemistry/covalent/capped-native/`: complete native inputs,
  AM1 optimization log, charged mol2, frcmod, prmtop, inpcrd and converted XML.
- `build/advanced-chemistry/covalent/native-link-terms.json`: link-term inventory
  and high-penalty assignments.
- `build/advanced-chemistry/covalent/capped-native/conversion-validation.json`:
  XML versus native-topology checks in the same OpenMM engine.
- `build/advanced-chemistry/covalent/capped-native/sander-parity.json`:
  independently executed Sander checks and electrostatic conventions.

## Protein interface: native mechanics checkpoint and remaining work

The first cap above is all GAFF2. A second, explicitly separate mechanics
probe replaces the 18 canonical cap/backbone atom types with their native
ff14SB types, keeps the connected sidechain–ligand GAFF2 types, and obtains
the mixed terms through native `parmchk2`/LEaP. Every term wholly within the
canonical protein atom-type set comes directly from loaded ff14SB; the
generated GAFF supplement only contains terms touching a GAFF atom type.
Analogy provenance remains visible. No crosslink constants were authored by hand. Its charges remain
the provisional joint AM1-BCC charges, so this probe is not yet a finished
ff14SB-compatible residue model.

The resulting topology is complete. An independent native LEaP reference,
`sequence {ACE CYS NME}`, establishes that **all 17 bonds, 27 angles and 48
dihedral entries wholly within the canonical cap/backbone subset match**.
This isolates one important interface invariant; it does not validate the
modified sidechain, its QM energy surface, or insertion into a larger peptide.
`covalent_interface_spike.py` reproduces the native generation and comparison.
Evidence is in `build/advanced-chemistry/covalent/interface-native-scope-v4/`.

The native parameter exporter additionally records the full 94-atom
inventory, 98 bonds, 173 angles, 305 torsions, and 497 nonbonded exceptions,
including zero exclusions and scaled 1–4 pairs. It reads parameter values
from the independent native prmtop with ParmEd, never from the OpenMM System
being tested. Its topology identities describe the isolated capped model;
the original protein/ligand source map is separately preserved.

The larger peptide assembly path has also been exercised with an explicitly
synthetic charge fixture: remove the 12 exactly neutral cap atoms, declare
central N/C connection points, and build native `ACE–ALA–COV–GLY–NME`.
The resulting 111-atom topology retains the connected 82-atom adduct and has
two peptide boundary bonds, eight boundary angles and 41 boundary torsion
entries. The independent full native parameter manifest matches the actual
OpenMM System, including all 585 exclusions/1–4 exceptions. This verifies
the assembly and parameter plumbing; `usable_parameters` remains `false`,
and no MD is run using the synthetic charges. Actual quantum-derived charges
must be tested through this same assembly path before model approval.
Evidence: `build/advanced-chemistry/covalent/polymer-native-scope-v4/`;
reproducer: `covalent_polymer_interface_spike.py`.

An actual Context/independent-engine test caught and corrected errors in
the initial isolated hybrid prototypes. Protein-database `parmchk2` output may
contain unmarked all-zero torsion placeholders with IDIVF 0 and periodicity
0. Treating every no-`ATTN`/no-penalty line as a valid protein override
replaced legitimate GAFF torsions; a static parameter comparison alone did
not catch it. The same protein output can emit unmarked default impropers
for GAFF-only types, incorrectly replacing valid GAFF carbonyl impropers.
The final native-type-scoped construction never uses these protein-database
overrides: canonical protein terms remain native ff14SB, while the GAFF
supplement only covers terms touching GAFF types. Invalid adduct/mixed records
are rejected, and native Context creation with finite energy/forces is
required. Earlier artifacts are retained as failed/superseded prototypes.
Both final capped and larger
peptide systems pass independent Sander/OpenMM parity at three poses after
the separately documented electrostatic-constant normalization: maximum
energy differences 6.6e-6 and 1.11e-5 kJ/mol; force-component differences
0.00325 and 0.00354 kJ/mol/nm. This is implementation consistency, not
chemical accuracy, especially for the synthetic-charge peptide fixture.

Still required before full polymer insertion:

1. Freeze the reviewed ff14SB backbone types and canonical charges at the
   peptide boundary. Fit the sidechain plus complete ligand jointly under
   the adduct's explicit total charge. Constrain ACE/NME to their canonical
   neutral charges so their removal does not require an arbitrary charge
   redistribution.
2. Derive a native LEaP residue template with explicit N/C connection points.
   Preserve ff14SB peptide/backbone terms and obtain the remaining local terms
   from the complete connected model. Native parmchk2 supports the ff14SB
   supplement and prints analogy provenance; missing or poorly supported
   mixed-type terms require review/QM refinement rather than zero-force or
   generic constants.
3. Test the residue in a larger capped peptide including both adjacent
   residues. Verify every cross-boundary bond, angle, proper/improper torsion
   and 1–4 interaction. Testing the isolated central residue cannot detect
   missing peptide-interface terms.
4. Preserve an authoritative native Amber topology, charge/atom-type source
   map and explicit connectivity. Validate any XML/GROMACS conversion against
   the complete native system using coordinate and torsional perturbations.
5. Only then integrate with solvation, saved preparation reload, measurement
   atom maps and both simulation engines.

OpenMM explicitly represents external residue bonds in templates. However,
DynaMol's `IdentityForceField` presently prohibits external bonds for its
native **noncovalent** ligand maps. Extending this needs a distinct covalent
component identity contract with validated external endpoints; removing that
check globally would weaken existing noncovalent safeguards.
[OpenMM force-field template documentation](https://docs.openmm.org/latest/userguide/application/06_creating_ffs.html)

The prepared-system loader currently accepts a GAFF2/AM1-BCC noncovalent bundle.
The covalent model needs a separate versioned method descriptor, exact bonded
component identities, named protonation/product state, model/cap constraints,
native parameter files, source hashes, and required validation reports. A
native-system factory may be simpler initially than generalizing template
generation; it still needs a demonstrated route through water-box creation
and subsequent atom additions.

## Product-state handling across reaction classes

| Class | Required product evidence | Important failure mode |
|---|---|---|
| Cys Michael thioether | Recorded S–C endpoint; reacted bond orders and H inventory | Reducing an already reacted CCD or retaining S–H |
| Cys nitrile/thioimidate | Carbon/nitrogen bond orders and N protonation | Treating a neutral thioimidate as an unreduced nitrile |
| Epoxide opening | Which C–O bond opened; linkage and OH state | Keeping the epoxide ring closed after attachment |
| Ser carbamate/acyl product | Retained acyl fragment, leaving atoms, catalytic O state | Keeping leaving-group atoms or adding an O–H across the link |
| Lys sulfonyl adduct | Explicit S–N connection, leaving group, amine state | Keeping fluoride/leaving atoms or overvalent nitrogen |
| Boronate/tetrahedral products | Coordination number, formal charge and geometry | Replacing tetrahedral/charged boron with a generic neutral boronic acid |
| Peptide or multiple-residue adduct | Complete connected component and all polymer endpoints | Splitting one covalent component into independent ligands |

These are implementation requirements, not claims that all classes have been
parameterized. The frozen 15-case covalent panel includes varied classes and
declared difficult cases; success must be reported per class and complete
complex, with unsuccessful inputs retained in the denominator.

### Frozen 15-case product-graph pressure test

`covalent_product_graph_audit.py` constructs capped graph candidates for all
target connections in **14 of the 15 cases**, retaining original source/CCD
hashes, all observed product heavy coordinates, explicit endpoint valence/H
counts, observed stereochemistry, cap identities and leaving-atom decisions.
The deposited CCDs already encode the reacted bond orders in these examples,
so the audit makes no heavy bond-order transformations. Only absent,
explicitly CCD-flagged leaving atoms are omitted; observed atoms are never
deleted. All four Cys Michael cases produce neutral capped graph candidates
with no S–H and two H on the attached terminal carbon.

This does **not** establish successful force-field preparation or correct pH
states. For example, the E64 CCD candidate has formal charge +1; relevant
acid/base states still require explicit treatment. Source complex components
outside the capped product are not modeled in this graph-only audit.
`chemical_state_accepted` and `full_preparation_accepted` remain false.

The generic peptide-polymer expansion now handles **2H5I** by retaining its
entire ACE–ASP–GLU–VAL–ASJ entity: 35 ligand heavy atoms, four explicitly
sequence-supported peptide joins, and both deposited target covalent bonds.
Original residue/atom identities survive the temporary connected graph. A
missing whole ligand residue, broken peptide join or additional unresolved
external link is rejected. The capped combined candidate contains 46 heavy
atoms and 88 atoms including hydrogens. It has not been assigned validated
charges or accepted as a physical model.

**5LF3** remains blocked because it requires a reviewed tetrahedral
boronate/N-terminal Thr charge model. Neutral trigonal boron is not reused.
All 15 cases remain in the fixed denominator. Evidence:
`build/advanced-chemistry/covalent/frozen15-graph-v4/summary.json`; 14 focused
graph contract tests pass, including the complete peptide expansion.

### Native atom-type and parameter-coverage screen

The 14 distinct constructed graph chemistries all complete native GAFF2 atom
typing and parameter lookup without `ATTN` entries or invalid zero-periodicity
placeholders. The screen uses atom-typing-only `antechamber -j 1`, preserves
the input explicit bond orders, and verifies atom/index identity and native
coordinate serialization. It does **not** fit charges: MOL2 charge columns
are explicitly unusable zero placeholders. No system or dynamics is launched.

This is coverage evidence, not accuracy. Native analogy penalty maxima range
from 6.0 to 780.0 across these candidates; every assigned term, score and
explicit default is retained. The largest scores occur in remote aromatic
chemistry as well as some other environments, so a complete term table is
insufficient to accept a model. The E64 candidate still carries the unreviewed
CCD +1 state. Protein ff14SB boundary precedence, parent QM charge fitting,
attachment profiles and full-complex checks remain separate requirements.
Evidence: `build/advanced-chemistry/covalent/frozen15-native-coverage-v2/summary.json`.
The first audit output is retained as a software-format failure: its initial
four-decimal tolerance overlooked the native AC intermediate's three-decimal
coordinates. The corrected tolerance follows the actual file format; original
SDF coordinates remain untouched, while native typed files round some heavy
coordinates by at most 0.0005 Å per component. Eighteen combined
graph and coverage contract tests pass, including deliberately changed
coordinates, bond orders and unmarked zero-periodicity rejection.

## Acceptance before user-facing support

- Product graph, original heavy inventory, stereochemistry, explicit charge,
  proton count and atom mapping validated independently.
- Native topology complete with no missing terms; analogy scores and source
  assignments exposed, with targeted QM tests for unsupported local chemistry.
- RESP cap/backbone constraints, total charge, convergence, fit residuals and
  held-out conformer behavior recorded. A converged fit is not automatically
  a good fit.
- Native/XML/engine parity covers crosslinks and peptide boundaries, including
  exclusions, 1–4 scaling and improper order.
- Geometry/finite-force checks, restrained equilibration and brief stability
  runs preserve the intended product and stereo configuration.
- Saved/reloaded systems and both engines reproduce the same atom identities
  and intended model. Persistent bond chemistry does not imply that a
  reactive mechanism or experimental binding affinity has been validated.

No central DynaMol preparation guards were changed by this prototype.

## RESP execution and local runtime evidence

The isolated RESP fitter fixes canonical ACE/NME charges and the central
N/H/CA/HA/C/O backbone charges, jointly fitting the connected sidechain and
ligand to the declared total charge. All 18 boundary/cap atoms are fixed;
cap removal therefore does not redistribute arbitrary charge. Stage 2 only
refines the four noncanonical methyl groups with equivalent methyl H charges;
potentially diastereotopic methylene H atoms are not forced equal. This
conservative two-stage policy is explicit and still requires charge-model
validation; it is not claimed as a universal residue fitting prescription.

A deliberately synthetic 4,000-point ESP fixture verifies native file format,
atom order and exact charge constraints. The original native RESP and the
local heap-array rebuild give **identical printed charges for all 94 atoms in
both stages**, with zero frozen-boundary error. This synthetic fixture is not
quantum evidence and its charges must not be used as production parameters.
`resp-synthetic-heap-v2/parity.json` preserves the comparison and executable
hash; the selected runtime is audit-only and has not replaced the installed
or distributed application binary.

The first actual capped-adduct RHF/6-31G* optimization was deliberately stopped
after 18 minutes and three completed initial SCF cycles, with its progress,
checkpoint and intentional-stop reason retained. Conventional SCF is
expensive on this laptop; the installed PySCF wheel also reports no OpenMP
support, so the requested two-thread setting does not parallelize its OpenMP
code. A separately labeled, fixed-initial-geometry **DF-RHF/6-31G*** comparison
uses the explicit `def2-universal-jkfit` auxiliary basis for Coulomb and
exchange. It completed energy and gradient in 430.4 seconds with 20 SCF
cycles at the unchanged strict tolerances. This is an additional
approximation; no optimization or fitting results from it are accepted as
equivalent without comparing converged energies/gradients to the conventional
calculation. The conventional comparison starts from the converged DF density
at identical geometry/basis/state, which changes only its initial guess.
The comparison helper verifies checkpoint provenance, orbital dimensions,
density electron count and identity; a deliberately changed geometry is
rejected. PySCF documents both
auxiliary-basis selection and the distinction between J and JK fitting.
[PySCF density-fitting documentation](https://pyscf.org/user/df.html)

The conventional same-geometry reference has now actually completed in
3,808.35 seconds. An independent reread verifies both accepted array hashes,
identical atom identities and bitwise-identical initial coordinates. DF minus
conventional energy is +0.0013747982 Hartree (+0.862699 kcal/mol); maximum and
RMS gradient-component differences are 5.53992×10⁻⁵ and 1.45855×10⁻⁵
Hartree/Bohr, respectively (0.065694 and 0.017296 kcal/mol/Å). The largest
component difference is at the temporary ACE oxygen. The absolute energy
offset at this one geometry is not a conformational-energy error, and this
comparison does not validate the entire DF potential-energy surface.
Evidence: `build/advanced-chemistry/covalent/df-rhf-initial-comparison-independent.json`.

A separate private PySCF runtime with actual OpenMP support was checked on
the same full 94-atom DF calculation: identical coordinates, energy difference
−2.18×10⁻¹¹ Hartree and maximum gradient difference 1.90×10⁻¹² Hartree/Bohr.
Its concurrent wall time (411 versus 430 seconds) does not establish a material
speedup. The full-adduct DF geometry candidate now runs there with two threads,
8 GB memory, five frozen temporary cap heavy atoms, a 100-step limit and a
six-hour deadline. This provisional calculation overlaps the exact reference;
it is not accepted as an accurate geometry before the reference comparison,
stereochemical/geometry checks and later attachment-profile tests. Live process
IDs, progress, results and required next steps are in `LIVE-JOBS.json`.

The durable `continue_covalent_qm.py` controller waits for that existing
optimization and a numerically completed conventional reference. Before any
new ESP calculation it verifies source/result hashes, exact atom order and
state, the five frozen temporary cap heavy atoms, recorded product stereo
centers, and gross bonded distances. It then requests **conventional RHF**
ESP at the accepted DF geometry, using the final DF checkpoint only as the
initial density, and runs the constrained native RESP stages. The DF geometry
approximation remains explicit. The ESP stage has a four-hour native wall
limit and an independent subprocess deadline; each native RESP stage has a
120-second limit. Existing started or completed stages are not overwritten
or duplicated, and a shared launch lock plus active-process inspection limits
concurrent quantum jobs to three. Twenty-four combined graph, coverage and
continuation contract tests pass; the controller itself is currently waiting
for real optimization completion. Its successful terminal status will be
`research_candidate`, with full-adduct torsion and complete-complex validation
still required.

The resident controller is now PID 71825, with a 15 GiB free-disk gate before
its future conventional ESP calculation. Its exact loaded source is archived
under its recorded SHA-256 in the continuation's `source-checkpoints/` folder.
The source on disk additionally guards against a second continuation of the
same parent; a real negative launch check rejected the resident process before
creating another workspace. A briefly overlapping *waiting* fallback controller
was stopped before any child/QM launch, with its evidence retained. The active
parent optimizer was never restarted or duplicated.

The six full-parent RHF profile jobs have not been launched and are now held
for the separately documented [reference-method feasibility check](covalent-reference-method-feasibility.md).
It proposes a small full-parent dispersion-corrected DFT / GFN2-xTB comparison
before committing to long scans, with explicit limits on applicability and
disk use. The current RHF/RESP charge calculation remains unchanged.
