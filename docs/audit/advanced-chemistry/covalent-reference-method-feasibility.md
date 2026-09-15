# A practical independent reference for attachment flexibility

The active RHF/RESP charge workflow continues unchanged. The six proposed
full-parent RHF torsion optimizations are **on hold** until a more useful
reference strategy is qualified. HF-only agreement would leave correlation
and dispersion untested in the large, flexible sotorasib adduct.

## What the external evidence supports

The BespokeFit paper compared several reference-data strategies on its TYK2
ligand set. DFT single-point energies on GFN2-xTB scan geometries reduced the
reported energy RMSE from about 1.10 to 0.60 kcal/mol; full
B3LYP-D3BJ/DZVP scans performed better, at about 0.29 kcal/mol. This supports
testing a cheaper geometry-plus-energy strategy. It does **not** establish
its accuracy for our 94-atom capped covalent adduct, for GAFF2/ff14SB, or for
the actual protein-bound conformational ensemble.
[BespokeFit paper, Table 3](https://pmc.ncbi.nlm.nih.gov/articles/PMC9709916/).

OpenFF recommends its full DFT reference by default and treats faster
semiempirical methods as needing validation for the particular chemistry.
Its SMIRNOFF torsion parameters also cannot simply replace Amber terms.
[OpenFF reference-method guidance](https://docs.openforcefield.org/projects/bespokefit/en/latest/getting-started/faq.html),
[BespokeFit scope](https://docs.openforcefield.org/projects/bespokefit/en/latest/).

## Proposed bounded feasibility experiment

1. After the current parent optimization and actual RESP charges finish,
   qualify a dispersion-corrected DFT calculation on **one unchanged complete
   parent geometry**. Record wall time, memory, scratch usage, SCF convergence,
   numerical-grid settings, expanded orbital/auxiliary bases, XC definition,
   and the separate electronic and D3 contributions. Set a four-hour outer
   budget for this first point; do not start six long optimizations first.
2. Generate only three complete-parent geometry candidates for the primary
   S–C attachment: the baseline and the two predeclared ±30° offsets. A
   separately qualified GFN2-xTB/geomeTRIC adapter could make this inexpensive.
   Keep the same five Cartesian cap freezes, actual reaction product, all 94
   atoms, electronic state and recorded stereochemistry. Do not replace the
   ligand with a truncated fragment or use GFN2 charges in the Amber model.
   Give each exploratory geometry a 15-minute limit; measured performance
   determines whether further work is practical.
3. Evaluate the same dispersion-corrected DFT single-point method at those
   geometries. Then perform **one** matching constrained DFT optimization at
   a predeclared displaced point, within a six-hour limit. Comparing its
   geometry and relative energy with the inexpensive geometry at the same
   target tests whether the mixed method is usable here. Record every point,
   including nonconvergence, stereo changes and constraint failures.
4. If this limited comparison is informative and affordable, extend the
   fixed six-target plan sequentially. Use actual joint RESP charges and
   verified native protein-interface terms for the classical comparison.
   Compare MM and QM relative energies at identical geometries before any
   optional MM relaxation. Investigate disagreements rather than tuning on
   the same points and presenting the fit as independent validation.

This is a staged proposal, not a result or an automatically accepted error
threshold. The three-point experiment cannot establish an entire torsion
surface, rotational barriers, reaction kinetics or protein-bound populations.
Coupled motion and unusual intra-adduct hydrogen bonds may require additional
held-out conformers. A gas-phase reference is also not a validation of the
full solvated complex.

## Implementation details that must be explicit

To reproduce the paper's nominal method, resolve the actual **DZVP** basis
and B3LYP variant; do not substitute def2-SVP or assume all B3LYP aliases
mean the same functional. Psi4 documents DZVP separately, with its own
auxiliary-basis associations. Pin the reference versions, expanded shell
data, integration grid, D3(BJ) damping parameters and whether the three-body
term is enabled, then check a small common-system energy/gradient before
using a second implementation as an equivalent reference.
[Psi4 basis-family definitions](https://psicode.org/psi4manual/master/basissets_byfamily.html).

PySCF supports dispersion through an additional extension, and Simple DFT-D3
provides independently inspectable energy and derivative interfaces. At the
time of this initial review, our isolated QM worker had no validated D3 request/provenance path;
passing a suggestive functional string is insufficient. The private runtime
inventory inspected today contains the PySCF interface source but no installed
`pyscf.dispersion`, `dftd3`, or xTB executable in the inspected runtimes.
No new dependency was installed and no new expensive calculation was launched
for the initial review. The later qualification below adds a separate provider.
[PySCF dispersion documentation](https://pyscf.org/user/dft.html),
[Simple DFT-D3 PySCF API](https://dftd3.readthedocs.io/en/latest/api/pyscf.html).

xTB's native scan controls use restraining potentials; that is distinct from
an exactly constrained target. A prospective adapter must independently
verify the final periodic angle and frozen atoms, preserve unbiased electronic
energy separately from constraint energy, and reject changed product graphs
or stereo. Prefer the already tested geomeTRIC constraint semantics after
qualifying an xTB energy/gradient adapter, rather than silently treating a
harmonic restraint as exact.
[xTB scan documentation](https://xtb-docs.readthedocs.io/en/latest/scan.html).

All large DF jobs stay sequential unless more than 30 GiB is free and the
existing memory/process limits also permit concurrency. The current parent
job has demonstrated roughly 12 GiB of temporary DF storage. Any future
adapter must measure peak scratch use and enforce the shared launch lock;
an atom-count or three-process limit alone is insufficient.

## Qualification checkpoint, 2026-09-13

The separate research runtime and provider are now numerically qualified for
closed-shell H/C/N/O/F/S systems. The exact Psi4 Godbout DZVP orbital basis,
LibXC 402 B3LYP convention, and explicit two-body D3(BJ) constants were checked
against pinned primary sources. Four small-molecule energy/gradient checks,
the upstream D3 fixture, finite differences, and actual launch/parent-state
rejections passed. The runtime adds about 12 MiB and leaves all active QM
workers unchanged. [Method, evidence and limitations](COVALENT-REFERENCE-PROVIDER.md)

The full94-atom request plan is ready but remains unlaunched. Actual coordinates
will be materialized only from an accepted, hash-verified parent optimization
after stereo, cap and gross bonded-geometry checks. No attachment flexibility,
DFT geometry convergence, or force-field accuracy has yet been established by
this qualification. The HF ESP/RESP workflow remains unchanged.
