# Isolated dispersion-corrected reference provider

This research provider passed the numerical checks below on 2026-09-13. It has **not** validated the covalent force field, attachment torsions, or a complete protein–ligand preparation. The production QM worker and the running RHF/RESP calculation were not changed. No full-adduct DFT calculation or torsion scan was launched during qualification.

The implementation is `covalent_reference_worker.py`; the native runtime is `.tools/covalent-reference/bin/python`. Final qualification evidence is `build/advanced-chemistry/covalent/reference-provider/qualification-v3/result.json`. The worker SHA256 is `e74a0ee5cae3fe954b918a824b4ea466aa99f46ac5311c0696eeea73e9557225`.

## Exact method

The basis is the **Godbout–Salahub–Andzelm–Wimmer DZVP DFT orbital basis**, read from the pinned Psi4 v1.9.1 `dzvp.gbs`, with spherical functions. It is not selected through a potentially ambiguous local alias. The source file cites the 1992 Canadian Journal of Chemistry paper; its SHA256 is `629b5f1d8b7eb52ab43dd84e7c0c5e3e4d0f25cff8ed858edf407788d79a968c`. H/C/N/O/F/S are the currently qualified elements. An independent reader reproduced all exponents and coefficients exactly; shell counts are 2/6/6/6/6/8 respectively. [Pinned Psi4 basis source](https://raw.githubusercontent.com/psi4/psi4/v1.9.1/psi4/share/psi4/basis/dzvp.gbs)

The electronic functional is PySCF `B3LYPG`, explicitly resolved to **LibXC 402**, the same `HYB_GGA_XC_B3LYP` identifier used in the pinned Psi4 definition. Independent mixture checks use 20% HF + 8% Slater + 72% B88 exchange, with 81% LYP + 19% VWN-RPA correlation. The VWN5 variant differs and is deliberately not substituted. Conventional RKS is used, with no density fitting; grid level 4 and grid-response terms in analytic gradients are explicit. This establishes the functional convention, not identical numerical quadratures or full-engine parity with Psi4. [Psi4 functional definition](https://raw.githubusercontent.com/psi4/psi4/v1.9.1/psi4/driver/procrouting/dft/libxc_functionals.py), [PySCF DFT documentation](https://pyscf.org/user/dft.html)

The added dispersion is **D3(BJ), two-body only**: `s6=1`, `s8=1.9889`, `a1=0.3981`, `a2=4.4211`, `s9=0`; ATM is disabled. The named parameters and these explicit values agree in s-dftd3 1.6.0. The primary parameter table associates the B3LYP BJ values with DOI 10.1002/jcc.21759. [Pinned parameter source](https://raw.githubusercontent.com/dftd3/simple-dftd3/v1.6.0/assets/parameters.toml), [Primary damping paper](https://doi.org/10.1002/jcc.21759)

## Numerical evidence

Qualification v3 completed in 9.73 seconds after package imports. Earlier v1 failed immediately because an explicit parser import was missing; that rejected attempt is retained. V2 passed before resource/state-reporting hardening; v3 is the final matching source check.

| Check | Actual result |
| --- | --- |
| Exact orbital basis, independent parser | All H/C/N/O/F/S shells and primitives match |
| Upstream 18-atom D3(BJ) energy fixture | Exactly −0.034889411100056264 Hartree |
| All 54 D3 derivative components | Maximum finite-difference error 1.22 × 10⁻¹¹ Hartree/Bohr |
| B3LYP LibXC 402 vs explicit mixture | Energy difference 4.27 × 10⁻¹⁴ Hartree; maximum gradient difference 2.00 × 10⁻¹⁴ Hartree/Bohr |
| Deliberately different VWN5 convention | Water energy changes by +0.0370760 Hartree |
| Water, ammonia, H₂S, fluoromethane | Independently constructed basis/RKS + upstream D3 adapter agree within 2.12 × 10⁻¹³ Hartree and 5.03 × 10⁻¹⁴ Hartree/Bohr |
| Total-energy derivative checks | All 9 water coordinates and 2 each on the other fixtures; maximum error 1.42 × 10⁻⁸ Hartree/Bohr |
| Invalid state, IDs, elements, units, resources, unknown settings | 11 cases rejected |

The D3 fixture comes from the upstream test suite. Direct, named, and PySCF-adapter dispersion paths share the same underlying library; they are not three independent implementations. The independent checks are the saved upstream expected value, separate basis/functional construction, and numerical differentiation. We have not run a second quantum chemistry engine on these fixtures. Numerical derivative thresholds were fixed before the runs, at 10⁻⁸ Hartree/Bohr for D3 and 2 × 10⁻⁶ Hartree/Bohr for total gradients; observed errors were smaller. [Upstream fixture](https://raw.githubusercontent.com/dftd3/simple-dftd3/v1.6.0/python/dftd3/test_pyscf.py)

`resource-and-parent-contract-tests.json` additionally records an actual rejected launch while the shared lock was held, a successful native CLI water calculation with the resource check, and rejection of materialization while the parent optimization had no final result. Rejected jobs are retained; no output is relabeled as successful.

## Runtime and result contract

The new environment occupies about **12 MiB**. It adds dftd3 1.6.0, cffi 2.1.1 and pycparser 3.0, and uses an explicit new `sitecustomize.py` to read the already verified PySCF 2.14.0/numpy/scipy runtime. The existing `.tools/qm-parallel` environment was not modified. This development environment is not a redistributable self-contained runtime. `runtime-manifest.json` pins the new D3 native binary, path setup, and previous PySCF build/parity reports. Source snapshots retain their upstream notices; their licenses are separate from DynaMol's license.

Inputs require stable unique atom IDs, one element per atom, an explicit integer charge and spin 0, exactly one coordinate array with declared units, and the fixed method name `B3LYP-D3BJ/Psi4-DZVP`. No implicit charge, basis, density-fitting or functional override is accepted. The scope is closed-shell H/C/N/O/F/S systems with at most 128 atoms.

Outputs preserve atom order and electronic state. `result.json` identifies the resolved basis, XC convention, D3 constants, versions, hashes and finite/converged status. `arrays.npz` contains coordinates in Bohr and Å, total/electronic/dispersion energies in **Hartree**, and separately named total/electronic/dispersion gradients in **Hartree/Bohr**. Total energy is the converged RKS energy plus D3; total gradient is their sum. Native parameter fitting must consume the coordinate-matched total energy, with an explicit Hartree-to-fitting-unit conversion. No synthetic qualification energies or charges enter biological models.

CLI execution holds `build/advanced-chemistry/QM-LAUNCH.lock`, counts active neutral/reference workers, and rejects a fourth quantum process. It limits each job to 2 threads and 8 GB, requiring at least 15 GiB free for models larger than 16 atoms. The current resident controller honors this lock even though its older process-name scanner does not know this new runtime. A caller must also enforce the external process-group timeout; Python alarm delivery alone cannot interrupt every native call immediately.

## Prepared full-parent request

`covalent-parent-reference-request-plan.json` contains the exact 94 atom IDs/elements, charge 0/spin 0, pinned parent input and cap graph hashes, final provider/qualification hashes, and a 4-hour single-point energy/gradient request. It intentionally has no provisional coordinates. `materialize_covalent_reference.py` creates the concrete request only after the **actual** parent optimization result exists, its hashes/state agree, and the original identity/stereochemistry/frozen-cap/bonded-geometry gates pass. It never launches a calculation.

After those gates and a new resource review, the first calculation would be **B3LYP-D3(BJ)/DZVP at the DF-RHF/6-31G* parent geometry**. That is a labeled mixed-level single point, not a DFT-optimized structure. The existing conventional HF ESP/RESP continuation remains the charge workflow. A durable outer launcher should allow 14,400 seconds native plus 120 seconds cleanup, retain rejected output, and record exact PIDs/source hashes. Further large jobs remain deferred while the present geometry, RESP and metal calculations use resources.

The first point estimates reference cost and supplies a distinct dispersion-corrected energy/gradient. It cannot validate an attachment torsion or support fitting by itself. Full-parent varied torsion geometries plus separate held-out structures are still required. The BespokeFit study supports considering B3LYP-D3(BJ)/DZVP single points on cheaper GFN2-xTB geometries as a feasibility direction; its benchmark does not establish accuracy for this 94-atom sotorasib adduct, its declared protonation state, or the GAFF2/ff14SB interface. [BespokeFit primary study](https://pmc.ncbi.nlm.nih.gov/articles/PMC9709916/)
