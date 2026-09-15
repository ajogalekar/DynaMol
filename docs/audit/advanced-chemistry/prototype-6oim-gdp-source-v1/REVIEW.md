# Retained GDP: source and native mechanics audit

The published GDP template loads with all 28 heavy atoms and the same heavy-atom
connectivity as the frozen GDP CCD. Its eight negative periodic-torsion
amplitudes exposed an overly restrictive DynaMol validator check. The validator
now preserves signed finite Fourier amplitudes and compares their exact values
against the independent native manifest. Bond/angle stiffness, finite values,
positive periodicities, constraints and complete force coverage remain checked.

The [contributor files](https://personalpages.manchester.ac.uk/staff/Richard.Bryce/amber/cof/phos_inf.html)
identify the Meagher–Redman–Carlson polyphosphate model. The maintainer's old
index hyperlinks returned 404, but the same PREP, FRCMOD and contributor page are
available on its personalpages host. Original bytes, URLs and hashes are in
`retrieval.json`; these files were not installed in the app or included in a release.

The [published FRCMOD](https://personalpages.manchester.ac.uk/staff/Richard.Bryce/amber/cof/frcmod.phos)
explicitly contains negative amplitudes. For the
[OpenMM periodic expression](https://docs.openmm.org/latest/userguide/theory/02_standard_forces.html#periodictorsionforce),
`k(1 + cos(nθ − phase))`, a negative coefficient remains bounded. Taking its
absolute value would change the surface. No coefficient, phase or energy offset
was modified to make the comparison pass.

Native AmberTools TLEAP loaded the unedited PREP/FRCMOD with their Amber99 base
in 0.39 seconds. The isolated template has 40 atoms and total charge approximately
−3 e. The deposited CCD's neutral reference charge is not silently substituted
for this model state. Name mapping allows only the declared star-to-prime
normalization. This is not a protonation-equilibrium determination.

The native specimen and original failed validator result remain in
`build/advanced-chemistry/covalent/gdp-source-native-v1`. After the code fix, all
40 atoms, 42 bonds, 73 angles, 135 torsion terms and 204 exclusions/scaled pairs
match the independently exported native manifest. The focused validator suite
passed 47 tests, including an analytic negative-amplitude energy, deliberate
sign mismatch, nonfinite coefficients and negative harmonic-stiffness rejection.

This establishes source-graph and parameter-transport consistency only. Native
PREP-generated coordinates produce a 0.373 Å nonbonded contact warning; they
are not the experimental 6OIM coordinates or a ready starting structure. No
MD or QM was run for GDP. Complete 6OIM assembly must retain its experimental
GDP/Mg/water identities, repair the declared protein loop, establish the selected
nucleotide/ion/water model's applicability, and validate unrestrained motion.
The covalent adduct still awaits its actual charges and conformational evidence.

Primary work on [GDP/GTP–Mg models](https://pubmed.ncbi.nlm.nih.gov/23280996/)
also demonstrates sensitivity to magnesium–oxygen interactions; availability
of a GDP template alone does not settle the coordination model. Download
permission is recorded separately from permission to redistribute a parameter
bundle. Nothing here promotes this candidate into app readiness.
