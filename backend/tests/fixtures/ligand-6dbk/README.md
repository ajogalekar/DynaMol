# 6DBK ligand identity regression fixture

These three small files were copied unchanged from the failed DynaMol preparation of PDB 6DBK. `provenance.json` records their source and SHA-256 hashes. They contain one 54-atom ligand and its native Amber parameters, coordinates, and converted OpenMM XML.

The regression compares the native Amber system with the XML system at the retained pose and two deterministic perturbations. The original graph-only template assignment selected a symmetry-related atom mapping, changing four nonzero improper torsion quartets. It failed the existing strict energy/force check. The tests retain that failure as a negative control, then check that an explicit validated atom identity map restores native terms, energies, and forces without changing these files or relaxing the tolerances.

This is a parameter-transfer and software regression fixture. Agreement does not validate the force field, protonation state, binding pose, or simulation sampling for scientific inference. Charge fitting is not repeated in this test suite.
