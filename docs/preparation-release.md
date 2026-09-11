# Preparation and solvent preview extension — September 10, 2026

Added independent plot-overlay visibility; Simulate source upload/fetch/SMILES; real background protein preparation; sequence-aware missing-region warnings and opt-in construction; and actual TIP3P preview with exact OpenMM state reuse. Preparation checks template handedness and rejects invalid models. Narrow panels keep the molecule above scrolling setup controls.

Validation: **48 backend tests**, **13 real Chrome tests**, and two additional camera/responsive-layout checks passed. Independent numerical and plot checks passed 8 and 6 groups. The latest user-visible ubiquitin example has 18,589 atoms including 17,358 water atoms; preparation checked 84 stereocenters. This is an unequilibrated starting system, not a scientific result. Source fetching was tested against actual RCSB/PubChem responses. Browser tests verified visible water removal using rendered pixels and preserved the original input.

Detailed records: [release validation](audit/preparation-release-validation.json), [preparation review](preparation-review.md), and [narrow-panel screenshot](dynamol-preparation-tablet.png). Scientific accuracy/convergence promotion remains outside the validated scope. Prepared GROMACS conversion and arbitrary ligand parameterization remain unsupported.
