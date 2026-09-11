# Third-party software and data

DynaMol's MIT license covers its original application code. It does not replace the licenses of bundled fonts or installed scientific/JavaScript dependencies. Full license texts are included in the respective installed packages and upstream repositories.

- [NGL](https://github.com/nglviewer/ngl) supplies the molecular renderer and depends on Three.js.
- [React](https://github.com/facebook/react), [Vite](https://github.com/vitejs/vite), [Lucide](https://github.com/lucide-icons/lucide), [TypeScript](https://github.com/microsoft/TypeScript), and [Playwright](https://github.com/microsoft/playwright) supply interface/build/test components.
- [FastAPI](https://github.com/fastapi/fastapi), [Uvicorn](https://github.com/encode/uvicorn), [NumPy](https://github.com/numpy/numpy), and [MDTraj](https://github.com/mdtraj/mdtraj) supply the local API and scientific processing. MDTraj is LGPL software; its license and source remain available separately.
- [OpenMM](https://github.com/openmm/openmm) and [GROMACS](https://gitlab.com/gromacs/gromacs) are independently maintained MD engines, with their own licenses, force-field attribution, scientific references, and distribution requirements. GROMACS binaries are installed separately and are not part of this source repository.
- [PDBFixer](https://github.com/openmm/pdbfixer) supplies missing-atom and sequence-supported residue reconstruction; [RDKit](https://github.com/rdkit/rdkit) supplies small-molecule parsing and conformer generation. Their upstream licenses and scientific references apply independently. Optional structure fetch uses [RCSB PDB](https://www.rcsb.org/) and [PubChem](https://pubchem.ncbi.nlm.nih.gov/).
- [DM Sans](https://github.com/googlefonts/dm-fonts) and [Manrope](https://github.com/sharanda/manrope), served locally through Fontsource, are licensed under the SIL Open Font License.
- Experimental coordinates for the bundled example come from [wwPDB / RCSB PDB entry 1UBQ](https://www.rcsb.org/structure/1UBQ). PDB archive data are made freely available under [wwPDB's data policy](https://www.wwpdb.org/about/usage-policies). The example's simulated coordinates were generated with OpenMM and are labeled separately from the experimental input.

Dependency versions are recorded by `uv.lock`, `frontend/package-lock.json`, and each simulation's provenance. Installing the application is not an endorsement by any upstream project or author.
