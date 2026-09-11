"""Regenerate small UI-only fixtures from the bundled demo and ethanol.

Run with .venv/bin/python docs/audit/ui-audit-fixtures/generate.py from DynaMol.
"""
import hashlib
import json
from pathlib import Path

import mdtraj as md
import numpy as np
from openmm import unit
from openmm.app import PDBxFile
from rdkit import Chem
from rdkit.Chem import AllChem

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
source = ROOT / "examples/demo/physical.npz"
raw = np.load(source)
topology = md.load_pdb(str(ROOT / "examples/demo/topology.pdb")).topology
indices = [atom.index for atom in topology.atoms if atom.residue.index < 6]
trajectory = md.Trajectory(raw["xyz"][:5], topology, time=np.arange(5) * 0.2).atom_slice(indices)
for extension in ["pdb", "gro", "h5", "xtc", "dcd", "trr", "nc", "xyz", "mdcrd"]:
    trajectory.save(str(HERE / f"short.{extension}"))
periodic = trajectory[:]
periodic.unitcell_lengths = np.full((periodic.n_frames, 3), 10.0)
periodic.unitcell_angles = np.full((periodic.n_frames, 3), 90.0)
periodic.save_lammpstrj(str(HERE / "short.lammpstrj"))
with (HERE / "short.cif").open("w") as output:
    PDBxFile.writeFile(trajectory.topology.to_openmm(), trajectory.xyz[0] * unit.nanometer, output)

mixed = trajectory.topology.copy()
positions = list(trajectory.xyz[0])
center = np.mean(trajectory.xyz[0], axis=0)
for name, elements, offset in [
    ("LIG", ["C", "C", "O"], [0.5, 0, 0]),
    ("HOH", ["O", "H", "H"], [0, 0.5, 0]),
    ("NA", ["Na"], [0, 0, 0.5]),
]:
    residue = mixed.add_residue(name, mixed.add_chain())
    atoms = []
    for index, symbol in enumerate(elements):
        atoms.append(mixed.add_atom(f"{symbol}{index + 1}" if name != "NA" else "NA", md.element.get_by_symbol(symbol), residue))
        positions.append(center + np.array(offset) + np.array([index * 0.12, (index % 2) * 0.07, 0]))
    if name == "LIG":
        mixed.add_bond(atoms[0], atoms[1])
        mixed.add_bond(atoms[1], atoms[2])
    if name == "HOH":
        mixed.add_bond(atoms[0], atoms[1])
        mixed.add_bond(atoms[0], atoms[2])
md.Trajectory(np.array(positions)[None], mixed).save_pdb(str(HERE / "mixed-groups.pdb"))

# A schematic phosphodiester backbone tests nucleic-acid category handling.
# It has no bases, no validated stereochemistry, and is not suitable for MD.
nucleic = md.Topology()
chain = nucleic.add_chain("N")
nucleic_positions = []
previous = None
for index, name in enumerate(["DA", "DT", "DG", "DC"]):
    residue = nucleic.add_residue(name, chain, resSeq=index + 1)
    local = []
    for j, (atom_name, symbol) in enumerate([("P", "P"), ("O5'", "O"), ("C5'", "C"), ("C4'", "C"), ("C3'", "C"), ("O3'", "O")]):
        atom = nucleic.add_atom(atom_name, md.element.get_by_symbol(symbol), residue)
        local.append(atom)
        nucleic_positions.append([index * 0.72 + j * 0.12, np.sin(index + j * 0.4) * 0.12, np.cos(index + j * 0.4) * 0.12])
        if j:
            nucleic.add_bond(local[-2], atom)
    if previous is not None:
        nucleic.add_bond(previous, local[0])
    previous = local[-1]
md.Trajectory(np.array(nucleic_positions)[None], nucleic).save_pdb(str(HERE / "schematic-nucleic.pdb"))

ethanol = Chem.AddHs(Chem.MolFromSmiles("CCO"))
AllChem.EmbedMolecule(ethanol, randomSeed=2026)
AllChem.UFFOptimizeMolecule(ethanol)
(HERE / "ethanol.mol").write_text(Chem.MolToMolBlock(ethanol))
with Chem.SDWriter(str(HERE / "ethanol.sdf")) as writer:
    writer.write(ethanol)
for extension in ["smi", "smiles"]:
    (HERE / f"ethanol.{extension}").write_text("CCO ethanol\n")
(HERE / "manifest.json").write_text(json.dumps({
    "source": "bundled examples/demo physical coordinates, first 5 frames and first 6 residues",
    "source_sha256": hashlib.sha256(source.read_bytes()).hexdigest(),
    "n_atoms": trajectory.n_atoms,
    "n_frames": trajectory.n_frames,
    "time_ps": trajectory.time.tolist(),
    "first_pair_distance_angstrom": (np.linalg.norm(trajectory.xyz[:, 0] - trajectory.xyz[:, 1], axis=1) * 10).tolist(),
    "mixed_fixture": "Visualization-only arbitrary placed ethanol, water and sodium; not an MD-ready complex.",
    "lammps_fixture": "Original demo coordinates in an artificial 10 nm cubic cell, for file decoding only.",
    "nucleic_fixture": "Schematic four-residue DNA backbone without bases or validated stereochemistry, for visibility only; not suitable for MD.",
    "ethanol": {"smiles": "CCO", "generator": "RDKit EmbedMolecule + UFFOptimizeMolecule", "seed": 2026},
}, indent=2) + "\n")
