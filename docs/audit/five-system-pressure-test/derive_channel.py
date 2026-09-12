"""Build the recorded 1K4C channel tetramer, without Fab or duplicate axis ions."""
from pathlib import Path
import hashlib
import json

import gemmi
import numpy as np


folder = Path(__file__).resolve().parent
source_path = folder / "review-source-1K4C.cif"
source = gemmi.read_structure(str(source_path))
model = gemmi.Model("1")
source_channel = source[0]["C"]
operators = source.assemblies[0].generators[0].operators
assert len(operators) == 4
protein_counts = []
water_candidates = []
for chain_id, operator in zip("ABCD", operators):
    chain = gemmi.Chain(chain_id)
    for residue in source_channel:
        if residue.entity_type == gemmi.EntityType.Polymer:
            copy = residue.clone()
            copy.subchain = chain_id
            copy.entity_id = "1"
            for atom in copy:
                position = operator.transform.apply(atom.pos)
                atom.pos = gemmi.Position(position.x, position.y, position.z)
            chain.add_residue(copy)
        elif residue.is_water():
            copy = residue.clone()
            for atom in copy:
                position = operator.transform.apply(atom.pos)
                atom.pos = gemmi.Position(position.x, position.y, position.z)
            water_candidates.append(copy)
    assert len(chain) == 103
    protein_counts.append(len(chain))
    model.add_chain(chain)

selected_ions = [3002, 3004, 3005]
ions = gemmi.Chain("I")
ion_xyz = []
for residue in source_channel:
    if residue.name == "K" and residue.seqid.num in selected_ions:
        copy = residue.clone()
        copy.subchain = "I"
        copy.entity_id = "2"
        for atom in copy:
            # One concrete particle is placed per selected site; raw occupancy stays in provenance.
            atom.occ = 1
            ion_xyz.append([atom.pos.x, atom.pos.y, atom.pos.z])
        ions.add_residue(copy)
assert len(ions) == 3
model.add_chain(ions)

# Retain only observed nearby waters, deduplicating special-position copies.
waters = gemmi.Chain("W")
water_xyz = []
for residue in water_candidates:
    oxygen = next((atom for atom in residue if atom.element.name == "O"), None)
    if oxygen is None:
        continue
    xyz = np.array([oxygen.pos.x, oxygen.pos.y, oxygen.pos.z])
    if np.linalg.norm(np.array(ion_xyz) - xyz, axis=1).min() > 4:
        continue
    if water_xyz and np.linalg.norm(np.array(water_xyz) - xyz, axis=1).min() < 1:
        continue
    residue.seqid = gemmi.SeqId(len(waters) + 1, " ")
    residue.subchain = "W"
    residue.entity_id = "3"
    for atom in residue:
        atom.occ = 1
    waters.add_residue(residue)
    water_xyz.append(xyz)
if len(waters):
    model.add_chain(waters)

output = gemmi.Structure()
output.name = "1K4C channel tetramer - water-only software test"
output.add_model(model)
entity = gemmi.Entity("1")
entity.entity_type = gemmi.EntityType.Polymer
entity.polymer_type = gemmi.PolymerType.PeptideL
entity.full_sequence = list(source.get_entity("3").full_sequence)
entity.subchains = list("ABCD")
output.entities.append(entity)
output.setup_entities()
target = folder / "channel-tetramer.cif"
output.make_mmcif_document().write_file(str(target))
record = {
    "source_url": "https://files.rcsb.org/download/1K4C.cif",
    "source_sha256": hashlib.sha256(source_path.read_bytes()).hexdigest(),
    "derived_sha256": hashlib.sha256(target.read_bytes()).hexdigest(),
    "assembly": "1", "operators": [op.name for op in operators],
    "source_protein_entity": "3", "source_chain": "C", "output_chains": list("ABCD"),
    "protein_residues_per_chain": protein_counts, "selected_ion_author_ids": selected_ions,
    "ions_copied_once": True, "ion_coordinates_angstrom": ion_xyz,
    "retained_observed_water_count": len(waters), "retained_observed_water_coordinates_angstrom": [p.tolist() for p in water_xyz],
    "excluded": "Fab entities 1/2; crystallization lipids; K3001/3003/3006/3007; bulk waters not within 4 Å of selected K; near-coincident symmetry duplicates",
    "interpretation": "Explicit discrete ion occupancy and water-only software-test scenario, not inferred experimental simultaneous occupancy or a realistic membrane model",
    "terminal_omissions": "Deposited channel construct has no observed residues 1–21; they remain omitted",
}
(folder / "channel-derivation.json").write_text(json.dumps(record, indent=2) + "\n")
print(json.dumps({"file": str(target), "protein_residues": sum(protein_counts), "ions": len(ions), "waters": len(waters)}))
