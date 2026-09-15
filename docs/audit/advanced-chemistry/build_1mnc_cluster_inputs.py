"""Construct explicitly declared 1MNC catalytic-Zn prototype cluster inputs.

No QM, dynamics, parameter acceptance or application backend changes. Native
MCPB supplies its published small/large capping construction. Inputs keep the
full PLH inhibitor; structural Zn/Ca sites are outside this isolated prototype.
"""
from __future__ import annotations

import hashlib
import io
import json
import os
import random
import re
import subprocess
import sys
from pathlib import Path

import gemmi
import numpy as np
from openmm import Platform, app, unit
from rdkit import Chem

ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT))
from backend.complex_topology import subset, residue_key
from backend.ligands import LigandModel, _add_hydrogens, _graph, _topology

BASE = Path(__file__).resolve().parent


def key(atom):
    r = atom.residue
    return f"{r.chain.id}:{r.id}:{r.insertionCode.strip()}:{r.name}:{atom.name}"


def pdb_atoms(path):
    return [{"serial": int(line[6:11]), "name": line[12:16].strip(),
             "resname": line[17:20].strip(), "chain": line[21:22],
             "resid": line[22:26].strip(), "xyz": [float(line[i:i+8]) for i in (30,38,46)]}
            for line in path.read_text().splitlines() if line.startswith(("ATOM  ", "HETATM"))]


def main():
    output = BASE / "prototype-1mnc-inputs"
    output.mkdir(exist_ok=False)
    cif = BASE / "inputs/1MNC.cif"
    source = app.PDBxFile(str(cif))
    source_xyz = np.asarray(source.positions.value_in_unit(unit.angstrom))
    flanks = {217,218,219,221,222,223,227,228,229}
    protein = [r for r in source.topology.residues() if r.chain.id == "A" and int(r.id) in flanks]
    if len(protein) != 9:
        raise ValueError("Expected complete source tripeptide neighborhoods are absent")
    original_protein = {key(a): source_xyz[a.index].copy() for r in protein for a in r.atoms()}
    modeller = subset(source.topology, source.positions, {residue_key(r) for r in protein}, remove_hydrogens=True)
    variants = ["HID" if r.name == "HIS" else None for r in modeller.topology.residues()]
    random.seed(2026)
    selected = modeller.addHydrogens(forcefield=None, pH=7.0, variants=variants,
                                     platform=Platform.getPlatformByName("Reference"))
    # Validate heavy coordinates before native capping/serialization. Hydrogens
    # are initial guesses and are slated for QM relaxation, not validated poses.
    hxyz = np.asarray(modeller.positions.value_in_unit(unit.angstrom))
    for atom in modeller.topology.atoms():
        if atom.element != app.element.hydrogen:
            if not np.array_equal(hxyz[atom.index], original_protein[key(atom)]):
                if not np.allclose(hxyz[atom.index], original_protein[key(atom)], atol=1e-12, rtol=0):
                    raise ValueError("Hydrogen placement moved observed protein heavy atoms")
    for r, variant in zip(modeller.topology.residues(), selected, strict=True):
        if r.name == "HIS":
            if variant != "HID" or "HE2" in {a.name for a in r.atoms()}:
                raise ValueError("Coordinating His NE2 was protonated")
            r.name = "HID"

    ligand = next(r for r in source.topology.residues() if r.name == "PLH" and r.id == "280")
    ccd = BASE / "inputs/ccd/PLH.cif"
    mol, ligand_atoms, names, stereo = _graph(ligand, source_xyz, block=gemmi.cif.read_file(str(ccd)).sole_block())
    oxygen = mol.GetAtomWithIdx(names.index("O1"))
    if oxygen.GetFormalCharge() != 0 or oxygen.GetTotalNumHs() != 1:
        raise ValueError("Expected hydroxamic-acid O1 identity/state is absent")
    oxygen.SetFormalCharge(-1)
    oxygen.SetNoImplicit(True)
    oxygen.SetNumExplicitHs(0)
    Chem.SanitizeMol(mol)
    if Chem.GetFormalCharge(mol) != -1:
        raise ValueError("Declared PLH hydroxamate charge was not obtained")
    model = LigandModel(mol, ligand, [a.index for a in ligand_atoms], names, {})
    molh, h_method = _add_hydrogens(model)
    if molh.GetAtomWithIdx(names.index("O1")).GetTotalNumHs(includeNeighbors=True) != 0:
        raise ValueError("Hydroxamate O1 regained a hydrogen")
    ligand_names = [a.GetProp("_TriposAtomName") for a in molh.GetAtoms()]
    ligand_top, ligand_pos = _topology(model, molh, ligand_names)
    (output / "PLH-declared-state.sdf").write_text(Chem.MolToMolBlock(molh)+"\n$$$$\n")
    # Intermediate mol2 charges are formal-state bookkeeping only. No BCC or
    # invented simulation charges replace the subsequent site RESP fit.
    (output / "formal-state.chg").write_text("\n".join(f"{float(a.GetFormalCharge()):.8f}" for a in molh.GetAtoms())+"\n")
    zn = next(r for r in source.topology.residues() if r.name == "ZN" and r.id == "281")
    zn_model = subset(source.topology, source.positions, {residue_key(zn)})
    modeller.add(zn_model.topology, zn_model.positions)
    modeller.add(ligand_top, ligand_pos)
    with (output / "site-input.pdb").open("w") as stream:
        app.PDBFile.writeFile(modeller.topology, modeller.positions, stream, keepIds=True)
    rows = pdb_atoms(output / "site-input.pdb")
    if len(rows) != len(list(modeller.topology.atoms())):
        raise ValueError("PDB atom count differs")
    original_map = {}
    for row, atom in zip(rows, modeller.topology.atoms(), strict=True):
        original_map[row["serial"]] = {**row, "stable_id": key(atom), "element": atom.element.symbol,
                                       "placement": "new hydrogen" if atom.element == app.element.hydrogen else "observed source heavy atom"}
    zn_serial = next(row["serial"] for row in rows if row["resname"] == "ZN")
    zn_xyz = next(row["xyz"] for row in rows if row["resname"] == "ZN")
    donor_serials = [row["serial"] for row in rows
                     if (row["resid"] in {"218","222","228"} and row["name"] == "NE2")
                     or (row["resname"] == "PLH" and row["name"] in {"O1","O2"})]
    if len(donor_serials) != 5:
        raise ValueError("Declared five-donor catalytic site is incomplete")
    (output / "ZN.mol2").write_text(
        "@<TRIPOS>MOLECULE\nZN\n1 0 1 0 0\nSMALL\nUSER_CHARGES\n\n@<TRIPOS>ATOM\n"
        +f"1 ZN {zn_xyz[0]:.4f} {zn_xyz[1]:.4f} {zn_xyz[2]:.4f} Zn 1 ZN 2.0\n"
        +"@<TRIPOS>BOND\n@<TRIPOS>SUBSTRUCTURE\n1 ZN 1 TEMP 0 **** **** 0 ROOT\n")
    amber = ROOT / ".tools/ambertools"
    env = dict(os.environ, AMBERHOME=str(amber), OMP_NUM_THREADS="2", OPENBLAS_NUM_THREADS="2")
    env["PATH"] = str(amber / "bin")+os.pathsep+env.get("PATH", "")
    def run(name, args):
        with (output/f"{name}.log").open("w") as stream:
            result = subprocess.run(args, cwd=output, env=env, stdout=stream,
                                    stderr=subprocess.STDOUT, timeout=120)
        if result.returncode:
            raise RuntimeError(f"{name} failed; see preserved log")
    run("antechamber", [str(amber/"bin/antechamber"), "-i", "PLH-declared-state.sdf", "-fi", "sdf",
        "-o", "PLH-native-types.mol2", "-fo", "mol2", "-c", "rc", "-cf", "formal-state.chg",
        "-nc", "-1", "-at", "gaff2", "-rn", "PLH", "-s", "2", "-pf", "n", "-dr", "yes", "-seq", "n"])
    # Restore the preserved atom names only after confirming native atom order,
    # coordinates, graph edges and formal charge sum.
    text = (output / "PLH-native-types.mol2").read_text()
    head, body = text.split("@<TRIPOS>ATOM\n",1)
    atomtext, tail = body.split("@<TRIPOS>BOND",1)
    typed_rows = [line.split() for line in atomtext.splitlines() if line.strip()]
    if len(typed_rows) != molh.GetNumAtoms():
        raise ValueError("Native atom typing changed ligand count")
    ligand_xyz = molh.GetConformer().GetPositions()
    for i, row in enumerate(typed_rows):
        # SDF serializes four decimals, native AC intermediates three; only
        # their combined rounding is allowed here, never coordinate relaxation.
        if int(row[0]) != i+1 or not np.allclose([float(s) for s in row[2:5]], ligand_xyz[i], atol=.00055+1e-10, rtol=0):
            raise ValueError("Native atom typing changed ligand coordinate/order")
        if row[5] in {"du", "DU"}:
            raise ValueError("Native GAFF2 lacks an atom type")
        row[1] = ligand_names[i]
    if abs(sum(float(row[8]) for row in typed_rows)+1) > 1e-6:
        raise ValueError("Native mol2 state-bookkeeping charge differs")
    native_edges = {tuple(sorted((int(row.split()[1])-1,int(row.split()[2])-1)))
                    for row in tail.split("@<TRIPOS>",1)[0].splitlines() if row.strip()}
    reference_edges = {tuple(sorted((b.GetBeginAtomIdx(),b.GetEndAtomIdx()))) for b in molh.GetBonds()}
    if native_edges != reference_edges:
        raise ValueError("Native typing changed ligand graph")
    (output / "PLH.mol2").write_text(head+"@<TRIPOS>ATOM\n"+"\n".join(" ".join(r) for r in typed_rows)+"\n@<TRIPOS>BOND"+tail)
    run("parmchk2", [str(amber/"bin/parmchk2"), "-i", "PLH.mol2", "-f", "mol2", "-o", "PLH.frcmod", "-s", "gaff2"])
    if "ATTN" in (output / "PLH.frcmod").read_text():
        raise ValueError("Native GAFF2 reported missing ligand terms")
    (output / "site.in").write_text(
        "original_pdb site-input.pdb\ngroup_name site\ncut_off 2.8\n"
        +f"ion_ids {zn_serial}\nion_mol2files ZN.mol2\nnaa_mol2files PLH.mol2\nfrcmod_files PLH.frcmod\n"
        +"force_field ff14SB\nwater_model tip3p\ngaff 2\nlarge_opt 1\n"
        +"smmodel_chg 1\nsmmodel_spin 1\nlgmodel_chg 1\nlgmodel_spin 1\n")
    run("mcpb-step1", [str(amber/"bin/python"), str(amber/"bin/MCPB.py"), "-i", "site.in", "-s", "1"])
    links = [line.split()[1:] for line in (output / "site_standard.fingerprint").read_text().splitlines() if line.startswith("LINK")]
    actual_edges = {tuple(sorted(int(value.split("-",1)[0]) for value in row)) for row in links}
    expected_edges = {tuple(sorted((zn_serial, donor))) for donor in donor_serials}
    if actual_edges != expected_edges:
        raise ValueError(f"Native model coordination graph differs: {actual_edges ^ expected_edges}")
    clusters = {}
    for size in ("small", "large"):
        cluster = pdb_atoms(output / f"site_{size}.pdb")
        atom_map = []
        for row in cluster:
            original = original_map.get(row["serial"])
            retained = original and row["name"] == original["name"] and row["resname"] == original["resname"]
            atom_map.append({**row, "role": "retained" if retained else "native MCPB cap",
                             "source_or_placed_id": original["stable_id"] if original else None,
                             "stable_id": original["stable_id"] if retained else f"cap:{size}:{row['serial']}:{row['name']}"})
            if retained and not np.array_equal(row["xyz"], original["xyz"]):
                raise ValueError("Native cluster construction moved retained atom")
        com = (output / f"site_{size}_{'opt' if size=='small' else 'mk'}.com").read_text()
        qrows = []
        for line in com.splitlines():
            pieces = line.split()
            if len(pieces) in {4,5} and re.fullmatch(r"[A-Z][a-z]?", pieces[0]):
                qrows.append((pieces[0], [float(v) for v in pieces[-3:]]))
        if len(qrows) != len(atom_map) or not np.array_equal([r[1] for r in qrows], [r["xyz"] for r in atom_map]):
            raise ValueError("Native QM/PDB order differs")
        request = {"schema_version": 1, "atom_ids": [r["stable_id"] for r in atom_map],
                   "elements": [r[0] for r in qrows], "coords_angstrom": [r[1] for r in qrows],
                   "charge": 1, "spin": 0, "method": "RKS", "functional": "B3LYPG",
                   "basis": "6-31g*", "grid_level": 4, "density_fit": False,
                   "threads": 2, "max_memory_mb": 8000, "operations": ["gradient"],
                   "max_wall_seconds": 7200}
        (output/f"{size}-screening-qm-request.json").write_text(json.dumps(request, indent=2)+"\n")
        (output/f"{size}-atom-map.json").write_text(json.dumps(atom_map, indent=2)+"\n")
        clusters[size] = {"atom_count": len(atom_map), "retained": sum(r["role"] == "retained" for r in atom_map),
                          "caps": sum(r["role"] != "retained" for r in atom_map)}
    report = {"status": "cluster_inputs_created_not_parameterized", "pdb": "1MNC",
              "source_sha256": hashlib.sha256(cif.read_bytes()).hexdigest(),
              "ccd_sha256": hashlib.sha256(ccd.read_bytes()).hexdigest(),
              "declared_state": {"Zn281": "Zn(II), closed-shell singlet", "His218/222/228": "HID; NE2 unprotonated",
                                 "PLH280": "hydroxamate candidate: O1 deprotonated (-1), N1 remains protonated",
                                 "basis": "explicit model hypothesis, not a state inferred from neutral CCD or distances",
                                 "cluster_charge": 1, "spin_2S": 0},
              "ligand_heavy_stereocenters_verified": stereo, "initial_ligand_hydrogens": h_method,
              "input_mol2_charges": "formal-state bookkeeping only; these are not accepted simulation partial charges",
              "native_coordination_edges": sorted(map(list, actual_edges)), "clusters": clusters,
              "excluded_from_isolated_prototype": ["structural Zn282", "Ca283", "nonlocal protein and solvent"],
              "not_claimed": "Full 1MNC preparation, QM convergence, site force-field accuracy or dynamics stability"}
    (output/"cluster-input-report.json").write_text(json.dumps(report, indent=2)+"\n")
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()
