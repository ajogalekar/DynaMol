"""Explicit, bounded completion of small missing acyclic CCD ligand fragments.

The authoritative reference determines atoms and bonds. Modeled positions are
low-confidence hypotheses; observed heavy-atom coordinates remain exact.
"""
from __future__ import annotations

import hashlib
import importlib.metadata
from pathlib import Path

import numpy as np
from openmm import app
from rdkit import Chem
from rdkit.Chem import AllChem, rdMolAlign
from rdkit.Geometry import Point3D
from scipy.spatial import cKDTree

MAX_MISSING = 8
ATTEMPTS = 8
MAX_ITERATIONS = 200


def _reference(residue, block):
    from .ligands import _graph
    table = block.get_mmcif_category("_chem_comp_atom.")
    rows = [i for i, element in enumerate(table.get("type_symbol", [])) if element.upper() not in {"H", "D"}]
    names = [table["atom_id"][i] for i in rows]
    if len(names) != len(set(names)):
        raise ValueError("The CCD reference contains ambiguous duplicate atom names.")
    top = app.Topology()
    full = top.addResidue(residue.name, top.addChain(residue.chain.id), residue.id, residue.insertionCode)
    positions = []
    for i in rows:
        top.addAtom(table["atom_id"][i], app.element.Element.getBySymbol(table["type_symbol"][i]), full)
        try:
            point = [float(table[f"pdbx_model_Cartn_{axis}_ideal"][i]) for axis in "xyz"]
        except (KeyError, ValueError, TypeError) as exc:
            raise ValueError("The full CCD ideal coordinates are unavailable; a stereochemically validated repair reference is required.") from exc
        if not np.isfinite(point).all():
            raise ValueError("The CCD ideal coordinates are not finite.")
        positions.append(point)
    try:
        mol, _, names, stereo = _graph(full, np.asarray(positions), block=block)
    except ValueError as exc:
        detail = str(exc).replace('Bound ligand stereochemistry disagrees with CCD', 'Reference coordinates disagree with CCD stereochemistry')
        raise ValueError('The CCD reference geometry could not be validated for missing-atom repair. '
                         + detail + ' Supply a validated complete ligand reference; the input structure was not modified.') from exc
    return mol, full, names, stereo


def inspect_repair(residue, coordinates, block):
    atoms = [atom for atom in residue.atoms() if atom.element != app.element.hydrogen]
    observed = [atom.name for atom in atoms]
    table = block.get_mmcif_category("_chem_comp_atom.")
    expected = [name for i, name in enumerate(table.get("atom_id", [])) if table["type_symbol"][i].upper() not in {"H", "D"}]
    missing, extra = sorted(set(expected) - set(observed)), sorted(set(observed) - set(expected))
    result = {"missing_heavy_atoms": missing, "extra_heavy_atoms": extra, "can_repair": False, "repair_reason": None}
    if not missing and not extra:
        return result
    try:
        if extra or len(observed) != len(set(observed)):
            raise ValueError("Atom identities disagree with the authoritative CCD graph; extra or duplicate atoms cannot be repaired by adding atoms.")
        if len(missing) > MAX_MISSING or len(missing) > len(expected) / 2:
            raise ValueError("This local repair supports at most eight missing heavy atoms and at most half of the ligand; supply a complete structure for a larger missing fragment.")
        if len(observed) < 3:
            raise ValueError("At least three observed heavy atoms are required to anchor ligand repair.")
        anchor_positions = np.asarray(coordinates)[[atom.index for atom in atoms]]
        if np.linalg.matrix_rank(anchor_positions - anchor_positions.mean(axis=0), tol=.05) < 2:
            raise ValueError("Observed ligand atoms are nearly collinear; their orientation does not constrain a reliable anchored repair.")
        mol, _, names, _ = _reference(residue, block)
        indices = {name: i for i, name in enumerate(names)}
        for atom in atoms:
            if mol.GetAtomWithIdx(indices[atom.name]).GetSymbol() != atom.element.symbol:
                raise ValueError(f"CCD element mismatch for observed atom {atom.name}.")
        if any(mol.GetAtomWithIdx(indices[name]).IsInRing() for name in missing):
            raise ValueError("A missing ring or ring-core atom needs a separately validated complete ligand model; local fragment repair does not reconstruct rings.")
        present = {indices[name] for name in observed}
        reached, pending = set(), [next(iter(present))]
        while pending:
            index = pending.pop()
            if index in reached:
                continue
            reached.add(index)
            pending += [atom.GetIdx() for atom in mol.GetAtomWithIdx(index).GetNeighbors() if atom.GetIdx() in present and atom.GetIdx() not in reached]
        if reached != present:
            raise ValueError("Missing atoms disconnect the observed ligand core; use a complete ligand model rather than joining uncertain fragments.")
        result.update(can_repair=True, repair_reason="A small acyclic fragment can be modeled from the full CCD graph while fixing every observed heavy atom. Newly modeled positions are low confidence and must be reviewed.")
    except ValueError as exc:
        result["repair_reason"] = str(exc)
    return result


def repair_ligand(residue, coordinates, block, *, seed=2026, on_progress=None, check_cancel=None, environment_indices=None):
    from .ligands import _graph
    inspection = inspect_repair(residue, coordinates, block)
    if not inspection["can_repair"]:
        raise ValueError(inspection["repair_reason"] or "No supported missing ligand atoms were identified.")
    reference, full_residue, names, _ = _reference(residue, block)
    observed = {atom.name: atom for atom in residue.atoms() if atom.element != app.element.hydrogen}
    fixed = {names.index(name): np.asarray(coordinates)[atom.index].copy() for name, atom in observed.items()}
    generated = [i for i, name in enumerate(names) if name not in observed]
    owned = {atom.index for atom in residue.atoms()}
    environment_indices = range(len(coordinates)) if environment_indices is None else environment_indices
    other = np.asarray([coordinates[i] for i in environment_indices if i not in owned])
    environment = cKDTree(other) if len(other) else None
    best, trials = None, []
    reference_h = Chem.AddHs(reference, addCoords=True)
    for attempt in range(ATTEMPTS):
        if check_cancel:
            check_cancel()
        if on_progress:
            on_progress(f"Modeling missing atoms in {residue.name} {residue.chain.id}:{residue.id}: candidate {attempt + 1}/{ATTEMPTS}; observed heavy atoms fixed.")
        mol = Chem.Mol(reference_h)
        parameters = AllChem.ETKDGv3()
        parameters.randomSeed = (seed + attempt) % 2_147_483_646 + 1
        parameters.numThreads = 1
        parameters.maxIterations = 200
        parameters.useRandomCoords = True
        parameters.SetCoordMap({i: Point3D(*map(float, point)) for i, point in fixed.items()})
        if AllChem.EmbedMolecule(mol, parameters) != 0:
            trials.append({"candidate": attempt + 1, "status": "embedding failed"})
            continue
        # CoordMap constrains relative geometry; align its global frame, then
        # impose exact experimental anchors as fixed optimization points.
        anchor = Chem.Mol(mol)
        for i, point in fixed.items():
            anchor.GetConformer().SetAtomPosition(i, Point3D(*map(float, point)))
        rdMolAlign.AlignMol(mol, anchor, atomMap=[(i, i) for i in fixed])
        for i, point in fixed.items():
            mol.GetConformer().SetAtomPosition(i, Point3D(*map(float, point)))
        if AllChem.MMFFHasAllMoleculeParams(mol):
            properties = AllChem.MMFFGetMoleculeProperties(mol, mmffVariant="MMFF94s")
            force = AllChem.MMFFGetMoleculeForceField(mol, properties)
            method = "MMFF94s"
        elif AllChem.UFFHasAllMoleculeParams(mol):
            force = AllChem.UFFGetMoleculeForceField(mol)
            method = "UFF"
        else:
            raise ValueError("The CCD ligand lacks local relaxation parameters; supply a complete ligand model.")
        for i in fixed:
            force.AddFixedPoint(i)
        force.Initialize()
        convergence = force.Minimize(maxIts=MAX_ITERATIONS)
        xyz = np.asarray(mol.GetConformer().GetPositions())[:len(names)]
        row = {"candidate": attempt + 1, "optimization_converged": convergence == 0, "method": method}
        trials.append(row)
        if convergence != 0 or not np.isfinite(xyz).all() or any(not np.array_equal(xyz[i], point) for i, point in fixed.items()):
            row["status"] = "nonconverged or changed observed coordinates"
            continue
        nearest = float(environment.query(xyz[generated])[0].min()) if environment else None
        if nearest is not None and nearest < 1.6:
            row.update(status="severe environment overlap", minimum_environment_distance_angstrom=nearest)
            continue
        try:
            checked, _, _, stereo = _graph(full_residue, xyz, block=block)
        except ValueError as exc:
            row.update(status="chemical or stereo validation failed", reason=str(exc))
            continue
        energy = float(force.CalcEnergy())
        row.update(status="accepted", local_energy=energy, minimum_environment_distance_angstrom=nearest)
        if best is None or energy < best[0]:
            best = energy, checked, xyz.copy(), stereo, row
    if best is None:
        raise ValueError("Constrained ligand repair could not produce a converged, stereo-consistent model without severe local overlap. Observed atoms were not moved. Supply a complete ligand structure or explicitly remove this residue.")
    _, mol, xyz, stereo, chosen = best
    report = {"method": "Full CCD graph and stereo; seeded constrained ETKDGv3 plus local fixed-anchor MMFF94s/UFF relaxation",
              "confidence": "low", "observed_heavy_atoms": list(observed), "modeled_heavy_atoms": inspection["missing_heavy_atoms"],
              "modeled_coordinates_angstrom": {names[i]: xyz[i].tolist() for i in generated}, "observed_heavy_max_displacement_angstrom": 0.0,
              "seed": seed, "attempts": trials, "selected_candidate": chosen["candidate"], "max_iterations_per_candidate": MAX_ITERATIONS,
              "minimum_environment_distance_angstrom": chosen["minimum_environment_distance_angstrom"], "gross_overlap_cutoff_angstrom": 1.6,
              "rdkit_version": importlib.metadata.version("rdkit"), "implementation_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
              "warnings": ["LOW-CONFIDENCE MODELED LIGAND ATOMS: " + ", ".join(inspection["missing_heavy_atoms"]) + ". Their identities and bonds come from the CCD, but these coordinates are modeled, not experimentally observed. The bounded conformer search and gross-overlap screen do not establish the native pose or binding energetics."]}
    return mol, names, [observed[name].index if name in observed else None for name in names], stereo, report
