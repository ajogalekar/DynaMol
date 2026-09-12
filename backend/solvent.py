"""Create an actual, bounded TIP3P solvent box for preview and MD reuse."""
import hashlib
import math
from functools import lru_cache
from pathlib import Path

import mdtraj as md
import numpy as np

from . import config
from .prepared_system import copy_ligand_parameters, ligand_parameter_files, load_prepared_forcefield
from .storage import atomic_json, dataset_dir, get_dataset, save_dataset


@lru_cache(maxsize=4)
def _tip3p_lattice(path, mtime_ns, size):
    """Read the same installed template as Modeller; cache by file revision."""
    from openmm import app, unit
    source = Path(path)
    pdb = app.PDBFile(str(source))
    coordinates = np.asarray(pdb.positions.value_in_unit(unit.nanometer))
    periods = np.asarray(pdb.topology.getUnitCellDimensions().value_in_unit(unit.nanometer))
    oxygens = []
    for residue in pdb.topology.residues():
        atoms = list(residue.atoms())
        if len(atoms) != 3 or sum(a.element == app.element.oxygen for a in atoms) != 1 or sum(a.element == app.element.hydrogen for a in atoms) != 2:
            raise ValueError("The installed TIP3P template is not a three-site water lattice; a safe allocation bound could not be established.")
        oxygens.append(coordinates[next(a.index for a in atoms if a.element == app.element.oxygen)])
    if not oxygens or not np.isfinite(oxygens).all() or not np.isfinite(periods).all() or np.any(periods <= 0):
        raise ValueError("The installed TIP3P solvent lattice has invalid coordinates or dimensions.")
    return np.asarray(oxygens), periods, hashlib.sha256(source.read_bytes()).hexdigest()


def _accepted_axis_tiles(oxygen, period, side):
    """Count native inclusive-boundary placements without allocating a tile grid."""
    total = math.ceil(side / period)
    lower, upper = 0, total
    # Binary search the literal native floating-point expression instead of
    # flooring (side-oxygen)/period, which can undercount exact boundaries.
    while lower < upper:
        middle = (lower + upper) // 2
        if oxygen + middle * period <= side:
            lower = middle + 1
        else:
            upper = middle
    return lower


def tip3p_allocation_bound(positions, padding_nm=1):
    """Upper bound for OpenMM's cubic TIP3P construction, before all deletions.

    Count every tiled template water accepted by the native box-boundary test.
    Solute exclusions, water-edge filtering and zero-molar monovalent ion
    replacement can only reduce this count. This is not a density estimate.
    """
    from openmm import Vec3, unit, version
    from openmm.app import modeller as implementation
    coordinates = np.asarray(positions.value_in_unit(unit.nanometer) if hasattr(positions, 'value_in_unit') else positions, dtype=float)
    if coordinates.ndim != 2 or coordinates.shape[1] != 3 or not np.isfinite(coordinates).all():
        raise ValueError("Finite N×3 nanometer coordinates are required to bound solvent allocation.")
    if isinstance(padding_nm, bool) or not math.isfinite(float(padding_nm)) or float(padding_nm) <= 0:
        raise ValueError("A finite positive solvent padding is required.")
    if len(coordinates):
        minimum = Vec3(*(float(v) for v in coordinates.min(axis=0)))
        maximum = Vec3(*(float(v) for v in coordinates.max(axis=0)))
        center = .5 * (minimum + maximum)
        # Mirror native Vec3/unit.norm arithmetic, including center definition.
        radius = max(unit.norm(center - Vec3(*(float(v) for v in point))) for point in coordinates)
    else:
        radius = 0
    side = max(2 * radius + float(padding_nm), 2 * float(padding_nm))
    if not math.isfinite(side):
        raise ValueError("The requested box dimensions are too large for bounded solvent allocation.")
    # Include neighboring representable boundaries conservatively. No change
    # is made to the native box width or physical solute coordinates.
    count_side = math.nextafter(math.nextafter(side, math.inf), math.inf)
    source = Path(implementation.__file__).parent / 'data' / 'tip3p.pdb'
    stat = source.stat()
    oxygens, periods, digest = _tip3p_lattice(str(source), stat.st_mtime_ns, stat.st_size)
    waters = sum(math.prod(_accepted_axis_tiles(float(oxygen[i]), float(periods[i]), count_side) for i in range(3)) for oxygen in oxygens)
    return {"box_side_nm": side, "count_boundary_nm": count_side, "candidate_water_molecules": waters,
            "maximum_added_water_atoms": 3 * waters, "maximum_total_atoms": len(coordinates) + 3 * waters,
            "template_box_nm": periods.tolist(), "template_sha256": digest, "openmm_version": version.version,
            "method": "Cubic TIP3P template oxygen lattice count before solute/edge exclusions or neutralizing monovalent-ion replacement; inclusive native upper boundaries plus two upward floating-point steps. No tiled grid is allocated; final actual atom cap is also enforced."}


def solvate_dataset(dataset_id, padding_nm=1, seed=2026, ph=7):
    """Executed in a dedicated preparation worker, never on the API event loop.

    Only an explicitly prepared protein or parameterized complex is accepted. Existing DynaMol previews
    are reused, or rebuilt from their parent when padding changes.  There is no
    second hydrogen/protonation pass and no automatic heterogen deletion.
    """
    import random
    from openmm import app, unit, version

    if not isinstance(padding_nm, (int, float)) or not math.isfinite(padding_nm) or not 1 <= padding_nm <= 3:
        raise ValueError("Choose solvent padding from 1 to 3 nm.")
    if not isinstance(seed, int) or isinstance(seed, bool) or not 1 <= seed <= 2_147_483_646:
        raise ValueError("Choose an integer random seed from 1 to 2147483646.")
    if not isinstance(ph, (int, float)) or not math.isfinite(ph) or not 0 <= ph <= 14:
        raise ValueError("Choose a pH from 0 to 14.")
    parent = get_dataset(dataset_id)
    ligand_parameter_files(dataset_dir(dataset_id), parent.get("preparation"))
    existing = parent.get("solvation")
    if existing:
        if existing.get("padding_nm") == padding_nm and existing.get("seed") == seed:
            return parent
        dataset_id = existing["parent_dataset_id"]
        parent = get_dataset(dataset_id)
    preparation = parent.get("preparation")
    source = dataset_dir(dataset_id) / "prepared.pdb"
    if not preparation or not source.is_file():
        raise ValueError("Prepare the protein first, then choose explicit water. This builds a real TIP3P box around the prepared hydrogen and ionization states.")
    if any(atom["category"] == "nucleic" for atom in parent["atoms"]):
        raise ValueError("Nucleic acids need a separately parameterized workflow; no atoms were removed.")
    if any(atom["category"] == "ligands" for atom in parent["atoms"]):
        ligand_parameter_files(dataset_dir(dataset_id), preparation, required=True)
    if not any(atom["category"] == "protein" for atom in parent["atoms"]):
        raise ValueError("Choose a prepared standard protein before adding explicit water.")
    from .modified_residues import register_topology_definitions
    register_topology_definitions()
    pdb = app.PDBFile(str(source))
    coordinates = np.asarray(pdb.positions.value_in_unit(unit.nanometer))
    if not len(coordinates) or not np.isfinite(coordinates).all():
        raise ValueError("Prepared coordinates are empty or non-finite.")
    allocation_bound = tip3p_allocation_bound(pdb.positions, padding_nm)
    if allocation_bound['maximum_total_atoms'] > config.MAX_ATOMS:
        raise ValueError(f"The requested solvent box may exceed the {config.MAX_ATOMS:,}-atom local limit. Reduce padding or prepare a smaller structure.")
    forcefield, forcefield_files = load_prepared_forcefield(dataset_dir(dataset_id), preparation)
    modeller = app.Modeller(pdb.topology, pdb.positions)
    try:
        unmatched = forcefield.getUnmatchedResidues(modeller.topology)
        if unmatched:
            names = ", ".join(sorted({r.name for r in unmatched}))
            raise ValueError(f"Prepared residues do not match the explicit-water force field: {names}. Prepare the protein again and inspect missing atoms or unsupported residues.")
        # The caller uses a dedicated subprocess. Restore the RNG even in direct
        # tests; the recorded seed determines counterion water replacement.
        previous_state = random.getstate()
        try:
            random.seed(seed)
            modeller.addSolvent(forcefield, model="tip3p", padding=padding_nm * unit.nanometer,
                                boxShape="cube", neutralize=True, ionicStrength=0 * unit.molar)
        finally:
            random.setstate(previous_state)
        if modeller.topology.getNumAtoms() > allocation_bound['maximum_total_atoms']:
            raise ValueError("Native solvent construction exceeded its template-derived allocation bound; the preview is rejected for inspection.")
        if modeller.topology.getNumAtoms() > config.MAX_ATOMS:
            raise ValueError(f"The generated box exceeds {config.MAX_ATOMS:,} atoms. Reduce solvent padding.")
        # Template validation is not an equilibration or a stability claim.
        forcefield.createSystem(modeller.topology, nonbondedMethod=app.PME,
                                nonbondedCutoff=1 * unit.nanometer, constraints=app.HBonds)
    except Exception as exc:
        raise ValueError(f"Could not build the explicit-water preview: {exc}") from exc
    top = md.Topology.from_openmm(modeller.topology)
    traj = md.Trajectory(np.asarray(modeller.positions.value_in_unit(unit.nanometer), dtype=np.float32)[None], top, time=[0])
    traj.unitcell_vectors = np.asarray(modeller.topology.getPeriodicBoxVectors().value_in_unit(unit.nanometer))[None]
    water_atoms = sum(a.residue.name in {"HOH", "WAT", "SOL"} for a in modeller.topology.atoms())
    old_water_atoms = sum(a["category"] == "water" for a in parent["atoms"])
    ions = sum(1 for r in modeller.topology.residues() if r.name in {"NA", "CL"})
    details = {"parent_dataset_id": dataset_id, "padding_nm": padding_nm, "water_model": "tip3p", "seed": seed,
               "added_water_atoms": water_atoms - old_water_atoms, "water_atoms": water_atoms, "neutralizing_ions": ions,
               "ionic_strength_molar": 0, "box_shape": "cube", "box_vectors_nm": traj.unitcell_vectors[0].tolist(),
               "openmm_version": version.version, "forcefield": forcefield_files,
               "prepared_pdb": "prepared.pdb", "input_sha256": hashlib.sha256(source.read_bytes()).hexdigest(),
               "allocation_bound": allocation_bound, "actual_atoms": modeller.topology.getNumAtoms(),
               "equilibrated": False, "preserves_prepared_protonation": True}
    provenance = {"operation": "explicit_water_preview", "parent_dataset_id": dataset_id, "preparation": preparation, "solvation": details}
    warnings = ["Real TIP3P coordinates and periodic box; this is an unequilibrated starting system. Run minimization and equilibration before interpreting dynamics.",
                "The prepared hydrogen and protonation states are retained. OpenMM simulations reuse this saved system without adding a second solvent box."]
    if old_water_atoms:
        warnings.append("Existing prepared water molecules were retained; new solvent was placed around the input without duplicating those waters.")
    if preparation.get("ph") is not None and preparation["ph"] != ph:
        warnings.append(f"The box retains the protein prepared at pH {preparation['ph']:g}. Change pH by preparing the protein again.")
    metadata = save_dataset(traj, parent["name"] + " · explicit water", "solvation",
                            "Prepared molecular system in a real TIP3P water box, ready for OpenMM minimization/equilibration.", warnings=warnings, provenance=provenance)
    folder = dataset_dir(metadata["id"])
    copy_ligand_parameters(dataset_dir(dataset_id), folder, preparation)
    with (folder / "prepared.pdb").open("w") as output:
        app.PDBFile.writeFile(modeller.topology, modeller.positions, output, keepIds=True)
    metadata.update(parent_dataset_id=dataset_id, preparation=preparation, solvation=details)
    atomic_json(folder / "metadata.json", metadata)
    return metadata
