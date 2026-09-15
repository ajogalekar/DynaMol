"""Canonical MDTraj storage; analysis always uses original physical coordinates."""
import json
import re
import uuid
from pathlib import Path

import mdtraj as md
import numpy as np

from . import config
from .ions import ION_STATES, SUPPORTED_IONS

WATERS = {"HOH", "WAT", "SOL", "TIP3", "TIP3P"}
IONS = SUPPORTED_IONS
TOPOLOGY_EXTENSIONS = {".pdb", ".pdbx", ".cif", ".mmcif", ".gro", ".h5", ".hdf5"}
TRAJECTORY_EXTENSIONS = {".dcd", ".xtc", ".trr", ".nc", ".netcdf", ".pdb", ".h5", ".hdf5", ".mdcrd", ".crd", ".lammpstrj", ".xyz"}


def safe_id(value: str) -> str:
    if not re.fullmatch(r"[A-Za-z0-9_-]{1,80}", value):
        raise ValueError("Invalid identifier.")
    return value


def atomic_json(path: Path, payload):
    temporary = path.with_suffix(path.suffix + "." + uuid.uuid4().hex + ".tmp")
    temporary.write_text(json.dumps(payload, indent=2, allow_nan=False))
    temporary.replace(path)


def dataset_dir(dataset_id: str) -> Path:
    return config.DATASETS_DIR / safe_id(dataset_id)


def get_dataset(dataset_id: str) -> dict:
    path = dataset_dir(dataset_id) / "metadata.json"
    if not path.exists():
        raise FileNotFoundError(f"Dataset '{dataset_id}' was not found.")
    return json.loads(path.read_text())


def list_datasets() -> list[dict]:
    results = []
    for path in sorted(config.DATASETS_DIR.glob("*/metadata.json"), key=lambda p: p.stat().st_mtime, reverse=True):
        try:
            meta = json.loads(path.read_text())
            results.append({**meta, "atoms": [], "bonds": []})
        except (ValueError, OSError):
            continue
    return results


def load_physical(dataset_id: str) -> md.Trajectory:
    folder = dataset_dir(dataset_id)
    get_dataset(dataset_id)
    with np.load(folder / "physical.npz", allow_pickle=False) as data:
        top = md.load_topology(str(folder / "topology.pdb"))
        chemistry_path = folder / "chemistry.json"
        if chemistry_path.is_file():
            chemistry = json.loads(chemistry_path.read_text())
            if chemistry.get("authoritative_bonds"):
                # PDB cannot retain bond orders and may infer extra bonds from
                # residue names. Rebuild from the uploaded chemical graph.
                if len(chemistry["atoms"]) != top.n_atoms:
                    raise ValueError("Stored chemical graph does not match the coordinate atom count.")
                restored = md.Topology()
                chains, residues, atoms = {}, {}, []
                for chain in top.chains:
                    chains[chain.index] = restored.add_chain(chain.chain_id)
                for residue in top.residues:
                    residues[residue.index] = restored.add_residue(residue.name, chains[residue.chain.index], residue.resSeq, residue.segment_id)
                for atom in top.atoms:
                    chemical = chemistry["atoms"][atom.index]
                    if chemical["index"] != atom.index:
                        raise ValueError("Stored chemical graph has inconsistent atom ordering.")
                    atoms.append(restored.add_atom(atom.name, atom.element, residues[atom.residue.index], atom.serial, chemical.get("formal_charge")))
                for bond in chemistry["bonds"]:
                    a, b = bond["atoms"]
                    if not (0 <= a < len(atoms) and 0 <= b < len(atoms)) or a == b:
                        raise ValueError("Stored chemical graph contains an invalid bond.")
                    order = bond.get("order")
                    restored.add_bond(atoms[a], atoms[b], type=md.core.topology.Aromatic if bond.get("aromatic") else None, order=int(order) if order is not None and float(order).is_integer() else None)
                top = restored
        traj = md.Trajectory(data["xyz"], top, time=data["time"])
        if "lengths" in data:
            traj.unitcell_lengths = data["lengths"]
            traj.unitcell_angles = data["angles"]
    return traj


def check_size(traj: md.Trajectory):
    if traj.n_atoms == 0 or traj.n_frames == 0:
        raise ValueError("The file contains no atoms or frames.")
    if traj.n_atoms > config.MAX_ATOMS:
        raise ValueError(f"Too many atoms: the local viewer supports at most {config.MAX_ATOMS:,} atoms.")
    if traj.n_frames > config.MAX_FRAMES or traj.xyz.nbytes > config.MAX_COORD_BYTES:
        raise ValueError("Trajectory is too large for this local viewer. Increase stride or trim the trajectory (10,000 frames and 256 MiB of coordinates maximum).")
    if not np.isfinite(traj.xyz).all():
        raise ValueError("The trajectory contains non-finite coordinates.")
    if not np.isfinite(traj.time).all():
        raise ValueError("The trajectory contains non-finite timestamps.")


def load_uploaded(topology: Path, trajectory: Path | None, stride: int, frame_interval_ps: float | None) -> tuple[md.Trajectory, list[str]]:
    if topology.suffix.lower() not in TOPOLOGY_EXTENSIONS:
        raise ValueError("Supported topology formats: PDB, mmCIF, GRO, and MDTraj HDF5.")
    if trajectory and trajectory.suffix.lower() not in TRAJECTORY_EXTENSIONS:
        raise ValueError("Supported trajectory formats: DCD, XTC, TRR, NetCDF, PDB, HDF5, MDCRD, LAMMPS, and XYZ.")
    warnings = []
    source = trajectory or topology
    kwargs = {} if source.suffix.lower() in {".pdb", ".cif", ".mmcif", ".pdbx", ".h5", ".hdf5"} else {"top": str(topology)}
    # Chunked decoding enforces caps before retaining the whole trajectory in memory.
    try:
        chunks = []
        frames = 0
        size = 0
        reader = [md.load(str(source))[::stride]] if source.suffix.lower() in {".cif", ".mmcif", ".pdbx"} else md.iterload(str(source), chunk=100, stride=stride, **kwargs)
        for chunk in reader:
            check_size(chunk)
            frames += chunk.n_frames
            size += chunk.xyz.nbytes
            if frames > config.MAX_FRAMES or size > config.MAX_COORD_BYTES:
                raise ValueError("Trajectory exceeds the local memory cap. Increase stride or trim the file (10,000 frames / 256 MiB coordinates).")
            chunks.append(chunk)
        if not chunks:
            raise ValueError("The selected file has no readable frames.")
        traj = md.join(chunks, check_topology=True)
        if trajectory:
            provided_topology = md.load_topology(str(topology))
            if provided_topology.n_atoms != traj.n_atoms:
                raise ValueError("Topology and trajectory have different atom counts; use matching files in the same atom order.")
            if source.suffix.lower() in {".pdb", ".cif", ".mmcif", ".pdbx", ".h5", ".hdf5"}:
                def signature(top):
                    return [(a.name, a.residue.name, a.residue.resSeq, a.residue.chain.chain_id or a.residue.chain.index) for a in top.atoms]
                if signature(provided_topology) != signature(traj.topology):
                    raise ValueError("Atom identities or ordering differ between the supplied topology and the self-describing trajectory.")
            traj.topology = provided_topology
    except Exception as exc:
        raise ValueError(f"Could not read the molecular files: {exc}") from exc
    if frame_interval_ps is not None:
        if not np.isfinite(frame_interval_ps) or frame_interval_ps <= 0:
            raise ValueError("Frame interval must be a positive number of picoseconds.")
        traj.time = np.arange(traj.n_frames, dtype=float) * frame_interval_ps * stride
        warnings.append(f"Timestamps supplied by user: {frame_interval_ps:g} ps per original frame (stride {stride}).")
    elif source.suffix.lower() in {".pdb", ".gro", ".cif", ".mmcif", ".pdbx", ".dcd", ".xyz", ".mdcrd", ".crd", ".lammpstrj"}:
        traj.time = np.arange(traj.n_frames, dtype=float) * stride
        warnings.append("Physical timestamps are unavailable from this file/reader; values are original frame indices, not measured picoseconds. Supply a frame interval for time-based interpretation.")
    if trajectory:
        warnings.append("Atom counts checked. Binary trajectories cannot verify atom identities: the topology must use the same atom ordering.")
    check_size(traj)
    return traj, warnings


def save_dataset(traj: md.Trajectory, name: str, source: str, description: str, warnings: list[str] | None = None, dataset_id: str | None = None, provenance: dict | None = None, *, exact_pdb: str | None = None) -> dict:
    check_size(traj)
    dataset_id = safe_id(dataset_id or uuid.uuid4().hex[:16])
    folder = dataset_dir(dataset_id)
    folder.mkdir(parents=True, exist_ok=True)
    warnings = list(warnings or [])
    physical = {"xyz": traj.xyz.astype(np.float32), "time": np.asarray(traj.time, dtype=np.float64)}
    has_cell = traj.unitcell_vectors is not None and bool(np.all(np.isfinite(traj.unitcell_vectors))) and bool(np.all(traj.unitcell_volumes > 0))
    if has_cell:
        physical.update(lengths=traj.unitcell_lengths, angles=traj.unitcell_angles)
    else:
        traj.unitcell_vectors = None
        warnings.append("No periodic box: measurements use ordinary Cartesian geometry.")
    np.savez_compressed(folder / "physical.npz", **physical)
    pdb_frame = traj[0]
    pdb_frame.topology = traj.topology.copy()
    # mmCIF readers can retain string atom IDs. The canonical display PDB needs
    # integer serials; renumber a copy without changing source order or bytes.
    for atom in pdb_frame.topology.atoms:
        atom.serial = atom.index + 1
    if exact_pdb is None:
        pdb_frame.save_pdb(str(folder / "topology.pdb"))
    else:
        # MDTraj cannot represent residue insertion codes. Native structure
        # writers supply the exact PDB before the canonical read so distinct
        # residues sharing an author number cannot collapse into one residue.
        (folder / "topology.pdb").write_text(exact_pdb)
    # PDB round trip avoids metadata and measurement connectivity disagreeing.
    canonical = md.load_topology(str(folder / "topology.pdb"))
    if canonical.n_atoms != traj.n_atoms:
        raise ValueError("Canonical topology serialization changed the atom inventory; exact source residue identities are required.")
    traj.topology = canonical
    display = traj.slice(slice(None), copy=True)
    if has_cell:
        try:
            molecules = display.topology.find_molecules()
            largest = max(molecules, key=len)
            display.image_molecules(anchor_molecules=[largest], inplace=True)
            warnings.append("Viewer molecules are made whole across periodic boundaries; analysis uses the original coordinates with minimum-image geometry.")
        except Exception:
            warnings.append("Periodic imaging unavailable for this topology; molecules may cross display box edges.")
    from .residue_identity import protein_residue_keys, residue_key
    protein_keys = protein_residue_keys(dataset_id, canonical, traj.xyz[0])
    protein_ca = np.array([a.index for a in canonical.atoms if residue_key(a.residue) in protein_keys and a.name == "CA"], dtype=int)
    if len(protein_ca) >= 3:
        display.superpose(display, 0, atom_indices=protein_ca, parallel=False)
        warnings.append("Display frames are aligned on protein alpha carbons; plotted geometry uses unaligned simulation frames.")
    else:
        display.xyz -= display.xyz[0].mean(axis=0)[None, None, :]
    (display.xyz * 10).astype("<f4").tofile(folder / "coordinates.bin")
    bonds = [[a.index, b.index] for a, b in canonical.bonds]
    adjacency: dict[int, list[int]] = {}
    for a, b in bonds:
        adjacency.setdefault(a, []).append(b)
        adjacency.setdefault(b, []).append(a)
    atom_list = list(canonical.atoms)
    atoms = []
    for atom in atom_list:
        element = atom.element.symbol if atom.element else "X"
        residue = atom.residue
        ion_state = ION_STATES.get(residue.name.upper())
        is_ion = residue.n_atoms == 1 and ion_state is not None and element == ion_state[0]
        # Some aliases (CAL) also appear in amino-acid name inventories. Exact
        # monatomic identity wins; a multiatom molecule never takes this path.
        category = "ions" if is_ion else "protein" if residue_key(residue) in protein_keys else "nucleic" if residue.is_nucleic else "water" if residue.name.upper() in WATERS else "ligands"
        atoms.append({"index": atom.index, "name": atom.name, "element": element, "residue": residue.name, "resid": residue.resSeq, "chain": residue.chain.chain_id or str(residue.chain.index + 1), "category": category,
                      "nonpolar_hydrogen": element == "H" and any(atom_list[b].element and atom_list[b].element.symbol == "C" for b in adjacency.get(atom.index, []))})
    metadata = {"id": dataset_id, "name": name, "n_atoms": traj.n_atoms, "n_residues": traj.n_residues, "n_frames": traj.n_frames, "times_ps": traj.time.tolist(), "time_unit": "frame" if any("Physical timestamps are unavailable" in w for w in warnings) else "ps", "source": source, "description": description,
                "topology_url": f"/api/datasets/{dataset_id}/topology", "coordinates_url": f"/api/datasets/{dataset_id}/coordinates", "atoms": atoms, "bonds": bonds, "has_unitcell": has_cell, "warnings": warnings}
    atomic_json(folder / "metadata.json", metadata)
    if provenance:
        atomic_json(folder / "provenance.json", provenance)
    return metadata
