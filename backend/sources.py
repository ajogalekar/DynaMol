"""Bounded structure import and fixed-provider retrieval with source provenance.

Small-molecule chemistry is kept in a sidecar because PDB is not a chemistry
exchange format.  Generated conformers are deliberately distinct from MD inputs.
"""
import datetime as dt
import hashlib
import io
import json
import re
import shutil
import tempfile
import urllib.error
import urllib.parse
import urllib.request
import uuid
from pathlib import Path

import mdtraj as md
import numpy as np

from . import config
from .storage import atomic_json, dataset_dir, save_dataset

MAX_SOURCE_BYTES = min(config.MAX_UPLOAD_BYTES, 25 * 1024 * 1024)
MAX_SMALL_MOLECULE_ATOMS = 500
SOURCE_EXTENSIONS = {".pdb", ".cif", ".mmcif", ".pdbx", ".mol2", ".mol", ".sdf", ".smi", ".smiles"}
STANDARD_RESIDUES = {"ALA", "ARG", "ASN", "ASP", "CYS", "GLN", "GLU", "GLY", "HIS", "ILE", "LEU", "LYS", "MET", "PHE", "PRO", "SER", "THR", "TRP", "TYR", "VAL", "HID", "HIE", "HIP", "CYX", "ASH", "GLH", "LYN"}


def _name(name, fallback):
    return str(name or fallback).strip()[:100] or fallback


def _seed(seed):
    if not isinstance(seed, int) or isinstance(seed, bool) or not 1 <= seed <= 2_147_483_646:
        raise ValueError("Choose an integer random seed from 1 to 2147483646.")
    return seed


def _source_metadata(data, suffix, provenance):
    return {**(provenance or {}), "source_file": "source" + suffix, "source_format": suffix.lstrip("."),
            "source_sha256": hashlib.sha256(data).hexdigest(), "retrieved_at": dt.datetime.now(dt.timezone.utc).isoformat()}


def _save(traj, name, data, suffix, provenance, warnings, chemistry=None, sdf=None, native=None):
    dataset_id = uuid.uuid4().hex[:16]
    try:
        return _save_files(traj, name, data, suffix, provenance, warnings, chemistry, sdf, dataset_id, native)
    except Exception:
        # A serializer or sidecar failure must never publish a partial import.
        shutil.rmtree(dataset_dir(dataset_id), ignore_errors=True)
        raise


def _save_files(traj, name, data, suffix, provenance, warnings, chemistry, sdf, dataset_id, native=None):
    provenance = _source_metadata(data, suffix, provenance)
    original_chains = [chain.chain_id for chain in traj.topology.chains]
    exact_pdb = None
    if native is not None:
        from openmm import app, unit
        output = io.StringIO()
        app.PDBFile.writeFile(native.topology, traj.xyz[0] * unit.nanometer, output, keepIds=True)
        exact_pdb = output.getvalue()
    metadata = save_dataset(traj, name, provenance.get("provider", "upload"),
                            "Imported molecular structure" if not provenance.get("conformer") else "Computed 3D conformer for viewing; MD parameters are not assigned.",
                            warnings=warnings, provenance=provenance, dataset_id=dataset_id, exact_pdb=exact_pdb)
    folder = dataset_dir(metadata["id"])
    (folder / provenance["source_file"]).write_bytes(data)
    metadata.update(source_file=provenance["source_file"], source_format=provenance["source_format"])
    canonical_chains = list(traj.topology.chains)
    mapping = [{"index": i, "original": original, "canonical": chain.chain_id} for i, (original, chain) in enumerate(zip(original_chains, canonical_chains))]
    metadata["chain_id_mapping"] = mapping
    provenance["chain_id_mapping"] = mapping
    atomic_json(folder / "provenance.json", provenance)
    if chemistry:
        atomic_json(folder / "chemistry.json", chemistry)
        metadata["bonds"] = [[b["atoms"][0], b["atoms"][1]] for b in chemistry["bonds"]]
        metadata["bond_orders"] = chemistry["bonds"]
        for atom, chemical in zip(metadata["atoms"], chemistry["atoms"]):
            atom.update({key: value for key, value in chemical.items() if key != "index"})
        adjacency = {}
        for a, b in metadata["bonds"]:
            adjacency.setdefault(a, []).append(b)
            adjacency.setdefault(b, []).append(a)
        for atom in metadata["atoms"]:
            atom["nonpolar_hydrogen"] = atom["element"] == "H" and any(metadata["atoms"][i]["element"] == "C" for i in adjacency.get(atom["index"], []))
        metadata["chemistry"] = {key: value for key, value in chemistry.items() if key not in {"atoms", "bonds"}}
    if sdf:
        (folder / "molecule.sdf").write_text(sdf)
    atomic_json(folder / "metadata.json", metadata)
    return metadata


def _rdkit_trajectory(mol):
    """Create a topology directly, retaining RDKit atom indices and bonds."""
    top = md.Topology()
    chain = top.add_chain("A")
    residue = top.add_residue("LIG", chain, resSeq=1)
    counts = {}
    atoms = []
    for atom in mol.GetAtoms():
        symbol = atom.GetSymbol()
        counts[symbol] = counts.get(symbol, 0) + 1
        atoms.append(top.add_atom(f"{symbol}{counts[symbol]}", md.element.get_by_symbol(symbol), residue))
    for bond in mol.GetBonds():
        top.add_bond(atoms[bond.GetBeginAtomIdx()], atoms[bond.GetEndAtomIdx()])
    return md.Trajectory(np.asarray(mol.GetConformer().GetPositions(), dtype=np.float32)[None] / 10, top, time=[0])


def _rdkit_chemistry(mol, extra=None):
    from rdkit import Chem, rdBase
    heavy = Chem.RemoveHs(mol)
    return {"schema_version": 1, "authoritative_bonds": True, "rdkit_version": rdBase.rdkitVersion,
            "canonical_isomeric_smiles": Chem.MolToSmiles(heavy, isomericSmiles=True),
            "total_formal_charge": Chem.GetFormalCharge(mol),
            "unassigned_stereocenters": sum(label == "?" for _, label in Chem.FindMolChiralCenters(heavy, includeUnassigned=True)),
            "atoms": [{"index": a.GetIdx(), "formal_charge": a.GetFormalCharge(), "aromatic": a.GetIsAromatic(), "isotope": a.GetIsotope()} for a in mol.GetAtoms()],
            "bonds": [{"atoms": [b.GetBeginAtomIdx(), b.GetEndAtomIdx()], "order": b.GetBondTypeAsDouble(), "aromatic": b.GetIsAromatic()} for b in mol.GetBonds()], **(extra or {})}


def _embed(mol, seed):
    from rdkit import Chem
    from rdkit.Chem import AllChem
    if mol.GetNumAtoms() > MAX_SMALL_MOLECULE_ATOMS or len(Chem.GetMolFrags(mol)) != 1:
        raise ValueError("3D generation supports one connected molecule with at most 500 atoms. Load separate molecules individually.")
    mol = Chem.AddHs(mol)
    if mol.GetNumAtoms() > MAX_SMALL_MOLECULE_ATOMS:
        raise ValueError("3D generation is limited to 500 atoms including hydrogens.")
    parameters = AllChem.ETKDGv3()
    parameters.randomSeed = _seed(seed)
    parameters.numThreads = min(2, config.CPU_THREADS)
    parameters.maxIterations = 1000
    if hasattr(parameters, "timeout"):
        parameters.timeout = 20
    parameters.enforceChirality = True
    if AllChem.EmbedMolecule(mol, parameters) != 0:
        raise ValueError("A 3D conformer could not be generated for this SMILES. Try a supplied 3D SDF or MOL2 structure.")
    if AllChem.MMFFHasAllMoleculeParams(mol):
        method = "MMFF94s"
        status = AllChem.MMFFOptimizeMolecule(mol, mmffVariant=method, maxIters=500)
    elif AllChem.UFFHasAllMoleculeParams(mol):
        method = "UFF"
        status = AllChem.UFFOptimizeMolecule(mol, maxIters=500)
    else:
        method, status = "none (parameters unavailable)", None
    return mol, {"method": "RDKit ETKDGv3", "seed": seed, "optimization": method, "optimization_converged": status == 0,
                 "interpretation": "Computed conformer; neither experimental coordinates nor a parameterized MD system."}


def create_smiles(smiles, name=None, seed=2026):
    from rdkit import Chem
    if not isinstance(smiles, str) or not smiles.strip() or len(smiles) > 10_000 or len(smiles.splitlines()) != 1:
        raise ValueError("Enter one SMILES string (maximum 10,000 characters).")
    smiles = smiles.strip()
    mol = Chem.MolFromSmiles(smiles)
    if mol is None or mol.GetNumAtoms() == 0:
        raise ValueError("The SMILES could not be parsed. Check atoms, bond symbols, rings, and stereochemistry.")
    mol, conformer = _embed(mol, seed)
    warnings = ["Computed 3D conformer from SMILES; protonation and tautomer state come from your input. Ligand MD parameters are not assigned."]
    chemistry = _rdkit_chemistry(mol, {"input_smiles": smiles, "conformer": conformer})
    if chemistry["unassigned_stereocenters"]:
        warnings.append("The input has unspecified stereocenters; the displayed conformer does not establish their experimental configuration.")
    if not conformer["optimization_converged"]:
        warnings.append("Conformer optimization did not converge or parameters were unavailable; inspect the generated geometry.")
    return _save(_rdkit_trajectory(mol), _name(name, "SMILES molecule"), smiles.encode(), ".smiles",
                 {"provider": "smiles", "conformer": conformer, "input_smiles": smiles}, warnings, chemistry,
                 Chem.MolToMolBlock(mol) + "\n$$$$\n")


def _mol2(data):
    """Read Tripos coordinates, residue identities, explicit connectivity/charges.

    Multiple molecules and underspecified connectivity are rejected, never guessed.
    Protein MOL2 requires residue records; PDB/mmCIF is better for sequence evidence.
    """
    text = data.decode("utf-8")
    if text.upper().count("@<TRIPOS>MOLECULE") != 1:
        raise ValueError("Load a MOL2 file containing exactly one molecular structure.")
    sections = {}
    current = None
    for line in text.splitlines():
        if line.upper().startswith("@<TRIPOS>"):
            current = line.strip().upper()[9:]
            sections[current] = []
        elif current and line.strip() and not line.lstrip().startswith("#"):
            sections[current].append(line.split())
    rows = sections.get("ATOM", [])
    if not rows or len(rows) > config.MAX_ATOMS:
        raise ValueError("The MOL2 contains no atoms or exceeds the local atom cap.")
    top = md.Topology()
    chains, residues, atoms, coords, chemistry_atoms = {}, {}, {}, [], []
    substructures = {int(r[0]): r for r in sections.get("SUBSTRUCTURE", [])}
    previous_residue, previous_chain = None, None
    atom_names = {}
    for row in rows:
        if len(row) < 6:
            raise ValueError("A MOL2 atom record is incomplete.")
        serial, atom_name = int(row[0]), row[1]
        if serial in atoms:
            raise ValueError("MOL2 atom IDs must be unique.")
        # Only unambiguous Tripos element labels; GAFF types must use SDF instead.
        raw_element = row[5].split(".")[0]
        element = md.element.get_by_symbol(raw_element)
        resid = int(row[6]) if len(row) > 6 else 1
        if not -999 <= resid <= 9999:
            raise ValueError("MOL2 residue IDs must fit the current PDB-based viewer (-999 through 9999).")
        original_residue = row[7] if len(row) > 7 else "LIG"
        match = re.fullmatch(r"([A-Za-z]{3})(?:\d+)?", original_residue)
        residue_name = match[1].upper() if match else original_residue[:3].upper()
        substructure = substructures.get(resid, [])
        chain_id = substructure[5] if len(substructure) > 5 and substructure[5] not in {"****", "0"} else "A"
        if len(chain_id) != 1 or len(atom_name) > 4:
            raise ValueError("This MOL2 has chain or atom names that cannot be preserved in the current viewer. Convert to PDB/mmCIF with stable atom identities.")
        key = (chain_id, resid)
        if key in residues and key != previous_residue:
            raise ValueError("MOL2 atoms from each residue must be contiguous; interleaved residues cannot preserve atom ordering. Export an ordered PDB/mmCIF or MOL2.")
        if chain_id in chains and previous_chain != chain_id:
            raise ValueError("MOL2 atoms from each chain must be contiguous to preserve atom ordering.")
        if atom_name in atom_names.setdefault(key, set()):
            raise ValueError("MOL2 atom names must be unique within each residue to preserve atom identities.")
        atom_names[key].add(atom_name)
        previous_residue, previous_chain = key, chain_id
        if key not in residues:
            if chain_id not in chains:
                chains[chain_id] = top.add_chain(chain_id)
            chain = chains[chain_id]
            residues[key] = top.add_residue(residue_name, chain, resSeq=resid)
        elif residues[key].name != residue_name:
            raise ValueError("MOL2 assigns conflicting names to one residue ID.")
        atom = top.add_atom(atom_name, element, residues[key], serial=serial)
        atoms[serial] = atom
        coords.append([float(value) for value in row[2:5]])
        partial_charge = float(row[8]) if len(row) > 8 else None
        if partial_charge is not None and not np.isfinite(partial_charge):
            raise ValueError("MOL2 contains non-finite partial charges.")
        chemistry_atoms.append({"index": atom.index, "input_atom_id": serial, "input_atom_type": row[5], "input_residue": original_residue, "partial_charge": partial_charge})
    chemistry_bonds, seen = [], set()
    orders = {"1": 1.0, "2": 2.0, "3": 3.0, "ar": 1.5, "am": 1.0}
    for row in sections.get("BOND", []):
        if len(row) < 4 or row[3].lower() not in orders:
            raise ValueError("MOL2 contains an unsupported or unspecified bond type. Supply explicit bonds in SDF, MOL2, or PDB.")
        a, b = atoms[int(row[1])], atoms[int(row[2])]
        pair = tuple(sorted([a.index, b.index]))
        if a is b or pair in seen:
            raise ValueError("MOL2 contains duplicated or self bonds.")
        seen.add(pair)
        top.add_bond(a, b)
        chemistry_bonds.append({"atoms": [a.index, b.index], "order": orders[row[3].lower()], "aromatic": row[3].lower() == "ar", "input_type": row[3]})
    if len(atoms) > 1 and not chemistry_bonds:
        raise ValueError("MOL2 has no bond records; connectivity cannot be preserved.")
    for residue in top.residues:
        if residue.name in STANDARD_RESIDUES and not {"N", "CA", "C"}.issubset({a.name for a in residue.atoms}):
            raise ValueError("A MOL2 amino-acid residue lacks backbone atom identities. Use a protein PDB/mmCIF; the structure was not relabeled as a ligand.")
    molecule_records = sections.get("MOLECULE", [])
    if len(molecule_records) > 2 and molecule_records[2][0].upper() in {"PROTEIN", "BIOPOLYMER"} and not any(residue.name in STANDARD_RESIDUES for residue in top.residues):
        raise ValueError("Protein MOL2 requires amino-acid residue identities. Supply a protein PDB/mmCIF instead; the protein was not relabeled as a ligand.")
    chemistry = {"schema_version": 1, "authoritative_bonds": True, "atoms": chemistry_atoms, "bonds": chemistry_bonds,
                 "input_format": "mol2", "charge_note": "MOL2 partial charges retained as source data; compatibility with an MD force field is not established."}
    return md.Trajectory(np.asarray(coords, dtype=np.float32)[None] / 10, top, time=[0]), chemistry


def import_structure(path, name=None, provenance=None):
    path = Path(path)
    suffix = path.suffix.lower()
    if suffix not in SOURCE_EXTENSIONS:
        raise ValueError("Choose PDB, mmCIF, MOL2, SDF, MOL, SMI, or SMILES.")
    if path.stat().st_size > MAX_SOURCE_BYTES:
        raise ValueError("Structure imports are limited to 25 MiB.")
    data = path.read_bytes()
    if not data:
        raise ValueError("The selected structure file is empty.")
    if suffix in {".smi", ".smiles"}:
        rows = [line.strip() for line in data.decode("utf-8").splitlines() if line.strip() and not line.startswith("#")]
        if len(rows) != 1:
            raise ValueError("Load a SMILES file with exactly one nonempty structure line.")
        fields = rows[0].split(maxsplit=1)
        metadata = create_smiles(fields[0], name=name or (fields[1] if len(fields) == 2 else path.stem))
        folder = dataset_dir(metadata["id"])
        old = json.loads((folder / "provenance.json").read_text())
        combined = _source_metadata(data, suffix, {**old, **(provenance or {})})
        (folder / combined["source_file"]).write_bytes(data)
        atomic_json(folder / "provenance.json", combined)
        metadata.update(source_file=combined["source_file"], source_format=combined["source_format"])
        atomic_json(folder / "metadata.json", metadata)
        return metadata
    chemistry, sdf, native = None, None, None
    warnings = []
    try:
        if suffix == ".mol2":
            traj, chemistry = _mol2(data)
            warnings.append("MOL2 connectivity, residue labels, bond orders, and input partial charges are retained. MOL2 has no reference polymer sequence; missing-residue identities cannot be recovered from numbering alone.")
        elif suffix in {".mol", ".sdf"}:
            from rdkit import Chem
            supplier = list(Chem.ForwardSDMolSupplier(io.BytesIO(data), sanitize=True, removeHs=False))
            if len(supplier) != 1 or supplier[0] is None:
                raise ValueError("Choose an SDF/MOL containing exactly one readable molecule.")
            mol = supplier[0]
            if mol.GetNumAtoms() > MAX_SMALL_MOLECULE_ATOMS:
                raise ValueError("Small-molecule imports are limited to 500 atoms.")
            if not mol.GetNumConformers():
                raise ValueError("This file has no molecular coordinates; use the SMILES field to generate a 3D conformer.")
            if not mol.GetConformer().Is3D():
                mol, conformer = _embed(Chem.RemoveHs(mol), 2026)
                provenance = {**(provenance or {}), "conformer": conformer}
                warnings.append("The input was 2D. A computed ETKDGv3 3D conformer was generated with seed 2026; it is not experimental geometry.")
            traj = _rdkit_trajectory(mol)
            chemistry = _rdkit_chemistry(mol)
            if (provenance or {}).get("conformer"):
                chemistry["conformer"] = provenance["conformer"]
                warnings.append("This is a computed 3D conformer, not an experimentally resolved molecular structure.")
            sdf = Chem.MolToMolBlock(mol) + "\n$$$$\n"
            warnings.append("Chemical graph, formal charges, and stereochemistry are retained. Ligand MD parameters are not assigned.")
        else:
            if suffix == ".pdb":
                traj = md.load_frame(str(path), 0)
                from openmm import app
                native = app.PDBFile(str(path))
            else:
                from openmm import app, unit
                import gemmi
                # OpenMM reads struct_conn bonds without applying their crystal
                # symmetry. A different-image endpoint is not an atom in this
                # loaded coordinate model; keep those records in the original
                # source while excluding them from a private parser view.
                document = gemmi.cif.read_file(str(path))
                structure = gemmi.make_structure_from_block(document.sole_block())
                from .sequence_evidence import different_image_connection_ids
                different_images = different_image_connection_ids(document.sole_block())
                omitted = [link.name for link in structure.connections if link.name in different_images]
                if omitted:
                    table = document.sole_block().find_mmcif_category("_struct_conn.")
                    id_column = list(table.tags).index("_struct_conn.id")
                    for index in reversed(range(len(table))):
                        if gemmi.cif.as_string(table[index][id_column]) in omitted:
                            table.remove_row(index)
                    provenance = {**(provenance or {}), "excluded_symmetry_connection_ids": omitted,
                                  "symmetry_connection_scope": "Different crystallographic images are absent from the loaded coordinate unit; original connection records remain in the unchanged source."}
                    warnings.append(f"{len(omitted)} source connection(s) to crystallographic symmetry mates are recorded but are not bonds between the loaded atoms.")
                cif = app.PDBxFile(io.StringIO(document.as_string()))
                native = cif
                traj = md.Trajectory(np.asarray(cif.positions.value_in_unit(unit.nanometer))[None], md.Topology.from_openmm(cif.topology), time=[0])
                vectors = cif.topology.getPeriodicBoxVectors()
                if vectors is not None:
                    traj.unitcell_vectors = np.asarray(vectors.value_in_unit(unit.nanometer))[None]
            warnings.append("The first coordinate model is loaded. Original source and polymer sequence records are retained for protein preparation.")
    except (ValueError, KeyError, IndexError, TypeError, UnicodeDecodeError) as exc:
        raise ValueError(f"Could not import this structure: {exc}") from exc
    if not any(a.residue.is_protein for a in traj.topology.atoms):
        warnings.append("This structure can be viewed and measured; the current automatic MD presets require a standard protein.")
    return _save(traj, _name(name, path.stem), data, suffix, provenance, warnings, chemistry, sdf, native)


class _NoRedirect(urllib.request.HTTPRedirectHandler):
    def redirect_request(self, req, fp, code, msg, headers, newurl):
        raise ValueError("The structure provider redirected the request; automatic redirects are disabled.")


def _download(url):
    parsed = urllib.parse.urlparse(url)
    if parsed.scheme != "https" or parsed.netloc not in {"files.rcsb.org", "pubchem.ncbi.nlm.nih.gov"}:
        raise ValueError("Only fixed RCSB PDB and PubChem HTTPS endpoints are supported.")
    request = urllib.request.Request(url, headers={"User-Agent": "DynaMol/0.1 (open-source local molecular viewer)", "Accept": "*/*"})
    try:
        with urllib.request.build_opener(_NoRedirect).open(request, timeout=20) as response:
            length = response.headers.get("Content-Length")
            if length and int(length) > MAX_SOURCE_BYTES:
                raise ValueError("The provider response exceeds the 25 MiB import limit.")
            data = response.read(MAX_SOURCE_BYTES + 1)
            if len(data) > MAX_SOURCE_BYTES:
                raise ValueError("The provider response exceeds the 25 MiB import limit.")
            return data
    except urllib.error.HTTPError as exc:
        if exc.code == 404:
            raise ValueError("No structure matched that identifier at the selected provider.") from exc
        raise ValueError(f"The structure provider returned HTTP {exc.code}. Retry later or upload a local file.") from exc
    except (urllib.error.URLError, TimeoutError, OSError) as exc:
        raise ValueError("Could not reach the structure provider. Retry later or upload a local file.") from exc


def fetch_structure(provider, identifier, name=None, seed=2026):
    identifier = str(identifier).strip()
    _seed(seed)
    if provider == "pdb":
        if not re.fullmatch(r"(?:[1-9][A-Za-z0-9]{3}|pdb_[A-Za-z0-9]{8})", identifier, flags=re.IGNORECASE):
            raise ValueError("Enter an RCSB PDB accession such as 1UBQ or pdb_00001ubq, not a URL.")
        accession = identifier.lower()
        url = f"https://files.rcsb.org/download/{accession}.cif"
        data = _download(url)
        provenance = {"provider": "rcsb_pdb", "identifier": identifier.upper(), "url": url, "assembly": "deposited coordinate model; no biological-assembly reconstruction"}
        suffix, default_name = ".cif", identifier.upper()
    elif provider == "pubchem":
        if not identifier or len(identifier) > 200 or any(ord(c) < 32 for c in identifier) or "://" in identifier:
            raise ValueError("Enter a PubChem CID or compound name (maximum 200 characters), not a URL.")
        cid = identifier
        lookup_url = None
        if not identifier.isdigit():
            lookup_url = "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/" + urllib.parse.quote(identifier, safe="") + "/cids/JSON"
            try:
                payload = json.loads(_download(lookup_url))
                cids = payload["IdentifierList"]["CID"]
                if len(cids) != 1:
                    raise ValueError("The name matches multiple PubChem records; enter one specific PubChem CID.")
                cid = str(cids[0])
            except (KeyError, json.JSONDecodeError) as exc:
                raise ValueError("PubChem returned no unambiguous compound identifier.") from exc
        if not re.fullmatch(r"[1-9][0-9]{0,11}", cid):
            raise ValueError("Enter a positive PubChem compound CID.")
        url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/{cid}/SDF?record_type=3d"
        data = _download(url)
        provenance = {"provider": "pubchem", "identifier": identifier, "cid": int(cid), "url": url,
                      "lookup_url": lookup_url, "conformer": {"method": "PubChem computed 3D conformer", "interpretation": "Computed PubChem conformer; not experimental coordinates or MD parameters."}}
        suffix, default_name = ".sdf", name or f"{identifier} · PubChem {cid}"
    else:
        raise ValueError("Choose the PDB or PubChem structure provider.")
    with tempfile.TemporaryDirectory(prefix="dynamol-fetch-") as temporary:
        path = Path(temporary) / ("structure" + suffix)
        path.write_bytes(data)
        return import_structure(path, _name(name, default_name), provenance)
