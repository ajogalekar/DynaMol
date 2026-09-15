"""Read-only inspection of retrieved model files; no preparation, QM, or MD."""
from pathlib import Path
import hashlib
import json
import re
import warnings

HERE = Path(__file__).resolve().parent


def residues(path):
    result = {}
    current = None
    for line in path.read_text().splitlines():
        fields = line.split("!", 1)[0].split()
        if not fields:
            continue
        if fields[0].upper() in {"RESI", "PRES"}:
            current = fields[1] if fields[0].upper() == "RESI" else None
            if current:
                result[current] = {"declared_charge": float(fields[2]), "atoms": {}}
        elif fields[0].upper() == "ATOM" and current:
            result[current]["atoms"][fields[1]] = {"type": fields[2], "charge": float(fields[3])}
    for residue in result.values():
        residue["summed_charge"] = sum(a["charge"] for a in residue["atoms"].values())
    return result


def psf(path):
    text = path.read_text()
    atom_match = re.search(r"^\s*(\d+)\s+!NATOM[^\n]*\n", text, re.M)
    atom_lines = text[atom_match.end():].strip().splitlines()[:int(atom_match[1])]
    atoms = {}
    for line in atom_lines:
        fields = line.split()
        atoms[int(fields[0])] = dict(segment=fields[1], resid=fields[2], residue=fields[3],
                                     name=fields[4], type=fields[5], charge=float(fields[6]))
    terms = {}
    for tag, width in [("NBOND", 2), ("NTHETA", 3), ("NPHI", 4), ("NIMPHI", 4)]:
        match = re.search(rf"^\s*(\d+)\s+!{tag}\b[^\n]*\n", text, re.M)
        count = int(match[1])
        tail = re.split(r"^.*!", text[match.end():], maxsplit=1, flags=re.M)[0]
        values = [int(v) for v in tail.split()][:count * width]
        assert len(values) == count * width
        terms[tag] = [values[i:i + width] for i in range(0, len(values), width)]
    return atoms, terms


def main():
    files = json.loads((HERE / "retrieval-dryad-browser.json").read_text())
    for f in files:
        assert hashlib.sha256((HERE / f["file"]).read_bytes()).hexdigest() == f["sha256"]
    folder = HERE / "dryad-oxy-myoglobin"
    old = residues(folder / "top_heme.inp")
    new = residues(folder / "top_heme_oxy_optimized_charges.inp")
    atoms, terms = psf(folder / "horse_oxy_only_wbi_10ns.psf")
    site = {i: a for i, a in atoms.items() if a["residue"] in {"HEM", "OXY"}}
    oxygen_ids = {i for i, a in site.items() if a["residue"] == "OXY"}
    target = next(c for c in json.loads((HERE / "target-sites.json").read_text()) if c["id"] == "1A6M")
    source = HERE.parent / target["source"]["path"]
    assert hashlib.sha256(source.read_bytes()).hexdigest() == target["source"]["sha256"]
    # Native CHARMM reader resolves the original parameter file; no System or Context is created.
    from openmm.app import CharmmParameterSet
    with warnings.catch_warnings(record=True) as parser_warnings:
        params = CharmmParameterSet(str(folder / "top_heme_oxy_optimized_charges.inp"),
                                   str(folder / "par_all27.inp"))
    bend = params.angle_types[("FE", "OM", "OM")]
    torsion = params.dihedral_types[("X", "FE", "OM", "X")]
    output = {
        "scope": "Static source inspection only; no new structure preparation or simulation",
        "source_hash_verified": target["source"]["sha256"],
        "downloaded_file_hashes_verified": len(files),
        "native_reader_warnings": [str(w.message) for w in parser_warnings],
        "published_topology": {"baseline": {k: old[k] for k in ["HEM", "OXY", "HSD", "TIP3"]},
                               "updated": {k: new[k] for k in ["HEM", "OXY", "HSD", "TIP3"]}},
        "published_psf_key_charges": [a for a in site.values() if a["name"] == "FE" or a["residue"] == "OXY"],
        "published_psf_uses_updated_charges": all(abs(a["charge"] - new[a["residue"]]["atoms"][a["name"]]["charge"]) < 1e-8 for a in site.values()),
        "published_psf_oxygen_terms": {k: [[dict(psf_id=i, **atoms[i]) for i in t] for t in v if oxygen_ids.intersection(t)] for k, v in terms.items()},
        "original_parameter_terms": {"O2_O1_Fe_bend": {"k_kcal_mol_rad2": bend.k, "equilibrium_degrees": bend.theteq},
                                     "X_Fe_OM_X_torsions": [{"k_kcal_mol": t.phi_k, "periodicity": t.per, "phase_degrees": t.phase} for t in torsion]},
        "sulfate_template_in_published_topology": "SO4" in new,
        "construction_script_loads": [s for s in (folder / "script_construction_oxy.pgn").read_text().splitlines() if s.startswith("topology ")],
        "component_inventory": target["retained_nonstandard_residues"],
        "validation_verdict": "Plausible published model candidate; not validated for 1A6M. Saved PSF charge mismatch and sulfate coverage remain explicit.",
    }
    assert output["published_psf_uses_updated_charges"] is False
    assert bend.k == 0
    assert all(t.phi_k == 0 and t.per > 0 for t in torsion)
    (HERE / "source-inspection.json").write_text(json.dumps(output, indent=2) + "\n")
    print(json.dumps({"inspection": "passed", "model_validated": False, "source_psf_charge_mismatch": True,
                      "zero_bend_preserved": True, "zero_torsion_preserved": True, "sulfate_template_present": "SO4" in new}))


if __name__ == "__main__":
    main()
