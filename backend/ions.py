"""Named monatomic ion states with verified Amber TIP3P nonbonded parameters.

A recognized element is not sufficient: the deposited component identity,
optional formal charge, and the actual force-field particle must agree.  These
12-6 Lennard-Jones models do not represent bonded metal coordination, charge
transfer, ligand-field energetics, or changing oxidation states.
"""
from __future__ import annotations

import copy
from functools import lru_cache
import hashlib
import math
from pathlib import Path
from xml.etree import ElementTree

# Canonical component identities and formal charges were checked against the
# wwPDB Chemical Component Dictionary (RCSB core/chemcomp records, 2026-09-11).
# Aliases retain the existing DynaMol/CHARMM input spellings.  Do not derive a
# charge merely from an element or from whichever OpenMM template happens to
# match: for example, CCD AG is Ag+ but Amber's Ag template is Ag2+; Fe has
# multiple non-identical automatic template matches.
_CANONICAL_STATES = {
    "LI": ("Li", 1), "NA": ("Na", 1), "K": ("K", 1),
    "RB": ("Rb", 1), "CS": ("Cs", 1),
    "F": ("F", -1), "CL": ("Cl", -1), "BR": ("Br", -1), "IOD": ("I", -1),
    "MG": ("Mg", 2), "CA": ("Ca", 2),
    "ZN": ("Zn", 2), "MN": ("Mn", 2), "CO": ("Co", 2),
    "NI": ("Ni", 2), "CU": ("Cu", 2),
}
_ALIASES = {"SOD": "NA", "CLA": "CL", "POT": "K", "CAL": "CA"}
ION_STATES = {**_CANONICAL_STATES, **{alias: _CANONICAL_STATES[name] for alias, name in _ALIASES.items()}}
SUPPORTED_IONS = frozenset(ION_STATES)
ION_FORCEFIELD = "amber14/tip3p.xml"
_ION_ELEMENTS = frozenset({
    "Li", "Na", "K", "Rb", "Cs", "Mg", "Ca", "Sr", "Ba", "Zn", "Fe", "Mn", "Cu", "Co", "Ni", "Cd", "Hg", "Ag",
    "Al", "Be", "Ce", "Cr", "Dy", "Eu", "Er", "Gd", "Hf", "In", "La", "Lu", "Nd", "Pb", "Pd", "Pr", "Pt",
    "Pu", "Ra", "Sm", "Sn", "Tb", "Th", "Tl", "Tm", "U", "V", "Y", "Yb", "Zr", "F", "Cl", "Br", "I",
})
_MODEL_LIMITS = "Amber TIP3P 12-6 nonbonded ion; not a bonded coordination, charge-transfer or oxidation-state model."


@lru_cache(maxsize=1)
def _reference_parameters():
    """Read the installed, selected water/ion model rather than another XML set."""
    from openmm import app
    path = Path(app.__file__).parent / "data" / ION_FORCEFIELD
    raw = path.read_bytes()
    root = ElementTree.fromstring(raw)
    types = {item.get("name"): item.attrib for item in root.findall("AtomTypes/Type")}
    nonbonded = {item.get("type"): item.attrib for item in root.findall("NonbondedForce/Atom")}
    result = {}
    for residue in root.findall("Residues/Residue"):
        atoms = residue.findall("Atom")
        if len(atoms) != 1 or residue.findall("Bond") or residue.findall("ExternalBond"):
            continue
        atom = atoms[0]
        atom_type = atom.get("type")
        if atom_type not in types or atom_type not in nonbonded:
            continue
        result[residue.get("name")] = {
            "template": residue.get("name"), "atom_type": atom_type,
            "element": types[atom_type]["element"], "charge_e": float(atom.get("charge")),
            "sigma_nm": float(nonbonded[atom_type]["sigma"]),
            "epsilon_kj_mol": float(nonbonded[atom_type]["epsilon"]),
            "forcefield_file": ION_FORCEFIELD, "forcefield_sha256": hashlib.sha256(raw).hexdigest(),
            "parameter_sources": [item.text for item in root.findall("Info/Source")],
            "references": [item.text for item in root.findall("Info/Reference")],
        }
    return result


def _particle_parameters(system, index):
    from openmm import NonbondedForce, unit
    forces = [force for force in system.getForces() if isinstance(force, NonbondedForce)]
    if len(forces) != 1:
        raise ValueError("the selected ion model requires one standard NonbondedForce")
    charge, sigma, epsilon = forces[0].getParticleParameters(index)
    return {"charge_e": charge.value_in_unit(unit.elementary_charge),
            "sigma_nm": sigma.value_in_unit(unit.nanometer),
            "epsilon_kj_mol": epsilon.value_in_unit(unit.kilojoule_per_mole)}


def _require_matching_parameters(actual, reference):
    for key in ("charge_e", "sigma_nm", "epsilon_kj_mol"):
        value = actual[key]
        if not math.isfinite(value) or not math.isclose(value, reference[key], abs_tol=1e-10, rel_tol=1e-8):
            raise ValueError(f"the actual ion {key} does not match the declared Amber TIP3P model")
    if actual["sigma_nm"] <= 0 or actual["epsilon_kj_mol"] <= 0:
        raise ValueError("the declared ion has no positive finite Lennard-Jones parameters")


def _verify_model(name, forcefield):
    from openmm import app
    canonical = _ALIASES.get(name, name)
    expected = ION_STATES[name]
    reference = _reference_parameters().get(canonical)
    if reference is None or reference["element"] != expected[0] or reference["charge_e"] != expected[1]:
        raise ValueError("the installed named ion template does not match its declared component element/charge")
    top = app.Topology()
    residue = top.addResidue(name, top.addChain("I"), "1")
    top.addAtom(name, app.element.Element.getBySymbol(expected[0]), residue)
    matches = forcefield.getMatchingTemplates(top)
    if len(matches) != 1 or matches[0].name != canonical:
        raise ValueError("the selected force field does not uniquely select the declared ion template")
    system = forcefield.createSystem(top, nonbondedMethod=app.NoCutoff, removeCMMotion=False)
    if system.getNumParticles() != 1:
        raise ValueError("the selected force field did not produce a monatomic ion")
    actual = _particle_parameters(system, 0)
    _require_matching_parameters(actual, reference)
    return {**reference, "identity_source": f"https://www.rcsb.org/ligand/{canonical}",
            "validation": "Named template and actual OpenMM NonbondedForce particle charge/Lennard-Jones parameters agree."}


@lru_cache(maxsize=len(ION_STATES))
def _default_model(name):
    from openmm import app
    return _verify_model(name, app.ForceField("amber14/protein.ff14SB.xml", ION_FORCEFIELD))


def inspect_ions(topology, forcefield=None):
    """Inspect named states without altering atoms, charges or coordination.

    Supplying the already loaded force field checks that specific model.  The
    default is the protein/TIP3P force field selected by the DynaMol workflow.
    Supported means parameter/identity compatibility, not site-specific accuracy.
    """
    records = []
    for residue in topology.residues():
        atoms = list(residue.atoms())
        name = residue.name.upper()
        expected = ION_STATES.get(name)
        element = atoms[0].element.symbol if len(atoms) == 1 and atoms[0].element else None
        deposited_charge = getattr(atoms[0], "formalCharge", None) if len(atoms) == 1 else None
        if expected is None and not (len(atoms) == 1 and element in _ION_ELEMENTS):
            continue
        key = ":".join((residue.chain.id, residue.id, (residue.insertionCode or "").strip(), residue.name))
        error, parameters = None, None
        if expected is None:
            error = f"{key}: this named ion identity/state has no verified DynaMol preparation model. Supply a compatible specialized model; no oxidation state is guessed."
        elif len(atoms) != 1 or element != expected[0]:
            error = f"{key}: ion name requires one {expected[0]} atom, but the deposited atom count/element disagrees. Correct the chemical identity before preparation; no ion was changed."
        elif deposited_charge is not None and (type(deposited_charge) is not int or deposited_charge != expected[1]):
            error = f"{key}: deposited ion charge {deposited_charge!r} disagrees with the supported {expected[0]} state ({expected[1]:+d}). A different oxidation/charge state requires an appropriate model; no charge was changed."
        else:
            try:
                parameters = copy.deepcopy(_default_model(name) if forcefield is None else _verify_model(name, forcefield))
            except Exception as exc:
                error = f"{key}: ion force-field verification failed: {exc}. No replacement parameters were assigned."
        records.append({"key": key, "residue": residue.name, "element": element,
                        "deposited_formal_charge": deposited_charge,
                        "formal_charge": expected[1] if expected and error is None else None,
                        "supported": error is None, "error": error, "model": _MODEL_LIMITS,
                        "parameters": parameters})
    seen, duplicates = set(), set()
    for record in records:
        if record["key"] in seen:
            duplicates.add(record["key"])
        seen.add(record["key"])
    for record in records:
        if record["key"] in duplicates:
            record.update(supported=False, formal_charge=None, parameters=None,
                          error=f"{record['key']}: ion residue identifiers are not unique. Resolve the deposited identities before preparation.")
    return records


def validate_ion_system(topology, system, expected_ions=None):
    """Validate retained identities and actual prepared-system q/LJ particles.

    Return provenance suitable for ``preparation['ions']``.  Passing the input
    inspection records also verifies that no retained ion disappeared or changed
    identity during topology reconstruction.
    """
    records = inspect_ions(topology)
    for record in records:
        if not record["supported"]:
            raise ValueError(record["error"])
    if expected_ions is not None:
        if any(not record.get("supported") for record in expected_ions):
            raise ValueError("The original ion inventory contains an unsupported identity/state.")
        identity = lambda record: (record["key"], record["element"], record["formal_charge"])
        if sorted(map(identity, expected_ions)) != sorted(map(identity, records)):
            raise ValueError("The retained ion inventory changed during preparation; no missing or substituted ion was accepted.")
    if system.getNumParticles() != topology.getNumAtoms():
        raise ValueError("The prepared ion system and topology atom counts disagree.")
    residues = {":".join((residue.chain.id, residue.id, (residue.insertionCode or "").strip(), residue.name)): residue
                for residue in topology.residues()}
    for record in records:
        atom = next(residues[record["key"]].atoms())
        try:
            _require_matching_parameters(_particle_parameters(system, atom.index), record["parameters"])
        except ValueError as exc:
            raise ValueError(f"{record['key']}: {exc}; the prepared ion model was not accepted.") from exc
        record["parameters"]["prepared_particle_index"] = atom.index
        record["parameters"]["prepared_system_verified"] = True
    return records
