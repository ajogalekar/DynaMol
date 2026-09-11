"""Named, monatomic ion states supported by the explicit TIP3P workflow."""
ION_STATES = {
    "NA": ("Na", 1), "SOD": ("Na", 1),
    "CL": ("Cl", -1), "CLA": ("Cl", -1),
    "K": ("K", 1), "POT": ("K", 1),
    "MG": ("Mg", 2), "CA": ("Ca", 2), "CAL": ("Ca", 2),
}
SUPPORTED_IONS = frozenset(ION_STATES)


def inspect_ions(topology):
    """Check atom count and element, rather than trusting residue names alone."""
    records = []
    metals = {"Na", "K", "Mg", "Ca", "Zn", "Fe", "Mn", "Cu", "Co", "Ni", "Cd", "Hg"}
    for residue in topology.residues():
        atoms = list(residue.atoms())
        expected = ION_STATES.get(residue.name.upper())
        element = atoms[0].element.symbol if len(atoms) == 1 and atoms[0].element else None
        deposited_charge = getattr(atoms[0], "formalCharge", None) if len(atoms) == 1 else None
        if expected is None and not (len(atoms) == 1 and element in metals | {"Cl"}):
            continue
        key = ":".join((residue.chain.id, residue.id, (residue.insertionCode or "").strip(), residue.name))
        error = None
        if expected is None:
            error = f"{key}: this ion identity/state has no validated DynaMol preparation model. Supply a compatible specialized model; no oxidation state is guessed."
        elif len(atoms) != 1 or element != expected[0]:
            error = f"{key}: ion name requires one {expected[0]} atom, but the deposited atom count/element disagrees. Correct the chemical identity before preparation; no ion was changed."
        elif deposited_charge is not None and deposited_charge != expected[1]:
            error = f"{key}: deposited ion charge {deposited_charge:+d} disagrees with the supported {expected[0]} state ({expected[1]:+d}). A different oxidation/charge state requires an appropriate model; no charge was changed."
        records.append({"key": key, "residue": residue.name, "element": element,
                        "deposited_formal_charge": deposited_charge,
                        "formal_charge": expected[1] if expected and error is None else None,
                        "supported": error is None, "error": error,
                        "model": "TIP3P-compatible Amber nonbonded ion; not a bonded coordination or oxidation-state model."})
    return records
