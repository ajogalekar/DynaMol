"""Exact mmCIF polymer positions, independent of gaps in author numbering."""
from __future__ import annotations

from collections import defaultdict

_ALIASES = {"HID": "HIS", "HIE": "HIS", "HIP": "HIS", "CYX": "CYS", "CYM": "CYS", "ASH": "ASP", "GLH": "GLU", "LYN": "LYS"}


def _text(value):
    return "" if value in {None, False, "", ".", "?"} else str(value).strip()


def load_scheme(path, label_to_author, original_to_canonical, chain_ids, sequences):
    """Read explicit sequence identities; no sequence/numbering interpolation."""
    import gemmi
    table = gemmi.cif.read_file(str(path)).sole_block().get_mmcif_category("_pdbx_poly_seq_scheme.")
    if not table:
        return {}
    groups = defaultdict(list)
    labels = defaultdict(set)
    for i, label in enumerate(table.get("asym_id", [])):
        author = _text(table.get("pdb_strand_id", [None] * len(table["asym_id"]))[i]) or label_to_author.get(label, label)
        chain = original_to_canonical.get(author, author)
        if chain not in chain_ids:
            continue
        if "," in author:
            raise ValueError("The polymer scheme has an ambiguous author chain ID")
        resid = _text(table.get("auth_seq_num", [None] * len(table["asym_id"]))[i]) or _text(table.get("pdb_seq_num", [None] * len(table["asym_id"]))[i])
        if not resid:
            raise ValueError(f"Polymer sequence position {table['seq_id'][i]} has no explicit author/PDB residue ID")
        labels[chain].add(label)
        groups[chain].append({"seq_id": int(table["seq_id"][i]), "resid": resid,
                              "insertion_code": _text(table.get("pdb_ins_code", [None] * len(table["asym_id"]))[i]),
                              "residue": table["mon_id"][i], "label_asym_id": label})
    for sequence in sequences:
        chain = sequence.chainId
        rows = sorted(groups.get(chain, []), key=lambda row: row["seq_id"])
        if len(labels[chain]) != 1 or [r["seq_id"] for r in rows] != list(range(1, len(sequence.residues) + 1)):
            raise ValueError(f"Polymer scheme for chain {chain} is absent, incomplete, or ambiguous")
        if [r["residue"] for r in rows] != sequence.residues:
            raise ValueError(f"Polymer scheme and declared sequence disagree for chain {chain}; residue identities will not be guessed")
        if len({(r["resid"], r["insertion_code"]) for r in rows}) != len(rows):
            raise ValueError(f"Polymer scheme has duplicate residue identities in chain {chain}")
        groups[chain] = rows
    return dict(groups)


def recover_observed_insertion_codes(fixer):
    """Recover omitted codes only from unique number/name matches and chain order.

    Validate the whole mapping before changing any in-memory identity. Repeated
    names at the same author number are ambiguous, even if an ordering-based
    guess could assign them. Source files and coordinates are never changed.
    """
    from .residue_identity import protein_residue_keys, residue_key
    if not getattr(fixer, "sequence_scheme", None):
        return []
    protein_keys = protein_residue_keys(None, fixer.topology, fixer.positions)
    assignments, visited = [], set()
    for chain in fixer.topology.chains():
        if chain.id not in fixer.sequence_scheme:
            continue
        residues = [r for r in chain.residues() if residue_key(r) in protein_keys]
        if not residues:
            continue
        if chain.id in visited:
            raise ValueError(f"Multiple protein chains share scheme identity {chain.id}")
        visited.add(chain.id)
        rows = fixer.sequence_scheme[chain.id]
        previous = 0
        for residue in residues:
            code = _text(residue.insertionCode)
            candidates = [row for row in rows if row['resid'] == residue.id
                          and _ALIASES.get(row['residue'], row['residue']) == _ALIASES.get(residue.name, residue.name)
                          and (not code or row['insertion_code'] == code)]
            if len(candidates) != 1:
                raise ValueError(f"Observed residue {chain.id}:{residue.id}{code} {residue.name} has no unique number/name/insertion-code match in the polymer scheme")
            row = candidates[0]
            if row['seq_id'] <= previous:
                raise ValueError(f"Observed residues in chain {chain.id} do not follow the polymer scheme uniquely")
            previous = row['seq_id']
            if not code and row['insertion_code']:
                assignments.append((residue, row))
    records = []
    for residue, row in assignments:
        records.append({'chain': residue.chain.id, 'resid': residue.id, 'residue': residue.name,
                        'previous_insertion_code': _text(residue.insertionCode), 'insertion_code': row['insertion_code'],
                        'seq_id': row['seq_id'], 'label_asym_id': row['label_asym_id'],
                        'method': 'Unique source residue number/name match with full observed chain order verified'})
        residue.insertionCode = row['insertion_code']
    return records


def missing_from_scheme(fixer):
    """Produce PDBFixer's insertion-index contract with exact source identities."""
    from .residue_identity import protein_residue_keys, residue_key
    protein_keys = protein_residue_keys(None, fixer.topology, fixer.positions)
    missing, identities, visited = {}, {}, set()
    for chain in fixer.topology.chains():
        if chain.id not in fixer.sequence_scheme:
            continue
        residues = list(chain.residues())
        observed = [(i, residue) for i, residue in enumerate(residues) if residue_key(residue) in protein_keys]
        if not observed:
            continue
        if chain.id in visited:
            raise ValueError(f"Multiple protein chains share scheme identity {chain.id}")
        visited.add(chain.id)
        rows = fixer.sequence_scheme[chain.id]
        lookup = {(row["resid"], row["insertion_code"]): row for row in rows}
        positions = []
        for index, residue in observed:
            row = lookup.get((residue.id, _text(residue.insertionCode)))
            if row is None or _ALIASES.get(row["residue"], row["residue"]) != _ALIASES.get(residue.name, residue.name):
                raise ValueError(f"Observed residue {chain.id}:{residue.id}{_text(residue.insertionCode)} {residue.name} does not match an exact polymer-scheme identity")
            positions.append((row["seq_id"], index))
        if any(a[0] >= b[0] for a, b in zip(positions, positions[1:])):
            raise ValueError(f"Observed residues in chain {chain.id} do not follow the polymer scheme uniquely")
        retained = {seq_id for seq_id, _ in positions}
        following = 0
        for row in rows:
            if row["seq_id"] in retained:
                continue
            while following < len(positions) and positions[following][0] < row["seq_id"]:
                following += 1
            insertion = positions[following][1] if following < len(positions) else len(residues)
            # Missing N-terminal sequence is always terminal, even if other
            # retained components precede the first protein residue.
            if row["seq_id"] < positions[0][0]:
                insertion = 0
            key = chain.index, insertion
            missing.setdefault(key, []).append(row["residue"])
            identities.setdefault(key, []).append(dict(row))
    if set(fixer.sequence_scheme) - visited:
        raise ValueError("A retained polymer scheme has no unambiguous observed protein chain")
    fixer.missingResidues = missing
    fixer.missingResidueIdentities = identities


def restoration_plan(fixer, selected_loops):
    """Capture exact retained and inserted IDs before PDBFixer replaces topology."""
    if not getattr(fixer, "sequence_scheme", None):
        return None
    plans = []
    for chain in fixer.topology.chains():
        residues = list(chain.residues())
        rows = []
        for index in range(len(residues) + 1):
            key = chain.index, index
            if key in selected_loops:
                additions = fixer.missingResidueIdentities.get(key)
                if not additions or [row["residue"] for row in additions] != selected_loops[key]:
                    raise ValueError("Selected loop lacks exact polymer-scheme identities")
                rows.extend({**row, "added": True} for row in additions)
            if index < len(residues):
                residue = residues[index]
                rows.append({"residue": residue.name, "resid": residue.id, "insertion_code": residue.insertionCode or "", "added": False})
        plans.append({"chain": chain.id, "residues": rows})
    return plans


def restore_residue_identities(topology, plans):
    """Validate all topology identities before assigning source IDs; no coordinates."""
    if plans is None:
        return []
    chains = list(topology.chains())
    if len(chains) != len(plans):
        raise ValueError("Loop construction changed the recorded chain count")
    assignments = []
    for chain, plan in zip(chains, plans):
        residues = list(chain.residues())
        if chain.id != plan["chain"] or [r.name for r in residues] != [r["residue"] for r in plan["residues"]]:
            raise ValueError("Loop construction changed the recorded residue identity/order")
        assignments.extend((residue, row) for residue, row in zip(residues, plan["residues"]))
    added = []
    for residue, row in assignments:
        residue.id, residue.insertionCode = row["resid"], row["insertion_code"]
        if row["added"]:
            added.append({"chain": residue.chain.id, **{key: row[key] for key in ("seq_id", "resid", "insertion_code", "residue", "label_asym_id")}})
    return added
