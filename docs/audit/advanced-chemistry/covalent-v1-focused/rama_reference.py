"""Research-only native CCTBX Ramachandran reference; no fitted cutoffs."""
import argparse
import importlib.metadata
import json
from pathlib import Path

import numpy as np
from openmm import app, unit
import mmtbx_validation_ramachandran_ext as rama_ext

from prepare_adduct import digest

STANDARD = set('ALA ARG ASN ASP CYS GLN GLU GLY HIS ILE LEU LYS MET PHE PRO SER THR TRP TYR VAL HID HIE HIP'.split())
CLASSES = ['general', 'glycine', 'cis-proline', 'trans-proline', 'pre-proline', 'isoleucine or valine']


def residue_key(residue):
    return (residue.chain.id, residue.id, (residue.insertionCode or '').strip(), residue.name)


def dihedral64(points):
    p0, p1, p2, p3 = np.asarray(points, dtype=np.float64)
    axis = p2-p1
    if np.linalg.norm(axis) < 1e-12:
        raise ValueError('Zero-length central bond in backbone dihedral')
    axis = axis/np.linalg.norm(axis)
    left, right = p0-p1, p3-p2
    left = left-axis*np.dot(left, axis)
    right = right-axis*np.dot(right, axis)
    if min(np.linalg.norm(left), np.linalg.norm(right)) < 1e-12:
        raise ValueError('Collinear backbone atoms cannot receive a torsion score')
    return float(np.degrees(np.arctan2(np.dot(np.cross(axis, left), right), np.dot(left, right))))


def rama_report(topology, positions, selected=None):
    xyz = np.asarray(positions.value_in_unit(unit.nanometer))
    if not np.isfinite(xyz).all():
        raise ValueError('Nonfinite coordinates cannot receive a Ramachandran score')
    residues = list(topology.residues())
    names = {r: {a.name: a for a in r.atoms()} for r in residues}
    previous, following = {}, {}
    for a, b in topology.bonds():
        if a.residue is b.residue or {a.name, b.name} != {'C', 'N'}:
            continue
        c, n = (a, b) if a.name == 'C' else (b, a)
        if c.residue in following or n.residue in previous:
            raise ValueError('Ambiguous peptide connectivity in reference scoring')
        following[c.residue] = n.residue
        previous[n.residue] = c.residue
    requested = set(map(tuple, selected)) if selected is not None else None
    found = set()
    rows, unavailable = [], []
    evaluator = rama_ext.rama_eval()
    for res in residues:
        key = residue_key(res)
        if requested is not None and key not in requested:
            continue
        found.add(key)
        if res.name not in STANDARD or res not in previous or res not in following:
            unavailable.append({'residue': key, 'reason': 'Nonstandard residue or incomplete peptide neighborhood'})
            continue
        before, after = previous[res], following[res]
        try:
            p = [names[before]['C'], names[res]['N'], names[res]['CA'], names[res]['C']]
            q = [names[res]['N'], names[res]['CA'], names[res]['C'], names[after]['N']]
            w = [names[before]['CA'], names[before]['C'], names[res]['N'], names[res]['CA']]
        except KeyError:
            unavailable.append({'residue': key, 'reason': 'Missing defining backbone atom'})
            continue
        phi, psi, omega = [dihedral64(xyz[[a.index for a in atoms]]) for atoms in [p, q, w]]
        if not np.isfinite([phi, psi, omega]).all():
            raise ValueError('Undefined backbone dihedral')
        cls = 1 if res.name == 'GLY' else (2 if -90 < omega < 90 else 3) if res.name == 'PRO' else 4 if after.name == 'PRO' else 5 if res.name in {'ILE', 'VAL'} else 0
        score = evaluator.get_score(cls, float(phi), float(psi))
        classification = evaluator.evaluate_score(cls, score)
        rows.append({'residue': key, 'phi_degrees': float(phi), 'psi_degrees': float(psi),
                     'preceding_omega_degrees': float(omega), 'rama_class': CLASSES[cls],
                     'score': score, 'classification': ['outlier', 'allowed', 'favored'][classification],
                     'phi_atoms': [a.index for a in p], 'psi_atoms': [a.index for a in q]})
    if requested is not None and requested != found:
        raise ValueError('Requested residue identities are missing or ambiguous')
    library = Path(rama_ext.__file__)
    return {'method': 'Native CCTBX rama_eval tables, residue classes and classification thresholds; double-precision vector-projection angles.',
            'cctbx_version': importlib.metadata.version('cctbx-base'),
            'native_library': str(library), 'native_library_sha256': digest(library),
            'rows': rows, 'unavailable': unavailable,
            'outliers': [row for row in rows if row['classification'] == 'outlier'],
            'all_selected_scored_without_outliers': bool(rows) and not unavailable and all(row['classification'] != 'outlier' for row in rows),
            'scope': 'Static structural reference quality, not proof of a native loop or an MD per-frame acceptance rule.'}


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--pdb', type=Path, required=True)
    parser.add_argument('--selection-result', type=Path)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    selected = None
    hashes = {str(args.pdb): digest(args.pdb)}
    if args.selection_result is not None:
        data = json.loads(args.selection_result.read_text())
        selected = data['modeled_residues']+data['native_changed_context_residues']
        hashes[str(args.selection_result)] = digest(args.selection_result)
    structure = app.PDBFile(str(args.pdb))
    report = rama_report(structure.topology, structure.positions, selected)
    report.update(stage='native_ramachandran_reference_audit', source_sha256=hashes,
                  app_ready=False, physical_model_validated=False)
    if hashes != {name: digest(Path(name)) for name in hashes}:
        raise ValueError('Source changed during reference audit')
    args.output.mkdir(exist_ok=False)
    (args.output/'result.json').write_text(json.dumps(report, indent=2)+'\n')
    print(json.dumps({k: v for k, v in report.items() if k not in ['rows', 'source_sha256']}, indent=2))
