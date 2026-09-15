"""Read-only Ramachandran scoring of fresh default-path panel outputs.

Scores requested modeled residues plus fixed boundary residues whose phi/psi
include a modeled atom, with the same native CCTBX tables as the Sep 14 audit.
"""
import json, sys, glob
from pathlib import Path
import numpy as np
from openmm import app
sys.path.insert(0, str(Path(__file__).resolve().parents[3]))
from backend.loop_backbone_reference import rama_report, residue_key

label = sys.argv[1]
summary = {}
for result_path in sorted(glob.glob(f'docs/audit/loop-fallback/runs/{label}/*/result.json')):
    r = json.load(open(result_path))
    if not r['passed']:
        summary[r['id']] = {'accepted': False}
        continue
    target = Path(r['workspace']) / 'datasets' / r['prepared_dataset_id']
    pdb = app.PDBFile(str(target / 'prepared.pdb'))
    loop = r['preparation']['loop_construction']
    modeled = {tuple(a['identity'][:4]) for a in loop['modeler']['worker_report']['atoms']}
    residues = list(pdb.topology.residues())
    keyed = {residue_key(res): res for res in residues}
    # boundary: observed residues bonded (peptide) to a modeled residue
    prev, foll = {}, {}
    for a, b in pdb.topology.bonds():
        if a.residue is b.residue or {a.name, b.name} != {'C', 'N'}:
            continue
        c, n = (a, b) if a.name == 'C' else (b, a)
        foll[c.residue] = n.residue; prev[n.residue] = c.residue
    boundary = set()
    for key in modeled:
        res = keyed[key]
        for nb in (prev.get(res), foll.get(res)):
            if nb is not None and residue_key(nb) not in modeled:
                boundary.add(residue_key(nb))
    rep_modeled = rama_report(pdb.topology, pdb.positions, selected=sorted(modeled))
    rep_boundary = rama_report(pdb.topology, pdb.positions, selected=sorted(boundary))
    fmt = lambda rows: [f"{row['residue'][0]}:{row['residue'][1]} {row['residue'][3]} ({row['classification']}, phi {row['phi_degrees']:.0f}, psi {row['psi_degrees']:.0f})" for row in rows]
    summary[r['id']] = {
        'accepted': True, 'modeled_residues': len(modeled), 'boundary_residues': len(boundary),
        'modeled_scored': len(rep_modeled['rows']), 'modeled_unavailable': len(rep_modeled['unavailable']),
        'modeled_outliers': fmt(rep_modeled['outliers']),
        'modeled_allowed': sum(row['classification'] == 'allowed' for row in rep_modeled['rows']),
        'boundary_outliers': fmt(rep_boundary['outliers']),
        'cctbx_version': rep_modeled['cctbx_version'],
    }
print(json.dumps(summary, indent=1))
