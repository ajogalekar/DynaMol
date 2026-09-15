"""Bounded native GAFF2 type/parameter coverage, never physical acceptance.

No charges are fitted. Generated MOL2 charge columns are zero placeholders and
must never be used for a simulation. All-GAFF2 caps here test chemical coverage;
the separate protein interface model still requires native ff14SB precedence.
"""
from __future__ import annotations
import argparse
from collections import Counter
import hashlib
import json
import os
from pathlib import Path
import re
import shutil
import subprocess
import time

import numpy as np
from rdkit import Chem

ROOT = Path(__file__).resolve().parents[3]
AMBER = ROOT / '.tools/ambertools'
SECTIONS = {'MASS': (1, 1), 'BOND': (2, 2), 'ANGLE': (3, 2),
            'DIHE': (4, 4), 'IMPROPER': (4, 3), 'NONBON': (1, 2)}


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def parse_frcmod(path):
    section = None
    entries = []
    for number, line in enumerate(Path(path).read_text().splitlines(), 1):
        if line.strip() in SECTIONS:
            section = line.strip()
            continue
        if not section or not line.strip():
            continue
        ntypes, nvalues = SECTIONS[section]
        if ntypes == 1:
            fields = line.split()
            atom_types, values = [fields[0]], [float(x) for x in fields[1:1+nvalues]]
        else:
            width = 3 * ntypes - 1
            atom_types = [x.strip() for x in line[:width].split('-')]
            values = [float(x) for x in line[width:].split()[:nvalues]]
        if len(atom_types) != ntypes or len(values) != nvalues or not np.isfinite(values).all():
            raise ValueError(f'Unparsed or nonfinite native parameter at line {number}')
        invalid = None
        if section in ('MASS', 'BOND') and any(x <= 0 for x in values):
            invalid = 'Nonpositive mass, bond stiffness or length'
        elif section == 'ANGLE' and (values[0] <= 0 or not 0 < values[1] <= 180):
            invalid = 'Invalid angle stiffness or equilibrium angle'
        elif section == 'DIHE' and (values[0] <= 0 or abs(values[3]) < 1 or not abs(values[3]).is_integer()):
            invalid = 'Invalid torsion divisor or periodicity; unmarked placeholders are not accepted'
        elif section == 'IMPROPER' and (abs(values[2]) < 1 or not abs(values[2]).is_integer()):
            invalid = 'Invalid improper periodicity'
        elif section == 'NONBON' and any(x <= 0 for x in values):
            invalid = 'Nonpositive Lennard-Jones radius or epsilon'
        score = re.search(r'penalty score\s*=\s*([\d.eE+\-]+)', line)
        entries.append({'line': number, 'section': section, 'atom_types': atom_types,
                        'native_values': values, 'raw': line,
                        'analogy_penalty': float(score.group(1)) if score else None,
                        'explicit_default': 'default' in line.lower(),
                        'attention': 'ATTN' in line, 'invalid': invalid})
    return {'entries': entries, 'counts_by_section': dict(Counter(x['section'] for x in entries)),
            'attention_terms': [x for x in entries if x['attention']],
            'invalid_terms': [x for x in entries if x['invalid']],
            'explicit_default_terms': [x for x in entries if x['explicit_default']],
            'analogy_terms': [x for x in entries if x['analogy_penalty'] is not None],
            'max_analogy_penalty': max((x['analogy_penalty'] for x in entries if x['analogy_penalty'] is not None), default=None),
            'interpretation': 'Penalty scores and defaults describe native assignment provenance, not a calibrated accuracy estimate. Absence of ATTN is not physical validation.'}


def verify_typing(source, typed):
    molecule = Chem.SDMolSupplier(str(source), sanitize=False, removeHs=False)[0]
    if molecule is None:
        raise ValueError('Source SDF cannot be parsed')
    sections = {}
    current = None
    for line in Path(typed).read_text().splitlines():
        if line.startswith('@<TRIPOS>'):
            current = line[len('@<TRIPOS>'):]
            sections[current] = []
        elif current and line.strip():
            sections[current].append(line.split())
    atoms = sections['ATOM']
    if len(atoms) != molecule.GetNumAtoms() or [int(x[0]) for x in atoms] != list(range(1, len(atoms)+1)):
        raise ValueError('Native typing changed atom count or index order')
    symbols = [re.sub(r'\d+$', '', x[1]).title() for x in atoms]
    if symbols != [x.GetSymbol() for x in molecule.GetAtoms()]:
        raise ValueError('Native typing changed element/index identity')
    xyz = np.array([[float(x) for x in row[2:5]] for row in atoms])
    max_component = float(np.max(np.abs(xyz - molecule.GetConformer().GetPositions())))
    # Native atomtype reads the intermediate AC format, which serializes three
    # decimals even though final MOL2 prints four. Original SDF remains intact.
    if max_component > 0.000501:
        raise ValueError('Native typing changed coordinates beyond intermediate AC three-decimal serialization')
    heavy = np.array([x.GetAtomicNum() != 1 for x in molecule.GetAtoms()])
    source_bonds = {tuple(sorted([b.GetBeginAtomIdx()+1, b.GetEndAtomIdx()+1])): b.GetBondTypeAsDouble()
                    for b in molecule.GetBonds()}
    native_bonds = {}
    for row in sections['BOND']:
        pair = tuple(sorted([int(row[1]), int(row[2])]))
        if pair in native_bonds:
            raise ValueError('Native typing duplicated a bond')
        native_bonds[pair] = float(row[3])
    if native_bonds != source_bonds:
        raise ValueError('Native typing changed original explicit connectivity or heavy bond orders')
    if any(not row[5] or row[5] in {'du', 'DU', 'unknown'} for row in atoms):
        raise ValueError('Native atom typing returned an unknown type')
    charges = np.array([float(row[8]) for row in atoms])
    if not np.isfinite(charges).all() or np.any(charges != 0):
        raise ValueError('Typing-only output unexpectedly contains nonzero charges')
    return {'atom_count': len(atoms), 'bond_count': len(native_bonds),
            'element_and_index_order_preserved': True, 'explicit_bonds_and_orders_preserved': True,
            'maximum_coordinate_component_serialization_error_angstrom': max_component,
            'maximum_heavy_coordinate_component_serialization_error_angstrom': float(np.max(np.abs(xyz[heavy] - molecule.GetConformer().GetPositions()[heavy]))),
            'native_coordinate_format': 'AC intermediate has three decimals; original SDF coordinates retained separately.',
            'gaff2_type_counts': dict(Counter(row[5] for row in atoms)),
            'charge_columns_are_unusable_zero_placeholders': True,
            'atom_types': [{'index': i, 'element': symbols[i], 'native_type': row[5]} for i, row in enumerate(atoms)]}


def execute(command, folder, log):
    env = dict(os.environ, AMBERHOME=str(AMBER))
    env['PATH'] = str(AMBER/'bin') + os.pathsep + env.get('PATH', '')
    began = time.monotonic()
    with (folder/log).open('w') as stream:
        done = subprocess.run(command, cwd=folder, env=env, stdout=stream,
                              stderr=subprocess.STDOUT, timeout=120)
    record = {'command': command, 'elapsed_seconds': time.monotonic()-began,
              'exit_code': done.returncode, 'log': log, 'timeout_seconds': 120}
    if done.returncode:
        raise ValueError(f'Native coverage tool failed; inspect {log}')
    return record


def run(graph_folder, output):
    graph_folder = graph_folder.resolve()
    output = output.resolve()
    output.mkdir(parents=True, exist_ok=False)
    summary_path = graph_folder/'summary.json'
    graphs = json.loads(summary_path.read_text())
    result = {'scope': 'Native GAFF2 typing/parameter coverage only. No charge fitting, no QM, no dynamics, no prepared complex or physical acceptance.',
              'full_preparation_accepted': False, 'chemical_state_accepted': False,
              'script_sha256': sha(__file__), 'graph_summary_sha256': sha(summary_path),
              'panel_sha256': graphs['panel_sha256'], 'original_case_count': graphs['case_count'],
              'runtime_sha256': {str(path.relative_to(AMBER)): sha(path) for path in
                [AMBER/'bin/antechamber', AMBER/'bin/atomtype', AMBER/'bin/parmchk2',
                 AMBER/'dat/antechamber/ATOMTYPE_GFF2.DEF', AMBER/'dat/antechamber/PARMCHK.DAT',
                 AMBER/'dat/leap/parm/gaff2.dat']}, 'cases': [], 'unique_chemistries': []}
    seen = {}
    for case in graphs['cases']:
        row = {'id': case['id'], 'adduct_class': case['adduct_class'], 'connections': []}
        for connection in case['connections']:
            name = connection['name']
            if connection['status'] != 'graph_candidate_constructed':
                row['connections'].append({'name': name, 'status': 'graph_blocked', 'reason': connection['reason']})
                continue
            graph = json.loads((graph_folder/case['id']/(name+'.json')).read_text())
            chemistry = graph['mapped_product_smiles']
            key = hashlib.sha256(chemistry.encode()).hexdigest()[:16]
            row['connections'].append({'name': name, 'chemistry_key': key})
            if chemistry in seen:
                seen[chemistry]['represented_connections'].append({'case': case['id'], 'connection': name})
                continue
            folder = output/(case['id']+'-'+name)
            folder.mkdir()
            source = graph_folder/case['id']/(name+'.sdf')
            shutil.copyfile(source, folder/'input.sdf')
            item = {'chemistry_key': key, 'folder': str(folder.relative_to(output)),
                    'represented_connections': [{'case': case['id'], 'connection': name}],
                    'canonical_isomeric_smiles': chemistry, 'source_sdf_sha256': sha(source),
                    'graph_report_sha256': sha(graph_folder/case['id']/(name+'.json')),
                    'candidate_formal_charge': graph['formal_charge_from_ccd_candidate'],
                    'protonation_assumption': 'Regenerated hydrogens at deposited CCD formal-charge state; no pH-state or tautomer validation.',
                    'full_preparation_accepted': False, 'usable_parameters': False, 'commands': []}
            seen[chemistry] = item
            try:
                item['commands'].append(execute([str(AMBER/'bin/antechamber'), '-i', 'input.sdf', '-fi', 'sdf',
                    '-o', 'typed.mol2', '-fo', 'mol2', '-at', 'gaff2', '-j', '1', '-nc', str(item['candidate_formal_charge']),
                    '-s', '2', '-seq', 'n', '-dr', 'yes', '-pf', 'n'], folder, 'antechamber.log'))
                item['typing_verification'] = verify_typing(folder/'input.sdf', folder/'typed.mol2')
                item['commands'].append(execute([str(AMBER/'bin/parmchk2'), '-i', 'typed.mol2', '-f', 'mol2',
                    '-o', 'coverage.frcmod', '-s', 'gaff2', '-a', 'Y'], folder, 'parmchk2.log'))
                coverage = parse_frcmod(folder/'coverage.frcmod')
                (folder/'parameter-assignments.json').write_text(json.dumps(coverage, indent=2)+'\n')
                item['parameter_assignment_summary'] = {k: v for k,v in coverage.items() if k not in ['entries', 'analogy_terms', 'explicit_default_terms']}
                item['analogy_assignment_count'] = len(coverage['analogy_terms'])
                item['explicit_default_count'] = len(coverage['explicit_default_terms'])
                item['native_parameter_coverage_complete'] = not coverage['attention_terms'] and not coverage['invalid_terms']
                item['status'] = 'coverage_screened'
                item['typed_mol2_sha256'] = sha(folder/'typed.mol2')
                item['coverage_frcmod_sha256'] = sha(folder/'coverage.frcmod')
            except Exception as error:
                item.update(status='coverage_blocked', reason=str(error), error_type=type(error).__name__, native_parameter_coverage_complete=False)
            result['unique_chemistries'].append(item)
            (folder/'result.json').write_text(json.dumps(item, indent=2)+'\n')
        result['cases'].append(row)
        (output/'progress.json').write_text(json.dumps({'completed_cases': len(result['cases']), 'total_cases': len(graphs['cases']), 'full_preparation_accepted': False}, indent=2)+'\n')
    for item in result['unique_chemistries']:
        (output/item['folder']/'result.json').write_text(json.dumps(item, indent=2)+'\n')
    result['unique_chemistry_count'] = len(seen)
    result['unique_chemistries_with_native_parameter_coverage'] = sum(x.get('native_parameter_coverage_complete', False) for x in seen.values())
    (output/'summary.json').write_text(json.dumps(result, indent=2)+'\n')
    return {k: v for k,v in result.items() if k not in ['cases', 'unique_chemistries', 'runtime_sha256']}


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('graphs', type=Path)
    parser.add_argument('output', type=Path)
    args = parser.parse_args()
    print(json.dumps(run(args.graphs, args.output), indent=2))
