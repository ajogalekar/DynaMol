"""Export a backbone-only seed from preserved native missing+context atoms."""
import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
from openmm import app, unit

from screen_kic_backbones import geometry


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--input', type=Path, required=True)
    parser.add_argument('--worker-output', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    args.output.mkdir(exist_ok=False)
    hashes = {str(p):sha(p) for p in [args.input, args.worker_output, Path(__file__)]}
    data = json.loads(args.worker_output.read_text())
    config = json.loads(args.input.read_text())
    if data['status'] != 'candidate':
        raise ValueError('Native worker did not produce a candidate')
    native_rows = data['atoms'] + data['research_context_atoms']
    lookup = {tuple(r['identity']):r['xyz_nm'] for r in native_rows}
    if len(lookup) != len(native_rows):
        raise ValueError('Native missing/context identities overlap or duplicate')
    topology = app.Topology()
    positions, rows = [], []
    for chainrow in config['chains']:
        chain = topology.addChain(chainrow['chain_id'])
        previous = None
        for record in chainrow['residues']:
            identity = [record['chain'], record['resid'], record['insertion_code'], record['residue']]
            atom_names = ['N', 'CA', 'C', 'O']
            if not any(tuple(identity+[a]) in lookup for a in atom_names):
                continue
            if not all(tuple(identity+[a]) in lookup for a in atom_names):
                raise ValueError('Native backbone incomplete for ' + str(identity))
            residue = topology.addResidue(record['residue'], chain, record['resid'], record['insertion_code'])
            atoms = {}
            for name in atom_names:
                element = app.element.nitrogen if name == 'N' else app.element.oxygen if name == 'O' else app.element.carbon
                atoms[name] = topology.addAtom(name, element, residue)
                positions.append(lookup[tuple(identity+[name])])
                rows.append({'identity':identity+[name], 'xyz_angstrom':(np.asarray(positions[-1])*10).tolist()})
            for first, last in [('N', 'CA'), ('CA', 'C'), ('C', 'O')]:
                topology.addBond(atoms[first], atoms[last])
            if previous is not None:
                topology.addBond(previous, atoms['N'])
            previous = atoms['C']
    seed = args.output/'coherent-native-BACKBONE-seed.pdb'
    with seed.open('w') as stream:
        app.PDBFile.writeFile(topology, np.asarray(positions)*unit.nanometer, stream, keepIds=True)
    local = {f'{first}-{last}':geometry([r for r in rows if r['identity'][0]=='A' and first<=int(r['identity'][1])<=last])
             for first,last in [(176,180),(175,181)]}
    if hashes != {name:sha(Path(name)) for name in hashes}:
        raise ValueError('Export input changed')
    (args.output/'implementation.py').write_bytes(Path(__file__).read_bytes())
    (args.output/'result.json').write_text(json.dumps({'source_sha256':hashes, 'seed_sha256':sha(seed),
        'backbone_atom_count':len(positions), 'local_geometry':local,
        'all_local_backbone_geometry_pass':all(g['accepted'] for g in local.values()),
        'observed_coordinates_restored_before_seed_export':False,
        'scope':'Native missing and observed backbone coordinates exported together for research. No full-atom or app admission.'},indent=2)+'\n')


if __name__ == '__main__':
    main()
