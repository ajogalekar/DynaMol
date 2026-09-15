"""Check mapped loop geometry in a native complete-complex topology.

No construction forces are added. The same geometry bounds are used before
MD and at sampled frames; this remains a gross screen, not native-loop proof.
"""
import argparse
import json
from pathlib import Path
import sys

import numpy as np
from openmm import app
from pdbfixer import PDBFixer

from prepare_adduct import digest


class LoopMonitor:
    def __init__(self, folder, protein, topology, record_peptide_excursions=False):
        sys.path.insert(0, str(Path(__file__).resolve().parents[4]))
        from backend.residue_identity import residue_key
        self.topology = topology
        self.protein = protein
        self.record_peptide_excursions = record_peptide_excursions
        self.recorded_excursions = []
        self.initial_peptide_states = None
        prep = json.loads((protein/'result.json').read_text())
        assembled = json.loads((folder/'result.json').read_text())
        if not prep.get('retained_environment_check', {}).get('passed'):
            raise ValueError('This monitored loop test requires a passing retained-environment audit')
        for name in ['protein.pdb', 'retained-environment.pdb', 'result.json']:
            if assembled['source_sha256'].get(str(protein/name)) != digest(protein/name):
                raise ValueError('Loop protein intermediate is not bound to this native assembly')
        for name in ['solvated.prmtop', 'source-mapping.json']:
            if assembled['outputs_sha256'][name] != digest(folder/name):
                raise ValueError('Native assembly changed before loop monitoring')
        wanted = {tuple(k) for k in prep['modeled_residues'] + prep.get('native_changed_context_residues', [])}
        atoms = list(topology.atoms())
        mapped = {key: set() for key in wanted}
        mapping = json.loads((folder/'source-mapping.json').read_text())
        for row in mapping:
            key = tuple(row.get('source', {}).get('residue', []))
            if key not in mapped:
                continue
            atom = atoms[row['native_index']]
            if atom.name != row['source']['atom']:
                raise ValueError('Loop native mapping changed an atom name')
            mapped[key].add(residue_key(atom.residue))
        if any(len(matches) != 1 for matches in mapped.values()):
            raise ValueError('Every source loop/context residue must map to one native residue')
        self.keys = {next(iter(matches)) for matches in mapped.values()}
        if len(self.keys) != len(wanted):
            raise ValueError('Different source residues map to the same native residue')
        self.remodeled_atom_indices = {atom.index for atom in atoms if residue_key(atom.residue) in self.keys}
        self.source_mapping = [{'source': list(key), 'native': list(next(iter(matches)))}
                               for key, matches in sorted(mapped.items())]
        self.templates = PDBFixer(filename=str(protein/'protein.pdb')).templates

    def check(self, positions, output, stage, time_ps):
        from backend.loop_geometry import loop_geometry_report
        from backend.preparation_worker import stereochemistry_report, require_valid_stereochemistry
        output.mkdir(exist_ok=True)
        geometry = loop_geometry_report(self.topology, positions, self.keys)
        stereo = stereochemistry_report(self.topology, positions, self.templates)
        error = None
        try:
            require_valid_stereochemistry(stereo, 'Full-complex loop test')
        except Exception as exc:
            error = str(exc)
        states = {tuple(row['atoms']): row['state'] for row in geometry['omega_checks']}
        if self.initial_peptide_states is None:
            if time_ps != 0 or not geometry['accepted'] or error is not None:
                raise ValueError('Loop monitor needs a passing initial geometry/stereo frame')
            self.initial_peptide_states = states
        basin_changes = [{'atoms': list(key), 'initial': self.initial_peptide_states.get(key), 'current': state}
                         for key, state in states.items() if self.initial_peptide_states.get(key) != state]
        failed_omega = [row for row in geometry['omega_checks'] if not row['accepted']]
        expected_omega_errors = {
            f"Peptide omega before {tuple(row['residues'][1])} is undefined or more than 35° from cis/trans."
            for row in failed_omega}
        # Observation mode never changes a failed geometry result into a pass.
        # It only permits additional diagnostic sampling for finite peptide
        # excursions, with no other errors or change of cis/trans basin.
        diagnostic_continuation = bool(self.record_peptide_excursions and time_ps > 0
            and error is None and not basin_changes and failed_omega
            and set(geometry['errors']) == expected_omega_errors
            and all(row['deviation_degrees'] is not None and np.isfinite(row['deviation_degrees']) for row in failed_omega))
        report = {'stage': stage, 'time_ps': time_ps,
            'passed': geometry['accepted'] and error is None and not basin_changes,
            'app_ready': False, 'physical_model_validated': False,
            'source_to_native_residues': self.source_mapping, 'geometry': geometry,
            'stereochemistry': stereo, 'stereo_error': error,
            'peptide_basin_changes': basin_changes,
            'diagnostic_continuation_permitted': diagnostic_continuation,
            'scope': 'Whole-coordinate, nonperiodic local geometry screen in the full native complex, plus standard protein stereochemistry. Temporary construction forces are absent; native loop conformation and equilibrium motion are not established.'}
        filename = f'{time_ps:03d}-'+stage.replace(' ', '-')+'.json'
        (output/filename).write_text(json.dumps(report, indent=2)+'\n')
        if diagnostic_continuation:
            self.recorded_excursions.append({'stage': stage, 'time_ps': time_ps, 'report': filename,
                                            'omega': failed_omega, 'geometry_accepted': False})
            (output/'recorded-peptide-excursions.json').write_text(json.dumps(self.recorded_excursions, indent=2)+'\n')
        if not report['passed'] and not diagnostic_continuation:
            raise ValueError('Full-complex loop geometry/stereochemistry failed: '+str(output/filename))
        return {'loop_geometry_accepted': geometry['accepted'], 'loop_stereochemistry_accepted': error is None,
                'loop_peptide_excursion_recorded': diagnostic_continuation,
                'loop_soft_overlap_count': len(geometry['soft_overlaps'])}


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    for name in ['input', 'protein', 'output']:
        parser.add_argument('--'+name, type=Path, required=True)
    args = parser.parse_args()
    args.output.mkdir(exist_ok=False)
    topology = app.AmberPrmtopFile(str(args.input/'solvated.prmtop')).topology
    positions = app.AmberInpcrdFile(str(args.input/'solvated.inpcrd')).positions
    monitor = LoopMonitor(args.input, args.protein, topology)
    print(json.dumps(monitor.check(positions, args.output, 'assembled', 0), indent=2))
