"""Export a screened research loop candidate for native complex assembly.

No coordinates are refined here; histidine names only express the hydrogens
already present. Temporary construction forces are not part of a PDB export.
"""
import argparse
import json
from pathlib import Path

import numpy as np
from openmm import app, unit

from prepare_adduct import digest


def run(candidate, audit, environment, output):
    local = json.loads((candidate/'result.json').read_text())
    checked = json.loads((audit/'result.json').read_text())
    if not local['passed_local_checks'] or not checked['passed']:
        raise ValueError('The candidate must pass local and retained-environment checks')
    for name in [candidate/'refined-candidate.pdb', candidate/'result.json', environment]:
        if checked['source_sha256'].get(str(name)) != digest(name):
            raise ValueError('Candidate/environment is not bound to the passing audit')
    for name, expected in checked['source_sha256'].items():
        if digest(name) != expected:
            raise ValueError('Audited source changed: '+name)
    molecule = app.PDBFile(str(candidate/'refined-candidate.pdb'))
    from rama_reference import rama_report
    rama = rama_report(molecule.topology, molecule.positions,
                       local['modeled_residues']+local['native_changed_context_residues'])
    if not rama['all_selected_scored_without_outliers']:
        rejected = [row['residue'] for row in rama['outliers']+rama['unavailable']]
        raise ValueError('Native Ramachandran reference rejects the loop/context before export: '+str(rejected))
    output.mkdir(exist_ok=False)
    (output/'loop-ramachandran.json').write_text(json.dumps(rama, indent=2)+'\n')
    variants = []
    for residue in molecule.topology.residues():
        if residue.name != 'HIS':
            continue
        names = {a.name for a in residue.atoms()}
        ring_hydrogens = frozenset(names & {'HD1', 'HE2'})
        variant = {frozenset({'HD1'}): 'HID', frozenset({'HE2'}): 'HIE',
                   frozenset({'HD1', 'HE2'}): 'HIP'}.get(ring_hydrogens)
        if variant is None:
            raise ValueError('Histidine hydrogen pattern has no unambiguous Amber variant')
        variants.append({'chain': residue.chain.id, 'resid': residue.id,
                         'variant': variant, 'basis': 'Existing generated hydrogen names'})
        residue.name = variant
    with (output/'protein.pdb').open('w') as handle:
        app.PDBFile.writeFile(molecule.topology, molecule.positions, handle, keepIds=True)
    saved = app.PDBFile(str(output/'protein.pdb'))
    if saved.topology.getNumAtoms() != molecule.topology.getNumAtoms():
        raise ValueError('Intermediate serialization changed the atom inventory')
    before = np.asarray(molecule.positions.value_in_unit(unit.nanometer))
    after = np.asarray(saved.positions.value_in_unit(unit.nanometer))
    if not np.array_equal(before, after):
        raise ValueError('Intermediate serialization moved coordinates')
    (output/'retained-environment.pdb').write_bytes(environment.read_bytes())
    report = {'stage': 'protein_intermediate_complete', 'full_complex_prepared': False,
        'app_ready': False, 'physical_model_validated': False,
        'atoms': molecule.topology.getNumAtoms(), 'modeled_residues': local['modeled_residues'],
        'native_changed_context_residues': local['native_changed_context_residues'],
        'environment_preserved_separately': True, 'retained_environment_check': checked,
        'loop_geometry': json.loads((audit/'geometry.json').read_text()),
        'loop_ramachandran': str(output/'loop-ramachandran.json'),
        'native_ramachandran_reference_passed': True,
        'loop_stereochemistry': str(audit/'stereochemistry.json'),
        'loop_refinement': str(candidate/'refinement.json'), 'histidine_variants': variants,
        'source_sha256': {str(p): digest(p) for p in [candidate/'refined-candidate.pdb',
             candidate/'result.json', candidate/'refinement.json', audit/'result.json', environment]},
        'outputs_sha256': {n: digest(output/n) for n in ['protein.pdb', 'retained-environment.pdb', 'loop-ramachandran.json']},
        'construction_forces_exported': False,
        'scope': 'Screened protein coordinate intermediate only. Complete-complex minimization and dynamics are required; no new terminal segments or source molecules were added or removed here.'}
    (output/'result.json').write_text(json.dumps(report, indent=2)+'\n')
    print(json.dumps({'stage': report['stage'], 'output': str(output), 'atoms': report['atoms'],
                      'histidine_variants': variants, 'app_ready': False}, indent=2))


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    for name in ['candidate', 'audit', 'environment', 'output']:
        parser.add_argument('--'+name, type=Path, required=True)
    args = parser.parse_args()
    run(args.candidate, args.audit, args.environment, args.output)
