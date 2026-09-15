"""Audit a saved loop candidate against unchanged retained molecules.

Read-only scientific screening: no minimization, parameter changes or app
promotion. Gross geometry does not establish a native loop conformation.
"""
import argparse
import json
from pathlib import Path
import sys

import numpy as np
from openmm import app, unit
from pdbfixer import PDBFixer

from prepare_adduct import digest


def run(repo, original, candidate, environment, output):
    output.mkdir(exist_ok=False)
    sys.path.insert(0, str(repo))
    from backend.loop_geometry import _RADII, loop_geometry_report
    from backend.preparation_worker import atom_key, stereochemistry_report, require_valid_stereochemistry

    paths = [original/'before-refinement.pdb', candidate/'refined-candidate.pdb',
             candidate/'result.json', environment]
    hashes = {str(p): digest(p) for p in paths}
    local = json.loads((candidate/'result.json').read_text())
    if not local['passed_local_checks']:
        raise ValueError('Candidate did not pass its local construction checks')
    modeled = {tuple(k) for k in local['modeled_residues']}
    context = {tuple(k) for k in local['native_changed_context_residues']}
    source = PDBFixer(filename=str(original/'before-refinement.pdb'))
    saved = app.PDBFile(str(candidate/'refined-candidate.pdb'))
    retained = app.PDBFile(str(environment))

    def heavy_inventory(molecule):
        atoms = list(molecule.topology.atoms())
        if any(a.element is None for a in atoms):
            raise ValueError('Unspecified element in source, candidate or retained environment')
        heavy = [a for a in atoms if a.element.symbol not in {'H', 'D'}]
        xyz = np.asarray(molecule.positions.value_in_unit(unit.nanometer))
        inventory = {atom_key(a): (a.element.symbol, xyz[a.index]) for a in heavy}
        if len(inventory) != len(heavy):
            raise ValueError('Ambiguous heavy-atom identity')
        return inventory

    before = heavy_inventory(source)
    after = heavy_inventory(saved)
    outside = heavy_inventory(retained)
    if before.keys() != after.keys() or any(before[k][0] != after[k][0] for k in before):
        raise ValueError('Saved candidate changed the source heavy-atom inventory')
    if after.keys() & outside.keys():
        raise ValueError('Retained environment duplicates protein atom identities')
    displacements = [dict(atom=list(k), displacement_nm=float(np.linalg.norm(after[k][1]-before[k][1])))
                     for k in before if k[:4] not in modeled]
    fixed = [row for row in displacements if tuple(row['atom'][:4]) not in context]
    max_fixed = max((row['displacement_nm'] for row in fixed), default=0)
    max_observed = max((row['displacement_nm'] for row in displacements), default=0)
    if max_fixed > .00018 or max_observed > .1:
        raise ValueError('Saved coordinates violate fixed-atom precision or 1 angstrom context limit')
    external = [(xyz, _RADII.get(element, .170)) for element, xyz in outside.values()]
    geometry = loop_geometry_report(saved.topology, saved.positions, modeled | context,
                                    environment_positions=external)
    stereo = stereochemistry_report(saved.topology, saved.positions, source.templates)
    stereo_error = None
    try:
        require_valid_stereochemistry(stereo, 'Saved loop candidate with retained environment')
    except Exception as exc:
        stereo_error = str(exc)
    for name, value in [('geometry.json', geometry), ('stereochemistry.json', stereo),
                        ('observed-displacements.json', displacements)]:
        (output/name).write_text(json.dumps(value, indent=2)+'\n')
    if hashes != {str(p): digest(p) for p in paths}:
        raise ValueError('Source changed during audit')
    report = dict(stage='saved_loop_retained_environment_screen_complete',
        passed=geometry['accepted'] and stereo_error is None,
        app_ready=False, physical_model_validated=False, full_environment_validated=False,
        retained_environment_geometry_checked=True, source_atom_inventory_preserved=True,
        protein_heavy_atoms=len(after), retained_environment_heavy_atoms=len(outside),
        maximum_fixed_atom_displacement_nm=max_fixed,
        maximum_observed_context_displacement_nm=max_observed,
        geometry_errors=geometry['errors'], stereo_error=stereo_error,
        requires_minimization=geometry['requires_minimization'],
        external_soft_overlap_count=sum(r['external_environment_index'] is not None for r in geometry['soft_overlaps']),
        source_sha256=hashes,
        scope='PDB-roundtrip identity, displacement, stereo and gross geometry screen of modeled and changed context residues against all retained nonprotein heavy atoms. Complete-complex mechanics and dynamics remain untested.')
    (output/'result.json').write_text(json.dumps(report, indent=2)+'\n')
    print(json.dumps(report, indent=2))


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    for name in ['repo', 'original', 'candidate', 'environment', 'output']:
        parser.add_argument('--'+name, type=Path, required=True)
    args = parser.parse_args()
    try:
        run(args.repo, args.original, args.candidate, args.environment, args.output)
    except Exception as exc:
        if args.output.is_dir():
            (args.output/'failure.json').write_text(json.dumps({'type': type(exc).__name__, 'error': str(exc)}, indent=2)+'\n')
        raise
