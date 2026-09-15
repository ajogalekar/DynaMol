"""Fail-closed research handoff: exact-parent RESP -> native adduct mechanics.

This adapter never supplies missing charges and never declares a force field
physically accepted. Its only input route is a completed, hash-bound parent
optimization / conventional-HF-ESP / constrained-RESP continuation. The
source plan is pinned before charge completion. Synthetic unit fixtures may
exercise individual validators; there is no synthetic bypass in the CLI.
"""
from __future__ import annotations

import argparse
from collections import Counter
import hashlib
import importlib.util
import json
import math
import os
from pathlib import Path
import shutil
import subprocess
import time

ROOT = Path(__file__).resolve().parents[3]
BASE = Path(__file__).resolve().parent


def read_bytes(path):
    path = Path(path)
    payload = path.read_bytes()
    if not payload or len(payload) != path.stat().st_size:
        raise ValueError('Artifact is empty or incompletely downloaded: '+str(path))
    return payload


def sha(path):
    return hashlib.sha256(read_bytes(path)).hexdigest()


def read(path):
    return json.loads(read_bytes(path))


def write(path, value):
    Path(path).write_text(json.dumps(value, indent=2, allow_nan=False)+'\n')


def module(path, name):
    spec = importlib.util.spec_from_file_location(name, path)
    result = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(result)
    return result


def require(condition, message):
    if not condition:
        raise ValueError(message)


def verify_artifacts(artifacts):
    require(bool(artifacts), 'No pinned artifacts')
    for path, expected in artifacts.items():
        require(sha(path) == expected, 'Pinned source artifact changed: '+path)


def validate_terminal_state(state, plan):
    require(state.get('status') == 'research_candidate' and
            state.get('stage') == 'awaiting_full_adduct_torsion_evidence',
            'Exact-parent RESP continuation is unfinished or failed')
    require(state.get('full_preparation_ready') is False,
            'Unexpected physical acceptance claim in research continuation')
    for key in ('cap_graph_sha256', 'cap_sdf_sha256', 'optimizer_input_sha256'):
        require(state.get(key) == plan[key], 'Continuation source identity changed: '+key)
    require(state.get('script_sha256') in plan['allowed_continuation_sha256'],
            'Continuation code is not a reviewed pinned version')
    require(Path(state['optimizer']).resolve() == Path(plan['parent']).resolve(),
            'Continuation belongs to another parent')
    require(state.get('stages', {}).get('resp', {}).get('fit', {}).get('constraint_checks_passed') is True,
            'RESP did not complete its explicit charge constraints')


def validate_charges(source, constraints, expected_constraints, fit, raw_charges):
    require(constraints == expected_constraints, 'Canonical charge constraints or equivalence policy changed')
    records = source['atom_map']
    count = len(records)
    require([r['index'] for r in records] == list(range(count)), 'Source atom indices are not contiguous')
    require(len({r['name'] for r in records}) == count, 'Source atom names are not unique')
    require(count == constraints['atom_count'] and source['formal_charge'] == constraints['formal_charge'],
            'Charge inventory or total state differs from the source graph')
    charges = fit.get('charges', [])
    require(fit.get('constraint_checks_passed') is True and len(charges) == count,
            'Missing or unfinished jointly fitted charge vector')
    require(all(type(q) in (int, float) and math.isfinite(q) for q in charges), 'Invalid fitted charge value')
    require(charges == raw_charges, 'Charge vector differs from native RESP stage2 output')
    require(abs(sum(charges)-source['formal_charge']) < 1e-4, 'Joint adduct charge changed')
    fixed = constraints['frozen_atoms']
    require(len(fixed) == 18 and len({r['index'] for r in fixed}) == 18,
            'Expected the reviewed 18 canonical cap/backbone charges')
    for row in fixed:
        require(records[row['index']]['name'] == row['capped_name'], 'Canonical charge atom identity changed')
        require(abs(charges[row['index']]-row['charge_e']) <= 1e-7, 'Canonical cap/backbone charge changed')
    cap_indices = [r['index'] for r in records if r['role'] == 'temporary-cap' or
                   (r['role'] == 'generated-hydrogen' and
                    records[r['parent_index']]['role'] == 'temporary-cap')]
    require(cap_indices == constraints['cap_indices'], 'Temporary cap removal inventory changed')
    require(abs(sum(charges[i] for i in cap_indices)) < 1e-6, 'Cap removal changes the formal charge')
    for group in constraints['stage2_methyl_equivalence_groups']:
        require(max(charges[i] for i in group[1:])-min(charges[i] for i in group[1:]) <= 1e-7,
                'Native RESP methyl equivalence changed')
    return {'charges_e': list(charges), 'removed_cap_indices': cap_indices,
            'atom_ids': [f'{r["index"]:03d}:{r["name"]}' for r in records],
            'canonical_fixed_count': len(fixed), 'charge_sum_e': sum(charges),
            'retained_charge_sum_e': sum(q for i, q in enumerate(charges) if i not in cap_indices)}


def verify_evidence(plan):
    # Refuse pending chemistry before importing expensive/offloaded libraries.
    continuation = Path(plan['continuation'])
    state = read(continuation/'result.json')
    validate_terminal_state(state, plan)
    verify_artifacts(plan['pinned_artifacts'])
    import numpy as np
    from rdkit import Chem

    cap, parent = Path(plan['cap']), Path(plan['parent'])
    source = read(cap/'capped-adduct.json')
    request, result = read(parent/'input.json'), read(parent/'result.json')
    require(result.get('accepted') is True and result.get('scf', {}).get('converged') is True,
            'Parent result is not an accepted converged quantum calculation')
    require(result['input_sha256'] == plan['optimizer_input_sha256'], 'Parent input binding changed')
    require(sha(parent/'result.json') == state['stages']['parent']['result_sha256'], 'Parent result hash changed')
    require(sha(parent/'arrays.npz') == result['arrays_sha256'] == state['stages']['parent']['arrays_sha256'],
            'Parent coordinate artifact changed')
    arrays = np.load(parent/'arrays.npz', allow_pickle=False)
    checker = module(BASE/'continue_covalent_qm.py', 'handoff_geometry_check')
    molecule, geometry = checker.check_geometry(cap, request, result, arrays)
    require(geometry['passed'] is True, 'Optimized adduct geometry failed identity/stereo/bond checks')
    require(read(continuation/'geometry.json')['passed'] is True, 'Original geometry check was not passed')

    espdir = continuation/'exact-esp'
    esp_request, esp_result = read(continuation/'exact-esp-input.json'), read(espdir/'result.json')
    require(sha(continuation/'exact-esp-input.json') == esp_result['input_sha256'] ==
            state['stages']['exact_esp_launch']['input_sha256'], 'Exact ESP input binding changed')
    require(esp_result.get('accepted') is True and esp_result.get('scf', {}).get('converged') is True,
            'Conventional HF ESP did not converge')
    require(esp_request.get('method') == 'RHF' and esp_request.get('density_fit') is False and
            esp_request.get('basis', '').lower() == '6-31g*' and not esp_request.get('ecp') and
            esp_request.get('charge') == source['formal_charge'] and esp_request.get('spin') == 0 and
            esp_request.get('operations') == ['esp'], 'ESP method or electronic state changed')
    require(esp_request['initial_checkpoint']['result_sha256'] == sha(parent/'result.json'),
            'Exact ESP was initialized from a different parent')
    require(esp_request['initial_checkpoint']['checkpoint_sha256'] == result['checkpoint_sha256'] ==
            sha(parent/'scf.chk'), 'Exact ESP parent checkpoint binding changed')
    for directory, outcome, density_fit in ((parent, result, True), (espdir, esp_result, False)):
        require(outcome.get('charge') == source['formal_charge'] and outcome.get('spin_2S') == 0,
                'Resolved quantum result state changed')
        resolved = read(directory/'resolved-method.json')
        require(sha(directory/'resolved-method.json') == outcome['resolved_method_file_sha256'] and
                resolved['method'] == 'RHF' and resolved['density_fit'] is density_fit and
                resolved['basis_requested'].lower() == '6-31g*' and not resolved['ecp_resolved'],
                'Resolved quantum method or basis changed')
    require(result['method']['basis_sha256'] == esp_result['method']['basis_sha256'],
            'Parent and ESP expanded basis differ')
    require(sha(espdir/'result.json') == state['stages']['exact_esp_result']['result_sha256'], 'ESP result hash changed')
    require(sha(espdir/'arrays.npz') == esp_result['arrays_sha256'] == state['stages']['exact_esp_result']['arrays_sha256'],
            'ESP array hash changed')
    esp_arrays = np.load(espdir/'arrays.npz', allow_pickle=False)
    ids = [f'{r["index"]:03d}:{r["name"]}' for r in source['atom_map']]
    for candidate in (request, result, esp_request, esp_result):
        require(candidate['atom_ids'] == ids and candidate['elements'] == request['elements'],
                'Parent/ESP atom identity or element order changed')
    require(np.array_equal(esp_arrays['coords_bohr'], arrays['coords_bohr']), 'ESP uses a different geometry')
    require(np.array_equal(np.asarray(esp_request['coords_bohr']), arrays['coords_bohr']), 'ESP request geometry changed')
    points, potentials = esp_arrays['esp_points_bohr'], esp_arrays['esp_hartree_per_e']
    require(points.ndim == 2 and points.shape[1] == 3 and potentials.shape == (len(points),) and
            len(points) > len(ids) and np.isfinite(points).all() and np.isfinite(potentials).all(), 'Invalid ESP evidence')

    fitfolder = continuation/'fitting'
    constraints, fit = read(fitfolder/'resp/constraints.json'), read(fitfolder/'resp/fit-constraints.json')
    require(fit == state['stages']['resp']['fit'], 'RESP report differs from the terminal continuation')
    require(sha(fitfolder/'capped-adduct.json') == plan['cap_graph_sha256'], 'RESP atom map changed')
    require(sha(fitfolder/'capped-adduct.sdf') == sha(continuation/'optimized-adduct.sdf'), 'RESP fitted a different molecular graph')
    fitted_mol = Chem.SDMolSupplier(str(fitfolder/'capped-adduct.sdf'), removeHs=False)[0]
    require(fitted_mol is not None and Chem.MolToSmiles(fitted_mol) == Chem.MolToSmiles(molecule),
            'RESP molecule graph or stereochemistry differs from accepted parent')
    require(np.allclose(fitted_mol.GetConformer().GetPositions(), arrays['coords_angstrom'], rtol=0, atol=5.1e-5),
            'RESP SDF coordinates differ beyond its recorded precision')
    require(fit['esp_sha256'] == sha(continuation/'hf-resp.esp') == sha(fitfolder/'resp/esp.dat'),
            'Native RESP ESP file changed')
    # Bind the rounded native RESP representation back to the exact QM arrays.
    lines = read_bytes(continuation/'hf-resp.esp').decode().splitlines()
    require([int(x) for x in lines[0].split()] == [len(ids), len(points), 0], 'Native RESP ESP inventory changed')
    nuclei = np.asarray([[float(x) for x in row.split()] for row in lines[1:1+len(ids)]])
    samples = np.asarray([[float(x) for x in row.split()] for row in lines[1+len(ids):]])
    require(np.allclose(nuclei, arrays['coords_bohr'], rtol=5.1e-8, atol=5.1e-8) and
            samples.shape == (len(points), 4) and
            np.allclose(samples[:, 0], potentials, rtol=5.1e-8, atol=5.1e-12) and
            np.allclose(samples[:, 1:], points, rtol=5.1e-8, atol=5.1e-8), 'RESP file does not reproduce accepted QM ESP')
    require(fit['resp_executable_sha256'] == plan['resp_executable_sha256'] and
            sha(fit['resp_executable']) == plan['resp_executable_sha256'], 'RESP executable changed')
    require(fit['fitter_source_sha256'] == plan['fitter_source_sha256'] == state['stages']['resp']['fitter_sha256'],
            'RESP fitter source changed')
    qout = [float(x) for x in read_bytes(fitfolder/'resp/stage2.qout').decode().split()]
    values = validate_charges(source, constraints, read(cap/'resp/constraints.json'), fit, qout)
    for relative, expected in constraints['source_sha256'].items():
        require(sha(ROOT/relative) == expected, 'Canonical Amber library changed')
    evidence = [parent/'input.json', parent/'result.json', parent/'arrays.npz', continuation/'result.json',
                continuation/'exact-esp-input.json', espdir/'result.json', espdir/'arrays.npz',
                continuation/'hf-resp.esp', fitfolder/'capped-adduct.sdf', fitfolder/'resp/constraints.json',
                fitfolder/'resp/fit-constraints.json', fitfolder/'resp/stage2.qout']
    return source, constraints, values, arrays['coords_angstrom'], {
        'artifact_sha256': {str(p): sha(p) for p in evidence}, 'geometry': geometry,
        'charge_validation': values, 'physical_parameter_acceptance': False,
        'scope': 'Completed quantum/RESP numerical evidence and identity transport only; torsion and held-out physical validation remain required.'}


def rewrite_mol2(text, source, charges, coords):
    require(len(charges) == len(coords) == len(source['atom_map']), 'Native charge/coordinate inventory differs')
    require(all(math.isfinite(float(x)) for row in coords for x in row), 'Native coordinates are nonfinite')
    mapping = source['native_parameters']['validation']['atom_mapping']['atoms']
    require([r['native_index'] for r in mapping] == list(range(len(charges))), 'Native source index order changed')
    lines, section, seen = text.splitlines(), None, []
    for i, line in enumerate(lines):
        if line.startswith('@<TRIPOS>'):
            section = line
            continue
        if section == '@<TRIPOS>ATOM' and line.strip():
            row = line.split()
            index = int(row[0])-1
            require(index == len(seen) and row[1] == mapping[index]['prepared_name'], 'Native MOL2 atom identity/order changed')
            require(len(row) >= 9, 'Incomplete native MOL2 atom record')
            row[2:5] = [f'{float(x):.10f}' for x in coords[index]]
            row[8] = f'{charges[index]:.10f}'
            lines[i] = ' '.join(row)
            seen.append(index)
    require(len(seen) == len(charges), 'Native MOL2 inventory is incomplete')
    # Explicitly replace the old AM1-BCC metadata as well as the actual charges.
    marker = lines.index('@<TRIPOS>MOLECULE')
    require(lines[marker+4].lower() in ('bcc', 'user_charges'), 'Unexpected MOL2 charge convention')
    lines[marker+4] = 'USER_CHARGES'
    return '\n'.join(lines)+'\n'


def run(plan_path, output, assemble=False):
    plan, output = read(plan_path), Path(output).resolve()
    output.mkdir(parents=True, exist_ok=False)
    report = {'started_unix': time.time(), 'status': 'checking', 'assembled': False,
              'simulation_ready': False, 'physical_parameter_acceptance': False,
              'adapter_sha256': sha(__file__), 'plan_sha256': sha(plan_path)}
    try:
        source, constraints, values, coords, evidence = verify_evidence(plan)
        report['evidence'] = evidence
        interface = Path(plan['interface'])
        updated = rewrite_mol2(read_bytes(interface/'hybrid.mol2').decode(), source, values['charges_e'], coords)
        (output/'hybrid.mol2').write_text(updated)
        shutil.copyfile(interface/'gaff-adduct.frcmod', output/'gaff-adduct.frcmod')
        write(output/'atom-map-and-charges.json', {'source_atom_map': source['atom_map'], **values})
        report['status'] = 'charge_transport_candidate'
        if assemble:
            native_assembly(output, plan, source, constraints, values, coords, report)
        # Recheck every input after work so concurrent edits cannot silently pass.
        verify_artifacts(plan['pinned_artifacts'])
        verify_artifacts(evidence['artifact_sha256'])
    except Exception as error:
        report.update(status='failed', error_type=type(error).__name__, error=str(error))
    report['elapsed_seconds'] = time.time()-report['started_unix']
    write(output/'result.json', report)
    return report


def native_assembly(output, plan, source, constraints, values, coords, report):
    import numpy as np
    import parmed

    amber = ROOT/'.tools/ambertools'
    env = dict(os.environ, AMBERHOME=str(amber), DYLD_FALLBACK_LIBRARY_PATH=str(amber/'lib'),
               OMP_NUM_THREADS='2', OPENBLAS_NUM_THREADS='2')
    leap = ('source leaprc.protein.ff14SB\nsource leaprc.gaff2\n'
            'loadamberparams gaff-adduct.frcmod\nadduct=loadmol2 hybrid.mol2\ncheck adduct\n'
            'saveamberparm adduct hybrid.prmtop hybrid.inpcrd\nquit\n')
    (output/'leap.in').write_text(leap)
    process = subprocess.run([str(amber/'bin/tleap'), '-f', 'leap.in'], cwd=output, env=env,
                             capture_output=True, text=True, timeout=120)
    (output/'tleap.log').write_text(process.stdout+'\n'+process.stderr)
    require(process.returncode == 0, 'Native capped adduct assembly failed')
    native = parmed.load_file(str(output/'hybrid.prmtop'), xyz=str(output/'hybrid.inpcrd'))
    expected_names = [r['prepared_name'] for r in source['native_parameters']['validation']['atom_mapping']['atoms']]
    require([a.name for a in native.atoms] == expected_names, 'Capped native atom order changed')
    require(np.allclose([a.charge for a in native.atoms], values['charges_e'], rtol=0, atol=1e-8),
            'Native assembly changed fitted charges')
    require(np.allclose(native.coordinates, coords, rtol=0, atol=5.1e-7), 'Native assembly moved the optimized geometry')
    baseline = parmed.load_file(str(Path(plan['interface'])/'hybrid.prmtop'))
    require(mechanical_inventory(native) == mechanical_inventory(baseline),
            'Native assembly changed graph, atom types or non-charge mechanical parameters')
    for name in ('reference.prmtop',):
        shutil.copyfile(Path(plan['interface'])/name, output/name)
    checker = module(BASE/'covalent_interface_spike.py', 'handoff_native_backbone')
    report['canonical_backbone'] = checker.compare_canonical_backbone(output, constraints)
    polymer = module(BASE/'covalent_polymer_interface_spike.py', 'handoff_native_polymer')
    result = polymer.run(Path(plan['cap']), output, Path(plan['continuation'])/'fitting',
                         output/'polymer', 'quantum_parameter_research')
    model = parmed.load_file(str(output/'polymer/peptide.prmtop'))
    for row in result['adduct_mapping']:
        require(abs(model.atoms[row['peptide_index']].charge-values['charges_e'][row['capped_index']]) < 1e-8,
                'Peptide assembly changed a jointly fitted adduct charge')
    report.update(assembled=True, status='native_research_candidate', polymer=result,
                  all_capped_noncharge_mechanics_match_pinned_v4=True,
                  parameter_precedence='Native ff14SB canonical terms, then pinned adduct-only GAFF supplement; no protein parmchk2 overrides, charge renormalization, or charge substitution.')


def mechanical_inventory(model):
    """Charge-independent native parameters; repeated torsion terms stay counted."""
    result = {'atoms': [(a.name, a.type, a.atomic_number, a.mass, a.rmin, a.epsilon,
                         a.rmin_14, a.epsilon_14) for a in model.atoms]}
    for group, count, fields in (('bonds', 2, ('req', 'k')), ('angles', 3, ('theteq', 'k')),
                                ('dihedrals', 4, ('per', 'phase', 'phi_k', 'scee', 'scnb'))):
        rows = []
        for term in getattr(model, group):
            atoms = tuple(getattr(term, f'atom{i+1}').idx for i in range(count))
            require(term.type is not None, 'Missing native mechanical term')
            values = tuple(float(getattr(term.type, field)) for field in fields)
            require(all(math.isfinite(v) for v in values), 'Nonfinite native mechanical term')
            if group == 'dihedrals': values += (term.improper, term.ignore_end)
            rows.append((min(atoms, atoms[::-1]), values))
        result[group] = Counter(rows)
    require(not model.adjusts and not model.cmaps and not model.impropers,
            'Unsupported extra native terms need their own fidelity comparison')
    return result


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('plan', type=Path)
    parser.add_argument('output', type=Path)
    parser.add_argument('--assemble', action='store_true')
    args = parser.parse_args()
    final = run(args.plan, args.output, args.assemble)
    print(json.dumps({k: v for k, v in final.items() if k != 'evidence'}, indent=2))
    raise SystemExit(final['status'] == 'failed')
