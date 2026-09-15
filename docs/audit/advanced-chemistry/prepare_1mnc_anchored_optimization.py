"""Prepare, never launch, a versioned 1MNC constrained-cluster research request.

No model acceptance, QM calculation, Hessian or production topology is generated.
All source inputs stay immutable; output creation refuses an existing directory.
"""
from __future__ import annotations

import hashlib
import importlib.util
import itertools
import json
import math
from pathlib import Path
import shutil
import time

ROOT = Path(__file__).resolve().parents[3]
AUDIT = ROOT / "docs/audit/advanced-chemistry"
OUT = ROOT / "build/advanced-chemistry/metal/qm-1mnc-anchored-optimization-v1"
DF = ROOT / "build/advanced-chemistry/metal/qm-1mnc-df-comparison-v2/calculation"
INPUTS = AUDIT / "prototype-1mnc-inputs"
EXPECTED_WORKER = "ddb67a1f8e9be6ed05c3c95e9864e20b79f0233604e27926e30e5ffc09295bc1"
EXPECTED_DF_RESULT = "de9be565922bb299c1e949f101dd2d5420cc9141de3ff9ca650fe4552ff3539c"


def read(path):
    raw = path.read_bytes()
    if not raw or len(raw) != path.stat().st_size:
        raise ValueError(f"Incomplete/empty source read: {path}")
    return raw


def sha(raw):
    return hashlib.sha256(raw).hexdigest()


def write(path, value):
    path.write_text(json.dumps(value, indent=2, allow_nan=False) + "\n")


def distance(a, b):
    return math.sqrt(sum((x-y)**2 for x, y in zip(a, b, strict=True)))


def volume(a, b, c, d):
    x, y, z = ([v-w for v, w in zip(p, d, strict=True)] for p in (a, b, c))
    return (x[0]*(y[1]*z[2]-y[2]*z[1]) - x[1]*(y[0]*z[2]-y[2]*z[0])
            + x[2]*(y[0]*z[1]-y[1]*z[0]))


def main():
    # Verify scientific provenance before creating anything.
    result_raw = read(DF / "result.json")
    if sha(result_raw) != EXPECTED_DF_RESULT:
        raise ValueError("Completed DF result changed")
    result = json.loads(result_raw)
    if not (result['accepted'] and result['scf']['converged'] and result['status'] == 'completed'):
        raise ValueError("Expected completed DF reference unavailable")
    for name, key in [('arrays.npz', 'arrays_sha256'), ('scf.chk', 'checkpoint_sha256'),
                      ('resolved-method.json', 'resolved_method_file_sha256')]:
        if sha(read(DF / name)) != result[key]:
            raise ValueError(f"DF source artifact changed: {name}")
    worker_raw = read(ROOT / "backend/qm_worker.py")
    if sha(worker_raw) != EXPECTED_WORKER:
        raise ValueError("Worker version requires explicit new review")
    request = json.loads(read(DF / "input.json"))
    original = json.loads(read(INPUTS / "small-screening-qm-request.json"))
    mapping = json.loads(read(INPUTS / "small-atom-map.json"))
    readiness = json.loads(read(AUDIT / "prototype-1mnc-df-independent-review/optimization-readiness.json"))
    cluster = json.loads(read(INPUTS / "cluster-input-report.json"))
    for key in ['atom_ids', 'elements', 'coords_angstrom', 'charge', 'spin', 'method', 'functional', 'basis', 'grid_level']:
        if request[key] != original[key]:
            raise ValueError(f"Reference model identity differs: {key}")
    ids, elements, xyz = request['atom_ids'], request['elements'], request['coords_angstrom']
    if len(ids) != 97 or sum(':PLH:' in a and e != 'H' for a, e in zip(ids, elements, strict=True)) != 25:
        raise ValueError("Whole inhibitor/cluster atom inventory changed")
    anchors = readiness['verified_cap_carbons_at_original_CA']
    frozen = [a['atom_id'] for a in anchors]
    if frozen != ['cap:small:15:CH3', 'cap:small:54:CH3', 'cap:small:92:CH3']:
        raise ValueError("Reviewed exact cap-carbon map changed")
    for a in anchors:
        i = ids.index(a['atom_id'])
        if elements[i] != 'C' or xyz[i] != a['coords_angstrom'] or mapping[i]['stable_id'] != ids[i]:
            raise ValueError("Cap identity/source coordinate verification failed")
    # DF checkpoint is a same-geometry initial density only, re-converged by worker.
    request.update(max_wall_seconds=21600, max_scf_cycles=100, scf_conv_tol=1e-10,
                   scf_conv_tol_grad=1e-7, threads=2, max_memory_mb=8000,
                   operations=['gradient'], density_fit=True, auxbasis='def2-universal-jkfit',
                   initial_checkpoint={'result_path': str(DF / 'result.json'),
                       'result_sha256': EXPECTED_DF_RESULT, 'checkpoint_path': str(DF / 'scf.chk'),
                       'checkpoint_sha256': result['checkpoint_sha256']},
                   optimization={'enabled': True, 'max_steps': 100, 'freeze_atom_ids': frozen,
                       'dihedral_constraints': [], 'convergence_energy': 1e-6,
                       'convergence_grms': 3e-4, 'convergence_gmax': 4.5e-4,
                       'convergence_drms': 1.2e-3, 'convergence_dmax': 1.8e-3})
    spec = importlib.util.spec_from_file_location('frozen_qm_validation', ROOT / 'backend/qm_worker.py')
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    resolved_input = module.validate_input(request)  # stdlib-only, no native import or calculation
    # Preserve explicit ligand bond order/parent-H records; never infer bonds from distances.
    mol2 = read(INPUTS / 'PLH.mol2').decode()
    atomtext = mol2.split('@<TRIPOS>ATOM\n', 1)[1].split('@<TRIPOS>BOND', 1)[0]
    atomnames = {int(r[0]): r[1] for r in (line.split() for line in atomtext.splitlines() if line.strip())}
    bondtext = mol2.split('@<TRIPOS>BOND\n', 1)[1].split('@<TRIPOS>', 1)[0]
    ligand_bonds = []
    neighbors = {n: [] for n in atomnames.values()}
    for row in (line.split() for line in bondtext.splitlines() if line.strip()):
        names = [atomnames[int(row[1])], atomnames[int(row[2])]]
        atom_ids = [f'E:280::PLH:{n}' for n in names]
        neighbors[names[0]].append(names[1]); neighbors[names[1]].append(names[0])
        ligand_bonds.append({'atom_ids': atom_ids, 'native_mol2_order': row[3],
                             'initial_distance_angstrom': distance(*(xyz[ids.index(a)] for a in atom_ids))})
    if len(atomnames) != 51 or len(ligand_bonds) != 51 or any(n.startswith('H') for n in neighbors['O1']) or 'H24' not in neighbors['N1']:
        raise ValueError("PLH graph or declared hydroxamate H inventory differs")
    donors = ['A:218::HID:NE2', 'A:222::HID:NE2', 'A:228::HID:NE2', 'E:280::PLH:O1', 'E:280::PLH:O2']
    zn = 'B:281::ZN:ZN'
    coord = {a: xyz[ids.index(a)] for a in ids}
    angles = []
    for a, b in itertools.combinations(donors, 2):
        u, v = ([x-y for x, y in zip(coord[d], coord[zn], strict=True)] for d in (a, b))
        cos = sum(x*y for x, y in zip(u, v, strict=True))/(distance(coord[a],coord[zn])*distance(coord[b],coord[zn]))
        angles.append({'atom_ids': [a, zn, b], 'initial_degrees': math.degrees(math.acos(max(-1,min(1,cos))))})
    chiral = []
    for name, cip in cluster['ligand_heavy_stereocenters_verified'].items():
        ns = sorted(neighbors[name])
        if len(ns) != 4:
            raise ValueError("Stereo center must retain four named neighbors")
        v = volume(*(coord[f'E:280::PLH:{n}'] for n in ns))
        if abs(v) < .1:
            raise ValueError("Source stereocenter is geometrically ambiguous")
        chiral.append({'atom_id': f'E:280::PLH:{name}', 'verified_source_CIP': cip,
                      'ordered_neighbor_ids': [f'E:280::PLH:{n}' for n in ns],
                      'initial_signed_determinant_angstrom3': v,
                      'meaning': 'Sign preservation for this fixed graph/order, independently of source CIP verification; not a bond-order inference.'})
    OUT.mkdir(parents=True, exist_ok=False)
    (OUT / 'qm_worker.py').write_bytes(worker_raw)
    write(OUT / 'input.json', request)
    write(OUT / 'validated-input.json', resolved_input)
    (OUT / 'native-unconstrained-small-opt.com').write_bytes(read(INPUTS / 'site_small_opt.com'))
    (OUT / 'native-unconstrained-small-opt.inp').write_bytes(read(INPUTS / 'site_small_opt.inp'))
    (OUT / 'small-atom-map.json').write_bytes(read(INPUTS / 'small-atom-map.json'))
    (OUT / 'PLH-declared-state.sdf').write_bytes(read(INPUTS / 'PLH-declared-state.sdf'))
    (OUT / 'PLH.mol2').write_bytes(read(INPUTS / 'PLH.mol2'))
    baseline = {'atom_count': 97, 'whole_PLH_atoms': 51, 'whole_PLH_heavy_atoms': 25,
        'frozen_atom_ids': frozen, 'free_atoms': 94, 'free_cartesian_degrees': 282,
        'anchors': anchors, 'declared_state': cluster['declared_state'],
        'coordination': [{'atom_ids': [zn, d], 'initial_distance_angstrom': distance(coord[zn],coord[d])} for d in donors],
        'zinc_centered_angles': angles, 'ligand_bond_graph': ligand_bonds,
        'stereochemistry': chiral,
        'all_other_N_O_distances_to_Zn': [{'atom_id': a, 'initial_distance_angstrom': distance(coord[zn],coord[a])}
             for a,e in zip(ids,elements,strict=True) if e in ['N','O'] and a not in donors],
        'cap_hydrogens_free': readiness['remaining_native_cap_atoms']}
    write(OUT / 'initial-geometry-diagnostics.json', baseline)
    sources = [
        {'id': 'amber_mcpb', 'url': 'https://ambermd.org/tutorials/advanced/tutorial20/mcpbpy.php',
         'role': 'Official implementation documentation',
         'supports': 'Optimize small model before same-level Hessian; inspect coordination after optimization; failure can require a different method or multistep strategy. Original preserved generated input is unconstrained.',
         'does_not_support': 'No claim that all five-coordinate hydroxamate sites work; no inherited cap-freeze policy.'},
        {'id': 'moubarak_2022', 'url': 'https://doi.org/10.3389/fmolb.2022.945415',
         'title': 'Structural and electronic properties of the active site of [ZnFe] SulE',
         'role': 'Primary original research, Methods: QM cluster optimization',
         'supports': 'Truncated source alpha carbons are a published location for coordinate anchoring in metal-site quantum cluster models.',
         'does_not_support': 'Different enzyme, metals, states and method. It does not validate this three-carbon-only 1MNC hypothesis or its gas-phase missing environment.'},
        {'id': 'geometric_constraints', 'url': 'https://geometric.readthedocs.io/en/latest/constraints.html',
         'role': 'Official optimizer documentation', 'supports': 'Explicit Cartesian freeze constraints with named atom-index mapping.'},
        {'id': 'geometric_convergence', 'url': 'https://geometric.readthedocs.io/en/latest/how-it-works.html',
         'role': 'Official optimizer documentation', 'supports': 'Five numerical convergence criteria; local installed source additionally inspected for actual per-atom gradient norms.'}]
    write(OUT / 'primary-source-policy.json', sources)
    policy = {'schema_version': 1, 'status': 'prepared_for_root_review_not_launched',
        'model_accuracy_validated': False, 'full_preparation_ready': False,
        'selected_hypothesis': 'Exact Cartesian anchors at three verified native cap carbons that occupy source His218/222/228 CA positions; all 94 remaining atoms free.',
        'choice_is_inference': True,
        'reason': 'Retain attachment locations imposed by the missing protein while allowing cap hydrogens, metal-donor geometry and full PLH to relax. Any coordination/state failure remains a failed hypothesis, never hidden by metal-distance constraints.',
        'not_selected': 'Native unconstrained small-cluster optimization is retained as a distinct possible sensitivity calculation, not silently overwritten or queued.',
        'omitted_environment': {'protein_and_solvent_embedding': 'None; gas-phase isolated first-shell model.',
            'missing_nonlocal_electrostatics': True, 'missing_steric_contacts': True,
            'full_PLH_relaxation': 'All 51 atoms free, including distal substituents; distortion can reflect missing pocket contacts. Record all ligand torsion and displacement diagnostics before parameter extraction.',
            'required_model_review': 'Large free-ligand distortion or coordination/proton transfer requires a revised environment/model; convergence cannot cure this limitation.',
            'other_sites_still_required': ['structural Zn282', 'Ca283']},
        'endpoint_checks': {
            'identity': 'Exact ordered 97 IDs/elements and +1/singlet state; original PLH graph/order retained, no deletion/mutation.',
            'numerical': 'Finite energy, coordinates and all derivatives; independently converged final DF SCF and geomeTRIC optimization.',
            'cap_displacement_max_bohr': 1e-5,
            'gradient_reporting': 'All 97 atomic vectors/components, the 94 free-atom vectors, and three cap vectors separately. Report RMS over atoms and components separately. Frozen-atom forces do not invalidate constrained convergence and must not be silently zeroed.',
            'independent_free_DF_gradient_max_vector_hartree_per_bohr': 4.5e-4,
            'independent_free_DF_gradient_rms_vector_hartree_per_bohr': 3e-4,
            'coordination_alarms_angstrom': {'lower': 1.7, 'upper': 2.8, 'scope': 'All five named donors; broad investigation bounds, not fitted equilibrium targets.'},
            'new_N_O_neighbor_alarm_angstrom': 2.8,
            'zinc_centered_angle_change_review_degrees': 30,
            'stereo': 'Both named tetrahedral determinant signs unchanged, absolute magnitude >=0.1 Angstrom^3, and independent graph-based final C3(R)/C9(S) assignment.',
            'proton_state': 'All H atoms retain their named parent connectivity by geometry review; O1 stays H-free; three NE2 donors stay H-free; HID ND1-H and PLH N1-H retained. Ambiguous proton transfer blocks extraction.',
            'C_N_O_covalent_bond_alarm_angstrom': [0.95,1.90],
            'hydrogen_parent_bond_alarm_angstrom': [0.70,1.35],
            'nonbonded_heavy_contact_alarm_angstrom': 1.0,
            'retained_heavy_frame_rmsd_review_angstrom': 1.0,
            'retained_heavy_max_displacement_review_angstrom': 2.0,
            'required_displacement_reports': 'Source-frame every-atom displacement, cap-pair distances, Zn/donor/core RMSD, whole/distal PLH heavy RMSD, ligand proper-torsion changes and contact list against original omitted protein after source-frame placement.',
            'alarm_policy': 'Prospective engineering diagnostics chosen before results. Any alarm stops automatic downstream extraction and requires documented model review; absence of alarms is not chemical acceptance.'},
        'conventional_endpoint_check': {
            'required_before_hessian': True,
            'method': 'Same final coordinates, same RKS/B3LYPG spherical 6-31g*, grid4, charge/spin, conventional integrals; DF final density only as verified initial guess.',
            'max_vector_DF_gradient_difference_hartree_per_bohr': 1e-4,
            'rms_vector_DF_gradient_difference_hartree_per_bohr': 3e-5,
            'free_conventional_max_gradient_vector_hartree_per_bohr': 5.5e-4,
            'free_conventional_rms_gradient_vector_hartree_per_bohr': 3.3e-4,
            'all_raw_energy_and_force_differences_retained': True,
            'basis_of_limits': 'Predeclared numerical screening limits, not literature chemical-accuracy guarantees; initial DF/conventional same-point differences do not establish endpoint equivalence.',
            'automatic_launch': False},
        'hessian_policy': 'No Hessian requested now. A constrained stationary point is not a full unconstrained minimum. Future native electronic Hessian must exclude artificial anchor-force terms, be analyzed in the correct free subspace, and separately inspect boundary coupling; a standard full-space six-zero-mode minimum claim is invalid. Native MCPB extraction remains a separate adapter/model validation task.',
        'resources': {'threads': 2, 'max_memory_mb_hint': 8000, 'hard_memory_limit': False,
            'native_wall_seconds': 21600, 'outer_wall_seconds': 21720,
            'max_optimization_steps': 100, 'max_all_QM_workers': 3, 'max_large_DF_workers': 2,
            'start_free_GiB_if_no_other_large_DF': 22, 'start_free_GiB_if_other_large_DF': 35,
            'running_free_GiB_reserve': 8,
            'launch_lock': str(ROOT / 'build/advanced-chemistry/QM-LAUNCH.lock'),
            'launch_rule': 'Root review plus real PID/resource inspection and lock-protected duplicate check. Future controller must stop only its own process on low disk/deadline, preserve all partial logs/results and never extend deadlines silently.',
            'memory_and_IO': 'Inspect memory pressure, swap and file-provider hydration before launch; thread/memory settings are not an OS cap.'},
        'launch_authorized_by_this_artifact': False, 'runtime': str(ROOT / '.tools/qm-parallel/bin/python'),
        'input_sha256': sha(read(OUT / 'input.json')), 'worker_sha256': EXPECTED_WORKER,
        'created_unix': time.time()}
    write(OUT / 'optimization-policy.json', policy)
    source_paths = [Path(__file__), ROOT/'backend/qm_worker.py', DF/'input.json', DF/'result.json',
        DF/'arrays.npz', DF/'scf.chk', DF/'resolved-method.json', INPUTS/'small-atom-map.json',
        INPUTS/'small-screening-qm-request.json', INPUTS/'cluster-input-report.json',
        INPUTS/'site_small_opt.com', INPUTS/'site_small_opt.inp', INPUTS/'PLH.mol2',
        INPUTS/'PLH-declared-state.sdf', AUDIT/'inputs/1MNC.cif', AUDIT/'inputs/ccd/PLH.cif',
        AUDIT/'prototype-1mnc-df-independent-review/optimization-readiness.json',
        ROOT/'.tools/qm-parallel/dynamol-qm-runtime.json',
        ROOT/'.tools/qm-parallel/lib/python3.12/site-packages/geometric/optimize.py']
    write(OUT / 'source-manifest.json', {str(p): {'sha256': sha(read(p)), 'bytes': p.stat().st_size} for p in source_paths})
    write(OUT / 'artifact-manifest.json', {p.name: {'sha256':sha(read(p)), 'bytes':p.stat().st_size} for p in OUT.iterdir() if p.is_file()})
    print(json.dumps({'output': str(OUT), 'status': policy['status'], 'input_sha256': policy['input_sha256'],
        'policy_sha256': sha(read(OUT/'optimization-policy.json')), 'QM_launched': False,
        'initial_distances_angstrom': baseline['coordination']}, indent=2))


if __name__ == '__main__':
    main()
