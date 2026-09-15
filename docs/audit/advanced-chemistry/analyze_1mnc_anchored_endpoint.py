"""Read-only endpoint diagnostics for the fixed 1MNC anchored research request.

No QM, Hessian, charge fitting, force-field extraction or model acceptance.
An optional conventional-gradient input is written only after every declared
endpoint gate passes. Missing endpoints are not read or treated as failures.
"""
from __future__ import annotations

import argparse
import hashlib
import itertools
import json
import math
from pathlib import Path
import time
import xml.etree.ElementTree as ET

import gemmi
import numpy as np
from rdkit import Chem, rdBase

ROOT = Path(__file__).resolve().parents[3]
BASE = ROOT / 'build/advanced-chemistry/metal/qm-1mnc-anchored-optimization-v1'
AUDIT = ROOT / 'docs/audit/advanced-chemistry'
INPUT_HASH = '3b6fc0c51200556bd78bf884760b77c7feeeb079d614ceee253480e9887a5cc5'
POLICY_HASH = '61c06fa5679fbaafd7db06afb225b57b021985b7976fa3e7aca8a94da3ebd93f'
WORKER_HASH = 'ddb67a1f8e9be6ed05c3c95e9864e20b79f0233604e27926e30e5ffc09295bc1'
ARTIFACT_MANIFEST_HASH = 'f38a43dc9f5bd492dba35186855e4767e05c8ffffe93d4251035653c00f66fd9'
SOURCE_MANIFEST_HASH = '2cfe10d1dfb0223767e7913f41be39a345387e98b9cacc05cb9dc1e892b1949f'
HID_SOURCE_HASH = 'd9f9779c09d67cd5f8bc657692f174ffab14c469dfd06d560ac1899fa7e976b8'
BOHR = 0.52917721092
Z = {'H': 1, 'C': 6, 'N': 7, 'O': 8, 'Zn': 30}


class EvidenceError(ValueError):
    """Incomplete or mismatched evidence cannot be interpreted chemically."""


def raw(path):
    path = Path(path)
    data = path.read_bytes()
    if not data or len(data) != path.stat().st_size:
        raise EvidenceError(f'Incomplete/empty read: {path}')
    return data


def sha(path):
    return hashlib.sha256(raw(path)).hexdigest()


def load(path, expected=None):
    data = raw(path)
    if expected is not None and hashlib.sha256(data).hexdigest() != expected:
        raise EvidenceError(f'Hash mismatch: {path}')
    return json.loads(data)


def write_new(path, report):
    path = Path(path)
    with path.open('x') as stream:
        stream.write(json.dumps(report, indent=2, allow_nan=False) + '\n')


def check(condition, message):
    if not condition:
        raise EvidenceError(message)


def dist(a, b):
    return float(np.linalg.norm(a-b))


def rms(vectors):
    return float(np.sqrt(np.mean(np.sum(vectors*vectors, axis=1))))


def angle(a, b, c):
    u, v = a-b, c-b
    n = np.linalg.norm(u)*np.linalg.norm(v)
    if n <= 1e-12:
        return None
    return float(np.degrees(np.arccos(np.clip(np.dot(u,v)/n,-1,1))))


def torsion(points):
    a, b, c, d = np.asarray(points)
    axis = c-b
    if np.linalg.norm(axis) <= 1e-12:
        return None
    axis /= np.linalg.norm(axis)
    x = a-b - np.dot(a-b,axis)*axis
    y = d-c - np.dot(d-c,axis)*axis
    if np.linalg.norm(x)*np.linalg.norm(y) <= 1e-12:
        return None
    return float(np.degrees(np.arctan2(np.dot(np.cross(axis,x),y),np.dot(x,y))))


def signature(mol):
    return {'elements': [a.GetSymbol() for a in mol.GetAtoms()],
            'charges': [a.GetFormalCharge() for a in mol.GetAtoms()],
            'bonds': sorted((min(b.GetBeginAtomIdx(),b.GetEndAtomIdx()),
                            max(b.GetBeginAtomIdx(),b.GetEndAtomIdx()),
                            str(b.GetBondType())) for b in mol.GetBonds())}


def load_context(base=BASE):
    base = Path(base)
    req = load(base/'input.json', INPUT_HASH)
    policy = load(base/'optimization-policy.json', POLICY_HASH)
    check(sha(base/'qm_worker.py') == WORKER_HASH, 'Prepared worker changed')
    manifest = load(base/'artifact-manifest.json', ARTIFACT_MANIFEST_HASH)
    for name, evidence in manifest.items():
        check(Path(name).name == name and sha(base/name) == evidence['sha256'], f'Prepared artifact changed: {name}')
    source_manifest = load(base/'source-manifest.json', SOURCE_MANIFEST_HASH)
    # Bind only source artifacts needed here. Do not read a live runtime/source
    # merely to recertify unrelated paths that the prepared worker already copied.
    source = AUDIT/'inputs/1MNC.cif'
    ccd = AUDIT/'inputs/ccd/PLH.cif'
    for path in [source, ccd]:
        check(sha(path) == source_manifest[str(path)]['sha256'], f'Frozen structure source changed: {path}')
    baseline = load(base/'initial-geometry-diagnostics.json')
    mapping = load(base/'small-atom-map.json')
    ids, elements = req['atom_ids'], req['elements']
    check(len(ids) == 97 and len(set(ids)) == 97, 'Expected 97 unique IDs')
    check([a['stable_id'] for a in mapping] == ids, 'Atom map order differs')
    check(baseline['frozen_atom_ids'] == req['optimization']['freeze_atom_ids'], 'Cap map differs')
    check(len(baseline['coordination']) == 5 and len(baseline['zinc_centered_angles']) == 10, 'Coordination manifest differs')
    index = {a:i for i,a in enumerate(ids)}
    bonds = {tuple(sorted(row['atom_ids'])) for row in baseline['ligand_bond_graph']}
    # Native standard HID connectivity plus exactly documented MCPB methyl cuts.
    hid_xml = ROOT/'.venv/lib/python3.12/site-packages/openmm/app/data/amber14/protein.ff14SB.xml'
    hid_source = raw(hid_xml)
    check(hashlib.sha256(hid_source).hexdigest() == HID_SOURCE_HASH, 'Reviewed native HID source changed')
    template = ET.fromstring(hid_source).find("./Residues/Residue[@name='HID']")
    check(template is not None, 'Native HID template missing')
    for resid in ['218','222','228']:
        atoms = [a for a in mapping if a['resid'] == resid and a['resname'] == 'HID']
        retained = {a['name']:a['stable_id'] for a in atoms if a['role'] == 'retained'}
        caps = {a['name']:a['stable_id'] for a in atoms if a['role'] == 'native MCPB cap'}
        check(set(caps) == {'H1','H2','H3','CH3'}, f'Native cap inventory differs for His{resid}')
        retained['CA'] = caps['CH3']
        for bond in template.findall('Bond'):
            a, b = bond.attrib['atomName1'], bond.attrib['atomName2']
            if a in retained and b in retained:
                bonds.add(tuple(sorted((retained[a],retained[b]))))
        for h in ['H1','H2','H3']:
            bonds.add(tuple(sorted((caps['CH3'],caps[h]))))
    check(len(bonds) == 96, 'Expected complete 96-bond nonmetal cluster graph')
    neighbors = {a:set() for a in ids}
    for a,b in bonds:
        neighbors[a].add(b);neighbors[b].add(a)
    parents = {}
    for a,e in zip(ids,elements,strict=True):
        if e == 'H':
            check(len(neighbors[a]) == 1, f'Hydrogen parent is not unique: {a}')
            parents[a] = next(iter(neighbors[a]))
    check(len(parents) == 50, 'Hydrogen inventory differs')
    check(not any(p == 'E:280::PLH:O1' or p.endswith(':HID:NE2') for p in parents.values()), 'Donor gained a covalent hydrogen')
    for n in ['218','222','228']:
        check(parents[f'A:{n}::HID:HD1'] == f'A:{n}::HID:ND1', 'HID tautomer map differs')
    check(parents['E:280::PLH:H24'] == 'E:280::PLH:N1', 'PLH N1 proton map differs')
    ligand_ids = [a for a in ids if ':PLH:' in a]
    ligand = Chem.MolFromMolBlock(raw(base/'PLH-declared-state.sdf').decode().split('$$$$')[0],removeHs=False,sanitize=True)
    check(ligand is not None and ligand.GetNumAtoms() == len(ligand_ids) == 51, 'PLH source graph cannot be read')
    mol2_text = raw(base/'PLH.mol2').decode().split('@<TRIPOS>ATOM\n')[1].split('@<TRIPOS>BOND')[0]
    mol2_names = [r[1] for r in (line.split() for line in mol2_text.splitlines() if line.strip())]
    check([a.rsplit(':',1)[1] for a in ligand_ids] == mol2_names, 'SDF/MOL2 name order is not verified')
    ligand_graph = signature(ligand)
    check(ligand_graph['elements'] == [elements[index[a]] for a in ligand_ids], 'Ligand element order differs')
    check(sum(ligand_graph['charges']) == -1, 'Ligand formal charge differs')
    source_pairs = {tuple(sorted((ligand_ids[b.GetBeginAtomIdx()],ligand_ids[b.GetEndAtomIdx()]))) for b in ligand.GetBonds()}
    check(source_pairs == {tuple(sorted(r['atom_ids'])) for r in baseline['ligand_bond_graph']}, 'Ligand source graph differs')
    check(np.max(np.abs(ligand.GetConformer().GetPositions()-np.asarray(req['coords_angstrom'])[[index[a] for a in ligand_ids]])) <= 0.00051,
          'Ligand SDF order/initial coordinate alignment differs')
    # Cross-check bond orders against the native MOL2 graph (SDF aromatic bonds
    # are sanitized to aromatic, so no arbitrary Kekule representation is used).
    for row in baseline['ligand_bond_graph']:
        a,b = [ligand_ids.index(x) for x in row['atom_ids']]
        bond = ligand.GetBondBetweenAtoms(a,b)
        actual = 'ar' if bond.GetIsAromatic() else str(int(bond.GetBondTypeAsDouble()))
        check(actual == row['native_mol2_order'], 'Ligand native bond order differs')
    context = {'base':base,'request':req,'policy':policy,'baseline':baseline,'mapping':mapping,
               'ids':ids,'elements':elements,'index':index,'bonds':bonds,'parents':parents,
               'ligand_ids':ligand_ids,'ligand':ligand,'ligand_signature':ligand_graph,
               'neighbors':neighbors,'hid_source_sha256':hashlib.sha256(hid_source).hexdigest(),
               'hid_source_path':str(hid_xml),'source_sha256':sha(source)}
    context['omitted_protein'] = omitted_protein(context, source)
    return context


def omitted_protein(ctx, source):
    block = gemmi.cif.read_file(str(source)).sole_block()
    columns = ['group_PDB','type_symbol','label_atom_id','label_comp_id','label_asym_id',
               'auth_seq_id','pdbx_PDB_ins_code','Cartn_x','Cartn_y','Cartn_z','label_alt_id',
               'pdbx_PDB_model_num','label_seq_id']
    table = block.find('_atom_site.',columns)
    check(len(table)>0, 'Frozen mmCIF atom table missing')
    represented = set()
    for atom,e in zip(ctx['mapping'],ctx['elements'],strict=True):
        if e == 'H':
            continue
        origin = atom['source_or_placed_id'] if atom['role'] == 'native MCPB cap' else atom['stable_id']
        represented.add(origin.replace(':HID:',':HIS:'))
    atoms, all_atoms = [], {}
    for row in table:
        group,e,name,res,chain,seq,icode,x,y,z,alt,model,polyseq = map(str,row)
        if model != '1' or e == 'H':
            continue
        check(alt in ['.','?',''], 'Frozen benchmark unexpectedly has alternate atoms; selection requires review')
        atom_id = f'{chain}:{seq}:{"" if icode in [".","?"] else icode}:{res}:{name}'
        check(atom_id not in all_atoms, 'Duplicate source heavy atom identity')
        entry = {'source_atom_id':atom_id,'element':e,'coords_angstrom':[float(x),float(y),float(z)]}
        all_atoms[atom_id] = entry
        if polyseq not in ['.','?'] and atom_id not in represented:
            atoms.append(entry)
    for atom,e in zip(ctx['mapping'],ctx['elements'],strict=True):
        if e == 'H':
            continue
        origin = atom['source_or_placed_id'] if atom['role'] == 'native MCPB cap' else atom['stable_id']
        origin = origin.replace(':HID:',':HIS:')
        check(origin in all_atoms, f'Retained heavy source identity missing: {origin}')
        check(np.max(np.abs(np.asarray(all_atoms[origin]['coords_angstrom'])-np.asarray(atom['xyz']))) < 0.00051,
              f'Retained heavy source coordinate mismatch: {origin}')
    check(len(atoms)>1000, 'Omitted protein inventory unexpectedly incomplete')
    return atoms


def read_endpoint(ctx, folder):
    folder = Path(folder)
    if not (folder/'result.json').is_file():
        return None
    result = load(folder/'result.json')
    check(result.get('schema_version') == 1 and result.get('status') == 'completed' and result.get('accepted') is True,
          'Endpoint worker did not report a completed numerical result')
    req = ctx['request']
    for k,expected in [('atom_ids',req['atom_ids']),('elements',req['elements']),('charge',1),('spin_2S',0),
                       ('explicit_electron_count',372),('nao',778),('input_sha256',INPUT_HASH),('worker_sha256',WORKER_HASH)]:
        check(result.get(k) == expected, f'Endpoint metadata differs: {k}')
    check(sha(folder/'input.json') == INPUT_HASH, 'Actual worker input differs')
    check(result.get('scf',{}).get('converged') is True, 'Final SCF not converged')
    check(result['scf'].get('conv_tol_hartree') == req['scf_conv_tol'] and result['scf'].get('conv_tol_grad') == req['scf_conv_tol_grad'],
          'Final SCF thresholds differ')
    opt = result.get('optimization',{})
    check(opt.get('requested') is True and opt.get('converged') is True, 'Geometry optimization not converged')
    expected_opt = load(ctx['base']/'validated-input.json')['optimization']
    check(opt.get('settings') == expected_opt, 'Actual optimizer settings differ')
    check(opt.get('constraint_atom_indices_zero_based') == [ctx['index'][a] for a in expected_opt['freeze_atom_ids']],
          'Actual Cartesian cap constraint indices differ')
    constraint_bytes=('$freeze\nxyz '+','.join(str(ctx['index'][a]+1) for a in expected_opt['freeze_atom_ids'])+'\n').encode()
    check(raw(folder/'constraints.txt') == constraint_bytes and opt.get('constraint_file_sha256') == hashlib.sha256(constraint_bytes).hexdigest(),
          'Actual constraint file differs from the three-carbon-only request')
    check(opt.get('dihedral_constraints_satisfied') is True and opt.get('final_dihedral_constraints') == [],
          'Unexpected or unverified dihedral constraints')
    check(opt.get('steps') and all(s.get('scf_converged') is True and math.isfinite(s.get('energy_hartree',float('nan'))) for s in opt['steps']),
          'Unconverged/nonfinite optimization step evidence')
    check(math.isfinite(result.get('energy_hartree',float('nan'))), 'Nonfinite endpoint energy')
    density_count=result.get('density_electron_count',float('nan'))
    check(math.isfinite(density_count) and abs(density_count-372)<1e-6, 'Endpoint density electron count differs')
    check(result.get('units',{}).get('bohr_to_angstrom') == BOHR and
          result['units'].get('gradient_hartree_per_bohr') == 'hartree/bohr', 'Endpoint derivative/coordinate units differ')
    for name,key in [('arrays.npz','arrays_sha256'),('scf.chk','checkpoint_sha256'),('resolved-method.json','resolved_method_file_sha256')]:
        check(result.get({'arrays.npz':'arrays_file','scf.chk':'checkpoint_file','resolved-method.json':'resolved_method_file'}[name]) == name,
              'Unexpected result artifact layout')
        check(sha(folder/name) == result.get(key), f'Endpoint artifact checksum failed: {name}')
    reference_path = Path(req['initial_checkpoint']['result_path']).parent/'resolved-method.json'
    reference_result = load(Path(req['initial_checkpoint']['result_path']),req['initial_checkpoint']['result_sha256'])
    reference = load(reference_path,reference_result['resolved_method_file_sha256'])
    method = load(folder/'resolved-method.json')
    for k in ['method','basis_sha256','ecp_sha256','spherical_basis','symmetry','density_fit',
              'functional_resolved_sha256','grid_level','auxbasis_sha256']:
        check(method.get(k) == reference.get(k), f'Resolved endpoint method changed: {k}')
    runtime = load(ctx['base']/'source-manifest.json',SOURCE_MANIFEST_HASH)[str(ROOT/'.tools/qm-parallel/dynamol-qm-runtime.json')]['sha256']
    check(result.get('runtime_manifest',{}).get('sha256') == runtime, 'Declared runtime manifest differs')
    check(sha(folder/'runtime-manifest.json') == runtime, 'Copied runtime manifest differs')
    check(result.get('actual_pyscf_threads') == 2 and result.get('resources',{}).get('max_wall_seconds') == 21600,
          'Actual resource settings differ')
    check(result['resources'].get('max_memory_mb') == 8000 and result.get('versions',{}).get('pyscf') == '2.14.0',
          'Requested memory/provider version differs')
    arrays = {}
    with np.load(folder/'arrays.npz', allow_pickle=False) as payload:
        expected_arrays=['coords_bohr','coords_angstrom','gradient_hartree_per_bohr','atomic_numbers','effective_nuclear_charges']
        check(set(payload.files) == set(expected_arrays), 'Missing or unrequested endpoint arrays')
        for key in expected_arrays:
            check(key in payload, f'Missing endpoint array: {key}')
            value = np.asarray(payload[key])
            shape = (97,3) if key in ['coords_bohr','coords_angstrom','gradient_hartree_per_bohr'] else (97,)
            check(value.shape == shape and np.issubdtype(value.dtype,np.number) and np.isfinite(value).all(), f'Malformed/nonfinite endpoint array: {key}')
            arrays[key] = value.copy()
    numbers = np.asarray([Z[e] for e in ctx['elements']])
    check(np.array_equal(arrays['atomic_numbers'],numbers) and np.array_equal(arrays['effective_nuclear_charges'],numbers),
          'Atomic number/effective nuclear charge map differs')
    check(np.max(np.abs(arrays['coords_bohr']*BOHR-arrays['coords_angstrom'])) <= 1e-12, 'Bohr/Angstrom coordinate arrays disagree')
    return {'folder':folder,'result':result,'arrays':arrays,'hashes':{
        name:sha(folder/name) for name in ['result.json','input.json','arrays.npz','scf.chk','resolved-method.json','runtime-manifest.json']}}


def chemistry_report(ctx, xyz, gradient):
    xyz, gradient = np.asarray(xyz,dtype=float), np.asarray(gradient,dtype=float)
    check(xyz.shape == (97,3) and gradient.shape == (97,3) and np.isfinite(xyz).all() and np.isfinite(gradient).all(),
          'Endpoint geometry/gradient shape or finite-value check failed')
    ids, idx, elements = ctx['ids'],ctx['index'],ctx['elements']
    initial = np.asarray(ctx['request']['coords_angstrom'])
    p = ctx['policy']['endpoint_checks']
    alarms = []
    def alarm(code, **detail):
        alarms.append({'code':code, **detail})
    coords = dict(zip(ids,xyz,strict=True))
    frozen = ctx['request']['optimization']['freeze_atom_ids']
    free = [a for a in ids if a not in frozen]
    gradient_rows = [{'atom_id':a,'role':'frozen_cap' if a in frozen else 'free',
                      'gradient_hartree_per_bohr':gradient[idx[a]].tolist(),
                      'force_hartree_per_bohr':(-gradient[idx[a]]).tolist(),
                      'gradient_norm_hartree_per_bohr':float(np.linalg.norm(gradient[idx[a]]))} for a in ids]
    gradient_stats = {}
    for name,atoms in [('full',ids),('free',free),('frozen_caps',frozen)]:
        g = gradient[[idx[a] for a in atoms]]
        gradient_stats[name] = {'atoms':len(atoms),'rms_atomic_vector':rms(g),
            'max_atomic_vector':float(np.max(np.linalg.norm(g,axis=1))),
            'rms_component':float(np.sqrt(np.mean(g*g))),'max_absolute_component':float(np.max(np.abs(g)))}
    if gradient_stats['free']['rms_atomic_vector'] > p['independent_free_DF_gradient_rms_vector_hartree_per_bohr']:
        alarm('free_gradient_rms_exceeded')
    if gradient_stats['free']['max_atomic_vector'] > p['independent_free_DF_gradient_max_vector_hartree_per_bohr']:
        alarm('free_gradient_max_exceeded')
    displacement = np.linalg.norm(xyz-initial,axis=1)
    cap_displacements = [{'atom_id':a,'angstrom':float(displacement[idx[a]]),'bohr':float(displacement[idx[a]]/BOHR)} for a in frozen]
    if max(r['bohr'] for r in cap_displacements) > p['cap_displacement_max_bohr']:
        alarm('frozen_cap_displacement_exceeded')
    cap_pairs = [{'atom_ids':[a,b],'initial_angstrom':dist(initial[idx[a]],initial[idx[b]]),
                  'final_angstrom':dist(coords[a],coords[b])} for a,b in itertools.combinations(frozen,2)]
    coordination = []
    donor_ids = []
    for row in ctx['baseline']['coordination']:
        a,b = row['atom_ids'];donor_ids.append(b)
        d = dist(coords[a],coords[b]);lo=p['coordination_alarms_angstrom']['lower'];hi=p['coordination_alarms_angstrom']['upper']
        coordination.append({**row,'final_distance_angstrom':d,'change_angstrom':d-row['initial_distance_angstrom']})
        if not lo <= d <= hi:
            alarm('named_coordination_distance_outside_bounds',atom_ids=[a,b],distance_angstrom=d)
    zn = ctx['baseline']['coordination'][0]['atom_ids'][0]
    new_neighbors = []
    for a,e in zip(ids,elements,strict=True):
        if e in ['N','O'] and a not in donor_ids:
            d=dist(coords[zn],coords[a]);entry={'atom_id':a,'initial_angstrom':dist(initial[idx[zn]],initial[idx[a]]),'final_angstrom':d}
            new_neighbors.append(entry)
            if d < p['new_N_O_neighbor_alarm_angstrom']:
                alarm('undeclared_N_O_near_zinc',**entry)
    angles = []
    for row in ctx['baseline']['zinc_centered_angles']:
        value = angle(*(coords[a] for a in row['atom_ids']))
        delta = None if value is None else value-row['initial_degrees']
        angles.append({**row,'final_degrees':value,'change_degrees':delta})
        if delta is None or abs(delta) > p['zinc_centered_angle_change_review_degrees']:
            alarm('coordination_angle_change',atom_ids=row['atom_ids'],change_degrees=delta)
    bond_rows = []
    for a,b in sorted(ctx['bonds']):
        is_h = elements[idx[a]] == 'H' or elements[idx[b]] == 'H'
        lo,hi = p['hydrogen_parent_bond_alarm_angstrom'] if is_h else p['C_N_O_covalent_bond_alarm_angstrom']
        d=dist(coords[a],coords[b]);bond_rows.append({'atom_ids':[a,b],'initial_angstrom':dist(initial[idx[a]],initial[idx[b]]),'final_angstrom':d,'hydrogen_bonded_atom':is_h})
        if not lo <= d <= hi:
            alarm('covalent_bond_geometry',atom_ids=[a,b],distance_angstrom=d)
    heavy = [a for a,e in zip(ids,elements,strict=True) if e != 'H']
    proton_rows = []
    for h,parent in sorted(ctx['parents'].items()):
        d=dist(coords[h],coords[parent])
        others=sorted((dist(coords[h],coords[a]),a) for a in heavy if a != parent and a != zn)
        near,other=others[0]
        competing=[{'atom_id':a,'distance_angstrom':v} for v,a in others if v <= p['hydrogen_parent_bond_alarm_angstrom'][1]]
        proton_rows.append({'hydrogen_id':h,'expected_parent_id':parent,'parent_distance_angstrom':d,
                            'nearest_other_heavy_atom_id':other,'nearest_other_distance_angstrom':near,'competing_covalent_range_neighbors':competing})
        if near <= d or competing:
            alarm('ambiguous_or_transferred_hydrogen',hydrogen_id=h,parent_id=parent,other_id=other)
    forbidden_donors = [a for a in donor_ids if a != 'E:280::PLH:O2']
    for donor in forbidden_donors:
        if any(dist(coords[h],coords[donor]) <= p['hydrogen_parent_bond_alarm_angstrom'][1] for h in ctx['parents']):
            alarm('unprotonated_donor_has_close_hydrogen',donor_id=donor)
    # Reassign stereochemistry from actual coordinates; erase inherited 3D tags
    # and properties so the source result cannot mask a geometry inversion.
    ligand = Chem.Mol(ctx['ligand'])
    for atom in ligand.GetAtoms():
        atom.SetChiralTag(Chem.ChiralType.CHI_UNSPECIFIED)
        if atom.HasProp('_CIPCode'):atom.ClearProp('_CIPCode')
    conformer = ligand.GetConformer()
    for i,a in enumerate(ctx['ligand_ids']):
        conformer.SetAtomPosition(i,tuple(float(v) for v in coords[a]))
    Chem.AssignStereochemistryFrom3D(ligand,replaceExistingTags=True)
    Chem.AssignStereochemistry(ligand,cleanIt=True,force=True)
    check(signature(ligand) == ctx['ligand_signature'], 'Stereo assignment unexpectedly changed ligand graph')
    stereo = []
    for row in ctx['baseline']['stereochemistry']:
        a,b,c,d = [coords[x] for x in row['ordered_neighbor_ids']]
        determinant=float(np.linalg.det(np.vstack((a-d,b-d,c-d))))
        atom=ligand.GetAtomWithIdx(ctx['ligand_ids'].index(row['atom_id']))
        cip=atom.GetProp('_CIPCode') if atom.HasProp('_CIPCode') else None
        same_sign=determinant*row['initial_signed_determinant_angstrom3'] > 0
        stereo.append({**row,'final_signed_determinant_angstrom3':determinant,'sign_preserved':same_sign,'final_graph_based_CIP':cip})
        if not same_sign or abs(determinant)<.1 or cip != row['verified_source_CIP']:
            alarm('ligand_stereochemistry_changed_or_ambiguous',atom_id=row['atom_id'],CIP=cip)
    contacts = []
    coordination_pairs = {tuple(sorted(r['atom_ids'])) for r in coordination}
    for a,b in itertools.combinations(heavy,2):
        pair=tuple(sorted((a,b)))
        if pair in ctx['bonds'] or pair in coordination_pairs:
            continue
        d=dist(coords[a],coords[b])
        if d < p['nonbonded_heavy_contact_alarm_angstrom']:
            contacts.append({'atom_ids':[a,b],'distance_angstrom':d})
            alarm('nonbonded_heavy_collision',atom_ids=[a,b],distance_angstrom=d)
    retained_heavy=[a for a in heavy if a not in frozen]
    ligand_heavy=[a for a in retained_heavy if ':PLH:' in a]
    ligand_core=['E:280::PLH:'+n for n in ['C1','N1','O1','O2']]
    core=[zn]+ligand_core+[a for a in retained_heavy if ':HID:' in a and a.rsplit(':',1)[1] in ['CG','ND1','CE1','NE2','CD2']]
    groups={'retained_heavy':retained_heavy,'Zn':[zn],'donors':donor_ids,'site_core':core,
            'whole_PLH_heavy':ligand_heavy,'distal_PLH_heavy':[a for a in ligand_heavy if a not in ligand_core]}
    movements={}
    for name,atoms in groups.items():
        v=(xyz-initial)[[idx[a] for a in atoms]]
        movements[name]={'atom_ids':atoms,'source_frame_rmsd_angstrom':rms(v),'max_displacement_angstrom':float(np.max(np.linalg.norm(v,axis=1)))}
    if movements['retained_heavy']['source_frame_rmsd_angstrom'] > p['retained_heavy_frame_rmsd_review_angstrom']:
        alarm('retained_heavy_source_frame_rmsd_exceeded')
    if movements['retained_heavy']['max_displacement_angstrom'] > p['retained_heavy_max_displacement_review_angstrom']:
        alarm('retained_heavy_max_displacement_exceeded')
    torsions=set()
    for b,c in ctx['bonds']:
        if b not in ligand_heavy or c not in ligand_heavy:
            continue
        for a in ctx['neighbors'][b]-{c}:
            for d in ctx['neighbors'][c]-{b}:
                path=(a,b,c,d)
                if len(set(path))==4 and set(path)<=set(ligand_heavy):
                    torsions.add(min(path,path[::-1]))
    torsion_rows=[]
    for path in sorted(torsions):
        start=torsion([initial[idx[a]] for a in path]);end=torsion([coords[a] for a in path])
        delta=None if start is None or end is None else (end-start+180)%360-180
        torsion_rows.append({'atom_ids':list(path),'initial_degrees':start,'final_degrees':end,'periodic_change_degrees':delta})
        if end is None:
            alarm('undefined_ligand_proper_torsion',atom_ids=list(path))
    pocket_coords=np.asarray([a['coords_angstrom'] for a in ctx['omitted_protein']])
    pocket_contacts=[];nearest_contacts=[]
    for a in retained_heavy:
        before=np.linalg.norm(pocket_coords-initial[idx[a]],axis=1)
        after=np.linalg.norm(pocket_coords-coords[a],axis=1)
        nearest=int(np.argmin(after))
        nearest_contacts.append({'cluster_atom_id':a,'source_protein_atom_id':ctx['omitted_protein'][nearest]['source_atom_id'],
                                  'initial_angstrom':float(before[nearest]),'final_angstrom':float(after[nearest])})
        for j in np.where((before<=4.0)|(after<=4.0))[0]:
            row={'cluster_atom_id':a,'source_protein_atom_id':ctx['omitted_protein'][int(j)]['source_atom_id'],
                 'initial_angstrom':float(before[j]),'final_angstrom':float(after[j])}
            pocket_contacts.append(row)
            if after[j] < p['nonbonded_heavy_contact_alarm_angstrom']:
                alarm('omitted_protein_heavy_collision',**row)
    return {'status':'endpoint_geometry_numerics_clear_for_conventional_check' if not alarms else 'endpoint_requires_model_review',
        'endpoint_checks_passed':not alarms,'model_accuracy_validated':False,'full_preparation_ready':False,
        'not_claimed':'No full unconstrained minimum, accepted force field, correct protonation equilibrium or dynamics validation.',
        'alarms':alarms,'gradient_units':'hartree/bohr','gradient_statistics':gradient_stats,'atomic_gradient_rows':gradient_rows,
        'cap_displacements':cap_displacements,'cap_pair_distances':cap_pairs,'coordination':coordination,'zinc_centered_angles':angles,
        'non_donor_N_O_distances':new_neighbors,'covalent_bond_distances':bond_rows,'hydrogen_parent_checks':proton_rows,
        'ligand_stereochemistry':stereo,'intracluster_nonbonded_heavy_collisions':contacts,
        'source_frame_displacements':[{'atom_id':a,'vector_angstrom':(xyz[idx[a]]-initial[idx[a]]).tolist(),'norm_angstrom':float(displacement[idx[a]])} for a in ids],
        'group_displacements':movements,'proper_ligand_heavy_torsions':torsion_rows,
        'omitted_protein_contacts':{'source_sha256':ctx['source_sha256'],'omitted_heavy_atoms':len(ctx['omitted_protein']),
            'coordinate_frame':'Original frame; no rigid-body alignment or pocket relaxation.',
            'contact_reporting_cutoff_angstrom':4.0,'reporting_cutoff_is_not_an_acceptance_threshold':True,
            'nearest_contact_per_retained_heavy_atom':nearest_contacts,'contacts_within_cutoff_at_either_geometry':pocket_contacts}}


def conventional_request(ctx, endpoint, analysis):
    check(analysis.get('endpoint_checks_passed') is True and not analysis.get('alarms'), 'Conventional request requires all endpoint checks')
    request=dict(ctx['request'])
    request.pop('coords_angstrom');request['coords_bohr']=endpoint['arrays']['coords_bohr'].tolist()
    request.update(density_fit=False,optimization={'enabled':False},operations=['gradient'],
        initial_checkpoint={'result_path':str(endpoint['folder']/'result.json'),'result_sha256':endpoint['hashes']['result.json'],
            'checkpoint_path':str(endpoint['folder']/'scf.chk'),'checkpoint_sha256':endpoint['hashes']['scf.chk']})
    request.pop('auxbasis',None)
    # No launch. The output explicitly preserves actual final Bohr coordinates,
    # and the worker independently verifies the same-state checkpoint at launch.
    return request


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--endpoint',type=Path,required=True)
    parser.add_argument('--output',type=Path,required=True)
    parser.add_argument('--write-conventional-input',action='store_true')
    args=parser.parse_args()
    args.output.mkdir(parents=True,exist_ok=False)
    report={'checked_unix':time.time(),'QM_launched':False,'model_accuracy_validated':False,'full_preparation_ready':False,
            'analyzer_sha256':sha(Path(__file__)),'versions':{'numpy':np.__version__,'gemmi':gemmi.__version__,'rdkit':rdBase.rdkitVersion}}
    try:
        ctx=load_context()
        endpoint=read_endpoint(ctx,args.endpoint)
        if endpoint is None:
            report['status']='endpoint_not_present_no_result_read'
        else:
            analysis=chemistry_report(ctx,endpoint['arrays']['coords_angstrom'],endpoint['arrays']['gradient_hartree_per_bohr'])
            report.update(analysis=analysis,status=analysis['status'],endpoint_hashes=endpoint['hashes'],
                endpoint_energy_hartree=endpoint['result']['energy_hartree'],input_sha256=INPUT_HASH,policy_sha256=POLICY_HASH,
                native_HID_connectivity_source={'path':ctx['hid_source_path'],'sha256':ctx['hid_source_sha256']})
            if args.write_conventional_input and analysis['endpoint_checks_passed']:
                path=args.output/'conventional-gradient-input.json'
                write_new(path,conventional_request(ctx,endpoint,analysis));report['conventional_input_sha256']=sha(path)
            elif args.write_conventional_input:
                report['conventional_input_written']=False
                report['reason']='Endpoint alarms require model review; no downstream request generated.'
    except Exception as exc:
        report.update(status='evidence_or_analysis_failed',error_type=type(exc).__name__,error=str(exc))
    write_new(args.output/'result.json',report)
    print(json.dumps({k:v for k,v in report.items() if k not in ['analysis','endpoint_hashes']},indent=2))
    if report['status']=='evidence_or_analysis_failed':
        raise SystemExit(1)


if __name__ == '__main__':
    main()
