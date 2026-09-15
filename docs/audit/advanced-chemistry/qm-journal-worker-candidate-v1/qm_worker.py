#!/usr/bin/env python3
"""Isolated, versioned PySCF worker; never imports into the app Python runtime.

Usage: private-python backend/qm_worker.py input.json FRESH_OUTPUT_DIRECTORY
An accepted result means the requested numerical calculation completed its
checks. It is not a validation of a force field, oxidation state, or QM model.
"""
import hashlib
import importlib.metadata
import json
import math
import os
from pathlib import Path
import platform
import signal
import sys
import time
import traceback

SCHEMA_VERSION = 1
OPT_DEFAULTS = {
    'convergence_energy': 1e-6, 'convergence_grms': 3e-4,
    'convergence_gmax': 4.5e-4, 'convergence_drms': 1.2e-3,
    'convergence_dmax': 1.8e-3,
}


def _json(path, value):
    path = Path(path)
    temporary = path.with_suffix(path.suffix + '.tmp')
    temporary.write_text(json.dumps(value, indent=2, allow_nan=False) + '\n')
    temporary.replace(path)


def _plain(value):
    if hasattr(value, 'tolist'):
        return value.tolist()
    if isinstance(value, dict):
        return {str(k): _plain(v) for k, v in value.items()}
    if isinstance(value, (list, tuple)):
        return [_plain(v) for v in value]
    return value


def _hash(value):
    return hashlib.sha256(json.dumps(_plain(value), sort_keys=True,
                                    separators=(',', ':'), allow_nan=False).encode()).hexdigest()


def _integer(value, name, minimum, maximum):
    if isinstance(value, bool) or not isinstance(value, int) or not minimum <= value <= maximum:
        raise ValueError(f'{name} must be an integer from {minimum} to {maximum}')
    return value


def _positive(value, name, maximum=1.0):
    if isinstance(value, bool) or not isinstance(value, (int, float)) or not math.isfinite(value) or not 0 < value <= maximum:
        raise ValueError(f'{name} must be finite and in (0, {maximum}]')
    return float(value)


def _sha256_string(value):
    return isinstance(value,str) and len(value)==64 and all(c in '0123456789abcdef' for c in value)


def dihedral_degrees(coordinates):
    """Signed A-B-C-D angle, with explicit rejection of undefined planes.

    Coordinates may be in any consistent length unit. A linear adjacent angle
    is singular; near-linearity is also rejected where geomeTRIC's derivative
    falls to zero (sine-squared <= 1e-6).
    """
    import numpy as np
    xyz=np.asarray(coordinates,dtype=float)
    if xyz.shape!=(4,3) or not np.isfinite(xyz).all():
        raise ValueError('Dihedral requires four finite Cartesian coordinates')
    left=xyz[0]-xyz[1];axis=xyz[2]-xyz[1];right=xyz[3]-xyz[2]
    lengths=np.array([np.linalg.norm(left),np.linalg.norm(axis),np.linalg.norm(right)])
    if lengths.min()<=1e-8:raise ValueError('Singular dihedral: coincident defining atoms')
    axis=axis/lengths[1]
    first=left-np.dot(left,axis)*axis
    second=right-np.dot(right,axis)*axis
    if np.linalg.norm(first)/lengths[0]<=1e-3 or np.linalg.norm(second)/lengths[2]<=1e-3:
        raise ValueError('Singular dihedral: a defining angle is linear or nearly linear')
    return math.degrees(math.atan2(float(np.dot(np.cross(axis,first),second)),float(np.dot(first,second))))


def _dihedral_rows(coordinates,atom_ids,constraints,tolerance):
    import numpy as np
    xyz=np.asarray(coordinates,dtype=float)
    rows=[]
    for constraint in constraints:
        indices=[atom_ids.index(value) for value in constraint['atom_ids']]
        measured=dihedral_degrees(xyz[indices])
        target=constraint['target_degrees']
        error=math.remainder(measured-target,360.0)
        rows.append({'atom_ids':constraint['atom_ids'],'atom_indices_zero_based':indices,
            'target_degrees':target,'canonical_target_degrees':math.remainder(target,360.0),
            'measured_degrees':measured,'periodic_error_degrees':error,
            'absolute_periodic_error_degrees':abs(error),'tolerance_degrees':tolerance,
            'target_satisfied':abs(error)<=tolerance})
    return rows


def capture_runtime_manifest(out):
    """Copy declared build provenance; this is not loaded-binary attestation."""
    source=Path(sys.prefix)/'dynamol-qm-runtime.json'
    if not source.exists():return {'present':False}
    if source.is_symlink() or not source.is_file() or source.stat().st_size>2*1024*1024:
        raise ValueError('Installed QM runtime manifest must be a regular JSON file below 2 MB')
    raw=source.read_bytes()
    if not isinstance(json.loads(raw),dict):raise ValueError('Installed QM runtime manifest must contain a JSON object')
    (out/'runtime-manifest.json').write_bytes(raw)
    return {'present':True,'source_path':str(source),'copied_file':'runtime-manifest.json',
        'sha256':hashlib.sha256(raw).hexdigest(),
        'meaning':'Declared installed build provenance copied verbatim; individual loaded binary hashes are not verified by this worker.'}


def validate_input(data):
    """Reject ambiguous units/methods and ignored options before native work."""
    if not isinstance(data, dict):
        raise ValueError('Input must be a JSON object')
    allowed = {'schema_version','atom_ids','elements','coords_angstrom','coords_bohr',
               'charge','spin','method','basis','ecp','functional','grid_level',
               'density_fit','auxbasis','threads','max_memory_mb','operations',
               'esp_points_bohr','esp_batch_size','esp_crosscheck_points',
               'max_scf_cycles','scf_conv_tol','scf_conv_tol_grad','max_cphf_cycles',
               'optimization','max_wall_seconds','initial_checkpoint'}
    if set(data) - allowed:
        raise ValueError(f'Unknown input options: {sorted(set(data) - allowed)}')
    if data.get('schema_version', 1) != SCHEMA_VERSION:
        raise ValueError('Unsupported schema_version')
    ids, elements = data.get('atom_ids'), data.get('elements')
    if not isinstance(ids, list) or not ids or len(ids) > 1000 or any(not isinstance(x, str) or not x for x in ids) or len(set(ids)) != len(ids):
        raise ValueError('atom_ids must be 1–1000 unique nonempty stable strings')
    if not isinstance(elements, list) or len(elements) != len(ids) or any(not isinstance(x, str) or not x.isalpha() or len(x)>2 for x in elements):
        raise ValueError('elements must contain one explicit chemical symbol per atom')
    if ('coords_angstrom' in data) == ('coords_bohr' in data):
        raise ValueError('Supply exactly one of coords_angstrom or coords_bohr')
    coords = data.get('coords_angstrom', data.get('coords_bohr'))
    if not isinstance(coords, list) or len(coords) != len(ids) or any(not isinstance(x,list) or len(x)!=3 for x in coords):
        raise ValueError('Coordinates must be an N×3 array in explicitly named units')
    for row in coords:
        if any(isinstance(x,bool) or not isinstance(x,(float,int)) or not math.isfinite(x) for x in row):
            raise ValueError('Coordinates must contain finite numbers')
    charge = _integer(data.get('charge'), 'charge', -100, 100)
    spin = _integer(data.get('spin'), 'spin', 0, 100)
    if spin != 0:
        raise ValueError('This worker currently supports closed-shell spin=0 only; no spin state is guessed')
    method = data.get('method')
    if method not in ('RHF','RKS'):
        raise ValueError('method must be explicitly RHF or RKS')
    basis = data.get('basis')
    for name, value in [('basis',basis),('ecp',data.get('ecp')),('auxbasis',data.get('auxbasis'))]:
        if value is None and name != 'basis':
            continue
        if not (isinstance(value,str) and value.strip() or isinstance(value,dict) and value and all(isinstance(k,str) and isinstance(v,str) and v.strip() for k,v in value.items())):
            raise ValueError(f'{name} must be a named basis/ECP string or element-to-name dictionary')
        if isinstance(value,dict):
            keys=[k.lower() for k in value]
            if len(set(keys))!=len(keys) or set(keys)-{x.lower() for x in elements}:
                raise ValueError(f'{name} element map contains duplicate or absent elements')
    functional = data.get('functional')
    if method=='RKS' and (not isinstance(functional,str) or not functional.strip()):
        raise ValueError('RKS requires an explicit functional; use B3LYPG for the Gaussian B3LYP reference')
    if method=='RKS' and functional.replace('-','').replace(' ','').upper()=='B3LYP':
        raise ValueError('Ambiguous B3LYP alias: use explicit B3LYPG or B3LYP5')
    if method=='RHF' and functional is not None:
        raise ValueError('functional is not applicable to RHF')
    operations = data.get('operations')
    if not isinstance(operations,list) or not operations or any(x not in ('energy','gradient','hessian','esp') for x in operations) or len(set(operations)) != len(operations):
        raise ValueError('operations must be a nonempty unique list of energy, gradient, hessian, esp')
    if ('esp' in operations) != ('esp_points_bohr' in data):
        raise ValueError('esp requires explicit esp_points_bohr, which is accepted only for the esp operation')
    if 'esp' in operations:
        points=data['esp_points_bohr']
        if not isinstance(points,list) or not 1<=len(points)<=1000000 or any(not isinstance(p,list) or len(p)!=3 or any(isinstance(x,bool) or not isinstance(x,(int,float)) or not math.isfinite(x) for x in p) for p in points):
            raise ValueError('esp_points_bohr must be a finite M×3 array')
    if not isinstance(data.get('density_fit',False),bool):
        raise ValueError('density_fit must be a boolean')
    if data.get('auxbasis') is not None and not data.get('density_fit',False):
        raise ValueError('auxbasis requires density_fit=true')
    checkpoint=data.get('initial_checkpoint')
    if checkpoint is not None:
        keys={'result_path','result_sha256','checkpoint_path','checkpoint_sha256'}
        if not isinstance(checkpoint,dict) or set(checkpoint)!=keys:
            raise ValueError('initial_checkpoint requires exactly result_path/result_sha256 and checkpoint_path/checkpoint_sha256')
        for key in ('result_path','checkpoint_path'):
            if not isinstance(checkpoint[key],str) or not Path(checkpoint[key]).is_absolute():
                raise ValueError('initial_checkpoint file paths must be absolute')
        if any(not _sha256_string(checkpoint[key]) for key in ('result_sha256','checkpoint_sha256')):
            raise ValueError('initial_checkpoint requires explicit lowercase SHA-256 hashes')
    result = dict(data)
    result.update(charge=charge,spin=spin,threads=_integer(data.get('threads',2),'threads',1,2),
        max_memory_mb=_integer(data.get('max_memory_mb',8000),'max_memory_mb',128,8000),
        max_scf_cycles=_integer(data.get('max_scf_cycles',100),'max_scf_cycles',1,500),
        max_cphf_cycles=_integer(data.get('max_cphf_cycles',100),'max_cphf_cycles',1,500),
        grid_level=_integer(data.get('grid_level',3),'grid_level',0,9),
        esp_batch_size=_integer(data.get('esp_batch_size',256),'esp_batch_size',1,4096),
        esp_crosscheck_points=_integer(data.get('esp_crosscheck_points',8),'esp_crosscheck_points',1,64),
        max_wall_seconds=_integer(data.get('max_wall_seconds',7200),'max_wall_seconds',1,172800),
        scf_conv_tol=_positive(data.get('scf_conv_tol',1e-10),'scf_conv_tol',1e-6),
        scf_conv_tol_grad=_positive(data.get('scf_conv_tol_grad',1e-7),'scf_conv_tol_grad',1e-4))
    opt=data.get('optimization',{})
    if not isinstance(opt,dict) or set(opt)-({'enabled','max_steps','freeze_atom_ids','dihedral_constraints','dihedral_tolerance_degrees'} | set(OPT_DEFAULTS)):
        raise ValueError('Unknown or malformed optimization options')
    if not isinstance(opt.get('enabled',False),bool):
        raise ValueError('optimization.enabled must be boolean')
    if opt and not opt.get('enabled',False) and set(opt)-{'enabled'}:
        raise ValueError('Optimization options supplied while optimization is disabled')
    freezes=opt.get('freeze_atom_ids',[])
    if not isinstance(freezes,list) or any(x not in ids for x in freezes) or len(set(freezes))!=len(freezes):
        raise ValueError('freeze_atom_ids must uniquely identify existing atoms')
    if opt.get('enabled',False) and len(freezes)==len(ids):
        raise ValueError('Optimization cannot freeze all atoms; request single-point operations instead')
    dihedrals=opt.get('dihedral_constraints',[])
    if not isinstance(dihedrals,list) or len(dihedrals)>32:
        raise ValueError('dihedral_constraints must be a list of at most 32 explicit constraints')
    if 'dihedral_tolerance_degrees' in opt and not dihedrals:
        raise ValueError('dihedral_tolerance_degrees requires at least one dihedral constraint')
    seen=set()
    for constraint in dihedrals:
        if not isinstance(constraint,dict) or set(constraint)!={'atom_ids','target_degrees'}:
            raise ValueError('Each dihedral constraint requires exactly atom_ids and target_degrees')
        selected=constraint['atom_ids']
        if (not isinstance(selected,list) or len(selected)!=4 or any(not isinstance(x,str) or x not in ids for x in selected) or
            len(set(selected))!=4):
            raise ValueError('A dihedral constraint needs four distinct existing stable atom IDs')
        target=constraint['target_degrees']
        if isinstance(target,bool) or not isinstance(target,(int,float)) or not math.isfinite(target) or abs(target)>360:
            raise ValueError('Dihedral target_degrees must be finite and between -360 and 360')
        key=min(tuple(selected),tuple(selected[::-1]))
        if key in seen:raise ValueError('Duplicate or conflicting dihedral constraints, including reversed atom order')
        seen.add(key)
        if set(selected)<=set(freezes):
            raise ValueError('A dihedral constraint whose four atoms are Cartesian-frozen is conflicting or redundant')
        dihedral_degrees([coords[ids.index(x)] for x in selected])
    result['optimization']={'enabled':opt.get('enabled',False),'max_steps':_integer(opt.get('max_steps',100),'optimization.max_steps',1,1000),
        'freeze_atom_ids':freezes,'dihedral_constraints':dihedrals,
        'dihedral_tolerance_degrees':_positive(opt.get('dihedral_tolerance_degrees',0.1),'dihedral_tolerance_degrees',0.5),
        **{k:_positive(opt.get(k,v),k,0.1) for k,v in OPT_DEFAULTS.items()}}
    return result


def flatten_hessian(native, natoms):
    import numpy as np
    native=np.asarray(native,dtype=float)
    if native.shape!=(natoms,natoms,3,3):
        raise ValueError(f'Unexpected native Hessian shape: {native.shape}')
    return native.transpose(0,2,1,3).reshape(3*natoms,3*natoms)


def electrostatic_potential(mol, density, points, batch_size=256, crosscheck_points=8):
    """031-MEP public API, with independent int1e_rinv spot checks.

    ECP convention: point effective ionic-core charge and explicit valence
    electron density. No claim to reconstruct an all-electron core potential.
    """
    import numpy as np
    from pyscf import df, gto
    points=np.asarray(points,dtype=float)
    coords=mol.atom_coords(unit='Bohr')
    charges=mol.atom_charges()
    if points.ndim!=2 or points.shape[1]!=3 or not len(points) or not np.isfinite(points).all():
        raise ValueError('Invalid ESP point array')
    # Bound each three-center integral batch to approximately 64 MB.
    batch_size=min(batch_size,max(1,64_000_000//(8*mol.nao_nr()**2)))
    potential=np.empty(len(points))
    for start in range(0,len(points),batch_size):
        grid=points[start:start+batch_size]
        distances=np.linalg.norm(grid[:,None,:]-coords[None,:,:],axis=2)
        if np.any(distances<1e-8):
            raise ValueError('An ESP point coincides with a nucleus; potential is singular')
        nuclear=(charges[None,:]/distances).sum(axis=1)
        fake=gto.fakemol_for_charges(grid)
        integrals=df.incore.aux_e2(mol,fake,intor='int3c2e',aosym='s1')
        electronic=np.einsum('ijp,ij->p',integrals,density)
        potential[start:start+len(grid)]=nuclear-electronic
    chosen=np.unique(np.linspace(0,len(points)-1,min(crosscheck_points,len(points)),dtype=int))
    errors=[]
    for index in chosen:
        point=points[index]
        nuclear=np.sum(charges/np.linalg.norm(coords-point,axis=1))
        with mol.with_rinv_origin(point):
            reference=nuclear-np.einsum('ij,ij->',mol.intor('int1e_rinv'),density)
        errors.append(abs(reference-potential[index]))
    largest=max(errors,default=0.0)
    if not np.isfinite(potential).all() or not math.isfinite(largest) or largest>1e-8:
        raise ValueError(f'ESP independent integral crosscheck failed ({largest} Hartree/e)')
    return potential,{'method':'Batched int3c2e/fakemol_for_charges; independent int1e_rinv',
        'batch_size':batch_size,'checked_point_indices':chosen.tolist(),
        'max_absolute_error_hartree_per_e':largest,'tolerance_hartree_per_e':1e-8,
        'passed':True,'ecp_convention':'Effective point ionic cores plus valence electron density; all-electron core ESP is not reconstructed.'}


def initial_density(data,mol,out,report):
    """Copy and verify an accepted same-state, same-geometry SCF checkpoint.

    A density-fit setting may differ: the density is only an initial guess and
    the requested new method must converge independently. No AO projection,
    atom reordering, charge change or implicit geometry transfer is performed.
    """
    import numpy as np
    import pyscf
    from pyscf import dft,scf
    request=data.get('initial_checkpoint')
    if request is None:return None
    destination=out/'initial-checkpoint';destination.mkdir(exist_ok=False)
    def copy_checked(source,name,expected):
        source=Path(source)
        if source.is_symlink() or not source.is_file() or source.stat().st_size>128*1024*1024:
            raise ValueError('Checkpoint inputs must be regular, non-symlink files below 128 MB')
        raw=source.read_bytes()
        if not _sha256_string(expected) or hashlib.sha256(raw).hexdigest()!=expected:
            raise ValueError(f'Initial checkpoint artifact hash mismatch: {name}')
        (destination/name).write_bytes(raw)
        return raw
    prior=json.loads(copy_checked(request['result_path'],'result.json',request['result_sha256']))
    if not (prior.get('schema_version')==1 and prior.get('accepted') is True and
            prior.get('status')=='completed' and prior.get('scf',{}).get('converged') is True):
        raise ValueError('Initial checkpoint must belong to a completed accepted worker SCF result')
    source_folder=Path(request['result_path']).parent
    if prior.get('arrays_file')!='arrays.npz' or prior.get('resolved_method_file')!='resolved-method.json':
        raise ValueError('Initial checkpoint result has an unsupported artifact layout')
    copy_checked(source_folder/'arrays.npz','arrays.npz',prior.get('arrays_sha256'))
    resolved=json.loads(copy_checked(source_folder/'resolved-method.json','resolved-method.json',prior.get('resolved_method_file_sha256')))
    copy_checked(request['checkpoint_path'],'scf.chk',request['checkpoint_sha256'])
    if prior.get('checkpoint_sha256') is not None and prior['checkpoint_sha256']!=request['checkpoint_sha256']:
        raise ValueError('Initial checkpoint hash differs from its accepted source result')
    elements=[mol.atom_symbol(i) for i in range(mol.natm)]
    same_state={'atom_ids':data['atom_ids'],'elements':elements,'charge':data['charge'],
        'spin_2S':data['spin'],'explicit_electron_count':mol.nelectron,'nao':mol.nao_nr()}
    if any(prior.get(key)!=value for key,value in same_state.items()):
        raise ValueError('Initial checkpoint atom order, charge, spin, electron count or AO count differs')
    if prior.get('versions',{}).get('pyscf')!=pyscf.__version__:
        raise ValueError('Initial checkpoint requires the same PySCF version and AO conventions')
    method_contract={'method':data['method'],'basis_sha256':_hash(mol._basis),
        'ecp_sha256':_hash(mol._ecp),'spherical_basis':True,'symmetry':False,
        'grid_level':data['grid_level'] if data['method']=='RKS' else None}
    if data['method']=='RKS':
        method_contract['functional_resolved_sha256']=_hash(dft.libxc.parse_xc(data['functional']))
    for key,value in method_contract.items():
        if resolved.get(key)!=value or prior.get('method',{}).get(key)!=value:
            raise ValueError(f'Initial checkpoint resolved method differs: {key}')
    if _hash(resolved.get('basis_resolved'))!=method_contract['basis_sha256'] or _hash(resolved.get('ecp_resolved'))!=method_contract['ecp_sha256']:
        raise ValueError('Initial checkpoint resolved basis/ECP payload differs from its hash')
    if data['method']=='RKS' and _hash(resolved.get('functional_resolved'))!=method_contract['functional_resolved_sha256']:
        raise ValueError('Initial checkpoint functional payload differs from its hash')
    checkpoint_mol,checkpoint_scf=scf.chkfile.load_scf(str(destination/'scf.chk'))
    if (checkpoint_mol.charge!=mol.charge or checkpoint_mol.spin!=mol.spin or
        checkpoint_mol.nelectron!=mol.nelectron or checkpoint_mol.nao_nr()!=mol.nao_nr() or
        bool(checkpoint_mol.cart)!=bool(mol.cart) or bool(checkpoint_mol.symmetry)!=bool(mol.symmetry) or
        [checkpoint_mol.atom_symbol(i) for i in range(checkpoint_mol.natm)]!=elements or
        _hash(checkpoint_mol._basis)!=_hash(mol._basis) or _hash(checkpoint_mol._ecp)!=_hash(mol._ecp)):
        raise ValueError('Checkpoint molecule state, AO convention, elements or resolved basis/ECP differs')
    with np.load(destination/'arrays.npz',allow_pickle=False) as arrays:
        source_coords=np.asarray(arrays['coords_bohr'])
        if not np.array_equal(arrays['effective_nuclear_charges'],mol.atom_charges()):
            raise ValueError('Initial checkpoint effective nuclear charges differ')
    coordinate_errors=[]
    for coords in (source_coords,checkpoint_mol.atom_coords(unit='Bohr')):
        if coords.shape!=(mol.natm,3) or not np.isfinite(coords).all():
            raise ValueError('Initial checkpoint geometry is malformed or nonfinite')
        coordinate_errors.append(float(np.max(np.abs(coords-mol.atom_coords(unit='Bohr')))))
    if max(coordinate_errors)>1e-10:
        raise ValueError('Initial checkpoint geometry differs; no density projection or coordinate transfer is permitted')
    energy=float(checkpoint_scf['e_tot']);source_energy=prior.get('energy_hartree')
    if (not isinstance(source_energy,(int,float)) or not math.isfinite(source_energy) or
        not math.isfinite(energy) or abs(energy-source_energy)>1e-9):
        raise ValueError('Initial checkpoint energy does not match its accepted source result')
    coefficients=np.asarray(checkpoint_scf['mo_coeff']);occupations=np.asarray(checkpoint_scf['mo_occ'])
    orbital_energies=np.asarray(checkpoint_scf['mo_energy'])
    if (coefficients.ndim!=2 or coefficients.shape[0]!=mol.nao_nr() or
        not 1<=coefficients.shape[1]<=mol.nao_nr() or occupations.shape!=(coefficients.shape[1],) or
        orbital_energies.shape!=occupations.shape or
        any(not np.isrealobj(x) or not np.isfinite(x).all() for x in (coefficients,occupations,orbital_energies)) or
        np.max(np.minimum(np.abs(occupations),np.abs(occupations-2)))>1e-8 or
        abs(float(occupations.sum())-mol.nelectron)>1e-8):
        raise ValueError('Initial checkpoint has invalid closed-shell orbitals/occupancies')
    overlap=mol.intor_symmetric('int1e_ovlp')
    orthogonality=float(np.max(np.abs(coefficients.T@overlap@coefficients-np.eye(len(occupations)))))
    if orthogonality>1e-7:raise ValueError('Initial checkpoint orbitals are not orthonormal in the requested AO basis')
    density=np.asarray(scf.hf.make_rdm1(coefficients,occupations))
    electrons=float(np.einsum('ij,ji->',density,overlap))
    if (density.shape!=(mol.nao_nr(),mol.nao_nr()) or not np.isfinite(density).all() or
        np.max(np.abs(density-density.T))>1e-10 or abs(electrons-mol.nelectron)>1e-6):
        raise ValueError('Initial checkpoint density symmetry, finiteness or electron count failed')
    report['initial_checkpoint']={'accepted_as_initial_guess':True,
        'source_result_sha256':request['result_sha256'],'source_checkpoint_sha256':request['checkpoint_sha256'],
        'source_arrays_sha256':prior['arrays_sha256'],'source_resolved_method_sha256':prior['resolved_method_file_sha256'],
        'source_worker_sha256':prior.get('worker_sha256'),'source_energy_hartree':energy,
        'energy_provenance':'Checkpoint e_tot matches the pinned accepted result; source energy is not recomputed here.',
        'source_density_fit':resolved.get('density_fit',False),'requested_density_fit':data.get('density_fit',False),
        'maximum_coordinate_difference_bohr':max(coordinate_errors),'coordinate_tolerance_bohr':1e-10,
        'orbital_orthogonality_max_error':orthogonality,'density_electron_count':electrons,
        'copied_folder':'initial-checkpoint','source_files_modified':False,
        'meaning':'Same-state same-geometry initial density only; requested SCF must independently converge.'}
    return density


def _create_evaluation_journal(data, elements, initial_bohr, mean_field, out, report,
                               bohr_to_angstrom):
    """Candidate-only observer: no SCF calls, coordinate edits or acceptance.

    Positional stable IDs come from the validated input. PySCF exposes element
    order, not independent stable IDs; same-element permutations cannot be
    independently detected here. Current engine/scanner/base coordinates and
    returned derivatives must agree exactly before an evaluation is committed.
    A completed evaluation can still be a rejected optimizer trial geometry.
    """
    import importlib.util
    helper_path=Path(__file__).with_name('qm_evaluation_journal.py')
    helper_sha=hashlib.sha256(helper_path.read_bytes()).hexdigest()
    if helper_sha!='b3a0543b5b2394373e7a787ad01ecd77865cc507f85bec8ec99efe6efed26528':
        raise ValueError('Candidate evaluation-journal helper hash mismatch')
    spec=importlib.util.spec_from_file_location('_candidate_qm_evaluation_journal',helper_path)
    journal_module=importlib.util.module_from_spec(spec);spec.loader.exec_module(journal_module)
    elements=tuple(elements)
    request=json.loads(json.dumps(_plain(data),allow_nan=False))
    original_resolved_request_sha=_hash(request)
    request['elements']=list(elements)  # Same canonical symbols selected by gto.M.
    binding=journal_module.make_binding(request,source_input_sha256=report['input_sha256'],
        worker_source_sha256=report['worker_sha256'])
    expected_identity=binding['identity']
    native_family=type(mean_field)

    def native_method(method):
        if not isinstance(method,native_family):
            raise ValueError('Current callback SCF method family differs from its initial method')
        molecule=method.mol
        density_fit=getattr(method,'with_df',None)
        return {'native_family':native_family.__module__+'.'+native_family.__qualname__,
            'basis_sha256':_hash(molecule._basis),'ecp_sha256':_hash(molecule._ecp),
            'cartesian_basis':bool(molecule.cart),'symmetry':bool(molecule.symmetry),
            'xc':getattr(method,'xc',None),
            'grid_level':getattr(method.grids,'level',None) if hasattr(method,'grids') else None,
            'density_fit':density_fit is not None,
            'auxbasis':_plain(getattr(density_fit,'auxbasis',None)),
            'max_scf_cycles':method.max_cycle,'scf_conv_tol':method.conv_tol,
            'scf_conv_tol_grad':method.conv_tol_grad,'max_memory_mb':method.max_memory}

    initial_method=native_method(mean_field)
    initial_method_sha=_hash(initial_method)
    conversion=journal_module._number(bohr_to_angstrom)
    if conversion<=0:raise ValueError('Invalid native Bohr to Angstrom conversion')
    binding['callback_contract']={
        'helper_source_sha256':helper_sha,
        'original_resolved_input_sha256':original_resolved_request_sha,
        'original_input_elements':list(data['elements']),
        'canonical_native_elements':list(elements),
        'input_coordinate_unit':'bohr' if 'coords_bohr' in data else 'angstrom',
        'native_initial_geometry_bohr_sha256':_hash(_plain(initial_bohr)),
        'native_initial_method':initial_method,'native_initial_method_sha256':initial_method_sha,
        'bohr_to_angstrom':conversion,'conversion_source':'pyscf.lib.param.BOHR',
        'callback_coordinate_unit':'bohr','callback_gradient_unit':'hartree/bohr',
        'gradient_convention':'positive dE/dx; no sign or unit conversion',
        'stable_id_origin':'validated positional input map; native element order checked',
        'geometry_meaning':'completed SCF/gradient evaluation, possibly an optimizer trial; not an accepted minimum'}
    journal=journal_module.EvaluationJournal(out/'evaluation-journal',binding,create=True)
    report['optimization']['evaluation_journal']={
        'directory':'evaluation-journal','manifest_sha256':journal.manifest_sha256,
        'helper_source_sha256':helper_sha,'accepted':False,'optimization_converged':False,
        'checkpoint_reuse_authorized':False,'meaning':binding['callback_contract']['geometry_meaning']}

    def record_evaluation(env):
        scanner=env['g_scanner']
        if scanner.converged is False:
            # Native PySCF asserts this immediately after the callback. Preserve
            # that behavior; there is no completed-SCF journal record to write.
            return None
        if scanner.converged is not True or scanner.base.converged is not True:
            raise ValueError('Callback requires an explicit native SCF convergence flag')
        if scanner.unit!='au':raise ValueError('Callback native gradient unit is not atomic units')
        atmlst=getattr(scanner,'atmlst',None)
        if atmlst is not None:
            indices=_plain(atmlst)
            if (not isinstance(indices,list) or any(type(i) is not int for i in indices)
                    or indices!=list(range(len(elements)))):
                raise ValueError('Callback gradient atom indices are not the complete input order')
        xyz=journal_module._matrix(env['coords'],len(elements),'callback coords',flat=True)
        gradient=journal_module._matrix(env['gradients'],len(elements),'callback gradient',flat=False)
        energy=journal_module._number(env['energy'])
        if energy!=journal_module._number(scanner.base.e_tot):
            raise ValueError('Callback energy differs from current SCF energy')
        if gradient!=journal_module._matrix(scanner.de,len(elements),'scanner gradient',flat=False):
            raise ValueError('Callback gradient differs from current scanner gradient')
        for molecule in (env['mol'],scanner.mol,scanner.base.mol):
            observed={'atom_ids':list(expected_identity['atom_ids']),
                'elements':[molecule.atom_symbol(i) for i in range(molecule.natm)],
                'charge':molecule.charge,'spin':molecule.spin}
            if journal_module._identity(observed)!=expected_identity:
                raise ValueError('Callback native element order or molecular state changed')
            if xyz!=journal_module._matrix(molecule.atom_coords(unit='Bohr'),len(elements),'native coords',flat=False):
                raise ValueError('Callback geometry differs from current native molecule')
        if _hash(native_method(scanner.base))!=initial_method_sha:
            raise ValueError('Callback current method or resolved basis/settings changed')
        result=journal.append(cycle=env['self'].cycle,coords=xyz,gradient=gradient,energy=energy,
            units=journal_module.UNITS,observed_identity=observed,scf_converged=True,
            requested_method_fingerprint_sha256=binding['requested_method_fingerprint_sha256'])
        return {'file':'evaluation-journal/'+Path(result['path']).name,
            'file_sha256':result['file_sha256'],'accepted':False}
    return record_evaluation


def calculate(data, out, report):
    import numpy as np
    import pyscf
    from pyscf import dft, gto, lib, scf
    lib.num_threads(data['threads'])
    actual_threads=int(lib.num_threads())
    report['resources'].update(requested_threads=data['threads'],actual_pyscf_threads=actual_threads)
    report['runtime_manifest']=capture_runtime_manifest(out)
    original=np.asarray(data.get('coords_bohr',data.get('coords_angstrom')),dtype=float)
    unit='Bohr' if 'coords_bohr' in data else 'Angstrom'
    elements=[gto.mole._std_symbol(x) for x in data['elements']]
    if any(gto.charge(x)<=0 for x in elements):
        raise ValueError('Ghost atoms and zero nuclear charges are not supported')
    log=out/'native.log'
    mol=gto.M(atom=list(zip(elements,original.tolist())),unit=unit,
        charge=data['charge'],spin=data['spin'],basis=data['basis'],ecp=data.get('ecp') or {},
        symmetry=False,cart=False,max_memory=data['max_memory_mb'],verbose=4,output=str(log))
    if mol.nelectron<=0 or mol.nelectron%2:
        raise ValueError('Closed-shell calculation requires a positive even number of explicit electrons')
    if len(mol._basis)!=len(set(elements)):
        raise ValueError('Resolved basis does not cover all requested elements')
    requested_ecp=data.get('ecp')
    if requested_ecp:
        resolved_elements={x.lower() for x in mol._ecp}
        if not resolved_elements or (isinstance(requested_ecp,dict) and
                set(x.lower() for x in requested_ecp)-resolved_elements):
            raise ValueError('A requested ECP was not resolved; silently reverting to an all-electron model is not permitted')
    initial_bohr=mol.atom_coords(unit='Bohr').copy()
    if len(elements)>1:
        distances=np.linalg.norm(initial_bohr[:,None,:]-initial_bohr[None,:,:],axis=2)
        np.fill_diagonal(distances,np.inf)
        if distances.min()<1e-6:
            raise ValueError('Coincident nuclei are not permitted')
    def progress(phase,**detail):
        report['phase']=phase
        _json(out/'progress.json',{'phase':phase,'accepted':False,**detail})
    def make_method(geometry):
        mean_field=scf.RHF(geometry) if data['method']=='RHF' else dft.RKS(geometry)
        if data['method']=='RKS':
            mean_field.xc=data['functional']
            mean_field.grids.level=data['grid_level']
        if data.get('density_fit',False):
            mean_field=mean_field.density_fit(auxbasis=data.get('auxbasis'))
        mean_field.max_cycle=data['max_scf_cycles']
        mean_field.conv_tol=data['scf_conv_tol']
        mean_field.conv_tol_grad=data['scf_conv_tol_grad']
        mean_field.max_memory=data['max_memory_mb']
        mean_field.chkfile=str(out/'scf.chk')
        return mean_field
    mean_field=make_method(mol)
    if data.get('initial_checkpoint'):progress('initial_checkpoint')
    guess=initial_density(data,mol,out,report)
    opt=data['optimization']
    report['optimization']={'requested':opt['enabled'],'converged':None,'settings':opt,'steps':[]}
    if opt['enabled']:
        from pyscf.geomopt import geometric_solver
        if guess is not None:
            progress('initial_checkpoint_scf')
            warm_energy=float(mean_field.kernel(dm0=guess))
            report['initial_checkpoint']['requested_method_initial_scf_converged']=bool(mean_field.converged)
            if not mean_field.converged or not math.isfinite(warm_energy):
                raise RuntimeError('Requested-method initial checkpoint SCF did not converge')
            # geomeTRIC's SCF scanner inherits these converged orbitals; its
            # subsequent geometry changes are explicit optimization steps.
            guess=None
        frozen_indices=[data['atom_ids'].index(x) for x in opt['freeze_atom_ids']]
        constraint_path=None;constraint_lines=[]
        if frozen_indices:
            constraint_lines.extend(['$freeze','xyz '+','.join(str(x+1) for x in frozen_indices)])
        dihedrals=opt['dihedral_constraints']
        initial_dihedrals=_dihedral_rows(initial_bohr,data['atom_ids'],dihedrals,opt['dihedral_tolerance_degrees'])
        report['optimization']['initial_dihedral_constraints']=initial_dihedrals
        if dihedrals:
            constraint_lines.append('$set')
            constraint_lines.extend('dihedral '+ ' '.join(str(x+1) for x in row['atom_indices_zero_based'])+
                ' '+format(row['canonical_target_degrees'],'.15g') for row in initial_dihedrals)
        if constraint_lines:
            constraint_path=out/'constraints.txt'
            constraint_path.write_text('\n'.join(constraint_lines)+'\n')
        report['optimization']['constraint_atom_indices_zero_based']=frozen_indices
        report['optimization']['constraint_file_sha256']=hashlib.sha256(constraint_path.read_bytes()).hexdigest() if constraint_path else None
        record_evaluation=_create_evaluation_journal(data,elements,initial_bohr,mean_field,out,report,lib.param.BOHR)
        def callback(env):
            evaluation_record=record_evaluation(env)
            energy=float(env['energy']);gradient=np.asarray(env['gradients'])
            if not math.isfinite(energy) or not np.isfinite(gradient).all():
                raise ValueError('Nonfinite energy/gradient during optimization')
            step={'cycle':int(env['self'].cycle),'energy_hartree':energy,
                  'gradient_norm_hartree_per_bohr':float(np.linalg.norm(gradient)),
                  'scf_converged':bool(env['g_scanner'].converged)}
            if evaluation_record is not None:step['evaluation_record']=evaluation_record
            if dihedrals:
                step['dihedrals']=_dihedral_rows(env['coords'],data['atom_ids'],dihedrals,opt['dihedral_tolerance_degrees'])
            report['optimization']['steps'].append(step)
            progress('optimization',**step)
        progress('optimization')
        converged,mol=geometric_solver.kernel(mean_field,assert_convergence=True,
            constraints=str(constraint_path) if constraint_path else None,callback=callback,
            maxsteps=opt['max_steps'],enforce=0.1,**{k:opt[k] for k in OPT_DEFAULTS})
        report['optimization']['converged']=bool(converged)
        final_dihedrals=_dihedral_rows(mol.atom_coords(unit='Bohr'),data['atom_ids'],dihedrals,opt['dihedral_tolerance_degrees'])
        report['optimization']['final_dihedral_constraints']=final_dihedrals
        report['optimization']['dihedral_constraints_satisfied']=all(row['target_satisfied'] for row in final_dihedrals)
        displacement=0.0
        if frozen_indices:
            displacement=float(np.linalg.norm(mol.atom_coords()-initial_bohr,axis=1)[frozen_indices].max())
            report['optimization']['max_frozen_displacement_bohr']=displacement
        _json(out/'optimization-report.json',report['optimization'])
        if not converged:
            raise RuntimeError('Geometry optimization did not converge within the step limit')
        if displacement>1e-5:
            raise RuntimeError(f'Frozen atom coordinate check failed ({displacement} bohr)')
        if not report['optimization']['dihedral_constraints_satisfied']:
            raise RuntimeError('Final actual dihedral does not satisfy the requested periodic target tolerance')
        # Recompute all requested properties consistently at the actual final
        # geometry; never report the input geometry next to optimized derivatives.
        mean_field=make_method(mol)
    if [mol.atom_symbol(i) for i in range(mol.natm)]!=elements:
        raise RuntimeError('Native optimization changed atom order/identity')
    progress('scf')
    energy=float(mean_field.kernel(dm0=guess))
    report['scf']={'converged':bool(mean_field.converged),'cycles':int(getattr(mean_field,'cycles',0)),
        'max_cycles':data['max_scf_cycles'],'conv_tol_hartree':data['scf_conv_tol'],
        'conv_tol_grad':data['scf_conv_tol_grad']}
    if not mean_field.converged or not math.isfinite(energy):
        raise RuntimeError('SCF did not converge to a finite energy within the cycle limit')
    density=np.asarray(mean_field.make_rdm1())
    if density.ndim!=2 or not np.isfinite(density).all():
        raise ValueError('Invalid restricted density matrix')
    density_electrons=float(np.einsum('ij,ji->',density,mean_field.get_ovlp()))
    if abs(density_electrons-mol.nelectron)>1e-6:
        raise ValueError('Density matrix electron count does not match the requested charge')
    arrays={'coords_bohr':mol.atom_coords(unit='Bohr'),'coords_angstrom':mol.atom_coords(unit='Angstrom'),
        'atomic_numbers':np.asarray([gto.charge(x) for x in elements],dtype=np.int64),
        'effective_nuclear_charges':np.asarray(mol.atom_charges(),dtype=float)}
    if 'gradient' in data['operations']:
        progress('gradient')
        gradient=mean_field.nuc_grad_method().kernel()
        if np.shape(gradient)!=(mol.natm,3):raise ValueError('Unexpected gradient shape')
        arrays['gradient_hartree_per_bohr']=np.asarray(gradient)
    if 'hessian' in data['operations']:
        progress('hessian')
        hessian_method=mean_field.Hessian()
        hessian_method.max_cycle=data['max_cphf_cycles']
        native=hessian_method.kernel()
        hessian=flatten_hessian(native,mol.natm)
        asymmetry=float(np.max(np.abs(hessian-hessian.T)))
        if asymmetry>1e-7:raise ValueError(f'Hessian symmetry check failed ({asymmetry})')
        arrays['hessian_hartree_per_bohr2']=hessian
        report['hessian']={'completed':True,'response_solver_completed':True,
            'max_cphf_cycles':data['max_cphf_cycles'],
            'native_order':'atom_i,atom_j,xyz_i,xyz_j','output_order':'atom_i/xyz_i,atom_j/xyz_j',
            'max_asymmetry_hartree_per_bohr2':asymmetry}
    if 'esp' in data['operations']:
        progress('esp',point_count=len(data['esp_points_bohr']))
        points=np.asarray(data['esp_points_bohr'],dtype=float)
        potential,check=electrostatic_potential(mol,density,points,data['esp_batch_size'],data['esp_crosscheck_points'])
        arrays.update(esp_points_bohr=points,esp_hartree_per_e=potential)
        report['esp']=check
    for key,value in arrays.items():
        if not np.isfinite(value).all():raise ValueError(f'Nonfinite output array: {key}')
    resolved={'method':data['method'],'basis_requested':data['basis'],'basis_resolved':_plain(mol._basis),
        'ecp_requested':data.get('ecp'),'ecp_resolved':_plain(mol._ecp),
        'basis_sha256':_hash(mol._basis),'ecp_sha256':_hash(mol._ecp),
        'spherical_basis':True,'symmetry':False,'density_fit':data.get('density_fit',False),
        'functional_requested':data.get('functional'),'grid_level':data['grid_level'] if data['method']=='RKS' else None,
        'requested_threads':data['threads'],'actual_pyscf_threads':actual_threads}
    if data['method']=='RKS':
        resolved['functional_resolved']=_plain(dft.libxc.parse_xc(data['functional']))
        resolved['functional_resolved_sha256']=_hash(resolved['functional_resolved'])
    if getattr(mean_field,'with_df',None):
        resolved['auxbasis_resolved']=_plain(mean_field.with_df.auxmol._basis)
        resolved['auxbasis_sha256']=_hash(resolved['auxbasis_resolved'])
    _json(out/'resolved-method.json',resolved)
    np.savez_compressed(out/'arrays.npz',**arrays)
    report.update(energy_hartree=energy,atom_ids=data['atom_ids'],elements=elements,
        requested_threads=data['threads'],actual_pyscf_threads=actual_threads,
        charge=data['charge'],spin_2S=data['spin'],explicit_electron_count=mol.nelectron,
        density_electron_count=density_electrons,nao=mol.nao_nr(),
        method={k:v for k,v in resolved.items() if not k.endswith('_resolved')},
        resolved_method_file='resolved-method.json',resolved_method_file_sha256=hashlib.sha256((out/'resolved-method.json').read_bytes()).hexdigest(),
        arrays_file='arrays.npz',arrays_sha256=hashlib.sha256((out/'arrays.npz').read_bytes()).hexdigest(),
        checkpoint_file='scf.chk',checkpoint_sha256=hashlib.sha256((out/'scf.chk').read_bytes()).hexdigest(),
        arrays={k:{'shape':list(v.shape),'dtype':str(v.dtype)} for k,v in arrays.items()},
        units={'energy':'hartree','coords_bohr':'bohr','coords_angstrom':'angstrom',
            'gradient_hartree_per_bohr':'hartree/bohr','hessian_hartree_per_bohr2':'hartree/bohr^2',
            'esp_points_bohr':'bohr','esp_hartree_per_e':'hartree/e','effective_nuclear_charges':'e',
            'bohr_to_angstrom':float(lib.param.BOHR),'bohr_constant_source':'pyscf.lib.param.BOHR'},
        versions={'python':platform.python_version(),'pyscf':pyscf.__version__,'numpy':np.__version__,
            'geometric':importlib.metadata.version('geometric'),'libxc':dft.libxc.__version__})
    mol.stdout.flush()
    return report


def run(input_path, output_directory):
    out=Path(output_directory).resolve()
    out.mkdir(parents=True,exist_ok=False)
    start=time.monotonic()
    report={'schema_version':SCHEMA_VERSION,'accepted':False,'status':'failed','phase':'input',
        'worker_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        'meaning_of_accepted':'Requested numerical calculation converged and passed stated checks; not force-field or chemical-model validation.'}
    old_handler=None
    try:
        raw=Path(input_path).read_bytes()
        report['input_sha256']=hashlib.sha256(raw).hexdigest()
        (out/'input.json').write_bytes(raw)
        data=validate_input(json.loads(raw))
        _json(out/'resolved-input.json',data)
        report['resources']={'threads':data['threads'],'max_memory_mb':data['max_memory_mb'],
            'memory_limit_kind':'PySCF allocation setting, not an operating-system hard limit',
            'max_wall_seconds':data['max_wall_seconds']}
        for key in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','VECLIB_MAXIMUM_THREADS','NUMEXPR_NUM_THREADS'):
            os.environ[key]=str(data['threads'])
        if hasattr(signal,'SIGALRM'):
            def timeout(signum,frame):
                raise TimeoutError('QM worker exceeded its configured wall-time limit')
            old_handler=signal.signal(signal.SIGALRM,timeout)
            signal.alarm(data['max_wall_seconds'])
        calculate(data,out,report)
        report.update(accepted=True,status='completed',phase='completed')
    except Exception as exc:
        report.update(accepted=False,status='failed',error={'type':type(exc).__name__,'message':str(exc)})
        (out/'error-traceback.txt').write_text(traceback.format_exc())
        # Never leave a result-array artifact next to a failed acceptance record.
        (out/'arrays.npz').unlink(missing_ok=True)
    finally:
        if old_handler is not None:
            signal.alarm(0);signal.signal(signal.SIGALRM,old_handler)
        report['elapsed_seconds']=round(time.monotonic()-start,6)
        _json(out/'result.json',report)
        _json(out/'progress.json',{'phase':report['phase'],'status':report['status'],'accepted':report['accepted']})
    return report


if __name__=='__main__':
    if len(sys.argv)!=3:
        raise SystemExit('Usage: qm_worker.py INPUT.json FRESH_OUTPUT_DIRECTORY')
    try:
        outcome=run(sys.argv[1],sys.argv[2])
    except Exception as exc:
        print(f'{type(exc).__name__}: {exc}',file=sys.stderr)
        raise SystemExit(2)
    print(json.dumps({'status':outcome['status'],'accepted':outcome['accepted'],'elapsed_seconds':outcome['elapsed_seconds']}))
    raise SystemExit(0 if outcome['accepted'] else 1)
