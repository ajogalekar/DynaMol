"""Conventional RHF reference from a validated same-geometry DF initial guess.

The converged calculation and gradient use the conventional RHF Hamiltonian.
DF supplies only the initial density; it never supplies reference properties.
"""
from __future__ import annotations
import argparse
import hashlib
import json
import signal
import time
from pathlib import Path

import numpy as np
from pyscf import gto,lib,scf
import pyscf


def digest(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def run(source,output):
    source=source.resolve();output=output.resolve();output.mkdir(parents=True,exist_ok=False)
    started=time.monotonic();report={'accepted':False,'method':'Conventional RHF/6-31G* with converged same-geometry DF density as initial guess only',
        'source':str(source),'source_result_sha256':digest(source/'result.json'),'initial_guess_checkpoint_sha256':digest(source/'scf.chk'),
        'comparison_script_sha256':digest(__file__),'pyscf_version':pyscf.__version__,'threads':2,'max_memory_mb':4000,'max_wall_seconds':7200}
    def progress(stage,**kwargs):
        (output/'progress.json').write_text(json.dumps({'accepted':False,'stage':stage,'elapsed_seconds':time.monotonic()-started,**kwargs},indent=2)+'\n')
    def timeout(_signum,_frame):raise TimeoutError('Conventional reference exceeded its two-hour wall limit')
    signal.signal(signal.SIGALRM,timeout);signal.alarm(7200)
    try:
        original=json.loads((source/'result.json').read_text());request=json.loads((source/'input.json').read_text())
        if not original.get('accepted') or request['method']!='RHF' or request.get('ecp') or not request.get('density_fit'):
            raise ValueError('Expected a successful all-electron DF-RHF source job.')
        if request.get('optimization',{}).get('enabled'):raise ValueError('Reference comparison requires the fixed-geometry DF specimen.')
        arrays=np.load(source/'arrays.npz');xyz=arrays['coords_angstrom']
        if digest(source/'arrays.npz')!=original['arrays_sha256']:raise ValueError('Accepted source array hash differs.')
        if not np.allclose(xyz,request['coords_angstrom'],rtol=0,atol=1e-10):raise ValueError('Accepted DF geometry differs from requested fixed geometry.')
        lib.num_threads(2)
        mol=gto.M(atom=list(zip(request['elements'],xyz)),unit='Angstrom',basis=request['basis'],charge=request['charge'],spin=request['spin'],
                  verbose=4,output=str(output/'native.log'),max_memory=4000)
        checkpoint_mol,checkpoint=scf.chkfile.load_scf(str(source/'scf.chk'))
        if (mol.natm!=checkpoint_mol.natm or mol.nelectron!=checkpoint_mol.nelectron or mol.spin!=checkpoint_mol.spin or mol._basis!=checkpoint_mol._basis
            or [mol.atom_symbol(i) for i in range(mol.natm)]!=[checkpoint_mol.atom_symbol(i) for i in range(mol.natm)]
            or not np.allclose(mol.atom_coords(),checkpoint_mol.atom_coords(),rtol=0,atol=1e-10)):
            raise ValueError('DF checkpoint atom order, state, basis or coordinates differ from the conventional reference.')
        coeff=np.asarray(checkpoint['mo_coeff']);occ=np.asarray(checkpoint['mo_occ'])
        if coeff.shape!=(mol.nao_nr(),mol.nao_nr()) or occ.shape!=(mol.nao_nr(),) or not np.isfinite(coeff).all() or not np.isfinite(occ).all():
            raise ValueError('Initial guess orbital dimensions or values are invalid.')
        density=scf.hf.make_rdm1(coeff,occ);overlap=mol.intor_symmetric('int1e_ovlp')
        electrons=float(np.einsum('ij,ji',density,overlap))
        if not np.isfinite(density).all() or not np.allclose(density,density.T,rtol=0,atol=1e-12) or abs(electrons-mol.nelectron)>1e-6:
            raise ValueError('Initial guess density is nonfinite, nonsymmetric or has the wrong electron count.')
        report['initial_guess_validation']={'same_atom_order_basis_state_coordinates':True,'density_shape':list(density.shape),
            'electron_count':electrons,'density_sha256':hashlib.sha256(np.ascontiguousarray(density).tobytes()).hexdigest(),
            'checkpoint_energy_hartree':float(checkpoint['e_tot'])}
        if abs(float(checkpoint['e_tot'])-original['energy_hartree'])>1e-8:raise ValueError('Checkpoint energy differs from accepted DF result.')
        reference=scf.RHF(mol);reference.chkfile=str(output/'scf.chk')
        reference.conv_tol=request.get('scf_conv_tol',1e-10);reference.conv_tol_grad=request.get('scf_conv_tol_grad',1e-7);reference.max_cycle=request.get('max_scf_cycles',100)
        reference.callback=lambda env:progress('conventional_scf',cycle=int(env['cycle']),energy_hartree=float(env['e_tot']))
        progress('conventional_scf');energy=float(reference.kernel(dm0=density))
        if not reference.converged or not np.isfinite(energy):raise RuntimeError('Conventional reference SCF did not converge.')
        progress('conventional_gradient');gradient=reference.nuc_grad_method().kernel()
        if gradient.shape!=(mol.natm,3) or not np.isfinite(gradient).all():raise RuntimeError('Conventional gradient is invalid.')
        np.savez(output/'arrays.npz',coords_angstrom=xyz,gradient_hartree_per_bohr=gradient)
        delta=arrays['gradient_hartree_per_bohr']-gradient
        report.update(accepted=True,energy_hartree=energy,scf_converged=True,scf_cycles=int(reference.cycles),
            atom_ids=request['atom_ids'],elements=request['elements'],
            df_minus_conventional_energy_hartree=float(original['energy_hartree']-energy),
            maximum_gradient_component_difference_hartree_per_bohr=float(np.max(np.abs(delta))),
            rms_gradient_component_difference_hartree_per_bohr=float(np.sqrt(np.mean(delta**2))),
            arrays_sha256=digest(output/'arrays.npz'),
            meaning_of_accepted='Conventional reference converged with validated initial-guess identity and finite derivatives; approximation error and parameter-model accuracy are reported separately.')
    except Exception as error:
        report.update(error=str(error),error_type=type(error).__name__)
    finally:
        signal.alarm(0);report['elapsed_seconds']=time.monotonic()-started
        (output/'result.json').write_text(json.dumps(report,indent=2)+'\n');progress('completed' if report['accepted'] else 'failed')
    if not report['accepted']:raise RuntimeError(report['error'])
    return report


if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('source',type=Path);parser.add_argument('output',type=Path)
    args=parser.parse_args();print(json.dumps(run(args.source,args.output),indent=2))
