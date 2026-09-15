"""Independent saved-trajectory, geometry and final-coordinate export checks."""
import argparse
import csv
import json
from pathlib import Path

import mdtraj as md
import numpy as np
import openmm as mm
from openmm import app,unit
import parmed

from prepare_adduct import digest


def analyze(folder,inputs):
    result=json.loads((folder/'result.json').read_text())
    diagnostic=bool(result.get('diagnostic_only'))
    if diagnostic and result.get('qualified_for_handoff') is not False:
        raise ValueError('Diagnostic source has inconsistent handoff qualification metadata')
    for name,expected in result['input_sha256'].items():
        if digest(inputs/name)!=expected:raise ValueError('Run input changed: '+name)
    native=parmed.load_file(str(inputs/'solvated.prmtop'))
    traj=md.load_dcd(str(folder/'trajectory.dcd'),top=str(inputs/'solvated.prmtop'))
    if len(traj)!=100 or traj.n_atoms!=len(native.atoms) or not np.isfinite(traj.xyz).all():
        raise ValueError('Saved trajectory inventory/coordinates failed')
    state=mm.XmlSerializer.deserialize((folder/'final-state.xml').read_text())
    last=np.asarray(state.getPositions(asNumpy=True).value_in_unit(unit.nanometer))
    last_box=np.asarray(state.getPeriodicBoxVectors(asNumpy=True).value_in_unit(unit.nanometer))
    coordinate_error=float(np.max(abs(last-traj.xyz[-1])))
    box_error=float(np.max(abs(last_box-traj.unitcell_vectors[-1])))
    if coordinate_error>2e-6 or box_error>2e-6:raise ValueError('Final saved frame differs from engine state')
    cov=next(r for r in native.residues if r.name=='COV')
    bonds=[b for b in native.bonds if b.atom1.residue is cov and b.atom2.residue is cov]
    pairs=[[b.atom1.idx,b.atom2.idx] for b in bonds]
    distances=md.compute_distances(traj,pairs,periodic=True)*10
    ideal=np.array([b.type.req for b in bonds])
    ratios=distances/ideal
    if ratios.min()<.6 or ratios.max()>1.5:raise ValueError('Saved adduct bond distortion')
    angles=[a for a in native.angles if any(x.residue is cov for x in [a.atom1,a.atom2,a.atom3])
            and all(x.atomic_number>1 for x in [a.atom1,a.atom2,a.atom3])]
    angle_values=np.degrees(md.compute_angles(traj,[[a.atom1.idx,a.atom2.idx,a.atom3.idx] for a in angles],periodic=True))
    angle_deviations=angle_values-np.array([a.type.theteq for a in angles])
    attachment=next(b for b in bonds if {b.atom1.atomic_number,b.atom2.atomic_number}=={6,16}
                    and {b.atom1.name,b.atom2.name}!={'CB','SG'})
    attached={attachment.atom1.idx,attachment.atom2.idx}
    torsions={}
    for d in native.dihedrals:
        if d.improper or {d.atom2.idx,d.atom3.idx}!=attached:continue
        indices=(d.atom1.idx,d.atom2.idx,d.atom3.idx,d.atom4.idx)
        if any(native.atoms[i].atomic_number==1 for i in indices):continue
        torsions[min(indices,indices[::-1])]=True
    torsion_rows=[]
    for indices in torsions:
        values=md.compute_dihedrals(traj,[indices],periodic=True)[:,0]
        unwrapped=np.degrees(np.unwrap(values))
        torsion_rows.append({'atoms':[native.atoms[i].name for i in indices],'indices':indices,
            'sampled_unwrapped_span_degrees':float(np.ptp(unwrapped)),
            'circular_mean_degrees':float(np.degrees(np.angle(np.mean(np.exp(1j*values))))),
            'circular_resultant_length':float(abs(np.mean(np.exp(1j*values))))})
    stereo=[]
    for center,neighbors in result['stereocenter_indices']:
        vectors=traj.xyz[:,neighbors,:]-traj.xyz[:,center,None,:]
        volume=np.linalg.det(vectors)
        if not (np.all(volume>0) or np.all(volume<0)):raise ValueError('Saved stereocenter inversion')
        stereo.append({'center':center,'minimum_absolute_volume_nm3':float(np.min(abs(volume)))})
    system=mm.XmlSerializer.deserialize((folder/'system.xml').read_text())
    constraints=[tuple(map(int,system.getConstraintParameters(i)[:2])) for i in range(system.getNumConstraints())]
    if any(set(p)==attached for p in constraints):raise ValueError('Covalent drug attachment was constrained')
    metal_rows=[]
    for site in result.get('initial_metal_sites',[]):
        ds=md.compute_distances(traj,[[site['metal_index'],a['index']] for a in site['donors']],periodic=True)*10
        metal_rows.append({'metal_index':site['metal_index'],'donors':[dict(a,
            sampled_min_angstrom=float(ds[:,i].min()),sampled_max_angstrom=float(ds[:,i].max()),
            sampled_mean_angstrom=float(ds[:,i].mean())) for i,a in enumerate(site['donors'])]})
    old_pdb=app.PDBFile(str(folder/'final.pdb'))
    old_box=np.asarray(old_pdb.topology.getPeriodicBoxVectors().value_in_unit(unit.nanometer))
    corrected=None
    if np.max(abs(old_box-last_box))>1e-4:
        corrected=folder/'final-with-current-box.pdb'
        if corrected.exists():raise ValueError('Preserve existing corrected export')
        topology=app.AmberPrmtopFile(str(inputs/'solvated.prmtop')).topology
        topology.setPeriodicBoxVectors(state.getPeriodicBoxVectors())
        with corrected.open('w') as handle:app.PDBFile.writeFile(topology,state.getPositions(),handle,keepIds=True)
        check=app.PDBFile(str(corrected))
        check_box=np.asarray(check.topology.getPeriodicBoxVectors().value_in_unit(unit.nanometer))
        if np.max(abs(check_box-last_box))>1e-4:raise ValueError('Corrected PDB box differs beyond export precision')
    report={'stage':'saved_run_integrity_and_geometry_checked','case':result['case'],
        'diagnostic_only':diagnostic,
        'qualified_for_handoff':False if diagnostic else result.get('qualified_for_handoff',True),
        'source_result_sha256':digest(folder/'result.json'),
        'loop_geometry_flags_present':bool(result.get('recorded_peptide_excursions')),
        'all_sampled_loop_geometry_checks_passed':result.get('all_sampled_loop_geometry_checks_passed'),
        'physical_model_validated':False,'app_ready':False,'frames':len(traj),'atoms':traj.n_atoms,
        'final_coordinate_error_nm':coordinate_error,'final_box_error_nm':box_error,
        'adduct_bond_ratio_range':[float(ratios.min()),float(ratios.max())],
        'adduct_heavy_angle_max_deviation_degrees':float(np.max(abs(angle_deviations))),
        'adduct_heavy_angle_rms_deviation_degrees':float(np.sqrt(np.mean(angle_deviations**2))),
        'attachment_constraint_absent':True,'attachment_torsions':torsion_rows,
        'sampled_stereocenters':stereo,'metal_sites':metal_rows,
        'corrected_final_pdb':str(corrected) if corrected else None,
        'original_final_pdb_preserved':True,
        'scope':'Independent saved-file integrity and short sampled adduct geometry check; this does not clear recorded loop failures, establish the correct torsion potential or equilibrium ensemble, or promote a diagnostic source to a qualified handoff.',
        'output_sha256':{name:digest(folder/name) for name in ['trajectory.dcd','final-state.xml','system.xml','final.pdb']}}
    if corrected:report['output_sha256'][corrected.name]=digest(corrected)
    (folder/'saved-run-audit.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps(report,indent=2))


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--run',type=Path,required=True);p.add_argument('--inputs',type=Path,required=True)
    a=p.parse_args();analyze(a.run,a.inputs)
