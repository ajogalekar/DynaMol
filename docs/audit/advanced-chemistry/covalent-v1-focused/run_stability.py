"""Bounded 100 ps numerical/site check of a complete native candidate.

This is a warm-up and stability test, not production MD, model-accuracy
validation, or a substitute for assessing uncertain attachment torsions.
"""
from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path
import time

import numpy as np
import openmm as mm
from openmm import app, unit
import parmed
import psutil

from prepare_adduct import digest


def run(folder,output,seed,loop_protein=None,free_remodeled_region=False,diagnostic_peptide_excursions=False,unrestrained_warmup=False):
    if free_remodeled_region and loop_protein is None:
        raise ValueError('Freeing a remodeled region requires its mapped, screened protein intermediate')
    if diagnostic_peptide_excursions and loop_protein is None:
        raise ValueError('Peptide excursion diagnostics require the complete mapped loop monitor')
    if unrestrained_warmup and not diagnostic_peptide_excursions:
        raise ValueError('The unrestrained-warm-up comparison is currently diagnostic-only')
    output.mkdir(parents=True,exist_ok=False)
    start=time.monotonic(); process=psutil.Process()
    evidence=json.loads((folder/'result.json').read_text())
    (output/'launch.json').write_text(json.dumps({'started_unix':time.time(),'pid':process.pid,
        'runner':str(Path(__file__).resolve()),'runner_sha256':digest(__file__),
        'input':str(folder),'seed':seed,'time_limit_seconds':7200,'rss_limit_bytes':4_000_000_000},indent=2)+'\n')
    native=parmed.load_file(str(folder/'solvated.prmtop'),xyz=str(folder/'solvated.inpcrd'))
    from native_elements import inventory
    inventory(folder/'solvated.prmtop')
    prmtop=app.AmberPrmtopFile(str(folder/'solvated.prmtop'))
    inpcrd=app.AmberInpcrdFile(str(folder/'solvated.inpcrd'))
    system=prmtop.createSystem(nonbondedMethod=app.PME,nonbondedCutoff=1.0*unit.nanometer,
                              constraints=app.HBonds,rigidWater=True,ewaldErrorTolerance=1e-5)
    loop_monitor=None
    if loop_protein is not None:
        from check_complex_loops import LoopMonitor
        loop_monitor=LoopMonitor(folder,loop_protein,prmtop.topology,diagnostic_peptide_excursions)
        allowed={'HarmonicBondForce','HarmonicAngleForce','PeriodicTorsionForce','NonbondedForce','CMMotionRemover'}
        if any(type(force).__name__ not in allowed for force in system.getForces()):
            raise ValueError('Unexpected force class in the native model before adding MD warm-up controls')
        (output/'loop-monitor-input.json').write_text(json.dumps({'protein':str(loop_protein),
            'protein_result_sha256':digest(loop_protein/'result.json'),
            'native_force_classes':[type(force).__name__ for force in system.getForces()],
            'free_remodeled_region_during_warmup':free_remodeled_region,
            'diagnostic_only':diagnostic_peptide_excursions,
            'warmup_restraints_released_before_heating':unrestrained_warmup,
            'remodeled_native_atom_indices':sorted(loop_monitor.remodeled_atom_indices),
            'source_to_native_residues':loop_monitor.source_mapping,
            'temporary_construction_forces_absent':True},indent=2)+'\n')
    xyz=np.asarray(inpcrd.positions.value_in_unit(unit.nanometer))
    restraint=mm.CustomExternalForce('0.5*position_k*((x-x0)^2+(y-y0)^2+(z-z0)^2)')
    restraint.addGlobalParameter('position_k',1000.)
    for name in ['x0','y0','z0']:restraint.addPerParticleParameter(name)
    for atom in native.atoms:
        if (atom.atomic_number>1 and atom.residue.name not in ['WAT','Na+','Cl-']
                and not (free_remodeled_region and atom.idx in loop_monitor.remodeled_atom_indices)):
            restraint.addParticle(atom.idx,xyz[atom.idx])
    restraint.setForceGroup(31);system.addForce(restraint)
    barostat=mm.MonteCarloBarostat(1*unit.bar,300*unit.kelvin,0)
    barostat.setRandomNumberSeed(seed+1);system.addForce(barostat)
    integrator=mm.LangevinMiddleIntegrator(50*unit.kelvin,1/unit.picosecond,.002*unit.picosecond)
    integrator.setRandomNumberSeed(seed)
    simulation=app.Simulation(prmtop.topology,system,integrator,mm.Platform.getPlatformByName('CPU'),
                              {'Threads':'2','DeterministicForces':'true'})
    simulation.context.setPositions(inpcrd.positions)
    simulation.context.setPeriodicBoxVectors(*inpcrd.boxVectors)
    cov=next(r for r in native.residues if r.name=='COV')
    edges=[(b.atom1.idx,b.atom2.idx,b.type.req*.1) for b in native.bonds
           if b.atom1.residue is cov and b.atom2.residue is cov]
    attachment=[(a,b,r) for a,b,r in edges if {native.atoms[a].atomic_number,native.atoms[b].atomic_number}=={6,16}
                and {native.atoms[a].name,native.atoms[b].name}!={'CB','SG'}]
    if len(attachment)!=1:raise ValueError('Expected one cysteine–drug carbon–sulfur bond')
    centers=[]
    ca=next(a for a in cov.atoms if a.name=='CA')
    centers.append((ca.idx,[next(a.idx for a in cov.atoms if a.name==n) for n in ['N','C','CB']]))
    # Ligand stereocenters are explicitly identified by the product-graph audit.
    mapping=json.loads((folder/'source-mapping.json').read_text())
    source_stereo={'5P9J':'CBE','6OIM':'C20'}
    if evidence['case'] in source_stereo:
        matched=[r for r in mapping if r['source'].get('atom')==source_stereo[evidence['case']]
                 and r['source'].get('resname') not in ['CYS',None]]
        if len(matched)!=1:raise ValueError('Ligand stereocenter mapping is ambiguous')
        atom=native.atoms[matched[0]['native_index']]
        neighbors=sorted(a.idx for a in atom.bond_partners if a.atomic_number>1)
        if len(neighbors)!=3:raise ValueError('Expected three heavy substituents at mapped stereocenter')
        centers.append((atom.idx,neighbors))
    def volumes(coords):
        return np.array([np.linalg.det(coords[neighbors]-coords[center]) for center,neighbors in centers])
    reference_volumes=volumes(xyz)
    if np.min(abs(reference_volumes))<1e-4:raise ValueError('Initial stereocenter is nearly planar')
    ca_indices=[a.idx for a in native.atoms if a.name=='CA' and a.residue.name!='WAT']
    reference_ca=xyz[ca_indices]-xyz[ca_indices].mean(axis=0)
    dof=3*system.getNumParticles()-system.getNumConstraints()-3
    mass=sum(a.mass for a in native.atoms)
    metals=[a for a in native.atoms if a.atomic_number==12]
    metal_sites=[]
    for metal in metals:
        donors=[a for a in native.atoms[:evidence['dry_atoms']] if a.atomic_number in [7,8]
                and np.linalg.norm(xyz[a.idx]-xyz[metal.idx])<.28]
        metal_sites.append({'metal_index':metal.idx,'donors':[{'index':a.idx,'name':a.name,
            'residue_index':a.residue.idx,'residue_name':a.residue.name,
            'initial_distance_angstrom':float(np.linalg.norm(xyz[a.idx]-xyz[metal.idx])*10)} for a in donors]})
    (output/'metal-sites.json').write_text(json.dumps(metal_sites,indent=2)+'\n')
    def observe(stage,ps):
        state=simulation.context.getState(getPositions=True,getEnergy=True,getForces=True)
        positions=np.asarray(state.getPositions(asNumpy=True).value_in_unit(unit.nanometer))
        force=np.asarray(state.getForces(asNumpy=True).value_in_unit(unit.kilojoule_per_mole/unit.nanometer))
        energy=float(state.getPotentialEnergy().value_in_unit(unit.kilojoule_per_mole))
        if not np.isfinite(positions).all() or not np.isfinite(force).all() or not np.isfinite(energy):
            raise ValueError('Nonfinite trajectory coordinates, energy or forces')
        ratios=[np.linalg.norm(positions[a]-positions[b])/r for a,b,r in edges]
        if min(ratios)<.6 or max(ratios)>1.5:raise ValueError('An adduct bond is grossly distorted')
        current_volumes=volumes(positions)
        if np.any(current_volumes*reference_volumes<=0):raise ValueError('A mapped stereocenter inverted')
        centered=positions[ca_indices]-positions[ca_indices].mean(axis=0)
        u,s,vt=np.linalg.svd(centered.T@reference_ca)
        correction=np.diag([1,1,np.linalg.det(u@vt)])
        rmsd=float(np.sqrt(np.mean(np.sum((centered@u@correction@vt-reference_ca)**2,axis=1)))*10)
        volume=float(state.getPeriodicBoxVolume().value_in_unit(unit.nanometer**3))
        kinetic=float(state.getKineticEnergy().value_in_unit(unit.kilojoule_per_mole))
        a,b,_=attachment[0]
        row={'stage':stage,'time_ps':ps,'potential_kj_mol':energy,
             'temperature_k':2*kinetic/(dof*.008314462618),
             'density_g_ml':mass*.00166053906660/volume,
             'attachment_angstrom':float(np.linalg.norm(positions[a]-positions[b])*10),
             'adduct_bond_min_ratio':min(ratios),'adduct_bond_max_ratio':max(ratios),
             'ca_rmsd_angstrom':rmsd,'max_force_kj_mol_nm':float(np.max(np.linalg.norm(force,axis=1))),
             'elapsed_seconds':time.monotonic()-start,'rss_bytes':process.memory_info().rss,
             'position_restraint_k':simulation.context.getParameter('position_k')}
        if metal_sites:
            donor_distances=[float(np.linalg.norm(positions[site['metal_index']]-positions[a['index']])*10)
                             for site in metal_sites for a in site['donors']]
            row['initial_metal_donor_min_angstrom']=min(donor_distances)
            row['initial_metal_donor_max_angstrom']=max(donor_distances)
        if loop_monitor is not None:
            try:
                row.update(loop_monitor.check(state.getPositions(),output/'loop-checks',stage,ps))
            except Exception:
                failed=simulation.context.getState(getPositions=True,getVelocities=True,getEnergy=True,getParameters=True)
                (output/'rejected-state.xml').write_text(mm.XmlSerializer.serialize(failed))
                (output/'rejected-system.xml').write_text(mm.XmlSerializer.serialize(system))
                simulation.saveCheckpoint(str(output/'rejected-checkpoint.chk'))
                raise
        if row['elapsed_seconds']>7200 or row['rss_bytes']>4_000_000_000:
            raise TimeoutError('Stability-run time or memory bound exceeded')
        (output/'progress.json').write_text(json.dumps(row,indent=2)+'\n')
        return row
    initial=observe('before minimization',0)
    print('Minimizing complete solvated complex',flush=True)
    simulation.minimizeEnergy(tolerance=10*unit.kilojoule_per_mole/unit.nanometer,maxIterations=1500)
    minimized=observe('minimized',0)
    if unrestrained_warmup:
        simulation.context.setParameter('position_k',0.)
    simulation.context.setVelocitiesToTemperature(50*unit.kelvin,seed)
    simulation.reporters.append(app.DCDReporter(str(output/'trajectory.dcd'),500,enforcePeriodicBox=False))
    rows=[]
    with (output/'observations.csv').open('w') as handle:
        writer=csv.DictWriter(handle,fieldnames=list(minimized));writer.writeheader()
        for picosecond in range(1,101):
            if picosecond<=10:
                temperature=50+250*(picosecond/10)
                integrator.setTemperature(temperature*unit.kelvin)
                stage='NVT temperature ramp'
            elif picosecond<=30:
                if picosecond==11:
                    if not unrestrained_warmup:simulation.context.setParameter('position_k',100.)
                    barostat.setFrequency(25)
                stage='unrestrained NPT equilibration' if unrestrained_warmup else 'restrained NPT equilibration'
            else:
                if picosecond==31:simulation.context.setParameter('position_k',0.)
                stage='unrestrained NPT stability'
            simulation.step(500)
            row=observe(stage,picosecond);rows.append(row);writer.writerow(row);handle.flush()
            if picosecond%10==0:
                print(f"{picosecond}/100 ps; T={row['temperature_k']:.1f} K; density={row['density_g_ml']:.3f}; attachment={row['attachment_angstrom']:.3f} A",flush=True)
                simulation.saveCheckpoint(str(output/'checkpoint.chk'))
    state=simulation.context.getState(getPositions=True,getVelocities=True,getEnergy=True,getParameters=True)
    # Context updates do not change a System's serialized default parameters.
    # A portable restart must not silently restore the warm-up restraint.
    restraint.setGlobalParameterDefaultValue(0,simulation.context.getParameter('position_k'))
    (output/'final-state.xml').write_text(mm.XmlSerializer.serialize(state))
    (output/'system.xml').write_text(mm.XmlSerializer.serialize(system))
    simulation.topology.setPeriodicBoxVectors(state.getPeriodicBoxVectors())
    with (output/'final.pdb').open('w') as handle:app.PDBFile.writeFile(simulation.topology,state.getPositions(),handle)
    report={'stage':'diagnostic_trajectory_complete' if diagnostic_peptide_excursions else 'short_stability_check_complete',
        'case':evidence['case'],'physical_model_validated':False,
        'app_ready':False,'production_simulation':False,'total_ps':100,'unrestrained_ps':100 if unrestrained_warmup else 70,
        'protocol':('10 ps 50–300 K NVT ramp; 20 ps unrestrained NPT equilibration; 70 ps unrestrained NPT. Positional restraints released after the same initial minimization, before heating. '
                    if unrestrained_warmup else '10 ps 50–300 K NVT ramp; 20 ps restrained NPT; 70 ps unrestrained NPT. ')
                   +'1 bar, 2 fs, PME, H-bond constraints. No constraint or bond restraint on the drug attachment.',
        'seed':seed,'initial':initial,'minimized':minimized,'final':rows[-1],
        'unrestrained_temperature_mean_k':float(np.mean([r['temperature_k'] for r in (rows if unrestrained_warmup else rows[30:])])),
        'unrestrained_temperature_mean_window_ps':[1 if unrestrained_warmup else 31,100],
        'stability_phase_temperature_mean_k':float(np.mean([r['temperature_k'] for r in rows[30:]])),
        'attachment_range_angstrom':[min(r['attachment_angstrom'] for r in rows),max(r['attachment_angstrom'] for r in rows)],
        'stereocenter_indices':centers,'stereochemistry_preserved_at_saved_frames':True,
        'initial_metal_sites':metal_sites,
        'full_complex_loop_checks':loop_monitor is not None,
        'free_remodeled_region_during_warmup':free_remodeled_region,
        'warmup_restraints_released_before_heating':unrestrained_warmup,
        'diagnostic_only':diagnostic_peptide_excursions,
        'qualified_for_handoff':not diagnostic_peptide_excursions,
        'recorded_peptide_excursions':loop_monitor.recorded_excursions if loop_monitor is not None else [],
        'all_sampled_loop_geometry_checks_passed':not loop_monitor.recorded_excursions if loop_monitor is not None else None,
        'maximum_ca_rmsd_angstrom':max(r['ca_rmsd_angstrom'] for r in rows),
        'input_sha256':{n:digest(folder/n) for n in ['solvated.prmtop','solvated.inpcrd','source-mapping.json']},
        'elapsed_seconds':time.monotonic()-start}
    (output/'result.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps(report,indent=2),flush=True)


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--input',type=Path,required=True)
    p.add_argument('--output',type=Path,required=True);p.add_argument('--seed',type=int,default=2026)
    p.add_argument('--loop-protein',type=Path)
    p.add_argument('--free-remodeled-region',action='store_true')
    p.add_argument('--diagnostic-peptide-excursions',action='store_true')
    p.add_argument('--unrestrained-warmup',action='store_true')
    a=p.parse_args()
    try:run(a.input,a.output,a.seed,a.loop_protein,a.free_remodeled_region,a.diagnostic_peptide_excursions,a.unrestrained_warmup)
    except Exception as exc:
        if a.output.is_dir():(a.output/'failure.json').write_text(json.dumps({'type':type(exc).__name__,'error':str(exc)},indent=2)+'\n')
        raise
