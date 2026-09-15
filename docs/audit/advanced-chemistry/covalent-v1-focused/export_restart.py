"""Export a portable physical restart without re-enabling warm-up restraints."""
import argparse
import json
from pathlib import Path

import numpy as np
import openmm as mm
from openmm import unit

from prepare_adduct import digest


def verify_roundtrip(output):
    report=json.loads((output/'manifest.json').read_text())
    for name,expected in report['exported_artifacts'].items():
        if digest(output/name)!=expected:raise ValueError('Portable restart artifact changed')
    system=mm.XmlSerializer.deserialize((output/'system.xml').read_text())
    state=mm.XmlSerializer.deserialize((output/'state.xml').read_text())
    integrator=mm.XmlSerializer.deserialize((output/'integrator.xml').read_text())
    context=mm.Context(system,integrator,mm.Platform.getPlatformByName('Reference'))
    context.setState(state)
    fresh=context.getState(getEnergy=True,getForces=True,getPositions=True,getParameters=True)
    if context.getParameter('position_k')!=0:raise ValueError('Reloaded restart has nonzero warm-up restraint')
    energy=float(fresh.getPotentialEnergy().value_in_unit(unit.kilojoule_per_mole))
    forces=np.asarray(fresh.getForces(asNumpy=True).value_in_unit(unit.kilojoule_per_mole/unit.nanometer))
    reference=np.asarray(state.getForces(asNumpy=True).value_in_unit(unit.kilojoule_per_mole/unit.nanometer))
    difference=float(np.max(abs(forces-reference)))
    delta=energy-report['comparisons'][-1]['energy_kj_mol']
    if abs(delta)>1e-5 or difference>1e-4:raise ValueError('Reloaded restart mechanics changed')
    audit={'stage':'serialized_restart_roundtrip_verified','warmup_restraint_k':context.getParameter('position_k'),
           'diagnostic_only':bool(report.get('diagnostic_only')),
           'qualified_for_handoff':False if report.get('diagnostic_only') else report.get('qualified_for_handoff',True),
           'energy_difference_kj_mol':delta,'maximum_force_component_difference_kj_mol_nm':difference,
           'artifact_sha256':report['exported_artifacts']}
    (output/'roundtrip.json').write_text(json.dumps(audit,indent=2)+'\n')
    del context,integrator
    return audit


def export(folder,output):
    result=json.loads((folder/'result.json').read_text())
    diagnostic=bool(result.get('diagnostic_only'))
    if diagnostic and result.get('qualified_for_handoff') is not False:
        raise ValueError('Diagnostic source has inconsistent handoff qualification metadata')
    value=result['final']['position_restraint_k']
    if result['final']['stage']!='unrestrained NPT stability' or value!=0:
        raise ValueError('This export requires a recorded unrestrained endpoint')
    output.mkdir(exist_ok=False)
    original_xml=(folder/'system.xml').read_text()
    old_state=mm.XmlSerializer.deserialize((folder/'final-state.xml').read_text())
    evidence=[];physical_xml=None;restart_state=None;reference_force=None;reference_energy=None
    for changed in [False,True]:
        system=mm.XmlSerializer.deserialize(original_xml)
        matched=[]
        for force in system.getForces():
            if not isinstance(force,mm.CustomExternalForce):continue
            for i in range(force.getNumGlobalParameters()):
                if force.getGlobalParameterName(i)=='position_k':
                    matched.append(force.getGlobalParameterDefaultValue(i))
                    if changed:force.setGlobalParameterDefaultValue(i,value)
        if len(matched)!=1:raise ValueError('Expected exactly one known warm-up restraint parameter')
        integrator=mm.LangevinMiddleIntegrator(300*unit.kelvin,1/unit.picosecond,.002*unit.picosecond)
        integrator.setRandomNumberSeed(2027)
        context=mm.Context(system,integrator,mm.Platform.getPlatformByName('Reference'))
        context.setState(old_state)
        if not changed:context.setParameter('position_k',value)
        state=context.getState(getPositions=True,getVelocities=True,getEnergy=True,getForces=True,getParameters=True)
        if context.getParameter('position_k')!=0:raise ValueError('Restart silently re-enabled a construction restraint')
        force=np.asarray(state.getForces(asNumpy=True).value_in_unit(unit.kilojoule_per_mole/unit.nanometer))
        energy=float(state.getPotentialEnergy().value_in_unit(unit.kilojoule_per_mole))
        if not np.isfinite(force).all() or not np.isfinite(energy):raise ValueError('Nonfinite restart mechanics')
        evidence.append({'corrected_system_defaults':changed,'original_system_default':matched[0],
                         'context_parameters':dict(state.getParameters()),'energy_kj_mol':energy})
        (output/'comparisons.json').write_text(json.dumps(evidence,indent=2)+'\n')
        if changed:
            (output/'differences.json').write_text(json.dumps({'energy_difference_kj_mol':energy-reference_energy,
                'maximum_force_component_difference_kj_mol_nm':float(np.max(abs(force-reference_force)))},indent=2)+'\n')
            if abs(energy-reference_energy)>1e-5 or np.max(abs(force-reference_force))>1e-4:
                raise ValueError('Changing the serialized default changed physical endpoint mechanics')
            position=np.asarray(state.getPositions(asNumpy=True).value_in_unit(unit.nanometer))
            original=np.asarray(old_state.getPositions(asNumpy=True).value_in_unit(unit.nanometer))
            if not np.array_equal(position,original):raise ValueError('Restart export changed coordinates')
            physical_xml=mm.XmlSerializer.serialize(system);restart_state=mm.XmlSerializer.serialize(state)
            (output/'integrator.xml').write_text(mm.XmlSerializer.serialize(integrator))
        else:reference_force=force;reference_energy=energy
        del context,integrator
    (output/'system.xml').write_text(physical_xml);(output/'state.xml').write_text(restart_state)
    report={'stage':'portable_physical_restart_verified','physical_model_validated':False,'app_ready':False,
        'diagnostic_only':diagnostic,
        'qualified_for_handoff':False if diagnostic else result.get('qualified_for_handoff',True),
        'parameters_included_in_state':True,'warmup_restraint_default_corrected':True,
        'coordinates_preserved_exactly':True,'comparisons':evidence,
        'comparison_platform':'OpenMM Reference',
        'scope':'Same physical endpoint and zero warm-up restraint, with a new explicitly seeded random stream. Not a bitwise RNG continuation; the original binary checkpoint is preserved separately. Export does not clear source geometry failures or diagnostic-only restrictions.',
        'original_artifacts':{n:digest(folder/n) for n in ['system.xml','final-state.xml','result.json']},
        'exported_artifacts':{n:digest(output/n) for n in ['system.xml','state.xml','integrator.xml']}}
    (output/'manifest.json').write_text(json.dumps(report,indent=2)+'\n')
    report['roundtrip']=verify_roundtrip(output)
    print(json.dumps(report,indent=2))


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--run',type=Path,required=True);p.add_argument('--output',type=Path,required=True)
    a=p.parse_args();export(a.run,a.output)
