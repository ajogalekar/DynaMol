"""Explain a native Sander/OpenMM near-pi convention difference independently.

Only a diagnostic in-memory System is changed. Native files, application
parameters, raw failed parity report and actual motion system stay unchanged.
"""
from pathlib import Path
import hashlib
import json
import math
import numpy as np
import openmm as mm
from openmm import app,unit

BASE=Path(__file__).resolve().parent
ROOT=BASE.parents[2]
folder=BASE/'prototype-1ca2-zaff-water-v3'
source=ROOT/'build/advanced-chemistry/amber-sander-source-check/set.F90'
code=source.read_text()
if 'if(abs(dum-pi) <= TEN_TO_MINUS3) dum = sign(pim,dum)' not in code:
    raise ValueError('Native phase-snap source evidence is missing')
initial=json.loads((folder/'sander-parity.json').read_text())
request=json.loads((folder/'sander-parity-input.json').read_text())
reference=json.loads((folder/'sander-parity-native.json').read_text())
native=app.AmberPrmtopFile(str(folder/'prepared.prmtop'))
system=native.createSystem(nonbondedMethod=app.NoCutoff,constraints=None,rigidWater=False)
ratio=initial['electrostatic_convention']['energy_scale_ratio']
nb=next(f for f in system.getForces() if isinstance(f,mm.NonbondedForce))
for i in range(nb.getNumParticles()):
    q,s,e=nb.getParticleParameters(i);nb.setParticleParameters(i,q*np.sqrt(ratio),s,e)
for i in range(nb.getNumExceptions()):
    a,b,q,s,e=nb.getExceptionParameters(i);nb.setExceptionParameters(i,a,b,q*ratio,s,e)
changed=[];bound=0.
torsions=next(f for f in system.getForces() if isinstance(f,mm.PeriodicTorsionForce))
for i in range(torsions.getNumTorsions()):
    a,b,c,d,n,phase,k=torsions.getTorsionParameters(i)
    p=phase.value_in_unit(unit.radian)
    if abs(p-math.pi)<=.001:
        bound+=abs(p-math.pi)*k.value_in_unit(unit.kilojoule_per_mole)
        changed.append({'term_index':i,'phase_original_radians':p,'phase_diagnostic_radians':math.pi})
        torsions.setTorsionParameters(i,a,b,c,d,n,math.pi,k)
    elif p!=0:
        raise ValueError('This bounded diagnostic is limited to actual zero/near-pi phases')
integrator=mm.VerletIntegrator(.001)
context=mm.Context(system,integrator,mm.Platform.getPlatformByName('Reference'))
records=[]
for pose,expected in zip(request['poses'],reference,strict=True):
    context.setPositions(np.asarray(pose['angstrom'])*unit.angstrom)
    state=context.getState(getEnergy=True,getForces=True)
    energy=state.getPotentialEnergy().value_in_unit(unit.kilojoule_per_mole)
    forces=np.asarray(state.getForces(asNumpy=True).value_in_unit(unit.kilojoule_per_mole/unit.nanometer))
    de=float(energy-expected['energy_kj_mol']);df=float(np.max(np.abs(forces-np.asarray(expected['forces_kj_mol_nm']))))
    records.append({'pose':pose['name'],'energy_difference_kj_mol':de,'maximum_force_difference_kj_mol_nm':df,
        'accepted_with_original_tolerances':bool(abs(de)<.001 and df<.01)})
del context,integrator
report={'scope':'Native implementation numerical-convention diagnostic only; production parameters and raw failed comparison unchanged',
    'native_source':{'path':str(source.relative_to(ROOT)),'sha256':hashlib.sha256(source.read_bytes()).hexdigest(),
        'archive':'build/loop-repair-review/new-source-materials/archives/52fb4fb3370a-AmberTools24_rc5.tar.bz2',
        'routine':'dihpar; lines433–437','behavior':'If phase is within 0.001 radians of pi, evaluate at exact mathematical pi; zero sine/cosine components whose absolute value is at most 1e-6'},
    'prior_report':{'path':'sander-parity.json','sha256':hashlib.sha256((folder/'sander-parity.json').read_bytes()).hexdigest(),
        'passed':initial['passed']},'phase_terms_matched':len(changed),'phase_change_radians':sorted(set(x['phase_original_radians']-x['phase_diagnostic_radians'] for x in changed)),
    'global_energy_difference_upper_bound_from_phase_rounding_kj_mol':bound,
    'bound_method':'Triangle inequality: sum_i amplitude_i * abs(delta_phase_i), using |cos(a)-cos(b)| <= |a-b|. Independent of tested coordinates.',
    'electrostatic_convention':initial['electrostatic_convention'],'conformations':records,
    'accepted_after_matching_documented_native_conventions':all(r['accepted_with_original_tolerances'] for r in records),
    'production_parameter_files_unchanged':{n:hashlib.sha256((folder/n).read_bytes()).hexdigest() for n in ('prepared.prmtop','prepared.inpcrd')}}
(folder/'sander-native-phase-diagnostic.json').write_text(json.dumps(report,indent=2)+'\n')
print(json.dumps(report,indent=2))
