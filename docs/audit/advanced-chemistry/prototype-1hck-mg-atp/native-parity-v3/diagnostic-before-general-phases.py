"""Static per-term/native-convention diagnosis; never edits model parameters."""
from pathlib import Path
import hashlib
import json
import math
import numpy as np
import openmm as mm
from openmm import app, unit

BASE=Path(__file__).resolve().parent
ROOT=BASE.parents[2]
MODEL=BASE/'prototype-1hck-mg-atp/native-v2'
OUT=BASE/'prototype-1hck-mg-atp/native-parity-v3'
SOURCE=ROOT/'build/advanced-chemistry/amber-sander-source-check/set.F90'

def sha(path): return hashlib.sha256(path.read_bytes()).hexdigest()
def scalar(value): return float(value.value_in_unit(unit.kilojoule_per_mole))

def summarize(diff):
    if not np.isfinite(diff).all(): raise ValueError('Nonfinite force difference')
    return {'maximum_component_kj_mol_nm':float(np.max(np.abs(diff))),
            'rms_component_kj_mol_nm':float(np.sqrt(np.mean(diff**2))),
            'maximum_vector_kj_mol_nm':float(np.max(np.linalg.norm(diff,axis=1))),
            'maximum_component_index':list(map(int,np.unravel_index(np.argmax(np.abs(diff)),diff.shape)))}

def components(system,poses):
    names=[]
    for i,force in enumerate(system.getForces()):
        force.setForceGroup(i);names.append(f'{i}:{type(force).__name__}')
    integrator=mm.VerletIntegrator(.001)
    context=mm.Context(system,integrator,mm.Platform.getPlatformByName('Reference'))
    rows=[];vectors={}
    try:
        for pose in poses:
            context.setPositions(np.asarray(pose['angstrom'])*unit.angstrom)
            values={};arrays={}
            for i,name in enumerate(names):
                state=context.getState(getEnergy=True,getForces=True,groups=1<<i)
                values[name]=scalar(state.getPotentialEnergy())
                arrays[name]=np.asarray(state.getForces(asNumpy=True).value_in_unit(unit.kilojoule_per_mole/unit.nanometer))
            full=context.getState(getEnergy=True,getForces=True)
            f=np.asarray(full.getForces(asNumpy=True).value_in_unit(unit.kilojoule_per_mole/unit.nanometer))
            e=scalar(full.getPotentialEnergy())
            if not np.isfinite(f).all() or not np.isfinite(e): raise ValueError('Nonfinite model')
            if abs(sum(values.values())-e)>1e-7: raise ValueError('Component sum differs from total')
            if np.max(np.abs(sum(arrays.values())-f))>1e-7: raise ValueError('Force components differ from total')
            rows.append({'pose':pose['name'],'energy_kj_mol':e,'components_kj_mol':values})
            vectors[pose['name']]={'total':f,**arrays}
    finally:
        del context,integrator
    return rows,vectors

def main():
    report_path=OUT/'component-convention-diagnostic.json'
    if report_path.exists():raise FileExistsError('Preserve prior diagnostic')
    before={p.name:sha(p) for p in [MODEL/'published.prmtop',MODEL/'published.inpcrd',MODEL/'published-openmm-system.xml']}
    code=SOURCE.read_text()
    if 'if(abs(dum-pi) <= TEN_TO_MINUS3) dum = sign(pim,dum)' not in code:raise ValueError('Missing documented native phase behavior')
    request=json.loads((OUT/'sander-parity-input.json').read_text())
    native=json.loads((OUT/'sander-parity-native.json').read_text())
    initial=json.loads((OUT/'sander-parity.json').read_text())
    raw=mm.XmlSerializer.deserialize((MODEL/'published-openmm-system.xml').read_text())
    ratio=initial['electrostatic_convention']['energy_scale_ratio']
    matched=mm.XmlSerializer.deserialize(mm.XmlSerializer.serialize(raw))
    nb=next(f for f in matched.getForces() if isinstance(f,mm.NonbondedForce))
    for i in range(nb.getNumParticles()):
        q,s,e=nb.getParticleParameters(i);nb.setParticleParameters(i,q*np.sqrt(ratio),s,e)
    for i in range(nb.getNumExceptions()):
        a,b,q,s,e=nb.getExceptionParameters(i);nb.setExceptionParameters(i,a,b,q*ratio,s,e)
    snapped=mm.XmlSerializer.deserialize(mm.XmlSerializer.serialize(matched))
    torsion=next(f for f in snapped.getForces() if isinstance(f,mm.PeriodicTorsionForce))
    changes=[];bound=0.
    for i in range(torsion.getNumTorsions()):
        a,b,c,d,n,phase,k=torsion.getTorsionParameters(i);p=phase.value_in_unit(unit.radian)
        if abs(p-math.pi)<=.001:
            changes.append({'index':i,'phase_radian':p,'diagnostic_phase_radian':math.pi})
            bound+=abs(p-math.pi)*scalar(k)
            torsion.setTorsionParameters(i,a,b,c,d,n,math.pi,k)
        elif p!=0:raise ValueError(f'Unreviewed nonzero/non-pi phase: {p}')
    arrays={};modes={}
    pairs={'HarmonicBondForce':['bond'],'HarmonicAngleForce':['angle'],
           'PeriodicTorsionForce':['dihedral'],'CMAPTorsionForce':['cmap'],
           'CustomNonbondedForce':['vdw'],'NonbondedForce':['elec','elec_14','vdw_14']}
    for mode,system in [('raw',raw),('native_electrostatics',matched),('native_electrostatics_and_phase',snapped)]:
        rows,forces=components(system,request['poses']);results=[]
        for row,expected in zip(rows,native,strict=True):
            if row['pose']!=expected['name']:raise ValueError('Pose identity mismatch')
            nf=np.asarray(expected['forces_kj_mol_nm']);f=forces[row['pose']]['total']
            d=row['energy_kj_mol']-expected['energy_kj_mol'];fd=summarize(f-nf)
            differences={}
            for key,value in row['components_kj_mol'].items():
                kind=key.split(':',1)[1]
                if kind in pairs:
                    n=sum(expected['components_kj_mol'][x] for x in pairs[kind])
                    differences[kind]={'openmm':value,'sander':n,'difference_kj_mol':value-n}
            results.append({**row,'sander_energy_kj_mol':expected['energy_kj_mol'],'difference_kj_mol':d,
                'force_difference':fd,'component_comparison':differences,
                'passed_original_tolerances':bool(abs(d)<.001 and fd['maximum_component_kj_mol_nm']<.01)})
            for kind,values in forces[row['pose']].items(): arrays[f'{mode}/{row["pose"]}/{kind}']=values
            arrays[f'sander/{row["pose"]}/total']=nf
        modes[mode]=results
    after={p.name:sha(p) for p in [MODEL/'published.prmtop',MODEL/'published.inpcrd',MODEL/'published-openmm-system.xml']}
    if before!=after:raise ValueError('Source model changed during diagnostic')
    force_path=OUT/'full-force-components.npz';np.savez_compressed(force_path,**arrays)
    report={'scope':'Identical complete static models, three input poses, native process per pose; no dynamics or chemistry accuracy inference.',
        'energy_tolerance_kj_mol':.001,'maximum_force_component_tolerance_kj_mol_nm':.01,
        'electrostatic_convention':initial['electrostatic_convention'],
        'phase_source':{'path':str(SOURCE.relative_to(ROOT)),'sha256':sha(SOURCE),'routine':'dihpar','lines':'433-437'},
        'changed_phase_terms':len(changes),'phase_deltas':sorted(set(x['phase_radian']-x['diagnostic_phase_radian'] for x in changes)),
        'coordinate_independent_phase_energy_difference_bound_kj_mol':bound,'modes':modes,
        'model_sha256_unchanged':before,'script_sha256':sha(Path(__file__)),
        'native_evidence_sha256':{name:sha(OUT/name) for name in ['sander-parity-input.json','sander-parity-native.json','sander-parity.json']},
        'full_force_components':{'path':force_path.name,'sha256':sha(force_path)}}
    report_path.write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps({mode:[{'pose':x['pose'],'energy_difference':x['difference_kj_mol'],'force_difference':x['force_difference'],'components':x['component_comparison'],'passed':x['passed_original_tolerances']} for x in values] for mode,values in modes.items()},indent=2))

if __name__=='__main__':main()
