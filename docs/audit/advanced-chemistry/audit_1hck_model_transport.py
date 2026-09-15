"""Inspect every mandatory Hu2024 added term without omitting CMAP/NBFIX.

Static transport evidence only. The ordinary Amber manifest exporter correctly
refuses this model family; this case-specific audit is not production support.
"""
import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
import openmm as mm
from openmm import app, unit
import parmed


def main():
    parser=argparse.ArgumentParser();parser.add_argument('folder',type=Path)
    folder=parser.parse_args().folder.resolve()
    native=parmed.load_file(str(folder/'published.prmtop'))
    raw=native.parm_data
    p=app.AmberPrmtopFile(str(folder/'published.prmtop'))
    system=p.createSystem(nonbondedMethod=app.NoCutoff,constraints=None,rigidWater=False)
    baseline=app.AmberPrmtopFile(str(folder/'baseline.prmtop')).createSystem(
        nonbondedMethod=app.NoCutoff,constraints=None,rigidWater=False)
    forces={type(f).__name__:f for f in system.getForces()}
    if len(forces)!=system.getNumForces():raise ValueError('Unexpected duplicate force class')
    if set(forces)!={'HarmonicBondForce','HarmonicAngleForce','PeriodicTorsionForce',
                     'CMAPTorsionForce','NonbondedForce','CustomNonbondedForce','CMMotionRemover'}:
        raise ValueError('Unreviewed force term')
    common=[]
    for f in baseline.getForces():
        name=type(f).__name__
        if name=='NonbondedForce':continue
        if mm.XmlSerializer.serialize(f)!=mm.XmlSerializer.serialize(forces[name]):
            raise ValueError('Mandatory postprocessor changed an unrelated imported force')
        common.append(name)
    # Native 24x24 map is tabulated from -pi; OpenMM stores 0..2pi with the
    # first angle changing fastest. Compare the entire table, not one energy.
    cmap=forces['CMAPTorsionForce']
    if cmap.getNumMaps()!=1 or cmap.getNumTorsions()!=1:raise ValueError('Missing/additional CMAP')
    size,values=cmap.getMapParameters(0)
    expected=np.roll(np.asarray(raw['CMAP_PARAMETER_01']).reshape(24,24),(-12,-12),axis=(0,1)).T.ravel()*4.184
    actual=np.asarray(values.value_in_unit(unit.kilojoule_per_mole))
    if size!=24 or not np.allclose(actual,expected,rtol=0,atol=1e-12):raise ValueError('CMAP values/order differ')
    a,b,c,d,e,map_index=raw['CMAP_INDEX']
    if list(cmap.getTorsionParameters(0))!=[map_index-1,a-1,b-1,c-1,d-1,b-1,c-1,d-1,e-1]:
        raise ValueError('CMAP uses wrong atom identities or torsion direction')
    custom=forces['CustomNonbondedForce']
    if ''.join(custom.getEnergyFunction().split())!='(a/r6)^2-b/r6;r6=r^6;a=acoef(type1,type2);b=bcoef(type1,type2);':
        raise ValueError('Unreviewed custom potential expression')
    if custom.getNumPerParticleParameters()!=1 or custom.getPerParticleParameterName(0)!='type':
        raise ValueError('Unknown pair type binding')
    if custom.getNumGlobalParameters()!=0 or custom.getNumInteractionGroups()!=0 or custom.getUseSwitchingFunction():
        raise ValueError('Additional unreviewed pair potential settings')
    n=native.ptr('ntypes');tables={}
    for i in range(custom.getNumTabulatedFunctions()):
        name=custom.getTabulatedFunctionName(i);tables[name]=custom.getTabulatedFunction(i).getFunctionParameters()
    if set(tables)!={'acoef','bcoef'}:raise ValueError('Missing/additional pair table')
    table_error={}
    for name,field,power,scale in [('acoef','LENNARD_JONES_ACOEF',.5,np.sqrt(4.184)*1e-6),
                                  ('bcoef','LENNARD_JONES_BCOEF',1.,4.184e-6)]:
        nx,ny,data=tables[name]
        if (nx,ny)!=(n,n):raise ValueError('Pair table dimensions differ')
        expected=np.asarray([raw[field][k-1]**power*scale for k in raw['NONBONDED_PARM_INDEX']])
        actual=np.asarray(data)
        error=float(np.max(np.abs(actual-expected)));table_error[name]=error
        if not np.allclose(actual,expected,rtol=2e-14,atol=1e-16):raise ValueError('Native pair coefficients differ')
    nb=forces['NonbondedForce']
    baseline_nb=next(f for f in baseline.getForces() if isinstance(f,mm.NonbondedForce))
    if custom.getNumParticles()!=len(native.atoms) or system.getNumParticles()!=len(native.atoms):raise ValueError('Particle inventory differs')
    for i,atom in enumerate(native.atoms):
        if list(custom.getParticleParameters(i))!=[float(atom.nb_idx-1)]:raise ValueError('Native pair type bound to wrong particle')
        charge,_,epsilon=nb.getParticleParameters(i)
        if abs(charge.value_in_unit(unit.elementary_charge)-atom.charge)>1e-12 or epsilon.value_in_unit(unit.kilojoule_per_mole)!=0:
            raise ValueError('Charge mismatch or double-counted particle LJ')
        if abs(system.getParticleMass(i).value_in_unit(unit.dalton)-atom.mass)>1e-12:raise ValueError('Native mass differs')
    if nb.getNumExceptions()!=baseline_nb.getNumExceptions():raise ValueError('Exception inventory differs')
    for i in range(nb.getNumExceptions()):
        if nb.getExceptionParameters(i)!=baseline_nb.getExceptionParameters(i):raise ValueError('Nonbonded exception parameters changed')
    expected_exclusions=set();offset=0
    for i,count in enumerate(raw['NUMBER_EXCLUDED_ATOMS']):
        for j in raw['EXCLUDED_ATOMS_LIST'][offset:offset+count]:
            if j:expected_exclusions.add(tuple(sorted((i,j-1))))
        offset+=count
    custom_exclusions={tuple(sorted(custom.getExclusionParticles(i))) for i in range(custom.getNumExclusions())}
    if custom_exclusions!=expected_exclusions:raise ValueError('Custom LJ exclusions differ from native Amber graph')
    # Identity comparison is by actual importer aliases, never just array size.
    source_map=json.loads((folder/'source-atom-map.json').read_text())['atom_map']
    for native_atom,atom,record in zip(native.atoms,p.topology.atoms(),source_map,strict=True):
        resname=app.PDBFile._residueNameReplacements.get(native_atom.residue.name,native_atom.residue.name)
        name=app.PDBFile._atomNameReplacements.get(resname,{}).get(native_atom.name,native_atom.name)
        if (resname,name,native_atom.atomic_number)!=(atom.residue.name,atom.name,atom.element.atomic_number):
            raise ValueError('Native/display identity mismatch')
        if record['native_index']!=native_atom.idx or record['native_atom_name']!=native_atom.name:
            raise ValueError('Source atom map is stale')
    (folder/'published-openmm-system.xml').write_text(mm.XmlSerializer.serialize(system))
    report={'status':'case-specific mandatory-term transport checks passed; not model accuracy or production lifecycle support',
            'native_atoms':len(native.atoms),'forces':list(forces),'common_imported_forces_unchanged':common,
            'cmap':{'maps':1,'size':24,'values_checked':576,'atom_order_checked':True,'maximum_grid_error_kj_mol':float(np.max(np.abs(np.asarray(values.value_in_unit(unit.kilojoule_per_mole))-np.roll(np.asarray(raw['CMAP_PARAMETER_01']).reshape(24,24),(-12,-12),axis=(0,1)).T.ravel()*4.184)))},
            'all_native_pair_table_entries_checked':n*n,'pair_coefficient_maximum_errors':table_error,
            'explicit_exclusion_pairs_checked':len(expected_exclusions),'native_masses_charges_types_identity_checked':True,
            'one_four_and_zero_exceptions_unchanged':nb.getNumExceptions(),'constraints':system.getNumConstraints(),
            'versions':{'openmm':mm.__version__,'parmed':parmed.__version__},
            'sha256':{name:hashlib.sha256((folder/name).read_bytes()).hexdigest() for name in ['baseline.prmtop','published.prmtop','source-atom-map.json','published-openmm-system.xml']},
            'audit_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
            'limitations':'Does not establish scientific accuracy of Mg-protein transfer, equilibrium coordination, phase-space sampling, or bundled availability. Native same-coordinate Sander results are recorded separately.'}
    (folder/'mandatory-term-transport.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps(report,indent=2))


if __name__=='__main__':main()
