"""Isolate the preserved published GDP model from protein parameter namespaces."""
import argparse
import copy
import json
from pathlib import Path

import numpy as np
import openmm as mm
from openmm import app,unit
import parmed

from prepare_adduct import native,digest
from assemble_complex import pdb_residues


def evaluate(prmtop,positions):
    topology=app.AmberPrmtopFile(str(prmtop))
    system=topology.createSystem(nonbondedMethod=app.NoCutoff,constraints=None,rigidWater=False)
    integrator=mm.VerletIntegrator(.001)
    context=mm.Context(system,integrator,mm.Platform.getPlatformByName('Reference'))
    context.setPositions(positions*unit.nanometer)
    state=context.getState(getEnergy=True,getForces=True)
    energy=float(state.getPotentialEnergy().value_in_unit(unit.kilojoule_per_mole))
    force=np.array(state.getForces(asNumpy=True).value_in_unit(unit.kilojoule_per_mole/unit.nanometer))
    del context,integrator
    return energy,force


def run(root,amber,out):
    out.mkdir(parents=True,exist_ok=False)
    source=root/'gdp-native-source'
    original=parmed.load_file(str(source/'gdp.prmtop'),xyz=str(source/'gdp.inpcrd'))
    scoped=copy.deepcopy(original)
    types={old:'Q'+chr(ord('A')+i) for i,old in enumerate(sorted({a.type for a in scoped.atoms}))}
    original_types=[a.type for a in scoped.atoms]
    for atom,old in zip(scoped.atoms,original_types):
        atom.type=types[old];atom.atom_type.name=types[old]
    parameters=parmed.amber.AmberParameterSet.from_structure(scoped)
    parameters.write(str(out/'gdp-scoped.frcmod'))
    # Preserve explicit improper ordering and source amplitude precision.
    # This source uses 0/180-degree phases. Native LEaP's angle conversion is
    # already present in the prmtop; reapplying its rounded-radian error as
    # an input-degree adjustment would count that conversion twice.
    exported=[];section=None
    for line in (out/'gdp-scoped.frcmod').read_text().splitlines():
        if line.strip() in ['MASS','BOND','ANGLE','DIHE','IMPROPER','NONB']:
            section=line.strip()
            if section=='IMPROPER':
                exported.append(line)
                for term in original.dihedrals:
                    if not term.improper:continue
                    key='-'.join(types[a.type] for a in [term.atom1,term.atom2,term.atom3,term.atom4])
                    if min(abs(term.type.phase),abs(term.type.phase-180))>.001:
                        raise ValueError('Noncanonical source phase needs explicit conversion qualification')
                    exported.append(f'{key}    {term.type.phi_k:.12f} {round(term.type.phase):.3f} {term.type.per:.1f}')
                exported.append('')
                continue
        elif section=='IMPROPER':
            continue
        elif section=='DIHE' and line.strip():
            fields=line.split();key=tuple(fields[0].split('-'))
            candidates=parameters.dihedral_types[key]
            if not isinstance(candidates,list):candidates=[candidates]
            matched=[t for t in candidates if t.per==abs(float(fields[4]))
                and abs(t.phi_k-float(fields[2])/float(fields[1]))<1e-7
                and abs(t.phase-float(fields[3]))<.001]
            if len(matched)!=1:raise ValueError('Ambiguous source proper torsion during exact export')
            fields[2]=f'{matched[0].phi_k*float(fields[1]):.12f}'
            if min(abs(matched[0].phase),abs(matched[0].phase-180))>.001:
                raise ValueError('Noncanonical source phase needs explicit conversion qualification')
            fields[3]=f'{round(matched[0].phase):.3f}';line=' '.join(fields)
        exported.append(line)
    (out/'gdp-scoped.frcmod').write_text('\n'.join(exported)+'\n')
    lines=[]
    for line in (root/'gdp-source/GDP.prep').read_text().splitlines():
        fields=line.split()
        if len(fields)==11 and fields[0].isdigit() and fields[2] in types:
            fields[2]=types[fields[2]];line=' '.join(fields)
        if line.strip()=='gdp  INT     1':line='GDP INT 1'
        lines.append(line)
    (out/'gdp-scoped.prep').write_text('\n'.join(lines)+'\n')
    from native_elements import inventory, write_unit_elements
    expected_elements = write_unit_elements(source/'gdp.prmtop', 'GDP', out/'gdp-elements.leap')
    (out/'leap.in').write_text('loadamberparams gdp-scoped.frcmod\nloadamberprep gdp-scoped.prep\nsource gdp-elements.leap\n'
        'check GDP\nsaveamberparm GDP gdp.prmtop gdp.inpcrd\nquit\n')
    commands=[native(amber,out,['tleap','-f','leap.in'])]
    if inventory(out/'gdp.prmtop') != expected_elements:
        raise ValueError('Native scoped GDP lost its explicit source elements')
    result=parmed.load_file(str(out/'gdp.prmtop'),xyz=str(out/'gdp.inpcrd'))
    if [(a.name,a.atomic_number) for a in result.atoms]!=[(a.name,a.atomic_number) for a in original.atoms]:
        raise ValueError('Scoping changed GDP atom inventory/elements')
    if not np.allclose([a.charge for a in result.atoms],[a.charge for a in original.atoms],rtol=0,atol=1e-9):
        raise ValueError('Scoping changed source GDP charges')
    environment=pdb_residues(root/'6oim-protein-v1/retained-environment.pdb')
    gdp=[lines for key,lines in environment.items() if key[3]=='GDP']
    if len(gdp)!=1:raise ValueError('Expected one retained GDP')
    observed={line[12:16].strip():np.array([float(line[a:b]) for a,b in [(30,38),(38,46),(46,54)]])*.1
              for line in gdp[0] if line[76:78].strip()!='H'}
    expected={a.name.replace('*',"'") for a in original.atoms if a.atomic_number>1}
    if expected!=set(observed):raise ValueError('Source GDP graph does not match every observed heavy atom')
    xyz=original.coordinates*.1
    for atom in original.atoms:
        if atom.atomic_number>1:xyz[atom.idx]=observed[atom.name.replace('*',"'")]
    for atom in original.atoms:
        if atom.atomic_number==1:
            parent=atom.bond_partners[0]
            xyz[atom.idx]=xyz[parent.idx]+(original.coordinates[atom.idx]-original.coordinates[parent.idx])*.1
    topology=app.AmberPrmtopFile(str(source/'gdp.prmtop'))
    system=topology.createSystem(nonbondedMethod=app.NoCutoff,constraints=None)
    for atom in original.atoms:
        if atom.atomic_number>1:system.setParticleMass(atom.idx,0)
    integrator=mm.VerletIntegrator(.001)
    context=mm.Context(system,integrator,mm.Platform.getPlatformByName('Reference'))
    context.setPositions(xyz*unit.nanometer)
    mm.LocalEnergyMinimizer.minimize(context,10,500)
    positions=np.array(context.getState(getPositions=True).getPositions(asNumpy=True).value_in_unit(unit.nanometer))
    del context,integrator
    heavy=[a.idx for a in original.atoms if a.atomic_number>1]
    if not np.array_equal(positions[heavy],xyz[heavy]):raise ValueError('GDP hydrogen construction moved observed heavy atoms')
    np.save(out/'bound-positions-nm.npy',positions)
    comparisons=[]
    rng=np.random.default_rng(2026)
    for i in range(3):
        probe=positions if i==0 else positions+rng.normal(0,.0002,positions.shape)
        e1,f1=evaluate(source/'gdp.prmtop',probe);e2,f2=evaluate(out/'gdp.prmtop',probe)
        maximum=float(np.max(abs(f1-f2)));delta=e2-e1
        (out/f'comparison-{i}.json').write_text(json.dumps({'source_energy_kj_mol':e1,'scoped_energy_kj_mol':e2,
            'energy_difference_kj_mol':delta,'maximum_source_force':float(np.max(abs(f1))),
            'maximum_force_difference':maximum},indent=2)+'\n')
        if not np.isfinite(e2) or not np.isfinite(f2).all() or abs(delta)>1e-4 or maximum>1e-3:
            raise ValueError(f'Scoped GDP energy/forces differ: {delta}, {maximum}')
        comparisons.append({'probe':i,'energy_difference_kj_mol':delta,'maximum_force_component_difference_kj_mol_nm':maximum})
    np.save(out/'bound-positions-nm.npy',positions)
    (out/'bound-atoms.json').write_text(json.dumps([{'name':a.name,'source_name':a.name.replace('*',"'"),'element':app.element.Element.getByAtomicNumber(a.atomic_number).symbol,
           'xyz_angstrom':(positions[a.idx]*10).tolist(),'source_heavy':a.atomic_number>1} for a in original.atoms],indent=2)+'\n')
    report={'stage':'scoped_source_GDP_candidate_complete','physical_model_validated':False,'app_ready':False,
        'source_model':'Preserved Meagher–Redman–Carlson GDP(-3) native Amber99-based source model, scoped without protein parameter overrides.',
        'source_parameters_unchanged':True,'original_heavy_coordinates_preserved':True,
        'protein_parameter_namespace_isolated':True,'type_renames':types,'comparisons':comparisons,
        'native_elements_explicitly_preserved':True,
        'atom_count':len(result.atoms),'charge_e':sum(a.charge for a in result.atoms),'commands':commands,
        'source_sha256':{str(p):digest(p) for p in [source/'gdp.prmtop',source/'gdp.inpcrd',root/'gdp-source/GDP.prep']}}
    (out/'result.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps(report,indent=2))


if __name__=='__main__':
    p=argparse.ArgumentParser()
    for n in ['root','amber','output']:p.add_argument('--'+n,type=Path,required=True)
    a=p.parse_args();run(a.root,a.amber,a.output)
