"""Compare declared native site parameters with an actual OpenMM System.

This validates implementation consistency, not chemical-state correctness or
force-field accuracy. Expected values must come from an independent native
parameter artifact; reconstructing them from the System under test is circular.
No bond or coordination partner is inferred from atom proximity.
"""
from __future__ import annotations

from collections import defaultdict
import hashlib
import json
import math
from pathlib import Path

import numpy as np

IDENTITY_FIELDS=('chain','resid','insertion','resname','atomname')


def _canonical_hash(value):
    return hashlib.sha256(json.dumps(value,sort_keys=True,separators=(',', ':'),allow_nan=False).encode()).hexdigest()


def atom_identity(atom):
    residue=atom.residue
    return (str(residue.chain.id),str(residue.id),
            str(getattr(residue,'insertionCode','') or '').strip(),
            str(residue.name),str(atom.name))


def _identity(record):
    if not isinstance(record,dict) or set(record)!=set(IDENTITY_FIELDS) or any(not isinstance(record[k],str) for k in IDENTITY_FIELDS):
        raise ValueError('Every atom identity must explicitly contain string chain,resid,insertion,resname,atomname')
    return tuple(record[k] for k in IDENTITY_FIELDS)


def _number(value,name,minimum=None):
    if isinstance(value,bool) or not isinstance(value,(float,int)) or not math.isfinite(value) or (minimum is not None and value<minimum):
        raise ValueError(f'{name} must be finite'+(f' and >= {minimum}' if minimum is not None else ''))
    return float(value)


def _key(indices):
    indices=tuple(indices)
    return min(indices,indices[::-1])


def _close(a,b,atol,rtol):
    return abs(a-b)<=atol+rtol*max(abs(a),abs(b))


def _phase_close(a,b,atol):
    return abs(math.remainder(a-b,2*math.pi))<=atol


def _manifest(manifest):
    if not isinstance(manifest,dict) or manifest.get('schema_version')!=1:
        raise ValueError('Expected a version-1 parameter manifest')
    if manifest.get('scope') not in ('full_system','site'):
        raise ValueError('scope must explicitly be full_system or site')
    if not isinstance(manifest.get('model_id'),str) or not manifest['model_id']:
        raise ValueError('A model_id is required')
    reference=manifest.get('native_reference')
    if not isinstance(reference,dict) or not isinstance(reference.get('sha256'),str) or len(reference['sha256'])!=64 or any(x not in '0123456789abcdef' for x in reference['sha256']):
        raise ValueError('native_reference must identify the independently generated source with a SHA-256 hash')
    atoms=manifest.get('atoms')
    if not isinstance(atoms,list) or not atoms:
        raise ValueError('A nonempty native atom inventory is required')
    ids=set();identities=set()
    for atom in atoms:
        if not isinstance(atom,dict) or not isinstance(atom.get('id'),str) or not atom['id'] or atom['id'] in ids:
            raise ValueError('Manifest atom IDs must be nonempty and unique')
        identity=_identity(atom.get('identity'))
        if identity in identities:raise ValueError('Manifest contains duplicate full atom identities')
        ids.add(atom['id']);identities.add(identity)
        if not isinstance(atom.get('element'),str) or not atom['element']:
            raise ValueError('Each atom must declare its element')
        for field,minimum in [('mass_da',0),('charge_e',None),('sigma_nm',0),('epsilon_kj_mol',0)]:
            _number(atom.get(field),field,minimum)
    definitions={'bonds':(2,[('length_nm',0),('k_kj_mol_nm2',0)]),
        'angles':(3,[('angle_radian',0),('k_kj_mol_rad2',0)]),
        # Periodic Fourier amplitudes can be signed (e.g. the published Amber
        # polyphosphate terms). Preserve their sign and phase exactly: replacing
        # k with abs(k) changes both the torsion surface and its energy origin.
        'torsions':(4,[('periodicity',0),('phase_radian',None),('k_kj_mol',None)]),
        'exceptions':(2,[('chargeprod_e2',None),('sigma_nm',0),('epsilon_kj_mol',0)]),
        'constraints':(2,[('length_nm',0)])}
    for group,(n,fields) in definitions.items():
        # Earlier unconstrained manifests omit this protocol field. Omission
        # means zero constraints, never permission to infer them from the System.
        terms=manifest.get(group,[]) if group=='constraints' else manifest.get(group)
        if not isinstance(terms,list):raise ValueError(f'{group} must explicitly be a list, including when empty')
        constraint_pairs=set()
        for term in terms:
            if not isinstance(term,dict) or not isinstance(term.get('atoms'),list) or len(term['atoms'])!=n or len(set(term['atoms']))!=n or any(x not in ids for x in term['atoms']):
                raise ValueError(f'Invalid or unmapped atom list in {group}')
            for field,minimum in fields:_number(term.get(field),field,minimum)
            if group=='constraints':
                pair=_key(term['atoms'])
                if term['length_nm']<=0 or pair in constraint_pairs:
                    raise ValueError('Declared constraints require positive lengths and unique atom pairs')
                constraint_pairs.add(pair)
            if group=='bonds' and (term['length_nm']<=0 or term['k_kj_mol_nm2']<=0):
                raise ValueError('A declared physical bond requires positive equilibrium length and force constant')
            if group=='angles' and not 0<term['angle_radian']<=math.pi:
                raise ValueError('Equilibrium angles must be in (0,pi] radians')
            if group=='torsions' and (isinstance(term['periodicity'],bool) or not isinstance(term['periodicity'],int)):
                raise ValueError('Torsion periodicity must be an integer')
            if group=='torsions' and term['periodicity']<1:
                raise ValueError('Torsion periodicity must be positive; even a zero-amplitude periodicity-zero native placeholder prevents an OpenMM Context')
            if group=='bonds' and term.get('topology_role','chemical_bond') not in ('chemical_bond','native_water_hh_geometry'):
                raise ValueError('Unknown bond topology_role')
    return atoms


def validate_bonded_site(topology,system,manifest,*,atol=1e-6,rtol=1e-7):
    """Validate all standard internal terms, masses, charges/LJ and exceptions.

    A site-scope atom inventory defines the exact coverage boundary. All actual
    standard terms entirely within it must match a declared native term. Terms
    touching atoms outside it are reported separately and are not validated.
    Topological bonds and physical HarmonicBondForce terms are both required.
    Constraints and custom restraint forces cannot satisfy that requirement.
    A full-system verdict also requires coverage of every potential-energy
    force. Unvalidated custom or CMAP terms cannot be silently accepted merely
    because the standard terms match. Constraints must match an independently
    declared protocol inventory, including when physical bond terms are present.
    Virtual-site geometry is currently unsupported within the validated scope.
    """
    from openmm import unit,XmlSerializer
    report={'schema_version':1,'accepted':False,'checks':[],'failures':[],
        'meaning_of_accepted':'Native-declared standard mechanical terms match the System within the explicit scope; not model accuracy, oxidation-state validation, or stability proof.',
        'tolerances':{'absolute':atol,'relative':rtol}}
    def fail(kind,**detail):report['failures'].append({'kind':kind,**detail})
    try:
        atoms=_manifest(manifest)
        _number(atol,'atol',0);_number(rtol,'rtol',0)
        report.update(model_id=manifest['model_id'],scope=manifest['scope'],native_reference=manifest['native_reference'],
            chemical_state_metadata=manifest.get('chemical_state'),chemical_state_validated=False)
        native_atoms=list(topology.atoms())
        report.update(validator_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
            system_xml_sha256=hashlib.sha256(XmlSerializer.serialize(system).encode()).hexdigest(),
            manifest_canonical_sha256=_canonical_hash(manifest),
            topology_identity_sha256=_canonical_hash({
                'atom_identities':[list(atom_identity(atom)) for atom in native_atoms],
                'bond_index_pairs':sorted([list(_key([b[0].index,b[1].index])) for b in topology.bonds()])}))
        if system.getNumParticles()!=len(native_atoms):
            fail('topology_system_particle_count_mismatch',topology=len(native_atoms),system=system.getNumParticles())
            return report
        lookup=defaultdict(list)
        for atom in native_atoms:lookup[atom_identity(atom)].append(atom)
        mapping={}
        for expected in atoms:
            matches=lookup[_identity(expected['identity'])]
            if len(matches)!=1:
                fail('atom_identity_not_unique' if matches else 'missing_atom',id=expected['id'],identity=expected['identity'],matches=len(matches))
            else:mapping[expected['id']]=matches[0].index
        if len(mapping)!=len(atoms):return report
        included=set(mapping.values());by_index={mapping[a['id']]:a['id'] for a in atoms}
        if manifest['scope']=='full_system' and len(included)!=len(native_atoms):
            fail('full_system_inventory_incomplete',declared=len(included),actual=len(native_atoms))
        report['atom_index_map']=mapping
        report['coverage']={'declared_atoms':len(included),'system_atoms':len(native_atoms),
            'rule':'Every standard term with all atoms inside the declared inventory; boundary terms are reported without validating them.'}
        virtual_sites=[]
        for index in range(system.getNumParticles()):
            if not system.isVirtualSite(index):continue
            site=system.getVirtualSite(index)
            virtual_sites.append({'particle_index':index,'atom_id':by_index.get(index),
                'identity':dict(zip(IDENTITY_FIELDS,atom_identity(native_atoms[index]))),
                'class':site.__class__.__name__,
                'parent_particle_indices':[site.getParticle(n) for n in range(site.getNumParticles())],
                'inside_declared_inventory':index in included})
        report['virtual_sites_not_validated']=virtual_sites
        scoped_virtual_sites=[site for site in virtual_sites
            if manifest['scope']=='full_system' or site['inside_declared_inventory']]
        if scoped_virtual_sites:
            fail('unvalidated_virtual_sites',sites=scoped_virtual_sites)
        nb=[force for force in system.getForces() if force.__class__.__name__=='NonbondedForce']
        if len(nb)!=1:
            fail('nonbonded_force_count',expected=1,actual=len(nb));return report
        nonbonded=nb[0]
        if nonbonded.getNumParticleParameterOffsets() or nonbonded.getNumExceptionParameterOffsets():
            fail('unsupported_nonbonded_parameter_offsets')
        for expected in atoms:
            index=mapping[expected['id']];actual_atom=native_atoms[index]
            element=actual_atom.element.symbol if actual_atom.element else None
            if element!=expected['element']:fail('element_mismatch',atom=expected['id'],expected=expected['element'],actual=element)
            charge,sigma,epsilon=nonbonded.getParticleParameters(index)
            values={'mass_da':system.getParticleMass(index).value_in_unit(unit.dalton),
                'charge_e':charge.value_in_unit(unit.elementary_charge),
                'sigma_nm':sigma.value_in_unit(unit.nanometer),
                'epsilon_kj_mol':epsilon.value_in_unit(unit.kilojoule_per_mole)}
            for key,actual in values.items():
                if key=='sigma_nm' and expected['epsilon_kj_mol']==0 and values['epsilon_kj_mol']==0:
                    continue  # Zero-epsilon particles have no LJ sigma-dependent interaction.
                if not math.isfinite(actual) or not _close(expected[key],actual,atol,rtol):
                    fail('atom_parameter_mismatch',atom=expected['id'],parameter=key,expected=expected[key],actual=actual)
        actual={name:defaultdict(list) for name in ('bonds','angles','torsions','exceptions','constraints')}
        boundary=[];additional=[]
        def add(group,indices,values,origin):
            indices=tuple(map(int,indices))
            if set(indices)<=included:actual[group][_key(indices)].append({'values':tuple(map(float,values)),'origin':origin})
            elif set(indices)&included:boundary.append({'group':group,'particle_indices':list(indices),'origin':origin})
        for force_index,force in enumerate(system.getForces()):
            kind=force.__class__.__name__;origin={'force_index':force_index,'class':kind,'name':force.getName()}
            if kind=='HarmonicBondForce':
                for n in range(force.getNumBonds()):
                    i,j,length,k=force.getBondParameters(n)
                    add('bonds',[i,j],[length.value_in_unit(unit.nanometer),k.value_in_unit(unit.kilojoule_per_mole/unit.nanometer**2)],dict(origin,term_index=n))
            elif kind=='HarmonicAngleForce':
                for n in range(force.getNumAngles()):
                    i,j,k,angle,constant=force.getAngleParameters(n)
                    add('angles',[i,j,k],[angle.value_in_unit(unit.radian),constant.value_in_unit(unit.kilojoule_per_mole/unit.radian**2)],dict(origin,term_index=n))
            elif kind=='PeriodicTorsionForce':
                for n in range(force.getNumTorsions()):
                    i,j,k,l,period,phase,constant=force.getTorsionParameters(n)
                    add('torsions',[i,j,k,l],[period,phase.value_in_unit(unit.radian),constant.value_in_unit(unit.kilojoule_per_mole)],dict(origin,term_index=n))
            elif kind=='NonbondedForce':
                for n in range(force.getNumExceptions()):
                    i,j,chargeprod,sigma,epsilon=force.getExceptionParameters(n)
                    add('exceptions',[i,j],[chargeprod.value_in_unit(unit.elementary_charge**2),sigma.value_in_unit(unit.nanometer),epsilon.value_in_unit(unit.kilojoule_per_mole)],dict(origin,term_index=n))
            elif kind!='CMMotionRemover':additional.append(origin)
        report['additional_forces_not_used_as_parameter_evidence']=additional
        report['boundary_terms_not_validated']=boundary
        report['all_potential_energy_forces_covered']=not additional
        if additional and manifest['scope']=='full_system':
            fail('unvalidated_additional_forces',forces=additional)
        constraints=[]
        for n in range(system.getNumConstraints()):
            first,second,length=system.getConstraintParameters(n)
            length_nm=length.value_in_unit(unit.nanometer)
            add('constraints',[first,second],[length_nm],{'class':'System','constraint_index':n})
            if {int(first),int(second)}&included:
                constraints.append({'particle_indices':[int(first),int(second)],'length_nm':length_nm})
        report['constraints_not_used_as_bond_evidence']=constraints
        report['constraint_inventory_policy']='Every constraint within the declared atom inventory must match the independent protocol; an omitted manifest inventory means no constraints.'
        topology_bonds={_key([b[0].index,b[1].index]) for b in topology.bonds()}
        expected_bonds={_key([mapping[x] for x in term['atoms']]) for term in manifest['bonds']}
        geometry_only=set();report['native_water_geometry_terms_without_topology_edges']=[]
        expected_by_id={atom['id']:atom for atom in atoms}
        neighbors=defaultdict(set)
        for first,second in topology_bonds:
            neighbors[first].add(second);neighbors[second].add(first)
        for term in manifest['bonds']:
            if term.get('topology_role')!='native_water_hh_geometry':continue
            first,second=[expected_by_id[x] for x in term['atoms']]
            pair=_key([mapping[x] for x in term['atoms']])
            shared_oxygen=[index for index in neighbors[pair[0]]&neighbors[pair[1]]
                if native_atoms[index].element and native_atoms[index].element.symbol=='O' and
                atom_identity(native_atoms[index])[:4]==_identity(first['identity'])[:4]]
            if (first['element']!='H' or second['element']!='H' or
                first.get('native_residue_name')!='WAT' or second.get('native_residue_name')!='WAT' or
                _identity(first['identity'])[:4]!=_identity(second['identity'])[:4] or len(shared_oxygen)!=1):
                fail('invalid_native_water_geometry_classification',atoms=term['atoms'])
            else:
                geometry_only.add(pair)
                if pair not in topology_bonds:
                    report['native_water_geometry_terms_without_topology_edges'].append(term['atoms'])
        for pair in expected_bonds-topology_bonds-geometry_only:
            fail('missing_topological_bond',atoms=[by_index[x] for x in pair])
        for pair in {p for p in topology_bonds if set(p)<=included}-expected_bonds:
            fail('undeclared_topological_bond',atoms=[by_index[x] for x in pair])
        names={'bonds':['length_nm','k_kj_mol_nm2'],'angles':['angle_radian','k_kj_mol_rad2'],
            'torsions':['periodicity','phase_radian','k_kj_mol'],
            'exceptions':['chargeprod_e2','sigma_nm','epsilon_kj_mol'],'constraints':['length_nm']}
        def matches(group,expected,observed):
            for index,(a,b) in enumerate(zip(expected,observed)):
                if not math.isfinite(b):return False
                if group=='exceptions' and index==1 and expected[2]==0 and observed[2]==0:
                    continue  # Sigma has no physical effect when epsilon is zero.
                if group=='torsions' and index==1:
                    if not _phase_close(a,b,atol):return False
                elif not _close(a,b,atol,rtol):return False
            return True
        report['zero_amplitude_torsions_omitted_by_reader']=[]
        for group,fields in names.items():
            for term in manifest.get(group,[]):
                indices=_key([mapping[x] for x in term['atoms']]);values=tuple(term[x] for x in fields)
                choices=actual[group].get(indices,[])
                match=next((j for j,x in enumerate(choices) if matches(group,values,x['values'])),None)
                if match is None:
                    if group=='torsions' and term['k_kj_mol']==0:
                        # GromacsTopFile omits exact-zero torsions. This is a
                        # null potential, not permission to alter nonzero terms
                        # or infer/exclude the independent native 1-4 graph.
                        report['zero_amplitude_torsions_omitted_by_reader'].append({
                            'atoms':term['atoms'],'periodicity':term['periodicity'],
                            'phase_radian':term['phase_radian'],'k_kj_mol':0.0})
                        continue
                    fail('missing_or_mismatched_constraint' if group=='constraints' else 'missing_or_mismatched_force_term',group=group,atoms=term['atoms'],expected=dict(zip(fields,values)),
                        actual_candidates=[dict(zip(fields,x['values'])) for x in choices])
                else:choices.pop(match)
            for indices,extra in actual[group].items():
                for term in extra:
                    fail('undeclared_constraint' if group=='constraints' else 'undeclared_force_term',group=group,atoms=[by_index[x] for x in indices],
                        values=dict(zip(fields,term['values'])),origin=term['origin'])
        report['checks']=['unique_full_atom_identity','particle_inventory','element_mass_charge_lj',
            'topological_bonds','harmonic_bond_terms','harmonic_angle_terms','ordered_periodic_torsion_terms',
            'nonbonded_exceptions','no_undeclared_internal_standard_terms','declared_constraint_inventory',
            'no_unvalidated_in_scope_virtual_sites']
        if manifest['scope']=='full_system':report['checks'].append('no_unvalidated_potential_energy_forces')
        report['accepted']=not report['failures']
    except (ValueError,TypeError,KeyError,AttributeError) as exc:
        fail('invalid_manifest_or_unsupported_system',message=str(exc))
    return report


def check_finite_forces(system,positions_nm,*,platform_name='Reference'):
    """One explicit configuration; finite energy/forces are necessary, not stability."""
    import openmm as mm
    from openmm import unit
    positions=np.asarray(positions_nm,dtype=float)
    if positions.shape!=(system.getNumParticles(),3) or not np.isfinite(positions).all():
        return {'accepted':False,'reason':'Positions must be a finite particle-count×3 array in nm'}
    context=None
    integrator=mm.VerletIntegrator(0.001*unit.picosecond)
    try:
        context=mm.Context(system,integrator,mm.Platform.getPlatformByName(platform_name))
        context.setPositions(positions*unit.nanometer)
        state=context.getState(getEnergy=True,getForces=True)
        energy=float(state.getPotentialEnergy().value_in_unit(unit.kilojoule_per_mole))
        forces=np.asarray(state.getForces(asNumpy=True).value_in_unit(unit.kilojoule_per_mole/unit.nanometer))
        return {'accepted':bool(math.isfinite(energy) and np.isfinite(forces).all()),
            'potential_energy_kj_mol':energy,'max_force_kj_mol_nm':float(np.max(np.linalg.norm(forces,axis=1))),
            'platform':platform_name,'meaning':'Finite single-configuration energy/forces; not stability or model-accuracy evidence.'}
    except Exception as exc:
        return {'accepted':False,'reason':f'{type(exc).__name__}: {exc}'}
    finally:
        del context
        del integrator


def _minimum_image(displacements,box):
    """Nearest images for the reduced triclinic convention used by OpenMM.

    Reject other cell conventions rather than silently returning a wrong image.
    The nearby lattice-image search also handles skewed reduced cells where
    component-wise fractional rounding alone is insufficient.
    """
    if box is None:return np.asarray(displacements,dtype=float)
    box=np.asarray(box,dtype=float)
    if box.shape!=(3,3) or not np.isfinite(box).all():
        raise ValueError('Periodic boxes must be finite 3×3 row-vector arrays in nm')
    tolerance=1e-9
    if (np.any(np.diag(box)<=0) or np.max(np.abs(np.triu(box,1)))>tolerance or
        abs(box[1,0])>box[0,0]/2+tolerance or abs(box[2,0])>box[0,0]/2+tolerance or
        abs(box[2,1])>box[1,1]/2+tolerance):
        raise ValueError('Periodic box must use a reduced OpenMM/GROMACS triclinic row-vector convention')
    values=np.asarray(displacements,dtype=float)
    fractional=values @ np.linalg.inv(box)
    fractional-=np.rint(fractional)
    shifts=np.array([(i,j,k) for i in (-1,0,1) for j in (-1,0,1) for k in (-1,0,1)])
    candidates=(fractional[:,None,:]+shifts[None,:,:]) @ box
    best=np.argmin(np.einsum('ijk,ijk->ij',candidates,candidates),axis=1)
    return candidates[np.arange(len(values)),best]


def audit_site_trajectory(frame_atom_ids,frames_nm,manifest,*,parameter_report,
                         bond_bounds,motion_atom_ids,boxes_nm=None,
                         forces_kj_mol_nm=None,minimum_internal_motion_nm=1e-4,
                         coordination_sites=None):
    """Assess explicit bond/coordination bounds, internal motion and finite forces.

    The caller supplies the exact saved trajectory atom-ID order, native model
    validation report, physically chosen bounds, and motion-selection IDs.
    Bounds are operational checks, not parameters inferred from distances.
    Forces must refer to the same saved frame/atom order. Missing force evidence
    produces an incomplete result rather than a claimed stability pass.
    """
    report={'accepted':False,'failures':[],
        'meaning_of_accepted':'The supplied frames satisfy the declared mechanical/geometry/motion checks; not force-field accuracy, chemical-state validation, or production stability.'}
    def fail(kind,**detail):report['failures'].append({'kind':kind,**detail})
    try:
        atoms=_manifest(manifest)
        ids={a['id'] for a in atoms}
        if (not isinstance(parameter_report,dict) or not parameter_report.get('accepted') or
            parameter_report.get('model_id')!=manifest['model_id'] or
            parameter_report.get('manifest_canonical_sha256')!=_canonical_hash(manifest) or
            parameter_report.get('native_reference',{}).get('sha256')!=manifest['native_reference']['sha256'] or
            set(parameter_report.get('atom_index_map',{}))!=ids):
            fail('missing_or_inconsistent_parameter_validation');return report
        if not isinstance(frame_atom_ids,(list,tuple)) or any(not isinstance(x,str) for x in frame_atom_ids) or len(set(frame_atom_ids))!=len(frame_atom_ids):
            raise ValueError('Saved frame atom IDs must be explicit and unique')
        mapping={value:index for index,value in enumerate(frame_atom_ids)}
        if ids-set(mapping):
            fail('trajectory_atom_map_missing_ids',ids=sorted(ids-set(mapping)));return report
        frames=np.asarray(frames_nm,dtype=float)
        if frames.ndim!=3 or frames.shape[1:]!=(len(frame_atom_ids),3) or len(frames)<2 or not np.isfinite(frames).all():
            raise ValueError('Trajectory requires at least two finite frames with the exact saved atom count and nm units')
        boxes=np.asarray(boxes_nm,dtype=float) if boxes_nm is not None else None
        if boxes is not None and boxes.shape!=(len(frames),3,3):
            raise ValueError('boxes_nm requires one explicit 3×3 row-vector box per frame')
        if not isinstance(bond_bounds,list) or not bond_bounds:
            raise ValueError('At least one explicitly chosen bond bound is required')
        native_bonds={_key(term['atoms']) for term in manifest['bonds']}
        by_pair={};bound_indices=[]
        for entry in bond_bounds:
            if not isinstance(entry,dict) or not isinstance(entry.get('atoms'),list) or len(entry['atoms'])!=2:
                raise ValueError('Malformed bond bound')
            pair=_key(entry['atoms'])
            if pair not in native_bonds or pair in by_pair:
                raise ValueError('Every monitored pair must be a unique explicitly parameterized native bond')
            lower=_number(entry.get('minimum_nm'),'minimum_nm',0)
            upper=_number(entry.get('maximum_nm'),'maximum_nm',0)
            if not 0<lower<upper:raise ValueError('Bond bounds must satisfy 0 < minimum_nm < maximum_nm')
            by_pair[pair]=len(bound_indices)
            bound_indices.append([mapping[x] for x in entry['atoms']])
        index=np.asarray(bound_indices,dtype=int)
        lengths=[]
        for frame_number,coordinates in enumerate(frames):
            delta=coordinates[index[:,1]]-coordinates[index[:,0]]
            delta=_minimum_image(delta,boxes[frame_number] if boxes is not None else None)
            lengths.append(np.linalg.norm(delta,axis=1))
        lengths=np.asarray(lengths)
        valid=np.empty(lengths.shape,dtype=bool)
        report['bonds']=[]
        for column,entry in enumerate(bond_bounds):
            series=lengths[:,column]
            valid[:,column]=(series>=entry['minimum_nm']) & (series<=entry['maximum_nm'])
            bad=np.flatnonzero(~valid[:,column])
            report['bonds'].append({'atoms':entry['atoms'],'minimum_observed_nm':float(series.min()),
                'maximum_observed_nm':float(series.max()),'mean_nm':float(series.mean()),
                'declared_minimum_nm':entry['minimum_nm'],'declared_maximum_nm':entry['maximum_nm'],
                'violating_frame_count':len(bad),'first_violating_frames':bad[:20].tolist()})
            if len(bad):fail('bond_distance_outside_declared_bounds',atoms=entry['atoms'],frames=bad[:20].tolist())
        report['coordination_sites']=[]
        for site in coordination_sites or []:
            metal=site['metal_atom_id'];donors=site['donor_atom_ids']
            if metal not in ids or not isinstance(donors,list) or not donors or len(set(donors))!=len(donors) or any(x not in ids for x in donors):
                raise ValueError('A coordination site needs explicit unique mapped metal/donor atom IDs')
            columns=[by_pair[_key([metal,donor])] for donor in donors]
            counts=valid[:,columns].sum(axis=1)
            report['coordination_sites'].append({'metal_atom_id':metal,'declared_donor_atom_ids':donors,
                'declared_donor_count':len(donors),'minimum_donors_within_declared_bounds':int(counts.min()),
                'maximum_donors_within_declared_bounds':int(counts.max()),
                'basis':'Only the supplied parameterized donor bonds; no new chemical bond assignment from proximity.'})
        if not isinstance(motion_atom_ids,list) or not 2<=len(motion_atom_ids)<=64 or len(set(motion_atom_ids))!=len(motion_atom_ids) or any(x not in ids for x in motion_atom_ids):
            raise ValueError('Choose 2–64 unique native atom IDs for the internal-motion diagnostic')
        threshold=_number(minimum_internal_motion_nm,'minimum_internal_motion_nm',0)
        motion_indices=[mapping[x] for x in motion_atom_ids]
        first,second=np.triu_indices(len(motion_indices),1)
        distances=[]
        for frame_number,coordinates in enumerate(frames):
            selected=coordinates[motion_indices]
            deltas=_minimum_image(selected[second]-selected[first],boxes[frame_number] if boxes is not None else None)
            distances.append(np.linalg.norm(deltas,axis=1))
        distances=np.asarray(distances)
        internal=np.sqrt(np.mean((distances-distances[0])**2,axis=1))
        report['motion']={'atom_ids':motion_atom_ids,'metric':'RMS change of selected pair distances relative to the first frame; rigid-body motion alone does not count.',
            'maximum_internal_motion_nm':float(internal.max()),'required_minimum_internal_motion_nm':threshold}
        if internal.max()<=threshold:fail('no_meaningful_internal_motion')
        if forces_kj_mol_nm is None:
            report['forces']={'checked':False};fail('force_evidence_not_supplied')
        else:
            forces=np.asarray(forces_kj_mol_nm,dtype=float)
            if forces.shape!=frames.shape or not np.isfinite(forces).all():
                report['forces']={'checked':True,'finite':False};fail('invalid_or_nonfinite_frame_forces')
            else:report['forces']={'checked':True,'finite':True,'maximum_norm_kj_mol_nm':float(np.linalg.norm(forces,axis=2).max())}
        report.update(frame_count=len(frames),saved_frame_atom_count=len(frame_atom_ids),
            monitored_native_bonds=len(bond_bounds),total_declared_native_bonds=len(manifest['bonds']),
            periodic_boxes_supplied=boxes is not None,model_id=manifest['model_id'],
            native_reference=manifest['native_reference'])
        report['parameter_validation_bindings']={key:parameter_report[key] for key in
            ('validator_sha256','system_xml_sha256','manifest_canonical_sha256','topology_identity_sha256')}
        report['constraints_not_accuracy_evidence']=parameter_report.get('constraints_not_used_as_bond_evidence',[])
        report['extra_forces_not_accuracy_evidence']=parameter_report.get('additional_forces_not_used_as_parameter_evidence',[])
        report['accepted']=not report['failures']
    except (ValueError,TypeError,KeyError,IndexError,np.linalg.LinAlgError) as exc:
        fail('invalid_trajectory_or_monitor_contract',message=str(exc))
    return report
