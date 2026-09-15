import copy
import hashlib
import json
import math

import numpy as np
import openmm as mm
import pytest
from openmm import app,unit

from backend.bonded_site_validation import validate_bonded_site,check_finite_forces,audit_site_trajectory,_minimum_image


def specimen():
    """Independent hand-declared four-atom mechanics control, not a fitted model."""
    topology=app.Topology();chain=topology.addChain('A')
    residue=topology.addResidue('LIG',chain,'42',insertionCode='B')
    elements=[app.element.carbon,app.element.sulfur,app.element.carbon,app.element.carbon]
    atoms=[topology.addAtom(name,element,residue) for name,element in zip(['C1','S1','C2','C3'],elements)]
    system=mm.System();nb=mm.NonbondedForce()
    masses=[12.01,32.06,12.01,12.01]
    charges=[-0.2,-0.1,0.2,0.1]
    for mass,q in zip(masses,charges):
        system.addParticle(mass);nb.addParticle(q,0.34,0.4)
    bonds=mm.HarmonicBondForce();angles=mm.HarmonicAngleForce();torsions=mm.PeriodicTorsionForce()
    for i,j in [(0,1),(1,2),(2,3)]:
        topology.addBond(atoms[i],atoms[j]);bonds.addBond(i,j,0.15,83680)
    for i,j,k in [(0,1,2),(1,2,3)]:angles.addAngle(i,j,k,math.radians(109.5),418.4)
    torsions.addTorsion(0,1,2,3,3,0.4,1.2552)
    for i,j in [(0,1),(1,2),(2,3),(0,2),(1,3)]:nb.addException(i,j,0,1,0)
    nb.addException(0,3,-0.2*0.1/1.2,0.34,0.2)
    for force in [bonds,angles,torsions,nb]:system.addForce(force)
    # These source-side numbers use explicit Amber-to-OpenMM conversions, not
    # inspection of the System under test.
    ids=['c1','s1','c2','c3']
    manifest={'schema_version':1,'scope':'full_system','model_id':'mechanics-control',
        'native_reference':{'sha256':'0'*64,'description':'Hand-declared test parameters'},
        'chemical_state':{'oxidation_state':None},
        'atoms':[{'id':ids[i],'identity':{'chain':'A','resid':'42','insertion':'B','resname':'LIG','atomname':atoms[i].name},
                  'element':elements[i].symbol,'mass_da':masses[i],'charge_e':charges[i],'sigma_nm':0.34,'epsilon_kj_mol':0.4} for i in range(4)],
        'bonds':[{'atoms':[ids[i],ids[j]],'length_nm':1.5*0.1,'k_kj_mol_nm2':100*2*4.184*100} for i,j in [(0,1),(1,2),(2,3)]],
        'angles':[{'atoms':[ids[i],ids[j],ids[k]],'angle_radian':math.radians(109.5),'k_kj_mol_rad2':50*2*4.184} for i,j,k in [(0,1,2),(1,2,3)]],
        'torsions':[{'atoms':ids,'periodicity':3,'phase_radian':0.4,'k_kj_mol':0.3*4.184}],
        'exceptions':[{'atoms':[ids[i],ids[j]],'chargeprod_e2':0,'sigma_nm':1,'epsilon_kj_mol':0} for i,j in [(0,1),(1,2),(2,3),(0,2),(1,3)]]+
            [{'atoms':['c1','c3'],'chargeprod_e2':-0.2*0.1/1.2,'sigma_nm':0.34,'epsilon_kj_mol':0.4/2}]}
    return topology,system,manifest


def test_full_native_term_inventory_matches_independent_atom_order():
    topology,system,manifest=specimen()
    manifest['atoms'].reverse()
    report=validate_bonded_site(topology,system,manifest)
    assert report['accepted'],report['failures']
    assert report['atom_index_map']=={'c3':3,'c2':2,'s1':1,'c1':0}
    assert report['chemical_state_validated'] is False


def test_bond_half_k_conversion_error_is_detected():
    topology,system,manifest=specimen()
    system.getForce(0).setBondParameters(0,0,1,0.15,41840)
    report=validate_bonded_site(topology,system,manifest)
    assert not report['accepted']
    assert any(x['kind']=='missing_or_mismatched_force_term' and x['group']=='bonds' for x in report['failures'])


def test_signed_native_fourier_amplitude_preserves_energy_and_parameter_identity():
    topology,system,manifest=specimen()
    # A negative Fourier amplitude is bounded, unlike a negative harmonic bond
    # stiffness. This is an independent analytic torsion control, not a fit.
    torsions=system.getForce(2)
    torsions.setTorsionParameters(0,0,1,2,3,3,0.,-1.2552)
    torsions.setForceGroup(1)
    manifest['torsions'][0].update(k_kj_mol=-0.3*4.184,phase_radian=0.)
    report=validate_bonded_site(topology,system,manifest)
    assert report['accepted'],report['failures']
    phi=.7
    xyz=np.array([[0,.1,0],[0,0,0],[.1,0,0],[.1,.1*math.cos(phi),.1*math.sin(phi)]])
    integrator=mm.VerletIntegrator(.001)
    context=mm.Context(system,integrator,mm.Platform.getPlatformByName('Reference'))
    try:
        context.setPositions(xyz*unit.nanometer)
        state=context.getState(getEnergy=True,getForces=True,groups={1})
        expected=-.3*4.184*(1+math.cos(3*phi))
        assert state.getPotentialEnergy().value_in_unit(unit.kilojoule_per_mole)==pytest.approx(expected,abs=1e-12)
        assert np.isfinite(state.getForces(asNumpy=True).value_in_unit(unit.kilojoule_per_mole/unit.nanometer)).all()
    finally:
        del context,integrator
    # A sign flip must still fail the independent native manifest comparison.
    torsions.setTorsionParameters(0,0,1,2,3,3,0.,1.2552)
    report=validate_bonded_site(topology,system,manifest)
    assert not report['accepted']
    assert any(x['kind']=='missing_or_mismatched_force_term' and x['group']=='torsions' for x in report['failures'])


@pytest.mark.parametrize('coefficient',[float('nan'),float('inf'),float('-inf'),True])
def test_signed_torsion_support_rejects_nonfinite_or_boolean_amplitudes(coefficient):
    topology,system,manifest=specimen()
    manifest['torsions'][0]['k_kj_mol']=coefficient
    assert not validate_bonded_site(topology,system,manifest)['accepted']


@pytest.mark.parametrize('group,field',[('bonds','k_kj_mol_nm2'),('angles','k_kj_mol_rad2')])
def test_signed_torsion_support_does_not_allow_negative_harmonic_stiffness(group,field):
    topology,system,manifest=specimen()
    manifest[group][0][field]=-1.
    assert not validate_bonded_site(topology,system,manifest)['accepted']


def test_constraint_and_custom_restraint_cannot_substitute_for_physical_bond():
    topology,system,manifest=specimen()
    system.removeForce(0)
    for i,j in [(0,1),(1,2),(2,3)]:system.addConstraint(i,j,0.15)
    restraint=mm.CustomBondForce('100000*(r-0.15)^2');restraint.setName('artificial site restraint')
    for i,j in [(0,1),(1,2),(2,3)]:restraint.addBond(i,j,[])
    system.addForce(restraint)
    report=validate_bonded_site(topology,system,manifest)
    assert not report['accepted']
    assert len(report['constraints_not_used_as_bond_evidence'])==3
    assert report['additional_forces_not_used_as_parameter_evidence'][0]['class']=='CustomBondForce'


def test_physical_bond_plus_unvalidated_restraint_is_not_full_parameter_coverage():
    topology,system,manifest=specimen()
    restraint=mm.CustomBondForce('100000*(r-0.15)^2')
    restraint.addBond(0,1,[])
    system.addForce(restraint)
    report=validate_bonded_site(topology,system,manifest)
    assert not report['accepted']
    assert report['all_potential_energy_forces_covered'] is False
    assert [failure['kind'] for failure in report['failures']]==['unvalidated_additional_forces']


def test_physical_bond_does_not_authorize_an_undeclared_constraint():
    topology,system,manifest=specimen()
    system.addConstraint(0,1,0.15)
    report=validate_bonded_site(topology,system,manifest)
    assert not report['accepted']
    assert [failure['kind'] for failure in report['failures']]==['undeclared_constraint']


def test_independently_declared_constraint_matches_by_identity_and_length():
    topology,system,manifest=specimen()
    system.addConstraint(0,1,0.15)
    manifest['constraints']=[{'atoms':['s1','c1'],'length_nm':0.15}]
    report=validate_bonded_site(topology,system,manifest)
    assert report['accepted'],report['failures']
    assert report['constraints_not_used_as_bond_evidence']==[
        {'particle_indices':[0,1],'length_nm':0.15}]
    assert 'declared_constraint_inventory' in report['checks']


@pytest.mark.parametrize('defect',['missing','length','pair','duplicate_actual'])
def test_changed_constraint_protocol_is_rejected(defect):
    topology,system,manifest=specimen()
    manifest['constraints']=[{'atoms':['c1','s1'],'length_nm':0.15}]
    if defect!='missing':
        system.addConstraint(0,2 if defect=='pair' else 1,0.16 if defect=='length' else 0.15)
        if defect=='duplicate_actual':system.addConstraint(0,1,0.15)
    report=validate_bonded_site(topology,system,manifest)
    assert not report['accepted']
    assert any(failure['kind'] in ('missing_or_mismatched_constraint','undeclared_constraint')
               for failure in report['failures'])


@pytest.mark.parametrize('defect',['duplicate_expected','zero','nonfinite','unknown_atom','not_list'])
def test_invalid_constraint_declarations_are_rejected(defect):
    topology,system,manifest=specimen()
    manifest['constraints']=[{'atoms':['c1','s1'],'length_nm':0.15}]
    if defect=='duplicate_expected':
        manifest['constraints'].append({'atoms':['s1','c1'],'length_nm':0.15})
    elif defect=='zero':manifest['constraints'][0]['length_nm']=0
    elif defect=='nonfinite':manifest['constraints'][0]['length_nm']=float('nan')
    elif defect=='unknown_atom':manifest['constraints'][0]['atoms'][1]='wrong'
    else:manifest['constraints']=None
    report=validate_bonded_site(topology,system,manifest)
    assert not report['accepted']
    assert report['failures'][0]['kind']=='invalid_manifest_or_unsupported_system'


def test_boundary_constraint_is_explicitly_unvalidated_in_partial_scope():
    topology,system,manifest=specimen()
    manifest['scope']='site';manifest['atoms']=manifest['atoms'][:2]
    manifest['bonds']=manifest['bonds'][:1];manifest['angles']=[];manifest['torsions']=[]
    manifest['exceptions']=manifest['exceptions'][:1]
    system.addConstraint(1,2,0.15)
    report=validate_bonded_site(topology,system,manifest)
    assert report['accepted'],report['failures']
    assert any(term['group']=='constraints' for term in report['boundary_terms_not_validated'])


def test_unvalidated_cmap_is_reported_even_when_standard_terms_match():
    topology,system,manifest=specimen()
    cmap=mm.CMAPTorsionForce()
    cmap.addMap(4,[float(i) for i in range(16)])
    cmap.addTorsion(0,0,1,2,3,0,1,2,3)
    system.addForce(cmap)
    report=validate_bonded_site(topology,system,manifest)
    assert not report['accepted']
    assert report['additional_forces_not_used_as_parameter_evidence'][0]['class']=='CMAPTorsionForce'
    manifest['scope']='site'
    partial=validate_bonded_site(topology,system,manifest)
    assert partial['accepted']
    assert partial['all_potential_energy_forces_covered'] is False
    assert partial['scope']=='site'


@pytest.mark.parametrize('scope',['full_system','site'])
def test_undeclared_virtual_site_is_rejected_despite_matching_mass_and_parameters(scope):
    topology,system,manifest=specimen()
    manifest['scope']=scope
    system.setParticleMass(3,0)
    manifest['atoms'][3]['mass_da']=0
    assert validate_bonded_site(topology,system,manifest)['accepted']
    system.setVirtualSite(3,mm.TwoParticleAverageSite(0,2,0.4,0.6))
    report=validate_bonded_site(topology,system,manifest)
    assert not report['accepted']
    assert [failure['kind'] for failure in report['failures']]==['unvalidated_virtual_sites']
    assert report['virtual_sites_not_validated'][0]['atom_id']=='c3'
    assert report['virtual_sites_not_validated'][0]['parent_particle_indices']==[0,2]
    assert report['virtual_sites_not_validated'][0]['class']=='TwoParticleAverageSite'


def test_outside_virtual_site_remains_visible_as_unvalidated_partial_scope_context():
    topology,system,manifest=specimen()
    manifest['scope']='site';manifest['atoms']=manifest['atoms'][:2]
    manifest['bonds']=manifest['bonds'][:1];manifest['angles']=[];manifest['torsions']=[]
    manifest['exceptions']=manifest['exceptions'][:1]
    system.setParticleMass(3,0)
    system.setVirtualSite(3,mm.TwoParticleAverageSite(0,2,0.4,0.6))
    report=validate_bonded_site(topology,system,manifest)
    assert report['accepted'],report['failures']
    assert report['virtual_sites_not_validated'][0]['inside_declared_inventory'] is False
    assert report['virtual_sites_not_validated'][0]['identity']['atomname']=='C3'


def test_changed_one_four_coulomb_or_lj_scaling_is_detected():
    for change in ['coulomb','lj']:
        topology,system,manifest=specimen();nb=system.getForce(3)
        i,j,q,sigma,epsilon=nb.getExceptionParameters(5)
        nb.setExceptionParameters(5,i,j,q*2 if change=='coulomb' else q,sigma,epsilon*2 if change=='lj' else epsilon)
        report=validate_bonded_site(topology,system,manifest)
        assert any(x.get('group')=='exceptions' for x in report['failures'])
        assert not report['accepted']


def test_missing_and_duplicate_torsions_are_detected():
    topology,system,manifest=specimen()
    system.getForce(2).addTorsion(0,1,2,3,3,0.4,1.2552)
    report=validate_bonded_site(topology,system,manifest)
    assert any(x['kind']=='undeclared_force_term' and x['group']=='torsions' for x in report['failures'])
    topology,system,manifest=specimen();system.removeForce(2)
    report=validate_bonded_site(topology,system,manifest)
    assert any(x['kind']=='missing_or_mismatched_force_term' and x['group']=='torsions' for x in report['failures'])


def test_ordered_torsion_and_equivalent_phase_are_distinguished():
    topology,system,manifest=specimen()
    manifest['torsions'][0]['phase_radian']+=2*math.pi
    assert validate_bonded_site(topology,system,manifest)['accepted']
    manifest['torsions'][0]['atoms']=['c1','c2','s1','c3']
    assert not validate_bonded_site(topology,system,manifest)['accepted']


def test_zero_amplitude_native_periodicity_zero_placeholder_is_unusable():
    topology,system,manifest=specimen()
    system.getForce(2).addTorsion(0,1,2,3,0,0,0)
    manifest['torsions'].append({'atoms':['c1','s1','c2','c3'],'periodicity':0,'phase_radian':0,'k_kj_mol':0})
    result=validate_bonded_site(topology,system,manifest)
    assert not result['accepted']
    assert 'prevents an OpenMM Context' in result['failures'][0]['message']


def test_reader_can_omit_exact_zero_torsion_without_relaxing_nonzero_terms_or_exceptions():
    topology,system,manifest=specimen()
    manifest['torsions'].append({'atoms':['c1','s1','c2','c3'],'periodicity':2,'phase_radian':0,'k_kj_mol':0})
    report=validate_bonded_site(topology,system,manifest)
    assert report['accepted'],report['failures']
    assert len(report['zero_amplitude_torsions_omitted_by_reader'])==1
    manifest['torsions'][-1]['k_kj_mol']=1e-9
    assert not validate_bonded_site(topology,system,manifest)['accepted']
    manifest['torsions'][-1]['k_kj_mol']=0
    nb=system.getForce(3);i,j,q,sigma,epsilon=nb.getExceptionParameters(5)
    nb.setExceptionParameters(5,i,j,0,sigma,0)
    report=validate_bonded_site(topology,system,manifest)
    assert not report['accepted']
    assert any(x.get('group')=='exceptions' for x in report['failures'])


def test_insertion_codes_and_duplicate_atom_identity_are_not_guessed():
    topology,system,manifest=specimen()
    manifest['atoms'][0]['identity']['insertion']='A'
    report=validate_bonded_site(topology,system,manifest)
    assert not report['accepted'] and report['failures'][0]['kind']=='missing_atom'
    topology,system,manifest=specimen()
    residue=list(topology.residues())[0]
    topology.addAtom('C1',app.element.carbon,residue);system.addParticle(12)
    system.getForce(3).addParticle(0,0.34,0.4)
    report=validate_bonded_site(topology,system,manifest)
    assert not report['accepted'] and report['failures'][0]['kind']=='atom_identity_not_unique'


def test_charge_mass_and_element_mismatches_are_independent_checks():
    topology,system,manifest=specimen()
    manifest['atoms'][1].update(charge_e=2,mass_da=65.38,element='Zn')
    report=validate_bonded_site(topology,system,manifest)
    assert not report['accepted']
    assert any(x['kind']=='element_mismatch' for x in report['failures'])
    assert {'mass_da','charge_e'} <= {x.get('parameter') for x in report['failures']}


def test_site_scope_reports_unvalidated_boundary_instead_of_claiming_full_coverage():
    topology,system,manifest=specimen()
    manifest['scope']='site';manifest['atoms']=manifest['atoms'][:2]
    manifest['bonds']=manifest['bonds'][:1];manifest['angles']=[];manifest['torsions']=[];manifest['exceptions']=manifest['exceptions'][:1]
    report=validate_bonded_site(topology,system,manifest)
    assert report['accepted'],report['failures']
    assert report['coverage']['declared_atoms']==2 and report['coverage']['system_atoms']==4
    assert report['boundary_terms_not_validated']
    manifest['scope']='full_system'
    assert not validate_bonded_site(topology,system,manifest)['accepted']


def test_undeclared_topological_link_and_parameter_offset_are_rejected():
    topology,system,manifest=specimen()
    atoms=list(topology.atoms());topology.addBond(atoms[0],atoms[3])
    nb=system.getForce(3);nb.addGlobalParameter('scale',1)
    nb.addParticleParameterOffset('scale',0,1,0,0)
    report=validate_bonded_site(topology,system,manifest)
    assert not report['accepted']
    assert {'undeclared_topological_bond','unsupported_nonbonded_parameter_offsets'} <= {x['kind'] for x in report['failures']}


def test_single_configuration_finite_forces_are_reported_without_stability_claim():
    _,system,_=specimen()
    coordinates=np.array([[0,0,0],[0.15,0,0],[0.2,0.14,0],[0.3,0.15,0.1]])
    report=check_finite_forces(system,coordinates)
    assert report['accepted'] and math.isfinite(report['potential_energy_kj_mol'])
    assert 'not stability' in report['meaning']
    coordinates[1,0]=np.nan
    assert not check_finite_forces(system,coordinates)['accepted']


def trajectory_control():
    topology,system,manifest=specimen()
    report=validate_bonded_site(topology,system,manifest)
    frame=np.array([[0,0,0],[0.15,0,0],[0.2,0.14,0],[0.3,0.15,0.1]])
    frames=np.stack([frame,frame.copy()]);frames[1,3,2]+=0.015
    kwargs={'parameter_report':report,'bond_bounds':[{'atoms':term['atoms'],'minimum_nm':0.13,'maximum_nm':0.18} for term in manifest['bonds']],
            'motion_atom_ids':['c1','s1','c2','c3'],'forces_kj_mol_nm':np.zeros_like(frames)}
    return manifest,frames,kwargs


def test_saved_atom_mapping_is_used_when_trajectory_order_changes():
    manifest,frames,kwargs=trajectory_control()
    result=audit_site_trajectory(['c3','c2','s1','c1'],frames[:,::-1],manifest,**kwargs)
    assert result['accepted'],result['failures']
    assert result['motion']['maximum_internal_motion_nm']>1e-4
    assert result['bonds'][0]['mean_nm']==0.15


def test_wrapped_periodic_site_is_not_reported_as_a_broken_bond():
    manifest,frames,kwargs=trajectory_control()
    frames[:,:,0]+=1.9;frames%=2
    result=audit_site_trajectory(['c1','s1','c2','c3'],frames,manifest,
                                 boxes_nm=np.repeat((np.eye(3)*2)[None,:,:],2,axis=0),**kwargs)
    assert result['accepted'],result['failures']
    assert abs(result['bonds'][0]['mean_nm']-0.15)<1e-12
    result=audit_site_trajectory(['c1','s1','c2','c3'],frames,manifest,**kwargs)
    assert any(x['kind']=='bond_distance_outside_declared_bounds' for x in result['failures'])


def test_reduced_triclinic_image_search_is_not_fractional_rounding_only():
    box=np.array([[2,0,0],[0.95,2,0],[0.2,0.3,2]])
    delta=np.array([[0.49,0.49,0]]) @ box
    image=_minimum_image(delta,box)
    assert np.linalg.norm(image)<np.linalg.norm(delta)-0.5
    np.testing.assert_allclose(image,delta-box[0],atol=1e-12)
    bad=box.copy();bad[1,0]=1.5
    import pytest
    with pytest.raises(ValueError,match='reduced'):_minimum_image(delta,bad)


def test_rigid_translation_and_repeated_frames_do_not_count_as_internal_motion():
    for translate in [False,True]:
        manifest,frames,kwargs=trajectory_control()
        frames[1]=frames[0]+(0.1 if translate else 0)
        result=audit_site_trajectory(['c1','s1','c2','c3'],frames,manifest,**kwargs)
        assert not result['accepted']
        assert any(x['kind']=='no_meaningful_internal_motion' for x in result['failures'])


def test_trajectory_requires_parameter_and_force_evidence():
    manifest,frames,kwargs=trajectory_control()
    kwargs['forces_kj_mol_nm']=None
    result=audit_site_trajectory(['c1','s1','c2','c3'],frames,manifest,**kwargs)
    assert any(x['kind']=='force_evidence_not_supplied' for x in result['failures'])
    kwargs['parameter_report']['accepted']=False
    result=audit_site_trajectory(['c1','s1','c2','c3'],frames,manifest,**kwargs)
    assert result['failures'][0]['kind']=='missing_or_inconsistent_parameter_validation'


def test_broken_declared_coordination_and_nonfinite_forces_are_visible():
    manifest,frames,kwargs=trajectory_control()
    frames[1,0,0]-=0.1
    kwargs['forces_kj_mol_nm'][1,1,1]=float('nan')
    result=audit_site_trajectory(['c1','s1','c2','c3'],frames,manifest,
        coordination_sites=[{'metal_atom_id':'s1','donor_atom_ids':['c1','c2']}],**kwargs)
    assert not result['accepted']
    assert result['coordination_sites'][0]['minimum_donors_within_declared_bounds']==1
    assert {'invalid_or_nonfinite_frame_forces','bond_distance_outside_declared_bounds'} <= {x['kind'] for x in result['failures']}


def test_missing_saved_atom_id_and_undeclared_monitor_bond_are_rejected():
    manifest,frames,kwargs=trajectory_control()
    result=audit_site_trajectory(['c1','s1','c2','unrelated'],frames,manifest,**kwargs)
    assert result['failures'][0]['kind']=='trajectory_atom_map_missing_ids'
    kwargs['bond_bounds'][0]['atoms']=['c1','c3']
    result=audit_site_trajectory(['c1','s1','c2','c3'],frames,manifest,**kwargs)
    assert result['failures'][0]['kind']=='invalid_trajectory_or_monitor_contract'


def test_validation_hashes_bind_system_and_exact_expected_manifest():
    topology,system,manifest=specimen()
    report=validate_bonded_site(topology,system,manifest)
    assert report['system_xml_sha256']==hashlib.sha256(mm.XmlSerializer.serialize(system).encode()).hexdigest()
    assert report['manifest_canonical_sha256']==hashlib.sha256(json.dumps(manifest,sort_keys=True,separators=(',', ':'),allow_nan=False).encode()).hexdigest()
    manifest,frames,kwargs=trajectory_control()
    manifest['atoms'][0]['charge_e']+=0.1
    result=audit_site_trajectory(['c1','s1','c2','c3'],frames,manifest,**kwargs)
    assert result['failures'][0]['kind']=='missing_or_inconsistent_parameter_validation'


def test_zero_epsilon_particle_sigma_is_physically_irrelevant_only_at_zero_epsilon():
    topology,system,manifest=specimen()
    system.getForce(3).setParticleParameters(0,-0.2,0.089,0)
    manifest['atoms'][0].update(sigma_nm=0,epsilon_kj_mol=0)
    assert validate_bonded_site(topology,system,manifest)['accepted']
    system.getForce(3).setParticleParameters(0,-0.2,0.089,0.001)
    assert not validate_bonded_site(topology,system,manifest)['accepted']


def test_water_geometry_term_requires_explicit_native_label_and_shared_oxygen():
    topology,system,manifest=specimen()
    # Recast the three bonded atoms into O-H-H with O-H topology edges only.
    atoms=list(topology.atoms());atoms[0].element=app.element.oxygen
    atoms[1].element=atoms[2].element=app.element.hydrogen
    for index,element in enumerate(['O','H','H']):manifest['atoms'][index]['element']=element
    topology._bonds=[bond for bond in topology._bonds if {bond[0].index,bond[1].index}!={1,2}]
    topology.addBond(atoms[0],atoms[2])
    system.getForce(0).addBond(0,2,0.15,83680)
    manifest['bonds'].append({'atoms':['c1','c2'],'length_nm':0.15,'k_kj_mol_nm2':83680})
    manifest['atoms'][1]['native_residue_name']=manifest['atoms'][2]['native_residue_name']='WAT'
    assert not validate_bonded_site(topology,system,manifest)['accepted']
    manifest['bonds'][1]['topology_role']='native_water_hh_geometry'
    result=validate_bonded_site(topology,system,manifest)
    assert result['accepted'],result['failures']
    assert result['native_water_geometry_terms_without_topology_edges']==[['s1','c2']]
    manifest['atoms'][1]['element']='Zn';atoms[1].element=app.element.zinc
    result=validate_bonded_site(topology,system,manifest)
    assert any(x['kind']=='invalid_native_water_geometry_classification' for x in result['failures'])
