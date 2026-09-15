import copy

import numpy as np
import pytest
from openmm import app

from backend.loop_context import (AnchoredConstructionForceField, admit_loop_proposal,
                                  affected_reference_keys, outside_observed_torsions,
                                  outer_torsion_anchor_indices, peptide_neighbors)
from backend.residue_identity import residue_key


SOURCE = 'a' * 64


def atom_key(atom):
    return (*residue_key(atom.residue), atom.name)


def fixture_protein(count=10):
    """Connected synthetic atoms; no physical conformation is claimed/tested."""
    top = app.Topology()
    chain = top.addChain('Z')
    residues, names, points, previous = [], {}, [], None
    elements = {'N': app.element.nitrogen, 'CA': app.element.carbon,
                'C': app.element.carbon, 'O': app.element.oxygen,
                'CB': app.element.carbon, 'H': app.element.hydrogen}
    offsets = {'N': [0, 0, 0], 'CA': [.145, .04, 0], 'C': [.24, .12, .05],
               'O': [.25, .23, .08], 'CB': [.17, -.07, .11], 'H': [-.03, .01, .05]}
    for i in range(1, count+1):
        # Nonsequential residue numbers/insertion codes ensure topology drives scope.
        residue = top.addResidue('ALA', chain, str(100+i*3), insertionCode='B' if i == 5 else '')
        residues.append(residue)
        local = {}
        for name, element in elements.items():
            local[name] = top.addAtom(name, element, residue)
            points.append(np.array(offsets[name]) + [i*.4, i*.1, 0])
        for a, b in [('N', 'CA'), ('CA', 'C'), ('C', 'O'), ('CA', 'CB'), ('N', 'H')]:
            top.addBond(local[a], local[b])
        if previous is not None:
            top.addBond(previous, local['N'])
        previous = local['C']
        names[i] = {n: a.index for n, a in local.items()}
    xyz = np.array(points)
    return top, residues, names, xyz


def proposal_for(top, xyz, modeled, context):
    def rows(keys, backbone_only):
        return [{'identity': list(atom_key(a)), 'element': a.element.symbol,
                 'xyz_nm': (xyz[a.index] + [.005, .008, -.002]).tolist()}
                for a in top.atoms() if residue_key(a.residue) in keys
                and a.element != app.element.hydrogen
                and (not backbone_only or a.name in {'N', 'CA', 'C', 'O'})]
    return {'schema_version': 1, 'status': 'candidate', 'source_sha256': SOURCE,
            'modeled_residue_keys': [list(k) for k in sorted(modeled)],
            'modeled_heavy_atoms': rows(modeled, False),
            'observed_context_atoms': rows(context, True)}


@pytest.fixture
def sample():
    top, residues, named, xyz = fixture_protein()
    modeled = {residue_key(residues[i-1]) for i in [4, 5]}
    flanks = {residue_key(residues[i-1]) for i in [3, 6]}
    return top, residues, named, xyz, modeled, flanks, proposal_for(top, xyz, modeled, flanks)


def admit(sample, proposal=None, **kwargs):
    top, residues, named, xyz, modeled, flanks, current = sample
    return admit_loop_proposal(current if proposal is None else proposal, top, xyz,
                               modeled, flanks, source_sha256=SOURCE, **kwargs)


def test_complete_seed_preserves_source_and_hydrogens_and_records_context(sample):
    top, residues, names, xyz, modeled, flanks, proposal = sample
    before, original = xyz.copy(), copy.deepcopy(proposal)
    result = admit(sample)
    assert result.modeled_keys == modeled and result.proposed_flank_keys == flanks
    assert result.mobile_flank_keys == flanks
    assert result.extra_fixed_indices == {names[3]['N'], names[3]['CA'], names[6]['CA'], names[6]['C']}
    assert np.array_equal(xyz, before) and proposal == original
    hydrogens = [a.index for a in top.atoms() if a.element == app.element.hydrogen]
    outside = [a.index for a in top.atoms() if residue_key(a.residue) not in modeled | flanks]
    assert np.array_equal(result.seeded_xyz_nm[hydrogens], xyz[hydrogens])
    assert np.array_equal(result.seeded_xyz_nm[outside], xyz[outside])
    for i in result.extra_fixed_indices:
        assert np.array_equal(result.seeded_xyz_nm[i], xyz[i])
    assert set(result.provenance['transferred_context_sidechain_atoms']) == {
        (*residue_key(residues[i-1]), 'CB') for i in [3, 6]}
    assert not result.provenance['complete_candidate_accepted']
    assert result.provenance['required_final_observed_heavy_cap_nm'] == .1
    with pytest.raises(ValueError):
        result.seeded_xyz_nm.setflags(write=True)
    with pytest.raises(TypeError):
        result.provenance['complete_candidate_accepted'] = True


def test_context_sidechain_transfer_preserves_source_residue_frame_geometry(sample):
    top, residues, names, xyz, modeled, flanks, proposal = sample
    result = admit(sample, observed_context_policy='mobile_flanks')
    for i in [3, 6]:
        delta = result.seeded_xyz_nm[names[i]['CB']] - xyz[names[i]['CB']]
        assert np.allclose(delta, [.005, .008, -.002], rtol=0, atol=1e-14)


def test_fixed_observed_policy_restores_entire_context_and_does_not_transfer(sample):
    top, residues, names, xyz, modeled, flanks, proposal = sample
    result = admit(sample, observed_context_policy='fixed_observed')
    assert not result.mobile_flank_keys and not result.extra_fixed_indices
    assert result.proposed_flank_keys == flanks
    observed = [a.index for a in top.atoms() if residue_key(a.residue) not in modeled]
    assert np.array_equal(xyz[observed], result.seeded_xyz_nm[observed])
    assert not result.provenance['transferred_context_sidechain_atoms']


def test_raw_observed_cap_is_not_an_admission_gate(sample):
    proposal = copy.deepcopy(sample[-1])
    for row in proposal['observed_context_atoms']:
        row['xyz_nm'][2] += .2
    result = admit(sample, proposal)
    assert result.provenance['maximum_proposed_context_displacement_nm'] > .1
    assert result.provenance['maximum_seeded_context_displacement_nm'] > .1
    assert result.provenance['required_final_observed_heavy_cap_nm'] == .1
    assert not result.provenance['complete_candidate_accepted']


@pytest.mark.parametrize('change,match', [
    ('subset', 'exact complete modeled'), ('extra_residue', 'exact complete modeled'),
    ('duplicate_residue', 'Duplicate proposal modeled'), ('missing_sidechain', 'every modeled heavy atom'),
    ('missing_backbone', 'every modeled heavy atom'), ('duplicate_atom', 'Duplicate proposal atom'),
    ('hydrogen', 'Unexpected hydrogen'), ('unknown_atom', 'Unknown proposal atom'),
    ('element', 'element differs'), ('wrong_group', 'groups disagree'),
    ('missing_context_oxygen', 'complete N/CA/C/O'), ('missing_context_field', 'explicit observed_context'),
    ('angstrom_only', 'exact identity, element, and xyz_nm'), ('nonfinite', 'Nonfinite'),
    ('wrong_shape', 'exact shape'), ('string_number', 'numeric nanometer'),
    ('wrong_source', 'source SHA256 differs'), ('schema_bool', 'schema_version=1'),
    ('status', 'fresh candidate'),
])
def test_malformed_incomplete_or_stale_proposal_fails(sample, change, match):
    top, residues, names, xyz, modeled, flanks, original = sample
    p = copy.deepcopy(original)
    if change == 'subset':
        p['modeled_residue_keys'].pop()
    elif change == 'extra_residue':
        p['modeled_residue_keys'].append(list(residue_key(residues[0])))
    elif change == 'duplicate_residue':
        p['modeled_residue_keys'].append(p['modeled_residue_keys'][0])
    elif change in {'missing_sidechain', 'missing_backbone'}:
        name = 'CB' if change == 'missing_sidechain' else 'O'
        p['modeled_heavy_atoms'].pop(next(i for i, r in enumerate(p['modeled_heavy_atoms']) if r['identity'][4] == name))
    elif change == 'duplicate_atom':
        p['modeled_heavy_atoms'].append(p['modeled_heavy_atoms'][0])
    elif change == 'hydrogen':
        p['modeled_heavy_atoms'].append({'identity': [*p['modeled_residue_keys'][0], 'H'], 'element': 'H', 'xyz_nm': [0., 0., 0.]})
    elif change == 'unknown_atom':
        p['modeled_heavy_atoms'][0]['identity'][4] = 'NX'
    elif change == 'element':
        p['modeled_heavy_atoms'][0]['element'] = 'C'
    elif change == 'wrong_group':
        p['observed_context_atoms'].append(p['modeled_heavy_atoms'].pop())
    elif change == 'missing_context_oxygen':
        p['observed_context_atoms'] = [r for r in p['observed_context_atoms'] if r['identity'][4] != 'O']
    elif change == 'missing_context_field':
        p.pop('observed_context_atoms')
    elif change == 'angstrom_only':
        p['modeled_heavy_atoms'][0]['xyz_angstrom'] = p['modeled_heavy_atoms'][0].pop('xyz_nm')
    elif change == 'nonfinite':
        p['modeled_heavy_atoms'][0]['xyz_nm'][0] = float('nan')
    elif change == 'wrong_shape':
        p['modeled_heavy_atoms'][0]['xyz_nm'].pop()
    elif change == 'string_number':
        p['modeled_heavy_atoms'][0]['xyz_nm'][0] = '0.1'
    elif change == 'wrong_source':
        p['source_sha256'] = 'b' * 64
    elif change == 'schema_bool':
        p['schema_version'] = True
    elif change == 'status':
        p['status'] = 'accepted'
    with pytest.raises(ValueError, match=match):
        admit(sample, p)


def test_unapproved_or_nonadjacent_context_is_rejected(sample):
    top, residues, names, xyz, modeled, flanks, proposal = sample
    wrong = {residue_key(residues[0])}
    p = proposal_for(top, xyz, modeled, wrong)
    with pytest.raises(ValueError, match='unapproved'):
        admit(sample, p)
    with pytest.raises(ValueError, match='not an immediate peptide flank'):
        admit_loop_proposal(p, top, xyz, modeled, wrong, source_sha256=SOURCE)


def test_modified_observed_context_rejected_even_when_explicitly_allowed(sample):
    top, residues, names, xyz, modeled, flanks, proposal = sample
    residues[2].name = 'TPO'
    new_flanks = {residue_key(residues[i-1]) for i in [3, 6]}
    p = proposal_for(top, xyz, modeled, new_flanks)
    with pytest.raises(ValueError, match='Unknown or modified'):
        admit_loop_proposal(p, top, xyz, modeled, new_flanks, source_sha256=SOURCE)


def test_protected_metal_contact_is_detected_from_source_before_proposal_moves_it(sample):
    top, residues, names, xyz, modeled, flanks, proposal = sample
    chain = top.addChain('M')
    metal = top.addResidue('ZN', chain, '1')
    top.addAtom('ZN', app.element.zinc, metal)
    expanded = np.vstack([xyz, xyz[names[3]['O']] + [0, 0, .2]])
    for r in proposal['observed_context_atoms']:
        r['xyz_nm'][2] += 10
    with pytest.raises(ValueError, match='protected metal'):
        admit_loop_proposal(proposal, top, expanded, modeled, flanks, source_sha256=SOURCE)


def test_explicit_specialized_protection_cannot_be_overridden(sample):
    with pytest.raises(ValueError, match='protected metal'):
        admit(sample, protected_residues={next(iter(sample[5]))})


@pytest.mark.parametrize('named_cn', [False, True])
def test_chemical_crosslink_and_crosschain_named_cn_are_protected(sample, named_cn):
    top, residues, names, xyz, modeled, flanks, proposal = sample
    atoms = list(top.atoms())
    if named_cn:
        chain = top.addChain('L')
        ligand = top.addResidue('LIG', chain, '1')
        extra = top.addAtom('N', app.element.nitrogen, ligand)
        top.addBond(atoms[names[3]['C']], extra)
        xyz = np.vstack([xyz, [30, 30, 30]])
    else:
        top.addBond(atoms[names[3]['CB']], atoms[names[8]['CB']])
    with pytest.raises(ValueError, match='chemical crosslink'):
        admit_loop_proposal(proposal, top, xyz, modeled, flanks, source_sha256=SOURCE)


@pytest.mark.parametrize('shared', [True, False])
def test_multi_gap_shared_or_adjacent_context_cannot_conflict(shared):
    top, residues, names, xyz = fixture_protein(12)
    modeled_indices, flank_indices = ([4, 6], [3, 5, 7]) if shared else ([4, 7], [3, 5, 6, 8])
    modeled = {residue_key(residues[i-1]) for i in modeled_indices}
    flanks = {residue_key(residues[i-1]) for i in flank_indices}
    p = proposal_for(top, xyz, modeled, flanks)
    with pytest.raises(ValueError, match='multi-gap'):
        admit_loop_proposal(p, top, xyz, modeled, flanks, source_sha256=SOURCE)


def test_disjoint_gaps_require_all_residues_and_can_share_one_complete_proposal():
    top, residues, names, xyz = fixture_protein(14)
    modeled = {residue_key(residues[i-1]) for i in [4, 10, 11]}
    flanks = {residue_key(residues[i-1]) for i in [3, 5, 9, 12]}
    p = proposal_for(top, xyz, modeled, flanks)
    result = admit_loop_proposal(p, top, xyz, modeled, flanks, source_sha256=SOURCE)
    assert sorted(map(len, result.provenance['modeled_components'])) == [1, 2]


def test_missing_actual_peptide_connection_cannot_be_inferred_from_residue_numbers(sample):
    top, residues, names, xyz, modeled, flanks, proposal = sample
    top._bonds = [b for b in top._bonds if {b[0].index, b[1].index} != {names[3]['C'], names[4]['N']}]
    with pytest.raises(ValueError, match='internal unambiguous peptide gap'):
        admit(sample)


def test_duplicate_topology_atom_identity_cannot_silently_overwrite_mapping(sample):
    top, residues, names, xyz, modeled, flanks, proposal = sample
    top.addAtom('CB', app.element.carbon, residues[-1])
    with pytest.raises(ValueError, match='ambiguous residue/atom'):
        admit_loop_proposal(proposal, top, np.vstack([xyz, [0, 0, 0]]), modeled, flanks, source_sha256=SOURCE)


def test_external_support_and_reference_selection_use_actual_changes(sample):
    top, residues, names, xyz, modeled, flanks, proposal = sample
    result = admit(sample)
    assert outside_observed_torsions(top, xyz, result.seeded_xyz_nm, modeled, flanks)['all_defining_coordinates_exactly_preserved']
    def row(i):
        return {'residue': list(residue_key(residues[i-1])),
                'phi_atoms': [names[i-1]['C'], names[i]['N'], names[i]['CA'], names[i]['C']],
                'psi_atoms': [names[i]['N'], names[i]['CA'], names[i]['C'], names[i+1]['N']]}
    report = {'rows': [row(i) for i in [3, 4, 5, 6, 7]], 'unavailable': []}
    expected = {residue_key(residues[i-1]) for i in [3, 4, 5, 6]}
    assert affected_reference_keys(top, xyz, xyz.copy(), modeled, report) == expected
    changed = result.seeded_xyz_nm.copy()
    changed[names[6]['CA'], 0] += 1e-9  # preceding omega only for residue seven
    assert residue_key(residues[6]) in affected_reference_keys(top, xyz, changed, modeled, report)
    assert not outside_observed_torsions(top, xyz, changed, modeled, flanks)['all_defining_coordinates_exactly_preserved']
    unavailable = {'rows': [], 'unavailable': [{'residue': list(residue_key(residues[6]))}]}
    assert residue_key(residues[6]) in affected_reference_keys(top, xyz, xyz.copy(), modeled, unavailable)


def test_ambiguous_peptide_branch_fails(sample):
    top, residues, names, xyz, modeled, flanks, proposal = sample
    atoms = list(top.atoms())
    top.addBond(atoms[names[6]['C']], atoms[names[9]['N']])
    with pytest.raises(ValueError, match='Ambiguous peptide'):
        peptide_neighbors(top)


def test_anchor_masses_do_not_change_the_physical_forcefield():
    class System:
        def __init__(self):
            self.masses = [12., 14., 16.]
        def setParticleMass(self, i, value):
            self.masses[i] = value
    class ForceField:
        def createSystem(self, **kwargs):
            return System()
    physical = ForceField()
    assert AnchoredConstructionForceField(physical, {1}).createSystem().masses == [12., 0., 16.]
    assert physical.createSystem().masses == [12., 14., 16.]


def test_unknown_policy_or_nonfinite_source_fails(sample):
    with pytest.raises(ValueError, match='Unknown observed-context'):
        admit(sample, observed_context_policy='soften_geometry')
    sample[3][0, 0] = np.inf
    with pytest.raises(ValueError, match='Nonfinite'):
        admit(sample)
