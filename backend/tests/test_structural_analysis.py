"""Analytic rigid-body/fluctuation/PBC fixtures, not simulated molecular evidence."""
import copy
import mdtraj as md
import numpy as np
import pytest
from backend import storage
from backend.structural_analysis import structural_analysis, StructuralAnalysisRequest

@pytest.fixture
def fixture(monkeypatch):
    top=md.Topology(); chain=top.add_chain('A')
    for i in range(3):
        residue=top.add_residue('ALA',chain,resSeq=i+1); top.add_atom('CA',md.element.carbon,residue)
    top.add_atom('CB',md.element.carbon,residue)
    for i in range(1,4): top.add_bond(top.atom(0),top.atom(i))
    ref=np.array([[0,0,0],[.1,0,0],[0,.1,0],[0,0,.2]],dtype=np.float32)
    def create(xyz, periodic=False):
        t=md.Trajectory(np.asarray(xyz,dtype=np.float32),top,time=np.arange(len(xyz))*.25)
        if periodic:
            t.unitcell_lengths=np.ones((len(xyz),3));t.unitcell_angles=np.ones((len(xyz),3))*90
        meta={'id':'analytic','n_atoms':4,'n_frames':len(xyz),'time_unit':'ps','atoms':[{'index':a.index,'name':a.name,'element':'C','category':'protein','residue':'ALA','resid':a.residue.resSeq,'chain':'A'} for a in top.atoms]}
        monkeypatch.setattr(storage,'load_physical',lambda _:t[:])
        monkeypatch.setattr(storage,'get_dataset',lambda _:copy.deepcopy(meta))
        return t
    return ref,create

def run(**kwargs): return structural_analysis('analytic',StructuralAnalysisRequest(**kwargs))

def test_rigid_body_fit_and_source_unchanged(fixture):
    ref,create=fixture; rot=np.array([[0,-1,0],[1,0,0],[0,0,1]])
    t=create([ref,ref@rot+[.8,-.2,.5]])
    before=t.xyz.copy()
    r=run(selection='protein-heavy')
    assert r['values']==pytest.approx([0,0],abs=2e-6)
    np.testing.assert_array_equal(t.xyz,before)
    assert r['unit']=='Å' and r['x']==[0,.25]

def test_cartesian_translation(fixture):
    ref,create=fixture;create([ref,ref+[.2,0,0]])
    assert run(alignment='none')['values']==pytest.approx([0,2],abs=1e-6)

def test_alignment_never_reflects_chirality(fixture):
    ref,create=fixture
    create([ref, ref * [-1, 1, 1]])
    result=run(selection='protein-heavy',alignment='selection')
    assert result['values'][0] == pytest.approx(0,abs=1e-6)
    assert result['values'][1] > .5

def test_nondegenerate_fit_matches_independent_mdtraj_rmsd(fixture):
    ref,create=fixture
    rng=np.random.default_rng(731)
    frames=ref + rng.normal(0,.009,(12,4,3))
    t=create(frames)
    expected=md.rmsd(t[:],t[:],frame=0,atom_indices=np.arange(4),parallel=False)*10
    result=run(selection='protein-heavy',alignment='selection')
    assert result['values'] == pytest.approx(expected,abs=2e-4)

def test_known_atom_rmsf_and_reference_outside_window(fixture):
    ref,create=fixture; lo=ref.copy(); hi=ref.copy();lo[3,2]=.1;hi[3,2]=.3
    create([ref,lo,hi])
    result=run(kind='rmsf',selection='custom',atoms=[3],start_frame=1,end_frame=2,reference_frame=0)
    assert result['values']==pytest.approx([1],abs=1e-6)
    assert result['frame_indices']==[1,2] and result['atom_groups']==[[3]]
    assert run(selection='custom',atoms=[3],start_frame=1,end_frame=2)['values']==pytest.approx([1,1],abs=1e-6)

def test_residue_rmsf_uses_mean_variance_not_mean_rmsf(fixture):
    ref,create=fixture; alt=ref.copy();alt[3,2]+=.2
    create([ref,alt])
    result=run(kind='rmsf',selection='custom',atoms=[2,3])
    assert result['values']==pytest.approx([np.sqrt(.5)],abs=1e-6)
    assert result['atom_groups']==[[2,3]]

def test_periodic_translation_is_removed_after_making_whole(fixture):
    ref,create=fixture;create([(ref+.91)%1,(ref+1.03)%1],periodic=True)
    result=run(selection='protein-heavy')
    assert result['values']==pytest.approx([0,0],abs=4e-6)
    assert 'periodic' in result['warnings'][0]
    assert max(run(selection='protein-heavy',alignment='none',periodic='cartesian')['values'])>2

@pytest.mark.parametrize('options',[
    {'atoms':[0,0], 'selection':'custom'}, {'atoms':[999], 'selection':'custom'},
    {'start_frame':1,'end_frame':0}, {'reference_frame':20}, {'end_frame':20},
    {'selection':'custom','atoms':[0,1],'alignment':'selection'},
    {'kind':'rmsf','start_frame':0,'end_frame':0},
])
def test_invalid_requests_are_rejected(fixture,options):
    ref,create=fixture;create([ref,ref])
    with pytest.raises(ValueError):run(**options)

def test_collinear_fit_is_rejected(fixture):
    ref,create=fixture;ref[2]=[.2,0,0];create([ref,ref])
    with pytest.raises(ValueError,match='collinear'):run()

def test_frame_stride_and_atom_output_mapping(fixture):
    ref,create=fixture;create([ref,ref,ref,ref,ref])
    result=run(kind='rmsf',selection='custom',atoms=[3,0],residue_average=False,stride=2)
    assert result['frame_indices']==[0,2,4]
    assert result['atom_groups']==[[3],[0]] and result['values']==pytest.approx([0,0])

@pytest.mark.parametrize('times',[[0,1,0],[0,0,1],[0,np.nan,2]])
def test_bad_timestamps_use_explicit_frame_axis(fixture,times):
    ref,create=fixture;t=create([ref,ref,ref]);t.time[:]=times
    result=run()
    assert result['x']==[0,1,2] and result['x_label']=='Saved frame index'
    assert any('timestamps' in w for w in result['warnings'])
