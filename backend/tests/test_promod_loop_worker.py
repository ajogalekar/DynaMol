"""Pure contract controls; native fragment/refinement evidence is recorded separately."""
import copy
import json
import pytest

from backend.promod_loop_worker import LoopConstructionError, construct_loops, main, validate_input


def request():
    rows=[{'chain':'A','resid':str(i),'insertion_code':'','residue':'ALA',
           'one_letter':'A','observed':i!=2} for i in range(1,4)]
    return {'input_pdb':'observed.pdb','chains':[{'chain_id':'A','residues':rows}],
            'max_res_extension':0}


def test_exact_internal_source_contract_preserved():
    config=request();before=copy.deepcopy(config)
    assert validate_input(config)==before
    assert config==before


@pytest.mark.parametrize('change,message',[
    (lambda c:c.update(max_res_extension=1),'extend'),
    (lambda c:c['chains'][0]['residues'][0].update(observed=False),'internal'),
    (lambda c:c['chains'][0]['residues'][1].update(residue='TPO',one_letter='T'),'standard'),
    (lambda c:c['chains'][0]['residues'][1].update(one_letter='G'),'identities'),
    (lambda c:c['chains'][0]['residues'][1].update(resid='1'),'ambiguous'),
    (lambda c:c['chains'][0].update(chain_id='A or all'),'canonical'),
    (lambda c:c['chains'][0]['residues'][1].update(observed=0),'observed flag'),
])
def test_unsafe_or_ambiguous_requests_rejected(change,message):
    config=request();change(config)
    with pytest.raises(ValueError,match=message):validate_input(config)


def test_missing_gap_limit_and_no_terminal_extension():
    config=request()
    config['chains'][0]['residues']=[{'chain':'A','resid':str(i),'insertion_code':'',
        'residue':'ALA','one_letter':'A','observed':i in (0,13)} for i in range(14)]
    validate_input(config)
    config['chains'][0]['residues'].insert(13,{'chain':'A','resid':'12A','insertion_code':'',
        'residue':'ALA','one_letter':'A','observed':False})
    with pytest.raises(ValueError,match='12 residues'):validate_input(config)


def test_failed_contract_publishes_actionable_json_without_native_imports(tmp_path):
    source=tmp_path/'request.json';out=tmp_path/'result.json'
    source.write_text(json.dumps({'input_pdb':'unused.pdb','chains':[]}))
    assert main([str(source),str(out)])==1
    assert json.loads(out.read_text())['status']=='failed'
    assert 'chain alignments' in json.loads(out.read_text())['error']
    assert not list(tmp_path.glob('result.json.*'))


class Search:
    def __init__(self, finish_at=None):
        self.finish_at=finish_at
        self.calls=[]
    def _search(self, method, model, kwargs):
        self.calls.append((method,kwargs))
        if len(self.calls)==self.finish_at:
            model.gaps=[]
    def FillLoopsByDatabase(self,model,*args,**kwargs):
        self._search('database',model,kwargs)
    def FillLoopsByMonteCarlo(self,model,*args,**kwargs):
        self._search('monte_carlo',model,kwargs)


@pytest.mark.parametrize('finish_at',[1,2,3,4,5])
def test_search_stops_on_success_and_preserves_bounded_fallbacks(finish_at):
    from types import SimpleNamespace
    model=SimpleNamespace(gaps=['A.ALA1-(G)-A.ALA3'])
    search=Search(finish_at); progress=[]
    report=construct_loops(model,search,None,None,None,lambda value:progress.append(copy.deepcopy(value)))
    assert len(search.calls)==finish_at
    assert report['gaps_after']==[]
    assert report['observed_context_changes_exported'] is False
    assert all(a['status']=='complete' for a in report['attempts'])
    assert all(a['context_extension_residues']<=2 for a in report['attempts'])
    assert search.calls[0][1]['min_loops_required']==4
    if finish_at>=2:assert search.calls[1][1]['min_loops_required']==1
    if finish_at==5:
        assert search.calls[-1][1]['max_extension']==2
        assert search.calls[-1][1]['mc_num_loops']==6
        assert search.calls[-1][1]['mc_steps']==5000
    assert all(kwargs['ring_punch_detection']==1 for _,kwargs in search.calls)
    assert len(progress)==2*finish_at


def test_exhausted_search_reports_all_attempts_without_accepting_candidate():
    from types import SimpleNamespace
    model=SimpleNamespace(gaps=['A.ALA1-(G)-A.ALA3'])
    with pytest.raises(LoopConstructionError,match='does not establish.*impossible') as error:
        construct_loops(model,Search(),None,None,None)
    assert len(error.value.report['attempts'])==5
    assert error.value.report['gaps_after']==model.gaps


def test_native_error_is_not_misrepresented_as_a_completed_search():
    from types import SimpleNamespace
    search=Search()
    def broken(*args,**kwargs):raise RuntimeError('native database corrupt')
    search.FillLoopsByDatabase=broken
    with pytest.raises(LoopConstructionError,match='database corrupt') as error:
        construct_loops(SimpleNamespace(gaps=['gap']),search,None,None,None)
    assert len(error.value.report['attempts'])==1
    assert error.value.report['attempts'][0]['status']=='failed'


def test_context_search_uses_observed_anchor_ranking_path_without_expanding_output_request():
    from types import SimpleNamespace
    model=SimpleNamespace(gaps=['bounded gap']);search=Search();calls=[]
    def rank(extension):
        calls.append(extension);model.gaps=[]
        return [{'selected_context_gap':'one extra context residue','candidate_count':4}]
    report=construct_loops(model,search,None,None,None,context_search=rank)
    assert calls==[1]
    assert len(search.calls)==2  # Original-stem database searches precede context search.
    assert report['attempts'][-1]['gap_searches'][0]['candidate_count']==4
    assert report['observed_context_changes_exported'] is False


def _ranked_native_search(candidates, **ranking_options):
    """Small native API adapter; geometry assertions are independent below."""
    from types import SimpleNamespace
    from backend.promod_loop_worker import context_fragment_search

    class Gap:
        length=1
        full_seq='AAA'
        before=SimpleNamespace(number=SimpleNamespace(num=1), IsValid=lambda:True)
        after=SimpleNamespace(IsValid=lambda:True)
        def IsTerminal(self): return False
        def GetChainIndex(self): return 0
        def Copy(self): return Gap()
        def __str__(self): return 'A.ALA1-(A)-A.ALA3'

    class Candidates(list):
        def ApplyCCD(self,*args): pass
        def CalculateBackboneScores(self,*args): pass

    selected=[]
    collection=Candidates(candidates)
    model=SimpleNamespace(gaps=[Gap()], model=SimpleNamespace(chains=[SimpleNamespace(name='A')]), seqres=['AAA'])
    methods=SimpleNamespace(
        IsBackboneScoringSetUp=lambda model:True,
        FullGapExtender=lambda *args:SimpleNamespace(Extend=lambda:False),
        CountEnclosedInsertions=lambda *args:1,
        LoopCandidates=SimpleNamespace(FillFromDatabase=lambda *args:collection),
        FilterCandidates=lambda *args:None,
        ScoreContainer=lambda:SimpleNamespace(LinearCombine=lambda weights:[c.score for c in collection]),
        ScoringWeights=SimpleNamespace(GetWeights=lambda:None),
        InsertLoopClearGaps=lambda model, candidate, gap:selected.append(candidate))
    observed={('A',resid,atom):(0.,0.,0.) for resid in (1,3) for atom in ('N','CA','C','O')}
    result=context_fragment_search(model, methods, SimpleNamespace(MaxFragLength=lambda:20), None, None, 1, observed, **ranking_options)
    return selected,result


class _BackboneCandidate:
    def __init__(self,label, backbone_delta=0., oxygen_delta=0.,score=0.):
        self.label,self.backbone_delta,self.oxygen_delta,self.score=label,backbone_delta,oxygen_delta,score
    def __len__(self): return 3
    def _point(self,delta):
        from types import SimpleNamespace
        return SimpleNamespace(x=delta,y=0.,z=0.)
    def GetN(self,index): return self._point(self.backbone_delta)
    GetCA=GetN
    GetC=GetN
    def GetO(self,index): return self._point(self.oxygen_delta)
    def Copy(self): return copy.copy(self)


def test_anchor_ranking_includes_observed_carbonyl_orientation():
    # The first fits N/CA/C exactly but puts the observed oxygen 2.1 Å away.
    # A small coherent displacement must outrank that misleading stem match.
    wrong_oxygen=_BackboneCandidate('wrong oxygen',oxygen_delta=2.1,score=-100.)
    coherent=_BackboneCandidate('coherent',backbone_delta=.15,oxygen_delta=.15)
    selected,rows=_ranked_native_search([wrong_oxygen,coherent],include_carbonyl_anchors=True)
    assert selected[0].label=='coherent'
    assert selected[0] is not coherent  # Native vector lifetime cannot alias selection.
    assert rows[0]['anchor_atoms']==['N','CA','C','O']
    assert rows[0]['anchor_backbone_atoms_compared']==8
    assert rows[0]['maximum_anchor_backbone_displacement_angstrom']==pytest.approx(.15)
    assert rows[0]['anchor_backbone_rmsd_angstrom']==pytest.approx(.15)


def test_candidate_budget_is_enforced_inside_each_native_collection():
    candidates=[_BackboneCandidate(str(i),oxygen_delta=1.) for i in range(40)]
    candidates.append(_BackboneCandidate('outside budget'))
    selected,rows=_ranked_native_search(candidates)
    assert rows[0]['candidate_count']==rows[0]['maximum_candidates']==40
    assert selected[0].label!='outside budget'


def test_incompatible_low_displacement_fragment_does_not_hide_compatible_alternative():
    candidates=[_BackboneCandidate('cis',score=-100),_BackboneCandidate('trans',backbone_delta=.15)]
    selected,rows=_ranked_native_search(candidates,candidate_policy=lambda candidate,gap:
                                       {'compatible':candidate.label=='trans','reason':'test peptide state'})
    assert selected[0].label=='trans'
    assert rows[0]['examined_candidates']==2
    assert rows[0]['candidate_count']==1
    assert rows[0]['policy_rejections'][0]['candidate_number']==1


def test_policy_rejections_still_consume_the_native_candidate_budget():
    candidates=[_BackboneCandidate('incompatible')for _ in range(40)]+[_BackboneCandidate('compatible')]
    selected,rows=_ranked_native_search(candidates,candidate_policy=lambda candidate,gap:
                                       {'compatible':candidate.label=='compatible'})
    assert selected==[]
    assert rows[0]['examined_candidates']==40
    assert len(rows[0]['policy_rejections'])==40
    assert rows[0]['candidate_count']==0


def test_nonfinite_observed_carbonyl_is_never_ranked_as_a_candidate():
    selected,rows=_ranked_native_search([_BackboneCandidate('nonfinite',oxygen_delta=float('nan'))],include_carbonyl_anchors=True)
    assert selected==[]
    assert rows[0]['candidate_count']==0


def test_experimental_ranking_remains_off_in_the_application_default():
    selected,rows=_ranked_native_search([
        _BackboneCandidate('original ranking',oxygen_delta=2.1),
        _BackboneCandidate('carbonyl ranking',backbone_delta=.15,oxygen_delta=.15)])
    assert selected[0].label=='original ranking'
    assert rows[0]['anchor_atoms']==['N','CA','C']
    config=request()
    assert 'research_include_carbonyl_anchors' not in validate_input(config)
    config['research_include_carbonyl_anchors']='true'
    with pytest.raises(ValueError,match='boolean'):validate_input(config)
