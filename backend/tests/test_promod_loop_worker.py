"""Pure contract controls; native fragment/refinement evidence is recorded separately."""
import copy
import json
import pytest

from backend.promod_loop_worker import main, validate_input


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
