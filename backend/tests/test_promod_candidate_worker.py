"""Software contracts only: these controls perform no native loop sampling."""
import copy
import hashlib
import json
import math
from pathlib import Path
import subprocess
import sys
from types import SimpleNamespace as NS

import pytest

from backend import promod_candidate_worker as worker


def peptide_candidate(omega, *, degenerate=False):
    points=[(0.,1.,0.),(0.,0.,0.),(1.,0.,0.),
            (1.,math.cos(math.radians(omega)),math.sin(math.radians(omega)))]
    if degenerate:
        points[2]=points[1]
    class Candidate:
        def __len__(self): return 2
        def GetCA(self,i): return NS(x=points[0 if i==0 else 3][0],y=points[0 if i==0 else 3][1],z=points[0 if i==0 else 3][2])
        def GetC(self,i): return NS(x=points[1][0],y=points[1][1],z=points[1][2])
        def GetN(self,i): return NS(x=points[2][0],y=points[2][1],z=points[2][2])
    return Candidate()


@pytest.mark.parametrize('omega,compatible',[(179,True),(-175,True),(5,False),(-4,False),(90,False),(-90,False),(90.0000005,False)])
def test_new_nonpro_fragment_links_match_construction_basin(omega,compatible):
    rows={1:dict(resid='275',insertion_code='A',residue='GLU',one_letter='E',observed=True),
          2:dict(resid='276',insertion_code='',residue='ASP',one_letter='D',observed=False)}
    result=worker.new_peptide_policy(peptide_candidate(omega),NS(before=NS(number=NS(num=1))),'C',rows)
    assert result['compatible'] is compatible
    row=result['new_nonpro_links'][0]
    assert row['omega_degrees']==pytest.approx(omega)
    assert row['atoms']==[['C','275','A','GLU','CA'],['C','275','A','GLU','C'],
                          ['C','276','','ASP','N'],['C','276','','ASP','CA']]


@pytest.mark.parametrize('observed,letter',[(True,'A'),(False,'P')])
def test_new_trans_policy_does_not_reinterpret_observed_or_pro_peptides(observed,letter):
    rows={1:dict(resid='1',insertion_code='',residue='ALA',one_letter='A',observed=True),
          2:dict(resid='2',insertion_code='',residue='PRO' if letter=='P' else 'ALA',one_letter=letter,observed=observed)}
    result=worker.new_peptide_policy(peptide_candidate(5),NS(before=NS(number=NS(num=1))),'A',rows)
    assert result['compatible'] is True
    assert result['new_nonpro_links']==[]


def test_undefined_new_peptide_is_explicitly_incompatible():
    rows={i:dict(resid=str(i),insertion_code='',residue='ALA',one_letter='A',observed=i==1)for i in (1,2)}
    result=worker.new_peptide_policy(peptide_candidate(180,degenerate=True),NS(before=NS(number=NS(num=1))),'A',rows)
    assert result['compatible'] is False
    assert result['new_nonpro_links'][0]['reason']=='undefined'
    assert result['new_nonpro_links'][0]['omega_degrees'] is None


def request(tmp_path, *, missing=(3,), count=7, strategy='fragment'):
    rows = [dict(chain='A', resid=str(i), insertion_code='', residue='ALA',
                 one_letter='A', observed=i not in missing) for i in range(1, count + 1)]
    config = dict(input_pdb=str(tmp_path / 'source.pdb'), chains=[dict(chain_id='A', residues=rows)],
                  max_res_extension=0, strategy=strategy, seed=23)
    write_source(config)
    return config


def write_source(config):
    lines, serial = [], 1
    for chain in config['chains']:
        for row in chain['residues']:
            if not row['observed']:
                continue
            for index, name in enumerate(('N', 'CA', 'C', 'O', 'CB')):
                x, y, z = float(row['resid']) * 4 + index, index % 2, index % 3
                lines.append(f"ATOM  {serial:5d} {name:>4s} {row['residue']:>3s} {row['chain']}{int(row['resid']):4d}{row['insertion_code'] or ' '}   {x:8.3f}{y:8.3f}{z:8.3f}{1.:6.2f}{20.:6.2f}          {name[0]:>2s}\n")
                serial += 1
    Path(config['input_pdb']).write_text(''.join(lines) + 'END\n')


def atom(name, xyz, element=None):
    return NS(name=name, element=element or name[0], pos=NS(x=xyz[0], y=xyz[1], z=xyz[2]))


def model_for(config, source):
    chains = []
    maps = {}
    for item in config['chains']:
        residues = []
        maps[item['chain_id']] = {}
        for pos, row in enumerate(item['residues'], 1):
            maps[item['chain_id']][pos] = row
            key = worker._key(item['chain_id'], row)
            aa = []
            for index, name in enumerate(('N', 'CA', 'C', 'O', 'CB')):
                xyz = source.get((*key, name), {}).get('xyz_angstrom', [pos * 4 + index, index % 2, index % 3])
                aa.append(atom(name, xyz))
            residues.append(NS(name=row['residue'], number=NS(num=pos), atoms=aa, IsValid=lambda: True))
        chains.append(NS(name=item['chain_id'], residues=residues))
    entity = NS(chains=chains)
    entity.FindResidue = lambda name, pos: next(r for c in chains if c.name == name for r in c.residues if r.number.num == pos)
    return NS(model=entity, gaps=['original gap']), maps


def fake_native(monkeypatch, config):
    source = worker.read_source(config['input_pdb'], config['chains'])
    models, calls = [], []
    def raw(conf, native):
        model, maps = model_for(conf, source)
        models.append(model)
        return model, maps, {key: list(value['xyz_angstrom']) for key, value in source.items()}
    def fragments(model, *args, **kwargs):
        calls.append(('fragment', args[-2], args[-1], kwargs))
        model.gaps = []
        return [{'candidate_count': 1}]
    native = NS(
        versions={'promod3': 'test', 'ost': 'test'},
        ost=NS(mol=NS(ResNum=lambda pos: pos)),
        modelling=NS(SetupDefaultBackboneScoring=lambda model: None,
                     ReconstructSidechains=lambda *args, **kwargs: calls.append(('sidechains', kwargs))),
        loop=NS(LoadFragDB=lambda: object(), LoadStructureDB=lambda: object(),
                LoadTorsionSamplerCoil=lambda **kwargs: calls.append(('torsions', kwargs)) or object()))
    monkeypatch.setattr(worker, '_load_native', lambda: native)
    monkeypatch.setattr(worker, '_raw_model', raw)
    monkeypatch.setattr(worker, '_native_gaps', lambda model, specs: {})
    monkeypatch.setattr(worker.legacy, 'context_fragment_search', fragments)
    return models, calls, native


@pytest.mark.parametrize('change,match', [
    (lambda c: c.update(seed=True), 'integer seed'),
    (lambda c: c.update(seed=-1), 'integer seed'),
    (lambda c: c.update(strategy='automatic'), 'explicit'),
    (lambda c: c.update(max_res_extension=3), 'extend'),
    (lambda c: c.update(max_res_extension=1, strategy='torsion_mc'), 'fragment strategy'),
    (lambda c: c.update(source_sha256='fake'), 'SHA256'),
    (lambda c: c['chains'][0]['residues'][1].update(resid='02'), 'exactly'),
])
def test_bad_requests_fail_before_native_loading(tmp_path, monkeypatch, change, match):
    config = request(tmp_path)
    monkeypatch.setattr(worker, '_load_native', lambda: pytest.fail('must not import native libraries'))
    change(config)
    with pytest.raises(ValueError, match=match):
        worker.run(config)


def test_source_identity_and_element_retained_with_author_insertion_code(tmp_path):
    config = request(tmp_path)
    config['chains'][0]['residues'][1]['insertion_code'] = 'A'
    write_source(config)
    atoms = worker.read_source(config['input_pdb'], config['chains'])
    assert atoms[('A', '2', 'A', 'ALA', 'CA')]['element'] == 'C'
    assert atoms[('A', '2', 'A', 'ALA', 'CA')]['xyz_angstrom'] == [9., 1., 1.]


@pytest.mark.parametrize('mutation,match', [
    (lambda s: s + s.splitlines(True)[0], 'duplicate'),
    (lambda s: s[:16] + 'A' + s[17:], 'alternates'),
    (lambda s: s[:30] + '     nan' + s[38:], 'nonfinite'),
    (lambda s: 'MODEL        1\n' + s + 'MODEL        2\n', 'one coordinate model'),
])
def test_malformed_observed_atoms_stop_the_attempt(tmp_path, mutation, match):
    config = request(tmp_path)
    path = Path(config['input_pdb'])
    path.write_text(mutation(path.read_text()))
    with pytest.raises(ValueError, match=match):
        worker.read_source(path, config['chains'])


def test_missing_source_row_is_integrity_failure(tmp_path):
    config = request(tmp_path)
    path = Path(config['input_pdb'])
    path.write_text(''.join(line for line in path.read_text().splitlines(True) if line[22:26].strip() != '2'))
    with pytest.raises(ValueError, match='identities/order'):
        worker.read_source(path, config['chains'])


def test_gap_specs_cover_every_original_stretch(tmp_path):
    config = request(tmp_path, missing=(3, 4, 9), count=12)
    gaps = worker.gap_specs(config['chains'])
    assert [(g['start'], g['end']) for g in gaps] == [(2, 5), (8, 10)]
    assert [k[1] for g in gaps for k in g['modeled']] == ['3', '4', '9']
    assert [k[1] for g in gaps for k in g['context']] == ['2', '5', '10', '8']


@pytest.mark.parametrize('missing', [(3, 5), (3, 6)])
def test_conflicting_multi_gap_context_is_unavailable(tmp_path, missing):
    config = request(tmp_path, missing=missing, count=9)
    with pytest.raises(worker.CandidateUnavailable, match='stem support'):
        worker.gap_specs(config['chains'])


def test_modified_stem_cannot_escape_native_standardization(tmp_path):
    config = request(tmp_path)
    config['chains'][0]['residues'][1].update(residue='TPO', one_letter='T')
    with pytest.raises(worker.CandidateUnavailable, match='modified'):
        worker.gap_specs(config['chains'])


def test_source_conditioning_has_no_hidden_terminal_defaults(tmp_path):
    config = request(tmp_path, missing=(2,))
    specs = worker.gap_specs(config['chains'])
    maps = {'A': {i + 1: r for i, r in enumerate(config['chains'][0]['residues'])}}
    with pytest.raises(worker.CandidateUnavailable, match='outside both'):
        worker.source_boundaries(specs[0], maps, {}, None)


def test_source_conditioning_uses_observed_atoms_and_neighbor_letters(tmp_path):
    config = request(tmp_path)
    source = worker.read_source(config['input_pdb'], config['chains'])
    maps = {'A': {i + 1: r for i, r in enumerate(config['chains'][0]['residues'])}}
    maps['A'][1]['one_letter'] = 'R'
    maps['A'][5]['one_letter'] = 'Y'
    calls = []
    geom = NS(Vec3=lambda *v: v, DihedralAngle=lambda *p: calls.append(p) or .25)
    values = worker.source_boundaries(worker.gap_specs(config['chains'])[0], maps, source, geom)
    assert values == dict(n_stem_phi=.25, c_stem_psi=.25, prev_aa='R', next_aa='Y')
    assert calls[0][0] == tuple(source[('A', '1', '', 'ALA', 'C')]['xyz_angstrom'])
    del source[('A', '1', '', 'ALA', 'C')]
    with pytest.raises(ValueError, match='support atom'):
        worker.source_boundaries(worker.gap_specs(config['chains'])[0], maps, source, geom)


def test_fresh_attempt_exports_all_missing_atoms_and_complete_context(tmp_path, monkeypatch):
    config = request(tmp_path, missing=(3, 4, 9), count=12)
    config['chains'][0]['residues'][2]['insertion_code'] = 'B'
    before = copy.deepcopy(config)
    models, calls, _ = fake_native(monkeypatch, config)
    first = worker.run(config)
    second = worker.run(config)
    assert models[0] is not models[1]
    assert config == before
    assert first['status'] == second['status'] == 'candidate'
    assert len(first['modeled_heavy_atoms']) == 15
    assert len(first['observed_context_atoms']) == 20
    assert ('A', '3', 'B', 'ALA') in map(tuple, first['modeled_residue_keys'])
    assert first['source_sha256'] == hashlib.sha256(Path(config['input_pdb']).read_bytes()).hexdigest()
    assert all(set(a) == {'identity', 'element', 'xyz_nm'} for a in first['modeled_heavy_atoms'] + first['observed_context_atoms'])
    assert all(call[1] == 0 for call in calls if call[0] == 'fragment')
    assert first['physical_model_validated'] is False


def test_extra_native_search_context_does_not_expand_export_scope(tmp_path, monkeypatch):
    config = request(tmp_path)
    config['max_res_extension'] = 2
    _, calls, _ = fake_native(monkeypatch, config)
    result = worker.run(config)
    assert result['max_res_extension'] == 2
    assert all(call[1] == 2 for call in calls if call[0] == 'fragment')
    assert {row['identity'][1] for row in result['modeled_heavy_atoms']} == {'3'}
    assert {row['identity'][1] for row in result['observed_context_atoms']} == {'2', '4'}


@pytest.mark.parametrize('failure', ['unclosed', 'native_error', 'incomplete_atoms'])
def test_unavailable_is_distinct_from_native_or_inventory_failure(tmp_path, monkeypatch, failure):
    config = request(tmp_path)
    _, _, native = fake_native(monkeypatch, config)
    if failure == 'unclosed':
        monkeypatch.setattr(worker.legacy, 'context_fragment_search', lambda *a, **k: [])
        expected = worker.CandidateUnavailable
    elif failure == 'native_error':
        def crash(*a, **k):
            raise RuntimeError('native database damaged')
        native.loop.LoadFragDB = crash
        expected = RuntimeError
    else:
        def remove_atom(model, **kwargs):
            model.FindResidue('A', 3).atoms.pop()
        native.modelling.ReconstructSidechains = remove_atom
        expected = ValueError
    with pytest.raises(expected):
        worker.run(config)


def test_source_change_overrides_a_candidate(tmp_path, monkeypatch):
    config = request(tmp_path)
    _, _, native = fake_native(monkeypatch, config)
    native.modelling.ReconstructSidechains = lambda *a, **k: Path(config['input_pdb']).write_text('END\n')
    with pytest.raises(ValueError, match='changed during'):
        worker.run(config)


def test_source_hash_mismatch_prevents_native_loading(tmp_path, monkeypatch):
    config = request(tmp_path)
    config['source_sha256'] = '0' * 64
    monkeypatch.setattr(worker, '_load_native', lambda: pytest.fail('native must not load'))
    with pytest.raises(ValueError, match='frozen request'):
        worker.run(config)


def test_export_never_includes_remote_modified_or_nonprotein_chemistry(tmp_path):
    config = request(tmp_path)
    source = worker.read_source(config['input_pdb'], config['chains'])
    model, maps = model_for(config, source)
    model.model.chains.append(NS(name='_', residues=[NS(name='ZN', atoms=[atom('ZN', [1, 2, 3], 'ZN')])]))
    maps['A'][7].update(residue='TPO', one_letter='T')
    model.model.FindResidue('A', 7).name = 'THR'
    modeled, context = worker.export_candidate(model, maps, worker.gap_specs(config['chains']), source)
    assert {a['identity'][3] for a in modeled + context} == {'ALA'}
    assert {a['identity'][1] for a in context} == {'2', '4'}


def test_context_atom_identity_change_is_integrity_failure(tmp_path):
    config = request(tmp_path)
    source = worker.read_source(config['input_pdb'], config['chains'])
    model, maps = model_for(config, source)
    model.model.FindResidue('A', 2).atoms[-1].element = 'S'
    with pytest.raises(ValueError, match='identity/element'):
        worker.export_candidate(model, maps, worker.gap_specs(config['chains']), source)


@pytest.mark.parametrize('exception,status,code', [
    (worker.CandidateUnavailable('no closed result'), 'unavailable', 2),
    (ValueError('bad source identity'), 'failed', 1),
    (RuntimeError('native failure'), 'failed', 1),
])
def test_cli_preserves_failure_classes(tmp_path, monkeypatch, exception, status, code):
    config = request(tmp_path)
    inp, out = tmp_path / 'request.json', tmp_path / 'result.json'
    inp.write_text(json.dumps(config))
    def fail(config):
        raise exception
    monkeypatch.setattr(worker, 'run', fail)
    assert worker.main([str(inp), str(out)]) == code
    assert json.loads(out.read_text())['status'] == status
    assert not list(tmp_path.glob('result.json.*'))


def test_standalone_isolated_bad_request_needs_no_native_imports(tmp_path):
    inp, out = tmp_path / 'request.json', tmp_path / 'result.json'
    inp.write_text('{}')
    process = subprocess.run([sys.executable, '-I', '-B', str(Path(worker.__file__).resolve()), str(inp), str(out)],
                             capture_output=True, text=True, timeout=10)
    assert process.returncode == 1
    assert json.loads(out.read_text())['error_type'] == 'ValueError'
    assert 'observed-coordinate' in json.loads(out.read_text())['error']


def test_existing_output_is_never_overwritten(tmp_path):
    out = tmp_path / 'result.json'
    out.write_text('previous result')
    with pytest.raises(ValueError, match='fresh attempt'):
        worker.main([str(tmp_path / 'absent.json'), str(out)])
    assert out.read_text() == 'previous result'


@pytest.mark.parametrize('strategy', ['torsion_mc', 'fragment_mc'])
def test_native_mc_budget_seed_closure_and_conditioning_contract(tmp_path, monkeypatch, strategy):
    config = request(tmp_path, strategy=strategy)
    source = worker.read_source(config['input_pdb'], config['chains'])
    model, maps = model_for(config, source)
    model.seqres = ['AAAAAAA']
    model.backbone_scorer, model.backbone_scorer_env = object(), object()
    spec = worker.gap_specs(config['chains'])[0]
    gap = NS(full_seq=spec['sequence'], before=object(), after=object(), GetChainIndex=lambda: 0)
    calls = {}
    class Backbone:
        def __len__(self):
            return len(gap.full_seq)
        def Copy(self):
            return Backbone()
    def record(name, result):
        def call(*args, **kwargs):
            calls[name] = (args, kwargs)
            return result
        return call
    methods = NS(
        PhiPsiSampler=record('phi_psi', object()),
        FraggerHandle=record('fragger_handle', NS(GetList=record('fragger_list', [[object()]]))),
        FragmentSampler=record('fragment', object()),
        CCDCloser=record('ccd', object()),
        ScoringWeights=NS(GetWeights=lambda: dict(reduced=1, cb_packing=2, clash=3, unused=4)),
        LinearScorer=record('scorer', NS(GetScore=lambda bb: -4.)),
        ExponentialCooler=record('cooler', object()),
        LoopCandidates=NS(FillFromMonteCarloSampler=record('sample', [Backbone()])),
        InsertLoopClearGaps=record('insert', None))
    native = NS(modelling=methods, geom=NS(Vec3=lambda *v: v, DihedralAngle=lambda *v: .5),
                loop=NS(LoadTorsionSamplerCoil=record('torsions', object()), LoadStructureDB=lambda: object()))
    result = worker._sample_gap(model, gap, spec, maps, source, native, strategy, 91)
    assert calls['torsions'][1] == {'seed': 91}
    assert calls['sample'][0][1:3] == (1, 5000)
    assert calls['sample'][0][-1] == 91
    assert calls['ccd'][0][:3] == (gap.before, gap.after, gap.full_seq)
    assert calls['ccd'][0][-1] == 91
    assert result['steps'] == 5000 and result['original_stems_only'] is True
    assert result['source_boundary_conditioning'] == dict(n_stem_phi=.5, c_stem_psi=.5, prev_aa='A', next_aa='A')
    if strategy == 'torsion_mc':
        assert calls['phi_psi'][1] == dict(seed=91, n_stem_phi=.5, c_stem_psi=.5, prev_aa='A', next_aa='A')
        assert 'fragment' not in calls
    else:
        # FragmentSampler has no boundary-torsion kwargs; do not claim otherwise.
        assert calls['fragment'][1] == dict(init_fragments=5, seed=91)
        assert calls['fragger_list'][0] == (spec['start'] - 1, spec['end'] - 1)
        assert 'no phi/psi keyword API' in result['sampler_conditioning']
    methods.LoopCandidates.FillFromMonteCarloSampler = lambda *a: []
    with pytest.raises(worker.CandidateUnavailable, match='no closed candidate'):
        worker._sample_gap(model, gap, spec, maps, source, native, strategy, 91)
    methods.LoopCandidates.FillFromMonteCarloSampler = lambda *a: [Backbone(), Backbone()]
    with pytest.raises(ValueError, match='candidate count'):
        worker._sample_gap(model, gap, spec, maps, source, native, strategy, 91)
    def unclosed(*args):
        raise RuntimeError('Failed to initialize monte carlo sampling protocol!')
    methods.LoopCandidates.FillFromMonteCarloSampler = unclosed
    with pytest.raises(worker.CandidateUnavailable, match='span 2-4, seed 91'):
        worker._sample_gap(model, gap, spec, maps, source, native, strategy, 91)
    def corrupt(*args):
        raise RuntimeError('Native scorer database is corrupt')
    methods.LoopCandidates.FillFromMonteCarloSampler = corrupt
    with pytest.raises(RuntimeError, match='corrupt'):
        worker._sample_gap(model, gap, spec, maps, source, native, strategy, 91)


def test_parsed_source_precision_is_distinct_from_actual_model_movement(tmp_path):
    config = request(tmp_path)
    source = worker.read_source(config['input_pdb'], config['chains'])
    model, maps = model_for(config, source)
    specs = worker.gap_specs(config['chains'])
    native = NS(ost=NS(mol=NS(ResNum=lambda pos: pos)))
    identity = ('A', '2', '', 'ALA', 'C')
    source[identity]['xyz_angstrom'] = [97.163, 98.664, 96.095]
    parsed = {key: list(value['xyz_angstrom']) for key, value in source.items()}
    # Measured native PDB parser output from the real frozen 1HCK input.
    parsed[identity] = [97.16299438476562, 98.66399383544922, 96.0949935913086]
    target = model.model.FindResidue('A', 2).atoms[2]
    target.pos = NS(**dict(zip(('x', 'y', 'z'), parsed[identity])))
    worker._assert_source_stems(model, specs, maps, source, native, parsed)
    target.pos.x += 1e-7
    with pytest.raises(ValueError, match='raw model changed'):
        worker._assert_source_stems(model, specs, maps, source, native, parsed)
    parsed[identity][0] += .001
    target.pos.x = parsed[identity][0]
    with pytest.raises(ValueError, match='parser changed'):
        worker._assert_source_stems(model, specs, maps, source, native, parsed)


def test_native_gap_inventory_and_sequence_must_match_exact_request(tmp_path):
    config = request(tmp_path)
    specs = worker.gap_specs(config['chains'])
    spec = specs[0]
    class Gap:
        full_seq = spec['sequence']
        before = NS(number=NS(num=spec['start']))
        after = NS(number=NS(num=spec['end']))
        def Copy(self):
            return copy.copy(self)
        def IsTerminal(self):
            return False
        def GetChainIndex(self):
            return 0
    original = Gap()
    model = NS(gaps=[original], model=NS(chains=[NS(name='A')]))
    result = worker._native_gaps(model, specs)
    assert result[('A', spec['start'], spec['end'])] is not original
    original.full_seq = 'GGG'
    with pytest.raises(ValueError, match='sequence'):
        worker._native_gaps(model, specs)
    model.gaps = []
    with pytest.raises(ValueError, match='spans'):
        worker._native_gaps(model, specs)
