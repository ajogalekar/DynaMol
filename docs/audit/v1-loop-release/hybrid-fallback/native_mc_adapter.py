"""Bounded native research generator. Each attempt starts from original gaps.

Raw backbone and coherent context are exported before any restoration. This
script never accepts a prepared model and never changes force-field chemistry.
"""
import argparse
import hashlib
import json
import math
from pathlib import Path
import time

import ost
from ost import io, seq, geom
from promod3 import loop, modelling
import promod3


def sha(p):
    return hashlib.sha256(Path(p).read_bytes()).hexdigest()


def save(path, data):
    path.write_text(json.dumps(data, indent=2, allow_nan=False)+'\n')


def point(v):
    return [v.x, v.y, v.z]


def exact_atoms(path):
    result = {}
    for line in Path(path).read_text().splitlines():
        if line[:6].strip() not in ('ATOM', 'HETATM'):
            continue
        key = (line[21].strip(), line[22:26].strip(), line[26].strip(), line[17:20].strip(), line[12:16].strip())
        if key in result:
            raise ValueError('Ambiguous original atom '+str(key))
        result[key] = [float(line[30:38]), float(line[38:46]), float(line[46:54])]
    return result


def raw_model(config):
    entity = io.LoadPDB(config['input_pdb'])
    alignments = seq.AlignmentList()
    maps = {}
    for chain in config['chains']:
        name, rows = chain['chain_id'], chain['residues']
        view = entity.Select('cname='+name)
        expected = [(r['resid'], r['insertion_code'], r['residue']) for r in rows if r['observed']]
        actual = [(str(r.number.num), str(r.number.ins_code).strip().replace('\x00',''), r.name) for r in view.residues]
        if actual != expected:
            raise ValueError('Exact observed source/alignment mapping differs')
        if [r.one_letter_code for r in view.residues] != [r['one_letter'] for r in rows if r['observed']]:
            raise ValueError('Observed one-letter sequence differs')
        aln = seq.CreateAlignment(seq.CreateSequence('target',''.join(r['one_letter'] for r in rows)),
                                  seq.CreateSequence('observed',''.join(r['one_letter'] if r['observed'] else '-' for r in rows)))
        aln.AttachView(1, view)
        alignments.append(aln)
        maps[name] = {i+1:r for i,r in enumerate(rows)}
    model = modelling.BuildRawModel(alignments, chain_names=[c['chain_id'] for c in config['chains']],
                                    include_ligands=True, aln_preprocessing=False)
    expected = {(name,i) for name,rows in maps.items() for i,r in rows.items() if not r['observed']}
    absent = {(name,i) for name,rows in maps.items() for i in rows if not model.model.FindResidue(name,ost.mol.ResNum(i)).IsValid()}
    if absent != expected:
        raise ValueError('Native raw model changed requested missing inventory')
    return model, maps


def residue_key(name, row):
    return [name,row['resid'],row['insertion_code'],row['residue']]


def choose_gap(model, maps, target, span):
    matches = []
    for gap in model.gaps:
        if gap.IsTerminal():
            continue
        name = model.model.chains[gap.GetChainIndex()].name
        keys = {tuple(residue_key(name,maps[name][i])) for i in range(gap.before.number.num+1,gap.after.number.num)}
        if keys == set(map(tuple,target)):
            matches.append(gap.Copy())
    if len(matches) != 1:
        raise ValueError('Target missing sequence does not identify exactly one native gap')
    original = matches[0]
    current = original.Copy()
    extender = modelling.FullGapExtender(current,model.seqres[current.GetChainIndex()],original.length+2)
    first = True
    while current.length <= original.length+2:
        if not first and not extender.Extend():
            break
        first = False
        name = model.model.chains[current.GetChainIndex()].name
        if [maps[name][current.before.number.num]['resid'],maps[name][current.after.number.num]['resid']] == [str(x) for x in span]:
            if modelling.CountEnclosedInsertions(model,current) != modelling.CountEnclosedInsertions(model,original):
                raise ValueError('Selected context would merge missing gaps')
            return current.Copy()
    raise ValueError('Frozen context span is not available from native bounded extender')


def collect_bb(bb, start, name, maps):
    return [{'identity':residue_key(name,maps[name][start+i])+[atom],
             'xyz_angstrom':point(getattr(bb,'Get'+atom)(i)),
             'xyz_nm':[v/10 for v in point(getattr(bb,'Get'+atom)(i))]}
            for i in range(len(bb)) for atom in ('N','CA','C','O')]


def source_boundaries(name, maps, start, end, source):
    """Condition on observed flanking torsions, with explicit terminal defaults.

    Missing atoms in an existing observed neighbor are integrity errors. Native
    defaults are retained only where a terminal torsion has no source value.
    """
    rows = maps[name]
    values = {'n_stem_phi': -1.0472, 'c_stem_psi': -.785398,
              'prev_aa': 'A', 'next_aa': 'A'}
    unavailable = []
    def atom(number, atom_name):
        if not rows[number]['observed']:
            raise ValueError('Boundary conditioning crosses an unresolved gap')
        return geom.Vec3(*source[tuple(residue_key(name, rows[number]) + [atom_name])])
    if start - 1 in rows:
        values['n_stem_phi'] = geom.DihedralAngle(atom(start - 1, 'C'), atom(start, 'N'), atom(start, 'CA'), atom(start, 'C'))
        values['prev_aa'] = rows[start - 1]['one_letter']
    else:
        unavailable.append('N-terminal stem has no observed preceding phi; explicit native default retained')
    if end + 1 in rows:
        values['c_stem_psi'] = geom.DihedralAngle(atom(end, 'N'), atom(end, 'CA'), atom(end, 'C'), atom(end + 1, 'N'))
        values['next_aa'] = rows[end + 1]['one_letter']
    else:
        unavailable.append('C-terminal stem has no observed following psi; explicit native default retained')
    if not all(math.isfinite(values[k]) for k in ('n_stem_phi', 'c_stem_psi')):
        raise ValueError('Nonfinite source boundary torsion')
    return values, unavailable


def native_attempt(case, mode, span, seed, output, deadline, boundary_protocol='native_default'):
    config = json.loads(Path(case['input']).read_text())
    original_atoms = exact_atoms(case['source'])
    model, maps = raw_model(config)  # No earlier candidate or cleared gap survives.
    gap = choose_gap(model,maps,case['modeled_residue_keys'],span)
    name = model.model.chains[gap.GetChainIndex()].name
    start, end = gap.before.number.num, gap.after.number.num
    wanted = set(map(tuple,case['modeled_residue_keys']))
    context = [residue_key(name,maps[name][i]) for i in range(start,end+1) if tuple(residue_key(name,maps[name][i])) not in wanted]
    inside = [k for k in context if k[1] not in (str(span[0]),str(span[1]))]
    if len(inside)>2 or any(maps[name][i]['residue'] not in {'ALA','ARG','ASN','ASP','CYS','GLN','GLU','GLY','HIS','ILE','LEU','LYS','MET','PHE','PRO','SER','THR','TRP','TYR','VAL'} for i in range(start,end+1)):
        raise ValueError('Candidate span exceeds canonical bounded context')
    modelling.SetupDefaultBackboneScoring(model)
    torsions = loop.LoadTorsionSamplerCoil(seed=seed)
    proposal = {'mode':mode,'seed':seed,'steps':5000,'native_gap':str(gap),'span':span,
                'original_model_reset':True,'closer':'DirtyCCDCloser','scoring':['reduced','cb_packing','clash']}
    if boundary_protocol not in {'native_default', 'source_dirty', 'source_ccd'}:
        raise ValueError('Unknown bounded boundary-conditioning protocol')
    proposal['boundary_protocol'] = boundary_protocol
    conditioning = {}
    if boundary_protocol != 'native_default':
        conditioning, unavailable = source_boundaries(name, maps, start, end, original_atoms)
        proposal['source_boundary_conditioning'] = dict(conditioning, terminal_defaults=unavailable,
            n_stem_phi_degrees=math.degrees(conditioning['n_stem_phi']),
            c_stem_psi_degrees=math.degrees(conditioning['c_stem_psi']))
    if mode == 'torsion_mc':
        sampler = modelling.PhiPsiSampler(gap.full_seq,torsions,seed=seed,**conditioning)
        proposal['sampler'] = 'PhiPsiSampler'
        proposal['sampler_boundary_conditioning'] = 'exact observed source' if conditioning else 'native defaults'
    else:
        handler = modelling.FraggerHandle(model.seqres[gap.GetChainIndex()], fragment_length=case['fragment_length'],
                                          fragments_per_position=100, structure_db=loop.LoadStructureDB())
        fraggers = handler.GetList(start-1,end-1)
        sampler = modelling.FragmentSampler(gap.full_seq,fraggers,init_fragments=5,seed=seed)
        proposal.update(sampler='FragmentSampler',fragment_length=case['fragment_length'],
                        fragments_per_position=100,fragger_count=len(fraggers),init_fragments=5,
                        ranking_input='Native BLOSUM62 sequence-only fragment search; no hidden control coordinates')
    if boundary_protocol == 'source_ccd':
        closer = modelling.CCDCloser(gap.before, gap.after, gap.full_seq, torsions, seed)
        proposal['closer'] = 'CCDCloser with native conditional torsion checks'
    else:
        closer = modelling.DirtyCCDCloser(gap.before,gap.after)
    weights = modelling.ScoringWeights.GetWeights()
    selected_weights = {k:weights[k] for k in ('reduced','cb_packing','clash')}
    scorer = modelling.LinearScorer(model.backbone_scorer,model.backbone_scorer_env,start,len(gap.full_seq),gap.GetChainIndex(),selected_weights)
    cooler = modelling.ExponentialCooler(round(5000/109),100,.9)
    save(output/'attempt-plan.json',proposal)
    if time.monotonic()>deadline:
        raise TimeoutError('Deadline reached before native MC')
    candidates = modelling.LoopCandidates.FillFromMonteCarloSampler(gap.full_seq,1,5000,sampler,closer,scorer,cooler,seed)
    if len(candidates)!=1:
        raise ValueError('Expected exactly one bounded native result')
    bb = candidates[0].Copy()
    atoms = collect_bb(bb,start,name,maps)
    if not all(math.isfinite(v) for a in atoms for v in a['xyz_nm']):
        raise ValueError('Nonfinite native coordinates')
    displacements = [{'identity':a['identity'],'displacement_angstrom':math.dist(a['xyz_angstrom'],original_atoms[tuple(a['identity'])])}
                     for a in atoms if tuple(a['identity'][:4]) not in wanted]
    baseline_links=[]
    for i in range(start,end):
        left,right=residue_key(name,maps[name][i]),residue_key(name,maps[name][i+1])
        if tuple(left) in wanted or tuple(right) in wanted:
            continue
        points=[original_atoms[tuple(k+[a])] for k,a in [(left,'CA'),(left,'C'),(right,'N'),(right,'CA')]]
        omega=geom.DihedralAngle(*[geom.Vec3(*p) for p in points])
        actual=bb.GetOmegaTorsion(i-start)
        baseline_links.append({'first':left,'second':right,'source_radians':omega,'candidate_radians':actual,
                               'basin_preserved':(-math.pi/2<omega<math.pi/2)==(-math.pi/2<actual<math.pi/2)})
    result={'schema_version':1,'case':case['case'],'generator':'promod3_'+mode,'settings':proposal,
            'sequence':gap.full_seq,'source':case['source'],'source_sha256':sha(case['source']),
            'modeled_residue_keys':case['modeled_residue_keys'],'context_residue_keys':context,
            'observed_internal_residue_keys':inside,'backbone_atoms':atoms,
            'outer_stem_residue_keys':[residue_key(name,maps[name][start]),residue_key(name,maps[name][end])],
            'observed_N_CA_C_O_displacements':displacements,
            'maximum_observed_backbone_displacement_angstrom':max(a['displacement_angstrom'] for a in displacements),
            'observed_omega_checks':baseline_links,'complete_source_atom_count':len(original_atoms),
            'source_environment_preserved_by_reference':True,
            'native_standardized_noncanonical_residues_are_not_exported':True,
            'native_backbone_is_before_observed_restoration':True,
            'native_backbone_score':scorer.GetScore(bb),'app_ready':False,'physical_model_validated':False}
    result['id']=hashlib.sha256(json.dumps([case['case'],proposal,atoms],sort_keys=True).encode()).hexdigest()
    save(output/'candidate.json',result)  # Preserve raw result even if sidechains fail.
    try:
        modelling.InsertLoopClearGaps(model,bb,gap)
        modelling.ReconstructSidechains(model.model,keep_sidechains=True,build_disulfids=False,consider_ligands=True)
        added=[]; moved=[]; other=[]
        for chain in model.model.chains:
            if chain.name not in maps:
                continue
            for residue in chain.residues:
                row=maps[chain.name][residue.number.num]
                key=residue_key(chain.name,row)
                if residue.name!=row['residue']:
                    other.append({'identity':key,'native_parent':residue.name,'exported':False})
                    continue
                for atom in residue.atoms:
                    if atom.element in ('H','D'):
                        continue
                    record={'identity':key+[atom.name],'xyz_nm':[v/10 for v in point(atom.pos)]}
                    if tuple(key) in wanted:
                        added.append(record)
                    elif tuple(record['identity']) in original_atoms:
                        displacement=math.dist(point(atom.pos),original_atoms[tuple(record['identity'])])
                        if displacement>1e-4:
                            moved.append(dict(record,displacement_angstrom=displacement))
        result.update(modeled_heavy_atoms=added,observed_context_atoms=moved,
                      native_parent_residue_omissions=other,sidechain_reconstruction_status='complete')
    except Exception as exc:
        result.update(sidechain_reconstruction_status='failed',sidechain_error=str(exc))
    save(output/'candidate.json',result)
    return result


def main():
    p=argparse.ArgumentParser();p.add_argument('--plan',type=Path,required=True);p.add_argument('--case',required=True)
    p.add_argument('--mode',choices=['torsion_mc','fragment_mc'],required=True);p.add_argument('--output',type=Path,required=True)
    p.add_argument('--boundary-protocol', choices=['native_default', 'source_dirty', 'source_ccd'], default='native_default')
    args=p.parse_args(); args.output.mkdir(parents=True,exist_ok=False)
    plan=json.loads(args.plan.read_text());case=next(c for c in plan['cases'] if c['case']==args.case)
    files=[args.plan,Path(__file__),Path(case['source']),Path(case['input']),Path(json.loads(Path(case['input']).read_text())['input_pdb'])]
    hashes={str(f):sha(f) for f in files};save(args.output/'provenance.json',{'source_sha256':hashes,'promod3':promod3.__version__,'ost':ost.__version__})
    started=time.monotonic();deadline=started+165;records=[]
    for span in case['spans']:
        for seed in case.get('seeds', plan['seeds']):
            folder=args.output/f'{span[0]}-{span[1]}-seed-{seed}';folder.mkdir()
            tick=time.monotonic()
            try:
                if time.monotonic()>deadline:raise TimeoutError('Batch deadline exhausted')
                candidate=native_attempt(case,args.mode,span,seed,folder,deadline,args.boundary_protocol)
                record={'span':span,'seed':seed,'status':'candidate','candidate':str(folder/'candidate.json'),'id':candidate['id']}
            except Exception as exc:
                record={'span':span,'seed':seed,'status':'failed','error':str(exc),'error_type':type(exc).__name__}
            record['elapsed_seconds']=time.monotonic()-tick;records.append(record)
            save(folder/'attempt-result.json',record);print(json.dumps(record),flush=True)
            save(args.output/'result.json',{'case':args.case,'mode':args.mode,'attempts':records,'app_ready':False})
    if hashes!={str(f):sha(f) for f in files}:raise ValueError('Frozen input changed')


if __name__=='__main__':main()
