"""Bounded coherent fragment pool and preregistered KIC subset for 8K5R."""
import argparse
import hashlib
import importlib.util
import itertools
import json
import math
from pathlib import Path
import sys
import time

import ost
from ost import io, seq
from promod3 import loop, modelling
import promod3

from sample_8k5r_kic import compare, residues, point

# Exact preserved source hash, not an input structure/model approval.
EXPECTED_WORKER = '400f2e4cf5860773ebb838c64671958d0427d8e48159c5ba9b5007a7aff28975'


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def save(path, value):
    path.write_text(json.dumps(value, indent=2, allow_nan=False)+'\n')


def collect(config, original, worker, deadline):
    worker.validate_input(config)
    entity = io.LoadPDB(config['input_pdb'])
    alignments = seq.AlignmentList()
    maps = {}
    for item in config['chains']:
        name, rows = item['chain_id'], item['residues']
        view = entity.Select('cname='+name)
        observed = [r for r in rows if r['observed']]
        identities = [(str(r.number.num), str(r.number.ins_code).strip().replace('\x00',''), r.name) for r in view.residues]
        if identities != [(r['resid'], r['insertion_code'], r['residue']) for r in observed]:
            raise ValueError('Frozen alignment and observed identity mapping disagree')
        aln = seq.CreateAlignment(seq.CreateSequence('target', ''.join(r['one_letter'] for r in rows)),
                                  seq.CreateSequence('observed', ''.join(r['one_letter'] if r['observed'] else '-' for r in rows)))
        aln.AttachView(1, view)
        alignments.append(aln)
        maps[name] = {i+1:r for i,r in enumerate(rows)}
    model = modelling.BuildRawModel(alignments, chain_names=[r['chain_id'] for r in config['chains']],
                                    include_ligands=True, aln_preprocessing=False)
    wanted = {('A','177','','ALA'), ('A','178','','LYS'), ('A','179','','ASN')}
    matching = []
    for gap in model.gaps:
        if gap.IsTerminal():
            continue
        name = model.model.chains[gap.GetChainIndex()].name
        keys = {(name, maps[name][i]['resid'], maps[name][i]['insertion_code'], maps[name][i]['residue'])
                for i in range(gap.before.number.num+1, gap.after.number.num)}
        if keys == wanted:
            matching.append(gap.Copy())
    if len(matching) != 1:
        raise ValueError('Exactly one native gap must map to source A177–179 AKN')
    initial = matching[0]
    current = initial.Copy()
    fdb, sdb = loop.LoadFragDB(), loop.LoadStructureDB()
    maximum = min(initial.length+2, fdb.MaxFragLength()-2)
    extender = modelling.FullGapExtender(current, model.seqres[initial.GetChainIndex()], maximum)
    insertions = modelling.CountEnclosedInsertions(model, initial)
    contexts, pool = [], []
    first = True
    while current.length <= maximum and len(pool) < 40:
        if time.monotonic() > deadline:
            raise TimeoutError('Pool enumeration exceeded its internal deadline')
        if not first and not extender.Extend():
            break
        first = False
        if not current.before.IsValid() or not current.after.IsValid():
            break
        if modelling.CountEnclosedInsertions(model, current) != insertions:
            continue
        name = model.model.chains[current.GetChainIndex()].name
        if name != 'A':
            raise ValueError('Unexpected target chain')
        native = modelling.LoopCandidates.FillFromDatabase(current.before, current.after, current.full_seq, fdb, sdb, True)
        count = min(len(native), 8, 40-len(pool))
        context = {'native_gap':str(current), 'raw_database_candidate_count':len(native),
                   'retained_first_native_candidates':count, 'candidate_ids':[],
                   'source_start':maps[name][current.before.number.num]['resid'],
                   'source_end':maps[name][current.after.number.num]['resid']}
        for ordinal, bb in enumerate(native):
            if ordinal >= count:
                break
            bb = bb.Copy()
            start = current.before.number.num
            atoms = []
            for i in range(len(bb)):
                row = maps[name][start+i]
                for atom in ('N','CA','C','O'):
                    atoms.append({'identity':[name,row['resid'],row['insertion_code'],row['residue'],atom],
                                  'xyz_angstrom':point(getattr(bb,'Get'+atom)(i))})
            if not all(math.isfinite(v) for a in atoms for v in a['xyz_angstrom']):
                raise ValueError('Native fragment contains nonfinite coordinates')
            geometry_hash = hashlib.sha256(json.dumps(atoms,sort_keys=True).encode()).hexdigest()
            identifier = hashlib.sha256(json.dumps([str(current),ordinal,atoms],sort_keys=True).encode()).hexdigest()
            displacement = []
            for a in atoms:
                ident = a['identity']
                if tuple(ident[:4]) in wanted:
                    continue
                residue = original[int(ident[1])]
                if residue.name != ident[3] or not residue.FindAtom(ident[4]).IsValid():
                    raise ValueError('Original full-anchor identity missing')
                d = math.dist(a['xyz_angstrom'],point(residue.FindAtom(ident[4]).pos))
                displacement.append({'identity':ident,'displacement_angstrom':d})
            observed_inside = {tuple(a['identity'][:4]) for a in atoms
                               if int(context['source_start']) < int(a['identity'][1]) < int(context['source_end'])
                               and tuple(a['identity'][:4]) not in wanted}
            if len(observed_inside) > 2:
                raise ValueError('Fragment exceeds the two-observed-internal-residue limit')
            record = {'id':identifier,'geometry_sha256':geometry_hash,'native_context':str(current),
                      'native_ordinal':ordinal,'source_start':int(context['source_start']),
                      'source_end':int(context['source_end']), 'backbone_atoms':atoms,
                      'observed_internal_residues':sorted(observed_inside),
                      'raw_observed_N_CA_C_O_displacements':displacement,
                      'raw_maximum_observed_backbone_displacement_angstrom':max(x['displacement_angstrom'] for x in displacement),
                      'raw_sum_squared_observed_backbone_displacement_angstrom2':sum(x['displacement_angstrom']**2 for x in displacement),
                      'raw_within_original_1A_cap':all(x['displacement_angstrom'] <= 1 for x in displacement),
                      'experimental_fragment_accession':'Not exposed by this native LoopCandidates API; coordinate hash and native ordinal retained.',
                      'no_observed_coordinate_restoration':True}
            pool.append(record)
            context['candidate_ids'].append(identifier)
        contexts.append(context)
    return {'native_target_gap':str(initial),'contexts':contexts,'candidates':pool,
            'candidate_count':len(pool),'unique_coordinate_sets':len({p['geometry_sha256'] for p in pool}),
            'raw_coordinates_before_CCD_KIC_or_restoration':True,
            'pool_reaches_overall_cap':len(pool)==40}


def export_overlay(base, candidate, target):
    requested = {tuple(a['identity']):a['xyz_angstrom'] for a in candidate['backbone_atoms']}
    found = set()
    lines = []
    maximum_rounding = 0
    for line in base.read_text().splitlines(keepends=True):
        if line.startswith(('ATOM  ','HETATM')):
            ident = (line[21].strip(),str(int(line[22:26])),line[26].strip(),line[17:20].strip(),line[12:16].strip())
            if ident in requested:
                if ident in found:
                    raise ValueError('Duplicate seed backbone atom identity')
                found.add(ident)
                xyz = requested[ident]
                formatted = ''.join(f'{v:8.3f}' for v in xyz)
                serialized = [float(formatted[i:i+8]) for i in range(0,24,8)]
                maximum_rounding = max(maximum_rounding,max(abs(a-b) for a,b in zip(xyz,serialized)))
                line = line[:30]+formatted+line[54:]
        lines.append(line)
    if found != set(requested) or maximum_rounding > .0005+1e-12:
        raise ValueError('Native candidate PDB mapping or rounding check failed')
    target.write_text(''.join(lines))
    return {'seed_sha256':sha(target),'maximum_component_rounding_angstrom':maximum_rounding,
            'replaced_backbone_atoms':len(found),
            'scope':'Only the candidate span is used for native KIC. Other seed coordinates supply topology context; original observed coordinates outside the span are restored for independent screening.'}


def main():
    parser = argparse.ArgumentParser()
    for name in ['input','worker','source','base-seed','control','output']:
        parser.add_argument('--'+name,type=Path,required=True)
    args = parser.parse_args()
    args.output.mkdir(exist_ok=False)
    config = json.loads(args.input.read_text())
    context_path = Path(config['input_pdb'])
    helper = Path(__file__).with_name('sample_8k5r_kic.py')
    paths = [args.input,args.worker,args.source,args.base_seed,args.control,context_path,Path(__file__),helper]
    hashes = {str(p):sha(p) for p in paths}
    if sha(args.worker) != EXPECTED_WORKER:
        raise ValueError('Preserved generic worker hash mismatch')
    for name,p in [('implementation.py',Path(__file__)),('kic-implementation.py',helper),('preserved-worker.py',args.worker)]:
        (args.output/name).write_bytes(p.read_bytes())
    save(args.output/'preregistered-plan.json',{
        'source_sha256':hashes,'target_source_residues':[['A','177','','ALA'],['A','178','','LYS'],['A','179','','ASN']],
        'maximum_observed_internal_context_residues':2,'maximum_raw_pool':40,
        'per_context_retention':'First up to 8 native database candidates, no outcome-based substitution.',
        'contexts':'Native FullGapExtender enumerates the original span and spans with at most two observed internal residues.',
        'closure_subset':'Up to 2 per enumerated context with lowest raw complete observed N/CA/C/O squared displacement; at most 12 total. Freeze IDs before closure/reference outcomes.',
        'closure_budget':'Up to 4 internal pivot sets per selected fragment, each unperturbed plus one independently sampled nonpivot conformer.',
        'sampling_seed':2026,'source_observed_backbone_cap_angstrom':1.0,
        'internal_deadline_seconds':160,'app_ready':False,'physical_model_validated':False})
    spec = importlib.util.spec_from_file_location('preserved_fragment_input_validator',args.worker)
    worker = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(worker)
    source_entity, original = residues(args.source)
    started = time.monotonic()
    deadline = started+160
    pool = collect(config,original,worker,deadline)
    save(args.output/'raw-pool.json',pool)
    selected = []
    for ctx in pool['contexts']:
        eligible = [c for c in pool['candidates'] if c['id'] in ctx['candidate_ids']]
        selected.extend(sorted(eligible,key=lambda c:(c['raw_sum_squared_observed_backbone_displacement_angstrom2'],c['id']))[:2])
    if len(selected) > 12:
        raise ValueError('Native context enumeration exceeded the frozen subset budget')
    plans = []
    for candidate in selected:
        n = candidate['source_end']-candidate['source_start']+1
        all_pivots = list(itertools.combinations(range(1,n-1),3))
        pivots = all_pivots if len(all_pivots)<=4 else [(1,3,5),(2,3,4),(1,2,5),(1,4,5)]
        if n > 7:
            raise ValueError('Source span unexpectedly exceeds two observed internal context residues')
        plans.append({'candidate_id':candidate['id'],'start':candidate['source_start'],'end':candidate['source_end'],
                      'internal_pivots':pivots,'nonpivot_conformer_draws_per_pivot':1 if n>5 else 0})
    save(args.output/'frozen-closure-plan.json',{'raw_pool_sha256':sha(args.output/'raw-pool.json'),'selected':plans,
                                               'selection_used_final_rama_or_geometry_outcomes':False})
    sampler = loop.LoadTorsionSamplerCoil(seed=2026)
    control_entity, control = residues(args.control)
    control_result = compare(control,control,1,6,set(),0,[(1,2,3),(1,2,4),(1,3,4),(2,3,4)],sampler,deadline)
    control_result['name'] = 'intact-control'
    cases = [control_result]
    for index,(candidate,plan) in enumerate(zip(selected,plans)):
        name = f'fragment-{index:02d}-{candidate["id"][:12]}'
        folder = args.output/name
        folder.mkdir()
        seed_path = folder/'coherent-BACKBONE-overlay.pdb'
        export = export_overlay(args.base_seed,candidate,seed_path)
        save(folder/'seed-export.json',export)
        seed_entity, seeded = residues(seed_path)
        result = compare(original,seeded,plan['start'],plan['end'],{177,178,179},plan['nonpivot_conformer_draws_per_pivot'],plan['internal_pivots'],sampler,deadline)
        result.update(name=name,candidate_id=candidate['id'],seed_path=str(seed_path),
                      raw_candidate= candidate,seed_export=export)
        save(folder/'result.json',result)
        cases.append(result)
        print(json.dumps({'case':name,'span':[plan['start'],plan['end']],
                          'raw_anchor_maximum_angstrom':candidate['raw_maximum_observed_backbone_displacement_angstrom'],
                          'solutions':result['raw_solutions'],'cap_and_omega_pass':result['solutions_within_observed_cap_and_omega_basins']}),flush=True)
    if hashes != {str(p):sha(p) for p in paths}:
        raise ValueError('Source or helper changed during the native comparison')
    save(args.output/'result.json',{'cases':cases,'pool':pool,'source_sha256':hashes,
        'frozen_closure_plan_sha256':sha(args.output/'frozen-closure-plan.json'),
        'promod3_version':promod3.__version__,'ost_version':ost.__version__,
        'elapsed_seconds':time.monotonic()-started,'app_ready':False,'physical_model_validated':False,
        'scope':'Multiple coherent raw fragments and bounded closure. No source-restored joins before KIC; no full-atom candidate is accepted.'})


if __name__ == '__main__':
    main()
