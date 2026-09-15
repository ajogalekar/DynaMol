#!/usr/bin/env python3
"""Private ProMod3 subprocess: output only requested missing-residue coordinates.

Run with the separate ProMod3 Python, passing input/output JSON paths. No
ProMod3/OpenStructure dependency is imported into the main application runtime.
The caller must enforce a process timeout, preserve observed atoms, and apply
fixed-heavy refinement plus independent full-environment geometry/stereo gates.
"""
import json, sys, time, os, math, re, tempfile
from pathlib import Path

_LETTERS=dict(zip('ALA ARG ASN ASP CYS GLN GLU GLY HIS ILE LEU LYS MET PHE PRO SER THR TRP TYR VAL'.split(),'ARNDCQEGHILKMFPSTWYV'))


class LoopConstructionError(RuntimeError):
    def __init__(self, message, report):
        super().__init__(message)
        self.report = report


def _write_json(path, value):
    path = Path(path)
    with tempfile.NamedTemporaryFile(mode='w', dir=path.parent, prefix=path.name+'.', delete=False) as output:
        json.dump(value, output, indent=2, allow_nan=False)
        output.write('\n')
        temporary = Path(output.name)
    temporary.replace(path)


def construct_loops(model, modelling, fdb, sdb, sampler, on_progress=lambda report: None):
    """Broaden construction context, never the caller's exported atom inventory.

    A rigid pair of crystal stems can defeat closure even for a short gap.
    Native ProMod3 normally extends those stems. We allow at most two context
    residues, then export only originally missing atoms. The caller restores
    every observed coordinate and independently refines/checks the joins.
    """
    stages = [
        ('fragment_database', 0, 'Searching loop fragments', dict(min_loops_required=4)),
        ('sparse_fragment_database', 0, 'Checking remaining fragment candidates', dict(min_loops_required=1)),
        ('context_fragment_database', 1, 'Searching with one neighboring residue of context', dict(min_loops_required=1)),
        ('context_fragment_database', 2, 'Searching with two neighboring residues of context', dict(min_loops_required=1)),
        ('monte_carlo', 2, 'Sampling conformations for remaining loops', dict(mc_num_loops=6, mc_steps=5000)),
    ]
    report = {'attempts': [], 'max_context_extension_residues': 2,
              'observed_context_changes_exported': False,
              'acceptance': 'Candidate only; original observed coordinates and final geometry gates are enforced by the caller.'}
    for method, extension, message, settings in stages:
        if not model.gaps:
            break
        attempt = dict(method=method, context_extension_residues=extension,
                       settings=settings, gaps_before=[str(g) for g in model.gaps], status='running')
        report['attempts'].append(attempt)
        report['message'] = f'{message} ({len(model.gaps)} gaps remaining)…'
        on_progress(report)
        started = time.monotonic()
        try:
            if method == 'monte_carlo':
                # This native high-level API has no seed parameter. Its sampler,
                # closer and acceptance generator use their documented defaults.
                modelling.FillLoopsByMonteCarlo(model, sampler, max_extension=extension,
                    max_loops_to_search=6, ring_punch_detection=1, **settings)
                attempt['native_default_random_seed'] = 0
            else:
                modelling.FillLoopsByDatabase(model, fdb, sdb, sampler,
                    max_res_extension=extension, max_loops_to_search=40,
                    ring_punch_detection=1, **settings)
            attempt['gaps_after'] = [str(g) for g in model.gaps]
            attempt['status'] = 'complete'
        except Exception as exc:
            attempt.update(status='failed', error=str(exc), error_type=type(exc).__name__)
            raise LoopConstructionError('Loop construction stopped: '+str(exc), report) from exc
        finally:
            attempt['elapsed_seconds'] = time.monotonic()-started
            on_progress(report)
    report['gaps_after'] = [str(g) for g in model.gaps]
    if model.gaps:
        raise LoopConstructionError(
            'No closed candidate was found within the bounded fragment and conformational search for: '
            + ', '.join(report['gaps_after'])
            + '. This does not establish that the gap is impossible to model. Supply a repaired structure; '
              'no observed coordinates were replaced and no prepared model was accepted.', report)
    return report


def validate_input(config):
    """Validate exact ordered identity/alignment metadata without native imports."""
    if not isinstance(config,dict) or not isinstance(config.get('input_pdb'),str):
        raise ValueError('An observed-coordinate PDB path is required.')
    if config.get('max_res_extension',0)!=0:
        raise ValueError('The requested missing-residue inventory cannot extend into observed residues.')
    chains=config.get('chains')
    if not isinstance(chains,list) or not 1<=len(chains)<=62:
        raise ValueError('Loop modelling requires explicit protein chain alignments.')
    seen_chains=set();missing=0
    for chain in chains:
        name=chain.get('chain_id') if isinstance(chain,dict) else None
        if not isinstance(name,str) or not re.fullmatch(r'[A-Za-z0-9]',name) or name in seen_chains:
            raise ValueError('Loop modelling needs unique canonical single-character protein chain IDs.')
        seen_chains.add(name);rows=chain.get('residues')
        if not isinstance(rows,list) or not 2<=len(rows)<=10000:
            raise ValueError('Each chain needs a bounded ordered sequence identity list.')
        seen=set();run_length=0
        for i,row in enumerate(rows):
            if not isinstance(row,dict) or row.get('chain')!=name or type(row.get('observed')) is not bool:
                raise ValueError('Each source residue requires its exact chain and observed flag.')
            fields=[row.get(field) for field in ['resid','insertion_code','residue','one_letter']]
            if any(not isinstance(value,str) for value in fields) or not fields[0] or len(fields[1])>1 or len(fields[3])!=1:
                raise ValueError('Residue identities and sequence letters are malformed.')
            key=tuple(fields[:3])
            if key in seen:raise ValueError('Source residue identities are ambiguous.')
            seen.add(key)
            if row['observed']:
                run_length=0
            else:
                if i in (0,len(rows)-1):raise ValueError('Only internal sequence-supported loops are modeled.')
                if _LETTERS.get(row['residue'])!=row['one_letter']:
                    raise ValueError('Missing residues must retain exact supported standard amino-acid identities.')
                run_length+=1;missing+=1
                if run_length>12:raise ValueError('Fragment loop modelling supports at most 12 residues per gap.')
    if not 1<=missing<=96:raise ValueError('Requested missing residues exceed the local modelling budget.')
    return config


def run(config):
    validate_input(config)
    started=time.monotonic()
    import ost, promod3
    from ost import io, seq
    from promod3 import modelling, loop
    # Preserve the native resolved-gap/stem identities in the worker log.
    ost.PushVerbosityLevel(2)
    entity=io.LoadPDB(config['input_pdb'])
    alignments=seq.AlignmentList()
    maps={}; warnings=[]
    for item in config['chains']:
        rows=item['residues']; name=item['chain_id']
        view=entity.Select('cname='+name)
        observed=[r for r in rows if r['observed']]
        actual=[(str(r.number.num),str(r.number.ins_code).strip().replace('\x00',''),r.name) for r in view.residues]
        expected=[(r['resid'],r['insertion_code'],r['residue']) for r in observed]
        if actual!=expected:raise ValueError('Observed source identities do not match input alignment.')
        if any(r.one_letter_code!=row['one_letter'] for r,row in zip(view.residues,observed)):
            raise ValueError('Observed residue sequence conflicts with the supplied alignment.')
        aln=seq.CreateAlignment(seq.CreateSequence('target',''.join(r['one_letter'] for r in rows)),
                                seq.CreateSequence('observed',''.join(r['one_letter'] if r['observed'] else '-' for r in rows)))
        aln.AttachView(1,view);alignments.append(aln)
        maps[name]={i+1:r for i,r in enumerate(rows)}
        for r in rows:
            if r['residue'] not in {'ALA','ARG','ASN','ASP','CYS','GLN','GLU','GLY','HIS','ILE','LEU','LYS','MET','PHE','PRO','SER','THR','TRP','TYR','VAL'}:
                warnings.append({'identity':[name,r['resid'],r['insertion_code'],r['residue']],
                                 'note':'ProMod3 may standardize this temporary scoring context to its parent. No observed atom is exported.'})
    m=modelling.BuildRawModel(alignments,chain_names=[x['chain_id'] for x in config['chains']],include_ligands=True,aln_preprocessing=False)
    selected_chains={item['chain_id'] for item in config['chains']}
    source_environment_atoms=sum(len(chain.atoms) for chain in entity.chains if chain.name not in selected_chains)
    imported_environment_atoms=sum(len(chain.atoms) for chain in m.model.chains if chain.name=='_')
    if source_environment_atoms>imported_environment_atoms:
        warnings.append({'note':'The fragment scoring context omits some nonprotein molecules from the input. Complete prepared-force-field refinement and full retained-environment geometry checks are required before acceptance.'})
    gaps_before=[str(g) for g in m.gaps]
    expected_missing={(name,pos) for name,rows in maps.items() for pos,r in rows.items() if not r['observed']}
    raw_missing={(name,pos) for name,rows in maps.items() for pos,r in rows.items()
                 if not m.model.FindResidue(name,ost.mol.ResNum(pos)).IsValid()}
    if raw_missing!=expected_missing:raise ValueError('ProMod3 introduced unexpected missing residues.')
    observed_context = {(c.name,r.number.num,a.name): (a.pos.x,a.pos.y,a.pos.z)
        for c in m.model.chains for r in c.residues for a in r.atoms
        if (c.name,r.number.num) not in expected_missing}
    db_start=time.monotonic(); fdb=loop.LoadFragDB();sdb=loop.LoadStructureDB();sampler=loop.LoadTorsionSamplerCoil()
    database_seconds=time.monotonic()-db_start
    db_max=fdb.MaxFragLength()
    model_start=time.monotonic()
    def progress(report):
        if config.get('_progress_path'):
            _write_json(config['_progress_path'], report)
    construction = construct_loops(m,modelling,fdb,sdb,sampler,progress)
    gaps_after=[str(g) for g in m.gaps]
    modelling.ReconstructSidechains(m.model,keep_sidechains=True,build_disulfids=False,consider_ligands=True)
    changed_context = {}
    for c in m.model.chains:
        for r in c.residues:
            for a in r.atoms:
                old = observed_context.get((c.name,r.number.num,a.name))
                if old is None:
                    continue
                displacement = math.dist(old,(a.pos.x,a.pos.y,a.pos.z))/10
                if displacement > 1e-8:
                    key = (c.name,r.number.num)
                    changed_context[key] = max(changed_context.get(key,0),displacement)
    construction['discarded_context_changes'] = [
        {'identity': [name, maps[name][pos]['resid'], maps[name][pos]['insertion_code'], maps[name][pos]['residue']]
         if name in maps and pos in maps[name] else [name,str(pos)],
         'max_heavy_atom_displacement_nm': displacement}
        for (name,pos),displacement in sorted(changed_context.items())]
    model_seconds=time.monotonic()-model_start
    atoms=[]
    for name,pos in sorted(expected_missing):
        r=maps[name][pos]; actual=m.model.FindResidue(name,ost.mol.ResNum(pos))
        if actual.name!=r['residue']:raise ValueError('Modeled residue identity differs from exact source sequence.')
        for a in actual.atoms:
            if a.element in ('H','D'):continue
            p=a.pos
            if not all(math.isfinite(value) for value in (p.x,p.y,p.z)):
                raise ValueError('Fragment modelling produced non-finite coordinates.')
            atoms.append({'identity':[name,r['resid'],r['insertion_code'],r['residue'],a.name],
                          'xyz_nm':[p.x/10,p.y/10,p.z/10]})
    return {'status':'candidate','atoms':atoms,'gaps_before':gaps_before,'gaps_after':gaps_after,
            'elapsed_seconds':time.monotonic()-started,'database_seconds':database_seconds,
            'model_seconds':model_seconds,'fragment_max_length_including_stems':db_max,
            'construction': construction,
            'context_warnings':warnings,'versions':{'promod3':promod3.__version__,'ost':ost.__version__},
            'environment_context':{'input_other_chain_atoms':source_environment_atoms,
                                   'imported_ligand_chain_atoms':imported_environment_atoms,
                                   'full_complex_scoring_claimed':False},
            'randomness':{'path':'Staged fragment search; native Monte Carlo fallback only for unresolved gaps.',
                          'requested_seed':config.get('requested_seed'),'seed_applied':False,
                          'native_monte_carlo_default_seed':0 if any(a['method']=='monte_carlo' for a in construction['attempts']) else None,
                          'note':'The high-level native APIs do not expose a requested-seed parameter. Native Monte Carlo defaults to seed 0; cross-platform bitwise reproducibility is not claimed.'},
            'method':'Staged fragment search with at most two temporary context residues, then bounded native Monte Carlo. Only originally missing atoms are exported; original observed coordinates, fixed-heavy refinement and independent geometry validation are required.'}

def main(argv=None):
    args=sys.argv[1:] if argv is None else argv
    if len(args)!=2:
        raise SystemExit("Usage: promod_loop_worker.py INPUT.json OUTPUT.json")
    source,destination=map(Path,args)
    os.environ.setdefault('PM3_OPENMM_CPU_THREADS',os.environ.get('DYNAMOL_CPU_THREADS','2'))
    try:
        if source.stat().st_size>16*1024*1024:
            raise ValueError('Loop modelling request exceeds its local input limit.')
        config=json.loads(source.read_text())
        config['_progress_path']=str(source.parent/'loop-model-progress.json')
        result=run(config)
        code=0
    except Exception as exc:
        result={'status':'failed','error':str(exc),'error_type':type(exc).__name__}
        if isinstance(exc,LoopConstructionError):
            result['construction']=exc.report
        code=1
    destination.parent.mkdir(parents=True,exist_ok=True)
    _write_json(destination,result)
    return code

if __name__=='__main__':
    raise SystemExit(main())
