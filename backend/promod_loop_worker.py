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


def validate_input(config):
    """Validate exact ordered identity/alignment metadata without native imports."""
    if not isinstance(config,dict) or not isinstance(config.get('input_pdb'),str):
        raise ValueError('An observed-coordinate PDB path is required.')
    if config.get('max_res_extension',0)!=0:
        raise ValueError('Loop modelling cannot extend into observed residues.')
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
    db_start=time.monotonic(); fdb=loop.LoadFragDB();sdb=loop.LoadStructureDB();sampler=loop.LoadTorsionSamplerCoil()
    database_seconds=time.monotonic()-db_start
    db_max=fdb.MaxFragLength()
    model_start=time.monotonic()
    modelling.FillLoopsByDatabase(m,fdb,sdb,sampler,max_res_extension=0,max_loops_to_search=40,min_loops_required=4)
    gaps_after=[str(g) for g in m.gaps]
    if gaps_after:raise RuntimeError('Fragment database did not close the requested gap: '+str(gaps_after))
    modelling.ReconstructSidechains(m.model,keep_sidechains=True,build_disulfids=False,consider_ligands=True)
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
            'context_warnings':warnings,'versions':{'promod3':promod3.__version__,'ost':ost.__version__},
            'environment_context':{'input_other_chain_atoms':source_environment_atoms,
                                   'imported_ligand_chain_atoms':imported_environment_atoms,
                                   'full_complex_scoring_claimed':False},
            'randomness':{'path':'Fragment database search and rotamer scoring; no Monte Carlo loop fallback.',
                          'requested_seed':config.get('seed'),'seed_applied':False,
                          'note':'This worker does not invoke a native random-seed API; it does not claim bitwise reproducibility across platforms.'},
            'method':'Fragment database only, extension0, keep complete sidechains, no Monte Carlo/global minimization. Candidate requires independent geometry validation.'}

def main(argv=None):
    args=sys.argv[1:] if argv is None else argv
    if len(args)!=2:
        raise SystemExit("Usage: promod_loop_worker.py INPUT.json OUTPUT.json")
    source,destination=map(Path,args)
    os.environ.setdefault('PM3_OPENMM_CPU_THREADS',os.environ.get('DYNAMOL_CPU_THREADS','2'))
    try:
        if source.stat().st_size>16*1024*1024:
            raise ValueError('Loop modelling request exceeds its local input limit.')
        result=run(json.loads(source.read_text()))
        code=0
    except Exception as exc:
        result={'status':'failed','error':str(exc),'error_type':type(exc).__name__}
        code=1
    destination.parent.mkdir(parents=True,exist_ok=True)
    with tempfile.NamedTemporaryFile(mode='w',dir=destination.parent,prefix=destination.name+'.',delete=False) as output:
        json.dump(result,output,indent=2,allow_nan=False);output.write('\n');temporary=Path(output.name)
    temporary.replace(destination)
    return code

if __name__=='__main__':
    raise SystemExit(main())
