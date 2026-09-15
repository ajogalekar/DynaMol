"""Single native shifted-anchor closure for each predeclared saved backbone.

No sampling, coordinate restoration, force field, or model admission occurs.
"""
import argparse,hashlib,importlib.util,json,math,time
from pathlib import Path
import ost
from ost import geom
from promod3 import loop,modelling
import promod3

def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def load(p):return json.loads(Path(p).read_text())
def save(p,x):Path(p).write_text(json.dumps(x,indent=2,allow_nan=False)+'\n')
def point(p):return [p.x,p.y,p.z]
def wrap(x):return (x+math.pi)%(2*math.pi)-math.pi

def run_one(attempt,case,output,helper):
    source_candidate=load(attempt['input_candidate']);config=load(case['input'])
    assert sha(attempt['input_candidate'])==attempt['input_candidate_sha256']
    source=helper.exact_atoms(case['source']);model,maps=helper.raw_model(config)
    span=source_candidate['settings']['span'];seed=source_candidate['settings']['seed']
    gap=helper.choose_gap(model,maps,case['modeled_residue_keys'],span)
    name=model.model.chains[gap.GetChainIndex()].name;start,end=gap.before.number.num,gap.after.number.num
    keys=[helper.residue_key(name,maps[name][i]) for i in range(start,end+1)]
    wanted=set(map(tuple,case['modeled_residue_keys']))
    assert source_candidate['source_sha256']==sha(case['source'])
    assert source_candidate['modeled_residue_keys']==case['modeled_residue_keys']
    assert keys[1:-1]==case['modeled_residue_keys']
    seq=''.join(maps[name][i]['one_letter'] for i in range(start,end+1))
    lookup={tuple(a['identity']):a for a in source_candidate['backbone_atoms']}
    assert len(lookup)==len(source_candidate['backbone_atoms'])==4*len(keys)
    bb=loop.BackboneList(seq)
    for i,key in enumerate(keys):
        xyz=[]
        for atom in ['N','CA','C','O']:
            row=lookup[tuple(key+[atom])]
            assert max(abs(row['xyz_angstrom'][j]/10-row['xyz_nm'][j]) for j in range(3))<1e-12
            xyz.append(geom.Vec3(*row['xyz_angstrom']))
        bb.Set(i,*xyz,seq[i])
    def collect(b):return helper.collect_bb(b,start,name,maps)
    def n_displacement(b):return {a:math.dist(point(getattr(b,'Get'+a)(0)),source[tuple(keys[0]+[a])]) for a in ['N','CA','C','O']}
    def joins(b):
        result=[]
        for i in range(len(b)-1):
            result.append({'first':keys[i],'second':keys[i+1],
              'C_N_A':geom.Distance(b.GetC(i),b.GetN(i+1)),
              'CA_C_N_degrees':math.degrees(geom.Angle(b.GetCA(i)-b.GetC(i),b.GetN(i+1)-b.GetC(i))),
              'O_C_N_degrees':math.degrees(geom.Angle(b.GetO(i)-b.GetC(i),b.GetN(i+1)-b.GetC(i))),
              'C_N_CA_degrees':math.degrees(geom.Angle(b.GetC(i)-b.GetN(i+1),b.GetCA(i+1)-b.GetN(i+1))),
              'omega_degrees':math.degrees(b.GetOmegaTorsion(i))})
        return result
    initial=bb.Copy();bb.ApplyTransform(bb.GetTransform(0,gap.before));aligned=bb.Copy()
    target=geom.DihedralAngle(*[geom.Vec3(*source[tuple(keys[0]+[a])]) for a in ['N','CA','C','O']])
    def carbon(b):return geom.DihedralAngle(b.GetN(0),b.GetCA(0),b.GetC(0),b.GetO(0))
    delta=wrap(target-carbon(bb));bb.SetPsiTorsion(0,bb.GetPsiTorsion(0)+delta,True)
    locked=bb.Copy();suffix=bb.Extract(1,len(bb));temp=bb.Extract(0,2).ToEntity();handles=list(temp.residues)
    assert ost.mol.InSequence(handles[0],handles[1])
    after=gap.after.GetNext();assert after.IsValid() and ost.mol.InSequence(gap.after,after)
    incoming_phi=geom.DihedralAngle(locked.GetC(0),locked.GetN(1),locked.GetCA(1),locked.GetC(1))
    settings={'mode':'shifted_anchor_ccd','seed':seed,'span':span,'native_gap':str(gap),
       'source_protocol':attempt['source_protocol'],'source_generator':source_candidate['generator'],
       'input_candidate':attempt['input_candidate'],'input_candidate_sha256':attempt['input_candidate_sha256'],
       'closure_calls':1,'native_maximum_iterations':1000,'native_c_stem_N_CA_C_rmsd_cutoff_A':.1,
       'carbonyl_rotation_sequential':True,'carbonyl_rotation_radians':delta,'temporary_phi_radians':incoming_phi,
       'suffix_sequence':suffix.GetSequence(),'native_closer':'CCDCloser',
       'source_observed_N_CA_C_O_cap_A':1,'resampling':False,'observed_coordinate_restoration':False}
    save(output/'attempt-plan.json',settings)
    diagnostics={'input':collect(initial),'N_CA_C_aligned':collect(aligned),'carbonyl_locked':collect(locked),
       'temporary_anchor_map':[{'temporary_chain':h.chain.name,'temporary_number':h.number.num,'source_identity':keys[i]} for i,h in enumerate(handles)],
       'source_n_stem_N_CA_C_O_displacement_A':{'input':n_displacement(initial),'aligned':n_displacement(aligned),'locked':n_displacement(locked)},
       'carbonyl_target_radians':target,'locked_orientation_error_radians':wrap(carbon(locked)-target),
       'joins_input':joins(initial),'joins_locked':joins(locked),
       'native_versions':{'promod3':promod3.__version__,'ost':ost.__version__},'native_closure_calls':0}
    save(output/'trial-coordinates.json',diagnostics)
    closer=modelling.CCDCloser(handles[1],gap.after,suffix.GetSequence(),loop.LoadTorsionSamplerCoil(seed=seed),seed)
    closed=suffix.Copy();diagnostics['native_closure_calls']=1
    success=closer.Close(suffix,closed)
    combined=locked.Copy();combined.ReplaceFragment(closed,1,False)
    diagnostics.update(native_closure_success=success,recombined=collect(combined),joins_recombined=joins(combined),
       prefix_displacements_during_closure_A={a:geom.Distance(getattr(locked,'Get'+a)(0),getattr(combined,'Get'+a)(0)) for a in ['N','CA','C','O']},
       virtual_anchor_N_CA_C_displacements_A={a:geom.Distance(getattr(locked,'Get'+a)(1),getattr(combined,'Get'+a)(1)) for a in ['N','CA','C']},
       incoming_phi_after_radians=geom.DihedralAngle(combined.GetC(0),combined.GetN(1),combined.GetCA(1),combined.GetC(1)))
    save(output/'trial-coordinates.json',diagnostics)
    assert all(v==0 for v in diagnostics['prefix_displacements_during_closure_A'].values())
    if not success:return {'status':'native_closure_failed','native_closure_success':False,'coordinates':str(output/'trial-coordinates.json')}
    atoms=collect(combined)
    assert all(math.isfinite(v) for a in atoms for v in a['xyz_nm'])
    context=[k for k in keys if tuple(k) not in wanted]
    displacements=[{'identity':a['identity'],'displacement_angstrom':math.dist(a['xyz_angstrom'],source[tuple(a['identity'])])} for a in atoms if tuple(a['identity'][:4]) not in wanted]
    result={'schema_version':1,'case':case['case'],'generator':'promod3_shifted_anchor_ccd','settings':settings,
        'sequence':seq,'source':case['source'],'source_sha256':sha(case['source']),
        'modeled_residue_keys':case['modeled_residue_keys'],'context_residue_keys':context,'observed_internal_residue_keys':[],
        'outer_stem_residue_keys':[keys[0],keys[-1]],'backbone_atoms':atoms,
        'observed_N_CA_C_O_displacements':displacements,'maximum_observed_backbone_displacement_angstrom':max(d['displacement_angstrom'] for d in displacements),
        'complete_source_atom_count':len(source),'source_environment_preserved_by_reference':True,
        'native_standardized_noncanonical_residues_are_not_exported':True,'native_backbone_is_before_observed_restoration':True,
        'native_closure_success':True,'native_closure_calls':1,'diagnostics':str(output/'trial-coordinates.json'),
        'app_ready':False,'physical_model_validated':False}
    result['id']=hashlib.sha256(json.dumps([case['case'],settings,atoms],sort_keys=True).encode()).hexdigest()
    save(output/'candidate.json',result)
    try:
        modelling.InsertLoopClearGaps(model,combined,gap)
        modelling.ReconstructSidechains(model.model,keep_sidechains=True,build_disulfids=False,consider_ligands=True)
        added=[];moved=[];other=[]
        for chain in model.model.chains:
            if chain.name not in maps:continue
            for residue in chain.residues:
                row=maps[chain.name][residue.number.num];key=helper.residue_key(chain.name,row)
                if residue.name!=row['residue']:
                    other.append({'identity':key,'native_parent':residue.name,'exported':False});continue
                for atom in residue.atoms:
                    if atom.element in ('H','D'):continue
                    record={'identity':key+[atom.name],'xyz_nm':[v/10 for v in point(atom.pos)]}
                    if tuple(key) in wanted:added.append(record)
                    elif tuple(record['identity']) in source:
                        displacement=math.dist(point(atom.pos),source[tuple(record['identity'])])
                        if displacement>1e-4:moved.append(dict(record,displacement_angstrom=displacement))
        result.update(modeled_heavy_atoms=added,observed_context_atoms=moved,native_parent_residue_omissions=other,sidechain_reconstruction_status='complete')
    except Exception as e:result.update(sidechain_reconstruction_status='failed',sidechain_error=str(e))
    save(output/'candidate.json',result)
    return {'status':'candidate','native_closure_success':True,'candidate':str(output/'candidate.json'),'id':result['id'],'sidechain_status':result['sidechain_reconstruction_status']}

def main():
    parser=argparse.ArgumentParser();parser.add_argument('--plan',type=Path,required=True);parser.add_argument('--case',required=True);parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args();plan=load(args.plan);case=next(c for c in plan['cases'] if c['case']==args.case)
    assert all(sha(p)==h for p,h in plan['frozen_sha256'].items())
    spec=importlib.util.spec_from_file_location('frozen_native_helper',plan['helper']);helper=importlib.util.module_from_spec(spec);spec.loader.exec_module(helper)
    rows=[];start=time.monotonic();args.output.mkdir(parents=True,exist_ok=False)
    for attempt in [a for a in plan['attempts'] if a['case']==args.case]:
        folder=args.output/attempt['name'];folder.mkdir();tick=time.monotonic()
        try:
            if time.monotonic()-start>110:raise TimeoutError('Frozen batch native work deadline reached')
            row=run_one(attempt,case,folder,helper)
        except Exception as e:row={'status':'error','error':str(e),'error_type':type(e).__name__}
        row.update(name=attempt['name'],case=args.case,elapsed_seconds=time.monotonic()-tick)
        save(folder/'attempt-result.json',row);rows.append(row);save(args.output/'result.json',rows);print(json.dumps(row),flush=True)
    assert all(sha(p)==h for p,h in plan['frozen_sha256'].items())
if __name__=='__main__':main()
