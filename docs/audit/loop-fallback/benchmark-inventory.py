"""Fetch actual experimental benign-host candidate PDBs and inventory internal gaps.

This pre-preparation inventory screens entry sequence/coordinate availability, not
native modeling success. Sources and all screened-out entries remain recorded.
"""
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor
import hashlib,json,time,urllib.request,datetime,sys
import gemmi
OUT=Path(__file__).resolve().parent
CANDIDATES=['8K5R','1UA2','1ATP','1HCK','1FIN','2SRC','1Y57','4HJO','2ITY','1M17','1IEP','2HYY','2CG9','1T46','1QCF','1PTY','1BZH','1OVA','1AKE','1HVY','1U72','1DHF','3MXF','4I23','1A52','2REN','1RNE','3HEG','6DBK','6A93','1YET','1YER','1BYQ','1AM1','1A4H','1AH3','1XU7','2BEL','3CPR','1GOS','1A9U','1AQ1','1FPU','2VV5','4HHY','3K8V']
def rows(block,category):
 table=block.find_mmcif_category('_'+category+'.')
 if not table:return []
 names=[str(tag).split('.',1)[1] for tag in table.tags]
 return [dict(zip(names,list(row))) for row in table]
def clean(v):
 return None if v in (None,'.','?') else gemmi.cif.as_string(v)
def one(pdb):
 begin=time.monotonic();p=OUT/'benchmark-inputs'/f'{pdb}.cif';url=f'https://files.rcsb.org/download/{pdb}.cif'
 r={'id':pdb,'source_url':url,'input_path':str(p),'retrieved_at':datetime.datetime.now(datetime.timezone.utc).isoformat()}
 try:
  if not p.exists():
   with urllib.request.urlopen(url,timeout=30) as reply:data=reply.read(15_000_000)
   if len(data)>=15_000_000:raise ValueError('Source exceeded inventory byte limit')
   p.write_bytes(data)
  r.update(sha256=hashlib.sha256(p.read_bytes()).hexdigest(),bytes=p.stat().st_size)
  b=gemmi.cif.read_file(str(p)).sole_block();r['title']=clean(b.find_value('_struct.title'))
  r['organisms']=sorted(set(clean(x.get(k)) for cat,k in [('entity_src_gen','pdbx_gene_src_scientific_name'),('entity_src_nat','pdbx_organism_scientific')] for x in rows(b,cat) if clean(x.get(k))))
  r['resolution_angstrom']=clean(b.find_value('_refine.ls_d_res_high'))
  entities={x['id']:clean(x.get('pdbx_description')) for x in rows(b,'entity')}
  polymer={x['entity_id']:x['type'] for x in rows(b,'entity_poly')}
  chains={x['id']:x['entity_id'] for x in rows(b,'struct_asym')};sites=rows(b,'atom_site');observed={};ligands={}
  for a in sites:
   if a.get('pdbx_PDB_model_num','1')!='1':continue
   label=a['label_asym_id'];seq=clean(a.get('label_seq_id'))
   if seq:observed.setdefault(label,set()).add(int(seq))
   else:ligands.setdefault(a['auth_asym_id'],set()).add(a['label_comp_id'])
  bychain={}
  for s in rows(b,'pdbx_poly_seq_scheme'):bychain.setdefault(s['asym_id'],[]).append(s)
  result=[]
  for label,scheme in bychain.items():
   entity=chains[label]
   if clean(polymer.get(entity))!='polypeptide(L)' or not observed.get(label):continue
   seen=observed[label];lo,hi=min(seen),max(seen);gaps=[];current=[]
   for s in sorted(scheme,key=lambda x:int(x['seq_id'])):
    n=int(s['seq_id'])
    if lo<n<hi and n not in seen:current.append(s)
    elif current:gaps.append(current);current=[]
   if current:gaps.append(current)
   gaprows=[{'length':len(g),'label_sequence_ids':[int(s['seq_id']) for s in g],'residues':[s['mon_id'] for s in g], 'author_residue_ids':[clean(s.get('auth_seq_num')) or clean(s.get('pdb_seq_num')) for s in g]} for g in gaps]
   author=scheme[0]['pdb_strand_id'];result.append({'label_chain':label,'author_chain':author,'entity_id':entity,'protein':entities.get(entity),'observed_residues':len(seen),'sequence_length':len(scheme),'internal_gaps':gaprows,'internal_missing_total':sum(g['length'] for g in gaprows),'ligand_codes_by_author_chain':sorted(ligands.get(author,set()))})
  r.update(chains=result,status='inventoried',elapsed_seconds=round(time.monotonic()-begin,2))
 except Exception as exc:r.update(status='error',error=str(exc))
 return r
if __name__=='__main__':
 with ThreadPoolExecutor(max_workers=6) as pool:
  rr=list(pool.map(one,CANDIDATES))
 (OUT/'benchmark-screening-inventory.json').write_text(json.dumps({'scope':'Pre-preparation sequence-gap inventory, actual source files; no modeling outcomes consulted for selection.','candidates':rr},indent=2)+'\n')
 for r in rr:
  print(r['id'],r.get('organisms'),r.get('title'),r.get('error',''))
  for c in r.get('chains',[]):
   print(' ',c['label_chain'],c['author_chain'],c['protein'],c['observed_residues'],'gaps',[g['length'] for g in c['internal_gaps']],'ligands',c['ligand_codes_by_author_chain'])
