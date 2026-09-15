from pathlib import Path
import os,sys,tempfile,json,time,hashlib,datetime
OUT=Path(__file__).resolve().parent;ROOT=OUT.parents[2]
DATA=Path(tempfile.mkdtemp(prefix='dynamol-loop-identity-reinspection-'))
os.environ['DYNAMOL_DATA_DIR']=str(DATA);os.environ['DYNAMOL_CPU_THREADS']='2';sys.path.insert(0,str(ROOT))
from backend import sources,storage,monomers,preparation
m=json.loads((OUT/'benchmark-manifest.json').read_text());report={'data_root':str(DATA),'started_at':datetime.datetime.now(datetime.timezone.utc).isoformat(),'cases':[]}
for case in m['cases']:
 if case['id'] not in {'6A93','2CG9','2REN','1RNE','4HHY'}:continue
 begin=time.monotonic();r={'id':case['id'],'chain_index':case['chain_index']}
 try:
  source=sources.import_structure(ROOT/case['input_path'],name='Identity check '+case['id']);r['source_dataset_id']=source['id']
  d=monomers.create_monomer(source['id'],monomers.MonomerRequest(chain_index=case['chain_index'],keep_associated_molecules=True));r['selected_dataset_id']=d['id'];r['inspection']=preparation.inspect_preparation(d['id']);r['monomer_selection']=d['monomer_selection']
  try:preparation.validate_preparation({**case['settings'],'dataset_id':d['id']});r['validation_ready']=True
  except Exception as e:r['validation_ready']=False;r['validation_error']=str(e)
  r['status']='inspected'
 except Exception as e:r['status']='error';r['error']=repr(e)
 r['elapsed_seconds']=round(time.monotonic()-begin,2);report['cases'].append(r)
 (OUT/'benchmark-identity-reinspection.json').write_text(json.dumps(report,indent=2)+'\n')
 print(r['id'],r['status'],r.get('validation_ready'),r.get('validation_error',r.get('error')),[x.get('count',len(x.get('residues',[]))) for x in r.get('inspection',{}).get('missing_residues',[])],flush=True)
