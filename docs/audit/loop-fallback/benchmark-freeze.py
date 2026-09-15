"""Freeze source/chain-selected real missing-loop cases before native repair runs."""
from pathlib import Path
import os,sys,tempfile,json,hashlib,time,datetime
OUT=Path(__file__).resolve().parent;ROOT=OUT.parents[2]
DATA=Path(tempfile.mkdtemp(prefix='dynamol-loop-benchmark-inventory-'))
os.environ['DYNAMOL_DATA_DIR']=str(DATA);os.environ['DYNAMOL_CPU_THREADS']='2'
sys.path.insert(0,str(ROOT))
from backend import sources,storage,monomers,preparation
PICKS=[('8K5R','A'),('1UA2','A'),('1HCK','A'),('4HJO','A'),('1M17','A'),('2HYY','C'),('3HEG','A'),('6A93','B'),('2CG9','X'),('2BEL','C'),('1FPU','B'),('4HHY','B'),('2REN','A'),('1RNE','A'),('2ITY','A')]
SCREEN=json.loads((OUT/'benchmark-screening-inventory.json').read_text())
INV={r['id']:r for r in SCREEN['candidates']}
FROZEN={'frozen_at':datetime.datetime.now(datetime.timezone.utc).isoformat(),'scope':'15 distinct experimental PDB entries with actual sequence-supported internal gaps; native preparation only. Selection considers protein/ligand categories and gap lengths before any native preparation outcomes. This convenience panel is not a random or universal success-rate estimate.','original_source_files_unchanged':True,'source_screening_inventory':'benchmark-screening-inventory.json','source_screening_sha256':hashlib.sha256((OUT/'benchmark-screening-inventory.json').read_bytes()).hexdigest(),'inventory_data_dir':str(DATA),'cases':[]}
FROZEN['excluded_candidates']=[{'id':r['id'],'reason':'Not selected from pre-preparation inventory: no internal gap in relevant chain, redundant chemistry, or panel-size limit. All actual screened source/gap inventories retained.'} for r in SCREEN['candidates'] if r['id'] not in {x[0] for x in PICKS}]
for accession,author in PICKS:
 r=INV[accession];chain=next(c for c in r['chains'] if c['author_chain']==author)
 expectation=('expected_unsupported_size' if accession=='2ITY' else 'expected_unsupported_covalent_glycan' if accession in {'2REN','1RNE'} else 'within_loop_size_budget')
 row={'id':accession,'pdb_id':accession,'input_path':str(Path(r['input_path']).relative_to(ROOT)),'input_sha256':r['sha256'],'source_url':r['source_url'],'protein_identity':chain['protein'],'organisms':r['organisms'],'author_chain':author,'label_chain':chain['label_chain'],'chain_index':None,'keep_associated_molecules':True,'source_gap_inventory':chain['internal_gaps'],'source_ligand_codes_by_author_chain':chain['ligand_codes_by_author_chain'],'scope_expectation':expectation,'settings':{'engine':'openmm','ph':7,'seed':42,'build_missing_residues':True,'add_missing_atoms':True,'optimize_sidechains':True,'remove_waters':True,'remove_heterogens':False,'ligand_actions':{}},'selection_note':'One observed protein chain with associated molecules retained by existing application selection; this does not assert a biological monomer.'}
 FROZEN['cases'].append(row)
(OUT/'benchmark-case-selection.json').write_text(json.dumps(FROZEN,indent=2)+'\n')
for row in FROZEN['cases']:
 begin=time.monotonic()
 try:
  meta=sources.import_structure(ROOT/row['input_path'],name='Loop benchmark '+row['id'],provenance={'source_url':row['source_url'],'source_sha256':row['input_sha256'],'audit':'loop-fallback inventory'})
  row['source_dataset_id']=meta['id']
  chain_options=monomers.list_monomers(meta['id']);row['monomer_choices']=chain_options
  selected=next((c for c in chain_options['chains'] if c['chain_id']==row['author_chain']),None)
  if selected is None:
   selected=next(c for c in chain_options['chains'] if c['chain_id']==row['label_chain'])
   row['chain_identity_note']='Source author chain differs from imported label-asym chain; both identifiers are retained. This is an observed import behavior, not a sequence remapping workaround.'
  row['chain_index']=selected['index'];row['imported_chain_id']=selected['chain_id']
  extracted=monomers.create_monomer(meta['id'],monomers.MonomerRequest(chain_index=row['chain_index'],keep_associated_molecules=True))
  row['inventory_dataset_id']=extracted['id'];row['monomer_selection']=extracted.get('monomer_selection');row['inspection']=preparation.inspect_preparation(extracted['id'])
  try:
   preparation.validate_preparation({**row['settings'],'dataset_id':extracted['id']});row['validation_ready']=True
  except Exception as e:row['validation_ready']=False;row['validation_error']=str(e)
  row['inventory_status']='inspected';row['source_unchanged']=hashlib.sha256((ROOT/row['input_path']).read_bytes()).hexdigest()==row['input_sha256']
 except Exception as e:row['inventory_status']='failed';row['inventory_error']=repr(e)
 row['source_unchanged']=hashlib.sha256((ROOT/row['input_path']).read_bytes()).hexdigest()==row['input_sha256']
 row['inventory_seconds']=round(time.monotonic()-begin,2)
 (OUT/'benchmark-manifest.json').write_text(json.dumps(FROZEN,indent=2,allow_nan=False)+'\n')
 print(row['id'],'chain_index',row['chain_index'],'inspection',row['inventory_status'],'ready',row.get('validation_ready'),'error',row.get('validation_error',row.get('inventory_error','')),flush=True)
FROZEN['inventory_complete_at']=datetime.datetime.now(datetime.timezone.utc).isoformat()
(OUT/'benchmark-manifest.json').write_text(json.dumps(FROZEN,indent=2,allow_nan=False)+'\n')
print('MANIFEST',OUT/'benchmark-manifest.json',flush=True)
