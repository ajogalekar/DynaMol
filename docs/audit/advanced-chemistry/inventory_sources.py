"""Preserve official inputs and inventory deposited advanced-chemistry evidence.

This performs no preparation, coordinate editing, oxidation-state inference, or
coordination assignment from proximity. Raw struct_conn rows remain available.
"""
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor
from collections import defaultdict
import datetime, hashlib, json, sys, urllib.request
import gemmi

OUT = Path(__file__).resolve().parent
METALS = {'LI','NA','K','RB','CS','MG','CA','SR','BA','AL','SC','TI','V','CR','MN','FE','CO','NI','CU','ZN','CD','HG','Y','ZR','MO','AG','AU','PT','PB'}
STANDARD = set('ALA ARG ASN ASP CYS GLN GLU GLY HIS ILE LEU LYS MET PHE PRO SER THR TRP TYR VAL'.split())

def clean(value):
    return None if value in (None, '.', '?') else gemmi.cif.as_string(value)

def rows(block, category):
    table = block.find_mmcif_category('_' + category + '.')
    names = [tag.split('.', 1)[1] for tag in table.tags]
    return [{key: clean(value) for key, value in zip(names, row)} for row in table]

def downloaded(path, url):
    path.parent.mkdir(parents=True, exist_ok=True)
    if not path.exists():
        with urllib.request.urlopen(url, timeout=60) as response:
            data = response.read(40_000_001)
        if len(data) > 40_000_000:
            raise ValueError('Input exceeded explicit 40 MB source limit')
        path.write_bytes(data)
    return {'path': str(path.relative_to(OUT)), 'source_url': url,
            'sha256': hashlib.sha256(path.read_bytes()).hexdigest(),
            'bytes': path.stat().st_size}

def endpoint(row, n):
    pre = f'ptnr{n}_'
    return {key: row.get(pre + field) for key, field in (
        ('label_chain', 'label_asym_id'), ('label_sequence_id', 'label_seq_id'),
        ('component', 'label_comp_id'), ('atom', 'label_atom_id'),
        ('author_chain', 'auth_asym_id'), ('author_residue_id', 'auth_seq_id'),
        ('symmetry', 'symmetry'))} | {
            'insertion_code': row.get(f'pdbx_ptnr{n}_PDB_ins_code'),
            'altloc': row.get(f'pdbx_ptnr{n}_label_alt_id')}

def inventory(code):
    result = {'id': code, 'recorded_at_utc': datetime.datetime.now(datetime.timezone.utc).isoformat()}
    try:
        path = OUT / 'inputs' / f'{code}.cif'
        result['input'] = downloaded(path, f'https://files.rcsb.org/download/{code}.cif')
        block = gemmi.cif.read_file(str(path)).sole_block()
        result['title'] = clean(block.find_value('_struct.title'))
        result['organisms'] = sorted({r[k] for cat,k in (
            ('entity_src_gen','pdbx_gene_src_scientific_name'),
            ('entity_src_nat','pdbx_organism_scientific')) for r in rows(block,cat) if r.get(k)})
        result['resolution_angstrom'] = clean(block.find_value('_refine.ls_d_res_high'))
        result['primary_citations'] = rows(block, 'citation')
        result['assemblies'] = rows(block, 'pdbx_struct_assembly')
        result['assembly_generators'] = rows(block, 'pdbx_struct_assembly_gen')
        result['assembly_operators'] = rows(block, 'pdbx_struct_oper_list')
        result['entities'] = rows(block, 'entity')
        result['components'] = rows(block, 'chem_comp')
        result['modified_residues'] = rows(block, 'pdbx_struct_mod_residue')
        result['raw_connections'] = rows(block, 'struct_conn')
        result['connections'] = [{'id':r['id'], 'type':r['conn_type_id'],
            'partner1':endpoint(r,1), 'partner2':endpoint(r,2),
            'deposited_distance_angstrom':r.get('pdbx_dist_value'),
            'deposited_order':r.get('pdbx_value_order'),
            'leaving_atom_flag':r.get('pdbx_leaving_atom_flag')}
            for r in result['raw_connections']]
        sites = [r for r in rows(block,'atom_site') if r.get('pdbx_PDB_model_num','1') == '1']
        observed, residues, metals = defaultdict(set), {}, []
        for r in sites:
            if r.get('label_seq_id'):
                observed[r['label_asym_id']].add(int(r['label_seq_id']))
            key = (r['label_asym_id'],r.get('label_seq_id'),r['auth_asym_id'],r['auth_seq_id'],r.get('pdbx_PDB_ins_code'),r['label_comp_id'])
            entry = residues.setdefault(key, {'label_chain':key[0], 'label_sequence_id':key[1], 'author_chain':key[2], 'author_residue_id':key[3], 'insertion_code':key[4], 'component':key[5], 'observed_heavy_atom_names':set()})
            if r['type_symbol'] not in ('H','D'):
                entry['observed_heavy_atom_names'].add(r['label_atom_id'])
            if r['type_symbol'].upper() in METALS:
                metals.append({'element':r['type_symbol'], 'component':r['label_comp_id'],
                    'label_chain':r['label_asym_id'],'author_chain':r['auth_asym_id'],
                    'author_residue_id':r['auth_seq_id'],'insertion_code':r.get('pdbx_PDB_ins_code'),
                    'atom':r['label_atom_id'],'altloc':r.get('label_alt_id'),
                    'occupancy':r.get('occupancy'),'source_formal_charge':r.get('pdbx_formal_charge'),
                    'cartesian_angstrom':[float(r[k]) for k in ('Cartn_x','Cartn_y','Cartn_z')]})
        result['metal_atoms'] = metals
        result['nonstandard_residues'] = [dict(r, observed_heavy_atom_names=sorted(r['observed_heavy_atom_names'])) for r in residues.values() if r['component'] not in STANDARD and r['component'] != 'HOH']
        result['water_residue_count'] = sum(r['component']=='HOH' for r in residues.values())
        result['atom_count_model1'] = len(sites)
        scheme = defaultdict(list)
        for r in rows(block,'pdbx_poly_seq_scheme'):
            scheme[r['asym_id']].append(r)
        result['chains'] = []
        for label, records in scheme.items():
            seen = observed[label]
            gaps, current = [], []
            if not seen:
                continue
            lo, hi = min(seen), max(seen)
            for r in sorted(records,key=lambda x:int(x['seq_id'])):
                n = int(r['seq_id'])
                if lo < n < hi and n not in seen:
                    current.append(r)
                elif current:
                    gaps.append(current); current=[]
            if current:gaps.append(current)
            result['chains'].append({'label_chain':label,'author_chain':records[0]['pdb_strand_id'],
                'entity_id':records[0]['entity_id'], 'sequence_length':len(records),'observed_residues':len(seen),
                'internal_gaps':[{'length':len(g),'residues':[r['mon_id'] for r in g],
                    'label_sequence_ids':[r['seq_id'] for r in g]} for g in gaps]})
        result['status'] = 'inventoried'
    except Exception as exc:
        result.update(status='error',error=f'{type(exc).__name__}: {exc}')
    target = OUT/'inventories'/f'{code}.json'; target.parent.mkdir(exist_ok=True)
    target.write_text(json.dumps(result,indent=2)+'\n')
    return result

def ccd(code):
    path=OUT/'inputs'/'ccd'/f'{code}.cif'
    source=downloaded(path,f'https://files.rcsb.org/ligands/download/{code}.cif')
    block=gemmi.cif.read_file(str(path)).sole_block()
    result={'id':code,'input':source,'chem_comp':rows(block,'chem_comp'),
        'atoms':rows(block,'chem_comp_atom'),'bonds':rows(block,'chem_comp_bond'),
        'descriptors':rows(block,'pdbx_chem_comp_descriptor')}
    target=OUT/'inventories'/'ccd'/f'{code}.json';target.parent.mkdir(parents=True,exist_ok=True)
    target.write_text(json.dumps(result,indent=2)+'\n')
    return result

if __name__=='__main__':
    with ThreadPoolExecutor(max_workers=6) as pool:
        results=list(pool.map(ccd if sys.argv[1]=='--ccd' else inventory,sys.argv[2:] if sys.argv[1]=='--ccd' else sys.argv[1:]))
    for r in results:
        print(r['id'],r.get('organisms'),r.get('resolution_angstrom'),r.get('title'),r.get('error',''))
        if 'chains' in r:
            print(' chains:', [(x['label_chain'],x['author_chain'],x['observed_residues'],[g['length'] for g in x['internal_gaps']]) for x in r['chains']])
            print(' chemistry:',sorted({x['component'] for x in r['nonstandard_residues']}))
            print(' covalent:',[(x['partner1'],x['partner2']) for x in r['connections'] if x['type']=='covale'])
            print(' metals:',len(r['metal_atoms']),'deposited metal links:',sum(x['type']=='metalc' for x in r['connections']))
