"""Read-only c6 and GDP source crosswalk; never assigns force-field parameters."""
import collections
import fnmatch
import hashlib
import itertools
import json
from pathlib import Path
import shlex

import parmed

HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[4]
AMBER=Path('/Users/ashujo/Library/Application Support/DynaMol/Runtimes/ambertools-b96b42253d7d46c5')
RESEARCH=Path('/Users/ashujo/.cache/dynamol-research/covalent-v1-focused')
CACHE=Path('/Users/ashujo/.cache/dynamol-research/complete-complex-track/c6-gdp-transfer-v1')
PRIOR=Path('/Users/ashujo/.cache/dynamol-research/complete-complex-track/panteva-applicability-v1')


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def off_atoms(path,residue):
    result=[]
    active=False
    for number,line in enumerate(path.read_text().splitlines(),1):
        if line.startswith('!'):
            active=line.startswith(f'!entry.{residue}.unit.atoms table')
            continue
        if active and line.strip():
            fields=shlex.split(line)
            result.append({'name':fields[0],'type':fields[1],'charge_e':float(fields[-1]),'atomic_number':int(fields[-2]),'source_line':number})
    if not result:
        raise ValueError(f'No exact OFF atom table for {residue}')
    return result


def polarizabilities(path):
    result={}
    for line_number,line in enumerate(path.read_text().splitlines(),1):
        fields=line.split(maxsplit=2)
        if not fields:
            continue
        if fields[0] in result:
            raise ValueError('Ambiguous duplicate polarizability key')
        result[fields[0]]={'alpha_angstrom3':float(fields[1]),'attribution':fields[2] if len(fields)>2 else '', 'line':line_number}
    return result


def mass_table(path):
    result={}
    for number,line in enumerate(path.read_text().splitlines()[1:],2):
        if not line.strip():
            break
        fields=line.split(maxsplit=3)
        result[fields[0]]={'mass':float(fields[1]),'mass_table_third_column':float(fields[2]),'definition':fields[3],'source_line':number}
    return result


def smallest_heavy_cycle(atom):
    choices=[]
    neighbors=[a for a in atom.bond_partners if a.atomic_number>1]
    for left,right in itertools.combinations(neighbors,2):
        queue=collections.deque([(left,[left.idx])])
        visited={atom.idx,left.idx}
        while queue:
            item,path=queue.popleft()
            if item is right:
                choices.append([atom.idx]+path)
                break
            for nxt in item.bond_partners:
                if nxt.atomic_number>1 and nxt.idx not in visited:
                    visited.add(nxt.idx)
                    queue.append((nxt,path+[nxt.idx]))
    return min(choices,key=len) if choices else None


def main():
    paths={
        'native_topology':RESEARCH/'6oim-complex-v3/solvated.prmtop',
        'mapping':RESEARCH/'6oim-complex-v3/source-mapping.json',
        'assembly':RESEARCH/'6oim-complex-v3/result.json',
        'gdp_scoping':RESEARCH/'gdp-scoped-v6/result.json',
        'gdp_source_prep':RESEARCH/'gdp-source/GDP.prep',
        'prior_inventory':HERE.parent/'panteva-applicability.json',
        'installed_pol':AMBER/'dat/leap/parm/lj_1264_pol.dat',
        'installed_gaff2':AMBER/'dat/leap/parm/gaff2.dat',
        'installed_nucleic':AMBER/'dat/leap/lib/nucleic12.lib',
        'installed_ff14sb_loader':AMBER/'dat/leap/cmd/oldff/leaprc.ff14SB',
        'installed_add1264':AMBER/'lib/python3.14/site-packages/parmed/tools/add1264.py',
        'installed_license':AMBER/'LICENSE',
        'paper':PRIOR/'paper.pdf', 'paper_text':PRIOR/'paper.txt',
        'SI':PRIOR/'supporting-information.pdf', 'SI_text':PRIOR/'supporting-information.txt',
        'retrieval':CACHE/'retrieval.json',
        'official_pol':CACHE/'lj_1264_pol.dat',
        'official_gaff2':CACHE/'gaff2.dat',
        'official_nucleic':CACHE/'nucleic12.lib',
        'official_ff14sb_loader':CACHE/'leaprc.ff14SB',
        'official_tutorial':CACHE/'amber-tutorial.html',
        'official_license':CACHE/'LICENSE',
        'metal_panel':ROOT/'docs/audit/advanced-chemistry/metal-panel.json',
        'covalent_panel':ROOT/'docs/audit/advanced-chemistry/covalent-panel.json',
        'implementation':Path(__file__)
    }
    hashes={name:{'path':str(path),'sha256':sha(path)} for name,path in paths.items()}
    scoped=json.loads(paths['gdp_scoping'].read_text())
    inverse={v:k for k,v in scoped['type_renames'].items()}
    prior=json.loads(paths['prior_inventory'].read_text())
    topology=parmed.load_file(str(paths['native_topology']))
    mapping={row['native_index']:row['source'] for row in json.loads(paths['mapping'].read_text())}
    installed_pol=polarizabilities(paths['installed_pol'])
    official_pol=polarizabilities(paths['official_pol'])
    installed_mass=mass_table(paths['installed_gaff2'])
    official_mass=mass_table(paths['official_gaff2'])
    c6=[]
    for atom in topology.atoms:
        if atom.type!='c6':
            continue
        cycle=smallest_heavy_cycle(atom)
        c6.append({'native_index':atom.idx,'name':atom.name,'charge_e':atom.charge,
                   'source_identity':mapping.get(atom.idx),'atomic_number':atom.atomic_number,
                   'native_nonbonded_type_index':int(topology.parm_data['ATOM_TYPE_INDEX'][atom.idx]),
                   'bond_partners':[{'index':other.idx,'name':other.name,'type':other.type,'atomic_number':other.atomic_number} for other in atom.bond_partners],
                   'smallest_heavy_cycle_indices':cycle,'smallest_heavy_cycle_size':len(cycle) if cycle else None,
                   'read_only_connectivity_scope':'Native bond graph only; no bond orders inferred or new valence assignment.',
                   'documented_polarizability_mapping_available':False,'assigned_alpha':None})
    gdp_residues=[r for r in topology.residues if r.name=='GDP']
    if len(gdp_residues)!=1:
        raise ValueError('Expected exactly one retained GDP')
    gdp=gdp_residues[0]
    native_gdp={a.name:a for a in gdp.atoms}
    gn={a['name']:a for a in off_atoms(paths['official_nucleic'],'GN')}
    # Both legacy hydrogen names are checked against their unique heavy parent.
    mapping_names={'H80':'H8','H1N':'H1'}
    hydrogen_parents={'H80':'C8','H1N':'N1'}
    base_names=['N9','C8','H80','N7','C5','C6','O6','N1','H1N','C2','N2','H21','H22','N3','C4']
    base=[]
    for name in base_names:
        atom=native_gdp[name]
        if name in hydrogen_parents:
            heavy=[p.name for p in atom.bond_partners if p.atomic_number>1]
            if heavy!=[hydrogen_parents[name]]:
                raise ValueError('Legacy GDP hydrogen identity mismatch')
        reference=gn[mapping_names.get(name,name)]
        original=inverse[atom.type]
        base.append({'native_index':atom.idx,'GDP_atom':name,'source_identity':mapping.get(atom.idx),
                     'GDP_scoped_type':atom.type,'GDP_original_type':original,'GDP_charge_e':atom.charge,
                     'GN_atom':reference['name'],'GN_type':reference['type'],'GN_charge_e':reference['charge_e'],
                     'GN_source_line':reference['source_line'],'charge_difference_e':atom.charge-reference['charge_e'],
                     'types_identical_after_recorded_scoping':original==reference['type'],
                     'GDP_default_alpha':installed_pol[original]['alpha_angstrom3'],
                     'GN_default_alpha':official_pol[reference['type']]['alpha_angstrom3']})
    phosphates=[]
    for name in ['O1A','O2A','O1B','O2B','O3B']:
        atom=native_gdp[name]
        original=inverse[atom.type]
        phosphates.append({'native_index':atom.idx,'name':name,'source_identity':mapping.get(atom.idx),
                            'scoped_type':atom.type,'original_type':original,'charge_e':atom.charge,
                            'bond_partners':[p.name for p in atom.bond_partners],
                            'documented_default_alpha':installed_pol[original]['alpha_angstrom3'],
                            'tutorial_OP_glob_matches_atom_name':fnmatch.fnmatchcase(name,'OP*'),
                            'DMP_fitted_correction_transfer_validated':False,'fitted_correction_assigned':False})
    assembly=json.loads(paths['assembly'].read_text())
    report={'schema_version':1,'stage':'read_only_c6_gdp_source_transfer_assessment',
            'case':prior['case'],'complete_complex_atom_counts':prior['complete_complex_atom_counts'],
            'source_hashes':hashes,
            'official_upstream':json.loads(paths['retrieval'].read_text()),
            'official_parameter_data_versions':{'installed_gaff2':paths['installed_gaff2'].read_text().splitlines()[0],
                                                'current_official_gaff2':paths['official_gaff2'].read_text().splitlines()[0],
                                                'polarizability_file_unchanged':sha(paths['installed_pol'])==sha(paths['official_pol']),
                                                'nucleic_template_file_unchanged':sha(paths['installed_nucleic'])==sha(paths['official_nucleic'])},
            'c6':{'atom_count':len(c6),'atoms':c6,'direct_installed_lookup':'c6' in installed_pol,
                  'direct_current_official_lookup':'c6' in official_pol,
                  'installed_mass_definition':installed_mass['c6'],'official_mass_definition':official_mass['c6'],
                  'source_code_behavior':'Exact-key lookup raises LJ12_6_4Error for any missing shared-class atom type; no c6 alias, hybridization, element, or equal-LJ fallback.',
                  'source_code_lines':[173,190],
                  'carbon_crosschecks':{name:{'GAFF2':installed_mass[name],'C4_source_alpha':installed_pol.get(name)} for name in ['c3','c5','c6','cx','cy']},
                  'status':'explicit_source_mapping_not_found',
                  'interpretation':'c6 is documented sp3 six-membered-ring carbon, but choosing another type\'s alpha would be a new crosswalk decision. GAFF2 mass-table numbers and C4-table alpha values differ; no replacement inferred.',
                  'parameters_assigned':False},
            'GDP':{'atom_count':len(gdp.atoms),'charge_e':sum(a.charge for a in gdp.atoms),
                   'retained_model':scoped['source_model'],'base_reference':'Amber ff14SB loader selects nucleic12.lib; GN is neutral guanosine. Current file matches installed source. Exact paper input topology was not recovered.',
                   'base_comparison':base,'all_15_base_charges_match_within_5e_9_e':all(abs(r['charge_difference_e'])<5e-9 for r in base),
                   'base_type_mismatches':[r for r in base if not r['types_identical_after_recorded_scoping']],
                   'purine_N7_transfer':'Local NB type/charge and surrounding base-charge agreement support a candidate transfer assessment, not validation of fitted binding thermodynamics for GDP(-3) or the complete protein site.',
                   'phosphate_atoms':phosphates,
                   'DMP_transfer':'Unestablished: DMP phosphodiester fit uses Dupradeau parameters; actual GDP uses Meagher diphosphate(-3) and O2/O3 site charges. Exact DMP atom/charge artifact not recovered; a donor-charge difference is not guessed.',
                   'documented_default_rule_coverage':True,'GDP_specific_fitted_validation_found':False,
                   'tutorial_selection_scope':'Example B is 1D23 DNA with TIP4P-Ew. G*/DG* N7 naming pattern would match GDP N7, but OP* does not match any of the five actual GDP nonbridging oxygen names. Pattern matches do not prove chemical transfer.'},
            'retrieval_limits':[{'url':'https://pubs.acs.org/doi/pdf/10.1021/ja00179a044','result':'Publisher full PDF unavailable; do not claim direct inspection of Miller full text.'},
                                {'url':'https://www.rsc.org/suppdata/cp/c0/c0cp00111b/c0cp00111b.pdf','result':'Web tool HTTP429; exact DMP dataset/charges not recovered.'},
                                {'url':'https://pmc.ncbi.nlm.nih.gov/articles/PMC2918240/','result':'Browser challenge; not bypassed.'},
                                {'path':str(CACHE/'meagher-2003.pdf'),'result':'Downloaded response is HTML, not PDF; excluded as a full-text scientific source. Publisher metadata and indexed author-repository abstract only.'}],
            'licenses':{'installed_Amber_dat_leap':'Public-domain grant already verified in installed LICENSE; no article/SI/code license inferred from it.',
                        'new_AmberClassic_root_license':'Root LICENSE discusses GPL msander code, not an independent new parameter-data clearance. Authoritative previous installed-data provenance retained.',
                        'new_papers_or_SI_copied_to_project':False},
            'claim_boundaries':['No full c6 mapping resolution achieved by name or LJ similarity.',
                                'Defined GDP O2/O3/NB defaults remain distinct from Panteva fitted site corrections.',
                                'Guanosine base-charge agreement is positive partial evidence, not an exact whole-GDP model match.',
                                'TIP4P-Ew conversion and exact native C4/GROMACS representation gaps remain unchanged.'],
            'next_bounded_step':'Obtain an explicitly sourced c6-to-polarizability crosswalk and the exact DMP/guanosine fit topology/charge artifacts; compare against retained GDP before any fitted override or complete-system C4 static test.',
            'parameters_assigned_or_changed':False,'topology_changed':False,'new_md_or_qm':False,'app_ready':False,'physical_model_validated':False}
    current={name:{'path':str(path),'sha256':sha(path)} for name,path in paths.items()}
    if current!=hashes:
        raise ValueError('Scientific source input changed during read-only audit')
    report['all_scientific_inputs_unchanged']=True
    checkpoint=ROOT/'docs/audit/advanced-chemistry/LIVE-JOBS.json'
    report['checkpoint_read']={'path':str(checkpoint),'sha256_at_report_time':sha(checkpoint),'shared_checkpoint_may_be_updated_by_root':True}
    (HERE/'evidence.json').write_text(json.dumps(report,indent=2,allow_nan=False)+'\n')
    print(json.dumps({'c6_count':len(c6),'c6_mapping_found':False,'base_charges_match':report['GDP']['all_15_base_charges_match_within_5e_9_e'],
                      'base_type_mismatch_atoms':[r['GDP_atom'] for r in report['GDP']['base_type_mismatches']],
                      'sources_unchanged':True,'output':str(HERE/'evidence.json')},indent=2))


if __name__=='__main__':
    main()
