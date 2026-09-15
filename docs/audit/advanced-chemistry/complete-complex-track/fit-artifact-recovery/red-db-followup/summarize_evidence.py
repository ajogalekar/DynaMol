"""Analyze cached public source artifacts; no parameter application or fitting."""
from datetime import datetime, timezone
from pathlib import Path
import hashlib
import html
import json
import math
import re

OUT = Path(__file__).resolve().parent
CACHE = Path('/Users/ashujo/.cache/dynamol-research/complete-complex-track/fit-artifact-recovery-red-db-v1')


def file_ref(path):
    path = Path(path)
    return {'path': str(path.resolve()), 'sha256': hashlib.sha256(path.read_bytes()).hexdigest()}


def clean(s):
    return ' '.join(html.unescape(re.sub('<[^>]*>', ' ', s)).split())


def project_fields(project):
    summary = (CACHE / project / 'summary.html').read_text()
    fields = {}
    for label in ('Firstname', 'Lastname', 'Institute', 'City', 'Country',
                  'Author(s)', 'Journal', 'Year', 'Volume', 'Page(s)'):
        match = re.search(r'<b>' + re.escape(label) + r'</b></span>(.*?)</p>', summary, re.S)
        fields[label] = clean(match.group(1)) if match else None
    text = clean(summary)
    for label in ('Upload', 'Update'):
        match = re.search(label + r' (\d{2}-\d{2}-\d{4} \d{2}:\d{2})', text)
        fields[label] = match.group(1) if match else None
    return fields


def main():
    inventory = json.loads((OUT / 'artifact-inventory.json').read_text())
    projects = []
    for p in inventory['projects']:
        project = p['project']
        mol = p['mol2']
        atoms = {a['id']: a for a in mol['atoms']}
        neighbors = {a: [] for a in atoms}
        for bond in mol['bonds']:
            neighbors[bond['a']].append(bond['b'])
            neighbors[bond['b']].append(bond['a'])
        nonbridging_oxygen = [a for a in atoms.values() if a['type'] == 'O'
                             and len(neighbors[a['id']]) == 1
                             and atoms[neighbors[a['id']][0]]['type'] == 'P']
        bridging_oxygen = [a for a in atoms.values() if a['type'] == 'O'
                          and sorted(atoms[i]['type'] for i in neighbors[a['id']]) == ['C', 'P']]
        inspected = CACHE / project / 'inspected-members'
        placeholder_files = [q.name for q in sorted(inspected.glob('script*.ff'))
                             if 'not yet provided by the author' in q.read_text()]
        resp = (inspected / 'input1.in').read_text()
        atom_lines = resp.split('Column not used by RESP', 1)[1].splitlines()[1:]
        constraints = []
        for line in atom_lines:
            f = line.split()
            if len(f) == 3:
                constraints.append({'atomic_number': int(f[0]), 'equivalence': int(f[1]), 'atom_id': int(f[2])})
        finite = all(math.isfinite(v) for a in atoms.values() for v in [*a['coordinates_A'], a['charge_e']])
        fields = project_fields(project)
        projects.append({
            'project': project, 'summary_url': p['summary_url'], 'metadata': fields,
            'metadata_time_zone_as_stated': 'Paris time (Day-Month-Year)',
            'submission_author_identity': {'first_name': fields['Firstname'], 'last_name': fields['Lastname'],
                                            'note': 'Blank live-page submission identity for W-13; do not infer the submitter from publication authors.' if project == 'W-13' else 'Named on live public project page.'},
            'publication_linkage': {'database_cites': 'Dupradeau et al., Phys. Chem. Chem. Phys. 2010, 12, 7821–7839',
                'doi_resolved_from_prior_primary_paper': '10.1039/C0CP00111B',
                'is_panteva_reference_51': True,
                'does_panteva_name_this_project': False,
                'does_project_name_panteva': False,
                'exact_panteva_model_identified': False},
            'source_method': {'molecule': 'dimethylphosphate, gauche/gauche conformer',
                'formal_total_charge_e': -1, 'method': 'HF/6-31G*; Connolly surface; two-stage RESP',
                'optimization_and_mep_program': 'GAMESS-US' if project in ('W-11','W-12') else 'Gaussian 1998',
                'molecular_orientation': 'program reorientation algorithm; one source conformer/orientation',
                'fit_strengths_source': {'stage_1_qwt': 0.0005, 'stage_2_qwt': 0.001},
                'fit_strengths_assigned': False},
            'integrity': {'atoms': mol['atom_count'], 'bonds': mol['bond_count'],
                'total_charge_e': mol['total_charge_e'], 'finite_coordinates_and_charges': finite,
                'atom_count_matches_resp': len(constraints) == mol['atom_count'] == 13,
                'net_charge_matches_resp': abs(mol['total_charge_e'] + 1) < 1e-7,
                'net_charge_residual_e': round(mol['total_charge_e'] + 1, 8),
                'charge_residual_interpretation': 'Four-decimal MOL2 charges sum to the reported value; no renormalization performed. The −0.0002 e residual for W-13/14 is compatible with output rounding, not proof of a different intended integer charge.',
                'mol2_types_are_element_labels': sorted(set(a['type'] for a in atoms.values())),
                'all_mol2_bond_orders_single': all(b['type'] == '1' for b in mol['bonds']),
                'complete_force_field_type_and_parameter_bundle': False},
            'nonbridging_oxygen_graph': [{'atom_id': a['id'], 'name': a['name'], 'charge_e': a['charge_e'],
                'neighbor': atoms[neighbors[a['id']][0]]['name']} for a in nonbridging_oxygen],
            'bridging_oxygen_graph': [{'atom_id': a['id'], 'name': a['name'], 'charge_e': a['charge_e'],
                'neighbors': [atoms[i]['name'] for i in neighbors[a['id']]]} for a in bridging_oxygen],
            'resp_stage_1_equivalence_records': [c for c in constraints if c['equivalence'] > 0],
            'placeholder_files': placeholder_files,
            'listed_converter_missing_from_archive': project == 'W-13',
            'redistribution_status': 'Public download verified; page says projects free but displays rights-reserved footer. No explicit redistribution license recovered from these archives.',
            'mol2_file': file_ref(inspected / 'tripos1.mol2'),
            'archive_file': file_ref(CACHE / project / (project + '.tar.bz2')),
        })
    # Comparison only, to frozen previously inspected GDP evidence; no model assignment.
    prior = OUT.parents[1] / 'c6-gdp-transfer' / 'evidence.json'
    assert prior.is_file(), prior
    prior_data = json.loads(prior.read_text())
    gdp_oxygen = {a['name']: a for a in prior_data['GDP']['phosphate_atoms']}
    assert all(round(gdp_oxygen[name]['charge_e'], 4) == -0.9474 for name in ('O1A','O2A'))
    assert all(round(gdp_oxygen[name]['charge_e'], 4) == -0.9552 for name in ('O1B','O2B','O3B'))
    result = {
        'schema': 'dynamol.reddb-dmp-recovery-assessment.v1',
        'created_utc': datetime.now(timezone.utc).isoformat(),
        'scope': 'Bounded public source recovery and read-only charge/graph comparison.',
        'lookup_budget': {'used': 2, 'maximum': 3, 'third_lookup_unneeded': True, 'guanosine_lookup_performed': False},
        'source_manifest': file_ref(OUT / 'artifact-inventory.json'),
        'inspection_script': file_ref(OUT / 'inspect_public_artifacts.py'),
        'prior_gdp_assessment': file_ref(prior),
        'diagnostic_corrections': [{'failure': 'Initial summary assertion searched floating-point JSON values as exact text substrings.',
            'evidence': 'Native Amber charge values contain normalization rounding, e.g. O1A −0.9473999989024436 e.',
            'correction': 'Use explicitly named phosphate atoms and numerical comparison rounded to the original four-decimal source precision.',
            'scientific_input_changed': False}],
        'projects': projects,
        'gdp_transfer_comparison': {'actual_gdp_charge_model': 'Meagher–Redman–Carlson GDP(−3), Amber99 based, preserved in corrected 6OIM',
            'native_nonbridging_oxygen_charge_groups_e': {'O1A/O2A': -0.9474, 'O1B/O2B/O3B': -0.9552},
            'actual_mg_donor': 'GDP O2B, original type O3, scoped type QR',
            'native_gdp_mg_donor_charge_e': -0.9552,
            'native_gdp_phosphate_atoms': [gdp_oxygen[name] for name in ('O1A','O2A','O1B','O2B','O3B')],
            'recovered_dmp_nonbridging_charge_options_e': [-0.7956, -0.7952],
            'gdp_mg_donor_minus_dmp_charge_options_e': [-0.1596, -0.1600],
            'interpretation': 'The recovered candidate DMP charges differ from native GDP. A charge difference alone neither proves failure nor validates transfer of the fitted correction. No exact Panteva source project or full force-field topology was identified.'},
        'remaining_gaps': ['Which exact DMP charge/orientation project, or other derivative, was used in Panteva 2015?',
            'What force-field atom types, bonded terms, and nonbonded terms completed the DMP simulation topology?',
            'Published or independent evidence for applying its fitted Mg correction to the preserved GDP(−3) charge model.',
            'Documented GAFF2 c6 polarizability crosswalk remains unrecovered by this DMP-only follow-up.',
            'Explicit redistribution terms for any recovered data intended to be bundled.'],
        'not_performed': ['parameter assignment', 'QM', 'MD', 'topology modification', 'application integration', 'authentication bypass', 'author contact'],
        'next_bounded_step': 'Source-only review of Dupradeau 2010 manuscript/SI for the W-11–W-14 orientation comparison and any exact DMP library identifiers; accept an exact Panteva identity only if a primary author artifact explicitly selects the charge/type bundle. Otherwise retain as candidates and require complete-system research validation before support.'
    }
    (OUT / 'evidence.json').write_text(json.dumps(result, indent=2) + '\n')
    print(json.dumps({'projects': [{'project':p['project'],'metadata':p['metadata'],'integrity':p['integrity']} for p in projects],
                      'evidence': file_ref(OUT / 'evidence.json')}, indent=2))


if __name__ == '__main__':
    main()
