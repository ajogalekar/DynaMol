"""Inventory one published Mg candidate without assigning or changing parameters."""
from __future__ import annotations

import argparse
from collections import Counter, defaultdict
from decimal import Decimal
import hashlib
import json
from pathlib import Path
import re
import shutil
import xml.etree.ElementTree as ET

import parmed

from audit_mg_coordination import digest


def sections(path):
    result = defaultdict(list)
    section = None
    for number, line in enumerate(path.read_text().splitlines(), 1):
        content = line.split(';', 1)[0].strip()
        if not content or content.startswith(('#', '*')):
            continue
        match = re.fullmatch(r'\[\s*(\S+)\s*\]', content)
        if match:
            section = match[1]
        elif section:
            result[section].append((number, content.split()))
    return result


def atomtypes(path):
    result = {}
    for number, row in sections(path)['atomtypes']:
        if len(row) != 7 or row[4] != 'A' or row[0] in result:
            raise ValueError('Unsupported/duplicate atomtype row')
        result[row[0]] = {'atomic_number': int(row[1]), 'mass_da': float(row[2]),
                         'type_charge_field_e': float(row[3]), 'sigma_nm': float(row[5]),
                         'epsilon_kj_mol': float(row[6]), 'sigma_text': row[5],
                         'epsilon_text': row[6], 'line': number}
    return result


def overrides(path):
    result = {}
    for number, row in sections(path)['nonbond_params']:
        if len(row) != 5 or float(row[2]) != 1:
            raise ValueError('Unsupported nonbonded pair row')
        if 'mMg' not in row[:2]:
            continue
        partner = row[1] if row[0] == 'mMg' else row[0]
        if partner in result:
            raise ValueError('Duplicate microMg pair')
        result[partner] = {'sigma_nm': float(row[3]), 'epsilon_kj_mol': float(row[4]),
                           'sigma_text': row[3], 'epsilon_text': row[4], 'line': number}
    return result


def rtp_atoms(path, residue):
    section = None
    selected = False
    rows = {}
    for number, line in enumerate(path.read_text().splitlines(), 1):
        content = line.split(';', 1)[0].strip()
        match = re.fullmatch(r'\[\s*(\S+)\s*\]', content)
        if match:
            name = match[1]
            if name not in ('atoms', 'bonds', 'impropers', 'dihedrals', 'exclusions', 'cmap'):
                selected = name == residue
            section = name
        elif selected and section == 'atoms' and content:
            row = content.split()
            rows[row[0]] = {'type': row[1], 'charge_e': float(row[2]), 'line': number}
    if not rows:
        raise ValueError('Residue source not found: ' + residue)
    return rows


def within_printed_precision(value, text):
    half = .5 * float(Decimal(10) ** Decimal(text).as_tuple().exponent)
    return abs(value - float(text)) <= half + 2e-9


def native_parameters(atom):
    return {'sigma_nm': atom.sigma / 10, 'epsilon_kj_mol': atom.epsilon * 4.184,
            'rmin_half_angstrom': atom.rmin, 'epsilon_kcal_mol': atom.epsilon}


def run(root, evidence):
    destination = evidence / 'compatibility-inventory.json'
    if destination.exists():
        raise ValueError('Existing inventory must be preserved')
    pinned = root / 'mg-gdp-model-assessment-v1'
    retrieval = json.loads((pinned / 'retrieval.json').read_text())
    source_manifest = json.loads((evidence / 'source-files.json').read_text())
    tree = json.loads((evidence / 'pinned-tree.json').read_text())
    if tree['truncated'] or tree['sha'] != retrieval['commit']:
        raise ValueError('Pinned recursive tree identity failed')
    blobs = {r['path']: r for r in tree['tree'] if r['type'] == 'blob'}
    used = [pinned / 'retrieval.json', evidence / 'retrieval.json', evidence / 'source-files.json',
            evidence / 'pinned-tree.json', evidence / 'published-micromg-paper.xml',
            evidence / 'paper-retrieval.json']
    for name, item in retrieval['files'].items():
        path = pinned / name
        if digest(path) != item['sha256']:
            raise ValueError('Original pinned source changed')
        raw = path.read_bytes()
        if hashlib.sha1(b'blob ' + str(len(raw)).encode() + b'\0' + raw).hexdigest() != blobs[name]['sha']:
            raise ValueError('Original source differs from pinned tree')
        used.append(path)
    for name, item in source_manifest.items():
        path = evidence / 'source' / name
        if digest(path) != item['sha256'] or item['git_blob_sha1'] != blobs[name]['sha']:
            raise ValueError('Additional pinned source changed')
        used.append(path)
    inputs = root / '6oim-complex-v3'
    used += [inputs / 'solvated.prmtop', inputs / 'result.json', inputs / 'source-mapping.json',
             root / 'gdp-scoped-v6/result.json']
    hashes = {str(p): digest(p) for p in used}
    assembly = json.loads((inputs / 'result.json').read_text())
    if digest(inputs / 'solvated.prmtop') != assembly['outputs_sha256']['solvated.prmtop']:
        raise ValueError('Corrected complex changed')
    native = parmed.load_file(str(inputs / 'solvated.prmtop'))
    mapping = {r['native_index']: r for r in json.loads((inputs / 'source-mapping.json').read_text())}
    scoped = json.loads((root / 'gdp-scoped-v6/result.json').read_text())
    reverse = {v: k for k, v in scoped['type_renames'].items()}
    published_atoms = atomtypes(pinned / 'ff_Mg.itp')
    published_pairs = overrides(pinned / 'ff_Mg.itp')
    source = evidence / 'source'
    base = source / 'example_MgCl2/amber14sb.ff'
    base_atoms = atomtypes(base / 'ffnonbonded.itp')
    for name in ('example_MgCl2/ions.ff/ff_Mg.itp', 'example_1y26withMg/ions.ff/ff_Mg.itp'):
        example = overrides(source / name)
        if {k: (v['sigma_nm'], v['epsilon_kj_mol']) for k, v in example.items()} != {
                k: (v['sigma_nm'], v['epsilon_kj_mol']) for k, v in published_pairs.items()}:
            raise ValueError('Pinned examples disagree with root microMg overrides')
    defaults = sections(base / 'forcefield.itp')['defaults']
    if defaults[0][1][:3] != ['1', '2', 'yes']:
        raise ValueError('Published Lorentz-Berthelot assumption changed')
    rna = rtp_atoms(source / 'example_1y26withMg/amber14sb.ff/rna.rtp', 'RG')
    ser = rtp_atoms(source / 'example_1y26withMg/amber14sb.ff/aminoacids.rtp', 'SER')
    paper = ET.fromstring((evidence / 'published-micromg-paper.xml').read_text())
    licenses = [ET.tostring(x, encoding='unicode') for x in paper.findall('.//license')]
    if not any('creativecommons.org/licenses/by/4.0' in x for x in licenses):
        raise ValueError('Article license evidence changed')
    table = next(x for x in paper.findall('.//table-wrap') if x.findtext('label') == 'Table 1')
    paper_values = []
    for row in table.findall('.//tbody/tr'):
        cells = [''.join(x.itertext()) for x in row.findall('td')]
        paper_values.append({'symbol': cells[0], 'unit': cells[1], 'microMg': cells[2]})
    groups = defaultdict(list)
    def category(atom):
        r = atom.residue.name
        if r == 'GDP': return 'GDP'
        if r == 'COV': return 'covalent_drug_residue'
        if r == 'WAT': return 'TIP3P_water'
        if atom.atomic_number == 12: return 'magnesium'
        if r == 'Na+': return 'sodium_counterion'
        return 'protein'
    for atom in native.atoms:
        groups[(category(atom), atom.type)].append(atom)
    records = []
    mg = published_atoms['mMg']
    for (component, kind), atoms in sorted(groups.items()):
        a = atoms[0]
        params = native_parameters(a)
        if any(native_parameters(x) != params for x in atoms):
            raise ValueError('Same native type has inconsistent Lennard-Jones parameters')
        original = reverse.get(kind, kind)
        direct = kind in published_pairs
        alias = component == 'GDP' and original in published_pairs
        partner = kind if direct else original if alias else None
        source_atom = published_atoms.get(partner or kind)
        signature = None if source_atom is None else {
            'sigma_matches_source_printed_precision': within_printed_precision(params['sigma_nm'], source_atom['sigma_text']),
            'epsilon_matches_source_printed_precision': within_printed_precision(params['epsilon_kj_mol'], source_atom['epsilon_text'])}
        if component == 'magnesium':
            applicability = 'candidate_replaces_current_mg_parameters_no_assignment_made'
        elif component == 'TIP3P_water':
            applicability = 'published_water_family_but_exact_oxygen_lj_variant_differs'
        elif component == 'GDP':
            applicability = 'not_validated_for_retained_diphosphate_model'
        elif direct:
            applicability = 'same_type_named_in_rna_overrides_application_to_this_component_unvalidated'
        else:
            applicability = 'no_component_specific_override_standard_mixing_would_still_define_a_pair'
        records.append({'component': component, 'native_type': kind, 'original_type_before_gdp_scoping': original,
                        'atom_count': len(atoms), 'atomic_numbers': sorted({x.atomic_number for x in atoms}),
                        'charge_range_e': [min(x.charge for x in atoms), max(x.charge for x in atoms)],
                        'native_lj': params, 'direct_name_override_exists': direct,
                        'override_exists_only_after_gdp_unscoping': alias,
                        'published_partner': partner,
                        'published_override': published_pairs.get(partner),
                        'source_lj_numeric_comparison': signature,
                        'nominal_unmodified_lb_pair_if_mg_alone_replaced': {
                            'sigma_nm': (mg['sigma_nm'] + params['sigma_nm']) / 2,
                            'epsilon_kj_mol': (mg['epsilon_kj_mol'] * params['epsilon_kj_mol']) ** .5},
                        'applicability': applicability,
                        'source_atom_examples': [{'native_index': x.idx, 'atom': x.name,
                                                 'source': mapping.get(x.idx, {}).get('source')} for x in atoms[:3]]})
    by_index = {a.idx: a for a in native.atoms}
    specific = {}
    for label, index in [('GDP_O2B', 2778), ('SER17_OG', 323), ('water_O', 2835), ('Mg', 2775)]:
        a = by_index[index]
        specific[label] = {'native_index': index, 'source': mapping[index]['source'], 'native_type': a.type,
                           'original_type': reverse.get(a.type, a.type), 'charge_e': a.charge,
                           **native_parameters(a)}
    oxygen = specific['water_O']
    published_water = base_atoms['OW']
    water_compare = {'current': oxygen, 'published_example_OW': published_water,
                     'sigma_difference_nm': oxygen['sigma_nm'] - published_water['sigma_nm'],
                     'epsilon_difference_kj_mol': oxygen['epsilon_kj_mol'] - published_water['epsilon_kj_mol'],
                     'published_example_hydrogen_and_oxygen_charges': [-.834, .417, .417],
                     'parameter_family': 'TIP3P', 'exact_lj_match': False,
                     'same_water_family_is_not_exact_reproduction': True}
    actual_gdp_types = [r for r in records if r['component'] == 'GDP']
    counts = Counter(category(a) for a in native.atoms)
    report = {'stage': 'pinned_micromg_complete_complex_applicability_inventory',
        'case': '6OIM', 'app_ready': False, 'physical_model_validated': False,
        'ready_for_direct_complete_complex_substitution': False,
        'parameters_assigned_or_changed': False, 'new_md_or_qm': False,
        'pinned_repository': {'url': 'https://github.com/bio-phys/Magnesium-FFs', 'commit': retrieval['commit'],
                              'source_blob_count': len(blobs), 'recursive_tree_complete': True,
                              'root_and_both_examples_microMg_overrides_equal': True},
        'paper': {'doi': '10.1021/acs.jctc.0c01281', 'full_text_source': 'https://www.ebi.ac.uk/europepmc/webservices/rest/PMC8047801/fullTextXML',
                  'table_1': paper_values, 'license_evidence': licenses,
                  'license_url': 'https://creativecommons.org/licenses/by/4.0/',
                  'calibration_scope': 'TIP3P hydration and exchange; Mamatkulov-Schwierz chloride; DMP with donor charge/type/angles changed to parmBSC0chiOL3 RNA, then RNA tests.',
                  'no_complete_KRAS_GDP_Ser17_covalent_validation_reported': True},
        'source_model': {'microMg_atomtype': mg, 'required_combination_rule': 'Lorentz-Berthelot with listed pair overrides',
                         'published_microMg_pairs': published_pairs, 'pair_count_including_self_and_chloride': len(published_pairs),
                         'nMg_excluded_from_this_inventory': True,
                         'published_type_charge_fields_are_not_residue_partial_charges': True},
        'complete_complex_atom_counts': dict(counts), 'native_unique_type_count': len({a.type for a in native.atoms}),
        'component_type_inventory': records, 'retained_site_atoms': specific,
        'GDP_applicability': {'actual_type_count': len(actual_gdp_types),
                             'types_with_no_override_even_after_unscoping': [r['original_type_before_gdp_scoping'] for r in actual_gdp_types if not r['override_exists_only_after_gdp_unscoping']],
                             'source_model': scoped['source_model'], 'RNA_RG_phosphate_oxygen': rna['O1P'],
                             'bound_GDP_oxygen_charge_difference_from_RNA_RG_e': specific['GDP_O2B']['charge_e'] - rna['O1P']['charge_e'],
                             'O3_to_O2_mapping_not_authorized_or_applied': True,
                             'diphosphate_is_not_the_published_phosphodiester_calibration': True},
        'SER17_applicability': {'source_SER_OG': ser['OG'], 'RNA_RG_ribose_OH': rna["O2'"],
                               'published_mMg_OH_override': published_pairs['OH'],
                               'shared_OH_type_does_not_establish_protein_site_validation': True},
        'water_compatibility': water_compare,
        'redistribution_evidence': {'license_named_files_in_pinned_tree': [name for name in blobs if re.search(r'(^|/)(licen[cs]e|copying|copyright|notice)([.\-_]|$)', name, re.I)],
                                    'repository_wide_explicit_grant_found': False,
                                    'scope': 'No license-named file in complete pinned tree and no explicit Mg code grant in retrieved model files/README. Third-party file headers were not all audited.',
                                    'article_CC_BY_4_0_verified_from_full_text': True,
                                    'possible_route': 'Independently implement equations/parameters from the licensed article with attribution and a record of changes; reconcile published precision and source differences before using it. Do not infer a license for copied repository files.'},
        'blocking_gaps_for_direct_substitution': [
            'Retained GDP O3 (including coordinating O2B) and CK have no published named pair override; scoped QA-QU names match no raw overrides.',
            'GDP charges and diphosphate chemistry differ from the RNA-matched DMP calibration; charge/type/angle equivalence is not established.',
            'Protein and GAFF covalent-adduct environments are outside the demonstrated calibration; sharing a type name does not supply that validation.',
            'Current TIP3P oxygen LJ values differ slightly from the pinned example.',
            'Native sodium counterions have no corresponding calibrated Mg-sodium pair in this source; published chloride correction does not apply to them.',
            'Repository redistribution grant is not established; licensed-paper reimplementation has a route but needs exact numerical provenance.'],
        'next_alternative': {'model': 'Panteva-Giambasu-York modified 12-6-4', 'doi': '10.1021/acs.jpcb.5b10423',
                             'primary_source': 'https://theory.rutgers.edu/resources/pdfs/Panteva_JChemPhysB_2015_v119_p15460.pdf',
                             'reason_to_inspect': 'Published site-specific phosphate/purine calibration addresses binding balance explicitly.',
                             'known_limits': 'TIP4P-Ew and defined nucleic-acid charge models; it is not a direct TIP3P/GDP/Ser replacement. Requires coherent water/cofactor model choice and full C4 pair transport.',
                             'accepted_or_implemented': False},
        'input_sha256': hashes}
    if sum(counts.values()) != 27477 or counts['GDP'] != 40:
        raise ValueError('Complete complex inventory changed')
    if any(digest(p) != sha for p, sha in hashes.items()):
        raise ValueError('Source changed during inventory')
    shutil.copy2(__file__, evidence / 'inventory-implementation.py')
    destination.write_text(json.dumps(report, indent=2) + '\n')
    print(json.dumps({'output': str(destination), 'native_types': report['native_unique_type_count'],
                      'components': dict(counts), 'published_pairs': len(published_pairs),
                      'GDP_missing_override_types': report['GDP_applicability']['types_with_no_override_even_after_unscoping'],
                      'water': water_compare}, indent=2))


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--root', type=Path, required=True)
    parser.add_argument('--evidence', type=Path, required=True)
    args = parser.parse_args()
    run(args.root, args.evidence)
