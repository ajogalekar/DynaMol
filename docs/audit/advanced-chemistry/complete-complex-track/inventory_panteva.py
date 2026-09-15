"""Read-only applicability inventory for one published modified 12-6-4 Mg model.

No topology assignment, model fitting, engine context, MD or QM is performed.
The computed C4 values are source-rule arithmetic, not approved parameters.
"""
from __future__ import annotations

import argparse
import ast
from collections import Counter, defaultdict
import hashlib
import json
from pathlib import Path
import re
import shlex
import shutil

import parmed


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def nonbond_rows(path):
    active = False
    result = {}
    for line_no, line in enumerate(path.read_text().splitlines(), 1):
        if line.strip() == 'NONBON':
            active = True
            continue
        if active and line.strip():
            row = line.split()
            result[row[0]] = {'rmin_half_angstrom': float(row[1]),
                              'epsilon_kcal_mol': float(row[2]), 'line': line_no}
    return result


def native_lj(atom):
    return {'rmin_half_angstrom': atom.rmin, 'epsilon_kcal_mol': atom.epsilon,
            'sigma_nm': atom.sigma / 10, 'epsilon_kj_mol': atom.epsilon * 4.184}


def component(atom):
    name = atom.residue.name
    if name == 'GDP': return 'GDP'
    if name == 'COV': return 'covalent_drug_residue'
    if name == 'WAT': return 'TIP3P_water'
    if atom.atomic_number == 12: return 'magnesium'
    if name == 'Na+': return 'sodium_counterion'
    return 'protein'


def run(root, amber, evidence, destination):
    if destination.exists():
        raise ValueError('Refusing to overwrite an existing inventory')
    project = Path(__file__).resolve().parents[4]
    fixed = {
        root / '6oim-complex-v3/solvated.prmtop':
            'd4a3c152a0a39d0c5a5dce8903178f0823a8da164d7a5bade438475c505f3f6c',
        project / 'docs/audit/advanced-chemistry/covalent-panel.json':
            '5fdb1dba8cfa43e90c3396c87bde1cdc2a1a8a7091ec813cb120255ad70a613a',
        project / 'docs/audit/advanced-chemistry/metal-panel.json':
            '28db88c10ce5a80dde8513d9c02542c2aa7ffda4910d8bbaef64ce9f312af694',
        evidence / 'paper.pdf':
            '5965ca3f43699d5ae5a9a237233ea336e52211cf7d9e186ee39ff25dbb29d923',
        evidence / 'supporting-information.pdf':
            'f25f3e1d064ecbf858b42490680d6605b41a10d9d663b4c36d1ffd59f229bf9b',
        evidence / 'gromacs-readir.cpp':
            'a32d477ac66f319dafe586921288940e7fb35f80f9d97cdf8e7f0aae2cf739ef',
        evidence / 'gromacs-mdp-options.rst':
            '49f9136b06c99b507018e3ede93d88e82f11dffc0b565db14e6e1710ed627faf',
    }
    for path, expected in fixed.items():
        if digest(path) != expected:
            raise ValueError(f'Pinned input changed: {path}')

    pythonlib = Path(parmed.__file__).parents[1]
    used = {
        'polarizabilities': amber / 'dat/leap/parm/lj_1264_pol.dat',
        'mg_tip4pew_1264': amber / 'dat/leap/parm/frcmod.ions234lm_1264_tip4pew',
        'water_parameters': amber / 'dat/leap/parm/frcmod.tip4pew',
        'water_template': amber / 'dat/leap/lib/tip4pewbox.off',
        'water_leaprc': amber / 'dat/leap/cmd/leaprc.water.tip4pew',
        'jc_sodium': amber / 'dat/leap/parm/frcmod.ionsjc_tip4pew',
        'lm_1264_sodium': amber / 'dat/leap/parm/frcmod.ions1lm_1264_tip4pew',
        'amber_license': amber / 'LICENSE',
        'amber_package': amber / 'conda-meta/ambertools-24.8-cuda_None_nompi_py314h5cdaa45_104.json',
        'parmed_package': amber / 'conda-meta/parmed-4.3.1-py314h4dda2b3_1.json',
        'amber_add1264': amber / 'lib/python3.14/site-packages/parmed/tools/add1264.py',
        'research_openmm_parser': pythonlib / 'openmm/app/internal/amber_file_parser.py',
        'research_parmed_gromacs_writer': pythonlib / 'parmed/gromacs/gromacstop.py',
        'research_parmed_amber': pythonlib / 'parmed/amber/_amberparm.py',
        'assembly': root / '6oim-complex-v3/result.json',
        'mapping': root / '6oim-complex-v3/source-mapping.json',
        'gdp_model': root / 'gdp-scoped-v6/result.json',
        'paper_text': evidence / 'paper.txt',
        'supporting_text': evidence / 'supporting-information.txt',
        'si_metadata': evidence / 'figshare-item.json',
        'si_retrieval': evidence / 'si-retrieval.json',
        'script': Path(__file__).resolve(),
    }
    originals = {str(p): digest(p) for p in set(used.values()) | set(fixed)}
    package = json.loads(used['amber_package'].read_text())
    paths = {x['_path']: x for x in package['paths_data']['paths']}
    # Verify parameter files against the installed conda package's own manifest.
    manifest_checks = {}
    for key in ('polarizabilities', 'mg_tip4pew_1264', 'water_parameters',
                'water_template', 'water_leaprc', 'jc_sodium', 'lm_1264_sodium'):
        relative = str(used[key].relative_to(amber))
        expected = paths[relative].get('sha256_in_prefix', paths[relative]['sha256'])
        manifest_checks[key] = digest(used[key]) == expected
        if not manifest_checks[key]:
            raise ValueError(f'Installed parameter differs from package manifest: {key}')

    pol = {}
    for n, line in enumerate(used['polarizabilities'].read_text().splitlines(), 1):
        row = line.split(maxsplit=2)
        if not row: continue
        if row[0] in pol: raise ValueError('Duplicate polarizability type')
        pol[row[0]] = {'alpha_angstrom3': float(row[1]), 'line': n,
                       'attribution': row[2] if len(row) == 3 else ''}
    tree = ast.parse(used['amber_add1264'].read_text())
    constants = {node.targets[0].id: ast.literal_eval(node.value)
                 for node in tree.body if isinstance(node, ast.Assign)
                 and isinstance(node.targets[0], ast.Name)
                 and node.targets[0].id in ('DEFAULT_C4_PARAMS', 'WATER_POL')}
    c4water = constants['DEFAULT_C4_PARAMS']['TIP4PEW']['Mg2']
    alpha_water = constants['WATER_POL']
    assert c4water == 180.5 and alpha_water == 1.444
    # Source values must agree independently between the SI, main table and Amber.
    supporting = used['supporting_text'].read_text()
    sirow = re.search(r'Mg2\+\s+0\.569\s+1\.090\s+1\.090\s+0\.170\s+1\.910\s+1\.925', supporting)
    if sirow is None: raise ValueError('SI Table S1 Mg row not reproduced')
    fitted = []
    paper = used['paper_text'].read_text()
    for site, label, key, alpha, radius, epsilon, old, new in [
        ('dimethyl_phosphate_nonbridging_O', 'dimethyl phosphate', 'OPMG', .170, 3.0972, 68.53801, 71.125, 21.25),
        ('adenosine_N7', 'adenosine', 'NAMG', 1.910, 3.2600, 61.66607, 136.25, 238.75),
        ('guanosine_N7', 'guanosine', 'NGMG', 1.925, 3.2600, 61.66607, 136.25, 240.625),
    ]:
        row = re.search(re.escape(label) + r'\s+Mg2\+\s+([\d.]+)\s+([\d.]+)\s+([\d.]+)\s+([\d.]+)', paper)
        if row is None or [float(v) for v in row.groups()] != [radius, epsilon, old, new]:
            raise ValueError('Main Table 2 source row changed')
        assert pol[key]['alpha_angstrom3'] == alpha
        assert abs(c4water / alpha_water * alpha - new) < 1e-10
        fitted.append({'site': site, 'amber_parameter_key': key, **pol[key],
                       'pair_rmin_angstrom': radius, 'pair_epsilon_cal_mol': epsilon,
                       'original_C4_kcal_angstrom4_mol': old,
                       'modified_C4_kcal_angstrom4_mol': new,
                       'modified_C4_kj_nm4_mol': new * 4.184e-4,
                       'source': 'paper Table 2 and SI Table S1', 'assigned_to_native_atoms': False})

    native = parmed.load_file(str(root / '6oim-complex-v3/solvated.prmtop'))
    mapping = {x['native_index']: x for x in json.loads(used['mapping'].read_text())}
    scoped = json.loads(used['gdp_model'].read_text())
    reverse = {v: k for k, v in scoped['type_renames'].items()}
    groups = defaultdict(list)
    classes = defaultdict(list)
    for atom in native.atoms:
        groups[(component(atom), atom.type)].append(atom)
        classes[atom.nb_idx].append(atom)
        if atom.type in reverse: assert component(atom) == 'GDP'
    records = []
    for (part, kind), atoms in sorted(groups.items()):
        source_kind = reverse.get(kind, kind)
        parameter = pol.get(source_kind)
        alpha = None if parameter is None else parameter['alpha_angstrom3']
        lj = native_lj(atoms[0])
        assert all(native_lj(a) == lj for a in atoms)
        records.append({
            'component': part, 'native_type': kind, 'original_type_before_gdp_scoping': source_kind,
            'atom_count': len(atoms), 'atomic_numbers': sorted({a.atomic_number for a in atoms}),
            'charge_range_e': [min(a.charge for a in atoms), max(a.charge for a in atoms)],
            'native_lj': lj, 'native_lj_type_indices': sorted({a.nb_idx for a in atoms}),
            'direct_polarizability_name_available': kind in pol,
            'polarizability_after_documented_alias_resolution': parameter,
            'lookup_status': 'missing_source_entry' if parameter is None else
                             'available_via_documented_GDP_alias' if source_kind != kind else 'available_directly',
            'default_Mg_selected_center_contribution_C4_kcal_angstrom4_mol':
                None if alpha is None else c4water / alpha_water * alpha,
            'default_Mg_selected_center_contribution_C4_kj_nm4_mol':
                None if alpha is None else c4water / alpha_water * alpha * 4.184e-4,
            'C4_scope': 'Unassigned base-rule arithmetic for a TIP4P-Ew Mg candidate; not a fitted site override or a validated transfer.',
            'source_atom_examples': [{'native_index': a.idx, 'atom': a.name,
                                      'source': mapping.get(a.idx, {}).get('source')} for a in atoms[:4]],
        })

    def describe(index):
        atom = native.atoms[index]
        return {'native_index': index, 'atom': atom.name, 'residue': atom.residue.name,
                'native_type': atom.type, 'original_type': reverse.get(atom.type, atom.type),
                'charge_e': atom.charge, 'native_lj_type_index': atom.nb_idx,
                'source': mapping.get(index, {}).get('source'), 'native_lj': native_lj(atom),
                'bonded_atom_names': [a.name for a in atom.bond_partners]}

    critical = {label: describe(index) for label, index in [
        ('Mg', 2775), ('SER17_OG', 323), ('GDP_O2B', 2778), ('GDP_N7', 2796),
        ('water_405', 2835), ('water_415', 2865), ('water_440', 2940), ('water_447', 2961)]}
    class_records = []
    for index, atoms in sorted(classes.items()):
        kinds = sorted({a.type for a in atoms})
        values = {pol[reverse.get(k, k)]['alpha_angstrom3'] for k in kinds if reverse.get(k, k) in pol}
        missing = [k for k in kinds if reverse.get(k, k) not in pol]
        class_records.append({'native_lj_type_index': index, 'native_types': kinds,
                              'atom_count': len(atoms), 'source_alpha_values': sorted(values),
                              'missing_lookup_types': missing,
                              'all_defaults_available_and_equal': not missing and len(values) == 1})
    selectivity = []
    for label, indices, fitted_key in [
        ('GDP_nonbridging_phosphate_analogy_only', [2776, 2778, 2779, 2782, 2783], 'OPMG'),
        ('GDP_guanosine_N7_analogy_only', [2796], 'NGMG')]:
        selected = set(indices)
        class_ids = {native.atoms[i].nb_idx for i in indices}
        affected = [a for c in class_ids for a in classes[c]]
        selectivity.append({'candidate_site': label, 'native_indices': indices,
                            'published_analogue_key': fitted_key,
                            'shared_native_lj_type_indices': sorted(class_ids),
                            'other_atoms_affected_if_class_modified': len([a for a in affected if a.idx not in selected]),
                            'shared_native_types': sorted({a.type for a in affected}),
                            'separate_nonbonded_type_index_needed_for_selective_pair': any(a.idx not in selected for a in affected),
                            'chemical_applicability_validated': False, 'assignment_made': False})

    waterlines = used['water_template'].read_text().splitlines()
    first = next(i for i, line in enumerate(waterlines) if '.unit.atoms table' in line) + 1
    water_atoms = []
    for line in waterlines[first:first + 4]:
        row = shlex.split(line)
        water_atoms.append({'name': row[0], 'type': row[1], 'atomic_number': int(row[6]), 'charge_e': float(row[7])})
    assert [a['charge_e'] for a in water_atoms] == [0.0, .524220, .524220, -1.048440]
    assert '0.1250' in used['water_parameters'].read_text()
    lic = used['amber_license'].read_text()
    assert 'force field parameter files' in lic and 'public domain by their authors' in lic
    metadata = json.loads(used['si_metadata'].read_text())
    assert metadata['doi'] == '10.1021/acs.jpcb.5b10423.s001'
    parser = used['research_openmm_parser'].read_text()
    assert "mm.CustomNonbondedForce('-c/r^4; c=ccoef(type1, type2)')" in parser
    assert 'cforce.setUseLongRangeCorrection(True)' in parser
    writer = used['research_parmed_gromacs_writer'].read_text()
    assert '[ nonbond_params ]' in writer
    assert not re.search(r'CCOEF|has_1264', writer)
    gmx_source = (evidence / 'gromacs-readir.cpp').read_text()
    assert 'With Verlet lists only cut-off and PME LJ interactions are supported' in gmx_source
    gmx_options = (evidence / 'gromacs-mdp-options.rst').read_text()
    assert re.search(r'mdp-value:: User\s+Currently unsupported\.', gmx_options)

    counts = Counter(component(a) for a in native.atoms)
    types = {a.type for a in native.atoms}
    missing_types = sorted(t for t in types if reverse.get(t, t) not in pol)
    assert len(native.atoms) == 27477 and len(types) == 72
    assert missing_types == ['c6'] and sum(a.type == 'c6' for a in native.atoms) == 4
    assert 'LENNARD_JONES_CCOEF' not in native.parm_data
    assembly = json.loads(used['assembly'].read_text())
    report = {
        'schema_version': 1, 'stage': 'panteva_modified_12_6_4_complete_complex_applicability',
        'case': 'corrected_6OIM_KRAS_sotorasib_GDP_Mg_water_counterions',
        'parameters_assigned_or_changed': False, 'new_md_or_qm': False,
        'app_ready': False, 'physical_model_validated': False,
        'complete_complex_atom_counts': dict(counts),
        'native_unique_type_count': len(types),
        'native_C4_term_present': False,
        'source_rule_coverage': {
            'direct_name_types': sum(t in pol for t in types),
            'GDP_aliases_explicitly_resolved_for_inspection_only': reverse,
            'types_after_documented_alias_resolution': len(types) - len(missing_types),
            'missing_types': missing_types, 'missing_type_atom_count': 4,
            'no_element_or_similar_type_guess': True,
            'missing_specific_fitted_override_does_not_mean_missing_default_rule': True,
        },
        'component_type_inventory': records, 'retained_coordination_site': critical,
        'all_native_nonbonded_classes': class_records,
        'site_specific_type_separation': selectivity,
        'published_model': {
            'doi': '10.1021/acs.jpcb.5b10423',
            'primary_paper_url': 'https://theory.rutgers.edu/resources/pdfs/Panteva_JChemPhysB_2015_v119_p15460.pdf',
            'si_doi': metadata['doi'], 'si_download': metadata['files'][0]['download_url'],
            'potential': 'A/r^12 - B/r^6 - C4/r^4 + Coulomb',
            'default_C4_rule': 'C4(Mg,j)=C4(Mg,water)*alpha(j)/alpha(water)',
            'Mg_water_C4_kcal_angstrom4_mol': c4water,
            'water_reference_alpha_angstrom3': alpha_water,
            'Mg_lj_from_Amber_TIP4PEw_1264_file': nonbond_rows(used['mg_tip4pew_1264'])['Mg2+'],
            'fitted_sites': fitted,
            'validation_scope': 'DMP nonbridging phosphate and adenosine/guanosine N7 binding, with TIP4P-Ew; RNA application. Not complete KRAS-GDP-protein-covalent-adduct validation.',
            'modified_TIP3P_site_coefficients_in_this_source': 'not_provided_or_validated',
            'base_TIP3P_Mg_water_C4_in_Amber_for_comparison_only': constants['DEFAULT_C4_PARAMS']['TIP3P']['Mg2'],
        },
        'water_applicability': {
            'current_family': 'TIP3P', 'current_water_count': counts['TIP3P_water'] // 3,
            'current_atoms_per_water': 3, 'current_extra_point_count': sum(a.atomic_number == 0 for a in native.atoms if component(a) == 'TIP3P_water'),
            'required_source_family': 'TIP4P-Ew', 'source_atoms_per_water': 4,
            'source_water_atoms': water_atoms, 'source_water_lj': nonbond_rows(used['water_parameters']),
            'source_O_extra_point_distance_angstrom': .125,
            'current_family_is_directly_compatible': False,
            'required_change': 'Rebuild/retype all retained and bulk waters with correct virtual sites, charges and geometry; Mg C4/radius replacement alone cannot reproduce this model.',
        },
        'GDP_applicability': {
            'current_model': scoped['source_model'],
            'net_charge_e': sum(a.charge for a in native.atoms if component(a) == 'GDP'),
            'default_O3_alpha_angstrom3': pol['O3']['alpha_angstrom3'],
            'default_O3_C4_kcal_angstrom4_mol': c4water / alpha_water * pol['O3']['alpha_angstrom3'],
            'DMP_fitted_C4_kcal_angstrom4_mol': 21.25,
            'GDP_diphosphate_specific_fit_in_inspected_source': 'not_provided',
            'whole_complex_transfer_validated': False,
            'required_decision': 'Do not substitute the DMP correction for the retained GDP(-3) diphosphate without an explicit transfer assessment and binding/coordination validation.',
        },
        'sodium_applicability': {
            'count': counts['sodium_counterion'], 'current': describe(next(a.idx for a in native.atoms if component(a) == 'sodium_counterion')),
            'TIP4PEw_Joung_Cheatham_source': nonbond_rows(used['jc_sodium'])['Na+'],
            'TIP4PEw_Li_Merz_1264_alternative_not_selected': nonbond_rows(used['lm_1264_sodium'])['Na+'],
            'Mg_only_selected_center_default_C4_kcal_angstrom4_mol': c4water / alpha_water * pol['Na+']['alpha_angstrom3'],
            'Na_as_second_1264_center_would_add_reverse_C4_kcal_angstrom4_mol': constants['DEFAULT_C4_PARAMS']['TIP4PEW']['Na1'] / alpha_water * pol['Mg2+']['alpha_angstrom3'],
            'counterion_policy_selected': False,
            'note': 'SI cites Joung-Cheatham NaCl. Keep a consistent source water/ion set; selecting Na as a second 12-6-4 center is a separate model choice, not implied by Mg selection.',
        },
        'engine_representation': {
            'OpenMM': {'inspected_research_version': '8.5.1', 'C4_parser_path_exists': True,
                       'required_prmtop_flag': 'LENNARD_JONES_CCOEF',
                       'C4_kcal_A4_to_kJ_nm4': 4.184e-4,
                       'custom_force': '-c/r^4; c=ccoef(type1, type2)',
                       'pair_matrix_indexing': 'NONBONDED_PARM_INDEX with native nonbonded type indices',
                       'exceptions': 'NonbondedForce exceptions copied as CustomNonbondedForce exclusions',
                       'periodic_C4': 'CutoffPeriodic plus long-range correction; not C4 PME',
                       'execution_test_in_this_assessment': False,
                       'validation_needed': 'Analytic pair energy/force and complete-topology mapping/exclusion/virtual-site/cutoff/tail/virial parity checks before MD.'},
            'GROMACS': {'inspected_version': '2025.4', 'official_source_tag': 'v2025.4',
                        'standard_nonbond_params_preserve_independent_C4': False,
                        'inspected_ParmEd_writer_version': parmed.__version__,
                        'inspected_writer_has_C4_path': False,
                        'User_nonbonded_tables_supported_in_this_version': False,
                        'current_native_conversion_supported': False,
                        'representation_requirement': 'A verified engine extension or supported exact additional pair-force implementation, preserving the C4 table, exclusions, periodic forces and virial. Ordinary C6/C12 .top conversion is insufficient.',
                        'general_impossibility_claim': False,
                        'stale_manual_caveat': 'Reference-manual table overview still describes User tables; the tagged 2025.4 option documentation marks them unsupported and preprocessing rejects them.'},
        },
        'licensing_evidence': {
            'Amber_parameter_files': 'Installed AmberTools LICENSE explicitly places dat/leap force-field parameter files in the public domain.',
            'Amber_code': 'GPL-3.0 with component exceptions; source code license is distinct from parameter-data notice.',
            'SI_license_from_current_Figshare_API': metadata['license'],
            'main_paper': 'ACS copyright 2015; no unrestricted article-copying grant established in this audit.',
            'redistribution_route': 'Use the separately public-domain Amber parameter files with citations and provenance. Do not relabel the CC BY-NC SI or ACS article as unrestricted open-source content.',
            'paper_or_SI_copied_into_project': False,
        },
        'package_provenance': {k: package.get(k) for k in ('name', 'version', 'build', 'url', 'sha256', 'license')},
        'parameter_file_package_manifest_checks': manifest_checks,
        'source_files': {key: {'path': str(path), 'sha256': digest(path)} for key, path in used.items()},
        'input_hashes': originals,
        'limits': ['Source/type inventory only; no parameter changes or physical validation.',
                   'Default polarizability coverage is not validation of GDP, Ser17 or the covalent adduct.',
                   'One missing source polarizability (c6); do not assign zero or borrow c3 silently.',
                   'Research OpenMM parser inspected; packaged application path has not been modified or verified for this candidate.',
                   'Current complete model and both frozen panels remain unchanged.'],
        'next_bounded_step': 'Resolve c6 polarizability provenance and explicitly assess DMP/guanosine transfer to the retained GDP(-3) charges; only then build a separate TIP4P-Ew reference topology for static analytic/native C4 energy-force validation. Keep GROMACS unsupported for this candidate until an exact implementation is demonstrated.',
    }
    assert sum(r['atom_count'] for r in records) == 27477
    # Every original is rehashed after the complete read/analysis pass.
    assert all(digest(Path(p)) == h for p, h in originals.items())
    report['all_original_inputs_unchanged_after_audit'] = True
    destination.write_text(json.dumps(report, indent=2, sort_keys=True) + '\n')
    snapshots = evidence / 'local-source-snapshots'
    snapshots.mkdir(exist_ok=True)
    for key, path in used.items():
        if path.is_relative_to(amber) or key.startswith('research_'):
            target = snapshots / (key + path.suffix)
            if target.exists():
                if digest(target) != digest(path): raise ValueError('Existing source snapshot differs')
            else:
                shutil.copyfile(path, target)
    print(json.dumps({'output': str(destination), 'sha256': digest(destination),
                      'atoms': len(native.atoms), 'native_types': len(types),
                      'covered_after_alias_resolution': len(types) - len(missing_types),
                      'missing_types': missing_types, 'shared_site_classes': selectivity,
                      'inputs_unchanged': True}, indent=2))


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--root', type=Path, required=True)
    parser.add_argument('--amber', type=Path, required=True)
    parser.add_argument('--evidence', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    run(args.root, args.amber, args.evidence, args.output)
