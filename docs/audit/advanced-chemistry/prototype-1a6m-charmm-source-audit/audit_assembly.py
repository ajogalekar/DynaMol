"""Audit published CHARMM PSF and prepare immutable-identity 1A6M inputs; no MD/QM."""
from pathlib import Path
from collections import Counter, defaultdict
import hashlib
import importlib.util
import json
import re
import warnings

import gemmi

HERE = Path(__file__).resolve().parent
AUDIT = HERE.parent
PUBLISHED = AUDIT / 'non-zinc-models/dryad-oxy-myoglobin'


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def save(name, value):
    (HERE / name).write_text(json.dumps(value, indent=2) + '\n')


def blocks(path):
    result = {}
    current = None
    for line in path.read_text().splitlines():
        fields = line.split('!', 1)[0].split()
        if not fields:
            continue
        if fields[0] in {'RESI', 'PRES'}:
            current = {'kind': fields[0], 'atoms': {}, 'delete_atoms': [], 'records': []}
            result[fields[1]] = current
        elif current is not None:
            current['records'].append(fields)
            if fields[0] == 'ATOM':
                current['atoms'][fields[1]] = {'type': fields[2], 'charge': float(fields[3])}
            elif fields[:2] == ['DELETE', 'ATOM']:
                current['delete_atoms'].extend(fields[2:])
    return result


case = next(c for c in json.loads((AUDIT / 'metal-panel.json').read_text())['cases'] if c['id'] == '1A6M')
source = AUDIT / case['source']['path']
assert sha(source) == case['source']['sha256']
old = blocks(PUBLISHED / 'top_heme.inp')
new = blocks(PUBLISHED / 'top_heme_oxy_optimized_charges.inp')
spec = importlib.util.spec_from_file_location('source_inspection', AUDIT / 'non-zinc-models/inspect_sources.py')
inspection = importlib.util.module_from_spec(spec)
spec.loader.exec_module(inspection)
atoms, terms = inspection.psf(PUBLISHED / 'horse_oxy_only_wbi_10ns.psf')
by_residue = defaultdict(dict)
for i, atom in atoms.items():
    by_residue[(atom['segment'], atom['resid'], atom['residue'])][atom['name']] = dict(atom, psf_index_1based=i)

# The deposited PSF declares these terminal patches explicitly in its title.
terminal_patches = {('A', '1'): 'GLYP', ('A', '153'): 'CTER'}
mismatches = []
for (segment, resid, resname), residue in by_residue.items():
    if resname not in old:
        mismatches.append({'identity': [segment, resid, resname], 'error': 'No supplied residue template'})
        continue
    expected = dict(old[resname]['atoms'])
    patch = terminal_patches.get((segment, resid))
    if patch:
        expected.update(old[patch]['atoms'])
        for name in old[patch]['delete_atoms']:
            expected.pop(name, None)
    missing = sorted(set(expected) - set(residue))
    extra = sorted(set(residue) - set(expected))
    changed = {name: {'expected': expected[name], 'actual': residue[name]} for name in set(expected) & set(residue)
               if expected[name]['type'] != residue[name]['type'] or abs(expected[name]['charge'] - residue[name]['charge']) > 1e-7}
    if missing or extra or changed:
        mismatches.append({'identity': [segment, resid, resname], 'missing': missing, 'extra': extra, 'changed': changed})

site_indices = {i for i, a in atoms.items() if a['residue'] in {'HEM', 'OXY'} or (a['segment'], a['resid']) == ('A', '93')}
site_terms = {kind: [[dict(index_1based=i, **atoms[i]) for i in t] for t in values if set(t) & site_indices]
              for kind, values in terms.items()}
save('published-site-terms.json', site_terms)
save('published-template-mapping.json', {
    'source_psf_sha256': sha(PUBLISHED / 'horse_oxy_only_wbi_10ns.psf'),
    'atoms': len(atoms), 'residues': len(by_residue),
    'residue_counts': dict(Counter(k[2] for k in by_residue)),
    'term_counts': {k: len(v) for k, v in terms.items()},
    'terminal_patches': [{'segment': k[0], 'resid': k[1], 'patch': v} for k, v in terminal_patches.items()],
    'template_mismatches': mismatches,
    'all_atoms_match_supplied_baseline_templates_after_declared_terminal_patches': not mismatches,
    'published_pdb_atom_count': sum(l.startswith(('ATOM ', 'HETATM')) for l in (PUBLISHED / 'horse_oxy.pdb').read_text().splitlines()),
    'caveat': 'The supplied horse PDB is unsolvated and is not a coordinate file for the full 27386-atom saved PSF. No coordinates are silently assigned by index.',
})

# Freeze whole-residue alternate choices, preserving rejected conformer identities.
table = gemmi.cif.read_file(str(source)).sole_block().find_mmcif_category('_atom_site.')
tags = [t.split('.', 1)[1] for t in table.tags]
rows = [dict(zip(tags, row)) for row in table]
groups = defaultdict(list)
for row in rows:
    assert row['pdbx_PDB_model_num'] == '1'
    assert row['label_asym_id'] in case['assembly_choice']['label_chains']
    key = (row['label_asym_id'], row['auth_asym_id'], row['auth_seq_id'], row['pdbx_PDB_ins_code'], row['label_comp_id'])
    groups[key].append(row)
selected, excluded, choices = [], [], []
for key, group in groups.items():
    scores = defaultdict(float)
    for r in group:
        if r['label_alt_id'] not in {'.', '?'}:
            scores[r['label_alt_id']] += float(r['occupancy'])
    choice = sorted(scores, key=lambda a: (-scores[a], a))[0] if scores else None
    kept = [r for r in group if r['label_alt_id'] in {'.', '?', choice}]
    dropped = [r for r in group if r not in kept]
    assert len({r['label_atom_id'] for r in kept}) == len(kept)
    selected.extend(kept)
    excluded.extend(dropped)
    if choice:
        choices.append({'residue_identity': key, 'selected_altloc': choice, 'summed_occupancies': dict(scores), 'excluded_atom_site_ids': [r['id'] for r in dropped]})
selected.sort(key=lambda r: int(r['id']))
identity_map = []
for index, r in enumerate(selected, 1):
    name, resname = r['auth_atom_id'], r['auth_comp_id']
    target_segment = {'A': 'PROA', 'B': 'SUL1', 'C': 'SUL2', 'D': 'HEME', 'E': 'OXY', 'F': 'WAT'}[r['label_asym_id']]
    target_name = 'OH2' if resname == 'HOH' else ('CD' if resname == 'ILE' and name == 'CD1' else name)
    if r['label_asym_id'] == 'A' and r['auth_seq_id'] == '151':
        target_name = {'O': 'OT1', 'OXT': 'OT2'}.get(name, target_name)
    target_resname = 'TIP3' if resname == 'HOH' else ('HSD' if resname == 'HIS' else resname)
    identity_map.append({'selected_index_1based': index, 'source_atom_site_id': r['id'],
        'source': r, 'proposed_charmm_identity': {'segment': target_segment, 'resid': r['auth_seq_id'], 'residue': target_resname, 'name': target_name},
        'histidine_state_status': ('proximal neutral ND1-protonated, NE2 donor as published PHEM model' if r['auth_seq_id'] == '93' else 'provisional published all-HIS-to-HSD policy; assess protonation before runnable 1A6M build') if resname == 'HIS' else None})
assert len(selected) == 1445 and len(excluded) == 139
assert Counter(r['label_comp_id'] for r in selected)['SO4'] == 10
assert Counter(r['label_comp_id'] for r in selected)['HOH'] == 186
save('1a6m-identity-map.json', identity_map)
save('1a6m-alternate-selection.json', {'policy': 'Whole-residue highest summed occupancy; lexical altloc tie-break; shared blank atoms retained. Alternate records remain in original CIF and this audit.', 'choices': choices, 'excluded_alternate_records': excluded})
save('1a6m-assembly-plan.json', {
    'status': 'Not runnable: additive all-atom sulfate model remains unresolved; no component removed.',
    'source': case['source'], 'assembly_choice': case['assembly_choice'],
    'source_atom_rows': len(rows), 'selected_heavy_atoms': len(selected), 'excluded_alternate_rows': len(excluded),
    'selected_residue_counts': dict(Counter(k[4] for k in groups)),
    'model_variants': [
        {'id': 'dryad-2015-baseline', 'topology': 'top_heme.inp', 'topology_sha256': sha(PUBLISHED / 'top_heme.inp'), 'published_saved_psf_matches': True},
        {'id': 'dryad-2015-updated-charges', 'topology': 'top_heme_oxy_optimized_charges.inp', 'topology_sha256': sha(PUBLISHED / 'top_heme_oxy_optimized_charges.inp'), 'published_saved_psf_matches': False}],
    'parameter_file': {'path': str(PUBLISHED.relative_to(AUDIT) / 'par_all27.inp'), 'sha256': sha(PUBLISHED / 'par_all27.inp')},
    'patch_order': ['PHEM PROA:93 HEME:154', 'PLO2 OXY:157 HEME:154 PROA:93'],
    'protein_terminal_patches': {'PROA:1:VAL': 'NTER', 'PROA:151:TYR': 'CTER'},
    'do_not_copy_horse_terminal_patch': 'Saved horse PSF uses GLYP on GLY1; 1A6M begins VAL1 and requires its own NTER template.',
    'additional_deposited_Fe_O2_contact': 'Retain as source coordination evidence; published PLO2 creates only Fe-O1. No second Fe-O2 bond is generated.',
    'water_model': 'Published additive CHARMM TIP3 including hydrogen Lennard-Jones terms; retain all 186 source water oxygens and account for every added H.',
    'sulfate_obligations': [r for r in case['retained_nonstandard_residues'] if r['component'] == 'SO4'],
    'protonation_obligations': ['Explicitly assess nonproximal His states; proposed all-HSD alias mirrors deposited script, not an independent pKa determination.', 'HEM propionates and two SO4 groups retain declared dianion model states; partial Fe charge is not oxidation/spin evidence.'],
    'validation_remaining': ['Complete matched sulfate topology, bonded terms and nonbonded parameters with provenance and redistribution terms.', 'Native full 1A6M PSF/coordinate build preserving every selected observed heavy atom.', 'Native/static energy-force and exception parity for entire CHARMM Hamiltonian.', 'Independent whole-site angular/torsion validation; no accuracy conclusion from zero-angle terms or intact bonds alone.'],
})

# Build mechanical definitions for the complete published horse PSF, never a
# sulfate-deleted 1A6M. This resolves parameters but evaluates no energy or motion.
from openmm import XmlSerializer, NonbondedForce, unit
from openmm.app import CharmmParameterSet, CharmmPsfFile, NoCutoff
summaries = []
for variant, topfile in [('baseline', 'top_heme.inp'), ('updated-charges', 'top_heme_oxy_optimized_charges.inp')]:
    text = (PUBLISHED / 'horse_oxy_only_wbi_10ns.psf').read_text()
    modifications = []
    if variant == 'updated-charges':
        lines = text.splitlines(keepends=True)
        start = next(i for i, line in enumerate(lines) if '!NATOM' in line) + 1
        for offset in range(len(atoms)):
            i = start + offset
            f = lines[i].split()
            if f[3] in {'HEM', 'OXY'}:
                expected = new[f[3]]['atoms'][f[4]]
                assert f[5] == expected['type']
                if float(f[6]) != expected['charge']:
                    modifications.append({'index_1based': int(f[0]), 'identity': f[1:5], 'old_charge_e': float(f[6]), 'new_charge_e': expected['charge']})
                    matches = list(re.finditer(r'\S+', lines[i])); m = matches[6]
                    lines[i] = lines[i][:m.start()] + f"{expected['charge']:.6f}" + lines[i][m.end():]
        text = ''.join(lines)
    output = HERE / ('published-horse-' + variant)
    output.mkdir(exist_ok=True)
    psf_path = output / 'model.psf'
    psf_path.write_text(text)
    with warnings.catch_warnings(record=True) as captured:
        params = CharmmParameterSet(str(PUBLISHED / topfile), str(PUBLISHED / 'par_all27.inp'))
        psf = CharmmPsfFile(str(psf_path))
        system = psf.createSystem(params, nonbondedMethod=NoCutoff, constraints=None, rigidWater=False)
    (output / 'system.xml').write_text(XmlSerializer.serialize(system))
    nb = next(f for f in system.getForces() if isinstance(f, NonbondedForce))
    charges = [nb.getParticleParameters(i)[0].value_in_unit(unit.elementary_charge) for i in range(nb.getNumParticles())]
    expected_atoms, expected_terms = inspection.psf(psf_path)
    assert expected_terms == terms
    assert all(abs(charges[i - 1] - a['charge']) < 1e-10 for i, a in expected_atoms.items())
    atom_section = re.search(r'^\s*(\d+)\s+!NATOM[^\n]*\n', text, re.M)
    native_masses = [float(line.split()[7]) for line in text[atom_section.end():].strip().splitlines()[:len(atoms)]]
    mass_error = max(abs(system.getParticleMass(i).value_in_unit(unit.dalton) - mass) for i, mass in enumerate(native_masses))
    assert mass_error < 1e-6
    summary = {'variant': variant, 'source_psf_sha256': sha(PUBLISHED / 'horse_oxy_only_wbi_10ns.psf'),
               'output_psf_sha256': sha(psf_path), 'system_xml_sha256': sha(output / 'system.xml'),
               'particle_count': system.getNumParticles(), 'force_classes': [type(f).__name__ for f in system.getForces()],
               'constraint_count': system.getNumConstraints(), 'nonbonded_exceptions': nb.getNumExceptions(),
               'total_partial_charge_e': sum(charges), 'particle_charges_match_psf': True, 'all_psf_term_indices_unchanged': True,
               'maximum_particle_mass_error_dalton': mass_error,
               'charge_modifications': modifications, 'reader_warnings': [str(w.message) for w in captured],
               'scope': 'Complete published horse PSF mechanical definition only; no coordinates, Context, energy, QM, MD or 1A6M parameterization acceptance.'}
    save(str(output.relative_to(HERE) / 'audit.json'), summary)
    summaries.append(summary)
save('published-native-system-audit.json', summaries)
print(json.dumps({'selected_1a6m_heavy_atoms': len(selected), 'sulfate_heavy_atoms_retained': 10,
                  'source_waters_retained': 186, 'template_mismatch_count': len(mismatches),
                  'published_psf_variants': len(summaries), '1a6m_runnable': False}))
