"""Audit source integrity and existing sampler artifacts without further sampling."""
import hashlib
import json
import os
from pathlib import Path
import sys

os.environ['PYTHONDONTWRITEBYTECODE']='1'
HERE=Path(__file__).resolve().parent
sys.path.insert(0,str(HERE))
from evaluate_private import SOURCE_CACHE, SOURCE, CASES, sha

CACHE=SOURCE_CACHE.parent/'disgro-evaluation-v2'
sys.path.insert(0,str(SOURCE))
from pydisgro import Structure
from pydisgro.structure import blank_loop


def main():
    manifest=json.loads((SOURCE_CACHE/'source-manifest.json').read_text())
    mismatches=[path for path,wanted in manifest.items() if sha(SOURCE_CACHE/path)!=wanted]
    if mismatches:
        raise ValueError(f'Third-party source changed: {mismatches}')
    source_lines=(CACHE/'1UA2/original-context.pdb').read_text().splitlines(keepends=True)
    probe=next(line for line in source_lines if line.startswith('ATOM  ') and line[12:16].strip()=='CA')
    inserted=probe[:26]+'A'+probe[27:]
    ordinary=Structure.readPdb([probe])
    insertion=Structure.readPdb([inserted])
    shorter_sequence=blank_loop(source_lines,43,56,'A','A')
    actual_placeholders=[line for line in shorter_sequence if line.startswith('ATOM  ') and 44<=int(line[22:26])<=55 and line[21]=='A']
    header='''Source inspection is not a C++/Python numerical parity validation. The synthetic input below exists only to expose parser behavior.'''
    report={'source_archive_unmodified':True,'source_manifest_sha256':sha(SOURCE_CACHE/'source-manifest.json'),
            'third_party_source_hash_mismatches':mismatches,'scope':header,
            'probe_insertion_code':{'ordinary_populated_atom_count':sum(bool(a._name) for r in ordinary._res[1:] for a in r._atom),
                                    'insertion_A_populated_atom_count':sum(bool(a._name) for r in insertion._res[1:] for a in r._atom),
                                    'insertion_silently_omitted':not any(a._name for r in insertion._res[1:] for a in r._atom),
                                    'empty_structure_default_numRes':insertion.numRes},
            'probe_short_sequence':{'requested_gap_length':12,'provided_sequence_length':1,'rejected':False,'placeholder_count':len(actual_placeholders)},
            'adapter_short_sequence_behavior':'Exact length required before sampler call',
            'upstream_provenance':json.loads((SOURCE_CACHE/'upstream-provenance/comparison.json').read_text()),
            'cases':[]}
    for case in CASES:
        folder=CACHE/case
        result=json.loads((folder/'result.json').read_text())
        parser=json.loads((folder/'parser-audit.json').read_text())
        raw=(folder/'candidate-000-RAW-LOSSY-PARSER-OUTPUT.pdb').read_text().splitlines()
        atomlines=[line for line in raw if line.startswith(('ATOM  ','HETATM'))]
        coordinates=lambda line:[float(line[p:p+8]) for p in (30,38,46)]
        raw_audit={'raw_atom_count':len(atomlines),'raw_heterogen_atom_count':sum(line.startswith('HETATM') for line in atomlines),
                   'chain_ids':sorted({line[21] for line in atomlines}),
                   'zero_coordinate_heavy_atom_count':sum(all(x==0 for x in coordinates(line)) and line[12:16].strip()!='H' for line in atomlines),
                   'CONECT_count':sum(line.startswith('CONECT') for line in raw),
                   'raw_output_can_replace_complex':False}
        candidates=result['candidates']
        summary={'case_id':case,'source_sha256':parser['source_sha256'],
                 'execution':json.loads((folder/'execution.json').read_text()),
                 'loop':CASES[case],'source_atom_count':parser['source_atom_count'],
                 'omitted_by_residue':parser['omitted_by_residue'],'sampler_environment_complete':parser['sampler_environment_complete'],
                 'trials':result['trial_counts'],'stored_proposals':result['sampler_stored'],'candidate_count':len(candidates),
                 'coherent_gross_geometry_pass':sum(c['coherent_gross_geometry_pass'] for c in candidates),
                 'original_anchor_restored_geometry_pass':sum(c['restored_gross_geometry_pass'] for c in candidates),
                 'coherent_rama_pass':sum(c['coherent_rama_pass'] for c in candidates),
                 'within_original_1A_anchor_cap':sum(c['maximum_anchor_displacement_angstrom']<=1 for c in candidates),
                 'anchor_shift_range_angstrom':[min(c['maximum_anchor_displacement_angstrom'] for c in candidates),max(c['maximum_anchor_displacement_angstrom'] for c in candidates)],
                 'candidate_paths':[c['path'] for c in candidates], 'raw_writer_audit':raw_audit,
                 'passed_whole_preparation':False}
        report['cases'].append(summary)
    report['prior_failed_attempt']={'path':str(SOURCE_CACHE),'stage':'Adapter import before any sampling','error':'SMC is in pydisgro.smc, not package __init__','retained':True,'corrected_only_in_private_adapter':True}
    report['app_ready']=False
    report['inherited_redistribution_clearance']=False
    report['new_sampling_in_this_audit']=False
    (HERE/'result-summary.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps({'source_unchanged':not mismatches,'probe_insertion':report['probe_insertion_code'],'probe_sequence':report['probe_short_sequence'],'cases':[{'id':c['case_id'],'raw_writer':c['raw_writer_audit']} for c in report['cases']]},indent=2))


if __name__=='__main__':
    main()
