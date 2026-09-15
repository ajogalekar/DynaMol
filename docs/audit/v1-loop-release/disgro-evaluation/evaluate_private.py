"""Private, bounded evaluation adapter; no third-party source or model vendored.

Coordinates are proposals only. Preserve the original complete complex verbatim;
never use the sampler's PDB round trip as a molecular identity authority.
"""
import argparse
import collections
import hashlib
import json
import os
from pathlib import Path
import resource
import subprocess
import sys
import time
import traceback

ROOT = Path(__file__).resolve().parents[4]
SOURCE_CACHE = Path('/Users/ashujo/.cache/dynamol-research/v1-loop-release/disgro-evaluation-v1')
CACHE = Path(os.environ.get('DYNAMOL_DISGRO_EVAL_CACHE', str(SOURCE_CACHE)))
SOURCE = SOURCE_CACHE/'source/pydisgro-0.1.0'
PYTHON = Path('/Users/ashujo/.cache/dynamol-runtimes/covalent-v1-focused/bin/python')
RAMA = Path('/Users/ashujo/.cache/dynamol-runtimes/cctbx-rama-reference-v1')
ARCHIVE = Path('/Users/ashujo/.cache/dynamol-research/v1-loop-release/alternative-methods-review-20260914/pydisgro-0.1.0.tar.gz')
CASES = {'1UA2': {'first': 44, 'last': 55, 'sequence': 'KLGHRSEAKDGI', 'seed': 91401},
         '8K5R': {'first': 177, 'last': 179, 'sequence': 'AKN', 'seed': 91402}}
BUDGET = {'wall_seconds': 180, 'resident_bytes': 4*1024**3, 'threads_per_worker': 1,
          'concurrent_workers': 1, 'trials': 256, 'distance_states': 32, 'confkeep': 20,
          'sample_sidechains': False}


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def write(path, value):
    Path(path).write_text(json.dumps(value, indent=2, allow_nan=False)+'\n')


def pdb_rows(lines):
    rows = []
    for line in lines:
        if line[:6] not in {'ATOM  ', 'HETATM'}:
            continue
        rows.append({'identity': [line[21], line[22:26].strip(), line[26].strip(), line[17:20].strip(), line[12:16].strip()],
                     'xyz_angstrom': [float(line[p:p+8]) for p in (30, 38, 46)],
                     'element': line[76:78].strip() or line[12:16].strip()[0],
                     'record': line[:6].strip(), 'serial': int(line[6:11]), 'altloc': line[16].strip()})
    if len({tuple(r['identity']) for r in rows}) != len(rows):
        raise ValueError('Input atom identities are ambiguous; adapter does not resolve alternate conformers')
    return rows


def atom_row(residue, slot):
    from pydisgro.residue import Residue
    name = str(Residue.cType[residue._type][slot])
    return {'identity': [residue._chainName, str(residue._pdbIndex), '', Residue.Name3[residue._type], name],
            'xyz_angstrom': [float(x) for x in residue._atom[slot].xyz], 'element': name[0]}


def geometry_screen(rows):
    # Reuse an existing independent gross-geometry screen unchanged.
    from screen_kic_backbones import geometry
    return geometry(rows)


def rama_screen(rows, first, last):
    import numpy as np
    from openmm import app, unit
    from rama_reference import rama_report
    topology = app.Topology()
    chain = topology.addChain('A')
    groups = collections.defaultdict(list)
    for row in rows:
        groups[tuple(row['identity'][:4])].append(row)
    xyz, previous, selected = [], None, []
    for identity, atoms in sorted(groups.items(), key=lambda kv:int(kv[0][1])):
        residue = topology.addResidue(identity[3], chain, identity[1], identity[2])
        lookup = {}
        for row in atoms:
            name = row['identity'][4]
            lookup[name] = topology.addAtom(name, app.element.get_by_symbol(row['element']), residue)
            xyz.append(row['xyz_angstrom'])
        for a, b in [('N', 'CA'), ('CA', 'C'), ('C', 'O')]:
            topology.addBond(lookup[a], lookup[b])
        if previous is not None:
            topology.addBond(previous['C'], lookup['N'])
        previous = lookup
        if first-1 <= int(identity[1]) <= last+1:
            selected.append(identity)
    return rama_report(topology, np.asarray(xyz)*unit.angstrom, selected)


def gross_contacts(generated, context, coherent):
    """Diagnostic heavy-atom overlaps, not a replacement for native FF exclusions."""
    import numpy as np
    scaffold = {tuple(r['identity']): r for r in context}
    scaffold.update({tuple(r['identity']): r for r in coherent})
    fixed = [r for r in scaffold.values() if r['element'] not in {'H', 'D'}]
    lookup = {tuple(r['identity']): np.asarray(r['xyz_angstrom']) for r in fixed+generated}
    # Graph-distance exclusions for canonical backbone and C-beta only.
    graph = collections.defaultdict(set)
    groups = collections.defaultdict(set)
    for key in lookup:
        groups[key[:4]].add(key[4])
    for group, names in groups.items():
        for a, b in [('N', 'CA'), ('CA', 'C'), ('C', 'O'), ('CA', 'CB')]:
            if a in names and b in names:
                x, y = (*group, a), (*group, b)
                graph[x].add(y); graph[y].add(x)
    by_number = {(k[0], int(k[1])): k for k in groups if k[1].lstrip('-').isdigit() and not k[2]}
    for (chain, number), group in by_number.items():
        other = by_number.get((chain, number+1))
        if other and 'C' in groups[group] and 'N' in groups[other]:
            a, b = (*group, 'C'), (*other, 'N')
            graph[a].add(b); graph[b].add(a)
    contacts = []
    fixedxyz = np.asarray([r['xyz_angstrom'] for r in fixed])
    for row in generated:
        key = tuple(row['identity'])
        excluded = {key}|graph[key]
        for neighbor in list(graph[key]):
            excluded |= graph[neighbor]
        distances = np.linalg.norm(fixedxyz-np.asarray(row['xyz_angstrom']), axis=1)
        for index in np.where(distances < 1.8)[0]:
            target = fixed[index]
            if tuple(target['identity']) not in excluded:
                contacts.append({'generated': row['identity'], 'source': target['identity'],
                                 'distance_angstrom': float(distances[index]), 'source_record': target.get('record')})
    return {'threshold_angstrom': 1.8, 'count': len(contacts), 'contacts': contacts,
            'scope': 'Gross overlaps against every retained source heavy atom; local backbone 1-2/1-3 excluded. Not a comprehensive van der Waals or chirality gate.'}


def worker(case_id):
    import numpy as np
    sys.path[:0] = [str(SOURCE), str(RAMA), str(ROOT/'docs/audit/v1-loop-release/native-sampling'),
                   str(ROOT/'docs/audit/advanced-chemistry/covalent-v1-focused')]
    from pydisgro import Structure
    from pydisgro.smc import SMC
    from pydisgro.structure import blank_loop, resnum_to_index
    from pydisgro.geom import seed
    case = CASES[case_id]
    out = CACHE/case_id
    lines = (out/'original-context.pdb').read_text().splitlines(keepends=True)
    source_rows = pdb_rows(lines)
    first, last = case['first'], case['last']
    if len(case['sequence']) != last-first+1:
        raise ValueError('Exact missing sequence/range mismatch')
    if any(r['identity'][0] == 'A' and first <= int(r['identity'][1]) <= last for r in source_rows):
        raise ValueError('Requested missing residue already exists in immutable context')
    blanked = blank_loop(lines, first-1, last+1, case['sequence'], 'A')
    (out/'sampler-input-NOT-FULL-COMPLEX-AUTHORITY.pdb').write_text(''.join(blanked))
    conf = Structure.readPdb(blanked)
    start, end = [resnum_to_index(conf, n, 'A') for n in (first-1, last+1)]
    if end-start != last-first+2:
        raise ValueError('Sampler range does not map to the exact intended sequence')
    parsed_keys = {tuple(atom_row(r, j)['identity']) for r in conf._res[1:] for j, a in enumerate(r._atom) if a._name}
    omitted = [r for r in source_rows if tuple(r['identity']) not in parsed_keys]
    parser_audit = {'source_atom_count': len(source_rows), 'source_heterogen_atom_count': sum(r['record']=='HETATM' for r in source_rows),
                    'source_sha256': sha(out/'original-context.pdb'), 'omitted_source_atoms': omitted,
                    'omitted_by_residue': dict(collections.Counter(':'.join(r['identity'][:4]) for r in omitted)),
                    'sampler_environment_complete': not omitted, 'sampler_loop_indices': [start,end],
                    'residue_mapping': [{'sampler_index':i, 'chain':r._chainName,'resid':r._pdbIndex} for i,r in enumerate(conf._res) if i]}
    write(out/'parser-audit.json', parser_audit)
    seed(case['seed'])
    sampler = SMC(conf, start, end, num_conf=BUDGET['trials'], num_distance_states=(32,),
                  confkeep=20, sample_sc=False, verbose=True)
    original_trial = sampler.trial
    counts = collections.Counter()
    progress = (out/'trial-events.jsonl').open('w', buffering=1)
    def logged_trial():
        began = time.monotonic()
        work, state = original_trial()
        counts['trials'] += 1
        counts['growth_returned'] += work is not None
        counts['closure_claimed'] += bool(state.closed)
        progress.write(json.dumps({'trial':counts['trials'], 'seconds':time.monotonic()-began,
                                   'growth_returned':work is not None,'closure_claimed':bool(state.closed)})+'\n')
        if counts['trials'] % 16 == 0:
            print(json.dumps(dict(counts)),flush=True)
        return work, state
    sampler.trial = logged_trial  # instrumentation only; no coordinate/scoring change
    results = sampler.run()
    progress.close()
    # Preserve all stored proposals, including ones the port's ranker rejected.
    np.savez_compressed(out/'all-stored-proposals.npz', xyz=np.asarray(sampler.LoopStore), energy=np.asarray(sampler.LoopEnergy))
    reports = []
    source_map = {tuple(r['identity']):r for r in source_rows}
    for number, result in enumerate(results):
        rebuilt = sampler.to_structure(result)
        generated, context = [], []
        for residue in rebuilt._res[start:end+1]:
            for slot in range(4):
                row = atom_row(residue, slot)
                if not np.isfinite(row['xyz_angstrom']).all() or np.all(np.asarray(row['xyz_angstrom']) == 0):
                    raise ValueError('Candidate has nonfinite or unplaced backbone atom')
                (generated if first <= residue._pdbIndex <= last else context).append(row)
        expected_count = 4*(last-first+1)
        if len(generated) != expected_count or len({tuple(r['identity']) for r in generated}) != expected_count:
            raise ValueError('Generated backbone identity count mismatch')
        margins = [r for r in source_rows if r['identity'][0]=='A' and int(r['identity'][1]) in (first-2,last+2) and r['identity'][4] in {'N','CA','C','O'}]
        coherent = margins+context+generated
        restored = margins+[source_map[tuple(r['identity'])] for r in context]+generated
        moves = [{'identity':r['identity'],'displacement_angstrom':float(np.linalg.norm(np.asarray(r['xyz_angstrom'])-source_map[tuple(r['identity'])]['xyz_angstrom']))} for r in context]
        candidate = {'schema_version':1, 'case_id':case_id, 'candidate_index':number,
                     'source_sha256':sha(out/'original-context.pdb'), 'immutable_scaffold':str(out/'original-context.pdb'),
                     'loop':dict(chain_id='A',start_residue=first,end_residue=last,sequence=case['sequence'],left_anchor=first-1,right_anchor=last+1),
                     'seed':case['seed'],'budget':BUDGET, 'generated_atoms':generated, 'research_context_atoms':context,
                     'sampler_energy':float(result.energy),'sampler_energy_is_not_force_field_energy':True,
                     'observed_anchor_displacements':moves,'maximum_observed_anchor_displacement_angstrom':max(r['displacement_angstrom'] for r in moves),
                     'coherent_backbone_geometry':geometry_screen(coherent), 'original_anchor_restored_geometry':geometry_screen(restored),
                     'coherent_rama':rama_screen(coherent,first,last), 'original_anchor_restored_rama':rama_screen(restored,first,last),
                     'coherent_complete_scaffold_gross_contacts':gross_contacts(generated,source_rows,context),
                     'original_anchor_restored_complete_scaffold_gross_contacts':gross_contacts(generated,source_rows,[]),
                     'preserved_source_atom_count':len(source_rows),'preserved_source_sha256':sha(out/'original-context.pdb'),
                     'source_scaffold_modified':False,'sampler_environment_complete':not omitted,
                     'scope':'Missing N/CA/C/O coordinate proposal with two mutable anchor proposals kept separately. All original complex atoms retained in immutable scaffold; missing sidechains absent. No full-complex refinement or admission.',
                     'app_ready':False,'whole_preparation_pass':False}
        path = out/f'candidate-{number:03d}.json'
        write(path,candidate)
        rebuilt.writePdb(str(out/f'candidate-{number:03d}-RAW-LOSSY-PARSER-OUTPUT.pdb'), l_start=start,l_end=end)
        reports.append({'path':str(path),'sha256':sha(path),'index':number,
                        'coherent_gross_geometry_pass':candidate['coherent_backbone_geometry']['accepted'],
                        'restored_gross_geometry_pass':candidate['original_anchor_restored_geometry']['accepted'],
                        'coherent_rama_pass':candidate['coherent_rama']['all_selected_scored_without_outliers'],
                        'restored_rama_pass':candidate['original_anchor_restored_rama']['all_selected_scored_without_outliers'],
                        'maximum_anchor_displacement_angstrom':candidate['maximum_observed_anchor_displacement_angstrom'],
                        'coherent_gross_contacts':candidate['coherent_complete_scaffold_gross_contacts']['count'],
                        'restored_gross_contacts':candidate['original_anchor_restored_complete_scaffold_gross_contacts']['count']})
    write(out/'result.json', {'case_id':case_id,'status':'completed','trial_counts':dict(counts),'sampler_closed':sampler.NumClosedconf,
                            'sampler_stored':len(sampler.LoopStore),'candidate_count':len(reports),'candidates':reports,
                            'parser_audit':str(out/'parser-audit.json'),'sampler_environment_complete':not omitted,
                            'app_ready':False,'whole_preparation_pass':False})
    print(json.dumps({'completed':case_id,'candidates':len(reports)}),flush=True)


def supervisor():
    frozen = {'archive':str(ARCHIVE),'archive_sha256':sha(ARCHIVE),'adapter_sha256':sha(__file__),
              'source_manifest_sha256':sha(SOURCE_CACHE/'source-manifest.json'),'cases':CASES,'budget':BUDGET,
              'scope':'Private source evaluation only; no software/model redistribution clearance.'}
    CACHE.mkdir(exist_ok=True)
    write(CACHE/'frozen-protocol.json',frozen)
    env = os.environ.copy()
    for key in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','NUMEXPR_NUM_THREADS','VECLIB_MAXIMUM_THREADS']:
        env[key]='1'
    env['PYTHONDONTWRITEBYTECODE']='1'
    env['PYTHONHASHSEED']='0'
    outputs=[]
    for case_id in CASES:
        out=CACHE/case_id
        out.mkdir(exist_ok=False)
        parent=CACHE.parent/'carbonyl-ranking-v2'/case_id
        for a,b in [('context.pdb','original-context.pdb'),('input.json','original-context-config.json')]:
            (out/b).write_bytes((parent/a).read_bytes())
        began=time.monotonic()
        peak=0
        reason=None
        with (out/'worker.log').open('w') as stream:
            proc=subprocess.Popen([str(PYTHON),'-u',str(Path(__file__).resolve()),'--worker',case_id],stdout=stream,stderr=subprocess.STDOUT,env=env,start_new_session=True)
            while proc.poll() is None:
                if time.monotonic()-began > BUDGET['wall_seconds']:
                    reason='wall_budget_exceeded'
                rss=subprocess.run(['ps','-o','rss=','-p',str(proc.pid)],capture_output=True,text=True)
                current=int(rss.stdout.strip() or '0')*1024
                peak=max(peak,current)
                if current > BUDGET['resident_bytes']:
                    reason='resident_memory_budget_exceeded'
                if reason:
                    import signal
                    os.killpg(proc.pid,signal.SIGKILL)
                    break
                time.sleep(0.2)
            proc.wait()
        record={'case_id':case_id,'exit_code':proc.returncode,'termination_reason':reason,
                'wall_seconds':time.monotonic()-began,'maximum_sampled_RSS_bytes':peak,
                'memory_enforcement':'Supervisor polls resident memory every 0.2 seconds and kills worker process group above 4 GiB; one native BLAS thread. No child jobs invoked by package.',
                'source_preserved':sha(parent/'context.pdb')==sha(out/'original-context.pdb')}
        write(out/'execution.json',record)
        outputs.append(record)
        print(json.dumps(record),flush=True)
    if frozen['adapter_sha256']!=sha(__file__) or frozen['archive_sha256']!=sha(ARCHIVE):
        raise ValueError('Frozen evaluation inputs changed')
    write(CACHE/'execution-summary.json',outputs)


if __name__=='__main__':
    parser=argparse.ArgumentParser()
    parser.add_argument('--worker',choices=CASES)
    args=parser.parse_args()
    if args.worker:
        try:
            worker(args.worker)
        except Exception as error:
            write(CACHE/args.worker/'failure.json',{'type':type(error).__name__,'error':str(error),'traceback':traceback.format_exc()})
            raise
    else:
        supervisor()
