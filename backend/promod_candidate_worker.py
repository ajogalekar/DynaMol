#!/usr/bin/env python3
"""Opt-in native proposals; complete-complex acceptance belongs to the caller.

This worker starts from the observed PDB and original sequence gaps on every
invocation. It never reads previous candidate coordinates. Native libraries are
loaded only in ``run``. The parent must supervise time, memory and cancellation.
CLI exit codes: 0 candidate, 2 explicitly unavailable, 1 malformed/native failure.
"""
from __future__ import annotations

import hashlib
import importlib.util
import json
import math
import os
from pathlib import Path
import re
import sys
import time
from types import SimpleNamespace

if __package__:
    from . import promod_loop_worker as legacy
    from .loop_search import CandidateUnavailable
else:
    # ``python -I /absolute/worker.py`` deliberately has no application imports.
    spec = importlib.util.spec_from_file_location('_dynamol_promod_contract', Path(__file__).with_name('promod_loop_worker.py'))
    legacy = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(legacy)

    class CandidateUnavailable(Exception):
        """A bounded native search produced no complete proposal."""


_SIDECHAINS = {
    'ALA': 'CB', 'ARG': 'CB CG CD NE CZ NH1 NH2', 'ASN': 'CB CG OD1 ND2',
    'ASP': 'CB CG OD1 OD2', 'CYS': 'CB SG', 'GLN': 'CB CG CD OE1 NE2',
    'GLU': 'CB CG CD OE1 OE2', 'GLY': '', 'HIS': 'CB CG ND1 CD2 CE1 NE2',
    'ILE': 'CB CG1 CG2 CD1', 'LEU': 'CB CG CD1 CD2', 'LYS': 'CB CG CD CE NZ',
    'MET': 'CB CG SD CE', 'PHE': 'CB CG CD1 CD2 CE1 CE2 CZ', 'PRO': 'CB CG CD',
    'SER': 'CB OG', 'THR': 'CB OG1 CG2', 'TRP': 'CB CG CD1 CD2 NE1 CE2 CE3 CZ2 CZ3 CH2',
    'TYR': 'CB CG CD1 CD2 CE1 CE2 CZ OH', 'VAL': 'CB CG1 CG2',
}
_BACKBONE = ('N', 'CA', 'C', 'O')
_STEPS = 5000
_MAX_SOURCE_BYTES = 64 * 1024 * 1024


def _key(name, row):
    return (name, row['resid'], row['insertion_code'], row['residue'])


def _sha(path):
    with Path(path).open('rb') as handle:
        return hashlib.file_digest(handle, 'sha256').hexdigest()


def validate_request(config):
    # The opt-in proposal worker can search a small amount of extra native
    # context; only the original missing residues and immediate flanks are
    # exported. The app's legacy construction policy is unchanged.
    extension = config.get('max_res_extension', 0)
    if type(extension) is not int or not 0 <= extension <= 2:
        raise ValueError('Native construction may extend by at most two context residues.')
    if extension and config.get('strategy') != 'fragment':
        raise ValueError('Extended native context currently requires the fragment strategy.')
    legacy.validate_input({**config, 'max_res_extension': 0})
    if config.get('strategy') not in {'fragment', 'torsion_mc', 'fragment_mc'}:
        raise ValueError('An explicit fragment, torsion_mc or fragment_mc strategy is required.')
    if type(config.get('seed')) is not int or not 0 <= config['seed'] <= 2**31 - 1:
        raise ValueError('An explicit nonnegative 31-bit integer seed is required.')
    if 'source_sha256' in config and not re.fullmatch(r'[0-9a-f]{64}', str(config['source_sha256'])):
        raise ValueError('Expected source SHA256 is malformed.')
    for chain in config['chains']:
        numbers = set()
        for row in chain['residues']:
            if not re.fullmatch(r'-?\d+', row['resid']) or str(int(row['resid'])) != row['resid'] or not -999 <= int(row['resid']) <= 9999:
                raise ValueError('Author residue numbers cannot be represented exactly in the source PDB.')
            if not re.fullmatch(r'[A-Za-z0-9]?', row['insertion_code']):
                raise ValueError('Author insertion code is not canonical.')
            pair = row['resid'], row['insertion_code']
            if pair in numbers:
                raise ValueError('Ambiguous author residue number/insertion code.')
            numbers.add(pair)
            if row['one_letter'] not in legacy._LETTERS.values():
                raise ValueError('Sequence requires a canonical parent amino-acid letter.')
            if row['residue'] in legacy._LETTERS and legacy._LETTERS[row['residue']] != row['one_letter']:
                raise ValueError('Standard residue and supplied sequence letter conflict.')
    return config


def read_source(path, chains):
    """Parse exact heavy-atom identities before any native standardization."""
    path = Path(path)
    if path.stat().st_size > _MAX_SOURCE_BYTES:
        raise ValueError('Observed source PDB exceeds the bounded input size.')
    atoms, identities, models = {}, {}, 0
    for line in path.read_text().splitlines():
        if line.startswith('MODEL '):
            models += 1
            if models > 1:
                raise ValueError('The observed source must contain exactly one coordinate model.')
        if line[:6].strip() not in {'ATOM', 'HETATM'}:
            continue
        if len(line) < 78 or line[16].strip():
            raise ValueError('Observed source has ambiguous alternates or missing element fields.')
        key = (line[21].strip(), line[22:26].strip(), line[26].strip(), line[17:20].strip(), line[12:16].strip())
        if any(not key[i] for i in (1, 3, 4)) or key in identities:
            raise ValueError('Observed source has missing or duplicate atom identities.')
        element = line[76:78].strip().upper()
        if not re.fullmatch(r'[A-Z]{1,2}', element):
            raise ValueError('Observed source atom element is malformed.')
        xyz = [float(line[start:start + 8]) for start in (30, 38, 46)]
        if not all(math.isfinite(v) for v in xyz):
            raise ValueError('Observed source coordinates are nonfinite.')
        identities[key] = element
        if element not in {'H', 'D'}:
            atoms[key] = {'element': element, 'xyz_angstrom': xyz}
    if not atoms:
        raise ValueError('Observed source contains no heavy atoms.')
    for chain in chains:
        actual = list(dict.fromkeys(k[:4] for k in atoms if k[0] == chain['chain_id']))
        expected = [_key(chain['chain_id'], row) for row in chain['residues'] if row['observed']]
        if actual != expected:
            raise ValueError('Observed source identities/order do not match requested alignment.')
    return atoms


def gap_specs(chains):
    """Identify original stretches and disallow ambiguous multi-gap composition."""
    gaps, occupied, outer_support = [], set(), set()
    for chain in chains:
        name, rows = chain['chain_id'], chain['residues']
        index = 0
        while index < len(rows):
            if rows[index]['observed']:
                index += 1
                continue
            first = index
            while index < len(rows) and not rows[index]['observed']:
                index += 1
            # Numbering below is the native 1-based alignment position.
            start, end = first, index + 1
            stems = {_key(name, rows[start - 1]), _key(name, rows[end - 1])}
            support = {_key(name, rows[i]) for i in (start - 2, end) if 0 <= i < len(rows)}
            if occupied & (stems | support) or outer_support & stems:
                raise CandidateUnavailable('Requested gaps share or conflict with observed stem support; no independent coherent composition is available.')
            if any(key[3] not in _SIDECHAINS for key in stems):
                raise CandidateUnavailable('An immediate loop stem is modified; standardized chemistry cannot be exported.')
            occupied.update(stems)
            outer_support.update(support)
            gaps.append({'chain': name, 'start': start, 'end': end,
                         'sequence': ''.join(row['one_letter'] for row in rows[start - 1:end]),
                         'modeled': [_key(name, rows[i]) for i in range(first, index)],
                         'context': sorted(stems)})
    return gaps


def _load_native():
    import ost
    import promod3
    from ost import io, seq, geom
    from promod3 import modelling, loop
    return SimpleNamespace(ost=ost, io=io, seq=seq, geom=geom, modelling=modelling, loop=loop,
                           versions={'promod3': promod3.__version__, 'ost': ost.__version__})


def _point(point):
    return [point.x, point.y, point.z]


def _omega_degrees(points):
    """Finite, nondegenerate vector-projection dihedral without native imports."""
    if any(not math.isfinite(v) for point in points for v in point):
        return None
    def dot(a,b): return sum(x*y for x,y in zip(a,b))
    def sub(a,b): return [x-y for x,y in zip(a,b)]
    axis=sub(points[2],points[1]); length=math.sqrt(dot(axis,axis))
    if length<1e-12:
        return None
    axis=[v/length for v in axis]
    left,right=sub(points[0],points[1]),sub(points[3],points[2])
    left=sub(left,[dot(left,axis)*v for v in axis])
    right=sub(right,[dot(right,axis)*v for v in axis])
    if min(dot(left,left),dot(right,right))<1e-24:
        return None
    cross=[axis[1]*left[2]-axis[2]*left[1],axis[2]*left[0]-axis[0]*left[2],axis[0]*left[1]-axis[1]*left[0]]
    return math.degrees(math.atan2(dot(cross,right),dot(left,right)))


def new_peptide_policy(candidate, gap, chain_name, rows):
    """Align fragment selection with the existing missing-link trans policy.

    Observed--observed peptides retain their source state elsewhere. X--Pro
    links keep their candidate basin. This does not establish an unknown
    native isomer state, and does not replace final geometry validation.
    """
    start=gap.before.number.num
    checked=[]
    for offset in range(1,len(candidate)):
        left,right=rows[start+offset-1],rows[start+offset]
        if (left['observed'] and right['observed']) or right['one_letter']=='P':
            continue
        points=[candidate.GetCA(offset-1),candidate.GetC(offset-1),
                candidate.GetN(offset),candidate.GetCA(offset)]
        omega=_omega_degrees([_point(p) for p in points])
        compatible=omega is not None and abs(omega)>90+1e-6
        checked.append({'atoms':[[*_key(chain_name,row),atom] for row,atom in
                                 ((left,'CA'),(left,'C'),(right,'N'),(right,'CA'))],
                        'omega_degrees':omega,'compatible':compatible,
                        'reason':None if compatible else 'undefined' if omega is None else
                                 'ambiguous_basin' if abs(abs(omega)-90)<=1e-6 else 'cis_basin'})
    return {'compatible':all(row['compatible'] for row in checked),'new_nonpro_links':checked,
            'policy':'New non-Pro links start in the trans basin; observed links and X-Pro candidate states are preserved separately.'}


def _raw_model(config, native):
    entity = native.io.LoadPDB(config['input_pdb'])
    alignments, maps, parsed_coordinates = native.seq.AlignmentList(), {}, {}
    for item in config['chains']:
        name, rows = item['chain_id'], item['residues']
        view = entity.Select('cname=' + name)
        observed = [row for row in rows if row['observed']]
        actual = [(str(r.number.num), str(r.number.ins_code).strip().replace('\x00', ''), r.name) for r in view.residues]
        if actual != [(r['resid'], r['insertion_code'], r['residue']) for r in observed]:
            raise ValueError('Native parser changed exact observed residue identities.')
        if [r.one_letter_code for r in view.residues] != [r['one_letter'] for r in observed]:
            raise ValueError('Native observed sequence differs from source alignment.')
        for residue, row in zip(view.residues, observed):
            key = _key(name, row)
            for atom in residue.atoms:
                identity = (*key, atom.name)
                if identity in parsed_coordinates:
                    raise ValueError('Native parser introduced duplicate observed atom identities.')
                parsed_coordinates[identity] = _point(atom.pos)
        alignment = native.seq.CreateAlignment(
            native.seq.CreateSequence('target', ''.join(r['one_letter'] for r in rows)),
            native.seq.CreateSequence('observed', ''.join(r['one_letter'] if r['observed'] else '-' for r in rows)))
        alignment.AttachView(1, view)
        alignments.append(alignment)
        maps[name] = {i + 1: row for i, row in enumerate(rows)}
    model = native.modelling.BuildRawModel(alignments, chain_names=list(maps), include_ligands=True, aln_preprocessing=False)
    expected = {(name, pos) for name, rows in maps.items() for pos, row in rows.items() if not row['observed']}
    absent = {(name, pos) for name, rows in maps.items() for pos in rows
              if not model.model.FindResidue(name, native.ost.mol.ResNum(pos)).IsValid()}
    if absent != expected:
        raise ValueError('Native raw model changed the requested missing-residue inventory.')
    return model, maps, parsed_coordinates


def source_boundaries(spec, maps, source, geom):
    """Read source phi/psi and neighbor identities, without terminal defaults."""
    name, start, end = spec['chain'], spec['start'], spec['end']
    rows = maps[name]
    if start - 1 not in rows or end + 1 not in rows:
        raise CandidateUnavailable('Source-conditioned sampling needs observed residues outside both original stems.')
    if any(not rows[pos]['observed'] for pos in (start - 1, start, end, end + 1)):
        raise CandidateUnavailable('Source boundary conditioning would cross another unresolved gap.')
    def point(pos, atom):
        key = (*_key(name, rows[pos]), atom)
        if key not in source:
            raise ValueError('Observed boundary support atom is missing: ' + repr(key))
        return geom.Vec3(*source[key]['xyz_angstrom'])
    values = {
        'n_stem_phi': geom.DihedralAngle(point(start - 1, 'C'), point(start, 'N'), point(start, 'CA'), point(start, 'C')),
        'c_stem_psi': geom.DihedralAngle(point(end, 'N'), point(end, 'CA'), point(end, 'C'), point(end + 1, 'N')),
        'prev_aa': rows[start - 1]['one_letter'], 'next_aa': rows[end + 1]['one_letter'],
    }
    if not all(math.isfinite(values[key]) for key in ('n_stem_phi', 'c_stem_psi')):
        raise ValueError('Observed boundary torsions are nonfinite.')
    return values


def _native_gaps(model, specs):
    found = {}
    for ref in model.gaps:
        gap = ref.Copy()  # Insertion invalidates native gap-vector references.
        if gap.IsTerminal():
            raise ValueError('Native raw model introduced a terminal gap.')
        key = model.model.chains[gap.GetChainIndex()].name, gap.before.number.num, gap.after.number.num
        if key in found:
            raise ValueError('Native model contains duplicate gaps.')
        found[key] = gap
    expected = {(s['chain'], s['start'], s['end']) for s in specs}
    if set(found) != expected:
        raise ValueError('Native gap spans differ from the original requested stretches.')
    for spec in specs:
        if found[(spec['chain'], spec['start'], spec['end'])].full_seq != spec['sequence']:
            raise ValueError('Native gap sequence differs from the exact requested alignment.')
    return found


def _assert_source_stems(model, specs, maps, source, native, parsed_coordinates):
    for spec in specs:
        for pos in (spec['start'], spec['end']):
            residue = model.model.FindResidue(spec['chain'], native.ost.mol.ResNum(pos))
            key = _key(spec['chain'], maps[spec['chain']][pos])
            if not residue.IsValid() or residue.name != key[3]:
                raise ValueError('Native scoring changed a canonical original stem identity.')
            by_name = {a.name: a for a in residue.atoms}
            for name in _BACKBONE:
                if (*key, name) not in source or name not in by_name:
                    raise ValueError('An observed original stem lacks required backbone support.')
                identity = (*key, name)
                parsed = parsed_coordinates.get(identity)
                original = source[identity]['xyz_angstrom']
                # Native PDB parsing uses finite precision; compare the original
                # three-decimal file values at file precision, then require exact
                # preservation of the parsed coordinates through BuildRawModel.
                # The separate final 1 A source-displacement cap is unchanged.
                if (parsed is None or not all(math.isfinite(v) for v in parsed)
                        or [f'{v:.3f}' for v in parsed] != [f'{v:.3f}' for v in original]):
                    raise ValueError('Native parser changed original observed PDB coordinates.')
                if _point(by_name[name].pos) != parsed:
                    raise ValueError('Native raw model changed original observed stem coordinates.')


def _sample_gap(model, gap, spec, maps, source, native, strategy, seed):
    modelling, loop = native.modelling, native.loop
    conditioning = source_boundaries(spec, maps, source, native.geom)
    torsions = loop.LoadTorsionSamplerCoil(seed=seed)
    if strategy == 'torsion_mc':
        sampler = modelling.PhiPsiSampler(gap.full_seq, torsions, seed=seed, **conditioning)
        details = {'sampler': 'PhiPsiSampler', 'sampler_conditioning': 'source phi/psi and neighbor letters'}
    else:
        length = min(3, len(gap.full_seq))
        handler = modelling.FraggerHandle(model.seqres[gap.GetChainIndex()], fragment_length=length,
                                          fragments_per_position=100, structure_db=loop.LoadStructureDB())
        fraggers = handler.GetList(spec['start'] - 1, spec['end'] - 1)
        if not fraggers or any(len(fragger) == 0 for fragger in fraggers):
            raise CandidateUnavailable('Sequence-based fragment sampling found no complete fragment inventory.')
        sampler = modelling.FragmentSampler(gap.full_seq, fraggers, init_fragments=5, seed=seed)
        details = {'sampler': 'FragmentSampler', 'fragment_length': length,
                   'sampler_conditioning': 'sequence-based fragments; no phi/psi keyword API',
                   'fragments_per_position': 100, 'init_fragments': 5}
    closer = modelling.CCDCloser(gap.before, gap.after, gap.full_seq, torsions, seed)
    weights = modelling.ScoringWeights.GetWeights()
    scorer = modelling.LinearScorer(model.backbone_scorer, model.backbone_scorer_env,
        spec['start'], len(gap.full_seq), gap.GetChainIndex(), {key: weights[key] for key in ('reduced', 'cb_packing', 'clash')})
    cooler = modelling.ExponentialCooler(round(_STEPS / 109), 100, .9)
    try:
        candidates = modelling.LoopCandidates.FillFromMonteCarloSampler(gap.full_seq, 1, _STEPS, sampler, closer, scorer, cooler, seed)
    except RuntimeError as exc:
        # ProMod3 documents this as a failure to find an initial closed loop.
        # It is an unavailable candidate, not a malformed input or native crash.
        # Only this exact error is recoverable; all other runtime errors remain fatal.
        if str(exc) != 'Failed to initialize monte carlo sampling protocol!':
            raise
        raise CandidateUnavailable(
            f'Native {strategy} could not initialize a closed loop for chain '
            f'{spec["chain"]}, span {spec["start"]}-{spec["end"]}, seed {seed}.') from exc
    if not len(candidates):
        raise CandidateUnavailable('Bounded native Monte Carlo returned no closed candidate.')
    if len(candidates) != 1:
        raise ValueError('Native Monte Carlo returned an unexpected candidate count.')
    backbone = candidates[0].Copy()
    if len(backbone) != len(gap.full_seq):
        raise ValueError('Native candidate changed the original gap span.')
    score = scorer.GetScore(backbone)
    if not math.isfinite(score):
        raise ValueError('Native candidate score is nonfinite.')
    modelling.InsertLoopClearGaps(model, backbone, gap)
    return dict(details, chain=spec['chain'], span=[spec['start'], spec['end']], seed=seed,
                steps=_STEPS, source_boundary_conditioning=conditioning, closer='CCDCloser',
                native_backbone_score=score, original_stems_only=True)


def export_candidate(model, maps, specs, source):
    """Export complete standard loops plus coherent immediate observed flanks."""
    wanted = {key for spec in specs for key in spec['modeled']}
    context = {key for spec in specs for key in spec['context']}
    modeled, observed, seen = [], [], set()
    for chain in model.model.chains:
        if chain.name not in maps:
            continue  # Native ligand/modified chemistry is never an authority.
        for residue in chain.residues:
            row = maps[chain.name].get(residue.number.num)
            if row is None:
                raise ValueError('Native model introduced an unmapped protein residue.')
            key = _key(chain.name, row)
            if key not in wanted | context:
                continue
            if key[3] not in _SIDECHAINS or residue.name != key[3]:
                raise ValueError('Native residue chemistry differs from a requested standard identity.')
            for atom in residue.atoms:
                element = atom.element.upper()
                if element in {'H', 'D'}:
                    continue
                identity = (*key, atom.name)
                if identity in seen:
                    raise ValueError('Native candidate has duplicate atom identities.')
                xyz = _point(atom.pos)
                if not all(math.isfinite(value) for value in xyz):
                    raise ValueError('Native candidate coordinates are nonfinite.')
                if key in context:
                    if identity not in source or source[identity]['element'] != element:
                        raise ValueError('Native context atom differs from exact observed identity/element.')
                elif element != atom.name[0]:
                    raise ValueError('Modeled standard atom element differs from its chemical identity.')
                seen.add(identity)
                record = {'identity': list(identity), 'element': element, 'xyz_nm': [v / 10 for v in xyz]}
                (modeled if key in wanted else observed).append(record)
    expected_modeled = {(*key, atom) for key in wanted for atom in (*_BACKBONE, *_SIDECHAINS[key[3]].split())}
    expected_context = {key for key in source if key[:4] in context}
    if {tuple(row['identity']) for row in modeled} != expected_modeled:
        raise ValueError('Native candidate does not contain every requested missing heavy atom.')
    if {tuple(row['identity']) for row in observed} != expected_context:
        raise ValueError('Native candidate does not contain complete coherent observed context.')
    for key in context:
        if not {(*key, name) for name in _BACKBONE} <= expected_context:
            raise ValueError('Observed context lacks a complete backbone.')
    return sorted(modeled, key=lambda row: row['identity']), sorted(observed, key=lambda row: row['identity'])


def run(config):
    validate_request(config)
    source_path = Path(config['input_pdb'])
    source_sha = _sha(source_path)
    if config.get('source_sha256', source_sha) != source_sha:
        raise ValueError('Observed source SHA256 differs from the frozen request.')
    started = time.monotonic()
    try:
        source = read_source(source_path, config['chains'])
        specs = gap_specs(config['chains'])
        native = _load_native()
        model, maps, parsed_coordinates = _raw_model(config, native)
        gaps = _native_gaps(model, specs)
        _assert_source_stems(model, specs, maps, source, native, parsed_coordinates)
        native.modelling.SetupDefaultBackboneScoring(model)
        if config['strategy'] == 'fragment':
            source_positions = {_key(name, row): (name, pos)
                                for name, rows in maps.items() for pos, row in rows.items() if row['observed']}
            observed = {(*source_positions[key[:4]], key[4]): value['xyz_angstrom']
                        for key, value in source.items() if key[:4] in source_positions}
            torsions = native.loop.LoadTorsionSamplerCoil(seed=config['seed'])
            attempts = legacy.context_fragment_search(model, native.modelling, native.loop.LoadFragDB(),
                native.loop.LoadStructureDB(), torsions, config.get('max_res_extension', 0), observed, include_carbonyl_anchors=True,
                candidate_policy=lambda candidate,gap: new_peptide_policy(candidate,gap,
                    model.model.chains[gap.GetChainIndex()].name,
                    maps[model.model.chains[gap.GetChainIndex()].name]))
            if model.gaps:
                unavailable=CandidateUnavailable('Bounded fragment search did not find a compatible closed candidate for every requested gap.')
                unavailable.attempts=attempts
                raise unavailable
        else:
            attempts = [_sample_gap(model, gaps[(spec['chain'], spec['start'], spec['end'])], spec,
                        maps, source, native, config['strategy'], config['seed']) for spec in specs]
            if model.gaps:
                raise ValueError('Native insertion left unexpected gaps after complete sampling.')
        native.modelling.ReconstructSidechains(model.model, keep_sidechains=True, build_disulfids=False, consider_ligands=True)
        modeled, context = export_candidate(model, maps, specs, source)
        return {'schema_version': 1, 'status': 'candidate', 'source_sha256': source_sha,
                'strategy': config['strategy'], 'seed': config['seed'], 'versions': native.versions,
                'modeled_residue_keys': sorted([list(key) for spec in specs for key in spec['modeled']]),
                'modeled_heavy_atoms': modeled, 'observed_context_atoms': context,
                'context_residue_keys': sorted([list(key) for spec in specs for key in spec['context']]),
                'attempts': attempts, 'original_model_reset': True,
                'max_res_extension': config.get('max_res_extension', 0),
                'export_scope': 'Original requested missing residues and immediate observed flanks only; additional search context is not transplanted.',
                'sampling_budget': {'gap_count': len(specs),
                                    'mc_steps_per_gap': _STEPS if config['strategy'] != 'fragment' else 0,
                                    'total_mc_steps_maximum': _STEPS * len(specs) if config['strategy'] != 'fragment' else 0,
                                    'external_process_supervision_required': True},
                'elapsed_seconds': time.monotonic() - started,
                'full_complex_scoring_claimed': False, 'physical_model_validated': False,
                'context_policy': 'Complete immediate standard flanks exported coherently; retained modified/nonprotein chemistry remains solely in the source snapshot.',
                'randomness': {'requested_seed': config['seed'], 'seed_applied_to_torsion_sampler': True,
                               'seed_applied_to_mc_sampler_closer_and_acceptance': config['strategy'] != 'fragment',
                               'bitwise_cross_platform_reproducibility_claimed': False}}
    finally:
        if _sha(source_path) != source_sha:
            raise ValueError('Observed source changed during native generation.')


def main(argv=None):
    args = sys.argv[1:] if argv is None else argv
    if len(args) != 2:
        raise SystemExit('Usage: promod_candidate_worker.py INPUT.json OUTPUT.json')
    request, output = map(Path, args)
    if output.exists():
        raise ValueError('Candidate output already exists; use a fresh attempt directory.')
    os.environ.setdefault('PM3_OPENMM_CPU_THREADS', '2')
    try:
        if request.stat().st_size > 16 * 1024 * 1024:
            raise ValueError('Candidate request exceeds its bounded input size.')
        result, code = run(json.loads(request.read_text())), 0
    except CandidateUnavailable as exc:
        result, code = {'schema_version': 1, 'status': 'unavailable', 'error': str(exc), 'error_type': type(exc).__name__}, 2
        if hasattr(exc,'attempts'):
            result['attempts']=exc.attempts
    except Exception as exc:
        result, code = {'schema_version': 1, 'status': 'failed', 'error': str(exc), 'error_type': type(exc).__name__}, 1
    output.parent.mkdir(parents=True, exist_ok=True)
    legacy._write_json(output, result)
    return code


if __name__ == '__main__':
    raise SystemExit(main())
