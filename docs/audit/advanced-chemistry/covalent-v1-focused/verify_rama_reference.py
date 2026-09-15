"""Verify the research adapter against the complete CCTBX PDB/ramalyze path."""
import argparse
import json
from pathlib import Path

import iotbx.pdb
from mmtbx.validation.ramalyze import ramalyze

from prepare_adduct import digest


def verify(reports, output):
    checks = []
    sources = {}
    for path in reports:
        report = json.loads(path.read_text())
        sources[str(path)] = digest(path)
        for name, expected in report['source_sha256'].items():
            if digest(Path(name)) != expected:
                raise ValueError('Adapter audit source changed')
        pdbs = [Path(name) for name in report['source_sha256'] if name.endswith('.pdb')]
        if len(pdbs) != 1:
            raise ValueError('Expected one unambiguous PDB source')
        pdb = pdbs[0]
        sources[str(pdb)] = digest(pdb)
        hierarchy = iotbx.pdb.input(file_name=str(pdb)).construct_hierarchy()
        full = ramalyze(pdb_hierarchy=hierarchy, outliers_only=False, quiet=True)
        rows = {}
        for row in full.results:
            key = (row.chain_id, row.resseq.strip(), row.icode.strip(), row.resname)
            if key in rows:
                raise ValueError('Multiple full-reference results for one residue')
            rows[key] = row
        angle_errors, score_errors = [], []
        for own in report['rows']:
            reference = rows[tuple(own['residue'])]
            if own['classification'] != ['outlier', 'allowed', 'favored'][reference.rama_type]:
                raise ValueError('Adapter classification differs from full CCTBX')
            for a, b in [(own['phi_degrees'], reference.phi), (own['psi_degrees'], reference.psi)]:
                angle_errors.append(abs((a-b+180) % 360-180))
            # The full reporting interface prints percentile units.
            score_errors.append(abs(own['score']-reference.score/100.))
        if max(angle_errors) > .001 or max(score_errors) > 1e-5:
            raise ValueError('Adapter angles or scores differ beyond coordinate precision')
        checks.append({'pdb': str(pdb), 'residues_compared': len(report['rows']),
                       'classifications_identical': True,
                       'maximum_angle_difference_degrees': max(angle_errors),
                       'maximum_score_difference': max(score_errors)})
    if sources != {name: digest(Path(name)) for name in sources}:
        raise ValueError('Source changed during reference verification')
    result = {'stage': 'full_cctbx_reference_parity_verified', 'checks': checks,
              'total_residues_compared': sum(row['residues_compared'] for row in checks),
              'source_sha256': sources, 'app_ready': False, 'physical_model_validated': False,
              'scope': 'Independent CCTBX PDB parser and backbone-angle calculation agree with the adapter on saved structures. Both use native CCTBX statistical tables; this verifies integration, not native loop accuracy.'}
    output.mkdir(exist_ok=False)
    (output/'result.json').write_text(json.dumps(result, indent=2)+'\n')
    print(json.dumps({k: v for k, v in result.items() if k != 'source_sha256'}, indent=2))


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--reports', nargs='+', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    verify(args.reports, args.output)
