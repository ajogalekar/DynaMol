"""Research candidate: native AM1-BCC/PREPGEN cysteine-adduct residue.

Uses the documented Amber modified-residue route. PREPGEN cap-charge
redistribution is recorded explicitly; these are not constrained RESP charges.
An assembled topology is not a force-field accuracy or app-readiness claim.
"""
from __future__ import annotations

import argparse
from decimal import Decimal
import hashlib
import json
import os
from pathlib import Path
import re
import signal
import subprocess
import time

import numpy as np
import openmm as mm
from openmm import app, unit
import parmed
import psutil


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def native(amber, folder, arguments, seconds=60, rss_limit=2_000_000_000):
    env = dict(os.environ, AMBERHOME=str(amber), OMP_NUM_THREADS='1',
               OPENBLAS_NUM_THREADS='1', OPENMM_CPU_THREADS='2')
    env['PATH'] = str(amber/'bin') + os.pathsep + env.get('PATH', '')
    log = folder/(arguments[0] + '.log')
    started = time.monotonic()
    peak = 0
    with log.open('w') as stream:
        proc = subprocess.Popen([str(amber/'bin'/arguments[0]), *arguments[1:]],
            cwd=folder, env=env, stdout=stream, stderr=subprocess.STDOUT,
            start_new_session=True)
        try:
            while proc.poll() is None:
                try:
                    parent = psutil.Process(proc.pid)
                    members = [parent, *parent.children(recursive=True)]
                    rss = 0
                    for member in members:
                        try:
                            rss += member.memory_info().rss
                        except psutil.NoSuchProcess:
                            pass
                    peak = max(peak, rss)
                except psutil.NoSuchProcess:
                    pass
                if peak > rss_limit or time.monotonic()-started > seconds:
                    raise TimeoutError('Native resource limit: ' + arguments[0])
                time.sleep(.1)
            if proc.returncode:
                raise RuntimeError('Native tool failed; inspect ' + str(log))
        finally:
            if proc.poll() is None:
                os.killpg(proc.pid, signal.SIGTERM)
                try:
                    proc.wait(timeout=3)
                except subprocess.TimeoutExpired:
                    os.killpg(proc.pid, signal.SIGKILL)
                    proc.wait(timeout=5)
    return {'arguments': arguments, 'seconds': time.monotonic()-started,
            'peak_sampled_rss_bytes': peak, 'rss_limit_bytes': rss_limit,
            'time_limit_seconds': seconds, 'log': log.name}


def mol2_sections(path):
    sections, current = {}, None
    for line in Path(path).read_text().splitlines():
        if line.startswith('@<TRIPOS>'):
            current = line[9:]
            sections[current] = []
        elif current and line.strip():
            sections[current].append(line.split())
    return sections


def validate_native_charge_total(baseline, charges, formal_charge):
    """Admit native SQM print precision only with its charge-conservation evidence.

    Some bundled SQM versions print individual Mulliken charges to 0.001 e.
    Antechamber carries that rounding into BCC output. PREPGEN subsequently
    enforces the declared residue charge; this function never changes charges.
    """
    q = np.asarray(charges, dtype=float)
    if not np.isfinite(q).all():
        raise ValueError('Nonfinite native atomic charge')
    error = float(q.sum()-formal_charge)
    report = {'input_total_e':float(q.sum()), 'declared_total_e':formal_charge,
              'input_error_e':error, 'charges_modified_by_admission':False}
    if abs(error) < 1e-5:
        return dict(report, basis='Total agrees within native output precision')
    folder = baseline/'capped-native'
    sqm = (folder/'sqm.out').read_text()
    if 'Calculation Completed' not in sqm:
        raise ValueError('Incomplete native SQM calculation')
    blocks = re.findall(r'Mulliken Charge\s*\n(.*?)Total Mulliken Charge\s*=\s*([-+\d.]+)',
                        sqm, re.S)
    if not blocks:
        raise ValueError('Missing native Mulliken charge evidence')
    table, total = blocks[-1]
    rows = [line.split() for line in table.splitlines() if line.strip()]
    if len(rows) != len(q) or [int(r[0]) for r in rows] != list(range(1,len(q)+1)):
        raise ValueError('Native Mulliken atom inventory differs')
    printed = [Decimal(r[-1]) for r in rows]
    bound = float(sum(Decimal('0.5')*Decimal(10)**x.as_tuple().exponent for x in printed))
    printed_sum = float(sum(printed))
    if abs(float(total)-formal_charge) > 1e-5 or abs(printed_sum-formal_charge) > bound+1e-8:
        raise ValueError('Native charge total disagrees with declared chemical state')
    intermediate = []
    for filename in ['ANTECHAMBER_AM1BCC_PRE.AC','ANTECHAMBER_AM1BCC.AC']:
        atoms = [line.split() for line in (folder/filename).read_text().splitlines()
                 if line.startswith('ATOM')]
        if len(atoms) != len(q):
            raise ValueError('Native BCC charge inventory differs')
        values = np.array([float(row[-2]) for row in atoms])
        if not np.isfinite(values).all() or abs(values.sum()-printed_sum) > len(q)*0.5e-6+1e-8:
            raise ValueError('Native BCC transformation did not conserve printed charge')
        intermediate.append(float(values.sum()))
    native_q = np.array([float(row[8]) for row in mol2_sections(folder/'charged.mol2')['ATOM']])
    if native_q.shape != q.shape or not np.allclose(native_q,q,atol=1e-9,rtol=0):
        raise ValueError('Candidate charges differ from native BCC output')
    if abs(error) > bound+len(q)*0.5e-6+1e-8 or abs(q.sum()-intermediate[-1]) > 1e-5:
        raise ValueError('Charge discrepancy exceeds evidenced native print precision')
    return dict(report, basis='Verified rounded SQM atomic charges; native BCC conserves their sum',
                sqm_reported_total_e=float(total), printed_atomic_sum_e=printed_sum,
                maximum_print_rounding_error_e=bound, bcc_intermediate_totals_e=intermediate)


def prepare(baseline, output, amber):
    output.mkdir(parents=True, exist_ok=False)
    constraints = json.loads((baseline/'resp/constraints.json').read_text())
    mapping = json.loads((baseline/'capped-adduct.json').read_text())['atom_map']
    source = baseline/'interface-native-scope-v4/hybrid.mol2'
    data = mol2_sections(source)
    atoms, bonds = data['ATOM'], data['BOND']
    assert len(atoms) == len(mapping) == constraints['atom_count']
    assert [int(a[0]) for a in atoms] == list(range(1, len(atoms)+1))
    caps = set(constraints['cap_indices'])
    canonical = {a['index']: a for a in constraints['frozen_atoms']}
    names = {}
    for i, record in enumerate(mapping):
        # Antechamber infers elements from atom names during format checking;
        # arbitrary L/Z labels are invalid even when MOL2 atom types are valid.
        element = re.sub(r'\d+$', '', atoms[i][1])
        if i in caps:
            names[i] = f'{element}{i:03d}'
        elif i in canonical:
            names[i] = canonical[i]['atom']
        elif record.get('resname') == 'CYS':
            names[i] = record['atom']
        else:
            names[i] = f'{element}{i:03d}'
    assert len(set(names.values())) == len(names)
    q = np.array([float(a[8]) for a in atoms])
    charge_admission = validate_native_charge_total(baseline,q,constraints['formal_charge'])
    lines = ['@<TRIPOS>MOLECULE', 'COV', f'{len(atoms)} {len(bonds)} 1 0 0',
             'SMALL', 'USER_CHARGES', '', '@<TRIPOS>ATOM']
    for i, row in enumerate(atoms):
        row = row.copy(); row[1] = names[i]; row[7] = 'COV'
        lines.append(' '.join(row))
    lines += ['@<TRIPOS>BOND', *[' '.join(row) for row in bonds],
              '@<TRIPOS>SUBSTRUCTURE', '1 COV 1 RESIDUE 1 A COV 1']
    (output/'capped.mol2').write_text('\n'.join(lines)+'\n')
    commands = [native(amber, output, ['antechamber', '-i', 'capped.mol2', '-fi',
        'mol2', '-o', 'capped.ac', '-fo', 'ac', '-j', '0', '-s', '0', '-pf', 'n'])]
    # No charge-generation flag: this invocation only converts the retained
    # native joint AM1-BCC calculation to PREPGEN's input format.
    mainchain = ['HEAD_NAME N', 'TAIL_NAME C', 'MAIN_CHAIN CA',
                 'PRE_HEAD_TYPE C', 'POST_TAIL_TYPE N',
                 f"CHARGE {constraints['formal_charge']:.1f}"]
    mainchain += ['OMIT_NAME ' + names[i] for i in sorted(caps)]
    (output/'residue.mc').write_text('\n'.join(mainchain)+'\n')
    commands.append(native(amber, output, ['prepgen', '-i', 'capped.ac', '-o',
        'covalent.prepi', '-m', 'residue.mc', '-rn', 'COV']))
    (output/'gaff-adduct.frcmod').write_bytes(
        (baseline/'interface-native-scope-v4/gaff-adduct.frcmod').read_bytes())
    (output/'leap.in').write_text('source leaprc.protein.ff14SB\n'
        'source leaprc.gaff2\nloadamberparams gaff-adduct.frcmod\n'
        'loadamberprep covalent.prepi\nmodel=sequence {ACE ALA COV GLY NME}\n'
        'check model\nsaveamberparm model peptide.prmtop peptide.inpcrd\n'
        'savepdb model peptide.pdb\nquit\n')
    commands.append(native(amber, output, ['tleap', '-f', 'leap.in']))
    model = parmed.load_file(str(output/'peptide.prmtop'), xyz=str(output/'peptide.inpcrd'))
    modified = [r for r in model.residues if r.name == 'COV']
    if len(modified) != 1:
        raise ValueError('Modified residue was not inserted exactly once')
    residue = modified[0]; actual = {a.name: a for a in residue.atoms}
    retained = set(range(len(atoms)))-caps
    if set(actual) != {names[i] for i in retained}:
        raise ValueError('PREPGEN changed retained atom inventory')
    expected_bonds = {tuple(sorted((names[int(b[1])-1], names[int(b[2])-1])))
        for b in bonds if int(b[1])-1 in retained and int(b[2])-1 in retained}
    actual_bonds = {tuple(sorted((b.atom1.name, b.atom2.name))) for b in model.bonds
        if b.atom1.residue is residue and b.atom2.residue is residue}
    if actual_bonds != expected_bonds:
        raise ValueError('Native assembly changed adduct connectivity')
    if abs(sum(a.charge for a in model.atoms)-constraints['formal_charge']) > 1e-4:
        raise ValueError('Peptide charge is inconsistent with declared state')
    boundary = {}
    for group, count in [('bonds', 2), ('angles', 3), ('dihedrals', 4)]:
        found = []
        for term in getattr(model, group):
            members = [getattr(term, 'atom'+str(i+1)) for i in range(count)]
            if any(a.residue is residue for a in members) and any(a.residue is not residue for a in members):
                if term.type is None:
                    raise ValueError('Missing peptide-boundary parameter')
                found.append([a.idx for a in members])
        boundary[group] = found
    if len(boundary['bonds']) != 2:
        raise ValueError('Expected exactly two peptide-boundary bonds')
    topology = app.AmberPrmtopFile(str(output/'peptide.prmtop'))
    system = topology.createSystem(nonbondedMethod=app.NoCutoff, constraints=None)
    integrator = mm.VerletIntegrator(.001)
    context = mm.Context(system, integrator, mm.Platform.getPlatformByName('Reference'))
    context.setPositions(model.coordinates*unit.angstrom)
    state = context.getState(getEnergy=True, getForces=True)
    energy = float(state.getPotentialEnergy().value_in_unit(unit.kilojoule_per_mole))
    force = state.getForces(asNumpy=True).value_in_unit(unit.kilojoule_per_mole/unit.nanometer)
    if not np.isfinite(energy) or not np.isfinite(force).all():
        raise ValueError('Nonfinite peptide energy or force')
    del context, integrator
    changes = [{'capped_index': i, 'name': names[i], 'source': mapping[i],
        'original_am1bcc_charge': float(q[i]), 'native_prepgen_charge': actual[names[i]].charge,
        'change_e': actual[names[i]].charge-float(q[i]),
        'canonical_ff14sb_charge': canonical[i]['charge_e'] if i in canonical else None}
        for i in sorted(retained)]
    report = {'stage': 'native_modified_residue_assembly_complete',
        'protocol': 'joint AM1-BCC, native PREPGEN cap removal/charge redistribution; ff14SB backbone mechanics and scoped GAFF2 adduct mechanics',
        'is_constrained_resp': False, 'physical_model_validated': False,
        'full_complex_prepared': False, 'app_ready': False,
        'source_sha256': {str(p.relative_to(baseline)): digest(p) for p in
            [source, baseline/'resp/constraints.json', baseline/'capped-adduct.json',
             baseline/'interface-native-scope-v4/gaff-adduct.frcmod']},
        'native_tools': {n: digest(amber/'bin'/n) for n in ['antechamber', 'prepgen', 'tleap']},
        'commands': commands, 'atom_count': len(model.atoms), 'adduct_atom_count': len(retained),
        'input_charge_precision':charge_admission,
        'removed_cap_charge_e': float(q[list(caps)].sum()),
        'total_peptide_charge_e': float(sum(a.charge for a in model.atoms)),
        'adduct_charge_e': float(sum(a.charge for a in residue.atoms)),
        'charge_changes': changes,
        'maximum_prepgen_charge_change_e': max(abs(x['change_e']) for x in changes),
        'canonical_backbone_charge_max_difference_e': max(abs(x['native_prepgen_charge']-x['canonical_ff14sb_charge'])
            for x in changes if x['canonical_ff14sb_charge'] is not None),
        'boundary_term_counts': {k:len(v) for k,v in boundary.items()},
        'initial_energy_kj_mol': energy,
        'initial_max_force_kj_mol_nm': float(np.max(np.linalg.norm(force, axis=1))),
        'outputs': {p.name: digest(p) for p in output.iterdir() if p.is_file()}}
    (output/'result.json').write_text(json.dumps(report, indent=2)+'\n')
    return {k:v for k,v in report.items() if k not in ['charge_changes', 'outputs', 'commands']}


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    for n in ['baseline', 'output', 'amber']:
        parser.add_argument('--'+n, type=Path, required=True)
    args = parser.parse_args()
    print(json.dumps(prepare(args.baseline, args.output, args.amber), indent=2))
