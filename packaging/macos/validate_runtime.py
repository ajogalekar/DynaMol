#!/usr/bin/env python3
"""Run native work using only the relocated bundled engines and Python."""
from __future__ import annotations
import argparse
import hashlib
import importlib.util
import json
import math
import os
from pathlib import Path
import subprocess
import sys
import tempfile
import time


def validate_promod3(prefix, expected_version, run, folder):
    """Exercise the loop runtime in its own interpreter, with native data lookup."""
    probe = '''import json, sys
from pathlib import Path
import ost, promod3, openmm
from ost import conop
prefix = Path(sys.prefix).resolve()
expected = Path(sys.argv[1]).resolve()
assert prefix == expected, f"Loop runtime escaped its private prefix: {prefix}"
assert promod3.__version__ == sys.argv[2], "ProMod3 version differs from manifest"
# This OpenStructure build links to the OpenMM 8.5.1 C++ ABI. Its OpenMM is
# private to loop modeling; the application's simulation engine is separate.
assert openmm.__version__ == "8.5.1", "Loop runtime OpenMM ABI pin changed"
paths = {"openstructure": Path(ost.GetSharedDataPath()).resolve(),
         "promod3": Path(promod3.GetProMod3SharedDataPath()).resolve()}
assert all(path.is_dir() and path.is_relative_to(prefix) for path in paths.values()), "Loop data escaped its private prefix"
assert conop.GetDefaultLib() is not None, "OpenStructure compound library was not loaded"
for module in (ost, promod3, openmm):
    assert Path(module.__file__).resolve().is_relative_to(prefix), "Loop module loaded outside its private prefix"
print(json.dumps({"version": promod3.__version__, "openstructure_version": ost.__version__,
                  "private_openmm_version": openmm.__version__, "prefix": str(prefix),
                  "data_paths": {name: str(path) for name, path in paths.items()},
                  "operation": "Private interpreter imports, pinned versions and bundled data lookup; no loop modeling"}))
'''
    output = run('promod3-library-check', [prefix / 'bin/python', '-I', '-B', '-c', probe, prefix, expected_version], folder)
    return json.loads(output.strip().splitlines()[-1])


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('resources', type=Path)
    parser.add_argument('--home', type=Path)
    parser.add_argument('--protein', type=Path)
    args = parser.parse_args()
    resources = args.resources.resolve()
    home = (args.home or Path(tempfile.mkdtemp(prefix='DynaMol Runtime Validation '))).resolve()
    bundle = resources.parent.parent
    if home.is_relative_to(bundle):
        raise SystemExit('Validation home must be outside the sealed application bundle.')
    home.mkdir(parents=True, exist_ok=True)
    spec = importlib.util.spec_from_file_location('dynamol_bootstrap', resources / 'bootstrap.py')
    bootstrap = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(bootstrap)
    manifest = json.loads((resources / 'manifest.json').read_text())
    with (home / 'unpack.log').open('a') as log:
        engines = bootstrap.ensure_engines(home, manifest, log)
    env = bootstrap.private_environment(home / 'Workspace', engines)
    env.update(OMP_NUM_THREADS='1', HTTP_PROXY='http://127.0.0.1:1', HTTPS_PROXY='http://127.0.0.1:1', ALL_PROXY='http://127.0.0.1:1', NO_PROXY='', DYLD_PRINT_LIBRARIES='1')
    operations = []
    def run(name, command, folder, input=None):
        started = time.monotonic()
        result = subprocess.run([str(arg) for arg in command], cwd=folder, env=env, input=input, text=True, capture_output=True, timeout=180)
        (folder / f'{name}.stdout.log').write_text(result.stdout)
        (folder / f'{name}.stderr.log').write_text(result.stderr)
        loaded = []
        for line in result.stderr.splitlines():
            if line.startswith('dyld[') and ' /' in line:
                path = line[line.index(' /') + 1:]
                loaded.append(path)
                if not path.startswith(('/usr/lib/', '/System/Library/', '/Library/Apple/usr/lib/', str(resources) + '/', str(home) + '/')):
                    raise AssertionError(f'{name} loaded an external non-system library: {path}')
        operations.append({'name': name, 'returncode': result.returncode, 'elapsed_seconds': round(time.monotonic()-started, 3), 'native_library_paths': sorted(set(loaded))})
        if result.returncode:
            raise RuntimeError(f'{name} failed. See {folder / (name + ".stderr.log")}')
        return result.stdout
    loop_runtime = None
    if 'promod3' in engines:
        loop_runtime = validate_promod3(engines['promod3'], manifest['engines']['promod3']['version'], run, home)
    assert str(Path(sys.prefix).resolve()).startswith(str(resources / 'python'))
    import openmm as mm
    from openmm import app, unit
    from rdkit import Chem
    from rdkit.Chem import AllChem
    import mdtraj, pdbfixer, gemmi, parmed, dimorphite_dl
    # Loading the actual installed XML resources verifies data lookup as well as
    # library import. CPU integration verifies the real native OpenMM platform.
    app.ForceField('amber14-all.xml', 'amber14/tip3p.xml')
    system = mm.System()
    system.addParticle(12)
    system.addParticle(12)
    force = mm.HarmonicBondForce()
    force.addBond(0, 1, 0.15, 1000)
    system.addForce(force)
    integrator = mm.VerletIntegrator(0.001)
    context = mm.Context(system, integrator, mm.Platform.getPlatformByName('CPU'), {'Threads': '1'})
    context.setPositions([[0,0,0],[0.18,0,0]])
    before = context.getState(getEnergy=True).getPotentialEnergy().value_in_unit(unit.kilojoules_per_mole)
    mm.LocalEnergyMinimizer.minimize(context, maxIterations=20)
    integrator.step(10)
    after = context.getState(getEnergy=True).getPotentialEnergy().value_in_unit(unit.kilojoules_per_mole)
    assert math.isfinite(after) and after < before
    del context, integrator
    amber_folder = home / 'Amber Ethanol Native Test'
    amber_folder.mkdir(exist_ok=True)
    mol = Chem.AddHs(Chem.MolFromSmiles('CCO'))
    assert AllChem.EmbedMolecule(mol, randomSeed=17) == 0
    AllChem.MMFFOptimizeMolecule(mol)
    writer = Chem.SDWriter(str(amber_folder / 'input.sdf')); writer.write(mol); writer.close()
    amber = engines['ambertools']
    run('antechamber', [amber/'bin/antechamber','-i','input.sdf','-fi','sdf','-o','charged.mol2','-fo','mol2','-c','bcc','-nc','0','-at','gaff2','-s','2','-pf','y'], amber_folder)
    run('parmchk2', [amber/'bin/parmchk2','-i','charged.mol2','-f','mol2','-o','ligand.frcmod','-s','gaff2'], amber_folder)
    (amber_folder/'leap.in').write_text('source leaprc.gaff2\nloadamberparams ligand.frcmod\nlig = loadmol2 charged.mol2\ncheck lig\nsaveamberparm lig ligand.prmtop ligand.inpcrd\nquit\n')
    run('tleap', [amber/'bin/tleap','-f','leap.in'], amber_folder)
    env['AMBERHOME'] = str(amber)
    run('sqm-library-check', [amber/'bin/sqm', '-O', '-i', 'sqm.in', '-o', 'sqm-library-check.out'], amber_folder)
    ligand = parmed.load_file(str(amber_folder/'ligand.prmtop'), str(amber_folder/'ligand.inpcrd'))
    assert len(ligand.atoms) == 9
    assert abs(sum(atom.charge for atom in ligand.atoms)) < 1e-5
    gmx_folder = home / 'GROMACS Native Test'
    gmx_folder.mkdir(exist_ok=True)
    gmx = engines['gromacs']/'bin/gmx'
    protein = (args.protein or resources/'app/examples/ubiquitin-start/topology.pdb').resolve()
    run('pdb2gmx', [gmx,'pdb2gmx','-f',protein,'-o','protein.gro','-p','topol.top','-i','posre.itp','-ff','amber99sb-ildn','-water','tip3p','-ignh'], gmx_folder)
    run('editconf', [gmx,'editconf','-f','protein.gro','-o','boxed.gro','-d','1.2','-bt','cubic'], gmx_folder)
    (gmx_folder/'minimize.mdp').write_text('integrator = steep\nnsteps = 5\nemtol = 1000\nemstep = 0.01\ncutoff-scheme = Verlet\nnstlist = 10\nrlist = 1.0\ncoulombtype = Cut-off\nrcoulomb = 1.0\nrvdw = 1.0\npbc = xyz\n')
    run('grompp', [gmx,'grompp','-f','minimize.mdp','-c','boxed.gro','-p','topol.top','-o','minimize.tpr'], gmx_folder)
    run('mdrun', [engines['gromacs']/'bin.ARM_NEON_ASIMD/gmx','mdrun','-deffnm','minimize','-ntmpi','1','-ntomp','1','-nb','cpu','-pme','cpu','-bonded','cpu','-update','cpu'], gmx_folder)
    assert (gmx_folder/'minimize.gro').stat().st_size > 1000
    run('python-library-check', [sys.executable, '-I', '-B', '-c', 'import openmm,rdkit,mdtraj; print(openmm.Platform.getNumPlatforms())'], home)
    for operation in operations:
        if operation['name'] in {'sqm-library-check', 'mdrun', 'python-library-check', 'promod3-library-check'}:
            assert operation['native_library_paths'], f"No loaded-library evidence for {operation['name']}"
    report = {'passed': True, 'python': sys.version, 'executable': sys.executable, 'prefix': sys.prefix, 'python_paths': sys.path, 'resources': str(resources), 'home': str(home), 'path': env['PATH'], 'network': 'No network operations used; all proxy endpoints disabled. Not an OS-level network sandbox.', 'openmm': {'version': mm.__version__, 'platform': 'CPU', 'steps': 10, 'energy_before_kj_mol': before, 'energy_after_kj_mol': after}, 'amber': {'atoms': len(ligand.atoms), 'charge': sum(atom.charge for atom in ligand.atoms), 'method': 'GAFF2/AM1-BCC with native SQM, Parmchk2 and LEaP'}, 'gromacs': {'system': 'bundled Ubiquitin', 'operation': 'pdb2gmx, editconf, grompp and 5-step vacuum minimization'}, 'operations': operations, 'limitation': 'Relocated paths with spaces and sanitized PATH on the development Mac; not a separate clean-machine or Gatekeeper test.'}
    if loop_runtime is not None:
        report['promod3'] = loop_runtime
    (home/'runtime-validation.json').write_text(json.dumps(report, indent=2)+'\n')
    print(json.dumps({'passed': True, 'report': str(home/'runtime-validation.json')},indent=2))


if __name__ == '__main__':
    main()
