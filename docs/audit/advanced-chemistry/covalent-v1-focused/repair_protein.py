"""Use DynaMol's existing protein/loop helpers for a focused research case.

The ligand/ion/water environment is kept separately for complete assembly.
This intermediate protein-only model must not be passed off as the complex.
"""
from __future__ import annotations

import argparse
import json
import os
from pathlib import Path
import sys
import time

import gemmi
import numpy as np
import openmm as mm
from openmm import app, unit
from pdbfixer import PDBFixer


def run(repo, source, output, extended_flanks=False):
    output.mkdir(parents=True,exist_ok=False)
    os.environ['DYNAMOL_DATA_DIR']=str(output.parent/'data')
    os.environ['DYNAMOL_PROMOD3']=str(Path.home()/'Library/Application Support/DynaMol/Runtimes/promod3-dfd2a559551cd2f2')
    os.environ['OPENMM_CPU_THREADS']='2'
    sys.path.insert(0,str(repo))
    from backend.preparation import find_missing_residues_preserving_identity, validate_loop_selection
    from backend.preparation_worker import (atom_key, proper_overlay_points,
        stereochemistry_report, require_valid_stereochemistry)
    from backend.sequence_evidence import (load_scheme, source_cif_chain_ids,
        restoration_plan, restore_residue_identities)
    from backend.residue_identity import STANDARD_PROTEINS, residue_key
    from backend.loop_modeling import generate_loop_model
    from backend.loop_refinement import minimize_modeled_loops
    from backend.loop_geometry import loop_geometry_report
    from backend.added_atom_stereo import correct_added_branches

    start=time.monotonic()
    def progress(stage, **kw):
        record={'stage':stage,'elapsed_seconds':time.monotonic()-start,**kw}
        (output/'progress.json').write_text(json.dumps(record,indent=2)+'\n')
        print(stage,flush=True)
    platform=mm.Platform.getPlatformByName('CPU');platform.setPropertyDefaultValue('Threads','2')
    progress('Loading deposited structure')
    fixer=PDBFixer(filename=str(source),platform=platform)
    original=app.Modeller(fixer.topology,fixer.positions)
    protected=set()
    for a,b in original.topology.bonds():
        if a.residue is not b.residue and {a.name,b.name}!={'C','N'}:
            protected.update((residue_key(a.residue),residue_key(b.residue)))
    label_map=source_cif_chain_ids(source,fixer.topology)
    protein_chains={r.chain.id for r in fixer.topology.residues() if r.name in STANDARD_PROTEINS}
    fixer.sequence_scheme=load_scheme(source,{}, {},protein_chains,fixer.sequences,label_map)
    modeller=app.Modeller(fixer.topology,fixer.positions)
    modeller.delete([r for r in modeller.topology.residues() if r.name not in STANDARD_PROTEINS])
    modeller.delete([a for a in modeller.topology.atoms() if a.element == app.element.hydrogen])
    environment=app.Modeller(original.topology,original.positions)
    environment.delete([r for r in environment.topology.residues() if r.name in STANDARD_PROTEINS])
    with (output/'retained-environment.pdb').open('w') as f:
        app.PDBFile.writeFile(environment.topology,environment.positions,f,keepIds=True)
    fixer.topology,fixer.positions=modeller.topology,modeller.positions
    find_missing_residues_preserving_identity(fixer)
    chains=list(fixer.topology.chains())
    terminal={str(k):v for k,v in fixer.missingResidues.items()
        if k[1] in [0,len(list(chains[k[0]].residues()))]}
    selected={k:v for k,v in fixer.missingResidues.items()
        if k[1] not in [0,len(list(chains[k[0]].residues()))]}
    validate_loop_selection([{'chain':chains[k[0]].id,'residues':v} for k,v in selected.items()],True)
    fixer.missingResidues=selected
    plan=restoration_plan(fixer,selected)
    observed={atom_key(a):np.array(fixer.positions[a.index].value_in_unit(unit.nanometer))
              for a in fixer.topology.atoms()}
    if len(observed)!=fixer.topology.getNumAtoms():
        raise ValueError('Ambiguous observed protein atom identity')
    fixer.findMissingAtoms()
    progress('Constructing missing atoms and internal loops',loops={str(k):v for k,v in selected.items()})
    import pdbfixer.pdbfixer as implementation
    saved=implementation._overlayPoints
    implementation._overlayPoints=proper_overlay_points
    try:
        fixer.addMissingAtoms(seed=2026)
    finally:
        implementation._overlayPoints=saved
    identities=restore_residue_identities(fixer.topology,plan)
    fixer.positions,stereo_repairs=correct_added_branches(fixer.topology,fixer.positions,fixer.templates,set(observed))
    modeled={residue_key(r) for r in fixer.topology.residues()
        if not any(atom_key(a) in observed for a in r.atoms())}
    if len(modeled)!=sum(map(len,selected.values())):
        raise ValueError('Loop identity count does not match declared missing sequence')
    loop_provenance=None
    if selected:
        progress('Rebuilding loops with DynaMol ProMod3 workflow')
        fixer.positions,loop_provenance=generate_loop_model(fixer.topology,fixer.positions,set(observed),
            output,2026,environment=environment,on_progress=lambda s:progress(s))
    require_valid_stereochemistry(stereochemistry_report(fixer.topology,fixer.positions,fixer.templates),'Constructed protein')
    progress('Assigning protein hydrogens')
    forcefield=app.ForceField('amber14/protein.ff14SB.xml','amber14/tip3p.xml')
    modeller=app.Modeller(fixer.topology,fixer.positions)
    variants=modeller.addHydrogens(forcefield,pH=7.4,platform=platform)
    refinement=None
    if modeled:
        progress('Refining modeled loops')
        reference_positions=modeller.positions
        with (output/'before-refinement.pdb').open('w') as handle:
            app.PDBFile.writeFile(modeller.topology,modeller.positions,handle,keepIds=True)
        modeller.positions,refinement=minimize_modeled_loops(modeller.topology,modeller.positions,forcefield,modeled,
                                                            protected_residues=protected)
        mobile_flanks={tuple(k) for k in refinement.get('flank_relaxation',{}).get('residues',[])}
        geometry=loop_geometry_report(modeller.topology,modeller.positions,modeled|mobile_flanks)
        if extended_flanks and not geometry['accepted'] and mobile_flanks:
            # Research fallback for fragments selected with two residues of
            # observed context: permit those extra local neighbors to relax
            # under the same bounded displacement/chirality restraints.
            from backend.loop_refinement import _minimize_attempt
            extra=set()
            for a,b in modeller.topology.bonds():
                if a.residue is b.residue or {a.name,b.name}!={'C','N'}:continue
                ka,kb=residue_key(a.residue),residue_key(b.residue)
                if ka in mobile_flanks:extra.add(kb)
                if kb in mobile_flanks:extra.add(ka)
            available={residue_key(r) for r in modeller.topology.residues() if r.name in STANDARD_PROTEINS}
            extra=(extra&available)-modeled-protected
            flanks=mobile_flanks|extra
            prior=refinement
            progress('Testing bounded relaxation of two neighboring residues')
            targets={tuple(row['atoms']):np.deg2rad(row['target_degrees'])
                     for row in prior['peptide_construction_restraints']['targets']}
            # Restart from the valid pre-relaxation candidate. A rejected
            # fixed-anchor minimum can contain collinear angles, which are a
            # poor starting point even after freeing the neighboring atoms.
            modeller.positions,refinement=_minimize_attempt(modeller.topology,reference_positions,forcefield,modeled,
                mobile_flank_keys=flanks,reference_positions=reference_positions,target_overrides=targets)
            refinement['prior_attempts']=prior
            refinement['research_extended_context_fallback']=True
            refinement['fallback_start']='Original pre-refinement candidate, not rejected fixed-anchor minimum'
            refinement['total_iteration_budget']=prior['total_iteration_budget']+1000
            mobile_flanks=flanks
            geometry=loop_geometry_report(modeller.topology,modeller.positions,modeled|mobile_flanks)
        if mobile_flanks and refinement['flank_relaxation']['maximum_displacement_nm']>.1:
            raise ValueError('Loop flank refinement exceeded its displacement bound')
        (output/'loop-refinement.json').write_text(json.dumps(refinement,indent=2)+'\n')
        with (output/'refined-candidate.pdb').open('w') as handle:
            app.PDBFile.writeFile(modeller.topology,modeller.positions,handle,keepIds=True)
        if not geometry['accepted']:
            (output/'failed-loop-geometry.json').write_text(json.dumps(geometry,indent=2)+'\n')
            raise ValueError('Modeled loop geometry did not pass existing DynaMol checks')
    else:
        geometry=None
    require_valid_stereochemistry(stereochemistry_report(modeller.topology,modeller.positions,fixer.templates),'Refined protein')
    for residue,variant in zip(modeller.topology.residues(),variants):
        if residue.name=='HIS' and variant in ['HID','HIE','HIP']:
            residue.name=variant
    with (output/'protein.pdb').open('w') as handle:
        app.PDBFile.writeFile(modeller.topology,modeller.positions,handle,keepIds=True)
    report={'stage':'protein_intermediate_complete','full_complex_prepared':False,'app_ready':False,
        'source':str(source),'atoms':modeller.topology.getNumAtoms(),
        'modeled_residues':[list(k) for k in sorted(modeled)],
        'unobserved_terminal_segments_not_built':terminal,
        'environment_preserved_separately':True,'retained_environment_atoms':environment.topology.getNumAtoms(),
        'residue_identity_restoration':identities,'stereo_repairs':stereo_repairs,
        'loop_provenance':loop_provenance,'loop_refinement':refinement,'loop_geometry':geometry,
        'elapsed_seconds':time.monotonic()-start}
    (output/'result.json').write_text(json.dumps(report,indent=2)+'\n')
    progress('Protein intermediate complete; covalent/cofactor assembly remains')


if __name__=='__main__':
    p=argparse.ArgumentParser()
    for n in ['repo','source','output']:p.add_argument('--'+n,type=Path,required=True)
    p.add_argument('--extended-flanks',action='store_true')
    a=p.parse_args()
    try:
        run(a.repo,a.source,a.output,a.extended_flanks)
    except Exception as exc:
        if a.output.is_dir():
            (a.output/'failure.json').write_text(json.dumps({'error_type':type(exc).__name__,'error':str(exc)},indent=2)+'\n')
        raise
