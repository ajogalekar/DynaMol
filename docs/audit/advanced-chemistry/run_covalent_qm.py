"""Bounded, durable isolated capped-adduct HF optimization → ESP → RESP.

Each stage is gated on the prior stage's actual acceptance. This does not
publish parameters or mark the full protein preparation as supported.
"""
from __future__ import annotations
import hashlib
import importlib.util
import json
import os
import signal
import subprocess
import sys
import time
from pathlib import Path
import numpy as np
from rdkit import Chem

ROOT=Path(__file__).resolve().parents[3]


def load_module(path,name):
    spec=importlib.util.spec_from_file_location(name,path);module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module);return module


def main(folder,workspace):
    folder=folder.resolve();workspace=workspace.resolve();workspace.mkdir(parents=True,exist_ok=False)
    started=time.time();state={'status':'running','stage':'initializing','pid':os.getpid(),'started_unix':started,
        'scope':'Quantum capped-adduct parameter model; not a full prepared protein or a validated MD force field.'}
    def update(**kwargs):
        state.update(kwargs);state['elapsed_seconds']=time.time()-started
        temporary=workspace/'progress.tmp';temporary.write_text(json.dumps(state,indent=2)+'\n');temporary.replace(workspace/'progress.json')
    def native_qm(stage,request,limit):
        update(stage=stage)
        path=workspace/(stage+'-input.json');path.write_text(json.dumps(request,indent=2)+'\n')
        output=workspace/stage
        with (workspace/(stage+'.log')).open('w') as log:
            child=subprocess.Popen([str(ROOT/'.tools/qm/bin/python'),str(ROOT/'backend/qm_worker.py'),str(path),str(output)],
                cwd=ROOT,stdout=log,stderr=subprocess.STDOUT,start_new_session=True)
            update(child_pid=child.pid,stage_deadline_unix=time.time()+limit)
            try:code=child.wait(timeout=limit)
            except subprocess.TimeoutExpired:
                os.killpg(child.pid,signal.SIGTERM)
                try:child.wait(timeout=5)
                except subprocess.TimeoutExpired:os.killpg(child.pid,signal.SIGKILL);child.wait()
                raise TimeoutError(f'{stage} exceeded its independent subprocess deadline')
        result=json.loads((output/'result.json').read_text())
        if code or not result['accepted']:raise RuntimeError(f'{stage} was not accepted: {result.get("error")}')
        if result['atom_ids']!=request['atom_ids'] or result['elements']!=request['elements']:
            raise ValueError(f'{stage} changed atom identities/order')
        return output,result,np.load(output/'arrays.npz')
    try:
        request=json.loads((folder/'qm-hf-optimization-input.json').read_text())
        source=json.loads((folder/'capped-adduct.json').read_text())
        original=Chem.SDMolSupplier(str(folder/'capped-adduct.sdf'),removeHs=False)[0]
        optdir,opt,arrays=native_qm('optimization',request,request['max_wall_seconds']+120)
        coords=arrays['coords_angstrom'];optimized=Chem.Mol(original)
        for i,row in enumerate(coords):optimized.GetConformer().SetAtomPosition(i,row)
        Chem.AssignStereochemistryFrom3D(optimized,replaceExistingTags=True)
        Chem.AssignStereochemistry(optimized,cleanIt=True,force=True)
        names={a['name']:a['index'] for a in source['atom_map']}
        for name,expected in source['stereochemistry'].items():
            atom=optimized.GetAtomWithIdx(names[name])
            if not atom.HasProp('_CIPCode') or atom.GetProp('_CIPCode')!=expected:
                raise ValueError(f'Quantum optimization changed recorded stereochemistry at {name}')
        table=Chem.GetPeriodicTable();geometry=[]
        for bond in optimized.GetBonds():
            a,b=bond.GetBeginAtom(),bond.GetEndAtom();distance=float(np.linalg.norm(coords[a.GetIdx()]-coords[b.GetIdx()]))
            scale=table.GetRcovalent(a.GetAtomicNum())+table.GetRcovalent(b.GetAtomicNum())
            geometry.append({'atoms':[a.GetIdx(),b.GetIdx()],'distance_angstrom':distance,'radii_ratio':distance/scale})
            if not .7<distance/scale<1.35:raise ValueError('Quantum model has an implausible bonded distance; no ESP fit started')
        (workspace/'optimized-adduct.sdf').write_text(Chem.MolToMolBlock(optimized)+'\n$$$$\n')
        (workspace/'geometry.json').write_text(json.dumps({'stereo_preserved':True,'bond_checks':geometry,
            'scope':'Gross graph-based distance/stereo checks; not an experimental or force-field validation.'},indent=2)+'\n')
        surface=load_module(ROOT/'docs/audit/advanced-chemistry/vendor/resp-surface/vdw_surface.py','dynamol_resp_surface')
        shells=[];elements=[x.upper() for x in request['elements']]
        for scale in (1.4,1.6,1.8,2.0):
            points,radii=surface.vdw_surface(coords,elements,scale,1.0,{})
            shells.append(points)
        bohr=opt['units']['bohr_to_angstrom'];points=np.concatenate(shells)/bohr
        esp_request={k:v for k,v in request.items() if k not in ('coords_angstrom','optimization')}
        esp_request.update(coords_bohr=arrays['coords_bohr'].tolist(),operations=['esp'],esp_points_bohr=points.tolist(),
                           esp_batch_size=64,max_wall_seconds=7200)
        espdir,esp,esp_arrays=native_qm('esp',esp_request,7320)
        update(stage='resp',child_pid=None)
        values=esp_arrays['esp_hartree_per_e'];nuclei=esp_arrays['coords_bohr'];grid=esp_arrays['esp_points_bohr']
        lines=[f'{len(nuclei):5d}{len(grid):5d}{0:5d}']
        lines += [' '*17+''.join(f'{x:16.7E}' for x in row) for row in nuclei]
        lines += [' '+''.join(f'{x:16.7E}' for x in [potential,*point]) for potential,point in zip(values,grid)]
        espfile=workspace/'hf-resp.esp';espfile.write_text('\n'.join(lines)+'\n')
        # Keep the actual QM fit separate from the earlier synthetic fixture.
        fitfolder=workspace/'fitting';fitfolder.mkdir()
        import shutil
        for name in ('capped-adduct.sdf','capped-adduct.json'):shutil.copyfile(folder/name,fitfolder/name)
        fitter=load_module(ROOT/'docs/audit/advanced-chemistry/covalent_resp_spike.py','dynamol_covalent_resp')
        fit=fitter.run(fitfolder,espfile)
        fitted=np.array(fit['charges']);predicted=(fitted/np.linalg.norm(grid[:,None,:]-nuclei[None,:,:],axis=2)).sum(axis=1)
        residual=predicted-values;rmse=float(np.sqrt(np.mean(residual**2)));baseline=float(np.sqrt(np.mean(values**2)))
        report={'scope':state['scope'],'full_preparation_ready':False,'status':'completed','numerical_fit_constraints_passed':fit['constraint_checks_passed'],
            'qm_optimization':{'path':str(optdir),'elapsed_seconds':opt['elapsed_seconds'],'energy_hartree':opt['energy_hartree'],'nao':opt['nao']},
            'qm_esp':{'path':str(espdir),'elapsed_seconds':esp['elapsed_seconds'],'point_count':len(grid),'surface_shells':[1.4,1.6,1.8,2.0],
                'surface_density_per_angstrom2':1.0,'surface_source':'BSD3 cdsgroup/resp vdw_surface.py, source hashes preserved in audit/vendor/resp-surface/source.json'},
            'fit':{'path':str(fitfolder/'resp'),'training_rmse_hartree_per_e':rmse,'training_relative_rmse':rmse/baseline,
                'maximum_fixed_charge_error_e':fit['maximum_fixed_charge_error_e'],'charge_sum_e':fit['charge_sum_e'],
                'accuracy_acceptance':'Not assessed; no held-out conformer, QM torsion benchmark or full protein-interface validation yet.'},
            'input_sha256':hashlib.sha256((folder/'qm-hf-optimization-input.json').read_bytes()).hexdigest(),
            'elapsed_seconds':time.time()-started}
        (workspace/'result.json').write_text(json.dumps(report,indent=2)+'\n');update(status='completed',stage='completed',child_pid=None)
    except Exception as exc:
        update(status='failed',error={'type':type(exc).__name__,'message':str(exc)},child_pid=None)
        (workspace/'result.json').write_text(json.dumps(state,indent=2)+'\n');raise


if __name__=='__main__':
    if len(sys.argv)!=3:raise SystemExit('Usage: run_covalent_qm.py CAP_FOLDER FRESH_OUTPUT_DIRECTORY')
    main(Path(sys.argv[1]),Path(sys.argv[2]))
