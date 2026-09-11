"""Expected rejections and actual ion-preparation continuity in an isolated workspace."""
import json
import os
from pathlib import Path
import sys
import time
import traceback

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT))
from check_chemistry_matrix import ccd_dataset
import numpy as np
from openmm import app, unit
from backend import config, jobs, preparation, sources, storage
from backend.prepared_system import load_prepared_forcefield
from backend.ligands import inspect_ligands

OUT=ROOT/'docs/audit/release-readiness/chemistry-boundaries.json'


def main():
    if config.DATA_ROOT==ROOT/'data': raise SystemExit('Use an isolated DYNAMOL_DATA_DIR.')
    report={'scope':'Expected unsupported chemistry remains preserved and visibly rejected; real retained-ion/capped-peptide preparation and solvation checks. No coordination accuracy claim.','checks':[],'passed':False}
    def save(): storage.atomic_json(OUT,report)
    def reject(name, meta, pattern):
        rows=inspect_ligands(meta['id'])
        assert any(pattern.lower() in str(r.get('error')).lower() for r in rows), rows
        report['checks'].append({'id':name,'expected':'unsupported','passed':True,'messages':[r['error'] for r in rows],'dataset_id':meta['id']});save()
    try:
        reject('boronic acid outside validated organic element set',sources.create_smiles('CB(O)O'),'elements outside')
        reject('radical methyl',sources.create_smiles('[CH3]'),'closed-shell')
        reject('heme with iron',ccd_dataset('HEM'),'elements outside')
        fad=ccd_dataset('FAD');folder=storage.dataset_dir(fad['id'])
        link=list(' '*80);link[:6]='LINK  ';link[17:20]='FAD';link[47:50]='CYS'
        with (folder/'source.pdb').open('a') as f: f.write(''.join(link)+'\n')
        reject('covalently linked FAD',fad,'covalent')
        atp=ccd_dataset('ATP');folder=storage.dataset_dir(atp['id'])
        pdb=app.PDBFile(str(folder/'topology.pdb'));model=app.Modeller(pdb.topology,pdb.positions)
        missing=next(a for a in model.topology.atoms() if a.name=='O1A');model.delete([missing])
        with (folder/'topology.pdb').open('w') as f:app.PDBFile.writeFile(model.topology,model.positions,f,keepIds=True)
        reject('ATP missing heavy atom',atp,'heavy-atom identity mismatch')
        # A complete capped alanine carries no deposited sequence gaps.
        source=config.DATA_ROOT/'chemistry-matrix/protein-native/ALA.pdb'
        peptide=app.PDBFile(str(source));model=app.Modeller(peptide.topology,peptide.positions)
        ions=app.Topology();xyz=[]
        for i,(name,symbol) in enumerate([('NA','Na'),('CL','Cl'),('K','K'),('MG','Mg'),('CA','Ca')]):
            r=ions.addResidue(name,ions.addChain(str(i+1)),'1');ions.addAtom(name,app.element.Element.getBySymbol(symbol),r)
            xyz.append([.5*(i-2),1.5,0])
        model.add(ions,np.array(xyz)*unit.nanometer)
        path=config.DATA_ROOT/'capped-alanine-ions.pdb'
        with path.open('w') as f:app.PDBFile.writeFile(model.topology,model.positions,f,keepIds=True)
        meta=sources.import_structure(path,name='ACE–ALA–NME with Na, Cl, K, Mg, Ca',provenance={'audit_fixture':'Deliberately separated nonbonded ions; not a coordination site'})
        settings,inspection=preparation.validate_preparation({'dataset_id':meta['id'],'optimize_sidechains':False})
        assert len(inspection['ions'])==5 and not inspection['ligands']
        def wait(job):
            start=time.monotonic()
            while True:
                state=jobs.get_job(job['id'])
                if state['status'] in {'completed','failed','cancelled','interrupted'}:
                    if state['status']!='completed': raise RuntimeError(state.get('error',state['status']))
                    return storage.get_dataset(state['dataset_id']),state
                if time.monotonic()-start>180:
                    jobs.cancel_job(job['id']);raise RuntimeError('Native ion continuity exceeded three-minute limit.')
                time.sleep(.2)
        prepared,prep_job=wait(preparation.submit_preparation(settings))
        assert prepared['preparation']['requires_explicit_solvent'] is True
        assert len(prepared['preparation']['ions'])==5
        exact=app.PDBFile(str(storage.dataset_dir(prepared['id'])/'prepared.pdb'))
        ff,_=load_prepared_forcefield(storage.dataset_dir(prepared['id']),prepared['preparation'],solvent='explicit')
        import openmm as mm
        system=ff.createSystem(exact.topology,nonbondedMethod=app.NoCutoff)
        force=next(f for f in system.getForces() if isinstance(f,mm.NonbondedForce))
        actual={a.residue.name:force.getParticleParameters(a.index)[0].value_in_unit(unit.elementary_charge) for a in exact.topology.atoms() if a.residue.name in {'NA','CL','K','MG','CA'}}
        assert actual=={'NA':1,'CL':-1,'K':1,'MG':2,'CA':2}
        solvated,solvent_job=wait(preparation.submit_solvation(prepared['id'],{'padding_nm':1,'ph':7,'seed':2026}))
        assert solvated['preparation']['ions']==prepared['preparation']['ions']
        report['checks'].append({'id':'capped protein with five retained ions','passed':True,'preparation_job':prep_job['id'],'solvation_job':solvent_job['id'],'prepared_dataset':prepared['id'],'solvated_dataset':solvated['id'],'charge_states_e':actual,'solvated_atoms':solvated['n_atoms'],'requires_explicit_solvent':True,'all_ion_provenance_retained':True})
        report['passed']=True
    except Exception as exc:
        report['error']=str(exc);report['traceback']=traceback.format_exc()
    save();print(json.dumps(report,indent=2))

if __name__=='__main__': main()
