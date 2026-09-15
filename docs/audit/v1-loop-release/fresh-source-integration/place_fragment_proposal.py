"""Properly place a complete native proposal before exact anchor locking."""
import copy
import importlib.util
import json
from pathlib import Path
import sys
import numpy as np
from backend.loop_context import admit_loop_proposal, _scope, peptide_neighbors, outer_torsion_anchor_indices
from backend.loop_snapshot import load_snapshot, _read, _json
from backend.residue_identity import residue_key

spec=importlib.util.spec_from_file_location('pool_selector',Path(__file__).with_name('select_fragment_pool.py'))
selector=importlib.util.module_from_spec(spec);spec.loader.exec_module(selector)


def place(request):
    if set(request)!={'snapshot','snapshot_sha256','proposal','proposal_sha256'}:
        raise ValueError('Exact snapshot and native proposal hash bindings are required')
    snapshot=load_snapshot(request['snapshot'],expected_sha256=request['snapshot_sha256'])
    proposal=_json(_read(Path(request['proposal']),request['proposal_sha256']))
    sources={f.name:f for f in snapshot.source_files}
    top,_=selector.validator.topology_from_snapshot(snapshot,sources['prepared-topology.cif'].path)
    atoms=list(top.atoms());residues={residue_key(r):r for r in top.residues()}
    modeled=set(map(tuple,snapshot.modeled_keys));previous,following=peptide_neighbors(top)
    flanks={residue_key(n)for key in modeled for n in (previous.get(residues[key]),following.get(residues[key]))
            if n is not None and residue_key(n)not in modeled}
    native=admit_loop_proposal(proposal,top,snapshot.prepared_xyz_nm,modeled,flanks,
        source_sha256=sources['loop-model-context.pdb'].sha256,observed_context_policy='mobile_flanks')
    actual_flanks=set(native.proposed_flank_keys)
    if actual_flanks!=flanks:raise ValueError('Coherent placement requires complete immediate observed flanks')
    xyz=native.seeded_xyz_nm.copy();fixed=outer_torsion_anchor_indices(top,modeled,flanks)
    components=_scope(top,modeled,flanks,residues);records=[]
    for component in components:
        neighbors={key for key in flanks if any(n is not None and residue_key(n)in component
                   for n in (previous.get(residues[key]),following.get(residues[key])))}
        selected=component|neighbors
        indices=[a.index for a in atoms if residue_key(a.residue)in selected and a.element.symbol not in {'H','D'}]
        anchors=sorted(set(indices)&fixed)
        rotation,translation=selector.proper_fit(xyz[anchors],snapshot.prepared_xyz_nm[anchors])
        residual=np.linalg.norm(xyz[anchors]@rotation+translation-snapshot.prepared_xyz_nm[anchors],axis=1)*10
        xyz[indices]=xyz[indices]@rotation+translation
        records.append({'component':sorted(component),'flanks':sorted(neighbors),
            'rotation':rotation.tolist(),'translation_nm':translation.tolist(),'determinant':float(np.linalg.det(rotation)),
            'anchor_residual_angstrom':residual.tolist(),'anchor_rms_angstrom':float(np.sqrt(np.mean(residual**2)))})
    result=copy.deepcopy(proposal)
    for group,keys in [('modeled_heavy_atoms',modeled),('observed_context_atoms',flanks)]:
        result[group]=sorted([{'identity':[*residue_key(a.residue),a.name],'element':a.element.symbol,'xyz_nm':xyz[a.index].tolist()}
            for a in atoms if residue_key(a.residue)in keys and a.element.symbol not in {'H','D'}],key=lambda r:r['identity'])
    result['coherent_placement']={'method':'Independent proper rigid component fit to outer torsion anchors before exact locking.',
        'parent_native_proposal_sha256':request['proposal_sha256'],'snapshot_sha256':snapshot.sha256,
        'components':records,'complete_candidate_accepted':False}
    admit_loop_proposal(result,top,snapshot.prepared_xyz_nm,modeled,flanks,
                        source_sha256=sources['loop-model-context.pdb'].sha256)
    return result


if __name__=='__main__':
    result=place(_json(_read(Path(sys.argv[1]))))
    with Path(sys.argv[2]).open('x')as handle:
        json.dump(result,handle,indent=2,allow_nan=False);handle.write('\n')
