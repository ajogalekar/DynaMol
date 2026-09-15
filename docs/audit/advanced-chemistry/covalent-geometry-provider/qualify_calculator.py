"""Small, bounded GFN2 numerical qualification; never optimizes a biomolecule."""
from pathlib import Path
import ast
import copy
import hashlib
import json
import os
import sys
import time
import traceback

for key in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','VECLIB_MAXIMUM_THREADS','NUMEXPR_NUM_THREADS']:os.environ[key]='1'
import numpy as np
import calculator as provider

HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[3]
OUT=ROOT/'build/advanced-chemistry/covalent/geometry-provider-qualification-v2'
BOHR_ANGSTROM=.529177210903


def save(path,data):path.write_text(json.dumps(data,indent=2)+'\n')
def request(elements,xyz,accuracy=.001):
    return {'schema_version':1,'method':'GFN2-xTB','atom_ids':[f'{i:03}:{e}' for i,e in enumerate(elements)],
        'elements':elements,'charge':0,'spin':0,'coords_bohr':np.asarray(xyz).tolist(),'accuracy':accuracy,'max_iter':250,'temperature_hartree':.00095}


def main():
    OUT.mkdir(parents=True,exist_ok=False);started=time.time();calls=0
    result={'status':'running','started_unix':started,'pid':os.getpid(),'tests':{},'physical_acceptance':False,
        'source_sha256':{p.name:provider.sha(p) for p in [HERE/'calculator.py',Path(__file__),HERE/'runtime-manifest.json',HERE/'sources/source-manifest.json']},
        'predeclared_limits':{'upstream_energy_hartree':1e-6,'upstream_gradient_hartree_bohr':1e-6,
            'finite_difference_max_error_hartree_bohr':2e-6,'finite_difference_steps_bohr':[.0002,.0001],
            'invariance_energy_hartree':1e-8,'invariance_gradient_hartree_bohr':1e-7,'outer_total_seconds':600,'threads':1}}
    def evaluate(r):
        nonlocal calls
        calls+=1;value=provider.evaluate(r);return value
    try:
        with provider.execution_guard(ROOT) as resources:
            result['resources']=resources;save(OUT/'progress.json',result)
            src=HERE/'sources/upstream-v0.7.0-test_interface.py';tree=ast.parse(src.read_text())
            func=next(n for n in tree.body if isinstance(n,ast.FunctionDef) and n.name=='get_ala')
            ns={'np':np};exec(compile(ast.Module(body=[func],type_ignores=[]),str(src),'exec'),ns)
            test=next(n for n in tree.body if isinstance(n,ast.FunctionDef) and n.name=='test_gfn2')
            gradnode=next(n.value for n in test.body if isinstance(n,ast.Assign) and any(isinstance(t,ast.Name) and t.id=='gradient' for t in n.targets))
            expected_grad=eval(compile(ast.Expression(gradnode),str(src),'eval'),{'np':np})
            elements_by_number={v:k for k,v in provider.ELEMENTS.items()};up=[]
            from tblite.interface import Calculator,Result
            numbers,xyz=ns['get_ala']('xac')
            native_calc=Calculator('GFN2-xTB',numbers,xyz,charge=0,uhf=0)
            native_calc.set('accuracy',1.0);native_calc.set('temperature',.00095);native_calc.set('verbosity',0)
            restart=Result()
            for conf,expected in [('xac',-32.97134042392298),('xag',-32.97132127543436)]:
                numbers,xyz=ns['get_ala'](conf)
                if conf=='xac':restart=native_calc.singlepoint(restart,copy=True)
                else:native_calc.update(xyz);native_calc.singlepoint(restart)
                calls+=1
                native_energy=float(restart.get('energy'));native_gradient=np.asarray(restart.get('gradient'))
                r=request([elements_by_number[int(n)] for n in numbers],xyz,1.0)
                save(OUT/f'upstream-{conf}-input.json',r)
                save(OUT/f'upstream-{conf}-result.json',{'energy_hartree':native_energy,'gradient_hartree_per_bohr':native_gradient.tolist(),
                    'protocol':'Exact upstream xac then updated xag with reused Result; SCC accuracy1.0, explicit charge0/spin0/temperature0.00095Ha.'})
                de=abs(native_energy-expected);assert de<1e-6
                row={'conformer':conf,'atoms':len(numbers),'expected_energy_hartree':expected,'actual_energy_hartree':native_energy,'absolute_error_hartree':de}
                if conf=='xag':
                    dg=float(np.max(np.abs(native_gradient-expected_grad)));assert dg<1e-6
                    row['all_66_gradient_components_max_error_hartree_bohr']=dg
                tight=request(r['elements'],xyz,.001);v=evaluate(tight);save(OUT/f'upstream-{conf}-tight-adapter-result.json',v)
                row['fresh_tight_adapter_energy_difference_hartree']=v['energy_hartree']-native_energy
                row['fresh_tight_adapter_gradient_difference_hartree_bohr']=float(np.max(np.abs(np.asarray(v['gradient_hartree_per_bohr'])-native_gradient)))
                up.append(row)
            result['tests']['upstream_gfn2_regression']={'source_sha256':provider.sha(src),'fixtures':up,'passed':True,
                'scope':'Exact upstream restart protocol. Difference to fresh/tighter SCC is separately reported, not hidden with a wider tolerance.',
                'preserved_first_attempt':'geometry-provider-qualification-v1: fresh default-accuracy second conformer failed stored restart-gradient tolerance.'}
            fixtures={
                'water':(['O','H','H'],[[0,0,0],[.96,0,0],[-.24,.93,0]]),
                'ammonia':(['N','H','H','H'],[[0,0,.11],[.94,0,-.22],[-.47,.82,-.23],[-.47,-.80,-.25]]),
                'hydrogen-sulfide':(['S','H','H'],[[0,0,0],[1.34,0,0],[-.06,1.34,0]]),
                'fluoromethane':(['C','F','H','H','H'],[[0,0,0],[1.38,0,0],[-.36,1.03,0],[-.36,-.51,.90],[-.36,-.51,-.89]]),
                'methanethiol':(['C','S','H','H','H','H'],[[0,0,0],[1.82,0,0],[-.36,1.03,0],[-.36,-.51,.90],[-.36,-.51,-.89],[2.23,1.27,.07]])}
            checks=[]
            for name,(elements,angstrom) in fixtures.items():
                r=request(elements,np.asarray(angstrom)/BOHR_ANGSTROM);v=evaluate(r);xyz=np.asarray(r['coords_bohr']);gradient=np.asarray(v['gradient_hartree_per_bohr'])
                save(OUT/f'{name}-input.json',r);save(OUT/f'{name}-result.json',v);errors=[];fd_values=[]
                for step in [.0002,.0001]:
                    fd=np.empty_like(xyz)
                    for index in np.ndindex(xyz.shape):
                        plus=copy.deepcopy(r);minus=copy.deepcopy(r)
                        xp=xyz.copy();xm=xyz.copy();xp[index]+=step;xm[index]-=step
                        plus['coords_bohr']=xp.tolist();minus['coords_bohr']=xm.tolist()
                        fd[index]=(evaluate(plus)['energy_hartree']-evaluate(minus)['energy_hartree'])/(2*step)
                    error=float(np.max(np.abs(fd-gradient)));errors.append(error);fd_values.append(fd)
                    if error>=2e-6:raise ValueError(f'{name} finite-difference failure at step {step}: {error}')
                arrays=OUT/f'{name}-derivatives.npz';np.savez(arrays,analytic=gradient,finite_difference=np.asarray(fd_values),coords_bohr=xyz)
                perm=np.arange(len(elements))[::-1];pr=copy.deepcopy(r)
                for field in ['atom_ids','elements','coords_bohr']:pr[field]=[r[field][i] for i in perm]
                pv=evaluate(pr);pe=abs(pv['energy_hartree']-v['energy_hartree']);pg=float(np.max(np.abs(np.asarray(pv['gradient_hartree_per_bohr'])-gradient[perm])))
                assert pe<1e-8 and pg<1e-7 and pv['atom_ids']==pr['atom_ids']
                q,_=np.linalg.qr(np.random.default_rng(310).normal(size=(3,3)))
                if np.linalg.det(q)<0:q[:,0]*=-1
                rr=copy.deepcopy(r);rr['coords_bohr']=(xyz@q+np.array([.25,-.60,.17])).tolist();rv=evaluate(rr)
                re=abs(rv['energy_hartree']-v['energy_hartree']);rg=float(np.max(np.abs(np.asarray(rv['gradient_hartree_per_bohr'])-gradient@q)))
                assert re<1e-8 and rg<1e-7
                checks.append({'fixture':name,'atoms':len(elements),'energy_hartree':v['energy_hartree'],'gradient_components_checked':gradient.size,
                    'finite_difference_max_errors_hartree_bohr':errors,'permutation_energy_error_hartree':pe,'permutation_gradient_error_hartree_bohr':pg,
                    'rigid_transform_energy_error_hartree':re,'rigid_transform_gradient_error_hartree_bohr':rg,'arrays_sha256':provider.sha(arrays)})
                result['tests']['small_molecules']=checks;result['stage']=name;save(OUT/'progress.json',result)
            bad=[];base=request(['O','H','H'],np.asarray(fixtures['water'][1])/BOHR_ANGSTROM)
            mutations={'wrong_method':{'method':'GFN1-xTB'},'unknown_setting':{'solvent':'water'},'implicit_spin':{'spin':None},'open_shell':{'spin':2},
                'fractional_charge':{'charge':.5},'odd_electron_count':{'charge':1},'duplicate_ids':{'atom_ids':['a','a','b']},
                'unsupported_element':{'elements':['Zn','H','H']},'missing_units':{'coords_bohr':None},'nonfinite':{'coords_bohr':[[float('nan'),0,0],[2,0,0],[0,2,0]]},
                'coincident':{'coords_bohr':[[0,0,0],[0,0,0],[0,2,0]]},'too_many_atoms':{'atom_ids':[str(i) for i in range(33)],'elements':['H']*33,'coords_bohr':[[i,0,0] for i in range(33)]},
                'bad_accuracy':{'accuracy':0},'bad_temperature_units':{'temperature_hartree':300},'bad_iterations':{'max_iter':0}}
            for name,mutation in mutations.items():
                r=copy.deepcopy(base);r.update(mutation)
                try:provider.validate(r)
                except (ValueError,TypeError) as error:bad.append({'case':name,'error':str(error)})
                else:raise AssertionError('Invalid input accepted: '+name)
            low=copy.deepcopy(base);low['max_iter']=1
            try:evaluate(low)
            except Exception as error:
                if type(error).__name__!='TBLiteRuntimeError':raise
                bad.append({'case':'native_nonconvergence','error_type':type(error).__name__,'error':str(error)})
            else:raise AssertionError('Native nonconvergence must not be accepted')
            result['tests']['rejections']={'cases':bad,'count':len(bad),'passed':True}
            result['status']='passed_small_molecule_numerical_qualification';result['native_evaluations_attempted']=calls
            result['elapsed_seconds']=time.time()-started;result['stage']='complete'
        try:provider.evaluate(base)
        except RuntimeError:result['tests']['unguarded_evaluation_rejected']=True
        else:raise AssertionError('Missing shared guard was accepted')
    except Exception as error:
        result.update(status='failed',error_type=type(error).__name__,error=str(error),elapsed_seconds=time.time()-started,native_evaluations_attempted=calls)
        (OUT/'traceback.log').write_text(traceback.format_exc())
    result['artifacts']={p.name:provider.sha(p) for p in sorted(OUT.iterdir()) if p.is_file() and p.name not in ['result.json','progress.json']}
    save(OUT/'result.json',result);save(OUT/'progress.json',result)
    print(json.dumps(result,indent=2))
    if result['status']=='failed':raise SystemExit(1)

if __name__=='__main__':main()
