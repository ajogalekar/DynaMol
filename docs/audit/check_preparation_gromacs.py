"""Bounded upstream compatibility probes, no MD; synthetic peptide excerpts only."""
from pathlib import Path
import json, os, random, subprocess
import numpy as np
from openmm import unit
from openmm.app import PDBFile, Modeller, ForceField
from pdbfixer import PDBFixer
ROOT=Path(__file__).resolve().parents[2]
OUT=ROOT/'docs/audit/preparation-gromacs-probes'
OUT.mkdir(exist_ok=True)
source=(ROOT/'data/sources/1UBQ.pdb').read_text().splitlines()
ff=ForceField('amber14/protein.ff14SB.xml','implicit/gbn2.xml')

def excerpt(name,start,end):
    atoms=[s for s in source if s.startswith('ATOM  ') and start<=int(s[22:26])<=end]
    res=[]
    for s in atoms:
        key=s[22:26]
        if not res or res[-1][0]!=key:res.append((key,s[17:20]))
    seq='SEQRES   1 A %4d  %s' % (len(res),' '.join(x[1] for x in res))
    p=OUT/f'{name}-heavy.pdb';p.write_text(seq+'\n'+'\n'.join(atoms)+'\nTER\nEND\n')
    f=PDBFixer(filename=str(p));f.findMissingResidues();f.findMissingAtoms();f.addMissingAtoms(seed=2026)
    return f

def atoms(p):
    # Match residue *ordinal*, name (residue labels change in mapped inputs).
    d={}; order=[]
    last=None; ordinal=-1
    for s in p.read_text().splitlines():
        if not s.startswith(('ATOM  ','HETATM')):continue
        r=(s[21],s[22:27]);a=s[12:16].strip()
        if r!=last:ordinal+=1;last=r
        key=f'{ordinal}:{a}'
        d[key]=np.array([float(s[30:38]),float(s[38:46]),float(s[46:54])])
        order.append(key)
    return d,order

# Upstream template bond lists let this probe expose the additional H-renaming work.
# Name assignment among H atoms on the same center is only a compatibility experiment.
rtps={}; active=None; section=None
for line in (ROOT/'.gromacs/share/gromacs/top/amber99sb-ildn.ff/aminoacids.rtp').read_text().splitlines():
    clean=line.split(';')[0].strip()
    if not clean:continue
    if clean.startswith('['):
        label=clean.strip('[] ').strip()
        if label in ('atoms','bonds','impropers','dihedrals','exclusions','cmap'):section=label
        else:active=label;section=None;rtps[active]=[]
    elif active and section=='bonds':rtps[active].append(clean.split()[:2])

def rename_hydrogens(m):
    changes=[]
    residues=list(m.topology.residues())
    bonds=list(m.topology.bonds())
    for r in residues:
        if r.name=='ILE':
            for a in r.atoms():
                if a.name=='CD1':
                    changes.append({'residue':r.id,'old':'CD1','new':'CD','element':'C'})
                    a.name='CD'
    for i,r in enumerate(residues):
        template=('N' if i==0 else 'C' if i==len(residues)-1 else '')+r.name
        desired={}
        for a,b in rtps[template]:
            if a.startswith('H'):desired.setdefault(b,[]).append(a)
            if b.startswith('H'):desired.setdefault(a,[]).append(b)
        existing={}
        for a,b in bonds:
            if a.residue==r and b.residue==r:
                if a.element.symbol=='H':existing.setdefault(b.name,[]).append(a)
                if b.element.symbol=='H':existing.setdefault(a.name,[]).append(b)
        assert set(desired)==set(existing),(template,desired,existing)
        for center,aa in existing.items():
            expected=desired[center]
            assert len(expected)==len(aa),(template,center,expected,[a.name for a in aa])
            # Preserve graph+coordinates; retain compatible names, then match others.
            common=set(expected)&{a.name for a in aa}
            remaining=sorted(set(expected)-common)
            for a,new in zip(sorted([a for a in aa if a.name not in common],key=lambda a:a.name),remaining):
                changes.append({'residue':r.id,'center':center,'old':a.name,'new':new})
                a.name=new
    return changes

results=[]
for name,start,end,ph in [('histidine-acid',66,72,2.0),('histidine-neutral',66,72,7.0),('acidic-acid',15,22,2.0),('acidic-neutral',15,22,7.0),('lysine-basic',24,30,12.0)]:
    f=excerpt(name,start,end)
    random.seed(2026);np.random.seed(2026)
    m=Modeller(f.topology,f.positions)
    try:
        variants=m.addHydrogens(ff,pH=ph)
    except Exception as exc:
        results.append({'fixture':name,'ph':ph,'stage':'addHydrogens','error':str(exc)})
        continue
    details=[{'id':r.id,'name':r.name,'variant':v,'hydrogens':[a.name for a in r.atoms() if a.element.symbol=='H']} for r,v in zip(m.topology.residues(),variants)]
    ff.createSystem(m.topology)
    for mode in ('plain','variants','variants_hmap'):
        d=OUT/f'{name}-{mode}';d.mkdir(exist_ok=True)
        if mode.startswith('variants'):
            for r,v in zip(m.topology.residues(),variants):
                if v:r.name=v
                names={a.name for a in r.atoms()}
                if r.name=='GLU' and 'HE2' in names:r.name='GLH'
                if r.name=='ASP' and 'HD2' in names:r.name='ASH'
                if r.name=='LYS' and 'HZ3' not in names:r.name='LYN'
        try:
            hmap=rename_hydrogens(m) if mode=='variants_hmap' else []
        except Exception as exc:
            results.append({'fixture':name,'ph':ph,'mode':mode,'stage':'H-name mapping','error':str(exc)})
            continue
        inp=d/'input.pdb'
        with inp.open('w') as fh:PDBFile.writeFile(m.topology,m.positions,fh,keepIds=True)
        cmd=[str(ROOT/'.gromacs/bin/gmx'),'pdb2gmx','-f','input.pdb','-o','processed.pdb','-p','topol.top','-i','posre.itp','-ff','amber99sb-ildn','-water','tip3p','-n','mapping.ndx']
        run=subprocess.run(cmd,cwd=d,stdout=subprocess.PIPE,stderr=subprocess.STDOUT,text=True,input='',timeout=20)
        (d/'pdb2gmx.log').write_text(run.stdout)
        row={'fixture':name,'ph':ph,'mode':mode,'variants':details,'hydrogen_name_mapping':hmap,'command':cmd[1:],'exit_code':run.returncode}
        if run.returncode==0:
            before,border=atoms(inp);after,aorder=atoms(d/'processed.pdb');common=before.keys()&after.keys()
            aliases={k:k.replace(':OXT',':OC1') if k.endswith(':OXT') else k.replace(':O',':OC2') if k.endswith(':O') and k not in after else k for k in before}
            all_identity_preserved=set(aliases.values())==set(after)
            all_displacements=[float(np.linalg.norm(before[k]-after[v])) for k,v in aliases.items() if v in after]
            row.update(all_atoms_preserved_after_terminal_aliases=all_identity_preserved,max_displacement_after_aliases_angstrom=max(all_displacements),input_atoms=len(before),output_atoms=len(after),missing=sorted(before.keys()-after.keys()),added=sorted(after.keys()-before.keys()),max_common_displacement_angstrom=max(float(np.linalg.norm(before[k]-after[k])) for k in common),reordered=border!=aorder)
        else:row['error_tail']=run.stdout[run.stdout.find('Fatal error:'):run.stdout.find('For more information and tips')] if 'Fatal error:' in run.stdout else run.stdout[-1500:]
        results.append(row)
        print(name,mode,run.returncode,flush=True)
(OUT.parent/'preparation-gromacs-compatibility.json').write_text(json.dumps({'scope':'Upstream parser compatibility on modeled tiny peptide fixtures, no MD. No general state preservation guarantee.','results':results},indent=2)+'\n')
