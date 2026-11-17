#!/usr/bin/env python3
"""Compare channel refinement fields independently of corrected mass-flux diagnostics."""
import argparse
import json
from pathlib import Path
import re
import numpy as np


def require(ok, message):
    if not ok:
        raise ValueError(message)


def load(prefix):
    prefix=Path(prefix)
    paths=sorted(prefix.parent.glob(prefix.name+'.rank*.csv'))
    require(bool(paths),'missing owned fields: '+str(prefix))
    a=np.concatenate([np.atleast_1d(np.genfromtxt(p,delimiter=',',names=True)) for p in paths])
    a=np.sort(a,order='node_id')
    require(np.array_equal(a['node_id'],np.arange(1,len(a)+1)),'incomplete or duplicate node identities')
    x=np.column_stack([a[k] for k in ('x','y','z')])
    u=np.column_stack([a[k] for k in ('velocity_x','velocity_y','velocity_z','pressure')])
    require(np.isfinite(x).all() and np.isfinite(u).all(),'non-finite fields')
    cells=np.loadtxt(str(prefix)+'.cells.csv',delimiter=',',skiprows=1,ndmin=2)
    require(cells.shape[1]==5 and np.isfinite(cells).all() and np.equal(cells,np.floor(cells)).all(),
            'invalid cell identities')
    require(np.array_equal(cells[:,0],np.arange(len(cells))),'duplicate or incomplete element identities')
    cells=cells[:,1:].astype(np.int64)-1
    require(cells.min()>=0 and cells.max()<len(x),'cell references absent node')
    require((np.diff(np.sort(cells,axis=1),axis=1)>0).all(),'repeated cell corner')
    return x,u,cells


def geometry(x,cells):
    edge=x[cells][:,1:]-x[cells][:,:1]
    volume=np.abs(np.linalg.det(edge))/6
    require(np.isfinite(volume).all() and (volume>0).all(),'invalid cell volume')
    return edge,volume


def metrics(x,u,cells):
    edge,volume=geometry(x,cells)
    # Each row of edge is a coordinate difference. Solve separately for each velocity component.
    gradient=np.linalg.solve(edge,u[cells][:,1:,:3]-u[cells][:,:1,:3])
    divergence=np.trace(gradient,axis1=1,axis2=2)
    return {'volume':float(volume.sum()),
            'raw_cell_divergence_l1':float(np.sum(volume*np.abs(divergence))),
            'raw_cell_divergence_l2':float(np.sqrt(np.sum(volume*divergence**2)))}


def difference(a,b):
    x,u,cells=a
    y,v,_=b
    # Procedural grids are nested. No nearest-node interpolation or interface averaging is used.
    table={tuple(np.round(c,12)):value for c,value in zip(y,v)}
    require(len(table)==len(y),'reference has duplicate geometric nodes')
    try: d=np.array([value-table[tuple(np.round(c,12))] for c,value in zip(x,u)])
    except KeyError as e: raise ValueError('reference does not contain a comparison node') from e
    _,volume=geometry(x,cells)
    dc=d[cells]
    # Exact P1 mass integral of the nodal difference interpolant, not of the fine-grid error field.
    squared=(np.sum(dc**2,axis=1)+np.sum(dc,axis=1)**2)/20
    l2=np.sqrt(np.sum(volume[:,None]*squared,axis=0))
    return {'velocity_max':float(np.max(np.linalg.norm(d[:,:3],axis=1))),
            'velocity_interpolant_l2':float(np.linalg.norm(l2[:3])),
            'pressure_max':float(np.max(abs(d[:,3]))),'pressure_interpolant_l2':float(l2[3])}


def inspect(root):
    root=Path(root)
    cases={}
    for mode,n in [(m,n) for n in (2,4,8) for m in ('shared','matched','nonmatching')]+[('shared',16)]:
        name=f'{mode}{n}'
        a=load(root/name)
        text=(root/(name+'.log')).read_text()
        require(f'PASS accuracy solve: n={n} mode={mode}' in text,'solve did not pass: '+name)
        cases[name]=a
    results={}
    for n in (2,4,8):
        for mode in ('shared','matched','nonmatching'):
            name=f'{mode}{n}'
            ref=f'shared{n if mode=="matched" else 2*n}'
            result=metrics(*cases[name])
            result.update(reference=ref,**difference(cases[name],cases[ref]))
            text=(root/(name+'.log')).read_text()
            result['reported_velocity_jump']=float(re.findall(r'velocity_jump=([\deE+.-]+)',text)[-1])
            results[name]=result
    checks={}
    for mode in ('matched','nonmatching'):
        for metric in ('velocity_max','velocity_interpolant_l2','pressure_max','pressure_interpolant_l2',
                       'reported_velocity_jump'):
            values=[results[f'{mode}{n}'][metric] for n in (2,4,8)]
            checks[f'{mode}_{metric}_decreases']=bool(values[2]<values[1]<values[0])
    checks['fine_nonmatching_raw_divergence_decreases']=bool(
        results['nonmatching8']['raw_cell_divergence_l2'] < results['nonmatching4']['raw_cell_divergence_l2'])
    return {'schema':'mars-nonconformal-flow-accuracy-v1','passed':all(checks.values()),
            'scope':'stationary_tet4_channel_refinement_not_general_accuracy_or_stability_proof',
            'comparison':'matched_versus_same_grid_shared_nodes; nonmatching_versus_finer_shared_nodes',
            'cases':results,'checks':checks}


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('directory',type=Path)
    p.add_argument('--output',type=Path,required=True)
    a=p.parse_args()
    try: report=inspect(a.directory)
    except (ValueError,OSError,IndexError) as e: p.exit(1,f'Accuracy inspection rejected: {e}\n')
    with a.output.open('x') as f: json.dump(report,f,indent=2); f.write('\n')
    print(json.dumps({'passed':report['passed'],'checks':report['checks']},indent=2))
    raise SystemExit(0 if report['passed'] else 1)

if __name__=='__main__': main()
