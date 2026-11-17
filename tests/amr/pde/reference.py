"""Independent dense P1 solve on host-exported public tet meshes; requires NumPy."""
import json, sys
from pathlib import Path
import numpy as np


def solve(mesh, quadratic=False, box=False):
    points=np.array(mesh['points']); cells=np.array(mesh['cells']); n=len(points)
    exact=np.sum(points**2,axis=1) if quadratic else 1+points@np.array([1.,2.,3.])
    stiffness=np.zeros((n,n)); load=np.zeros(n); volumes=[]
    fixed=np.zeros(n,dtype=bool)
    for e,f in mesh['boundary']:
        fixed[np.delete(cells[e],f)]=True
    missed=int(np.count_nonzero(fixed & ~np.any(np.isclose(points,points.min(0))|np.isclose(points,points.max(0)),axis=1)))
    if box:fixed=np.any(np.isclose(points,points.min(0))|np.isclose(points,points.max(0)),axis=1)
    for cell in cells:
        x=points[cell]; det=np.linalg.det((x[1:]-x[0]).T); assert det>0
        volume=det/6; volumes.append(volume)
        gradients=np.linalg.inv(np.column_stack((np.ones(4),x)))[1:,:].T
        stiffness[np.ix_(cell,cell)]+=volume*(gradients@gradients.T)
        load[cell]+=(-6 if quadratic else 0)*volume/4
    edges={tuple(sorted((int(c[a]),int(c[b])))) for c in cells for a in range(4) for b in range(a+1,4)}
    parents={node:(a,b) for node,a,b in mesh.get('ancestry',[]) if tuple(sorted((a,b))) in edges}
    used=set(map(int,cells.flat))
    masters=sorted(used-set(parents)); index={node:i for i,node in enumerate(masters)}
    P=np.zeros((n,len(masters)))
    def expand(node,seen=()):
        assert node not in seen, 'Cyclic constraints'
        if node in index:return {index[node]:1.}
        a,b=parents[node]; row={}
        for parent in (a,b):
            for col,w in expand(parent,seen+(node,)).items():row[col]=row.get(col,0)+.5*w
        return row
    for node in used:
        for col,w in expand(node).items():P[node,col]=w
    fixed_dof=(P.T@fixed.astype(float))>0
    exact_dof=exact[masters]; free=~fixed_dof; dof_solution=exact_dof.copy()
    constrained=P.T@stiffness@P; rhs=P.T@load
    if np.any(free):
        A=constrained[np.ix_(free,free)]; b=rhs[free]-constrained[np.ix_(free,fixed_dof)]@exact_dof[fixed_dof]
        # Cholesky independently checks the eliminated matrix is positive definite.
        factor=np.linalg.cholesky(A)
        dof_solution[free]=np.linalg.solve(factor.T,np.linalg.solve(factor,b))
        assert np.linalg.norm(A@dof_solution[free]-b) <= 1e-11*max(1.,np.linalg.norm(b))
    solution=P@dof_solution
    q,w=np.polynomial.legendre.leggauss(4);q=(q+1)/2;w=w/2
    phi=[];weights=[]
    for i in range(4):
        for j in range(4):
            for k in range(4):
                b,c,d=q[i],(1-q[i])*q[j],(1-q[i])*(1-q[j])*q[k]
                phi.append([1-b-c-d,b,c,d]);weights.append(w[i]*w[j]*w[k]*(1-q[i])**2*(1-q[j]))
    phi=np.array(phi);weights=np.array(weights)
    error=0
    for cell,volume in zip(cells,volumes):
        x=phi@points[cell]; target=np.sum(x*x,axis=1) if quadratic else 1+x@np.array([1.,2.,3.])
        error+=6*volume*np.dot(weights,(phi@solution[cell]-target)**2)
    return dict(elements=len(cells),nodes=len(used),dofs=len(masters),hanging=len(used & set(parents)),free=int(sum(free)),l2=float(np.sqrt(error)),
                maximum=float(np.max(np.abs(solution[list(used)]-exact[list(used)]))),volume=float(sum(volumes)),non_box_boundary=missed)


def main(prefix):
    records=[]
    for level in range(5):
        mesh=json.loads(Path(prefix+str(level)+'.json').read_text())
        affine=solve(mesh); quadratic=solve(mesh,True); old=solve(mesh,False,True)
        assert affine['maximum']<1e-10 and affine['non_box_boundary']>0
        assert old['maximum']>1e-3, 'Old box detector should fail this non-box control'
        assert abs(quadratic['volume']-solve(json.loads(Path(prefix+'0.json').read_text()),True)['volume'])<1e-12
        records.append(dict(level=level,affine=affine,quadratic=quadratic,old_box=old))
    assert records[-1]['quadratic']['free']>0
    assert records[-1]['quadratic']['l2']<.5*records[0]['quadratic']['l2']
    print(json.dumps(records,indent=2))
    return records

def red_reference(prefix):
    # These two fixture rounds are 1-irregular. The later boundary-only fixture
    # deliberately omits generation closure and is not a PDE reference.
    records=[]
    for level in range(2):
        mesh=json.loads(Path(prefix+str(level)+'.json').read_text())
        affine=solve(mesh); quadratic=solve(mesh,True)
        assert affine['maximum']<1e-10 and affine['hanging']>0
        assert abs(affine['volume']-6.431416353311024)<1e-10
        records.append(dict(level=level,affine=affine,quadratic=quadratic))
    assert records[-1]['quadratic']['free']>0
    assert records[-1]['quadratic']['l2']<.5*records[0]['quadratic']['l2']
    print(json.dumps(records,indent=2))
    return records

if __name__=='__main__':
    (red_reference if '--red' in sys.argv[2:] else main)(sys.argv[1])
