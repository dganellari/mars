"""Independent Q3 conduction checks with mixed boundary conditions and rotated hexes."""
import itertools
import json
import numpy as np

BITS=np.array([[0,0,0],[1,0,0],[1,1,0],[0,1,0],[0,0,1],[1,0,1],[1,1,1],[0,1,1]])
GLL=np.array([-1,-1/np.sqrt(5),1/np.sqrt(5),1])


def rotations():
    out=[]
    for axes in itertools.permutations(range(3)):
        for signs in itertools.product((-1,1),repeat=3):
            a=np.eye(3,dtype=int)[list(axes)]*np.array(signs)[:,None]
            if round(np.linalg.det(a))==1:out.append(a)
    return out


def basis(points):
    b=np.ones((len(points),4));d=np.zeros_like(b)
    for i in range(4):
        for j in range(4):
            if j!=i:b[:,i]*=(points-GLL[j])/(GLL[i]-GLL[j])
        for j in range(4):
            if j==i:continue
            term=np.full(len(points),1/(GLL[i]-GLL[j]))
            for k in range(4):
                if k not in (i,j):term*=(points-GLL[k])/(GLL[i]-GLL[k])
            d[:,i]+=term
    return b,d


def mapped(corners, points):
    signs=2*BITS-1
    factors=(1+points[:,None,:]*signs[None,:,:])/2
    n=factors.prod(2)
    grad=np.empty((len(points),8,3))
    for a in range(3):grad[:,:,a]=signs[None,:,a]/2*np.prod(factors[:,:,[d for d in range(3) if d!=a]],axis=2)
    return n@corners,np.einsum('ka,qkb->qab',corners,grad)


def assemble(warp,offset):
    q,w=np.polynomial.legendre.leggauss(5);b,d=basis(q)
    xyzq=np.array(list(itertools.product(q,repeat=3)))
    weights=np.array([w[i]*w[j]*w[k] for i,j,k in itertools.product(range(5),repeat=3)])
    grad=np.stack([np.einsum('ia,jb,kc->ijkabc',*(d if a==axis else b for a in range(3))).reshape(125,64) for axis in range(3)],axis=2)
    points=np.array(list(itertools.product(GLL,repeat=3)))
    n=7;matrix=np.zeros((n**3,n**3));coord=np.full((n**3,3),np.nan);locals=[]
    for e,ijk in enumerate(itertools.product(range(2),repeat=3)):
        rotation=rotations()[(offset+e)%24]
        bits=((2*BITS-1)@rotation.T+1)//2
        c=(bits+ijk)/2
        bump=warp*np.prod(np.sin(np.pi*c),axis=1)
        c=c+bump[:,None]*np.array([1,.4,-.3])
        x,jac=mapped(c,xyzq);det=np.linalg.det(jac)
        assert np.all(det>0)
        physical_grad=np.einsum('qia,qab->qib',grad,np.linalg.inv(jac))
        local=np.einsum('qia,qja,q->ij',physical_grad,physical_grad,weights*det)
        indices=[]
        for node in itertools.product(range(4),repeat=3):
            original=((2*np.array(node)-3)@rotation.T+3)//2
            g=np.array(ijk)*3+original;indices.append((g[0]*n+g[1])*n+g[2])
        indices=np.array(indices);px,_=mapped(c,points)
        old=coord[indices];present=np.isfinite(old[:,0]);np.testing.assert_allclose(old[present],px[present],atol=1e-14)
        coord[indices]=px
        matrix[np.ix_(indices,indices)]+=local;locals.append((indices,local))
    return matrix,coord,locals


def solve_case(warp,offset,wrong_sign=False):
    a,x,local=assemble(warp,offset)
    hot=np.isclose(x[:,0],0,atol=1e-13);cold=np.isclose(x[:,0],1,atol=1e-13);free=~(hot|cold)
    lift=hot.astype(float);rhs=(-1 if not wrong_sign else 1)*(a@lift)
    correction=np.zeros(len(x));correction[free]=np.linalg.solve(a[np.ix_(free,free)],rhs[free])
    theta=correction+lift;r=a@theta
    # The element-local pairings must produce the same reactions and energy as assembled rows.
    pair=np.zeros(3)
    for ids,k in local:
        y=k@theta[ids];pair += [hot[ids]@y,cold[ids]@y,theta[ids]@y]
    np.testing.assert_allclose(pair,[hot@r,cold@r,theta@r],atol=1e-13)
    error=float(np.max(np.abs(theta-(1-x[:,0]))))
    if wrong_sign:
        assert error>.1
        return dict(negative_control_error=error)
    np.testing.assert_allclose(theta,1-x[:,0],atol=2e-12,rtol=0)
    np.testing.assert_allclose(pair,[1,-1,1],atol=3e-12,rtol=0)
    assert np.max(np.abs(r[free]))<1e-13
    assert theta.min()>-1e-12 and theta.max()<1+1e-12
    return dict(warp=warp,rotation_offset=offset,nodes=len(x),max_error=error,reactions=pair.tolist(),free_residual=float(np.max(np.abs(r[free]))))


if __name__=='__main__':
    results=[solve_case(w,o) for w in (0,.05) for o in (0,8,16)]
    results.append(solve_case(.05,0,True))
    print(json.dumps(dict(status='PASS',cases=results),indent=2))
