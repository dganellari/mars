"""GPT/Codex, 2026-09-11: host audit replica, public procedural channel only.

Reconstruct frozen outlet operators and the linear (advection-disabled) time map.
Pressure is stored as z = dt_eff*p/rho to avoid extreme dimensional scaling.
This independently assembled NumPy model does not execute production CUDA kernels.
See docs/design/gpt_outlet_temporal_audit_2026-09-11.md for source correspondence.
"""
import argparse
import numpy as np


def mesh():
    nx, ny, nz = 16, 4, 4
    xyz = np.array([(i/4, j/4, k/4) for k in range(5) for j in range(5) for i in range(17)])
    plane = 85
    offsets = np.array([0, 1, 18, 17, plane, plane+1, plane+18, plane+17])
    kuhn = np.array([[0,1,2,6],[0,2,3,6],[0,3,7,6],[0,7,4,6],[0,4,5,6],[0,5,1,6]])
    bases = np.array([k*plane+j*17+i for k in range(nz) for j in range(ny) for i in range(nx)])
    cells = (bases[:,None,None]+offsets[kuhn][None,:,:]).reshape(-1,4)
    coords = xyz[cells]
    jac = np.stack([coords[:,k]-coords[:,0] for k in (1,2,3)],axis=2)
    det = np.linalg.det(jac)
    grad = np.zeros((len(cells),4,3))
    grad[:,1:,:] = np.linalg.inv(jac)
    grad[:,0,:] = -grad[:,1:,:].sum(axis=1)
    assert np.all(det > 0)
    return xyz, cells, det, grad


class Operators:
    def __init__(self, dt_eff, nu, beta, trace_mode='lagged', wall_precedence=False):
        xyz, cells, det, grad = mesh()
        n = len(xyz)
        self.n, self.xyz, self.cells, self.grad = n, xyz, cells, grad
        self.dt_eff = dt_eff
        self.inlet = xyz[:,0] == 0
        self.outlet = xyz[:,0] == 4
        wall = np.any((xyz[:,1:] == 0) | (xyz[:,1:] == 1),axis=1)
        fixed = (wall | self.inlet) if wall_precedence else (wall | self.inlet) & ~self.outlet
        self.fixed = fixed
        self.free = np.flatnonzero(~fixed)
        self.free3 = np.flatnonzero(np.repeat(~fixed,3))
        self.fixed3 = np.flatnonzero(np.repeat(fixed,3))
        mass, stiffness = np.zeros(n), np.zeros((n,n))
        bi, bo, binlet = (np.zeros((n,3*n)) for _ in range(3))
        edges = []
        faces = []
        weights = np.zeros(n)
        for e, nodes in enumerate(cells):
            vol = det[e]/6
            mass[nodes] += vol/4
            stiffness[np.ix_(nodes,nodes)] += vol*grad[e]@grad[e].T
            for l,r in ((0,1),(1,2),(0,2),(0,3),(1,3),(2,3)):
                a = det[e]*(grad[e,r]-grad[e,l])/24
                left,right = nodes[l],nodes[r]
                edges.append((e,left,right,a))
                for row,sign in ((left,1),(right,-1)):
                    for col in (left,right):
                        bi[row,3*col:3*col+3] += sign*.5*a
            for face in ((0,1,3),(1,2,3),(0,3,2),(0,2,1)):
                fn = nodes[list(face)]
                if not (np.all(xyz[fn,0] == 0) or np.all(xyz[fn,0] == 4)):
                    continue
                a = .5*np.cross(xyz[fn[1]]-xyz[fn[0]],xyz[fn[2]]-xyz[fn[0]])
                is_out = xyz[fn[0],0] == 4
                for row in fn:
                    (bo if is_out else binlet)[row,3*row:3*row+3] += a/3
                    if is_out:
                        weights[row] += np.linalg.norm(a)/3
                if is_out:
                    opposite = next(i for i in range(4) if i not in face)
                    faces.append((e,fn,face,opposite,a))
        self.mass, self.stiffness = mass, stiffness
        self.edges, self.faces = edges, faces
        self.b = bi+bo+binlet
        self.bfree = self.b[:,self.free3]
        self.gv = -bi.T/np.repeat(mass,3)[:,None]
        self.gt = bo.T/np.repeat(mass,3)[:,None]
        self.gt[self.fixed3,:] = 0
        self.gv -= self.gt
        self.gfree = self.gv[self.free3,:]
        f = np.diag(self.outlet.astype(float))-np.outer(self.outlet.astype(float),weights/weights.sum())
        self.trace = (1-beta)*f if trace_mode == 'lagged' else np.zeros_like(f)
        self.gpred = self.gv+self.gt@self.trace
        # tau/h, including the velocity-Dirichlet fallback and relax_u=0.3.
        c = .3/(1+dt_eff*nu*np.diag(stiffness)/mass)
        c[fixed] = .3
        cp, ct, smooth = np.zeros((n,n)),np.zeros((n,n)),np.zeros((n,3*n))
        for e,left,right,a in edges:
            coeff = .5*(c[left]+c[right])
            compact = -coeff*grad[e]@a
            nodes = cells[e]
            for row,sign in ((left,1),(right,-1)):
                cp[row,nodes] += sign*compact
                for col in (left,right):
                    if not fixed[col]:
                        smooth[row,3*col:3*col+3] += sign*coeff*.5*a
        for e,fn,face,opposite,a in faces:
            coeff = np.mean(c[fn])
            opp = cells[e,opposite]
            for row in fn:
                cp[row,opp] -= coeff*(grad[e,opposite]@a)/3
                for local,col in zip(face,fn):
                    ct[row,col] -= coeff*(grad[e,local]@a)/3
                    smooth[row,3*col:3*col+3] += coeff*a/18
                smooth[row,3*opp:3*opp+3] += coeff*a/6
        self.cp, self.ct, self.smooth = cp,ct,smooth
        self.implicit_trace = trace_mode == 'implicit'
        derivative = cp.copy()
        if self.implicit_trace:
            # Counterfactual audit: both gradient and flux get delta trace=F*phi.
            self.trace = (1-beta)*f
            self.gpred = self.gv+self.gt@self.trace
            self.gfree = self.gpred[self.free3,:]
            derivative += ct@self.trace
        self.j = -self.bfree@self.gfree+derivative
        self.jinv = np.linalg.inv(self.j)
        self.history_pressure = smooth@self.gpred+cp+ct@self.trace
        free = self.free
        fixed_nodes = np.flatnonzero(fixed)
        self.fixed_nodes = fixed_nodes
        a = np.diag(mass/dt_eff)+nu*stiffness
        aff = a[np.ix_(free,free)]
        self.hfree = np.linalg.solve(aff,np.diag(mass[free]/dt_eff))
        self.lift = np.linalg.solve(aff,-a[np.ix_(free,fixed_nodes)])
        self.inverse_error = np.max(np.abs(self.j@self.jinv-np.eye(n)))

    def step(self, u, um, z, target, bdf):
        w = u.copy() if bdf == 1 else (4*u-um)/3
        star = w-(self.gpred@z).reshape(self.n,3)
        star[self.fixed_nodes] = target[self.fixed_nodes]
        diff = star.copy()
        diff[self.free] = self.hfree@star[self.free]+self.lift@target[self.fixed_nodes]
        rhs = self.b@diff.ravel()+self.history_pressure@z
        phi = -self.jinv@rhs
        corrected = diff.copy()
        corrected[self.free] -= (self.gfree@phi).reshape(-1,3)
        pn = z+phi
        residual = self.b@corrected.ravel()+self.cp@pn+self.smooth@self.gpred@z
        residual += self.ct@self.trace@(pn if self.implicit_trace else z)
        return corrected,pn,dict(star=float(np.max(np.abs(star))),diff=float(np.max(np.abs(diff))),
                                 corrected=float(np.max(np.abs(corrected))),
                                 speed_max=float(np.max(np.linalg.norm(corrected,axis=1))),
                                 scaled_pressure=float(np.max(np.abs(pn))),
                                 rms=float(np.sqrt(np.sum(residual**2/self.mass)/np.sum(self.mass))))

    def amplification(self):
        nf=len(self.free3)
        h=np.kron(self.hfree,np.eye(3))
        jbi=self.jinv@self.bfree
        pressure_old=np.eye(self.n)-self.jinv@self.history_pressure
        velocity_projector=np.eye(nf)+self.gfree@jbi
        velocity_old=self.gfree@self.jinv@self.history_pressure
        pred_p=self.gpred[self.free3,:]
        return np.block([
            [4/3*velocity_projector@h, -1/3*velocity_projector@h,
             -velocity_projector@h@pred_p+velocity_old],
            [np.eye(nf),np.zeros((nf,nf)),np.zeros((nf,self.n))],
            [-4/3*jbi@h,1/3*jbi@h,jbi@h@pred_p+pressure_old]])


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument('--dt',type=float,default=2e-6)
    parser.add_argument('--nu',type=float,default=1e-4)
    parser.add_argument('--rho',type=float,default=1000.)
    parser.add_argument('--beta',type=float,default=.05)
    parser.add_argument('--trace',choices=['lagged','fixed','implicit'],default='lagged')
    parser.add_argument('--steps',type=int,default=200)
    parser.add_argument('--spectrum',action='store_true')
    parser.add_argument('--wall-precedence', action='store_true',
                        help='host counterfactual: keep wall/outlet intersections velocity-fixed')
    args=parser.parse_args()
    first=Operators(args.dt,args.nu,args.beta,args.trace,args.wall_precedence)
    later=Operators(2*args.dt/3,args.nu,args.beta,args.trace,args.wall_precedence)
    print('REPLICA: linear momentum, advection disabled, NumPy direct pressure solve',flush=True)
    print('nodes',first.n,'free',len(first.free),'volume',sum(first.mass),
          'nu*K_trace',args.nu*np.trace(first.stiffness),'inverse_errors',first.inverse_error,later.inverse_error,flush=True)
    if args.spectrum:
        mat=later.amplification()
        # Check the derived matrix against the independently executed timestep.
        rng=np.random.default_rng(19)
        u=np.zeros((first.n,3)); um=u.copy(); z=rng.normal(size=first.n)
        u[later.free]=rng.normal(size=(len(later.free),3))
        um[later.free]=rng.normal(size=(len(later.free),3))
        un,zn,_=later.step(u,um,z,np.zeros_like(u),2)
        action=mat@np.r_[u[later.free].ravel(),um[later.free].ravel(),z]
        direct=np.r_[un[later.free].ravel(),u[later.free].ravel(),zn]
        action_error=np.linalg.norm(action-direct)/np.linalg.norm(direct)
        if not np.isfinite(action_error) or action_error > 1e-12:
            raise RuntimeError('amplification matrix does not reproduce the timestep')
        print('time_map_action_error',action_error,flush=True)
        eig,vectors=np.linalg.eig(mat)
        order=np.argsort(np.abs(eig))[::-1]
        leading=vectors[:,order[0]]
        eigen_error=np.linalg.norm(mat@leading-eig[order[0]]*leading)/np.linalg.norm(mat@leading)
        if not np.isfinite(eigen_error) or eigen_error > 1e-10:
            raise RuntimeError('leading eigenpair failed its residual check')
        print('leading_eigenpair_residual',eigen_error,flush=True)
        print('spectral_radius',abs(eig[order[0]]),'unstable_modes',sum(abs(eig)>1+1e-8),
              'largest',eig[order[:8]],flush=True)
    u=np.zeros((first.n,3)); um=u.copy(); z=np.zeros(first.n)
    for step in range(1,args.steps+1):
        op=first if step==1 else later
        if step==2: z*=2/3
        target=np.zeros_like(u); target[op.inlet,0]=.5*min(1,step/100)
        un,zn,stats=op.step(u,um,z,target,1 if step==1 else 2)
        stats['p_max']=stats['scaled_pressure']*args.rho/op.dt_eff
        if step<=10 or step%10==0:
            print(step,stats,flush=True)
        um,u,z=u,un,zn
        if not np.all(np.isfinite(u)) or np.max(np.abs(u))>1e100: break


if __name__=='__main__': main()
