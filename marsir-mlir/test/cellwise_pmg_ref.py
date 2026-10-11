# Reference for p-multigrid in the cell-wise CVFEM solve (design:
# docs/design/cellwise_unstructured_dss.md, "Multigrid on unstructured meshes"):
# levels p = 7 -> 3 -> 1 on the same mesh, the Knaus CVFEM operator rediscretized at each
# degree, Chebyshev-Jacobi smoothing through the DSS, element-local GLL interpolation
# between degrees (restriction = its transpose on the raw residual), and at p = 1 either
# an exact solve of the assembled system (an ideal AMG) or k Chebyshev steps.
#
# Run: python3 cellwise_pmg_ref.py ne:degrees:pre:post:coarse ...
#   e.g. 8:7,3,1:0:3:exact   16:7,3,1:0:3:10   (coarse = exact, or Chebyshev steps)
import sys, time
import numpy as np
from numpy.polynomial import legendre as Lg

def lag(nodes, x):
    m = len(nodes); v = np.ones(m); d = np.zeros(m)
    for j in range(m):
        others = [k for k in range(m) if k != j]
        den = np.prod([nodes[j] - nodes[k] for k in others])
        v[j] = np.prod([x - nodes[k] for k in others]) / den
        d[j] = sum(np.prod([x - nodes[k] for k in others if k != i]) for i in others) / den
    return v, d

class Basis:
    def __init__(self, p):
        self.p = p; n = self.n = p + 1; self.P = p
        c = np.zeros(p + 1); c[p] = 1
        self.zeta = np.sort(np.concatenate(([-1.0], Lg.legroots(Lg.legder(c)), [1.0])))
        self.xi, _ = Lg.leggauss(p)
        self.Btil = np.array([lag(self.zeta, x)[0] for x in self.xi])
        self.Dtil = np.array([lag(self.zeta, x)[1] for x in self.xi])
        self.D = np.array([lag(self.zeta, x)[1] for x in self.zeta])
        pad = np.concatenate(([-1.0], self.xi, [1.0])); Winv = np.zeros((n, n))
        for i in range(n):
            _, d = lag(pad, self.zeta[i]); Winv[i] = -np.cumsum(d[:n])
        self.W = np.linalg.inv(Winv)

def apply(u, G, B):
    n, P, nn = B.n, B.P, B.n * B.n
    u = u.reshape(-1, n, n, n); y = np.zeros_like(u)
    for d in range(3):
        U = np.moveaxis(u, 1 + d, 1).reshape(-1, n, nn); Y = np.moveaxis(y, 1 + d, 1)
        interp = (B.Btil @ U).reshape(-1, P, n, n); deriv = (B.Dtil @ U).reshape(-1, P, n, n)
        g = G[:, d]
        flux = g[:, :, 2] * deriv + g[:, :, 0] * (interp @ B.D.T) + g[:, :, 1] * (B.D @ interp)
        intf = B.W @ (flux @ B.W.T)
        Y[:, :P] -= intf; Y[:, 1:] += intf
    return y.reshape(-1, n ** 3)

t1 = (1, 0, 0); t2 = (2, 2, 1)
corner_ref = np.array([[-1, -1, -1], [1, -1, -1], [1, 1, -1], [-1, 1, -1],
                       [-1, -1, 1], [1, -1, 1], [1, 1, 1], [-1, 1, 1]])

def shape(rf):
    f = 1 + corner_ref * rf[..., None, :]; dN = np.empty(f.shape)
    for k in range(3):
        g = f.copy(); g[..., k] = corner_ref[:, k]; dN[..., k] = 0.125 * g.prod(-1)
    return 0.125 * f.prod(-1), dN

def metric(C, B):
    n, P = B.n, B.P; G = np.empty((len(C), 3, P, 3, n, n))
    for d in range(3):
        rf = np.zeros((P, n, n, 3)); rf[..., d] = B.xi[:, None, None]
        rf[..., t1[d]] = B.zeta[None, :, None]; rf[..., t2[d]] = B.zeta[None, None, :]
        J = np.einsum("eca,lsrcb->elsrab", C, shape(rf)[1]); Ji = np.linalg.inv(J)
        gv = np.linalg.det(J)[..., None] * np.einsum("...ak,...k->...a", Ji, Ji[..., d, :])
        G[:, d, :, 0] = gv[..., t2[d]]; G[:, d, :, 1] = gv[..., t1[d]]; G[:, d, :, 2] = gv[..., d]
    return G

def diagonal(G, B):
    n, P, nn, W, D, Bt, Dt = B.n, B.P, B.n * B.n, B.W, B.D, B.Btil, B.Dtil
    dg = np.zeros((len(G), n ** 3))
    for node in range(n ** 3):
        a, b, cc = node // nn, (node // n) % n, node % n
        for d in range(3):
            A = (a, b, cc)[d]; Bq = b if d == 0 else a; C = b if d == 2 else cc
            for l, sign in ((A - 1, 1.0), (A, -1.0)):
                if not 0 <= l < P: continue
                g0, g1, g2 = G[:, d, l, 0], G[:, d, l, 1], G[:, d, l, 2]
                term = W[Bq, Bq] * W[C, C] * g2[:, Bq, C] * Dt[l, A]
                term = term + W[Bq, Bq] * Bt[l, A] * (g0[:, Bq, :] @ (W[C] * D[:, C]))
                term = term + W[C, C] * Bt[l, A] * (g1[:, :, C] @ (W[Bq] * D[:, Bq]))
                dg[:, node] += sign * term
    return dg

class Level:
    def __init__(self, ne, deform, B):
        self.ne, self.E, self.B = ne, ne ** 3, B; n, p = B.n, B.p
        lat = np.arange(ne + 1) / ne
        V = np.stack(np.meshgrid(lat, lat, lat, indexing="ij"), -1)
        V = V + deform * np.prod(np.sin(np.pi * V), -1, keepdims=True)
        eidx = np.indices((ne, ne, ne)).reshape(3, -1).T
        self.corners = np.stack([V[tuple((eidx + (s > 0)).T)] for s in corner_ref], 1)
        self.G = metric(self.corners, B)
        la, lb, lc = np.indices((n, n, n)).reshape(3, -1)
        outer = lambda e, l: ((e[:, None] == 0) & (l[None, :] == 0)) | ((e[:, None] == ne - 1) & (l[None, :] == p))
        self.mask = (~(outer(eidx[:, 0], la) | outer(eidx[:, 1], lb) | outer(eidx[:, 2], lc))).astype(float)
        self.mult = self.dss(np.ones((self.E, n ** 3)))
        self.diag = self.dss(diagonal(self.G, B))
        rf = np.stack([B.zeta[la], B.zeta[lb], B.zeta[lc]], -1)
        self.x = np.einsum("jc,eca->eja", shape(rf)[0], self.corners)

    def dss(self, fc):
        ne, n = self.ne, self.B.n
        v = fc.reshape(ne, ne, ne, n, n, n).copy()
        for ax in range(3):
            lo = [slice(None)] * 6; hi = [slice(None)] * 6
            lo[ax] = slice(0, ne - 1); hi[ax] = slice(1, ne); lo[3 + ax] = n - 1; hi[3 + ax] = 0
            s = v[tuple(lo)] + v[tuple(hi)]; v[tuple(lo)] = s; v[tuple(hi)] = s
        return v.reshape(self.E, n ** 3)

    def A(self, u): return apply(u, self.G, self.B)
    def PJ(self, r): return self.mask * self.dss(r) / self.diag
    def wdot(self, u, v): return float(np.sum(u * v / self.mult))

def transfer(Bf, Bc): return np.array([lag(Bc.zeta, z)[0] for z in Bf.zeta])   # n_f x n_c

def prolong(uc, I):
    nf, nc = I.shape
    return np.einsum("aj,bk,cl,ejkl->eabc", I, I, I, uc.reshape(-1, nc, nc, nc), optimize=True).reshape(-1, nf ** 3)

def restrict(rf, I):
    nf, nc = I.shape
    return np.einsum("aj,bk,cl,eabc->ejkl", I, I, I, rf.reshape(-1, nf, nf, nf), optimize=True).reshape(-1, nc ** 3)

class CoarseP1:
    """Exact solve of the assembled p = 1 system (an ideal AMG)."""
    def __init__(self, lev):
        ne, E = lev.ne, lev.E; V1 = ne + 1
        Ae = np.zeros((E, 8, 8))
        for j in range(8):
            u = np.zeros((E, 8)); u[:, j] = 1.0; Ae[:, :, j] = lev.A(u)
        eidx = np.indices((ne, ne, ne)).reshape(3, -1).T
        la, lb, lc = np.indices((2, 2, 2)).reshape(3, -1)
        self.gid = ((eidx[:, 0, None] + la) * V1 + eidx[:, 1, None] + lb) * V1 + eidx[:, 2, None] + lc
        self.NV = V1 ** 3
        A = np.zeros((self.NV, self.NV)); np.add.at(A, (self.gid[:, :, None], self.gid[:, None, :]), Ae)
        vi = np.indices((V1, V1, V1)).reshape(3, -1).T
        self.idx = np.where(np.all((vi > 0) & (vi < ne), axis=1))[0]
        self.Ainv = np.linalg.inv(A[np.ix_(self.idx, self.idx)])
    def __call__(self, r):
        rhs = np.zeros(self.NV); np.add.at(rhs, self.gid, r)
        x = np.zeros(self.NV); x[self.idx] = self.Ainv @ rhs[self.idx]
        return x[self.gid]

def hash_vector(size):
    t = np.arange(size, dtype=np.uint64)
    return ((t * np.uint64(2654435761)) & np.uint64(0xffffffff)).astype(float) / 4294967296.0 - 0.5

def power_iteration(lev, steps=30):
    v = lev.PJ(hash_vector(lev.mask.size).reshape(lev.mask.shape))
    for _ in range(steps):
        w = lev.PJ(lev.A(v)); nw = np.sqrt(lev.wdot(w, w)); lam = nw / np.sqrt(lev.wdot(v, v)); v = w / nw
    return lam

class PMultigrid:
    def __init__(self, ne, deform, degrees, pre, post, coarse="exact", rng=15.0):
        self.B = [Basis(p) for p in degrees]
        self.levels = [Level(ne, deform, b) for b in self.B]
        self.I = [transfer(self.B[l], self.B[l + 1]) for l in range(len(degrees) - 1)]
        self.pre, self.post, self.range = pre, post, rng
        self.lmax = [power_iteration(l) for l in self.levels]
        self.coarse = CoarseP1(self.levels[-1]) if coarse == "exact" else None
        self.coarse_steps = 0 if coarse == "exact" else int(coarse)

    def smooth(self, l, b, steps, x=None):
        lev = self.levels[l]; lmax = 1.1 * self.lmax[l]; lmin = lmax / self.range
        theta, delta = 0.5 * (lmax + lmin), 0.5 * (lmax - lmin); sigma = theta / delta; rho = 1.0 / sigma
        d = lev.PJ(b if x is None else b - lev.A(x)) / theta; x = d if x is None else x + d
        for _ in range(steps - 1):
            rho_new = 1.0 / (2.0 * sigma - rho)
            d = rho_new * rho * d + (2.0 * rho_new / delta) * lev.PJ(b - lev.A(x)); x = x + d; rho = rho_new
        return x

    def vcycle(self, b, l=0):
        lev = self.levels[l]
        if l == len(self.levels) - 1:
            return self.coarse(lev.mask * b) if self.coarse else self.smooth(l, b, self.coarse_steps)
        if self.pre:
            x = self.smooth(l, b, self.pre); r = lev.mask * (b - lev.A(x))
        else:
            r = lev.mask * b
        e = prolong(self.vcycle(restrict(r, self.I[l]), l + 1), self.I[l])
        x = x + e if self.pre else e
        return self.smooth(l, b, self.post, x) if self.post else x

def bicgstab(lev, P, b, tol=1e-10, maxit=300):
    A, dot = lev.A, lev.wdot
    x = np.zeros_like(b); r = P(b); rh = r.copy(); r0 = np.sqrt(dot(r, r))
    rho0 = al = om = 1.0; v = np.zeros_like(b); pp = np.zeros_like(b)
    for it in range(1, maxit + 1):
        rho = dot(rh, r); be = (rho / rho0) * (al / om); rho0 = rho
        pp = r + be * (pp - om * v); v = P(A(pp)); al = rho / dot(rh, v)
        s = r - al * v; t = P(A(s)); om = dot(t, s) / dot(t, t)
        x = x + al * pp + om * s; r = s - om * t
        if np.sqrt(dot(r, r)) <= tol * r0: return x, it
    return x, maxit

def run(ne, deform, degrees, pre, post, coarse="exact"):
    t = time.time()
    mg = PMultigrid(ne, deform, degrees, pre, post, coarse)
    f = mg.levels[0]
    u = f.mask * np.prod(np.sin(np.pi * f.x), -1)
    x, it = bicgstab(f, lambda r: mg.vcycle(r), f.A(u))
    err = np.sqrt(f.wdot(x - u, x - u) / f.wdot(u, u))
    print(f"ne={ne:3d} deform={deform} degrees={degrees} ({pre},{post}) coarse={coarse}: "
          f"{it:3d} it, err {err:.1e}, lmax {[round(l, 3) for l in mg.lmax[:-1]]} ({time.time() - t:.0f} s)", flush=True)

if __name__ == "__main__":
    for spec in sys.argv[1:]:
        ne, deg, pre, post, coarse = spec.split(":")
        run(int(ne), 0.05, [int(d) for d in deg.split(",")], int(pre), int(post), coarse)
