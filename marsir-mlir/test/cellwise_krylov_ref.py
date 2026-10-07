#!/usr/bin/env python3
"""Reference for solving on element-local (cell-wise) vectors, after Wichrowski,
"Coalesced Matrix-Free Finite Elements in Cell-Wise Storage" (internal-notes/FlexibleCG.pdf).

The paper proves flexible CG on cell-wise data reproduces the assembled solve, but
CG needs a SYMMETRIC operator, and the Knaus CVFEM operator MARSIR generates is
not: ||K - K^T|| / ||K|| is about 11% even on a Cartesian element (printed below).
The same idea still works for a nonsymmetric operator. Run the Krylov method on
the left-preconditioned operator P A, where P = mask * DSS / diag maps an
unassembled residual to a continuous field. Every vector is then continuous, and
an inner product that weights each element-local copy by 1/(number of copies)
equals the assembled inner product. So BiCGStab (or GMRES) runs on element-local
data only, communication stays inside P, and the iterates equal those of the
assembled solve. This script checks:
  * the dimensionally split DSS cascade equals G G^T;
  * cell-wise BiCGStab and assembled BiCGStab give the same iterates;
  * both converge to a manufactured solution.
The GPU implementations (mixed MARS + MARSIR, then all-MARSIR) must reproduce it.

Run: python3 test/cellwise_krylov_ref.py
"""
import numpy as np
from numpy.polynomial import legendre as L

p = 7; n = p + 1; P = p
# GLL nodes: -1, roots of P'_p, +1
c = np.zeros(p + 1); c[p] = 1
zeta = np.sort(np.concatenate(([-1.0], L.legroots(L.legder(c)), [1.0])))
xi, _ = L.leggauss(P)
def lag(nodes, x):
    m = len(nodes); v = np.ones(m); d = np.zeros(m)
    for j in range(m):
        others = [k for k in range(m) if k != j]
        den = np.prod([nodes[j] - nodes[k] for k in others])
        v[j] = np.prod([x - nodes[k] for k in others]) / den
        d[j] = sum(np.prod([x - nodes[k] for k in others if k != i]) for i in others) / den
    return v, d
Btil = np.array([lag(zeta, x)[0] for x in xi]); Dtil = np.array([lag(zeta, x)[1] for x in xi])
D = np.array([lag(zeta, x)[1] for x in zeta])
pad = np.concatenate(([-1.0], xi, [1.0]))
Winv = np.zeros((n, n))
for i in range(n):
    _, d = lag(pad, zeta[i]); Winv[i] = -np.cumsum(d[:n])
W = np.linalg.inv(Winv)
nn = n * n
def apply(u, G):        # G[d][l] = (g0, g1, g2) planes of shape (n, n), as knaus_oracle.h
    u = u.reshape(n, n, n); y = np.zeros((n, n, n))
    for d in range(3):
        U = np.moveaxis(u, d, 0); Y = np.moveaxis(y, d, 0)
        for l in range(P):
            interp = np.einsum("q,qsr->sr", Btil[l], U); deriv = np.einsum("q,qsr->sr", Dtil[l], U)
            dt2 = np.einsum("rq,sq->sr", D, interp); dt1 = np.einsum("sq,qr->sr", D, interp)
            g0, g1, g2 = G[d][l]
            flux = g2 * deriv + g0 * dt2 + g1 * dt1
            tmp = np.einsum("rq,sq->sr", W, flux); intf = np.einsum("sq,qr->sr", W, tmp)
            Y[l] -= intf; Y[l + 1] += intf
    return y.ravel()
def matrix(G):
    return np.array([apply(np.eye(n**3)[j], G) for j in range(n**3)]).T
t1 = {0: 1, 1: 0, 2: 0}; t2 = {0: 2, 1: 2, 2: 1}
def metric_from_J(Jfun):
    G = {}
    for d in range(3):
        G[d] = {}
        for l in range(P):
            g = np.zeros((3, n, n))
            for s in range(n):
                for r in range(n):
                    rf = [0.0] * 3; rf[d] = xi[l]; rf[t1[d]] = zeta[s]; rf[t2[d]] = zeta[r]
                    J = Jfun(rf); det = np.linalg.det(J); Ji = np.linalg.inv(J)
                    Gv = det * (Ji @ Ji.T)[:, d]
                    g[0, s, r] = Gv[t2[d]]; g[1, s, r] = Gv[t1[d]]; g[2, s, r] = Gv[d]
            G[d][l] = g
    return G

K0 = matrix(metric_from_J(lambda rf: np.eye(3)))
print("CVFEM element operator: ||K - K^T|| / ||K|| =", "%.3e" % (np.linalg.norm(K0 - K0.T) / np.linalg.norm(K0)))
ne = 3; h = 1.0 / ne
Ke = matrix(metric_from_J(lambda rf: (h / 2) * np.eye(3)))     # same for every element
NG = ne * p + 1                                                  # global nodes per axis
E = ne ** 3
# gather G: element (ex,ey,ez), local (a,b,c) -> global (ex p + a, ey p + b, ez p + c)
gidx = np.zeros((E, n, n, n), dtype=np.int64)
for ex in range(ne):
    for ey in range(ne):
        for ez in range(ne):
            e = (ex * ne + ey) * ne + ez
            a = np.arange(n)
            gx, gy, gz = np.meshgrid(ex * p + a, ey * p + a, ez * p + a, indexing="ij")
            gidx[e] = (gx * NG + gy) * NG + gz
gidx = gidx.reshape(E, n ** 3)
def gather(ug): return ug[gidx]
def scatter(fc):
    out = np.zeros(NG ** 3); np.add.at(out, gidx.ravel(), fc.ravel()); return out
# DSS by the dimensionally split cascade: 3 axis passes of one-to-one face sums
def dss(fc):
    v = fc.reshape(ne, ne, ne, n, n, n).copy()
    for ax in range(3):                         # element axis ax <-> local axis ax
        lo = [slice(None)] * 6; hi = [slice(None)] * 6
        lo[ax] = slice(0, ne - 1); hi[ax] = slice(1, ne)
        lo[3 + ax] = n - 1; hi[3 + ax] = 0       # left element's last plane, right's first
        s = v[tuple(lo)] + v[tuple(hi)]
        v[tuple(lo)] = s; v[tuple(hi)] = s
    return v.reshape(E, n ** 3)
# checks: S == G G^T on random data
rng = np.random.default_rng(1)
f = rng.standard_normal((E, n ** 3))
print("DSS cascade vs G G^T:", np.abs(dss(f) - gather(scatter(f))).max())
mult = dss(np.ones((E, n ** 3)))                 # copies of each node
coords = np.arange(NG) / (NG - 1)
X, Y, Z = np.meshgrid(coords, coords, coords, indexing="ij")
bnd = ((X == 0) | (X == 1) | (Y == 0) | (Y == 1) | (Z == 0) | (Z == 1)).ravel()
mask_c = (~bnd)[gidx].astype(float)              # Dirichlet: zero boundary copies
def A_cell(uc): return uc @ Ke.T                 # local operator, primal -> dual
def wdot(u, v): return float(np.sum(u * v / mult))
diag_c = dss(np.tile(np.diag(Ke), (E, 1)))       # assembled diagonal, on every copy
def P_cell(rc): return mask_c * dss(rc) / diag_c # Jacobi-DSS: dual -> continuous primal
# assembled reference on the interior global DOFs
inner = np.where(~bnd)[0]
def A_glob(ug):
    full = np.zeros(NG ** 3); full[inner] = ug
    return scatter(A_cell(gather(full)))[inner]
dg = scatter(np.tile(np.diag(Ke), (E, 1)))[inner]
def P_glob(rg): return rg / dg
uex = (np.sin(np.pi * X) * np.sin(np.pi * Y) * np.sin(np.pi * Z)).ravel()
b_cell = A_cell(gather(uex))                     # unassembled RHS (dual)
b_glob = scatter(b_cell)[inner]
def bicgstab(A, P, b, dot, iters):
    x = np.zeros_like(b); r = P(b - A(x)); rh = r.copy()
    rho0 = al = om = 1.0; v = np.zeros_like(b); pp = np.zeros_like(b); hist = []
    for k in range(iters):
        rho = dot(rh, r); be = (rho / rho0) * (al / om); rho0 = rho
        pp = r + be * (pp - om * v); v = P(A(pp)); al = rho / dot(rh, v)
        s = r - al * v; t = P(A(s)); om = dot(t, s) / dot(t, t)
        x = x + al * pp + om * s; r = s - om * t
        hist.append((x.copy(), np.sqrt(dot(r, r))))
    return hist
hc = bicgstab(A_cell, P_cell, b_cell, wdot, 40)
hg = bicgstab(A_glob, P_glob, b_glob, lambda a, c: float(a @ c), 40)
print(" it   ||Pr||_w (cell-wise)   max|x_cell - G x_glob|   max|x - u_exact|")
for k in [0, 4, 9, 19, 29, 39]:
    xc, rc = hc[k]; xg, rg = hg[k]
    full = np.zeros(NG ** 3); full[inner] = xg
    print(f"{k+1:3d}   {rc:.3e}              {np.abs(xc - gather(full)).max():.2e}"
          f"             {np.abs(xc - gather(uex)).max():.2e}")

def to_cell(xg):   # interior global vector -> element-local copies
    full = np.zeros(NG ** 3); full[inner] = xg
    return gather(full)
ok = (np.abs(dss(f) - gather(scatter(f))).max() < 1e-12
      and max(np.abs(xc - to_cell(xg)).max() for (xc, _), (xg, _) in zip(hc, hg)) < 1e-10
      and np.abs(hc[-1][0] - gather(uex)).max() < 1e-10)
print("CELL-WISE KRYLOV REFERENCE:", "PASS" if ok else "FAIL")
