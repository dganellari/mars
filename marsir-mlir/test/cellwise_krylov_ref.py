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
  * the closed-form diagonal of the GPU driver equals the probed one;
  * cell-wise BiCGStab and assembled BiCGStab give the same iterates;
  * the solve converges to the manufactured solution.

Mesh, solution and solver are those of
examples/distributed/unstructured/mars_marsir_cellwise_solve.cu: the unit cube as
ne^3 trilinear hexahedra, optionally deformed by x += a sin(pi x) sin(pi y) sin(pi z) (1,1,1),
u = sin(pi x) sin(pi y) sin(pi z) at the GLL nodes, b = A u, stop at
||P r||_w <= tol ||P r0||_w. Its residual history must match the GPU run with the
same --ne and --deform.

Run: python3 test/cellwise_krylov_ref.py [--ne 3] [--deform 0] [--tol 1e-10] [--maxit 500]
"""
import argparse

import numpy as np
from numpy.polynomial import legendre as L

ap = argparse.ArgumentParser()
ap.add_argument("--ne", type=int, default=3)
ap.add_argument("--deform", type=float, default=0.0)
ap.add_argument("--tol", type=float, default=1e-10)
ap.add_argument("--maxit", type=int, default=500)
args = ap.parse_args()

p = 7; n = p + 1; P = p; nn = n * n; n3 = nn * n
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

# G (E, 3, P, 3, n, n) = [dir][face][g0 g1 g2][s][r], the layout of knaus_oracle.h
def apply(u, G):        # y = A u on every element; u, y are (E, n^3)
    u = u.reshape(-1, n, n, n); y = np.zeros_like(u)
    for d in range(3):
        U = np.moveaxis(u, 1 + d, 1); Y = np.moveaxis(y, 1 + d, 1)
        interp = np.einsum("lq,eqsr->elsr", Btil, U); deriv = np.einsum("lq,eqsr->elsr", Dtil, U)
        dt2 = np.einsum("rq,elsq->elsr", D, interp); dt1 = np.einsum("sq,elqr->elsr", D, interp)
        g = G[:, d]
        flux = g[:, :, 2] * deriv + g[:, :, 0] * dt2 + g[:, :, 1] * dt1
        tmp = np.einsum("rq,elsq->elsr", W, flux); intf = np.einsum("sq,elqr->elsr", W, tmp)
        Y[:, :P] -= intf; Y[:, 1:] += intf
    return y.reshape(-1, n3)

t1 = (1, 0, 0); t2 = (2, 2, 1)
corner_ref = np.array([[-1, -1, -1], [1, -1, -1], [1, 1, -1], [-1, 1, -1],
                       [-1, -1, 1], [1, -1, 1], [1, 1, 1], [-1, 1, 1]])   # c_hexCornerRef
def shape(rf):          # trilinear shape functions and reference gradients at points rf (..., 3)
    f = 1 + corner_ref * rf[..., None, :]
    dN = np.empty(f.shape)
    for k in range(3):
        g = f.copy(); g[..., k] = corner_ref[:, k]
        dN[..., k] = 0.125 * g.prod(-1)
    return 0.125 * f.prod(-1), dN
def metric(C):          # det J (J^-1 J^-T)[:, dir], as ho_cvfem_metric_point
    G = np.empty((len(C), 3, P, 3, n, n))
    for d in range(3):
        rf = np.zeros((P, n, n, 3))
        rf[..., d] = xi[:, None, None]; rf[..., t1[d]] = zeta[None, :, None]
        rf[..., t2[d]] = zeta[None, None, :]
        J = np.einsum("eca,lsrcb->elsrab", C, shape(rf)[1])      # J[a][b] = dx_a / dxi_b
        Ji = np.linalg.inv(J)
        gv = np.linalg.det(J)[..., None] * np.einsum("...ak,...k->...a", Ji, Ji[..., d, :])
        G[:, d, :, 0] = gv[..., t2[d]]; G[:, d, :, 1] = gv[..., t1[d]]; G[:, d, :, 2] = gv[..., d]
    return G
def diagonal(G):        # the closed form of diagonal_kernel in mars_marsir_cellwise_solve.cu
    dg = np.zeros((len(G), n3))
    for node in range(n3):
        a, b, cc = node // nn, (node // n) % n, node % n
        for d in range(3):
            A = (a, b, cc)[d]; B = b if d == 0 else a; C = b if d == 2 else cc
            for l, sign in ((A - 1, 1.0), (A, -1.0)):
                if not 0 <= l < P:
                    continue
                g0, g1, g2 = G[:, d, l, 0], G[:, d, l, 1], G[:, d, l, 2]
                term = W[B, B] * W[C, C] * g2[:, B, C] * Dtil[l, A]
                term = term + W[B, B] * Btil[l, A] * (g0[:, B, :] @ (W[C] * D[:, C]))
                term = term + W[C, C] * Btil[l, A] * (g1[:, :, C] @ (W[B] * D[:, B]))
                dg[:, node] += sign * term
    return dg

K0 = np.stack([apply(np.eye(n3)[j][None], metric(1.0 * corner_ref[None]))[0]
               for j in range(n3)], 1)
print("CVFEM element operator: ||K - K^T|| / ||K|| =",
      "%.3e" % (np.linalg.norm(K0 - K0.T) / np.linalg.norm(K0)))

ne = args.ne; E = ne ** 3
lat = np.arange(ne + 1) / ne
V = np.stack(np.meshgrid(lat, lat, lat, indexing="ij"), -1)
V = V + args.deform * np.prod(np.sin(np.pi * V), -1, keepdims=True)
eidx = np.indices((ne, ne, ne)).reshape(3, -1).T            # e = (ex ne + ey) ne + ez
corners = np.stack([V[tuple((eidx + (s > 0)).T)] for s in corner_ref], 1)
G = metric(corners)
la, lb, lc = np.indices((n, n, n)).reshape(3, -1)
rf_nodes = np.stack([zeta[la], zeta[lb], zeta[lc]], -1)
xnodes = np.einsum("jc,eca->eja", shape(rf_nodes)[0], corners)
uex = np.prod(np.sin(np.pi * xnodes), -1)                    # at every element-local node

# DSS by the dimensionally split cascade: 3 axis passes of one-to-one face sums
def dss(fc):
    v = fc.reshape(ne, ne, ne, n, n, n).copy()
    for ax in range(3):                         # element axis ax <-> local axis ax
        lo = [slice(None)] * 6; hi = [slice(None)] * 6
        lo[ax] = slice(0, ne - 1); hi[ax] = slice(1, ne)
        lo[3 + ax] = n - 1; hi[3 + ax] = 0       # left element's last plane, right's first
        s = v[tuple(lo)] + v[tuple(hi)]
        v[tuple(lo)] = s; v[tuple(hi)] = s
    return v.reshape(E, n3)
NG = ne * p + 1
gidx = ((eidx[:, 0:1] * p + la) * NG + eidx[:, 1:2] * p + lb) * NG + eidx[:, 2:3] * p + lc
def gather(ug): return ug[gidx]
def scatter(fc):
    out = np.zeros(NG ** 3); np.add.at(out, gidx.ravel(), fc.ravel()); return out
rng = np.random.default_rng(1)
f = rng.standard_normal((E, n3))
dss_err = np.abs(dss(f) - gather(scatter(f))).max()
print("DSS cascade vs G G^T:", "%.1e" % dss_err)

dg_cell = diagonal(G)
Ec = min(E, 27)
probe = np.stack([apply(np.tile(np.eye(n3)[j], (Ec, 1)), G[:Ec])[:, j] for j in range(n3)], 1)
diag_err = np.abs(dg_cell[:Ec] - probe).max() / np.abs(probe).max()
print("closed-form diagonal vs probed:", "%.1e" % diag_err)

mult = dss(np.ones((E, n3)))                     # copies of each node
I3 = np.indices((NG, NG, NG)); bnd = ((I3 == 0) | (I3 == NG - 1)).any(axis=0).ravel()
mask_c = (~bnd)[gidx].astype(float)              # Dirichlet: zero boundary copies
def A_cell(uc): return apply(uc, G)              # local operator, primal -> dual
def wdot(u, v): return float(np.sum(u * v / mult))
diag_c = dss(dg_cell)                            # assembled diagonal, on every copy
def P_cell(rc): return mask_c * dss(rc) / diag_c # Jacobi-DSS: dual -> continuous primal
# assembled reference on the interior global DOFs
inner = np.where(~bnd)[0]
def A_glob(ug):
    full = np.zeros(NG ** 3); full[inner] = ug
    return scatter(A_cell(gather(full)))[inner]
dg = scatter(dg_cell)[inner]
def P_glob(rg): return rg / dg
b_cell = A_cell(uex)                             # unassembled RHS (dual)
b_glob = scatter(b_cell)[inner]
def bicgstab(A, P, b, dot, tol, iters, keep=False):
    x = np.zeros_like(b); r = P(b - A(x)); rh = r.copy(); r0 = np.sqrt(dot(r, r))
    rho0 = al = om = 1.0; v = np.zeros_like(b); pp = np.zeros_like(b); hist = []; xs = []
    for k in range(iters):
        rho = dot(rh, r); be = (rho / rho0) * (al / om); rho0 = rho
        pp = r + be * (pp - om * v); v = P(A(pp)); al = rho / dot(rh, v)
        s = r - al * v; t = P(A(s)); om = dot(t, s) / dot(t, t)
        x = x + al * pp + om * s; r = s - om * t
        hist.append(np.sqrt(dot(r, r)))
        if keep:
            xs.append(x.copy())
        if hist[-1] <= tol * r0:
            break
    return x, r0, hist, xs

ncmp = 40
_, _, _, xc = bicgstab(A_cell, P_cell, b_cell, wdot, 0.0, ncmp, keep=True)
_, _, _, xg = bicgstab(A_glob, P_glob, b_glob, lambda a, c: float(a @ c), 0.0, ncmp, keep=True)
def to_cell(xg):   # interior global vector -> element-local copies
    full = np.zeros(NG ** 3); full[inner] = xg
    return gather(full)
same_err = max(np.abs(a - to_cell(b)).max() for a, b in zip(xc, xg)) / np.abs(uex).max()
print(f"cell-wise vs assembled iterates, {ncmp} iterations: {same_err:.1e}")

x, r0, hist, _ = bicgstab(A_cell, P_cell, b_cell, wdot, args.tol, args.maxit)
err = np.sqrt(wdot(x - uex, x - uex) / wdot(uex, uex))
print(f"cell-wise BiCGStab, p={p}, {ne}^3 elements, {NG ** 3} unique DoFs, deform {args.deform:.3f}")
print(f"  iterations {len(hist)}, ||P r0||_w = {r0:.3e}, ||P r||_w = {hist[-1]:.3e}")
for k in (1, 5, 10, 20, 30, 40):
    if k <= len(hist):
        print(f"  it {k:3d}  ||P r||_w = {hist[k - 1]:.3e}")
print(f"  ||x - u||_w / ||u||_w = {err:.3e}")
ok = dss_err < 1e-12 and diag_err < 1e-12 and same_err < 1e-6 and err < 1e-8
print("CELL-WISE KRYLOV REFERENCE:", "PASS" if ok else "FAIL")
