#!/usr/bin/env python3
"""Reference for the CEED benchmark operators BP1 (mass) and BP3 (Laplacian) on
element-local (cell-wise) vectors, solved by the communication-free flexible CG of
Wichrowski, "Coalesced Matrix-Free Finite Elements in Cell-Wise Storage" (Algorithm 1,
internal-notes/FlexibleCG.pdf).

The operators are those of marsir-compiler/specs/bp1.op and bp3.op: the GLL basis
(n = p + 1 nodes per axis) at q = p + 2 Gauss points, B[a][i] = l_i(g_a), D[a][b] the
derivative of the Gauss-point Lagrange basis; geometry per Gauss point g0 = w det J
(BP1) and the symmetric w det J J^-1 J^-T, entries 00 01 02 11 12 22 (BP3).

FCG on cell-wise data: the iterate x, the search direction p and the preconditioned
residual z are continuous (primal) element-local vectors, the residual r and q = A p
unassembled (dual) ones, and every scalar pairs one of each, so it is a plain sum
over all copies (Lemma 1: <f, u> = <C^T f, u> for continuous u). The preconditioner
z = mask * DSS(r) / DSS(diag) is the only place copies are summed. This script checks:
  * BP1 applied to u = 1, summed, gives the domain's volume (1: the deformation keeps
    the cube's boundary); BP3 applied to u = 1 gives 0; the BP3 element matrix is
    symmetric (FCG needs it);
  * the closed-form diagonal (the GPU driver's setup kernel) equals the probed one;
  * cell-wise FCG and FCG on the assembled system give the same iterates;
  * the solve converges to the manufactured solution: BP3 -lap u = f, BP1 the L2
    projection M u = (g, phi), with u = g = sin(pi x) sin(pi y) sin(pi z).

Mesh and data are those of examples/distributed/unstructured/mars_marsir_bp.cu (the
unit cube as ne^3 trilinear hexahedra, deformed by x += a sin(pi x) sin(pi y) sin(pi z)
(1,1,1)); its residual history must match.

Run: python3 test/cellwise_bp_ref.py --spec bp3 [--p 3] [--ne 3] [--deform 0.05]
"""
import argparse

import numpy as np
from numpy.polynomial import legendre as L

ap = argparse.ArgumentParser()
ap.add_argument("--spec", choices=("bp1", "bp3"), default="bp3")
ap.add_argument("--p", type=int, default=3)
ap.add_argument("--q", type=int, default=0, help="Gauss points per axis (default p + 2)")
ap.add_argument("--ne", type=int, default=3)
ap.add_argument("--deform", type=float, default=0.05)
ap.add_argument("--tol", type=float, default=1e-10)
ap.add_argument("--maxit", type=int, default=1000)
args = ap.parse_args()

p = args.p; n = p + 1; n3 = n ** 3; Q = args.q or p + 2; q3 = Q ** 3
c = np.zeros(p + 1); c[p] = 1
zeta = np.sort(np.concatenate(([-1.0], L.legroots(L.legder(c)), [1.0])))   # GLL nodes
gq, wq = L.leggauss(Q)                                                     # Gauss points


def lag(nodes, x):
    m = len(nodes); v = np.ones(m); d = np.zeros(m)
    for j in range(m):
        others = [k for k in range(m) if k != j]
        den = np.prod([nodes[j] - nodes[k] for k in others])
        v[j] = np.prod([x - nodes[k] for k in others]) / den
        d[j] = sum(np.prod([x - nodes[k] for k in others if k != i]) for i in others) / den
    return v, d


B = np.array([lag(zeta, x)[0] for x in gq])     # q x n
D = np.array([lag(gq, x)[1] for x in gq])       # q x q
corner_ref = np.array([[-1, -1, -1], [1, -1, -1], [1, 1, -1], [-1, 1, -1],
                       [-1, -1, 1], [1, -1, 1], [1, 1, 1], [-1, 1, 1]])


def shape(rf):   # trilinear shape functions and reference gradients at points rf (..., 3)
    f = 1 + corner_ref * rf[..., None, :]
    dN = np.empty(f.shape)
    for k in range(3):
        g = f.copy(); g[..., k] = corner_ref[:, k]
        dN[..., k] = 0.125 * g.prod(-1)
    return 0.125 * f.prod(-1), dN


ne = args.ne; E = ne ** 3
lat = np.arange(ne + 1) / ne
V = np.stack(np.meshgrid(lat, lat, lat, indexing="ij"), -1)
V = V + args.deform * np.prod(np.sin(np.pi * V), -1, keepdims=True)
eidx = np.indices((ne, ne, ne)).reshape(3, -1).T            # e = (ex ne + ey) ne + ez
corners = np.stack([V[tuple((eidx + (s > 0)).T)] for s in corner_ref], 1)

ga, gb, gc = np.indices((Q, Q, Q)).reshape(3, -1)
rf_q = np.stack([gq[ga], gq[gb], gq[gc]], -1)
J = np.einsum("eca,qcb->eqab", corners, shape(rf_q)[1])     # J[a][b] = dx_a / dxi_b
detJw = np.linalg.det(J) * (wq[ga] * wq[gb] * wq[gc])
Ji = np.linalg.inv(J)
M = detJw[..., None, None] * np.einsum("eqdr,eqfr->eqdf", Ji, Ji)
SYM = ((0, 0), (0, 1), (0, 2), (1, 1), (1, 2), (2, 2))
G = detJw[:, None] if args.spec == "bp1" else np.stack([M[..., a, b] for a, b in SYM], 1)   # (E, K, q^3)
xq = np.einsum("qc,eca->eqa", shape(rf_q)[0], corners)      # physical Gauss points
la, lb, lc = np.indices((n, n, n)).reshape(3, -1)
xn = np.einsum("jc,eca->eja", shape(np.stack([zeta[la], zeta[lb], zeta[lc]], -1))[0], corners)


def interp(u):   # (E, n^3) -> (E, q, q, q)
    return np.einsum("ai,bj,ck,eijk->eabc", B, B, B, u.reshape(-1, n, n, n), optimize=True)


def integrate(w):   # (E, q, q, q) -> (E, n^3), the transpose of interp
    return np.einsum("ai,bj,ck,eabc->eijk", B, B, B, w, optimize=True).reshape(-1, n3)


def apply(u, G):
    uq = interp(u)
    if args.spec == "bp1":
        return integrate(G[:, 0].reshape(uq.shape) * uq)
    d = [np.einsum("am,embc->eabc", D, uq), np.einsum("bm,eamc->eabc", D, uq), np.einsum("cm,eabm->eabc", D, uq)]
    g = G.reshape(len(G), 6, Q, Q, Q)
    sym = [[0, 1, 2], [1, 3, 4], [2, 4, 5]]
    dv = [sum(g[:, sym[r][s]] * d[s] for s in range(3)) for r in range(3)]
    w = (np.einsum("ma,embc->eabc", D, dv[0]) + np.einsum("mb,eamc->eabc", D, dv[1])
         + np.einsum("mc,eabm->eabc", D, dv[2]))
    return integrate(w)


def diagonal(G):   # the closed form of the driver's setup kernel: entry (i,j,k) of A_e
    i, j, k = la, lb, lc
    Bi, Bj, Bk = B[:, i], B[:, j], B[:, k]                  # (q, n^3)
    if args.spec == "bp1":
        phi2 = (Bi[ga] * Bj[gb] * Bk[gc]) ** 2               # (q^3, n^3)
        return G[:, 0] @ phi2
    DB = D @ B                                               # derivative of l_i at the Gauss points
    grad = np.stack([DB[ga][:, i] * Bj[gb] * Bk[gc], Bi[ga] * DB[gb][:, j] * Bk[gc],
                     Bi[ga] * Bj[gb] * DB[gc][:, k]])        # (3, q^3, n^3)
    out = np.zeros((len(G), n3))
    for c, (a, b) in enumerate(SYM):
        out += (1 if a == b else 2) * (G[:, c] @ (grad[a] * grad[b]))
    return out


def dss(fc):   # the dimensionally split cascade: 3 axis passes of one-to-one face sums
    v = fc.reshape(ne, ne, ne, n, n, n).copy()
    for ax in range(3):
        lo = [slice(None)] * 6; hi = [slice(None)] * 6
        lo[ax] = slice(0, ne - 1); hi[ax] = slice(1, ne)
        lo[3 + ax] = n - 1; hi[3 + ax] = 0
        s = v[tuple(lo)] + v[tuple(hi)]
        v[tuple(lo)] = s; v[tuple(hi)] = s
    return v.reshape(E, n3)


NG = ne * p + 1
gidx = ((eidx[:, 0:1] * p + la) * NG + eidx[:, 1:2] * p + lb) * NG + eidx[:, 2:3] * p + lc


def scatter(fc):
    out = np.zeros(NG ** 3); np.add.at(out, gidx.ravel(), fc.ravel()); return out


# Checks of the operator itself.
ones = np.ones((E, n3))
if args.spec == "bp1":
    print(f"BP1 applied to 1, summed: {apply(ones, G).sum():.15f} (the volume: 1)")
else:
    print(f"BP3 applied to 1: max {np.abs(apply(ones, G)).max():.1e} (of max {np.abs(apply(xn[..., 0], G)).max():.1e} for u = x)")
K0 = np.stack([apply(np.eye(n3)[j][None], G[:1])[0] for j in range(n3)], 1)
print(f"element matrix: ||K - K^T|| / ||K|| = {np.linalg.norm(K0 - K0.T) / np.linalg.norm(K0):.1e}")
Ec = min(E, 8)
probe = np.stack([apply(np.tile(np.eye(n3)[j], (Ec, 1)), G[:Ec])[:, j] for j in range(n3)], 1)
dg_cell = diagonal(G)
print(f"closed-form diagonal vs probed: {np.abs(dg_cell[:Ec] - probe).max() / np.abs(probe).max():.1e}")

# Dirichlet (BP3): zero on the boundary copies; BP1 has none.
I3 = np.indices((NG, NG, NG))
bnd = ((I3 == 0) | (I3 == NG - 1)).any(axis=0).ravel() if args.spec == "bp3" else np.zeros(NG ** 3, bool)
mask = (~bnd)[gidx].astype(float)
diag_c = dss(dg_cell)
uex = np.prod(np.sin(np.pi * xn), -1)
fq = np.prod(np.sin(np.pi * xq), -1) * (3 * np.pi ** 2 if args.spec == "bp3" else 1.0)
b_cell = integrate((detJw * fq).reshape(E, Q, Q, Q))     # (f, phi) per element: dual


def A_cell(u): return apply(u, G)
def P_cell(r): return mask * dss(r) / diag_c
def pair(f, u): return float(np.sum(f * u))             # dual . primal, no weights


inner = np.where(~bnd)[0]
def gather(xg):
    full = np.zeros(NG ** 3); full[inner] = xg; return full[gidx]
def A_glob(xg): return scatter(A_cell(gather(xg)))[inner]
dg_glob = scatter(dg_cell)[inner]
def P_glob(r): return r / dg_glob
b_glob = scatter(b_cell)[inner]


def fcg(A, P, b, tol, iters, keep=False):   # Algorithm 1, flexible beta
    x = np.zeros_like(b); r = b - A(x); z = P(r); p_ = z.copy()
    rz = pair(r, z); rz0 = rz; hist = []; xs = []
    for _ in range(iters):
        q = A(p_); qp = pair(q, p_); al = rz / qp
        x = x + al * p_; r = r - al * q
        z = P(r); rz = pair(r, z); be = -pair(q, z) / qp
        p_ = z + be * p_
        hist.append(np.sqrt(max(rz, 0.0) / rz0))
        if keep:
            xs.append(x.copy())
        if hist[-1] <= tol:
            break
    return x, hist, xs


ncmp = 25
_, _, xc = fcg(A_cell, P_cell, b_cell, 0.0, ncmp, keep=True)
_, _, xg = fcg(A_glob, P_glob, b_glob, 0.0, ncmp, keep=True)
same = max(np.abs(a - gather(g)).max() for a, g in zip(xc, xg)) / np.abs(uex).max()
print(f"cell-wise vs assembled FCG iterates, {ncmp} iterations: {same:.1e}")

x, hist, _ = fcg(A_cell, P_cell, b_cell, args.tol, args.maxit)
err = np.abs(x - mask * uex).max() / np.abs(uex).max()
print(f"cell-wise FCG, {args.spec}, p={p}, q={Q}, {ne}^3 elements, {NG ** 3} unique DoFs, deform {args.deform:.3f}")
print(f"  iterations {len(hist)}, sqrt(<r,z>/<r0,z0>) = {hist[-1]:.3e}")
for k in (1, 5, 10, 20, 50, 100):
    if k <= len(hist):
        print(f"  it {k:3d}  {hist[k - 1]:.6e}")
print(f"  max |x - u| / max |u| at the nodes = {err:.3e}")

# The GPU driver's solve: b = A u with u = mask * sin sin sin at the nodes, so FCG returns
# u; its history must match these lines.
xd, hd, _ = fcg(A_cell, P_cell, A_cell(mask * uex), args.tol, args.maxit)
print(f"GPU-comparable (b = A u): iterations {len(hd)}, sqrt(<r,z>/<r0,z0>) = {hd[-1]:.3e}, "
      f"max |x - u| / max |u| = {np.abs(xd - mask * uex).max() / np.abs(uex).max():.3e}")
for k in (1, 5, 10, 20, 50, 100):
    if k <= len(hd):
        print(f"  it {k:3d}  {hd[k - 1]:.6e}")
