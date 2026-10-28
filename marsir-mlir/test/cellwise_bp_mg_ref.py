#!/usr/bin/env python3
"""Reference for the p-multigrid preconditioner of BP3 (the CEED Laplacian of
marsir-compiler/specs/bp3.op) inside the cell-wise flexible CG of
test/cellwise_bp_ref.py, before the GPU version (examples/.../mars_marsir_bp.cu --mg).

Levels p = 7 -> 3 -> 1 (or 3 -> 1) on the same mesh, the BP3 Galerkin form
rediscretized at every degree (q = p + 2 Gauss points, the geometry of the trilinear
map). One V-cycle per FCG preconditioner call, as PMultigrid::vcycle in
backend/distributed/unstructured/solvers/mars_cellwise_pmg.hpp:
  pre Chebyshev-Jacobi steps from x = 0 (DSS of the raw residual / assembled diagonal,
  zero on Dirichlet copies), r = b - A x, restrict r with the transposed GLL
  interpolation per element (raw in, raw out: no copy counted twice), the coarse
  level, x += interpolated coarse correction (continuous), post steps.
The smoother is the same polynomial (chebyshev_steps: [hi / 15, hi], hi = 1.1 lmax of
D^-1 A by 30 power steps) before and after, so the cycle is symmetric; on the p = 1
level an exact solve of the assembled p = 1 matrix stands in for BoomerAMG.

Checks:
  * restriction is the adjoint of prolongation: <r_f, P u_c> = <P^T r_f, u_c>;
  * a prolonged continuous vector is continuous (all copies of a node equal);
  * the fixed V-cycle, as an assembled matrix on a small mesh, is symmetric and
    positive definite (FCG needs an SPD preconditioner to keep its theory);
  * FCG iterations, Jacobi against the V-cycle, for growing ne: compared at the same
    TRUE assembled residual ||b - A x|| (each free node once) and solution error, not
    at the preconditioner-dependent sqrt(<r,z>/<r0,z0>).

--gpu: the run mars_marsir_bp.cu --mg --coarse k makes on its default (lexicographic,
unrotated) mesh, for its history to match: k Chebyshev steps on p = 1 instead of the
exact solve, the power iteration from the GPU's hash vector, b = A u, FCG stopped at
sqrt(<r,z>/<r0,z0>) <= --tol.

Run: python3 test/cellwise_bp_mg_ref.py [--levels 7,3,1] [--ne 2,4,8] [--deform 0.05]
     python3 test/cellwise_bp_mg_ref.py --gpu --levels 3,1 --ne 4 [--coarse 100] [--tol 1e-10]
"""
import argparse

import numpy as np
from numpy.polynomial import legendre as L

ap = argparse.ArgumentParser()
ap.add_argument("--levels", default="3,1")
ap.add_argument("--ne", default="2,4,8")
ap.add_argument("--deform", type=float, default=0.05)
ap.add_argument("--pre", type=int, default=2)
ap.add_argument("--post", type=int, default=2)
ap.add_argument("--rtol", type=float, default=1e-10, help="true relative assembled residual to stop at")
ap.add_argument("--maxit", type=int, default=3000)
ap.add_argument("--gpu", action="store_true", help="the GPU driver's run (see above)")
ap.add_argument("--coarse", type=int, default=100, help="--gpu: Chebyshev steps on p = 1")
ap.add_argument("--tol", type=float, default=1e-10, help="--gpu: FCG tolerance on sqrt(<r,z>/<r0,z0>)")
args = ap.parse_args()
degrees = [int(x) for x in args.levels.split(",")]

corner_ref = np.array([[-1, -1, -1], [1, -1, -1], [1, 1, -1], [-1, 1, -1],
                       [-1, -1, 1], [1, -1, 1], [1, 1, 1], [-1, 1, 1]])
SYM = ((0, 0), (0, 1), (0, 2), (1, 1), (1, 2), (2, 2))


def gll(p):
    c = np.zeros(p + 1); c[p] = 1
    return np.sort(np.concatenate(([-1.0], L.legroots(L.legder(c)), [1.0])))


def lag(nodes, x):
    m = len(nodes); v = np.ones(m); d = np.zeros(m)
    for j in range(m):
        others = [k for k in range(m) if k != j]
        den = np.prod([nodes[j] - nodes[k] for k in others])
        v[j] = np.prod([x - nodes[k] for k in others]) / den
        d[j] = sum(np.prod([x - nodes[k] for k in others if k != i]) for i in others) / den
    return v, d


def shape(rf):
    f = 1 + corner_ref * rf[..., None, :]
    dN = np.empty(f.shape)
    for k in range(3):
        g = f.copy(); g[..., k] = corner_ref[:, k]
        dN[..., k] = 0.125 * g.prod(-1)
    return 0.125 * f.prod(-1), dN


class Level:
    """BP3 at degree p on the ne^3 mesh: operator, diagonal, DSS, Dirichlet mask."""

    def __init__(self, p, ne, corners, eidx):
        self.p = p; n = self.n = p + 1; self.n3 = n ** 3; Q = self.Q = p + 2
        self.zeta = gll(p); gq, wq = L.leggauss(Q)
        self.B = np.array([lag(self.zeta, x)[0] for x in gq])
        self.D = np.array([lag(gq, x)[1] for x in gq])
        ga, gb, gc = np.indices((Q, Q, Q)).reshape(3, -1)
        rf = np.stack([gq[ga], gq[gb], gq[gc]], -1)
        J = np.einsum("eca,qcb->eqab", corners, shape(rf)[1])
        w = np.linalg.det(J) * (wq[ga] * wq[gb] * wq[gc])
        Ji = np.linalg.inv(J)
        M = w[..., None, None] * np.einsum("eqdr,eqfr->eqdf", Ji, Ji)
        self.G = np.stack([M[..., a, b] for a, b in SYM], 1)
        self.E = len(corners)
        la, lb, lc = np.indices((n, n, n)).reshape(3, -1)
        self.xn = np.einsum("jc,eca->eja", shape(np.stack([self.zeta[la], self.zeta[lb], self.zeta[lc]], -1))[0],
                            corners)
        NG = self.NG = ne * p + 1
        self.gidx = ((eidx[:, 0:1] * p + la) * NG + eidx[:, 1:2] * p + lb) * NG + eidx[:, 2:3] * p + lc
        I3 = np.indices((NG, NG, NG))
        bnd = ((I3 == 0) | (I3 == NG - 1)).any(axis=0).ravel()
        self.inner = np.where(~bnd)[0]
        self.mask = (~bnd)[self.gidx].astype(float)
        self.diag = self.dss(self.diagonal())

    def apply(self, u):
        n, Q, B, D = self.n, self.Q, self.B, self.D
        uq = np.einsum("ai,bj,ck,eijk->eabc", B, B, B, u.reshape(-1, n, n, n), optimize=True)
        d = [np.einsum("am,embc->eabc", D, uq), np.einsum("bm,eamc->eabc", D, uq), np.einsum("cm,eabm->eabc", D, uq)]
        g = self.G.reshape(self.E, 6, Q, Q, Q)
        sym = [[0, 1, 2], [1, 3, 4], [2, 4, 5]]
        dv = [sum(g[:, sym[r][s]] * d[s] for s in range(3)) for r in range(3)]
        w = (np.einsum("ma,embc->eabc", D, dv[0]) + np.einsum("mb,eamc->eabc", D, dv[1])
             + np.einsum("mc,eabm->eabc", D, dv[2]))
        return np.einsum("ai,bj,ck,eabc->eijk", B, B, B, w, optimize=True).reshape(-1, self.n3)

    def diagonal(self):
        n, B, D, Q = self.n, self.B, self.D, self.Q
        la, lb, lc = np.indices((n, n, n)).reshape(3, -1)
        ga, gb, gc = np.indices((Q, Q, Q)).reshape(3, -1)
        Bi, Bj, Bk = B[:, la], B[:, lb], B[:, lc]
        DB = D @ B
        grad = np.stack([DB[ga][:, la] * Bj[gb] * Bk[gc], Bi[ga] * DB[gb][:, lb] * Bk[gc],
                         Bi[ga] * Bj[gb] * DB[gc][:, lc]])
        out = np.zeros((self.E, self.n3))
        for c, (a, b) in enumerate(SYM):
            out += (1 if a == b else 2) * (self.G[:, c] @ (grad[a] * grad[b]))
        return out

    def scatter(self, fc):
        out = np.zeros(self.NG ** 3); np.add.at(out, self.gidx.ravel(), fc.ravel()); return out

    def dss(self, fc):
        return self.scatter(fc)[self.gidx]

    def nodes(self, u):   # a continuous vector's node values (each node once)
        out = np.zeros(self.NG ** 3); out[self.gidx.ravel()] = u.ravel(); return out

    def jacobi(self, r):   # z = Jacobi(DSS(r)), zero on Dirichlet copies
        return self.mask * self.dss(r) / self.diag

    def true_residual(self, b, x):   # ||b - A x|| assembled, each free node once
        return np.linalg.norm(self.scatter(b - self.apply(x))[self.inner])


def chebyshev_steps(lmax, steps):
    hi = 1.1 * lmax; lo = hi / 15.0
    theta, delta = 0.5 * (hi + lo), 0.5 * (hi - lo); sigma = theta / delta
    rho = 1.0 / sigma; out = [(0.0, 1.0 / theta)]
    for _ in range(1, steps):
        rho_new = 1.0 / (2.0 * sigma - rho)
        out.append((rho_new * rho, 2.0 * rho_new / delta)); rho = rho_new
    return out


def interpolation(fine, coarse):   # I[a][j] = coarse basis j at fine GLL node a
    return np.array([lag(coarse.zeta, x)[0] for x in fine.zeta])


class VCycle:
    def __init__(self, levels, pre, post, coarse_steps=0):
        self.levels = levels; self.pre = pre; self.post = post; self.coarse_steps = coarse_steps
        self.I = [interpolation(levels[k], levels[k + 1]) for k in range(len(levels) - 1)]
        self.lmax = [self.power(lv) for lv in levels]
        c = levels[-1]
        # The coarse p = 1 matrix, assembled from its element matrices (BoomerAMG's input).
        K = np.stack([c.apply(np.tile(np.eye(c.n3)[j], (c.E, 1))) for j in range(c.n3)], 2)   # (E, i, j)
        A = np.zeros((c.NG ** 3, c.NG ** 3))
        np.add.at(A, (c.gidx[:, :, None], c.gidx[:, None, :]), K)
        self.Ac = A[np.ix_(c.inner, c.inner)]

    def power(self, lv):   # 30 steps from the GPU's hash vector (pmg_detail::hash_kernel)
        t = np.arange(lv.E * lv.n3, dtype=np.uint64)
        h = ((t * np.uint64(2654435761)) & np.uint64(0xffffffff)).astype(float) / 4294967296.0 - 0.5
        v = lv.jacobi(h.reshape(lv.E, lv.n3))
        lam = 0.0
        for _ in range(30):
            w = lv.jacobi(lv.apply(v))
            nw = np.linalg.norm(lv.nodes(w)); nv = np.linalg.norm(lv.nodes(v))   # as the GPU's counted norm
            lam = nw / nv; v = w / nw
        return lam

    def prolong(self, k, uc):
        f, c, I = self.levels[k], self.levels[k + 1], self.I[k]
        return np.einsum("ai,bj,ck,eijk->eabc", I, I, I, uc.reshape(-1, c.n, c.n, c.n), optimize=True).reshape(-1, f.n3)

    def restrict(self, k, rf):
        f, c, I = self.levels[k], self.levels[k + 1], self.I[k]
        return np.einsum("ai,bj,ck,eabc->eijk", I, I, I, rf.reshape(-1, f.n, f.n, f.n), optimize=True).reshape(-1, c.n3)

    def smooth(self, k, b, x, steps, zero_x):
        lv = self.levels[k]; d = None
        for i, (c1, c2) in enumerate(chebyshev_steps(self.lmax[k], steps)):
            y = None if (i == 0 and zero_x) else lv.apply(x)
            z = lv.jacobi(b - (y if y is not None else 0.0))
            d = c2 * z if i == 0 else c1 * d + c2 * z
            x = d.copy() if (i == 0 and zero_x) else x + d
        return x

    def __call__(self, b, k=0):
        lv = self.levels[k]
        if k == len(self.levels) - 1:
            if self.coarse_steps:
                return self.smooth(k, b, np.zeros_like(b), self.coarse_steps, True)
            u = np.zeros(lv.NG ** 3)
            u[lv.inner] = np.linalg.solve(self.Ac, lv.scatter(b)[lv.inner])
            return u[lv.gidx]
        x = self.smooth(k, b, np.zeros_like(b), self.pre, True) if self.pre else np.zeros_like(b)
        r = (b - lv.apply(x)) * lv.mask
        xc = self(self.restrict(k, r), k + 1)
        x = x + self.prolong(k, xc)
        return self.smooth(k, b, x, self.post, False) if self.post else x


def fcg(lv, P, b, rtol, maxit):
    x = np.zeros_like(b); r = b.copy(); z = P(r); p_ = z.copy(); rz = np.sum(r * z)
    r0 = lv.true_residual(b, x)
    for it in range(1, maxit + 1):
        q = lv.apply(p_); qp = np.sum(q * p_); al = rz / qp
        x = x + al * p_; r = r - al * q
        z = P(r); rz = np.sum(r * z); be = -np.sum(q * z) / qp
        p_ = z + be * p_
        if lv.true_residual(b, x) <= rtol * r0:
            return x, it
    return x, maxit


def fcg_history(lv, P, b, tol, maxit):   # the driver's FCG: history sqrt(<r,z>/<r0,z0>)
    x = np.zeros_like(b); r = b.copy(); z = P(r); p_ = z.copy(); rz = np.sum(r * z); rz0 = rz; hist = []
    for _ in range(maxit):
        q = lv.apply(p_); qp = np.sum(q * p_); al = rz / qp
        x = x + al * p_; r = r - al * q
        z = P(r); rz = np.sum(r * z); be = -np.sum(q * z) / qp
        p_ = z + be * p_
        hist.append(np.sqrt(max(rz, 0.0) / rz0))
        if hist[-1] <= tol:
            break
    return x, hist


def mesh(ne):
    lat = np.arange(ne + 1) / ne
    V = np.stack(np.meshgrid(lat, lat, lat, indexing="ij"), -1)
    V = V + args.deform * np.prod(np.sin(np.pi * V), -1, keepdims=True)
    eidx = np.indices((ne, ne, ne)).reshape(3, -1).T
    return np.stack([V[tuple((eidx + (s > 0)).T)] for s in corner_ref], 1), eidx


if args.gpu:
    ne = int(args.ne.split(",")[0])
    corners, eidx = mesh(ne)
    levels = [Level(p, ne, corners, eidx) for p in degrees]
    mg = VCycle(levels, args.pre, args.post, args.coarse)
    f = levels[0]
    u = f.mask * np.prod(np.sin(np.pi * f.xn), -1)
    b = f.apply(u)
    x, hist = fcg_history(f, mg, b, args.tol, args.maxit)
    true = f.true_residual(b, x) / f.true_residual(b, np.zeros_like(b))
    print(f"GPU-comparable: bp3, p={degrees[0]}, levels {degrees}, Chebyshev ({args.pre},{args.post}), coarse "
          f"{args.coarse} Chebyshev steps, {ne}^3 elements, deform {args.deform}")
    print(f"  lmax {' '.join('%.4f' % l for l in mg.lmax)}")
    print(f"  FCG + V-cycle: {len(hist)} iterations, sqrt(<r,z>/<r0,z0>) = {hist[-1]:.3e}, true assembled residual "
          f"{true:.3e}")
    for k in (1, 5, 10, 20, 50, 100):
        if k <= len(hist):
            print(f"  it {k:3d}  {hist[k - 1]:.6e}")
    print(f"  max |x - u| / max |u| = {np.abs(x - u).max() / np.abs(u).max():.3e}")
    raise SystemExit(0)

print(f"BP3 p-multigrid {degrees}, Chebyshev ({args.pre},{args.post}), exact p=1 coarse solve, deform {args.deform}")
for ne in [int(x) for x in args.ne.split(",")]:
    corners, eidx = mesh(ne)
    levels = [Level(p, ne, corners, eidx) for p in degrees]
    mg = VCycle(levels, args.pre, args.post)
    f, c = levels[0], levels[1]
    rng = np.random.default_rng(ne)
    if ne == min(int(x) for x in args.ne.split(",")):
        # Adjointness of the transfers and continuity of a prolonged continuous vector.
        uc = c.mask * c.dss(rng.standard_normal((c.E, c.n3))) / c.dss(np.ones((c.E, c.n3)))
        rf = rng.standard_normal((f.E, f.n3))
        adj = abs(np.sum(rf * mg.prolong(0, uc)) - np.sum(mg.restrict(0, rf) * uc)) / abs(np.sum(rf * mg.prolong(0, uc)))
        pu = mg.prolong(0, uc); cont = np.abs(pu - f.dss(pu) / f.dss(np.ones_like(pu))).max()
        print(f"  transfers: adjoint to {adj:.1e}, prolonged continuous vector continuous to {cont:.1e}")
        # The fixed V-cycle as an assembled matrix: symmetric, positive definite.
        m = len(f.inner); Mmat = np.zeros((m, m))
        for j in range(m):
            e = np.zeros(f.NG ** 3); e[f.inner[j]] = 1.0
            raw = np.zeros((f.E, f.n3))   # a raw residual whose assembly is e: put it on one copy
            first = {}
            for idx, g in np.ndenumerate(f.gidx):
                if g == f.inner[j] and g not in first:
                    first[g] = idx
            raw[first[f.inner[j]]] = 1.0
            Mmat[:, j] = f.scatter(mg(raw) / f.dss(np.ones_like(raw)))[f.inner]
        sym = np.linalg.norm(Mmat - Mmat.T) / np.linalg.norm(Mmat)
        ev = np.linalg.eigvalsh(0.5 * (Mmat + Mmat.T))
        print(f"  V-cycle on ne={ne}: ||M - M^T|| / ||M|| = {sym:.1e}, eigenvalues of (M + M^T)/2 in "
              f"[{ev.min():.3e}, {ev.max():.3e}]")
    u = f.mask * np.prod(np.sin(np.pi * f.xn), -1)
    b = f.apply(u)
    xj, itj = fcg(f, f.jacobi, b, args.rtol, args.maxit)
    xm, itm = fcg(f, mg, b, args.rtol, args.maxit)
    ej = np.abs(xj - u).max() / np.abs(u).max(); em = np.abs(xm - u).max() / np.abs(u).max()
    print(f"  ne={ne:2d} ({f.NG ** 3} DoFs): lmax {['%.3f' % l for l in mg.lmax]}; to true residual {args.rtol:.0e}: "
          f"Jacobi {itj} it (err {ej:.1e}), V-cycle {itm} it (err {em:.1e})")
