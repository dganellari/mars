#!/usr/bin/env python3
"""Reference for the geometric multigrid preconditioner on element-local vectors
(backend/distributed/unstructured/solvers/mars_cellwise_multigrid.hpp), after
Wichrowski, "Coalesced Matrix-Free Geometric Multigrid on Persistent Cell-Wise
Storage" (arXiv:2607.03413), here at p = 7 for the nonsymmetric CVFEM operator.

Levels: the ne^3 block (ne = 2^k) coarsened 2:1 down to one element, all at p = 7.
Prolongation evaluates the coarse polynomial at the child's GLL nodes; restriction
is its exact transpose on the raw unassembled residual. Smoother: Chebyshev on P_J A
(P_J = DSS + Jacobi), interval [lmax / 15, lmax], lmax = 1.1 x a 30-step power
iteration from a hash vector, --pre steps before and --post steps after the coarse
correction (pre = 0 restricts the right-hand side directly). The one-element level
is solved exactly. One V-cycle per preconditioner call of the left-preconditioned
BiCGStab.

Checks: prolongation keeps continuity, restriction is the transpose, degree-7
polynomials are reproduced, the V-cycle is linear with continuous output, and the
solve converges to the manufactured solution. The GPU run with the same --ne and
--deform, --pre and --post (and --mg) must reproduce the history to about 4 digits.

Run: python3 test/cellwise_multigrid_ref.py [--ne 4] [--deform 0.05] [--pre 0] [--post 3]
"""
import argparse

import numpy as np

import cellwise_fem as fem
from cellwise_fem import n, n3, zeta

ap = argparse.ArgumentParser()
ap.add_argument("--ne", type=int, default=4)
ap.add_argument("--deform", type=float, default=0.05)
ap.add_argument("--tol", type=float, default=1e-10)
ap.add_argument("--maxit", type=int, default=200)
ap.add_argument("--pre", type=int, default=0)
ap.add_argument("--post", type=int, default=3)
args = ap.parse_args()

I = [np.array([fem.lag(zeta, (z - 1) / 2)[0] for z in zeta]),
     np.array([fem.lag(zeta, (z + 1) / 2)[0] for z in zeta])]

def prolong(uc, ne_c):
    U = uc.reshape(ne_c, ne_c, ne_c, n, n, n); ne_f = 2 * ne_c
    out = np.empty((ne_f, ne_f, ne_f, n, n, n))
    for cx in (0, 1):
        for cy in (0, 1):
            for cz in (0, 1):
                out[cx::2, cy::2, cz::2] = np.einsum("aj,bk,cl,XYZjkl->XYZabc",
                                                     I[cx], I[cy], I[cz], U, optimize=True)
    return out.reshape(-1, n3)

def restrict(rf, ne_c):
    ne_f = 2 * ne_c; R = rf.reshape(ne_f, ne_f, ne_f, n, n, n)
    out = np.zeros((ne_c, ne_c, ne_c, n, n, n))
    for cx in (0, 1):
        for cy in (0, 1):
            for cz in (0, 1):
                out += np.einsum("aj,bk,cl,XYZabc->XYZjkl", I[cx], I[cy], I[cz],
                                 R[cx::2, cy::2, cz::2], optimize=True)
    return out.reshape(-1, n3)

def hash_vector(size):   # hash_kernel in mars_cellwise_multigrid.hpp
    t = np.arange(size, dtype=np.uint64)
    return ((t * np.uint64(2654435761)) & np.uint64(0xffffffff)).astype(float) / 4294967296.0 - 0.5

def power_iteration(lev, steps=30):
    v = lev.PJ(hash_vector(lev.E * n3).reshape(lev.E, n3))
    for _ in range(steps):
        w = lev.PJ(lev.A(v)); nw = np.sqrt(lev.wdot(w, w))
        lam = nw / np.sqrt(lev.wdot(v, v)); v = w / nw
    return lam

class Multigrid:
    def __init__(self, ne, deform, pre, post, smoothing_range=15.0):
        self.levels = [fem.Level(m, deform) for m in [ne >> k for k in range(ne.bit_length())]]
        assert self.levels[-1].E == 1 and ne & (ne - 1) == 0, "ne must be a power of two"
        self.pre, self.post, self.range = pre, post, smoothing_range
        self.lmax = [power_iteration(l) for l in self.levels[:-1]]
        c = self.levels[-1]
        self.inner = np.where(c.mask[0] > 0)[0]
        Aint = np.empty((len(self.inner),) * 2)
        for k, j in enumerate(self.inner):
            e = np.zeros((1, n3)); e[0, j] = 1.0
            Aint[:, k] = c.A(e)[0, self.inner]
        self.Ainv = np.linalg.inv(Aint)
        self.vcycles = 0

    def smooth(self, l, b, steps, x=None):
        lev = self.levels[l]
        lmax = 1.1 * self.lmax[l]; lmin = lmax / self.range
        theta, delta = 0.5 * (lmax + lmin), 0.5 * (lmax - lmin)
        sigma = theta / delta; rho = 1.0 / sigma
        d = lev.PJ(b if x is None else b - lev.A(x)) / theta
        x = d if x is None else x + d
        for _ in range(steps - 1):
            rho_new = 1.0 / (2.0 * sigma - rho)
            d = rho_new * rho * d + (2.0 * rho_new / delta) * lev.PJ(b - lev.A(x))
            x = x + d; rho = rho_new
        return x

    def vcycle(self, b, l=0):
        if l == len(self.levels) - 1:
            x = np.zeros((1, n3)); x[0, self.inner] = self.Ainv @ b[0, self.inner]
            return x
        lev, coarse = self.levels[l], self.levels[l + 1]
        if self.pre:
            x = self.smooth(l, b, self.pre)
            r = lev.mask * (b - lev.A(x))
        else:
            r = lev.mask * b
        e = prolong(self.vcycle(restrict(r, coarse.ne), l + 1), coarse.ne)
        x = x + e if self.pre else e
        return self.smooth(l, b, self.post, x) if self.post else x

    def __call__(self, r):
        self.vcycles += 1
        return self.vcycle(r)

def bicgstab(lev, P, b, tol, maxit):
    A, dot = lev.A, lev.wdot
    x = np.zeros_like(b); r = P(b); rh = r.copy(); r0 = np.sqrt(dot(r, r))
    rho0 = al = om = 1.0; v = np.zeros_like(b); pp = np.zeros_like(b); hist = []
    for _ in range(maxit):
        rho = dot(rh, r); be = (rho / rho0) * (al / om); rho0 = rho
        pp = r + be * (pp - om * v); v = P(A(pp)); al = rho / dot(rh, v)
        s = r - al * v; t = P(A(s)); om = dot(t, s) / dot(t, t)
        x = x + al * pp + om * s; r = s - om * t
        hist.append(np.sqrt(dot(r, r)))
        if hist[-1] <= tol * r0:
            break
    return x, r0, hist

rng = np.random.default_rng(3)
f, c = fem.Level(4, 0.05), fem.Level(2, 0.05)
uc = c.mask * c.dss(rng.standard_normal((c.E, n3))) / c.mult
uf = prolong(uc, 2)
cont = np.abs(f.dss(uf) / f.mult - uf).max()
rf = rng.standard_normal((f.E, n3))
transp = abs(np.sum(uf * rf) - np.sum(uc * restrict(rf, 2))) / abs(np.sum(uf * rf))
f0, c0 = fem.Level(4, 0.0), fem.Level(2, 0.0)
poly = lambda x: (x[..., 0] + 2 * x[..., 1] - x[..., 2]) ** 7 + x[..., 0] ** 3 * x[..., 2] ** 4
repro = np.abs(prolong(poly(c0.x), 2) - poly(f0.x)).max() / np.abs(poly(f0.x)).max()
print(f"prolongation continuity {cont:.1e}, restriction = transpose {transp:.1e}, "
      f"degree-7 reproduction {repro:.1e}")

M = Multigrid(args.ne, args.deform, args.pre, args.post)
lev = M.levels[0]
r1, r2 = rng.standard_normal((lev.E, n3)), rng.standard_normal((lev.E, n3))
z = M.vcycle(2.0 * r1 - 3.0 * r2)
lin = np.abs(z - (2.0 * M.vcycle(r1) - 3.0 * M.vcycle(r2))).max() / np.abs(z).max()
zc = np.abs(lev.dss(z) / lev.mult - z).max() / np.abs(z).max()
print(f"V-cycle linear {lin:.1e}, continuous {zc:.1e}")
print("lambda_max(P_J A) per level:", " ".join(f"{l:.4f}" for l in M.lmax))

u = np.prod(np.sin(np.pi * lev.x), -1); b = lev.A(u)
x, r0, hist = bicgstab(lev, M, b, args.tol, args.maxit)
err = np.sqrt(lev.wdot(x - u, x - u) / lev.wdot(u, u))
print(f"multigrid ({args.pre},{args.post}) BiCGStab, p=7, {args.ne}^3 elements, deform {args.deform:.3f}")
print(f"  iterations {len(hist)}, ||P r0||_w = {r0:.3e}, ||P r||_w = {hist[-1]:.3e}")
for k, h in enumerate(hist, 1):
    print(f"  it {k:3d}  ||P r||_w = {h:.3e}")
print(f"  ||x - u||_w / ||u||_w = {err:.3e}")
ok = cont < 1e-12 and transp < 1e-12 and repro < 1e-12 and lin < 1e-12 and zc < 1e-12 and err < 1e-8
print("CELL-WISE MULTIGRID REFERENCE:", "PASS" if ok else "FAIL")
