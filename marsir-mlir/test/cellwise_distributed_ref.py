#!/usr/bin/env python3
"""Reference for the element-local solver on several MPI ranks
(backend/distributed/unstructured/solvers/mars_cellwise_layout.hpp and the multi-rank
paths of mars_cellwise_krylov.hpp and mars_cellwise_multigrid.hpp).

Each rank owns an equal sub-block of the global block; the ranks form a grid, rank
(rx * PY + ry) * PZ + rz, here the grids MPI_Dims_create returns (2: 2x1x1, 4: 2x2x1,
8: 2x2x2). Checks:
  * DSS with the ghost exchange: every rank receives from each of its up to 26
    neighbours the raw shared layer (a face plane, an edge line or a corner node per
    element) into a ghost shell, and gathers in the cascade order with global sharing.
    The result equals the one-domain DSS bit for bit. Ghost entries that no message
    filled are NaN, so reading anything the exchange did not deliver would show.
  * The multigrid solve on several ranks: each rank coarsens its own sub-block while all
    its dimensions are even; the level it cannot halve is gathered onto rank 0, which
    runs the rest of the V-cycle as one domain and scatters the correction back.
    BiCGStab reproduces the one-domain history up to the rounding of the per-rank
    partial sums.

Run: python3 test/cellwise_distributed_ref.py
"""
import itertools

import numpy as np

import cellwise_fem as fem
import cellwise_multigrid_ref as mgref
from cellwise_fem import n3

kN = fem.n
DIRS = [d for d in itertools.product((-1, 0, 1), repeat=3) if d != (0, 0, 0)]


class Decomp:
    def __init__(self, N, Pg):
        self.N, self.Pg = tuple(N), tuple(Pg)
        self.n = tuple(N[k] // Pg[k] for k in range(3))
        assert all(self.n[k] * Pg[k] == N[k] for k in range(3))
        self.ranks = list(itertools.product(*(range(p) for p in Pg)))   # rank order (rx, ry, rz)

    def offset(self, r):
        return tuple(r[k] * self.n[k] for k in range(3))

    def neighbour(self, r, d):
        q = tuple(r[k] + d[k] for k in range(3))
        return q if all(0 <= q[k] < self.Pg[k] for k in range(3)) else None

    def split(self, g):
        return {r: g[tuple(slice(o, o + m) for o, m in zip(self.offset(r), self.n))].copy() for r in self.ranks}

    def join(self, parts):
        g = np.empty(self.N + (kN, kN, kN))
        for r, v in parts.items():
            g[tuple(slice(o, o + m) for o, m in zip(self.offset(r), self.n))] = v
        return g


def cascade(v, N):
    """One-domain DSS of a global (NX, NY, NZ, 8, 8, 8) array."""
    v = v.copy()
    for ax in range(3):
        lo = [slice(None)] * 6; hi = [slice(None)] * 6
        lo[ax] = slice(0, N[ax] - 1); hi[ax] = slice(1, N[ax])
        lo[3 + ax] = kN - 1; hi[3 + ax] = 0
        s = v[tuple(lo)] + v[tuple(hi)]
        v[tuple(lo)] = s; v[tuple(hi)] = s
    return v


def layer(n, d, send):
    """Element and node slices of the shared layer toward direction d: what a rank sends
    to its neighbour at d (send), or where the neighbour's layer lands in the padded
    (n+2)^3 array (receive)."""
    es, ns = [], []
    for k in range(3):
        if d[k] > 0:
            es.append(slice(n[k] - 1, n[k]) if send else slice(n[k] + 1, n[k] + 2))
            ns.append(slice(kN - 1, kN) if send else slice(0, 1))
        elif d[k] < 0:
            es.append(slice(0, 1))
            ns.append(slice(0, 1) if send else slice(kN - 1, kN))
        else:
            es.append(slice(0, n[k]) if send else slice(1, n[k] + 1))
            ns.append(slice(0, kN))
    return tuple(es) + tuple(ns)


def exchange(dec, parts):
    """Every rank's padded array: own values inside, NaN shell except what arrived."""
    n = dec.n
    padded = {}
    for r in dec.ranks:
        p = np.full(tuple(m + 2 for m in n) + (kN, kN, kN), np.nan)
        p[1:n[0] + 1, 1:n[1] + 1, 1:n[2] + 1] = parts[r]
        padded[r] = p
    for r in dec.ranks:
        for d in DIRS:
            q = dec.neighbour(r, d)
            if q is not None:   # q sends its layer toward -d, which is r
                padded[r][layer(n, d, False)] = parts[q][layer(n, tuple(-x for x in d), True)]
    return padded


def share_side(a, g, Ng):
    return np.where(a == 0, np.where(g > 0, -1, 0), np.where(a == kN - 1, np.where(g < Ng - 1, 1, 0), 0))


def gather(dec, r, p):
    """Cascade-order gather on a rank's padded array, sharing from global coordinates."""
    n, o, N = dec.n, dec.offset(r), dec.N
    ex, ey, ez = (x.ravel() for x in np.indices(n))
    out = np.empty((ex.size, kN, kN, kN))
    gx, gy, gz = ex + o[0], ey + o[1], ez + o[2]
    for a, b, c in itertools.product(range(kN), repeat=3):
        sx, sy, sz = share_side(a, gx, N[0]), share_side(b, gy, N[1]), share_side(c, gz, N[2])
        lx, ly, lz = ex + 1 + np.minimum(sx, 0), ey + 1 + np.minimum(sy, 0), ez + 1 + np.minimum(sz, 0)
        a0, b0, c0 = np.where(sx != 0, kN - 1, a), np.where(sy != 0, kN - 1, b), np.where(sz != 0, kN - 1, c)
        zsum = None
        for k in (0, 1):
            ysum = None
            for j in (0, 1):
                bb, cc = np.where(j, 0, b0), np.where(k, 0, c0)
                xsum = p[lx, ly + j, lz + k, a0, bb, cc]
                xsum = np.where(sx != 0, xsum + p[lx + 1, ly + j, lz + k, 0, bb, cc], xsum)
                ysum = xsum if j == 0 else np.where(sy != 0, ysum + xsum, ysum)
            zsum = ysum if k == 0 else np.where(sz != 0, zsum + ysum, zsum)
        out[:, a, b, c] = zsum
    return out.reshape(n + (kN, kN, kN))


class DistLevel:
    """The ne^3 block on a rank grid. Vectors are (ranks * E_local, 512) arrays, rank by
    rank, each rank's elements in its local order."""
    def __init__(self, ne, Pg, deform):
        self.ne, self.dec = ne, Decomp((ne,) * 3, Pg)
        self.n = self.dec.n
        self.El = int(np.prod(self.n))
        V = fem.deform_vertices(ne, deform)
        loc = np.indices(self.n).reshape(3, -1).T
        g = np.concatenate([loc + np.array(self.dec.offset(r)) for r in self.dec.ranks])
        self.gidx = (g[:, 0] * ne + g[:, 1]) * ne + g[:, 2]
        corners = np.stack([V[tuple((g + (s > 0)).T)] for s in fem.corner_ref], 1)
        self.G = fem.metric(corners)
        la, lb, lc = np.indices((kN,) * 3).reshape(3, -1)
        outer = lambda e, l: ((e[:, None] == 0) & (l[None, :] == 0)) | ((e[:, None] == ne - 1) & (l[None, :] == kN - 1))
        self.mask = (~(outer(g[:, 0], la) | outer(g[:, 1], lb) | outer(g[:, 2], lc))).astype(float)
        sh = lambda e, l: ((l[None, :] == 0) & (e[:, None] > 0)) | ((l[None, :] == kN - 1) & (e[:, None] < ne - 1))
        self.w = 0.5 ** (sh(g[:, 0], la).astype(int) + sh(g[:, 1], lb) + sh(g[:, 2], lc))
        self.diag = self.dss(fem.diagonal(self.G))
        rf = np.stack([fem.zeta[la], fem.zeta[lb], fem.zeta[lc]], -1)
        self.x = np.einsum("jc,eca->eja", fem.shape(rf)[0], corners)

    def rows(self, i):
        return slice(i * self.El, (i + 1) * self.El)

    def dss(self, v):
        parts = {r: v[self.rows(i)].reshape(self.n + (kN,) * 3) for i, r in enumerate(self.dec.ranks)}
        padded = exchange(self.dec, parts)
        return np.concatenate([gather(self.dec, r, padded[r]).reshape(self.El, n3) for r in self.dec.ranks])

    def A(self, u): return fem.apply(u, self.G)
    def PJ(self, r): return self.mask * self.dss(r) / self.diag
    def wdot(self, u, v):   # per-rank partial sums, then the sum over ranks
        return float(sum(np.sum((u * v * self.w)[self.rows(i)]) for i in range(len(self.dec.ranks))))
    def to_global(self, v):
        out = np.empty_like(v); out[self.gidx] = v; return out
    def from_global(self, v): return v[self.gidx]


def halvable(n):
    return all(m % 2 == 0 for m in n)


class DistMultigrid:
    """The V-cycle of mars_cellwise_multigrid.hpp on a rank grid."""
    def __init__(self, ne, Pg, deform, pre, post, smoothing_range=15.0):
        self.levels = [DistLevel(ne, Pg, deform)]
        while halvable(self.levels[-1].n) and np.prod(self.levels[-1].n) > 1:
            self.levels.append(DistLevel(self.levels[-1].ne // 2, Pg, deform))
            if not (halvable(self.levels[-1].n) and np.prod(self.levels[-1].n) > 1):
                break
        self.pre, self.post, self.range = pre, post, smoothing_range
        self.lmax = [self.power(l) for l in self.levels[:-1]]
        self.root = mgref.Multigrid(self.levels[-1].ne, deform, pre, post)   # rank 0, one domain

    def power(self, lev, steps=30):   # hash of the global value index, as hash_kernel
        t = lev.gidx[:, None].astype(np.uint64) * np.uint64(n3) + np.arange(n3, dtype=np.uint64)[None, :]
        v = lev.PJ(((t * np.uint64(2654435761)) & np.uint64(0xffffffff)).astype(float) / 4294967296.0 - 0.5)
        for _ in range(steps):
            w = lev.PJ(lev.A(v)); nw = np.sqrt(lev.wdot(w, w))
            lam = nw / np.sqrt(lev.wdot(v, v)); v = w / nw
        return lam

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
        lev = self.levels[l]
        if l == len(self.levels) - 1:   # gathered onto rank 0, solved, scattered back
            return lev.from_global(self.root.vcycle(lev.to_global(b)))
        coarse = self.levels[l + 1]
        if self.pre:
            x = self.smooth(l, b, self.pre)
            r = lev.mask * (b - lev.A(x))
        else:
            r = lev.mask * b
        bc = np.concatenate([mgref.restrict(r[lev.rows(i)], coarse.n) for i in range(len(lev.dec.ranks))])
        xc = self.vcycle(bc, l + 1)
        e = np.concatenate([mgref.prolong(xc[coarse.rows(i)], coarse.n) for i in range(len(lev.dec.ranks))])
        x = x + e if self.pre else e
        return self.smooth(l, b, self.post, x) if self.post else x

    def __call__(self, r):
        return self.vcycle(r)


def main():
    rng = np.random.default_rng(11)
    ok = True
    print("DSS with the ghost exchange vs the one-domain DSS:")
    for N, Pg in [((4, 4, 4), (2, 2, 2)), ((8, 4, 2), (2, 2, 1)), ((4, 6, 2), (4, 3, 2)),
                  ((2, 2, 2), (2, 2, 2)), ((8, 8, 8), (2, 4, 1)), ((3, 3, 3), (3, 3, 3))]:
        dec = Decomp(N, Pg)
        g = rng.standard_normal(N + (kN, kN, kN))
        padded = exchange(dec, dec.split(g))
        res = dec.join({r: gather(dec, r, padded[r]) for r in dec.ranks})
        equal, nans = np.array_equal(res, cascade(g, N)), int(np.isnan(res).sum())
        ok &= equal and nans == 0
        print(f"  global {N}, ranks {Pg}, sub-block {dec.n}: bit-identical {equal}, unexchanged reads {nans}")

    print("multigrid BiCGStab on a rank grid vs one domain:")
    for ne, a, Pg in [(4, 0.05, (2, 1, 1)), (4, 0.05, (2, 2, 2)), (8, 0.1, (2, 2, 1))]:
        D = DistMultigrid(ne, Pg, a, 0, 3)
        lev = D.levels[0]
        u = np.prod(np.sin(np.pi * lev.x), -1)
        x, r0, h = mgref.bicgstab(lev, D, lev.A(u), 1e-10, 100)
        S = mgref.Multigrid(ne, a, 0, 3)
        one = S.levels[0]
        us = np.prod(np.sin(np.pi * one.x), -1)
        xs, r0s, hs = mgref.bicgstab(one, S, one.A(us), 1e-10, 100)
        dh = max(abs(p - q) / q for p, q in zip(h, hs)) if len(h) == len(hs) else np.inf
        dx = np.abs(lev.to_global(x) - xs).max()
        lam = D.lmax + D.root.lmax
        ok &= len(h) == len(hs) and dh < 1e-8 and dx < 1e-12
        print(f"  ne {ne}, deform {a}, ranks {Pg}: {len(h)} vs {len(hs)} iterations, "
              f"history rel diff {dh:.1e}, max |x - x_one| {dx:.1e}, "
              f"levels per rank {len(D.levels) - 1} then rank 0 ({D.levels[-1].ne}^3)")
        print("    lambda per level:", " ".join(f"{v:.4f}" for v in lam))
        print("    history:", " ".join(f"{v:.3e}" for v in h))
    print("CELL-WISE DISTRIBUTED REFERENCE:", "PASS" if ok else "FAIL")


if __name__ == "__main__":
    main()
