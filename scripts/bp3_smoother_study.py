#!/usr/bin/env python3
"""Host study: which smoother keeps the BP3 p-multigrid iterations flat on stretched elements.

On the GPU (docs/design/cellwise_bp_fcg.md, "Stretched cube") point Chebyshev-Jacobi p-multigrid
took 14 / 91 / 401 / 648 iterations at element aspect ratio 1 / 10 / 100 / 1000, in MARS and MFEM
alike. This runs the V-cycle of marsir-mlir/test/cellwise_bp_mg_ref.py (p = 3 -> 1, Chebyshev
(2,2) on [hi / 15, hi], an exact p = 1 coarse solve, so only the smoother varies) on the cube with
z divided by s after the deformation, and changes what the Chebyshev polynomial is applied to:

  jacobi : point Jacobi, DSS(r) / assembled diagonal (the GPU solver's smoother);
  element: additive Schwarz on element closures, each block the assembled operator restricted to
           the element's free nodes, solved exactly (dense here);
  column : additive Schwarz on the columns of elements along z (the stretched direction), solved
           exactly: a line smoother at element granularity;
  xcolumn: the same along x, across the stretch (control).

Reported: FCG iterations to a true assembled residual of --rtol, with b = A u as the GPU driver.
Run: python3 scripts/bp3_smoother_study.py [--ne 8] [--stretch 1,10,100,1000] [--smoothers jacobi,element,column]
"""
import argparse
import time

import numpy as np
from numpy.polynomial import legendre as L

ap = argparse.ArgumentParser()
ap.add_argument("--ne", type=int, default=8)
ap.add_argument("--p", type=int, default=3)
ap.add_argument("--deform", type=float, default=0.05)
ap.add_argument("--stretch", default="1,10,100,1000")
ap.add_argument("--smoothers", default="jacobi,element,column")
ap.add_argument("--pre", type=int, default=2)
ap.add_argument("--post", type=int, default=2)
ap.add_argument("--rtol", type=float, default=1e-10)
ap.add_argument("--maxit", type=int, default=3000)
ap.add_argument("--line-cut", type=int, default=0, help="nodeline: split element lines every k elements")
args = ap.parse_args()

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
    """BP3 at degree p on the ne^3 mesh, as cellwise_bp_mg_ref.Level."""

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
        NG = self.NG = ne * p + 1; self.N = NG ** 3
        self.gidx = ((eidx[:, 0:1] * p + la) * NG + eidx[:, 1:2] * p + lb) * NG + eidx[:, 2:3] * p + lc
        I3 = np.indices((NG, NG, NG))
        bnd = ((I3 == 0) | (I3 == NG - 1)).any(axis=0).ravel()
        self.free = ~bnd
        self.inner = np.where(self.free)[0]
        self.mask = self.free[self.gidx].astype(float)
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
        out = np.zeros(self.N); np.add.at(out, self.gidx.ravel(), fc.ravel()); return out

    def dss(self, fc):
        return self.scatter(fc)[self.gidx]

    def nodes(self, u):
        out = np.zeros(self.N); out[self.gidx.ravel()] = u.ravel(); return out

    def jacobi(self, r):
        return self.mask * self.dss(r) / self.diag

    def true_residual(self, b, x):
        return np.linalg.norm(self.scatter(b - self.apply(x))[self.inner])

    def assembled(self):   # the assembled operator: sorted keys row * N + col and their summed values
        K = np.stack([self.apply(np.tile(np.eye(self.n3)[j], (self.E, 1))) for j in range(self.n3)], 2)
        keys = (self.gidx[:, :, None] * self.N + self.gidx[:, None, :]).ravel()
        u, inv = np.unique(keys, return_inverse=True)
        return u, np.bincount(inv, weights=K.ravel())


class Schwarz:
    """z = sum over node sets S of R_S^T A_SS^-1 R_S DSS(r): A_SS the assembled operator on the free
    nodes of S, solved exactly; the result is continuous and zero on Dirichlet nodes."""

    def __init__(self, lv, sets):
        self.lv = lv; self.sets = sets
        keys, vals = lv.assembled()
        kk = sets[:, :, None] * lv.N + sets[:, None, :]
        pos = np.minimum(np.searchsorted(keys, kk), len(keys) - 1)
        A = np.where(keys[pos] == kk, vals[pos], 0.0)
        f = lv.free[sets].astype(float)
        A = A * f[:, :, None] * f[:, None, :] + (1.0 - f)[:, :, None] * np.eye(sets.shape[1])
        self.inv = np.linalg.inv(A); self.f = f

    def __call__(self, r):
        lv = self.lv
        z = np.einsum("sij,sj->si", self.inv, lv.scatter(r)[self.sets] * self.f) * self.f
        Z = np.zeros(lv.N); np.add.at(Z, self.sets.ravel(), z.ravel())
        return Z[lv.gidx]


def element_lines(lv, corners, eidx, ne, cut):
    """Chains of elements along their thin axis (the local axis with the shortest mean edge): e and
    its neighbour across a face are chained when the face lies across the thin axis of both. Each
    line is a list of (element, axis, entry side), longest first; cut > 0 splits lines into pieces
    of at most cut elements (as a rank boundary would)."""
    E = lv.E
    vid = ((eidx[:, None, 0] + (corner_ref[None, :, 0] > 0)) * (ne + 1) + eidx[:, None, 1]
           + (corner_ref[None, :, 1] > 0)) * (ne + 1) + eidx[:, None, 2] + (corner_ref[None, :, 2] > 0)
    length = np.zeros((E, 3))
    for k in range(3):
        pairs = [(i, j) for i in range(8) for j in range(i + 1, 8)
                 if np.sum(corner_ref[i] != corner_ref[j]) == 1 and corner_ref[i][k] != corner_ref[j][k]]
        length[:, k] = np.mean([np.linalg.norm(corners[:, i] - corners[:, j], axis=1) for i, j in pairs], 0)
    axis = np.argmin(length, 1)
    faces = {}
    for e in range(E):
        for k in range(3):
            for t in (-1, 1):
                faces.setdefault(tuple(sorted(vid[e, corner_ref[:, k] == t])), []).append((e, k, t))
    nxt = {}
    for fl in faces.values():
        if len(fl) == 2 and all(axis[e] == k for e, k, _ in fl):
            (e1, _, t1), (e2, _, t2) = fl
            nxt[(e1, t1)] = (e2, t2); nxt[(e2, t2)] = (e1, t1)
    seen = np.zeros(E, bool); out = []

    def walk(e, entry):
        line = []
        while not seen[e]:
            seen[e] = True; line.append((e, axis[e], entry))
            if (e, -entry) not in nxt:
                break
            e, entry = nxt[(e, -entry)]
        return line

    for e in range(E):
        for t in (-1, 1):
            if not seen[e] and (e, t) not in nxt:
                out.append(walk(e, t))
    for e in range(E):   # closed loops: open them anywhere
        if not seen[e]:
            out.append(walk(e, -1))
    if cut:
        out = [line[i:i + cut] for line in out for i in range(0, len(line), cut)]
    return sorted(out, key=len, reverse=True)


class NodeLines:
    """Block Jacobi on lines of nodes: every free node belongs to exactly one block, a run of nodes
    along the thin axis of its element line, solved exactly with the assembled operator restricted
    to it. Lines are claimed longest first; a node already claimed ends the run."""

    def __init__(self, lv, lines):
        n, p = lv.n, lv.p
        loc = [dict(zip(row, range(lv.n3))) for row in lv.gidx]
        abc = np.stack(np.indices((n, n, n)).reshape(3, -1), 1)
        claimed = np.zeros(lv.N, bool); blocks = []
        for line in lines:
            e0, k0, t0 = line[0]
            k1, k2 = [a for a in range(3) if a != k0]
            for i in range(n):
                for j in range(n):
                    c = np.zeros(3, int); c[k1] = i; c[k2] = j
                    chain = []
                    for m, (e, k, t) in enumerate(line):
                        if m > 0:   # continue at the node the previous element ended on
                            c = abc[loc[e][chain[-1]]].copy()
                            assert c[k] == (0 if t == -1 else p)
                        ks = range(n) if t == -1 else range(p, -1, -1)
                        for s, kk in enumerate(ks):
                            if m > 0 and s == 0:
                                continue
                            c[k] = kk; chain.append(lv.gidx[e, (c[0] * n + c[1]) * n + c[2]])
                    run = []
                    for g in chain:
                        if lv.free[g] and not claimed[g]:
                            claimed[g] = True; run.append(g)
                        elif run:
                            blocks.append(run); run = []
                    if run:
                        blocks.append(run)
        assert claimed[lv.inner].all()
        keys, vals = lv.assembled()
        self.lv = lv; self.groups = []
        for m in sorted({len(b) for b in blocks}):
            S = np.array([b for b in blocks if len(b) == m])
            kk = S[:, :, None] * lv.N + S[:, None, :]
            pos = np.minimum(np.searchsorted(keys, kk), len(keys) - 1)
            self.groups.append((S, np.linalg.inv(np.where(keys[pos] == kk, vals[pos], 0.0))))
        sizes = np.array([len(b) for b in blocks])
        self.stats = f"{len(lines)} element lines (longest {len(lines[0])}), {len(blocks)} node blocks, " \
                     f"mean {sizes.mean():.1f} nodes"

    def __call__(self, r):
        lv = self.lv; R = lv.scatter(r); Z = np.zeros(lv.N)
        for S, inv in self.groups:
            Z[S] = np.einsum("sij,sj->si", inv, R[S])
        return Z[lv.gidx]


def chebyshev_steps(lmax, steps):
    hi = 1.1 * lmax; lo = hi / 15.0
    theta, delta = 0.5 * (hi + lo), 0.5 * (hi - lo); sigma = theta / delta
    rho = 1.0 / sigma; out = [(0.0, 1.0 / theta)]
    for _ in range(1, steps):
        rho_new = 1.0 / (2.0 * sigma - rho)
        out.append((rho_new * rho, 2.0 * rho_new / delta)); rho = rho_new
    return out


class VCycle:
    """p = 3 -> 1 with Chebyshev on prec(A) at p = 3 and an exact p = 1 solve."""

    def __init__(self, fine, coarse, prec, pre, post):
        self.f, self.c, self.prec, self.pre, self.post = fine, coarse, prec, pre, post
        self.I = np.array([lag(coarse.zeta, x)[0] for x in fine.zeta])
        self.lmax = self.power()
        K = np.stack([coarse.apply(np.tile(np.eye(coarse.n3)[j], (coarse.E, 1))) for j in range(coarse.n3)], 2)
        A = np.zeros((coarse.N, coarse.N))
        np.add.at(A, (coarse.gidx[:, :, None], coarse.gidx[:, None, :]), K)
        self.Ac = A[np.ix_(coarse.inner, coarse.inner)]

    def power(self):
        lv = self.f
        t = np.arange(lv.E * lv.n3, dtype=np.uint64)
        h = ((t * np.uint64(2654435761)) & np.uint64(0xffffffff)).astype(float) / 4294967296.0 - 0.5
        v = self.prec(h.reshape(lv.E, lv.n3)); lam = 0.0
        for _ in range(30):
            w = self.prec(lv.apply(v))
            nw = np.linalg.norm(lv.nodes(w)); lam = nw / np.linalg.norm(lv.nodes(v)); v = w / nw
        return lam

    def smooth(self, b, x, steps, zero_x):
        d = None
        for i, (c1, c2) in enumerate(chebyshev_steps(self.lmax, steps)):
            z = self.prec(b if (i == 0 and zero_x) else b - self.f.apply(x))
            d = c2 * z if i == 0 else c1 * d + c2 * z
            x = d.copy() if (i == 0 and zero_x) else x + d
        return x

    def __call__(self, b):
        f, c, I = self.f, self.c, self.I
        x = self.smooth(b, np.zeros_like(b), self.pre, True)
        r = (b - f.apply(x)) * f.mask
        rc = np.einsum("ai,bj,ck,eabc->eijk", I, I, I, r.reshape(-1, f.n, f.n, f.n), optimize=True).reshape(-1, c.n3)
        u = np.zeros(c.N); u[c.inner] = np.linalg.solve(self.Ac, c.scatter(rc)[c.inner])
        xc = u[c.gidx]
        x = x + np.einsum("ai,bj,ck,eijk->eabc", I, I, I, xc.reshape(-1, c.n, c.n, c.n), optimize=True).reshape(-1, f.n3)
        return self.smooth(b, x, self.post, False)


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


def mesh(ne, s):
    lat = np.arange(ne + 1) / ne
    V = np.stack(np.meshgrid(lat, lat, lat, indexing="ij"), -1)
    V = V + args.deform * np.prod(np.sin(np.pi * V), -1, keepdims=True)
    V[..., 2] /= s
    eidx = np.indices((ne, ne, ne)).reshape(3, -1).T   # e = (ex ne + ey) ne + ez
    return np.stack([V[tuple((eidx + (c > 0)).T)] for c in corner_ref], 1), eidx


ne = args.ne
print(f"BP3 p = {args.p} -> 1, Chebyshev ({args.pre},{args.post}) on [hi / 15, hi], exact p = 1 solve, "
      f"{ne}^3 elements, deform {args.deform}, FCG to true residual {args.rtol:.0e}")
for s in [float(x) for x in args.stretch.split(",")]:
    corners, eidx = mesh(ne, s)
    f, c = Level(args.p, ne, corners, eidx), Level(1, ne, corners, eidx)
    xs = f.xn.copy(); xs[..., 2] *= s
    u = f.mask * np.prod(np.sin(np.pi * xs), -1)
    b = f.apply(u)
    row = []
    for name in args.smoothers.split(","):
        t0 = time.time()
        if name == "jacobi":
            prec = f.jacobi
        elif name == "element":
            prec = Schwarz(f, f.gidx)
        elif name == "column":
            cols = f.gidx.reshape(ne * ne, ne * f.n3)
            prec = Schwarz(f, np.stack([np.unique(g) for g in cols]))
        elif name == "xcolumn":   # control: the same lines, but along x, across the stretch
            cols = f.gidx.reshape(ne, ne * ne, f.n3).transpose(1, 0, 2).reshape(ne * ne, ne * f.n3)
            prec = Schwarz(f, np.stack([np.unique(g) for g in cols]))
        elif name == "nodeline":
            prec = NodeLines(f, element_lines(f, corners, eidx, ne, args.line_cut))
            name = f"nodeline [{prec.stats}]"
        else:
            raise SystemExit(f"unknown smoother {name}")
        mg = VCycle(f, c, prec, args.pre, args.post)
        x, it = fcg(f, mg, b, args.rtol, args.maxit)
        err = np.abs(x - u).max() / np.abs(u).max()
        row.append(f"{name} {it} it (err {err:.1e}, lmax {mg.lmax:.2f}, {time.time() - t0:.0f} s)")
    print(f"  s = {s:g}: " + "; ".join(row), flush=True)
