#!/usr/bin/env python3
"""Host study: which smoother keeps the BP3 p-multigrid iterations flat on stretched elements.

On the GPU (docs/design/cellwise_bp_fcg.md, "Stretched cube") point Chebyshev-Jacobi p-multigrid
took 14 / 91 / 401 / 648 iterations at --stretch 1 / 10 / 100 / 1000, in MARS and MFEM alike; with
deform 0.05 those elements are sheared as well as stretched. This runs the V-cycle of
marsir-mlir/test/cellwise_bp_mg_ref.py (p = 3 -> 1, Chebyshev (2,2) on [hi / 15, hi], an exact
p = 1 coarse solve, so only the smoother varies) on the cube with z divided by s after the
deformation, and changes what the Chebyshev polynomial is applied to:

  jacobi  : point Jacobi, DSS(r) / assembled diagonal (the GPU solver's smoother);
  element : additive Schwarz on element closures, each block the assembled operator restricted to
            the element's free nodes, solved exactly (dense here);
  column  : additive Schwarz on the columns of elements along z (the thin direction), solved
            exactly: a line smoother at element granularity (cube only);
  xcolumn : the same along x, across the thin direction (control, cube only);
  nodeline: block Jacobi on lines of nodes along each element line's thin axis (NodeLines).

The mesh code is generic (corner vertex ids, topological node numbering, Dirichlet on the faces
only one element has), so scripts/bp3_wing_block_study.py runs the same on blocks of real meshes.
Reported: FCG iterations to a true assembled residual of --rtol, with b = A u as the GPU driver.
Run: python3 scripts/bp3_smoother_study.py [--ne 8] [--deform 0.05] [--stretch 1,10,100,1000]
                                           [--smoothers jacobi,element,column,nodeline] [--line-cut k] [--grade g]
"""
import argparse
import time

import numpy as np
from numpy.polynomial import legendre as L

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


def faces_of(vid):   # {sorted corner ids of a face: [(element, axis, side), ...]}
    out = {}
    for e in range(len(vid)):
        for k in range(3):
            for t in (-1, 1):
                out.setdefault(tuple(sorted(vid[e, corner_ref[:, k] == t])), []).append((e, k, t))
    return out


def numbering(vid, p):
    """Global ids of the degree-p GLL nodes of hexahedra with corner vertex ids vid (E, 8, in
    corner_ref order): vertices by id, edge and face nodes by their vertices and a position in a
    canonical frame (so both elements of a face agree whatever their orientation), interior nodes
    per element. Returns gidx (E, (p + 1)^3) and free (N,), False on the faces only one element has."""
    n = p + 1; E = len(vid)
    corner = {tuple((corner_ref[c] > 0) * p): c for c in range(8)}
    boundary = {key for key, fl in faces_of(vid).items() if len(fl) == 1}
    ids = {}; gidx = np.empty((E, n ** 3), np.int64); fixed = set()
    for e in range(E):
        def V(x): return vid[e, corner[tuple(x)]]
        on_boundary = [(k, t) for k in range(3) for t in (-1, 1)
                       if tuple(sorted(vid[e, corner_ref[:, k] == t])) in boundary]
        for local, x in enumerate(np.ndindex(n, n, n)):
            ext = [i for i in range(3) if x[i] in (0, p)]
            if len(ext) == 3:
                key = ("v", V(x))
            elif len(ext) == 2:
                f = [i for i in range(3) if i not in ext][0]
                lo = list(x); lo[f] = 0; hi = list(x); hi[f] = p
                v0, v1 = V(lo), V(hi)
                key = ("e", min(v0, v1), max(v0, v1), x[f] if v0 < v1 else p - x[f])
            elif len(ext) == 1:
                k = ext[0]; k1, k2 = [i for i in range(3) if i != k]
                cs = {}
                for u in (0, p):
                    for w in (0, p):
                        y = list(x); y[k1] = u; y[k2] = w; cs[(u, w)] = V(y)
                ou, ow = min(cs, key=cs.get)
                u1 = p - x[k1] if ou == p else x[k1]; w1 = p - x[k2] if ow == p else x[k2]
                pos = (u1, w1) if cs[(p - ou, ow)] < cs[(ou, p - ow)] else (w1, u1)
                key = ("f", tuple(sorted(cs.values())), pos)
            else:
                key = ("i", e, local)
            g = ids.setdefault(key, len(ids)); gidx[e, local] = g
            if any(x[k] == (0 if t == -1 else p) for k, t in on_boundary):
                fixed.add(g)
    free = np.ones(len(ids), bool); free[list(fixed)] = False
    return gidx, free


class Level:
    """BP3 at degree p on the hexahedra `corners` (E, 8, 3) with node ids gidx, as
    cellwise_bp_mg_ref.Level; free is False on Dirichlet nodes."""

    def __init__(self, p, corners, gidx, free):
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
        self.gidx = gidx; self.free = free; self.N = len(free)
        self.inner = np.where(free)[0]
        self.mask = free[gidx].astype(float)
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

    def matrices(self):   # element matrices (E, n3, n3), probed
        return np.stack([self.apply(np.tile(np.eye(self.n3)[j], (self.E, 1))) for j in range(self.n3)], 2)

    def assembled(self):   # the assembled operator: sorted keys row * N + col and their summed values
        keys = (self.gidx[:, :, None] * self.N + self.gidx[:, None, :]).ravel()
        u, inv = np.unique(keys, return_inverse=True)
        return u, np.bincount(inv, weights=self.matrices().ravel())


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


def element_lines(corners, vid, cut=0):
    """Chains of elements along their thin axis (the local axis with the shortest mean edge): an
    element and its neighbour across a face are chained when the face lies across the thin axis of
    both. Each line is a list of (element, axis, entry side), longest first; cut > 0 splits lines
    into pieces of at most cut elements (as a rank boundary would)."""
    E = len(vid)
    length = np.zeros((E, 3))
    for k in range(3):
        pairs = [(i, j) for i in range(8) for j in range(i + 1, 8)
                 if np.sum(corner_ref[i] != corner_ref[j]) == 1 and corner_ref[i][k] != corner_ref[j][k]]
        length[:, k] = np.mean([np.linalg.norm(corners[:, i] - corners[:, j], axis=1) for i, j in pairs], 0)
    axis = np.argmin(length, 1)
    nxt = {}
    for fl in faces_of(vid).values():
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
    """p -> 1 with Chebyshev on prec(A) at p and an exact solve of the assembled p = 1 matrix."""

    def __init__(self, fine, coarse, prec, pre, post):
        self.f, self.c, self.prec, self.pre, self.post = fine, coarse, prec, pre, post
        self.I = np.array([lag(coarse.zeta, x)[0] for x in fine.zeta])
        self.lmax = self.power()
        A = np.zeros((coarse.N, coarse.N))
        np.add.at(A, (coarse.gidx[:, :, None], coarse.gidx[:, None, :]), coarse.matrices())
        self.Ac_inv = np.linalg.inv(A[np.ix_(coarse.inner, coarse.inner)])

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
        u = np.zeros(c.N); u[c.inner] = self.Ac_inv @ c.scatter(rc)[c.inner]
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


def levels(corners, vid, p):
    return Level(p, corners, *numbering(vid, p)), Level(1, corners, *numbering(vid, 1))


def solve(f, c, prec, u, pre, post, rtol, maxit):
    """FCG with the V-cycle on b = A u: iterations, max error relative to max |u|, lmax."""
    mg = VCycle(f, c, prec, pre, post)
    x, it = fcg(f, mg, f.apply(u), rtol, maxit)
    return it, np.abs(x - u).max() / np.abs(u).max(), mg.lmax


def smoother(name, f, corners, vid, cut=0, ne=0):
    if name == "jacobi":
        return f.jacobi, name
    if name == "element":
        return Schwarz(f, f.gidx), name
    if name in ("column", "xcolumn"):   # the cube's lattice order, e = (ex ne + ey) ne + ez
        g = f.gidx.reshape(ne, ne, ne, f.n3)
        cols = (g if name == "column" else g.transpose(1, 2, 0, 3)).reshape(ne * ne, ne * f.n3)
        return Schwarz(f, np.stack([np.unique(r) for r in cols])), name
    if name == "nodeline":
        prec = NodeLines(f, element_lines(corners, vid, cut))
        return prec, f"nodeline [{prec.stats}]"
    raise SystemExit(f"unknown smoother {name}")


def cube(ne, deform, s, grade=0.0):
    lat = np.arange(ne + 1) / ne
    V = np.stack(np.meshgrid(lat, lat, lat, indexing="ij"), -1)
    V = V + deform * np.prod(np.sin(np.pi * V), -1, keepdims=True)
    if grade:   # element layers growing by e^(grade / ne) in z: the stretch varies over the cube
        V[..., 2] = np.expm1(grade * V[..., 2]) / np.expm1(grade)
    V[..., 2] /= s
    eidx = np.indices((ne, ne, ne)).reshape(3, -1).T   # e = (ex ne + ey) ne + ez
    off = (corner_ref > 0).astype(int)
    corners = np.stack([V[tuple((eidx + o).T)] for o in off], 1)
    vid = ((eidx[:, None, 0] + off[:, 0]) * (ne + 1) + eidx[:, None, 1] + off[:, 1]) * (ne + 1) + eidx[:, None, 2] + off[:, 2]
    return corners, vid


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--ne", type=int, default=8)
    ap.add_argument("--p", type=int, default=3)
    ap.add_argument("--deform", type=float, default=0.05)
    ap.add_argument("--stretch", default="1,10,100,1000")
    ap.add_argument("--smoothers", default="jacobi,element,column,nodeline")
    ap.add_argument("--pre", type=int, default=2)
    ap.add_argument("--post", type=int, default=2)
    ap.add_argument("--rtol", type=float, default=1e-10)
    ap.add_argument("--maxit", type=int, default=3000)
    ap.add_argument("--line-cut", type=int, default=0, help="nodeline: split element lines every k elements")
    ap.add_argument("--grade", type=float, default=0.0, help="z layers grow by e^(grade / ne): thickest / thinnest = e^grade")
    args = ap.parse_args()
    ne = args.ne
    print(f"BP3 p = {args.p} -> 1, Chebyshev ({args.pre},{args.post}) on [hi / 15, hi], exact p = 1 solve, "
          f"{ne}^3 elements, deform {args.deform}, grade {args.grade}, FCG to true residual {args.rtol:.0e}")
    for s in [float(x) for x in args.stretch.split(",")]:
        corners, vid = cube(ne, args.deform, s, args.grade)
        f, c = levels(corners, vid, args.p)
        xs = f.xn.copy(); xs[..., 2] *= s
        if args.grade:
            xs[..., 2] = np.log1p(xs[..., 2] * np.expm1(args.grade)) / args.grade
        u = f.mask * np.prod(np.sin(np.pi * xs), -1)
        row = []
        for name in args.smoothers.split(","):
            t0 = time.time()
            prec, label = smoother(name, f, corners, vid, args.line_cut, ne)
            it, err, lmax = solve(f, c, prec, u, args.pre, args.post, args.rtol, args.maxit)
            row.append(f"{label} {it} it (err {err:.1e}, lmax {lmax:.2f}, {time.time() - t0:.0f} s)")
        print(f"  s = {s:g}: " + "; ".join(row), flush=True)


if __name__ == "__main__":
    main()
