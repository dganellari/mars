#!/usr/bin/env python3
"""ACM Stage 0 (host, NO-REGRET): directional aggregation + Galerkin-equivalence gate.

Proves the load-bearing algebra before any GPU code: that the paper's Eq.27-28 sum-of-fine
coarse operator (which the GPU builds by segmented-reduce over fine nonzeros grouped by
parent-aggregate) equals P^T A P bit-for-bit for a 0/1 block aggregate-membership P. If these
disagree, the GPU coarse-operator design (Stage 3) is wrong.

Builds the directional-P from the paper's criterion (Eqs 19-20, beta=0.5) with a GEOMETRIC
node-dual face weight |S_ij| = sum over triangles sharing edge (i,j) of |centroid - edge_mid|
(the median-dual SCS face length -- the legitimate node-centered analog of the paper's
cell-centered |face normal|, honestly labeled; NOT a bare FEM-edge proxy).

Reuses coupled_cavity_reference (skew_mesh) + phase1_precond_study (build = assembled coupled
[u,v,p] Stokes operator). Run: /tmp/p1env/bin/python scripts/acm_stage0_galerkin_gate.py
"""
import importlib.util
import numpy as np
import scipy.sparse as sp

spec = importlib.util.spec_from_file_location("p1", "scripts/phase1_precond_study.py")
p1 = importlib.util.module_from_spec(spec); spec.loader.exec_module(p1)
cav = p1.cav

BETA = 0.5
KMAX = 4         # paper fuses up to k_max-1 neighbors into a seed -> ~4 nodes/aggregate
NDOF = 3         # cavity is [u,v,p]


def dual_face_weights(P, tris):
    """omega[(i,j)] = median-dual SCS face length between nodes i,j (node-dual |S| analog)."""
    w = {}
    for t in tris:
        c = P[t].mean(axis=0)                      # triangle centroid
        for a, b in ((t[0], t[1]), (t[1], t[2]), (t[2], t[0])):
            mid = 0.5 * (P[a] + P[b])
            seg = np.linalg.norm(c - mid)          # SCS dual-face contribution from this triangle
            key = (min(a, b), max(a, b))
            w[key] = w.get(key, 0.0) + seg
    return w


def adjacency(w, N):
    nb = [[] for _ in range(N)]
    for (i, j) in w:
        nb[i].append(j); nb[j].append(i)
    return nb


def directional_aggregate(P, tris):
    """Faithful seed-growth (Eqs 19-20, beta=0.5, k_max), node-granular. Returns agg[node]."""
    N = len(P)
    w = dual_face_weights(P, tris)
    nb = adjacency(w, N)
    wij = lambda i, j: w[(min(i, j), max(i, j))]
    wmax = np.array([max((wij(i, j) for j in nb[i]), default=0.0) for i in range(N)])
    strong = lambda i, j: (wij(i, j) > BETA * wmax[i]) and (wij(i, j) > BETA * wmax[j])

    agg = np.full(N, -1, dtype=np.int64)
    naggr = 0
    for seed in range(N):
        if agg[seed] != -1:
            continue
        agg[seed] = naggr
        count = 1
        # fuse strong, non-aggregated neighbors of the seed up to k_max-1
        for j in nb[seed]:
            if count >= KMAX:
                break
            if agg[j] == -1 and strong(seed, j):
                agg[j] = naggr; count += 1
        naggr += 1

    # cleanup: leftover nodes join a neighboring aggregate (strong if possible, else any)
    for i in range(N):
        if agg[i] != -1:
            continue
        cand = [agg[j] for j in nb[i] if agg[j] != -1 and strong(i, j)]
        if not cand:
            cand = [agg[j] for j in nb[i] if agg[j] != -1]
        agg[i] = cand[0] if cand else naggr
        if not cand:
            naggr += 1
    return agg, naggr


def block_P(agg, naggr, ndof=NDOF):
    """0/1 block aggregate-membership: fine dof (ndof*node+c) -> coarse dof (ndof*I+c)."""
    N = len(agg)
    rows = np.arange(ndof * N)
    cols = np.array([ndof * agg[r // ndof] + (r % ndof) for r in rows])
    return sp.csr_matrix((np.ones(ndof * N), (rows, cols)), shape=(ndof * N, ndof * naggr))


def sumfine_coarse(A, agg, naggr, ndof=NDOF):
    """Eq.27-28 the GPU way: segmented reduce of fine nonzeros by (parent(row),parent(col))."""
    Ac = A.tocoo()
    dof2cdof = lambda d: ndof * agg[d // ndof] + (d % ndof)
    cr = np.array([dof2cdof(r) for r in Ac.row])
    cc = np.array([dof2cdof(c) for c in Ac.col])
    nC = ndof * naggr
    return sp.coo_matrix((Ac.data, (cr, cc)), shape=(nC, nC)).tocsr()


def run():
    for n in [16, 24, 32]:
        Pc, tris = cav.skew_mesh(n, 45.0)
        A, rhs = p1.build(n)                       # assembled coupled Stokes operator + RHS
        agg, naggr = directional_aggregate(Pc, tris)

        P = block_P(agg, naggr)
        A_ptap = (P.T @ A @ P).tocsr()
        A_sum = sumfine_coarse(A, agg, naggr)
        d_A = abs((A_ptap - A_sum)).max()

        r = rhs - A @ np.zeros(A.shape[0])         # residual at x=0 = rhs
        B_ptap = P.T @ r
        B_sum = np.zeros(NDOF * naggr)
        np.add.at(B_sum, np.array([NDOF * agg[d // NDOF] + (d % NDOF) for d in range(len(r))]), r)
        d_B = abs(B_ptap - B_sum).max()

        ratio = len(Pc) / naggr
        ok = (d_A < 1e-12) and (d_B < 1e-12)
        print(f"n={n:<3} nodes={len(Pc):<5} aggregates={naggr:<5} coarsening={ratio:4.2f}x  "
              f"|PtAP - sumfine|={d_A:.2e}  |Ptr - sumfine_r|={d_B:.2e}  -> {'PASS' if ok else 'FAIL'}")
    print("[gate] sum-of-fine (Eq.27-28, the GPU segmented-reduce recipe) == P^T A P for 0/1 block P.")
    print("[note] directional edge is NOT tested here (isotropic cavity) -- that is Stage 2 on a stretched mesh.")


if __name__ == "__main__":
    run()
