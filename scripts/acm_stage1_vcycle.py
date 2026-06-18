#!/usr/bin/env python3
"""ACM Stage 1 (host): full paper V-cycle as a preconditioner; faithfulness check.

Builds the multilevel hierarchy faithfully (directional aggregation -> sum-of-fine Galerkin
coarse [proven == P^T A P in Stage 0] -> order-zero injection / sum transfers), smooths with a
TRUE per-node block-Jacobi (the GPU-forced replacement for ILU(0)), and runs it as a GMRES
preconditioner. (The fixed-smoother V-cycle is a stationary linear operator, so plain
preconditioned GMRES is the correct host check; FlexGMRES is only needed on the GPU where the
cycle becomes non-stationary.)

GATE (this is a faithfulness check, NOT a win):
  - block-Jacobi V-cycle must REPRODUCE the known isotropic unsmoothed growth (~x1.7-2.0) --
    the replica is faithful, not magically flat;
  - point-Jacobi must be much worse / diverge (confirms block-Jacobi is mandatory on the saddle
    block where the pressure diagonal is ~0);
  - coarse per-node diagonal blocks must be invertible (singular coarse blocks = a known
    unsmoothed-aggregation failure mode).

Reuses Stage 0 (aggregation, block_P) + phase1_precond_study (build). Levels >=1 use the paper's
algebraic 'mutual coefficients' strength (block Frobenius norm) since no geometry survives.
Run: /tmp/p1env/bin/python scripts/acm_stage1_vcycle.py
"""
import importlib.util
import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spla

s0spec = importlib.util.spec_from_file_location("s0", "scripts/acm_stage0_galerkin_gate.py")
s0 = importlib.util.module_from_spec(s0spec); s0spec.loader.exec_module(s0)
p1, cav = s0.p1, s0.cav
NDOF, BETA, KMAX = s0.NDOF, s0.BETA, s0.KMAX
MAX_COARSE = 60          # direct solve below this many DOFs (Algorithm-A coarsest)


def algebraic_weights(A, ndof):
    """Level>=1 strength: omega_IJ = ||block_IJ||_F (paper's 'mutual coefficients' route)."""
    Ac = A.tocoo()
    nB = A.shape[0] // ndof
    acc = {}
    for r, c, v in zip(Ac.row, Ac.col, Ac.data):
        I, J = r // ndof, c // ndof
        if I == J:
            continue
        key = (min(I, J), max(I, J))
        acc[key] = acc.get(key, 0.0) + v * v
    return {k: np.sqrt(s) for k, s in acc.items()}


def aggregate(weights, N):
    """Faithful seed-growth on a generic strength graph (Eqs 19-20, beta=0.5, k_max)."""
    nb = s0.adjacency(weights, N)
    wij = lambda i, j: weights[(min(i, j), max(i, j))]
    wmax = np.array([max((wij(i, j) for j in nb[i]), default=0.0) for i in range(N)])
    strong = lambda i, j: (wij(i, j) > BETA * wmax[i]) and (wij(i, j) > BETA * wmax[j])
    agg = np.full(N, -1, np.int64); k = 0
    for seed in range(N):
        if agg[seed] != -1:
            continue
        agg[seed] = k; cnt = 1
        for j in nb[seed]:
            if cnt >= KMAX:
                break
            if agg[j] == -1 and strong(seed, j):
                agg[j] = k; cnt += 1
        k += 1
    for i in range(N):
        if agg[i] != -1:
            continue
        cand = [agg[j] for j in nb[i] if agg[j] != -1 and strong(i, j)] or \
               [agg[j] for j in nb[i] if agg[j] != -1]
        agg[i] = cand[0] if cand else k
        if not cand:
            k += 1
    return agg, k


def build_hierarchy(A0, Pcoord, tris):
    levels = [A0]; Ps = []
    A = A0; lvl = 0
    while A.shape[0] > MAX_COARSE:
        if lvl == 0:
            agg, na = s0.directional_aggregate(Pcoord, tris)
        else:
            N = A.shape[0] // NDOF
            w = algebraic_weights(A, NDOF)
            if not w:
                break
            agg, na = aggregate(w, N)
        if na <= 1 or na * NDOF >= A.shape[0]:      # no real coarsening -> stop
            break
        P = s0.block_P(agg, na, NDOF)
        A = (P.T @ A @ P).tocsr()
        Ps.append(P); levels.append(A); lvl += 1
    return levels, Ps


def block_dinv(A, ndof, reg=1e-10):
    """Per-node ndof x ndof diagonal-block inverse (regularized for the ~0 pressure diagonal)."""
    n = A.shape[0] // ndof
    D = A.toarray() if A.shape[0] < 4000 else None
    blocks = np.empty((n, ndof, ndof))
    Ad = A.tocsr()
    worst = np.inf
    for i in range(n):
        sl = slice(ndof * i, ndof * i + ndof)
        B = (D[sl, sl] if D is not None else Ad[sl, sl].toarray())
        worst = min(worst, abs(np.linalg.det(B)))
        blocks[i] = np.linalg.inv(B + reg * np.eye(ndof))
    return blocks, worst


def smooth_block(A, b, x, Dinv, ndof, sweeps, omega=0.7):
    n = A.shape[0] // ndof
    for _ in range(sweeps):
        r = (b - A @ x).reshape(n, ndof)
        x = x + omega * np.einsum('nij,nj->ni', Dinv, r).reshape(-1)
    return x


def smooth_point(A, b, x, dinv, sweeps, omega=0.7):
    for _ in range(sweeps):
        x = x + omega * dinv * (b - A @ x)
    return x


def make_vcycle(levels, Ps, mode):
    facs = {}                                       # cached coarsest LU + smoother data
    Acoarse = levels[-1]
    facs['lu'] = spla.splu(Acoarse.tocsc())
    sm = []
    worst_det = np.inf
    for A in levels[:-1]:
        if mode == 'block':
            Dinv, wd = block_dinv(A, NDOF); worst_det = min(worst_det, wd)
            sm.append(('block', Dinv))
        else:
            d = A.diagonal(); d[np.abs(d) < 1e-30] = 1e-30
            sm.append(('point', 1.0 / d))
    def vcycle(lvl, b):
        if lvl == len(levels) - 1:
            return facs['lu'].solve(b)
        A = levels[lvl]; kind, data = sm[lvl]
        sweep = (smooth_block if kind == 'block' else smooth_point)
        args = (data, NDOF) if kind == 'block' else (data,)
        x = sweep(A, b, np.zeros_like(b), *args, 2)         # 2 pre-sweeps (paper)
        r = b - A @ x
        ec = vcycle(lvl + 1, Ps[lvl].T @ r)                 # restrict (sum) -> recurse
        x = x + Ps[lvl] @ ec                                # prolong (injection) correct
        x = sweep(A, b, x, *args, 1)                        # 1 post-sweep (paper)
        return x
    M = spla.LinearOperator(levels[0].shape, matvec=lambda b: vcycle(0, b))
    return M, worst_det


def run():
    print(f"{'n':<4}{'DOFs':<7}{'lvls':<6}{'blockJac':<10}{'pointJac':<10}{'worst|det|':<12}")
    for n in [16, 24, 32]:
        Pc, tris = cav.skew_mesh(n, 45.0)
        A, b = p1.build(n)
        levels, Ps = build_hierarchy(A, Pc, tris)
        Mb, wd = make_vcycle(levels, Ps, 'block')
        Mp, _ = make_vcycle(levels, Ps, 'point')
        itb = p1.gmres_iters(A, b, M=Mb)
        itp = p1.gmres_iters(A, b, M=Mp)
        print(f"{n:<4}{A.shape[0]:<7}{len(levels):<6}{itb:<10}{itp:<10}{wd:<12.2e}")
    print(f"[gate] block-Jacobi V-cycle should converge + reproduce unsmoothed growth (~x1.7-2.0,")
    print(f"       faithful, NOT magically flat); point-Jacobi should be much worse (saddle block);")
    print(f"       worst|det| > 0 confirms coarse diagonal blocks are invertible. vs 1a: 13/19/39.")


if __name__ == "__main__":
    run()
