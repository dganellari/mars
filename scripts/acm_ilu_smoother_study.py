#!/usr/bin/env python3
"""ACM smoother host study -- ILU(0) vs damped-Jacobi on a convection-dominated operator.

Validates (BEFORE any GPU code) that ILU(0) is the right smoother for the high-Re stall: on a
convection-dominated operator the multigrid SMOOTHING FACTOR mu (high-frequency spectral radius of the
smoother iteration matrix S = I - M^-1 A) must stay bounded for ILU(0) while damped-Jacobi degrades
toward 1. mu -> 1 is exactly what makes the production V-cycle non-contractive and stalls FlexGMRES at
the restart cap at Re=1000.

Operator: 2D convection-diffusion on the unit square, FIRST-ORDER UPWIND convection along c=(1,1),
5-point stencil, Dirichlet. Upwind (NOT central) is the M-matrix darwish actually solves; central
makes the ILU(0) factorization itself blow up at high Pe (a misleading test). Interior stencil with
h=1/n, df=nu/h^2, ax=ay=1/h:  diag = 4 df + ax + ay,  west=-(df+ax), south=-(df+ay), east=-df, north=-df.

mu = spectral radius of S restricted to the HIGH-frequency eigenspace of sym(A) (the upper half of the
spectrum -- the part the coarse grid cannot fix, so the smoother must). Jacobi uses the optimal damping
omega = 4/(3 rho(D^-1 sym(A))).  numpy only.
"""
import numpy as np
from numpy.linalg import eigvals, eigh, inv


def build_2d_upwind_cd(n, nu, cx=1.0, cy=1.0):
    h = 1.0 / n
    df = nu / (h * h); ax = cx / h; ay = cy / h
    N = n * n
    A = np.zeros((N, N))
    ix = lambda i, j: i * n + j
    for i in range(n):
        for j in range(n):
            p = ix(i, j)
            A[p, p] = 4 * df + ax + ay
            if i > 0:     A[p, ix(i - 1, j)] = -(df + ay)   # south (upwind from i-1, cy>0)
            if i < n - 1: A[p, ix(i + 1, j)] = -df          # north
            if j > 0:     A[p, ix(i, j - 1)] = -(df + ax)   # west (upwind from j-1, cx>0)
            if j < n - 1: A[p, ix(i, j + 1)] = -df          # east
    return A, ax * h / (2 * nu) if nu > 0 else np.inf       # A, element Peclet ~ |c|h/(2nu)


def ilu0(A):
    """Dense IKJ ILU(0): factor on A's nonzero pattern, drop all fill. Returns L (unit lower), U."""
    n = A.shape[0]
    M = A.copy()
    nz = A != 0.0
    for i in range(1, n):
        for k in range(i):
            if nz[i, k]:
                M[i, k] /= M[k, k]
                for j in range(k + 1, n):
                    if nz[i, j] and nz[k, j]:     # ILU(0) drop rule: update only existing entries
                        M[i, j] -= M[i, k] * M[k, j]
    L = np.tril(M, -1) + np.eye(n)
    U = np.triu(M)
    return L, U


def smoothing_factor(S, Asym):
    """High-frequency spectral radius of S: restrict to the upper half of sym(A)'s spectrum."""
    w, V = eigh(Asym)
    hi = w > np.median(w)
    Vhf = V[:, hi]
    Shf = Vhf.T @ S @ Vhf
    return float(np.max(np.abs(eigvals(Shf))))


def jacobi_mu(A, Asym, I):
    d = np.diag(A).copy()
    Dinv = np.diag(1.0 / d)
    rho = float(np.max(np.abs(eigvals(Dinv @ Asym))))
    omega = 4.0 / (3.0 * rho)                     # optimal damped-Jacobi smoother
    S = I - omega * (Dinv @ A)
    return smoothing_factor(S, Asym), omega


def ilu_mu(A, Asym, I):
    L, U = ilu0(A)
    Minv = inv(U) @ inv(L)
    S = I - Minv @ A
    return smoothing_factor(S, Asym)


def greedy_color(A):
    """First-fit greedy coloring of the matrix adjacency graph (what the GPU would do, node-wise)."""
    n = A.shape[0]
    nz = A != 0.0
    color = -np.ones(n, dtype=int)
    for i in range(n):
        used = {color[j] for j in range(n) if nz[i, j] and j != i and color[j] >= 0}
        c = 0
        while c in used:
            c += 1
        color[i] = c
    return color, int(color.max()) + 1


def multicolor_ilu_mu(A, Asym, I):
    """ILU(0) on the COLOR-permuted matrix -- the GPU-parallel ordering. Tests whether the reordering
    keeps the convective smoothing strength or collapses it toward Jacobi (Duff-Meurant)."""
    color, nc = greedy_color(A)
    perm = np.argsort(color, kind="stable")           # group rows by color
    Aperm = A[np.ix_(perm, perm)]
    L, U = ilu0(Aperm)
    Minv_p = inv(U) @ inv(L)
    Minv = np.zeros_like(A)
    Minv[np.ix_(perm, perm)] = Minv_p                 # M^-1 back in natural order
    S = I - Minv @ A
    return smoothing_factor(S, Asym), nc


def dilu(A):
    """DILU (AMGX's default smoother): off-diagonals straight from A, only the diagonal D is factored.
    M = (D+L) D^-1 (D+U),  D_i = A_ii - sum_{k<i} A_ik D_k^-1 A_ki (over k with both A_ik,A_ki nonzero)."""
    n = A.shape[0]
    d = np.diag(A).astype(float).copy()
    for i in range(n):
        s = 0.0
        for k in range(i):
            if A[i, k] != 0.0 and A[k, i] != 0.0:
                s += A[i, k] * A[k, i] / d[k]
        d[i] = A[i, i] - s
    D = np.diag(d)
    return (D + np.tril(A, -1)) @ np.diag(1.0 / d) @ (D + np.triu(A, 1))


def multicolor_dilu_mu(A, Asym, I):
    """DILU on the color-permuted matrix -- AMGX's actual smoother. Cheapest to build/apply."""
    color, nc = greedy_color(A)
    perm = np.argsort(color, kind="stable")
    Mp = dilu(A[np.ix_(perm, perm)])
    Minv = np.zeros_like(A)
    Minv[np.ix_(perm, perm)] = inv(Mp)
    return smoothing_factor(I - Minv @ A, Asym), nc


if __name__ == "__main__":
    n = 22
    I = np.eye(n * n)
    print(f"2D upwind convection-diffusion, n={n} (N={n*n}), c=(1,1). Smoothing factor mu (lower=better).")
    print("mu -> 1 means the smoother fails on high-frequency convective error (the Re=1000 stall).")
    print("nat = natural order (serial-strong), mc = multicolor (GPU-parallel). DILU mc = what AMGX runs.\n")
    print(f"{'nu':>9} {'Pe_e':>8} | {'Jacobi':>8} | {'ILU nat':>8} {'ILU mc':>8} | {'DILU mc':>8} {'colors':>7}")
    okall = True
    for nu in [0.5, 0.1, 0.02, 0.005, 0.001, 0.0002]:
        A, pee = build_2d_upwind_cd(n, nu)
        Asym = 0.5 * (A + A.T)
        muj, _ = jacobi_mu(A, Asym, I)
        mui = ilu_mu(A, Asym, I)
        mumc, nc = multicolor_ilu_mu(A, Asym, I)
        mudmc, _ = multicolor_dilu_mu(A, Asym, I)
        okall = okall and (mumc < 0.7 and mumc < muj)
        print(f"{nu:>9.4f} {pee:>8.2f} | {muj:>8.3f} | {mui:>8.3f} {mumc:>8.3f} | {mudmc:>8.3f} {nc:>7d}")
    print(f"\nGATE: multicolor ILU mu stays <0.7 and < Jacobi across the Pe sweep (reordering keeps a real")
    print(f"smoothing advantage over Jacobi -- the GPU-parallel ordering is worth building).")
    print(f"-> {'PASS (multicolor retains the convective advantage)' if okall else 'CHECK (multicolor collapses toward Jacobi)'}")
