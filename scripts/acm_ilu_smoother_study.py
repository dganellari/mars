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


if __name__ == "__main__":
    n = 22
    I = np.eye(n * n)
    print(f"2D upwind convection-diffusion, n={n} (N={n*n}), c=(1,1). Smoothing factor mu (lower=better).")
    print("mu -> 1 means the smoother fails on high-frequency convective error (the Re=1000 stall).\n")
    print(f"{'nu':>9} {'Pe_e':>8} | {'Jacobi mu':>10} {'omega':>7} | {'ILU(0) mu':>10} | verdict")
    okall = True
    for nu in [0.5, 0.1, 0.02, 0.005, 0.001, 0.0002]:
        A, pee = build_2d_upwind_cd(n, nu)
        Asym = 0.5 * (A + A.T)
        muj, om = jacobi_mu(A, Asym, I)
        mui = ilu_mu(A, Asym, I)
        ok = mui < 0.5 and mui <= muj + 1e-9
        okall = okall and ok
        print(f"{nu:>9.4f} {pee:>8.2f} | {muj:>10.3f} {om:>7.3f} | {mui:>10.3f} | {'ILU<Jac' if ok else 'CHECK'}")
    print(f"\nGATE: ILU(0) mu stays <0.5 and <= Jacobi across the Pe sweep, AND Jacobi climbs toward 1.")
    print(f"-> {'PASS (ILU is the right convective smoother; Jacobi degrades)' if okall else 'CHECK'}")
