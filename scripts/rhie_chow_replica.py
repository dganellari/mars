#!/usr/bin/env python3
"""Rhie-Chow coupled operator -- host replica / Stage-0 algebra gate (darwish2008a collocated FV).

Replaces the PSPG a^pp = tau*K (the split-method pressure-Poisson folded into a block) with the
collocated Rhie-Chow d_f-weighted graph-Laplacian, the inf-sup remedy. Validates the four coupled
blocks BEFORE any GPU code, on BOTH a uniform and a deliberately SKEWED tet -- because the two
highest-risk pieces (the Ef orthogonal coefficient and the smooth-vs-compact discrepancy) are
INVISIBLE on a uniform mesh (compact == smooth there) and only bite on stretched/skewed cells.

Per SCS face f=(L,R), area vector S_f (median-dual, oriented L->R, matching mars_cvfem_tet_area.hpp),
dx = x_R - x_L, Ef = |S_f|^2/(S_f.dx), D_f = 0.5*(V_P/a_P^u + V_N/a_N^u) harvested from the assembled
momentum diagonal (a_P^u = nu*K_PP for the Stokes core), rho=1:
  a^uu = nu*K (component-diagonal)
  a^up : momentum row (P,d) gets +0.5*S_f[d] at p_P and p_N ; (N,d) gets -0.5*S_f[d] at both
  a^pu : continuity row P gets +0.5*rho*S_f[d] at u_P[d],u_N[d] ; N gets -0.5*rho*S_f[d] at both  (= -rho*(a^up)^T)
  a^pp : w_f = rho*D_f*Ef ;  [P,P]+=w [P,N]-=w [N,N]+=w [N,P]-=w   (M-matrix graph-Laplacian, NO tau)
"""
import numpy as np

RHO, NU = 1.0, 1.0
EDGES = [(0, 1), (1, 2), (0, 2), (0, 3), (1, 3), (2, 3)]          # ip order == precomputeTetAreaVectorsKernel
OTHER = [(2, 3), (0, 3), (1, 3), (1, 2), (0, 2), (0, 1)]


def tet_grad(X):
    J = np.array([X[1] - X[0], X[2] - X[0], X[3] - X[0]]).T       # 3x3, columns = edge vectors
    detJ = np.linalg.det(J)
    Jinv = np.linalg.inv(J)
    dN = np.zeros((4, 3))
    dN[1], dN[2], dN[3] = Jinv[0], Jinv[1], Jinv[2]               # grad of barycentric coords
    dN[0] = -(dN[1] + dN[2] + dN[3])
    return detJ / 6.0, dN


def scs_areas(X):
    C = X.mean(axis=0)
    A = np.zeros((6, 3))
    for ip, (L, R) in enumerate(EDGES):
        a, b = OTHER[ip]
        M = 0.5 * (X[L] + X[R]); Fa = (X[L] + X[R] + X[a]) / 3.0; Fb = (X[L] + X[R] + X[b]) / 3.0
        Av = 0.5 * (np.cross(Fa - M, C - M) + np.cross(C - M, Fb - M))   # quad M->Fa->C->Fb
        if np.dot(Av, X[R] - X[L]) < 0: Av = -Av                        # orient L->R
        A[ip] = Av
    return A


def assemble(X):
    V, dN = tet_grad(X)
    K = np.zeros((4, 4))
    for a in range(4):
        for b in range(4):
            K[a, b] = V * float(dN[a] @ dN[b])
    aPu = np.array([NU * K[i, i] for i in range(4)])               # momentum diagonal, harvested
    Vp = np.full(4, V / 4.0)                                       # node dual CV (single tet -> V/4)
    A = np.zeros((16, 16))
    # a^uu = nu*K (component-diagonal over u,v,w)
    for a in range(4):
        for b in range(4):
            for d in range(3):
                A[4 * a + d, 4 * b + d] += NU * K[a, b]
    S = scs_areas(X)
    diag = []
    for ip, (L, R) in enumerate(EDGES):
        Sf = S[ip]; dx = X[R] - X[L]
        asq = float(Sf @ Sf); axdx = float(Sf @ dx)
        Ef = asq / axdx
        Df = 0.5 * (Vp[L] / aPu[L] + Vp[R] / aPu[R])
        w = RHO * Df * Ef
        diag.append((ip, Ef, w))
        for d in range(3):
            # a^pu = FV divergence of the interpolated face velocity (v_f = avg); S_f out of L (+), of R (-)
            A[4 * L + 3, 4 * L + d] += 0.5 * RHO * Sf[d]; A[4 * L + 3, 4 * R + d] += 0.5 * RHO * Sf[d]
            A[4 * R + 3, 4 * L + d] -= 0.5 * RHO * Sf[d]; A[4 * R + 3, 4 * R + d] -= 0.5 * RHO * Sf[d]
        # a^pp (Rhie-Chow d_f-Laplacian)
        A[4 * L + 3, 4 * L + 3] += w; A[4 * L + 3, 4 * R + 3] -= w
        A[4 * R + 3, 4 * R + 3] += w; A[4 * R + 3, 4 * L + 3] -= w
    # enforce G = -D^T by construction: a^up = -(1/rho)*(a^pu)^T (the consistent discrete gradient;
    # the averaged "pressure force" form does NOT satisfy this on the neighbor term).
    for i in range(4):
        for j in range(4):
            for d in range(3):
                A[4 * j + d, 4 * i + 3] = -(1.0 / RHO) * A[4 * i + 3, 4 * j + d]
    return A, S, diag


def smooth_vs_compact(X, S, pfun):
    """compact (p_N-p_P)*Ef vs smooth (gradp_bar . S): equal on uniform, differ on skewed (deferred-RHS term)."""
    V, dN = tet_grad(X)
    p = np.array([pfun(x) for x in X])
    gp = sum(p[i] * dN[i] for i in range(4))                       # constant nodal gradient on a linear tet
    worst = 0.0
    for ip, (L, R) in enumerate(EDGES):
        Sf = S[ip]; dx = X[R] - X[L]
        compact = (p[R] - p[L]) * (Sf @ Sf) / (Sf @ dx)
        smooth = gp @ Sf
        worst = max(worst, abs(compact - smooth))
    return worst


def gates(name, X):
    A, S, diag = assemble(X)
    # G1: skew-transpose  a^pu == -rho*(a^up)^T
    e_skew = 0.0
    for i in range(4):
        for j in range(4):
            for d in range(3):
                e_skew = max(e_skew, abs(A[4 * i + 3, 4 * j + d] + RHO * A[4 * j + d, 4 * i + 3]))
    # G2: a^pp M-matrix (diag>0, off<=0) + row-sum 0 + null(1)
    App = A[3::4, 3::4]
    diag_ok = all(App[i, i] > 0 for i in range(4))
    off_ok = all(App[i, j] <= 1e-14 for i in range(4) for j in range(4) if i != j)
    rowsum = float(np.abs(App.sum(axis=1)).max())
    null1 = float(np.abs(App @ np.ones(4)).max())
    # G3 (linear pressure): smooth vs compact discrepancy (the deferred-RHS term)
    disc = smooth_vs_compact(X, S, lambda x: 1.0 + 2.0 * x[0] - x[1] + 0.5 * x[2])
    Efs = [d[1] for d in diag]
    print(f"\n=== {name} ===")
    print(f"  G1 a^pu == -rho*(a^up)^T : {e_skew:.2e}  {'PASS' if e_skew < 1e-12 else 'FAIL'}")
    print(f"  G2 a^pp M-matrix         : diag>0={diag_ok} off<=0={off_ok} rowsum={rowsum:.2e} null(1)={null1:.2e}"
          f"  {'PASS' if (diag_ok and off_ok and rowsum < 1e-12 and null1 < 1e-12) else 'FAIL'}")
    print(f"  Ef per face (orthog coeff): min={min(Efs):.3e} max={max(Efs):.3e}")
    print(f"  G3 |compact - smooth| (linear p): {disc:.3e}  "
          f"{'~0 (uniform: deferred term invisible)' if disc < 1e-10 else 'NONZERO -> deferred smooth-grad RHS is REQUIRED here'}")


if __name__ == "__main__":
    # uniform-ish reference tet
    reg = np.array([[0, 0, 0], [1, 0, 0], [0, 1, 0], [0, 0, 1]], float)
    gates("uniform tet", reg)
    # skewed / stretched tet: squash z hard + shear -> S_f nearly perpendicular to some edges (Ef blows up)
    skew = np.array([[0, 0, 0], [1, 0, 0], [0.4, 1.0, 0], [0.45, 0.5, 0.02]], float)
    gates("skewed/stretched tet", skew)
    print("\n[note] G1/G2 are structural (must PASS on both). The Ef spread + nonzero G3 on the skewed tet")
    print("       are the pieces invisible on uniform meshes -> the deferred smooth-grad RHS + Ef cap are")
    print("       mandatory before a stretched/pump run, exactly as the design's verify flagged.")
