#!/usr/bin/env python3
"""Host replica of the Phase-0 coupled [u,v,w,p] block assembly (Stokes, advection OFF).

Validates the COEFFICIENT ALGEBRA of internal-notes/phase0_coupled_assembly.md before
any GPU code is written -- the same host-replica methodology MARS used for the pressure
operator. Assembles the 4 sub-blocks the spec defines and runs the matrix-consistency
gates (no solver involved):

  a^up   = +G : momentum row (node a, comp d) <- pressure col (node b) : +(V/4) dNdx[b][d]
  a^pu   = D  : continuity row (node a)        <- momentum col (node b, comp d) : -(V/4) dNdx[a][d]
  a^uu   = nu*K (component-diagonal) [+ M/dt, dropped for steady Stokes]
  a^pp   = tau*K  (Galerkin stiffness scaled by tau -- NOT compact Rhie-Chow)
with K_ab = V*(dNdx[a].dNdx[b]).

Gates: (1) a^pu == -(a^up)^T, (2) a^pp.1 == 0 (const-p null mode), (3) D.(const u) == 0,
(4) single-tet 16x16 structure. Convention-robust: any correct linear-tet dNdx works,
since D and G share the same dNdx and sum_i dNdx_i = 0.
"""
import numpy as np

NU = 1.0e-3
TAU = 1.0 / 24.0   # tau ~ h^2/24; absolute scale irrelevant to the gates


def tet_dndx(r):
    """r: (4,3) node coords -> (V, dNdx (4,3))."""
    T = np.array([r[1] - r[0], r[2] - r[0], r[3] - r[0]])   # rows = edge vectors
    detT = np.linalg.det(T)
    V = detT / 6.0
    G = np.linalg.inv(T)                 # grad of barycentric coords 1,2,3 = columns of inv(T)
    dN = np.zeros((4, 3))
    dN[1] = G[:, 0]
    dN[2] = G[:, 1]
    dN[3] = G[:, 2]
    dN[0] = -(dN[1] + dN[2] + dN[3])
    return V, dN


def assemble(points, tets, steady=True, dt=1.0, avel=None):
    """Assemble the dense 4N x 4N coupled block matrix.

    avel=None -> Stokes (advection OFF). avel=(N,3) nodal velocity -> add the
    linearized (frozen-a) consistent-Galerkin convection into a^uu:
      C_ij = integral_e N_i (a.grad N_j) = (V/4) (a_e . dNdx[j]),  a_e = elem-mean(a).
    Same FEM/dNdx primitive as the other blocks; non-symmetric; row-sums to 0
    (advection conserves)."""
    N = len(points)
    A = np.zeros((4 * N, 4 * N))
    Mnode = np.zeros(N)
    for t in tets:
        r = points[t]
        V, dN = tet_dndx(r)
        if not (V > 0):
            raise SystemExit("non-positive tet volume")
        ae = avel[t].mean(axis=0) if avel is not None else None
        for a in range(4):
            Mnode[t[a]] += V / 4.0
        for a in range(4):
            ia = t[a]
            for b in range(4):
                ib = t[b]
                Kab = V * float(dN[a] @ dN[b])            # scalar stiffness
                for d in range(3):
                    A[4 * ia + d, 4 * ib + d] += NU * Kab   # a^uu diffusion (comp-diagonal)
                A[4 * ia + 3, 4 * ib + 3] += TAU * Kab        # a^pp = tau*K
                for d in range(3):
                    A[4 * ia + d, 4 * ib + 3] += (V / 4.0) * dN[b][d]   # a^up = +G
                    A[4 * ia + 3, 4 * ib + d] += -(V / 4.0) * dN[a][d]  # a^pu = -G^T
                if ae is not None:
                    Cab = (V / 4.0) * float(ae @ dN[b])    # convection C_ij = (V/4)(a.gradN_j)
                    for d in range(3):
                        A[4 * ia + d, 4 * ib + d] += Cab    # a^uu advection (comp-diagonal)
    if not steady:
        for i in range(N):
            for d in range(3):
                A[4 * i + d, 4 * i + d] += Mnode[i] / dt
    return A, Mnode


def cube_6tet():
    """Unit cube, 8 nodes, split into 6 tets (shared 0-6 diagonal) -> interior coupling."""
    P = np.array([[0, 0, 0], [1, 0, 0], [1, 1, 0], [0, 1, 0],
                  [0, 0, 1], [1, 0, 1], [1, 1, 1], [0, 1, 1]], dtype=float)
    K = [(0, 1, 2, 6), (0, 2, 3, 6), (0, 3, 7, 6), (0, 7, 4, 6), (0, 4, 5, 6), (0, 5, 1, 6)]
    tets = []
    for t in K:                                  # orient positive
        r = P[list(t)]
        if np.linalg.det(np.array([r[1] - r[0], r[2] - r[0], r[3] - r[0]])) < 0:
            t = (t[0], t[2], t[1], t[3])
        tets.append(t)
    return P, np.array(tets)


def gate(name, ok, detail):
    print(f"[replica][{name}] {'PASS' if ok else 'FAIL'}  ({detail})")
    return ok


def run():
    P, tets = cube_6tet()
    N = len(P)
    A, _ = assemble(P, tets, steady=True)
    TOL = 1e-12
    allok = True

    # extract blocks
    pres = [4 * i + 3 for i in range(N)]
    vel = [4 * i + d for i in range(N) for d in range(3)]
    A_pu = A[np.ix_(pres, vel)]
    A_up = A[np.ix_(vel, pres)]
    A_pp = A[np.ix_(pres, pres)]

    # (1) a^pu == -(a^up)^T
    e1 = np.abs(A_pu + A_up.T).max()
    allok &= gate("a_pu==-a_up^T", e1 < TOL, f"max|a_pu + a_up^T| = {e1:.3e}")

    # (2) a^pp . 1 == 0  (constant-pressure null mode of tau*K)
    e2 = np.abs(A_pp @ np.ones(N)).max()
    allok &= gate("a_pp.1==0", e2 < TOL, f"max|a_pp @ 1| = {e2:.3e}")

    # (3) D . (constant u) == 0 on interior: continuity rows of A applied to u=1,p=0
    x = np.zeros(4 * N)
    for i in range(N):
        for d in range(3):
            x[4 * i + d] = 1.0
    contin = A[np.ix_(pres, range(4 * N))] @ x      # = a^pu @ (const u)
    # cube has no interior node (all 8 are boundary), so this is the boundary-flux sum;
    # the true interior check is the global sum == 0 (closed divergence theorem)
    e3 = abs(contin.sum())
    allok &= gate("sum(D.const)==0", e3 < TOL, f"|sum_i (a_pu u)_i| = {e3:.3e} (global closure)")

    # (4) single-tet 16x16 structure: reference tet, print block norms
    Pt = np.array([[0, 0, 0], [1, 0, 0], [0, 1, 0], [0, 0, 1]], dtype=float)
    At, _ = assemble(Pt, np.array([(0, 1, 2, 3)]), steady=True)
    pres1 = [4 * i + 3 for i in range(4)]
    vel1 = [4 * i + d for i in range(4) for d in range(3)]
    s_pu = At[np.ix_(pres1, vel1)]
    s_up = At[np.ix_(vel1, pres1)]
    e4 = np.abs(s_pu + s_up.T).max()
    sym_uu = np.abs(At[np.ix_(vel1, vel1)] - At[np.ix_(vel1, vel1)].T).max()
    allok &= gate("single-tet 16x16", e4 < TOL and sym_uu < TOL,
                  f"max|a_pu+a_up^T|={e4:.3e}, a_uu symmetric residual={sym_uu:.3e} (Stokes a_uu must be SPD-symmetric)")

    print(f"\n[replica] OVERALL: {'ALL GATES PASS' if allok else 'FAILURE'}")
    return 0 if allok else 1


if __name__ == "__main__":
    raise SystemExit(run())
