#!/usr/bin/env python3
"""Increment-3 host replica: does the coupled collocated formulation reproduce the
published lid-driven SKEWED-cavity benchmark?

2D P1-triangle coupled [u,v,p] solver (same block algebra as the validated 3D tet
assembly, in the benchmark's native 2D), PSPG-stabilized equal-order, Picard outer
loop for the convection. Compares the centerline u-profile (line A-B) against the
exact Erturk-Dursun (arXiv:physics/0505121) reference for beta=45 deg, Re=100.

Blocks (per element, area A, integral_e N_i = A/3):
  a^uu = nu*K (comp-diagonal) + convection C_ij=(A/3)(a_e.gradN_j)
  a^up = +(A/3) gradN[b][d]     a^pu = -(A/3) gradN[a][d]     a^pp = tau*K
  K_ab = A*(gradN_a . gradN_b)
This is a *direct reference solve* (dense LU), the ground truth for the formulation.
"""
import numpy as np

# Erturk-Dursun Table 5, alpha=45 col: u along A-B, Re=100 (grid index 0..512; pos=idx/512)
ED_IDX = np.array([0,32,64,96,128,160,192,224,256,288,320,352,384,416,448,480,512]) / 512.0
ED_U45_RE100 = np.array([0.0,-2.389e-3,-8.901e-3,-1.949e-2,-3.426e-2,-5.373e-2,-7.846e-2,
                         -1.080e-1,-1.390e-1,-1.634e-1,-1.669e-1,-1.315e-1,-3.945e-2,
                         1.217e-1,3.544e-1,6.477e-1,1.0])


def tri_grad(r):
    """r:(3,2) -> (area, gradN (3,2)). Positive area."""
    T = np.array([r[1] - r[0], r[2] - r[0]])    # 2x2, rows = edges
    det = np.linalg.det(T)
    if det < 0:
        r = r[[0, 2, 1]]                         # flip to CCW
        T = np.array([r[1] - r[0], r[2] - r[0]])
        det = np.linalg.det(T)
        flipped = True
    else:
        flipped = False
    A = det / 2.0
    G = np.linalg.inv(T)
    dN = np.zeros((3, 2))
    dN[1] = G[:, 0]; dN[2] = G[:, 1]; dN[0] = -(dN[1] + dN[2])
    return A, dN, flipped


def skew_mesh(n, beta_deg):
    b = np.radians(beta_deg); cb, sb = np.cos(b), np.sin(b)
    P = np.zeros(((n + 1) * (n + 1), 2))
    for j in range(n + 1):
        for i in range(n + 1):
            x0, y0 = i / n, j / n
            P[j * (n + 1) + i] = [x0 + y0 * cb, y0 * sb]
    tris = []
    for j in range(n):
        for i in range(n):
            a = j * (n + 1) + i; bb = a + 1; c = a + (n + 1) + 1; d = a + (n + 1)
            tris.append([a, bb, c]); tris.append([a, c, d])
    return P, np.array(tris)


def assemble(P, tris, nu, tau, avel):
    N = len(P); ND = 3 * N
    A = np.zeros((ND, ND)); rhs = np.zeros(ND)
    for t in tris:
        r = P[t]
        Ae, dN, flip = tri_grad(r)
        tl = t[[0, 2, 1]] if flip else t          # node order matching dN
        ae = avel[tl].mean(axis=0)
        for a in range(3):
            ia = tl[a]
            for b in range(3):
                ib = tl[b]
                Kab = Ae * float(dN[a] @ dN[b])
                Cab = (Ae / 3.0) * float(ae @ dN[b])
                for d in range(2):
                    A[3 * ia + d, 3 * ib + d] += nu * Kab + Cab      # a^uu
                    A[3 * ia + d, 3 * ib + 2] += (Ae / 3.0) * dN[b][d]   # a^up
                    A[3 * ia + 2, 3 * ib + d] += -(Ae / 3.0) * dN[a][d]  # a^pu
                A[3 * ia + 2, 3 * ib + 2] += tau * Kab               # a^pp
    return A, rhs


def apply_bcs(A, rhs, P, n, ulid):
    """lid (top) u=ulid,v=0; other walls u=v=0; pin one pressure."""
    tol = 1e-9
    ymax = P[:, 1].max()
    for nd in range((n + 1) * (n + 1)):
        i = nd % (n + 1); j = nd // (n + 1)
        is_lid = (j == n)
        is_wall = (j == 0 or i == 0 or i == n) and not is_lid
        if is_lid or is_wall:
            uval = ulid if is_lid else 0.0
            for d, val in ((0, uval), (1, 0.0)):
                dof = 3 * nd + d
                A[dof, :] = 0.0; A[dof, dof] = 1.0; rhs[dof] = val
    pin = 3 * 0 + 2                                  # pin pressure at node 0
    A[pin, :] = 0.0; A[pin, pin] = 1.0; rhs[pin] = 0.0
    return A, rhs


def run(n=40, beta=45.0, Re=100.0, picard=25):
    P, tris = skew_mesh(n, beta)
    N = len(P)
    nu = 1.0 / Re                                    # U=1, L=1
    h = 1.0 / n
    tau = h * h / (4.0 * nu)                          # PSPG (Stokes-limit scaling)
    avel = np.zeros((N, 2))
    x = None
    for it in range(picard):
        A, rhs = assemble(P, tris, nu, tau, avel)
        A, rhs = apply_bcs(A, rhs, P, n, ulid=1.0)
        x = np.linalg.solve(A, rhs)
        unew = np.column_stack([x[0::3], x[1::3]])
        du = np.abs(unew - avel).max()
        avel = unew
        if du < 1e-8:
            break
    print(f"[cavity] n={n} beta={beta} Re={Re} picard={it+1} d(u)={du:.2e} "
          f"u_range=[{avel[:,0].min():.3f},{avel[:,0].max():.3f}]")

    # centerline A-B: column i=n/2, j=0..n; position s=j/n
    ic = n // 2
    s = np.array([j / n for j in range(n + 1)])
    ucl = np.array([avel[j * (n + 1) + ic, 0] for j in range(n + 1)])
    ed = np.interp(s, ED_IDX, ED_U45_RE100)
    err = np.abs(ucl - ed)
    print(f"[cavity] centerline u vs Erturk-Dursun(b45,Re100): max|err|={err.max():.3e}  "
          f"L2={np.sqrt(np.mean((ucl-ed)**2)):.3e}")
    print(f"[cavity] sample (s, u_mine, u_ED):")
    for k in range(0, n + 1, max(1, n // 8)):
        print(f"    s={s[k]:.2f}  mine={ucl[k]:+.4f}  ED={ed[k]:+.4f}")
    # qualitative shape check: min in lower-half, rises to ~1 at lid
    shape_ok = (ucl[:n//2].min() < -0.05) and (ucl[-1] > 0.9) and (ucl.argmin() < 0.8 * n)
    print(f"[cavity] qualitative shape (neg dip in lower half, ->1 at lid): "
          f"{'OK' if shape_ok else 'NO'}")


if __name__ == "__main__":
    run()
