#!/usr/bin/env python3
"""ACM Stage 2 (host): the DECISIVE stretched-mesh GO/NO-GO. De-rigged per review.

The whole justification for building the paper's faithful (unsmoothed/injection) ACM. The host
studies showed it grows x1.7-2.0 on ISOTROPIC meshes while smoothed-SA is flat (x1.4), so the
paper's form is only worth building if its DIRECTIONAL (Mavriplis) agglomeration beats smoothed-SA
on a STRETCHED/anisotropic (boundary-layer) mesh. We test three hierarchies in the SAME V-cycle
with the SAME point-Jacobi smoother (so the ONLY variable is aggregation + prolongation):

  dir_inj : directional aggregation + order-zero injection P   (THE PAPER)
  iso_inj : isotropic distance-1 aggregation + injection P     (strawman)
  iso_sa  : isotropic aggregation + SMOOTHED P                 (THE REAL COMPETITOR)

GO  = paper (dir_inj) growth g<1.6 AND it beats iso_sa on the stretched mesh -> build faithful ACM.
NO-GO = paper not flat, OR smoothed-SA matches/beats it -> smoothed-SA is the 'something better';
        escalate to the user (do NOT silently switch -- 'paper form first' was the instruction).

Also tracks aggregate ANISOTROPY (do directional aggregates elongate along the stretch direction,
as Mavriplis intends). Run: /tmp/p1env/bin/python scripts/acm_stage2_stretched.py
"""
import importlib.util
import numpy as np
import scipy.sparse as sp

s1spec = importlib.util.spec_from_file_location("s1", "scripts/acm_stage1_vcycle.py")
s1 = importlib.util.module_from_spec(s1spec); s1spec.loader.exec_module(s1)
s0, p1, cav = s1.s0, s1.p1, s1.cav
NDOF, KMAX = s1.NDOF, s1.KMAX


def stretched_mesh(n, p):
    """structured (n+1)^2, uniform x, POWER-graded y=(j/n)^p clustered at y=0 (wall aspect=n^(p-1))."""
    P = np.zeros(((n + 1) ** 2, 2))
    for j in range(n + 1):
        y = (j / n) ** p
        for i in range(n + 1):
            P[j * (n + 1) + i] = [i / n, y]
    tris = []
    for j in range(n):
        for i in range(n):
            a = j * (n + 1) + i; b = a + 1; c = a + (n + 1) + 1; d = a + (n + 1)
            tris.append([a, b, c]); tris.append([a, c, d])
    dy_min = P[(n + 1)][1] - P[0][1]                   # first cell height at the wall
    return P, np.array(tris), (1.0 / n) / dy_min       # aspect ratio dx/dy_min


def p_for_aspect(n, aspect):
    return 1.0 + np.log(aspect) / np.log(n)            # aspect = n^(p-1)


def build_stretched(n, p, Re=100.0):
    P, tris, ar = stretched_mesh(n, p)
    N = len(P); nu = 1.0 / Re; h = 1.0 / n; tau = h * h / (4.0 * nu)
    A, rhs = cav.assemble(P, tris, nu, tau, np.zeros((N, 2)))
    A, rhs = cav.apply_bcs(A, rhs, P, n, ulid=1.0)
    return sp.csr_matrix(A), rhs.copy(), P, tris, ar


def block_adjacency(A, ndof):
    Ac = A.tocoo(); N = A.shape[0] // ndof
    nb = [set() for _ in range(N)]
    for r, c in zip(Ac.row, Ac.col):
        I, J = r // ndof, c // ndof
        if I != J:
            nb[I].add(J)
    return [list(x) for x in nb]


def isotropic_aggregate(nb, N):
    """distance-1 greedy (no strength test) -- the isotropic baseline aggregation."""
    agg = np.full(N, -1, np.int64); k = 0
    for seed in range(N):
        if agg[seed] != -1:
            continue
        agg[seed] = k; cnt = 1
        for j in nb[seed]:
            if cnt >= KMAX:
                break
            if agg[j] == -1:
                agg[j] = k; cnt += 1
        k += 1
    for i in range(N):
        if agg[i] != -1:
            continue
        cand = [agg[j] for j in nb[i] if agg[j] != -1]
        agg[i] = cand[0] if cand else k
        if not cand:
            k += 1
    return agg, k


def smooth_prolongation(A, Pt, ndof, niter=10):
    """smoothed-aggregation P = (I - omega D^-1 A) P_tent, omega = 4/(3 rho(D^-1 A))."""
    d = A.diagonal().copy(); d[np.abs(d) < 1e-30] = 1e-30
    Dinv = sp.diags(1.0 / d); DA = (Dinv @ A).tocsr()
    x = np.random.RandomState(0).randn(A.shape[0])
    for _ in range(niter):
        x = DA @ x; x /= np.linalg.norm(x)
    rho = np.linalg.norm(DA @ x)
    return (Pt - (4.0 / (3.0 * rho)) * (DA @ Pt)).tocsr()


def build_hierarchy(A0, Pcoord, tris, mode):
    directional = mode.startswith('dir'); sa = mode.endswith('sa')
    levels = [A0]; Ps = []; A = A0; lvl = 0
    while A.shape[0] > s1.MAX_COARSE:
        N = A.shape[0] // NDOF
        if directional:
            if lvl == 0:
                agg, na = s0.directional_aggregate(Pcoord, tris)
            else:
                w = s1.algebraic_weights(A, NDOF)
                if not w:
                    break
                agg, na = s1.aggregate(w, N)
        else:
            agg, na = isotropic_aggregate(block_adjacency(A, NDOF), N)
        if na <= 1 or na * NDOF >= A.shape[0]:
            break
        Pt = s0.block_P(agg, na, NDOF)
        P = smooth_prolongation(A, Pt, NDOF) if sa else Pt
        A = (P.T @ A @ P).tocsr()
        Ps.append(P); levels.append(A); lvl += 1
        if lvl == 0:
            pass
    return levels, Ps


def aggregate_stats(agg, Pcoord):
    """mean nodes/aggregate and mean bbox aspect (y-extent/x-extent, capped) -- do aggregates align with the stretch?"""
    sizes, bbox = [], []
    for I in range(int(agg.max()) + 1):
        pts = Pcoord[agg == I]
        if len(pts) < 2:
            sizes.append(len(pts)); continue
        sizes.append(len(pts))
        dx = np.ptp(pts[:, 0]); dy = np.ptp(pts[:, 1])
        bbox.append(min((dy + 1e-12) / (dx + 1e-12), 100.0))   # >1 = elongated along stretch (y)
    return float(np.mean(sizes)), float(np.mean(bbox)) if bbox else 1.0


MODES = [('dir_inj', 'directional+inj (PAPER)'),
         ('iso_inj', 'isotropic+inj  (strawman)'),
         ('iso_sa',  'isotropic+SA   (COMPETITOR)')]


def solve_iters(A, b, Pc, tris, mode):
    levels, Ps = build_hierarchy(A, Pc, tris, mode)
    return p1.gmres_iters(A, b, M=s1.make_vcycle(levels, Ps, 'point')[0])


def run():
    # ---- Sweep 1: fixed n, INCREASING aspect -- does smoothed-SA blow up as cells stretch? ----
    n = 32
    print(f"[sweep 1] fixed n={n}, increasing wall aspect ratio (iters; lower=better):")
    print(f"  {'aspect':<9}{'dir_inj(PAPER)':<16}{'iso_inj':<10}{'iso_sa(SA)':<12}")
    for aspect in [1, 10, 100, 1000]:
        p = p_for_aspect(n, aspect)
        A, b, Pc, tris, ar = build_stretched(n, p)
        its = {m: solve_iters(A, b, Pc, tris, m) for m, _ in MODES}
        print(f"  {ar:<9.0f}{its['dir_inj']:<16}{its['iso_inj']:<10}{its['iso_sa']:<12}")

    # ---- Sweep 2: fixed high aspect ~200, INCREASING n -- mesh-independence (flatness) ----
    print(f"\n[sweep 2] fixed wall aspect ~=200, increasing n (growth = flatness):")
    iters = {m: [] for m, _ in MODES}
    for nn in [16, 24, 32, 40]:
        p = p_for_aspect(nn, 200.0)
        A, b, Pc, tris, ar = build_stretched(nn, p)
        line = f"  n={nn:<3} DOFs={A.shape[0]:<5} aspect={ar:6.0f}  "
        for m, _ in MODES:
            it = solve_iters(A, b, Pc, tris, m); iters[m].append(it)
            line += f"{m.split('_')[0]}+{m.split('_')[1]}={it:<5}"
        agg_d, _ = s0.directional_aggregate(Pc, tris)
        agg_i, _ = isotropic_aggregate(block_adjacency(A, NDOF), len(Pc))
        sd, bd = aggregate_stats(agg_d, Pc); si, bi = aggregate_stats(agg_i, Pc)
        line += f" | aggsize dir={sd:.1f}/iso={si:.1f} bbox-aspect dir={bd:.1f}/iso={bi:.1f}"
        print(line)

    print(f"\n  {'method':<28}{'iters':<22}{'growth'}")
    g = {}
    for m, name in MODES:
        v = iters[m]; ok = all(isinstance(x, int) and x > 0 for x in v)
        g[m] = (v[-1] / v[0]) if ok else float('inf')
        print(f"  {name:<28}{str(v):<22}{'x%.2f' % g[m] if ok else 'DIVERGED'}")

    # ---- honest multi-signal verdict (compare at the FINEST mesh + growth, not the n=16 tie) ----
    paper, iso, sa = iters['dir_inj'], iters['iso_inj'], iters['iso_sa']
    flat = g['dir_inj'] < 1.6
    beats_sa = isinstance(paper[-1], int) and isinstance(sa[-1], int) and paper[-1] < 0.9 * sa[-1]
    beats_iso = isinstance(paper[-1], int) and isinstance(iso[-1], int) and paper[-1] < 0.9 * iso[-1]
    print()
    print(f"[signals @ finest mesh] paper flat: {flat} (g={g['dir_inj']:.2f}) | "
          f"paper beats smoothed-SA: {beats_sa} ({paper[-1]} vs {sa[-1]}) | "
          f"directional beats isotropic-inj: {beats_iso} ({paper[-1]} vs {iso[-1]})")
    if flat and beats_sa and beats_iso:
        print("[VERDICT] GO -- the paper's DIRECTIONAL ACM is flat AND beats both smoothed-SA and")
        print("          isotropic injection on the stretched mesh -> build the faithful ACM (Stages 3-6).")
    elif beats_sa and not beats_iso:
        print("[VERDICT] PARTIAL -- the UNSMOOTHED (injection) approach beats smoothed-SA on stretching")
        print("          (so smoothed-SA is NOT the better choice), but directional aggregation does not")
        print("          yet separate from plain isotropic injection. Escalate: build unsmoothed ACM, but")
        print("          the directional-vs-isotropic value needs a sharper test (3D/tet, real |S|).")
    else:
        why = []
        if not flat: why.append(f"paper not flat (g={g['dir_inj']:.2f})")
        if not beats_sa: why.append("smoothed-SA matches/beats the paper")
        print(f"[VERDICT] NO-GO -- {'; '.join(why)}. Escalate (settled-tradeoff: 'paper form first').")


if __name__ == "__main__":
    run()
