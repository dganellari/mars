#!/usr/bin/env python3
"""Host checks of the multi-rank periodic projection (no GPU, no MPI).

    python3 tests/periodic/check_periodic_space.py            # all checks, 8^3 mesh
    python3 tests/periodic/check_periodic_space.py --reference --n 16 --steps 300

Part 1 transcribes the pre-fix multi-rank periodic operators of
mars_ns_solver.hpp and shows where they stop being consistent.
Part 2 checks the invariants of the reduced periodic space
(mars_periodic_space.hpp): rank invariance, symmetry, exact projection.
Part 3 runs the PeriodicNavierStokes time step on 1, 2 and 4 simulated ranks.
--reference prints the kinetic energy of the GPU validation case for comparison.
"""
import argparse
import numpy as np
from periodic_space_model import Mesh, partition, Dist, PeriodicNS, seam_mismatch, tgv_initial, smooth_field


def rel(a, b):
    return np.max(np.abs(a - b)) / max(np.max(np.abs(b)), 1e-300)


def group_field(mesh, seed, ncomp=None):
    if ncomp is None:
        return smooth_field(mesh, seed)[mesh.master]
    return np.stack([smooth_field(mesh, seed + c) for c in range(ncomp)], axis=1)[mesh.master]


# ---------------------------------------------------------------------------
# Part 1: pre-fix kernels (mars_ns_solver.hpp / mars_periodic_bc.hpp)
# ---------------------------------------------------------------------------
def bcast_ungated(d, F):
    """periodicBroadcastKernel: every slot with a partner copies the partner slot."""
    for rk in d.ranks:
        s = np.nonzero(rk.partner >= 0)[0]
        F[rk.r][s] = F[rk.r][rk.partner[s]]


def bcast_same_rank(d, F):
    """periodicBroadcastSameRankKernel: every slot whose partner is owned here."""
    for rk in d.ranks:
        s = np.nonzero(rk.partner >= 0)[0]
        s = s[rk.own[rk.partner[s]]]
        F[rk.r][s] = F[rk.r][rk.partner[s]]


def bcast_cross_rank(d, F):
    for rk in d.ranks:
        for q in np.unique(rk.crossQ):
            sel = rk.crossQ == q
            F[rk.r][rk.crossS[sel]] = F[q][rk.crossM[sel]]


def fold_ungated(d, F):
    """periodicFoldToMasterKernel, applied before the reverse halo."""
    for rk in d.ranks:
        s = np.nonzero(rk.partner >= 0)[0]
        np.add.at(F[rk.r], rk.partner[s], F[rk.r][s])
        F[rk.r][s] = 0


def normalize_owned(d, F, massNode):
    """normalizeGradientPerNodeKernel: owned slots / mass, 0 for ghosts or m == 0."""
    out = []
    for rk, f, m in zip(d.ranks, F, massNode):
        g = np.zeros_like(f)
        ok = rk.own & (m != 0)
        g[ok] = (f[ok].T / m[ok]).T
        out.append(g)
    return out


def legacy_mass(d):
    """setupNSStepper: fold, halo, then mirror onto slaves with the ungated broadcast."""
    m = d.scatter_mass()
    d.restrict(m)
    d.exchange(m)
    bcast_ungated(d, m)
    d.exchange(m)
    return m


def assembled_velocity_rows(mesh, d):
    """Owner-migration assembly over owned+halo elements into owned rows only.

    A cross-rank slave's row maps to its master's ghost DOF and is dropped; on
    the master's rank the slave is a ghost, whose rows are skipped too."""
    X = group_field(mesh, 7)
    trueK = d.scatter_K(d.from_global(X))
    d.restrict(trueK)
    yTrue = d.to_global(trueK)
    yAsm = np.zeros(mesh.nNodes)
    seam = np.zeros(mesh.nNodes, bool)
    for rk in d.ranks:
        for le in rk.g2l[mesh.conn[rk.localElems]]:
            for a in range(8):
                i = le[a]
                if not rk.own[i]:
                    continue
                g = rk.nodes[i]
                if mesh.isSlave[g] and (rk.partner[i] < 0 or not rk.own[rk.partner[i]]):
                    continue
                yAsm[mesh.master[g]] += mesh.Ke[a] @ X[rk.nodes[le]]
        seam[mesh.master[rk.nodes[rk.crossS]]] = True
    red = ~mesh.isSlave
    err = np.abs(yAsm - yTrue)
    scale = np.max(np.abs(yTrue[red]))
    return (np.max(err[red & seam]) / scale if seam.any() else 0.0), np.max(err[red & ~seam]) / scale


def full_fold_operator(d, massNode, haloAfterBroadcast):
    """applyDDTPerNode(applyPeriodic=true): P on phi, fold+normalize g, copy g to slaves, D, fold."""
    def op(P):
        bcast_same_rank(d, P)
        bcast_cross_rank(d, P)
        d.exchange(P)
        g = d.scatter_divT(P)
        d.restrict(g)
        g = normalize_owned(d, g, massNode)
        if haloAfterBroadcast:
            bcast_same_rank(d, g)
            bcast_cross_rank(d, g)
            d.exchange(g)
        else:
            d.exchange(g)
            bcast_same_rank(d, g)
            bcast_cross_rank(d, g)
        out = d.scatter_div(g)
        d.restrict(out)
        return out
    return op


def symmetry(d, op, mesh):
    u, v = d.from_global(group_field(mesh, 1)), d.from_global(group_field(mesh, 2))
    Av, Au = op([x.copy() for x in v]), op([x.copy() for x in u])
    a, b = d.rdot(u, Av), d.rdot(v, Au)
    return abs(a - b) / max(abs(a), abs(b))


def reduced_path_projection(mesh, d, rho, dtEff, slaveWritesLast, tol=1e-10):
    """Pre-fix multi-rank default: bare reduced operator, per-slot corrector,
    then the u[slave] := u[master] copy. Returns |D u|/|D u**| as PROJ-P3 measures it."""
    massNode = legacy_mass(d)
    scratch = d.zeros()

    def P_reduced(X):
        for rk, s in zip(d.ranks, scratch):
            s[rk.reduced] = X[rk.r][rk.reduced]
            s[rk.sameS] = X[rk.r][rk.sameM]
        d.exchange(scratch)
        bcast_ungated(d, scratch)
        return [s.copy() for s in scratch]

    def reduce_legacy(F):
        fold_ungated(d, F)
        d.reverse_add(F)

    def A_red(X):
        g = d.scatter_divT(P_reduced(X))
        d.reverse_add(g)
        g = normalize_owned(d, g, massNode)
        d.exchange(g)
        out = d.scatter_div(g)
        reduce_legacy(out)
        Y = d.zeros()
        for rk, y, o in zip(d.ranks, Y, out):
            y[rk.reduced] = o[rk.reduced]
            np.add.at(y, rk.sameM, o[rk.sameS])
        return Y

    def proj_p3_div(U):
        acc = d.scatter_div(U)
        reduce_legacy(acc)
        n = normalize_owned(d, acc, massNode)
        return max(np.max(np.abs(x[rk.own])) for rk, x in zip(d.ranks, n))

    U = d.from_global(group_field(mesh, 11, 3))
    acc = d.scatter_div(U)
    reduce_legacy(acc)
    b = d.zeros()
    for rk, bb, a in zip(d.ranks, b, acc):
        bb[rk.reduced] = -(rho / dtEff) * a[rk.reduced]
        if slaveWritesLast:
            # buildPressureRhsKernel assigns rhs[dof] from BOTH slots of a
            # same-rank pair; which write survives is not defined on a GPU
            bb[rk.sameM] = -(rho / dtEff) * a[rk.sameS]
    mean = d.rsum(b) / sum(int(np.sum(rk.reduced)) for rk in d.ranks)
    b = [bb - mean for bb in b]
    x, r = d.zeros(), [bb.copy() for bb in b]
    p = [rr.copy() for rr in r]
    rr = d.rdot(r, r)
    r0 = np.sqrt(rr)
    its = 0
    while np.sqrt(rr) > tol * r0 and its < 3000:
        Ap = A_red(p)
        alpha = rr / d.rdot(p, Ap)
        x = [xi + alpha * pi for xi, pi in zip(x, p)]
        r = [ri - alpha * ai for ri, ai in zip(r, Ap)]
        rrn = d.rdot(r, r)
        p = [ri + (rrn / rr) * pi for ri, pi in zip(r, p)]
        rr, its = rrn, its + 1
    phi = P_reduced(x)
    g = d.scatter_divT(phi)
    d.reverse_add(g)
    g = normalize_owned(d, g, massNode)
    d.exchange(g)
    Un = []
    for rk, u, gg in zip(d.ranks, U, g):
        un = u.copy()
        un[rk.own] = u[rk.own] + (dtEff / rho) * gg[rk.own]
        Un.append(un)
    d.exchange(Un)
    dIn, perSlot = proj_p3_div(U), proj_p3_div(Un)
    bcast_same_rank(d, Un)
    bcast_cross_rank(d, Un)
    d.exchange(Un)
    return perSlot / dIn, proj_p3_div(Un) / dIn, its


def part1(mesh):
    print("== Part 1: pre-fix multi-rank periodic operators (every ghost refreshed by the halo) ==")
    for R in (2, 4):
        d = Dist(mesh, partition(mesh, R, fullTopo=True))
        seamErr, bulkErr = assembled_velocity_rows(mesh, d)
        mass = PeriodicNS(d, 0.05, 1.0, 1e-3, True).mass
        symBad = symmetry(d, full_fold_operator(d, mass, False), mesh)
        symGood = symmetry(d, full_fold_operator(d, mass, True), mesh)
        print(f"R={R}: assembled velocity rows vs P^T K P: max rel error {seamErr:.2e} at cross-rank periodic "
              f"points, {bulkErr:.1e} elsewhere")
        print(f"      full-fold D M^-1 D^T symmetry: halo then slave copy {symBad:.1e}, slave copy then halo "
              f"{symGood:.1e}")
        for last in (False, True):
            perSlot, p3, its = reduced_path_projection(mesh, d, 1.0, 2e-3 / 3, last)
            print(f"      reduced operator, rhs race {'slave' if last else 'master'} wins: cg={its} "
                  f"per-slot |Du|/|Du**|={perSlot:.1e}, after u[S]:=u[M] PROJ-P3={p3:.2f}")


def part2(mesh):
    print("\n== Part 2: reduced periodic space, halo refreshes only corners of owned elements ==")
    ref = None
    for R in (1, 2, 3, 4):
        d = Dist(mesh, partition(mesh, R))
        ns = PeriodicNS(d, 0.05, 1.0, 1e-3, True)
        f, V = group_field(mesh, 3), group_field(mesh, 4, 3)
        out = dict(grad=d.to_global(ns.gradient(d.from_global(f)), 3),
                   div=d.to_global(ns.divergence(d.from_global(V))),
                   ddt=d.to_global(ns.pressure_op(d.from_global(f))),
                   visc=d.to_global(ns.visc_op(1e3)(d.from_global(f))))
        ref = ref or out
        inv = " ".join(f"{k}={rel(out[k], ref[k]):.0e}" for k in out)
        symP, symK = symmetry(d, ns.pressure_op, mesh), symmetry(d, ns.visc_op(1e3), mesh)
        U = d.from_global(V)
        Un, _, its = ns.project(U, 1e-3)
        ratio = d.rmax(ns.divergence(Un)) / d.rmax(ns.divergence(U))
        d.prolong(Un)
        adv = d.scatter_adv(Un, True)
        d.restrict(adv)
        upw = d.scatter_adv(Un, False)
        d.restrict(upw)
        ke = abs(sum(d.rdot([u[:, c] for u in Un], [a[:, c] for a in adv]) for c in range(3)))
        ku = abs(sum(d.rdot([u[:, c] for u in Un], [a[:, c] for a in upw]) for c in range(3)))
        print(f"R={R}: rank invariance {inv} | symmetry D M^-1 G {symP:.0e}, M+nuK {symK:.0e}")
        print(f"      projection cg={its} |D u^(n+1)|/|D u**|={ratio:.1e} slave-master mismatch "
              f"{seam_mismatch(d, Un):.0e} | skew KE production / upwind {ke / ku:.0e}")


def run_tgv(mesh, R, steps, nu, dt, skew, legacy=False):
    d = Dist(mesh, partition(mesh, R))
    ns = PeriodicNS(d, nu, 1.0, dt, skew, legacyPressureCoef=legacy)
    U, P = tgv_initial(mesh)
    ns.start(d.from_global(U), d.from_global(P))
    rows = [(0, ns.kinetic_energy(), ns.max_divergence(), 0.0, (0, 0, 0, 0), d.to_global(ns.U, 3))]
    for n in range(1, steps + 1):
        ns.step()
        rows.append((n, ns.kinetic_energy(), ns.max_divergence(), seam_mismatch(d, ns.U), ns.its,
                     d.to_global(ns.U, 3)))
    return rows


def part3(mesh, steps):
    nu, dt = 0.05, 1e-3
    print(f"\n== Part 3: TGV on {mesh.N}^3, nu={nu}, dt={dt}, skew advection, BDF1 then BDF2 ==")
    runs = {R: run_tgv(mesh, R, steps, nu, dt, True) for R in (1, 2, 4)}
    for n in range(steps + 1):
        k1 = runs[1][n][1]
        line = f"step {n}: KE={k1:.12f}"
        for R in (2, 4):
            line += f"  R={R}: dKE/KE={abs(runs[R][n][1] - k1) / k1:.0e} max|du|={np.max(np.abs(runs[R][n][5] - runs[1][n][5])):.0e}"
        line += f"  div={max(runs[R][n][2] for R in runs):.0e} slave-master={max(runs[R][n][3] for R in runs):.0e}"
        line += f"  cg(u,v,w,p)={'/'.join(map(str, runs[4][n][4]))}"
        print(line)


def reference(mesh, steps, every):
    nu, dt = 0.05, 1e-4
    k = 2 * np.pi
    print(f"== Reference: {mesh.N}^3 unit box, nu={nu}, dt={dt}, V0=1, rho=1 (mars_tgv --box-lo=0 --box-hi=1) ==")
    for skew in (True, False):
        for legacy in ((False, True) if skew else (False,)):
            rows = run_tgv(mesh, 1, steps, nu, dt, skew, legacy)
            ke0 = rows[0][1]
            label = ("skew" if skew else "upwind") + (", pre-fix BDF2 pressure coefficient rho/dt" if legacy else "")
            print(f"-- advection {label}")
            for n, ke, div, _, its, _ in rows:
                if n % every == 0:
                    t = n * dt
                    print(f"   step {n:5d} t={t:.4f} KE={ke:.10e} KE/KE_Stokes={ke / (ke0 * np.exp(-6 * nu * k * k * t)):.8f}"
                          f" div={div:.1e} cg(u,v,w,p)={'/'.join(map(str, its))}")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--n", type=int, default=8, help="elements per box edge")
    ap.add_argument("--steps", type=int, default=6)
    ap.add_argument("--reference", action="store_true", help="print the GPU validation reference only")
    ap.add_argument("--every", type=int, default=50, help="reference print interval")
    args = ap.parse_args()
    mesh = Mesh(args.n)
    if args.reference:
        reference(mesh, args.steps, args.every)
        return
    part1(mesh)
    part2(mesh)
    part3(mesh, args.steps)


if __name__ == "__main__":
    main()
