#!/usr/bin/env python3
"""Phase-1 preconditioner study (host signal, NOT the definitive Hypre answer).

The existential question for the coupled solver: does ANY AMG-class preconditioner give
near mesh-INDEPENDENT GMRES iteration counts on the coupled collocated saddle operator?
Built on the same coupled cavity operator the host replica validated (Stokes, PSPG-
stabilized [u,v,p], lid-driven BCs). Tests at several mesh sizes:

  - none      : GMRES, no preconditioner (baseline; iters should grow ~1/h)
  - ILU       : scipy spilu
  - AMG-RS    : pyamg classical Ruge-Stuben (scalar)
  - AMG-sys   : pyamg smoothed-aggregation with the 3-DOF block near-null-space
                (the closest host analog to Hypre point-block / nodal-systems AMG)

CAVEAT: pyamg != Hypre point-block BoomerAMG; this is a SIGNAL to decide whether 1a
(Hypre point-block) is worth building or we go straight to 1c (agglomeration ACM).
Run with the throwaway venv that has scipy+pyamg:  /tmp/p1env/bin/python scripts/phase1_precond_study.py
"""
import sys, importlib.util
import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spla
import pyamg

spec = importlib.util.spec_from_file_location("cav", "scripts/coupled_cavity_reference.py")
cav = importlib.util.module_from_spec(spec); spec.loader.exec_module(cav)


def build(n, Re=100.0):
    P, tris = cav.skew_mesh(n, 45.0)
    N = len(P)
    nu = 1.0 / Re
    h = 1.0 / n
    tau = h * h / (4.0 * nu)
    A, rhs = cav.assemble(P, tris, nu, tau, np.zeros((N, 2)))   # Stokes (advection off)
    A, rhs = cav.apply_bcs(A, rhs, P, n, ulid=1.0)
    return sp.csr_matrix(A), rhs.copy()


def gmres_iters(A, b, M=None, rtol=1e-8, restart=300, maxiter=3000):
    cnt = [0]
    def cb(_): cnt[0] += 1
    x, info = spla.gmres(A, b, M=M, rtol=rtol, atol=0.0, restart=restart,
                         maxiter=maxiter, callback=cb, callback_type='pr_norm')
    return (cnt[0] if info == 0 else -1)   # -1 = did not converge


def _block_B(ND):
    B = np.zeros((ND, 3))                     # component-constant near-null-space (u,v,p blocks)
    for c in range(3):
        B[c::3, c] = 1.0
    return B


def systems_amg(A, smooth='default'):
    # smooth=None -> UNSMOOTHED aggregation = injection prolongation + sum-of-fine
    # coarse operator (R A R^T) = the additive-correction agglomeration MG (ACM analog).
    kw = dict(B=_block_B(A.shape[0]), max_coarse=50)
    if smooth is None:
        kw['smooth'] = None
    return pyamg.smoothed_aggregation_solver(A, **kw).aspreconditioner('V')


def run():
    ns = [16, 24, 32]
    print(f"{'precond':<10} " + "  ".join(f"n={n}(DOFs={3*(n+1)**2})" for n in ns))
    rows = {}
    mats = {n: build(n) for n in ns}
    for name in ["none", "ILU", "AMG-RS", "AMG-sys", "ACM-agg"]:
        out = []
        for n in ns:
            A, b = mats[n]
            try:
                if name == "none":
                    M = None
                elif name == "ILU":
                    M = spla.LinearOperator(A.shape, spla.spilu(A.tocsc()).solve)
                elif name == "AMG-RS":
                    M = pyamg.ruge_stuben_solver(A, max_coarse=50).aspreconditioner('V')
                elif name == "AMG-sys":
                    M = systems_amg(A)                       # smoothed aggregation (systems)
                else:
                    M = systems_amg(A, smooth=None)          # ACM analog: unsmoothed agglomeration
                it = gmres_iters(A, b, M=M)
            except Exception as e:
                it = f"ERR:{type(e).__name__}"
            out.append(it)
        rows[name] = out
        print(f"{name:<10} " + "  ".join(f"{str(v):>12}" for v in out))

    # verdict: an AMG-class row that stays ~flat (mesh-independent) answers the question YES
    print()
    for name in ["AMG-RS", "AMG-sys", "ACM-agg"]:
        v = rows[name]
        if all(isinstance(x, int) and x > 0 for x in v):
            growth = v[-1] / max(1, v[0])
            verdict = "MESH-INDEPENDENT (AMG works)" if growth < 2.0 else "iters GROW with mesh"
            print(f"[verdict] {name}: iters {v} -> {verdict} (x{growth:.1f})")
        else:
            print(f"[verdict] {name}: {v} -> did NOT converge / errored")
    print("[note] pyamg signal only; definitive answer is Hypre point-block BoomerAMG on GPU (1a).")


if __name__ == "__main__":
    run()
