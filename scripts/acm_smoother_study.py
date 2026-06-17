#!/usr/bin/env python3
"""ACM host-prototype (Phase 1c de-risk): can a GPU-friendly polynomial smoother give the
mesh-INDEPENDENCE that l1-Jacobi cannot?

1a showed Hypre point-block BoomerAMG converges but grows ~x3 with mesh, because the only
GPU-viable Hypre smoother is weak l1-Jacobi (strong GS is a no-op on GPU). The GPU ACM's
plan is unsmoothed agglomeration (sum-of-fine + injection, the darwish2008 form) + a
GPU-friendly smoother. This tests, on the SAME coupled cavity operator, the ACM family
(pyamg smooth=None = injection) under two smoothers:

  ACM+Jacobi  : weighted Jacobi (mimics the GPU BoomerAMG l1-Jacobi ceiling) -> expect growth
  ACM+Cheby   : Chebyshev polynomial (runs fine on GPU; matvec-only) -> expect flat

If ACM+Cheby is flat where ACM+Jacobi grows, the GPU ACM's smoother choice is justified and
it will beat the 1a baseline. (pyamg signal; definitive answer is the GPU ACM build.)
Run: /tmp/p1env/bin/python scripts/acm_smoother_study.py
"""
import importlib.util
import numpy as np
import pyamg

spec = importlib.util.spec_from_file_location("p1", "scripts/phase1_precond_study.py")
p1 = importlib.util.module_from_spec(spec); spec.loader.exec_module(p1)   # reuses build()/gmres_iters()/_block_B()

JAC  = ('jacobi',    {'iterations': 2, 'omega': 4.0 / 3.0})   # weak, ~ GPU l1-Jacobi
CHEB = ('chebyshev', {'degree': 3, 'iterations': 1})          # polynomial, GPU-friendly


def acm(A, smoother):
    # smooth=None -> unsmoothed aggregation = injection prolongation + sum-of-fine coarse
    # = additive-correction agglomeration MG (the darwish2008 ACM form), block near-null-space.
    return pyamg.smoothed_aggregation_solver(
        A, B=p1._block_B(A.shape[0]), smooth=None, max_coarse=50,
        presmoother=smoother, postsmoother=smoother).aspreconditioner('V')


def run():
    ns = [16, 24, 32, 40]
    mats = {n: p1.build(n) for n in ns}
    print(f"{'precond':<14} " + "  ".join(f"n={n}" for n in ns))
    rows = {}
    for name, sm in [("ACM+Jacobi", JAC), ("ACM+Cheby", CHEB)]:
        out = []
        for n in ns:
            A, b = mats[n]
            try:
                out.append(p1.gmres_iters(A, b, M=acm(A, sm)))
            except Exception as e:
                out.append(f"ERR:{type(e).__name__}")
        rows[name] = out
        print(f"{name:<14} " + "  ".join(f"{str(v):>5}" for v in out))

    print()
    for name in ["ACM+Jacobi", "ACM+Cheby"]:
        v = rows[name]
        if all(isinstance(x, int) and x > 0 for x in v):
            g = v[-1] / max(1, v[0])
            print(f"[verdict] {name}: {v} -> {'MESH-INDEPENDENT' if g < 1.6 else 'grows'} (x{g:.1f})")
        else:
            print(f"[verdict] {name}: {v} -> errored")
    print("[note] If ACM+Cheby is flat where ACM+Jacobi grows, the GPU ACM's polynomial-smoother")
    print("       choice recovers the mesh-independence stock GPU BoomerAMG (l1-Jacobi) cannot.")


if __name__ == "__main__":
    run()
