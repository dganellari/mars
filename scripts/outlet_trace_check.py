#!/usr/bin/env python3
"""Acceptance gates for the average-pressure outlet trace.

Checks the algebraic properties the review requires of

    p_trace = p_ref + (1-beta) * (p_sample - mean_A(p_sample))

with mean_A the AREA-WEIGHTED mean over the outlet. These are the gates that need no GPU; the
CUDA/MPI ones (frozen-state Jacobian action, owned-continuity vs exterior flux, 1/2/4-rank
agreement) still have to run on Alps.

This mirrors mars_ns_pump_solver.hpp `outletTraceValue` + `buildOutletPressureTraceKernel`; it is a
replica of the algebra, not a call into the shipped kernel.
"""

import sys
import numpy as np

FAIL = []


def check(name, ok, detail=""):
    print(f"  {'PASS' if ok else 'FAIL'}  {name}{('  ' + detail) if detail else ''}")
    if not ok:
        FAIL.append(name)


def trace(p, area, p_ref, beta):
    tot = area.sum()
    mean = (area * p).sum() / tot if tot > 0 else p_ref
    return p_ref + (1.0 - beta) * (p - mean), mean


def main():
    rng = np.random.default_rng(20260909)
    n = 400
    area = rng.uniform(0.1, 3.0, n)          # deliberately NON-uniform: an unweighted mean would pass
    p = rng.normal(50.0, 12.0, n)            # a genuinely nonuniform outlet pressure
    beta, p_ref = 0.05, 7.5

    t, mean = trace(p, area, p_ref, beta)

    print("gate 1: prescribed weighted mean, retained variation")
    got = (area * t).sum() / area.sum()
    check("area-weighted mean of the trace equals p_ref", abs(got - p_ref) < 1e-10,
          f"got {got:.12f} vs {p_ref}")
    # beta=0.05 must RETAIN 95% of the fluctuation, not force 95% of it to the mean.
    ratio = np.std(t) / np.std(p)
    check("fluctuation retained = (1-beta)", abs(ratio - (1.0 - beta)) < 1e-12,
          f"std ratio {ratio:.12f} vs {1.0-beta}")
    check("trace is NOT uniform", np.std(t) > 0.5 * np.std(p), f"std {np.std(t):.4f}")

    print("gate 2: gauge invariance of the fluctuation")
    for c in (0.0, 1e3, -2.5e6):
        tc, _ = trace(p + c, area, p_ref, beta)
        check(f"adding {c:g} to volume pressure leaves trace-p_ref unchanged",
              np.allclose(tc - p_ref, t - p_ref, rtol=0, atol=1e-8))

    print("gate 3: area weighting actually matters")
    unweighted = p.mean()
    check("weighted mean differs from the unweighted one",
          abs(mean - unweighted) > 1e-6, f"{mean:.6f} vs {unweighted:.6f}")

    print("gate 4: partition invariance of the moment reduction")
    for nranks in (2, 3, 4):
        parts = np.array_split(np.arange(n), nranks)
        num = sum((area[q] * p[q]).sum() for q in parts)
        den = sum(area[q].sum() for q in parts)
        check(f"{nranks}-way split reproduces the global mean",
              abs(num / den - mean) < 1e-9)
    # A rank owning no outlet nodes must contribute (0,0) and still reduce.
    parts = [np.arange(n), np.array([], dtype=int)]
    num = sum((area[q] * p[q]).sum() for q in parts)
    den = sum(area[q].sum() for q in parts)
    check("a rank with zero outlet nodes changes nothing", abs(num / den - mean) < 1e-12)

    print("gate 5: empty / zero-area patch")
    empty = np.array([])
    t0, mean0 = trace(empty, empty, p_ref, beta)
    check("empty patch falls back to p_ref with no division by zero",
          mean0 == p_ref and t0.size == 0)
    zero_area = np.zeros(5)
    _, mz = trace(np.ones(5) * 33.0, zero_area, p_ref, beta)
    check("zero total area falls back to p_ref", mz == p_ref)

    print("gate 6: the degenerate case the review warns about")
    # Feeding a HARD p=0 face back through the map returns p_ref everywhere -- the face goes
    # uniform. This is why the trace must live in its own field and the unknowns stay free.
    clamped = np.zeros(n)
    tc, _ = trace(clamped, area, p_ref, beta)
    check("clamped p=0 face collapses the trace to a uniform p_ref",
          np.allclose(tc, p_ref) and np.std(tc) < 1e-12)

    print()
    if FAIL:
        print(f"{len(FAIL)} FAILED: {', '.join(FAIL)}")
        return 1
    print("all host gates passed; CUDA/MPI gates still pending on Alps")
    return 0


if __name__ == "__main__":
    sys.exit(main())
