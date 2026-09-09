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

    print("gate 7: exact values from the GPT spec (docs/design/gpt_outlet_boundary_spec_2026-09-09.md)")
    sp = np.array([2.0, 4.0, 8.0])
    sa = np.array([1.0, 2.0, 1.0])
    for b, want in ((0.05, [7.625, 9.525, 13.325]),
                    (0.0,  [7.5, 9.5, 13.5]),
                    (1.0,  [10.0, 10.0, 10.0])):
        t7, m7 = trace(sp, sa, 10.0, b)
        check(f"beta={b}: trace matches the spec", np.allclose(t7, want, atol=1e-12),
              f"{np.round(t7,6).tolist()} vs {want}")
    check("spec mean is 4.5", abs(m7 - 4.5) < 1e-12)
    t7a, _ = trace(sp + 100.0, sa, 10.0, 0.05)
    t7b, _ = trace(sp, sa, 10.0, 0.05)
    check("adding 100 Pa to every sample leaves the trace unchanged",
          np.allclose(t7a, t7b, atol=1e-12))

    print("gate 8: BENT surface -- scalar area vs the norm of the summed area vector")
    # Node 0 carries two PERPENDICULAR unit-area facet contributions, node 1 carries one.
    # Scalar area (correct):     |A|sum = 1+1 = 2   and 1
    # Norm of the summed vector: |(1,0,0)+(0,1,0)| = sqrt(2)  and 1   <- what the old code used
    pb = np.array([0.0, 3.0])
    a_scalar = np.array([2.0, 1.0])
    a_norm = np.array([np.sqrt(2.0), 1.0])
    _, m_true = trace(pb, a_scalar, 10.0, 0.05)
    _, m_bad = trace(pb, a_norm, 10.0, 0.05)
    check("true scalar-area mean is 1", abs(m_true - 1.0) < 1e-12, f"got {m_true}")
    check("vector-norm mean is the spec's wrong value",
          abs(m_bad - 1.242640687119285) < 1e-12, f"got {m_bad:.15f}")
    # The consequence: the trace no longer has the prescribed physical mean.
    t_bad = 10.0 + 0.95 * (pb - m_bad)
    phys = (a_scalar * t_bad).sum() / a_scalar.sum()
    check("vector-norm weights break the prescribed mean by the spec's amount",
          abs(phys - (10.0 - 0.230508652763321)) < 1e-12, f"physical mean {phys:.15f}")
    check("the two weightings genuinely disagree", abs(m_true - m_bad) > 1e-3)

    print("gate 9: the GEOMETRY builder, from triangle coordinates")
    # Replicates perNodeAreaVec: per triangle A_f = 0.5*cross(B-A, C-A), then each vertex takes
    # A_f/3 (vector) and |A_f|/3 (scalar). Gate 8 supplied weights directly and so could not see
    # this; here the weights are BUILT.
    def build(tris):
        vec, sca = {}, {}
        for (A, B, C) in tris:
            A, B, C = map(np.asarray, (A, B, C))
            Af = 0.5 * np.cross(B - A, C - A)
            for v in (A, B, C):
                k = tuple(np.round(v, 12))
                vec[k] = vec.get(k, np.zeros(3)) + Af / 3.0
                sca[k] = sca.get(k, 0.0) + np.linalg.norm(Af) / 3.0
        return vec, sca

    # Two PERPENDICULAR unit-area triangles sharing the origin.
    t1 = ((0, 0, 0), (2, 0, 0), (0, 1, 0))     # A=(0,0,1),  |A|=1
    t2 = ((0, 0, 0), (2, 0, 0), (0, 0, 1))     # A=(0,-1,0), |A|=1
    vec, sca = build([t1, t2])
    o = (0.0, 0.0, 0.0)
    check("builder: shared node scalar area is 2/3", abs(sca[o] - 2.0 / 3.0) < 1e-12,
          f"got {sca[o]:.12f}")
    check("builder: |summed vector| is sqrt(2)/3, NOT 2/3",
          abs(np.linalg.norm(vec[o]) - np.sqrt(2.0) / 3.0) < 1e-12,
          f"got {np.linalg.norm(vec[o]):.12f}")
    check("builder: the two disagree by the bend factor sqrt(2)/2",
          abs(np.linalg.norm(vec[o]) / sca[o] - np.sqrt(2.0) / 2.0) < 1e-12)

    # Opposed normals: the vector sum CANCELS on a node that has real area.
    t3 = ((0, 0, 0), (0, 1, 0), (2, 0, 0))     # t1 with reversed winding -> A=(0,0,-1)
    vec2, sca2 = build([t1, t3])
    check("cancelling normals: scalar area is still 2/3", abs(sca2[o] - 2.0 / 3.0) < 1e-12)
    check("cancelling normals: |summed vector| is 0 -- the node would be DROPPED",
          np.linalg.norm(vec2[o]) < 1e-14, f"got {np.linalg.norm(vec2[o]):.3e}")

    # And the patch-level gate: |sum_f A_f| ~ 0 while sum_f |A_f| = 2.
    patch_vec = np.zeros(3)
    patch_sca = 0.0
    for (A, B, C) in (t1, t3):
        A, B, C = map(np.asarray, (A, B, C))
        Af = 0.5 * np.cross(B - A, C - A)
        patch_vec += Af
        patch_sca += np.linalg.norm(Af)
    check("patch: |sum A_f| ~ 0 would skip setup entirely",
          np.linalg.norm(patch_vec) < 1e-14)
    check("patch: sum |A_f| = 2 is the real area to gate on",
          abs(patch_sca - 2.0) < 1e-12)

    print()
    if FAIL:
        print(f"{len(FAIL)} FAILED: {', '.join(FAIL)}")
        return 1
    print("all host gates passed; CUDA/MPI gates still pending on Alps")
    return 0


if __name__ == "__main__":
    sys.exit(main())
