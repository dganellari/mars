#!/usr/bin/env python3
"""Acceptance gates for the outlet boundary flux sample.

Fixtures from docs/design/gpt_outlet_boundary_spec_2026-09-09.md section 6. Replicates

    q_{f,r} = u_r . (A_f/3) + D_f (gbar_f - grad p_mix) . (A_f/3)
    grad p_mix = p_o grad N_o + sum_r t_r grad N_r
    D_f        = mean of the THREE FACE-node coefficients      <- no opposing node

This is a host replica of the algebra in `boundaryMassFluxKernel` /
`boundaryFaceCoefficient`; running the actual CUDA evaluator is still pending, because
mars_cvfem_kernel.hpp pulls in cuda_runtime.h and will not build on a CPU-only host.
"""

import sys
import numpy as np

FAIL = []
# Unit tetrahedron; outlet face is (1,2,3), opposite node 0.
GRAD_N = np.array([[-1.0, -1.0, -1.0], [1.0, 0, 0], [0, 1.0, 0], [0, 0, 1.0]])
AREA = np.array([0.5, 0.5, 0.5])      # outward area vector of face (1,2,3)


def check(name, ok, detail=""):
    print(f"  {'PASS' if ok else 'FAIL'}  {name}{('  ' + detail) if detail else ''}")
    if not ok:
        FAIL.append(name)


def face_coeff(d):
    """boundaryFaceCoefficient: face nodes only."""
    return sum(d) / 3.0


def triangle_flux(u_face, gbar, p_opp, trace, d_face, coeff=face_coeff):
    """Total volumetric flux over the triangle: the three samples summed."""
    grad_mix = p_opp * GRAD_N[0] + sum(trace[r] * GRAD_N[r + 1] for r in range(3))
    Df = coeff(d_face)
    total = 0.0
    for r in range(3):
        total += np.dot(u_face[r], AREA / 3.0)
        total += Df * np.dot(gbar - grad_mix, AREA / 3.0)
    return total


def main():
    zero3 = [np.zeros(3)] * 3

    print("gate 1: the coefficient is face-only, invariant to the opposite node")
    for d_opp in (2.0, 20.0):
        q = triangle_flux(zero3, np.zeros(3), 1.0, [0.0, 0.0, 0.0], (2.0, 4.0, 6.0))
        check(f"d_opp={d_opp:g}: total flux is 6", abs(q - 6.0) < 1e-12, f"got {q:.12f}")
    # And the blend it replaces really did depend on it, so this gate has teeth.
    blend = lambda d, o: 0.5 * (sum(d) / 3.0 + o)
    for d_opp, want in ((2.0, 4.5), (20.0, 18.0)):
        qb = triangle_flux(zero3, np.zeros(3), 1.0, [0.0] * 3, (2.0, 4.0, 6.0),
                           coeff=lambda d, o=d_opp: blend(d, o))
        check(f"the OLD face/opposite blend gave {want} at d_opp={d_opp:g}",
              abs(qb - want) < 1e-12, f"got {qb:.12f}")

    print("gate 2: affine pressure reproduced by the mixed element gradient")
    # p(x) = 7 + 2x - 3y + 5z  ->  nodal (7, 9, 4, 12); feed the manufactured trace directly.
    p_nodal = [7.0, 9.0, 4.0, 12.0]
    trace = p_nodal[1:]
    gbar = np.array([2.0, -3.0, 5.0])
    grad_mix = p_nodal[0] * GRAD_N[0] + sum(trace[r] * GRAD_N[r + 1] for r in range(3))
    check("mixed element gradient equals the analytic (2,-3,5)",
          np.allclose(grad_mix, gbar, atol=1e-12), f"got {grad_mix}")
    u = [np.array([1.0, 2.0, 3.0])] * 3
    q = triangle_flux(u, gbar, p_nodal[0], trace, (2.0, 2.0, 2.0))
    check("uniform velocity (1,2,3) gives advective flux 3 with zero RC difference",
          abs(q - 3.0) < 1e-12, f"got {q:.12f}")

    print("gate 3: sensitivities, D=2, delta=0.2")
    D, delta = (2.0, 2.0, 2.0), 0.2
    base = triangle_flux(zero3, np.zeros(3), 0.0, [0.0] * 3, D)
    dq_opp = triangle_flux(zero3, np.zeros(3), delta, [0.0] * 3, D) - base
    check("opposite-node pressure perturbation gives +1.5*D*delta = +0.6",
          abs(dq_opp - 0.6) < 1e-12, f"got {dq_opp:.12f}")
    dq_tr = triangle_flux(zero3, np.zeros(3), 0.0, [delta] * 3, D) - base
    check("common trace perturbation gives -1.5*D*delta = -0.6",
          abs(dq_tr + 0.6) < 1e-12, f"got {dq_tr:.12f}")

    print()
    if FAIL:
        print(f"{len(FAIL)} FAILED: {', '.join(FAIL)}")
        return 1
    print("all host gates passed; the CUDA evaluator gate is still pending")
    return 0


if __name__ == "__main__":
    sys.exit(main())
