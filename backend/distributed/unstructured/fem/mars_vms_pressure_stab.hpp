#pragma once

// Nalu-Wind-style VMS (variational multiscale) pressure stabilization for the
// CVFEM incompressible solver. SEPARATE, OPT-IN path -- does NOT touch the
// existing computeDivergenceRhieChowKernel or the DDT projection. Enable it
// explicitly (a flag in the driver); off by default, so the validated pump/TGV
// behavior is unchanged until this is proven and switched on deliberately.
//
// WHY this exists: the literal compact Rhie-Chow term tau*(p_R-p_L)*|A|^2/(A.dx)
// blows up on skewed tets (A.dx -> 0). Nalu-Wind avoids the A.dx quotient
// entirely: it adds a residual-based pressure stabilization built from the
// difference between the ELEMENT ip pressure gradient and the PROJECTED NODAL
// (Green-Gauss) gradient. That difference is a 4th-order pressure dissipation --
// it vanishes for smooth p (no accuracy loss) and damps the checkerboard mode.
//
// Stabilized SCS mass flux (Nalu theory manual, Domino CTR-2014):
//   mdot_ip = ( rho*uhat_i + tau*G_i(p) - tau*dp/dx_i|_ip ) * A_i
// where
//   uhat_i        provisional face velocity = 0.5*(u_L + u_R)
//   dp/dx_i|_ip   element ip gradient = sum_k dN_k/dx_i * p_k     (element basis)
//   G_i(p)        projected nodal gradient (Green-Gauss), at ip = 0.5*(G_L+G_R)
//   tau           stabilization time scale ~ dt
// The stabilization is the bracket tau*(G_i(p) - dp/dx_i)*A_i. NO A.dx anywhere.
//
// This header provides ONLY the tet kernel (the pump is tets). Hex is separate
// later work (hex needs per-ip shape derivatives).

// IMPORTANT — INCLUSION ORDER:
// This header is meant to be #include'd INSIDE mars_ns_pump_solver.hpp (or
// mars_ns_solver.hpp), AFTER its `using namespace mars; using namespace
// mars::fem;` and AFTER its file-local `scsLR<>` definition (~line 591). It
// deliberately lives in the GLOBAL namespace to match the solver's other
// kernels and to see the file-local scsLR<> and the unqualified
// Tet4CVFEM / ElemTraits brought in by those usings. Do not include it
// standalone.

// Tet VMS-stabilized divergence: accumulates the stabilized SCS mass flux into
// per-node divergence. Mirrors computeDivergenceRhieChowKernel's I/O so it can
// be dropped into the same call site, but uses the Nalu residual instead of the
// compact A.dx term.
//
//   gradPx/y/z : PROJECTED NODAL pressure gradient (Green-Gauss), per node,
//                MUST be halo-complete (computeGradientPerNodeKernel +
//                normalizeGradientPerNodeKernel on p, then halo exchange).
//   tau        : stabilization time scale (~dt). Pass the same value the
//                existing path uses for rhieChowTau.
template<typename KeyType, typename RealType>
__global__ void computeDivergenceVMSTetKernel(
    const KeyType* c0, const KeyType* c1, const KeyType* c2, const KeyType* c3,
    const RealType* vx, const RealType* vy, const RealType* vz,
    const RealType* p,
    const RealType* gradPx, const RealType* gradPy, const RealType* gradPz,
    const RealType* nodeX, const RealType* nodeY, const RealType* nodeZ,
    const RealType* areaVecX, const RealType* areaVecY, const RealType* areaVecZ,
    RealType tau,
    bool keepSmooth,          // true=full Nalu G-dpdx (blows up on Chorin u**); false=compact -dpdx only (default)
    RealType* divAccNode,
    size_t startElem, size_t numLocal)
{
    size_t k = blockIdx.x * blockDim.x + threadIdx.x;
    if (k >= numLocal) return;
    size_t e = startElem + k;

    constexpr int NPE  = ElemTraits<TetTag>::NodesPerElem;   // 4
    constexpr int NSCS = ElemTraits<TetTag>::ScsPerElem;     // 6

    KeyType n[NPE];
    {
        const KeyType* cc[4] = {c0, c1, c2, c3};
        for (int i = 0; i < NPE; ++i) n[i] = cc[i][e];
    }

    // Element ip pressure gradient dp/dx_i = sum_k dN_k/dx_i * p_k.
    // For a linear tet, dN/dx is CONSTANT over the element, so dp/dx is the same
    // for all 6 SCS -> compute it ONCE per element (hoisted out of the ip loop).
    RealType coords[4][3];
    for (int i = 0; i < NPE; ++i) {
        coords[i][0] = nodeX[n[i]];
        coords[i][1] = nodeY[n[i]];
        coords[i][2] = nodeZ[n[i]];
    }
    RealType det;
    RealType dNdx[4][3];
    Tet4CVFEM::jacobian_and_dNdx<RealType>(coords, det, dNdx);

    RealType dpdx = 0, dpdy = 0, dpdz = 0;
    for (int kk = 0; kk < NPE; ++kk) {
        RealType pk = p[n[kk]];
        dpdx += dNdx[kk][0] * pk;
        dpdy += dNdx[kk][1] * pk;
        dpdz += dNdx[kk][2] * pk;
    }

    #pragma unroll
    for (int ip = 0; ip < NSCS; ++ip)
    {
        int nodeL, nodeR; scsLR<TetTag>(ip, nodeL, nodeR);
        KeyType iL = n[nodeL];
        KeyType iR = n[nodeR];

        // Provisional (post-predictor) face velocity.
        RealType vfx = RealType(0.5) * (vx[iL] + vx[iR]);
        RealType vfy = RealType(0.5) * (vy[iL] + vy[iR]);
        RealType vfz = RealType(0.5) * (vz[iL] + vz[iR]);

        size_t off = e * NSCS + ip;
        RealType Ax = areaVecX[off];
        RealType Ay = areaVecY[off];
        RealType Az = areaVecZ[off];

        // base advective flux
        RealType flow = vfx * Ax + vfy * Ay + vfz * Az;

        // Projected nodal gradient interpolated to the ip = 0.5*(G_L + G_R).
        RealType Gx = RealType(0.5) * (gradPx[iL] + gradPx[iR]);
        RealType Gy = RealType(0.5) * (gradPy[iL] + gradPy[iR]);
        RealType Gz = RealType(0.5) * (gradPz[iL] + gradPz[iR]);

        // Stabilization flux. The FULL Nalu term is tau*(G - dp/dx).A. But that
        // KEEPS the smooth +tau*G.A half, and on a Chorin POST-PREDICTOR u**
        // (which already subtracted -dt/rho*grad p^n) the smooth half DOUBLE-COUNTS
        // grad p -> pressure ramps unbounded (validated empirically on cube16: |p|
        // climbs 44->626+ over 200 steps, same failure the RC kernel documents).
        //
        // FIX (keepSmooth=false, default): drop the smooth half, keep only the
        // compact element-ip-gradient term tau*(-dp/dx).A. This is the
        // checkerboard-mode amplitude (the ip gradient sees the sign-alternating p
        // strongly; subtracting it damps the mode) WITHOUT re-adding the smooth
        // grad p the predictor already applied. Mirrors what Rhie-Chow did for Chorin.
        RealType stab;
        if (keepSmooth)
            stab = tau * ((Gx - dpdx) * Ax + (Gy - dpdy) * Ay + (Gz - dpdz) * Az);
        else
            stab = tau * ((-dpdx) * Ax + (-dpdy) * Ay + (-dpdz) * Az);
        flow += stab;

        atomicAdd(&divAccNode[iL], +flow);
        atomicAdd(&divAccNode[iR], -flow);
    }
}
