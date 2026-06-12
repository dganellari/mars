#pragma once

// FEM-consistent (weak-form) divergence + adjoint gradient pair for the
// pressure projection. OPT-IN (s.useFemProjection, set by --pressure-k); off
// by default so the validated SCS behavior is byte-identical without the flag.
//
// WHY: the SCS divergence D (computeDivergencePerNodeKernel) and its corrector
// gradient M^-1 D^T are built from SIGNED area-vector sums that cancel on
// perfectly good tets, so D M^-1 D^T has near-zero modes the Galerkin K solve
// never sees. The per-step projection (I - D M^-1 D^T K^-1) then leaves those
// modes growing ~1.2x/step -> blowup by step ~15 even with a fully converged
// pressure solve. The weak-form pair below is a sum of squares: nothing
// cancels, ker(D_fem^T) = constants = ker(K), and D_fem M^-1 D_fem^T is
// spectrally equivalent to K -> the projection contracts divergence on ALL
// modes.
//
// Weak divergence (P1 tet, constant dNdx):
//   b_i = integral(N_i div u) = -integral(grad N_i . u) + surface(N_i u.n)
// Interior part per element e (u linear, integral_e u = (V/4) sum_j u_j):
//   b_i += -(V/4) * (dNdx_i . (u_0+u_1+u_2+u_3)),  V = det/6.
// The surface term at the inlet/outlet openings is the opening-flux source
// (--opening-flux-source). On THIS path it must be the CONSISTENT P1 face
// quadrature oint(N_i u.n dA) = sum_f (A_f/12)(2u_i+u_j+u_k).n_f -- the lumped
// per-node u_i.areaVec_i agrees in total but not node-by-node, and that fixed
// mismatch is a spurious mass source that blows up geometrically. The solver
// switches to the consistent weights (s.d_femFluxWin/Wout, built by the pump
// driver) when useFemProjection is on. Walls give 0 (u=0). Same positive-
// outflow divAccNode convention as the SCS scatter, so buildPressureRhs
// (rhs = -coef*divAcc) and everything downstream are untouched.
//
// Adjoint gradient (corrector): gradPhi_j = [sum_e (V/4)(grad phi)_e] / M_j,
// with M_j the existing lumped mass d_massNode (exactly sum_e V/4, the tet
// branch of computeLumpedMassPerNodeKernel). The accumulator is -D_fem^T, so
// M^-1 of it estimates +grad(phi): feed straight to applyCorrectorPerNodeKernel
// (q = q** - dt/rho * gradPhi), NO sign flip (unlike the divT path).

// IMPORTANT — INCLUSION ORDER:
// Like mars_vms_pressure_stab.hpp, this header is meant to be #include'd
// INSIDE mars_ns_pump_solver.hpp (or mars_ns_solver.hpp), AFTER its
// `using namespace mars; using namespace mars::fem;`. It deliberately lives in
// the GLOBAL namespace to match the solver's other kernels and to see the
// unqualified Tet4CVFEM / ElemTraits brought in by those usings. Do not
// include it standalone.

// Weak-form divergence volume term: one thread per element, atomicAdd scatter
// to all 4 nodes (ghosts included; the caller's reverseExchangeNodeHaloAdd
// folds ghost contributions to owners, exactly like the SCS divergence path).
template<typename KeyType, typename RealType>
__global__ void computeFemDivergenceTetKernel(
    const KeyType* c0, const KeyType* c1, const KeyType* c2, const KeyType* c3,
    const RealType* vx, const RealType* vy, const RealType* vz,
    const RealType* nodeX, const RealType* nodeY, const RealType* nodeZ,
    RealType* divAccNode,
    size_t startElem, size_t numLocal)
{
    size_t k = blockIdx.x * blockDim.x + threadIdx.x;
    if (k >= numLocal) return;
    size_t e = startElem + k;

    constexpr int NPE = ElemTraits<TetTag>::NodesPerElem;   // 4
    const KeyType* cc[4] = {c0, c1, c2, c3};
    KeyType n[NPE];
    for (int i = 0; i < NPE; ++i) n[i] = cc[i][e];

    RealType coords[4][3];
    for (int i = 0; i < NPE; ++i) {
        coords[i][0] = nodeX[n[i]];
        coords[i][1] = nodeY[n[i]];
        coords[i][2] = nodeZ[n[i]];
    }
    RealType det, dNdx[4][3];
    Tet4CVFEM::jacobian_and_dNdx<RealType>(coords, det, dNdx);
    RealType V = det / RealType(6);
    // Degenerate or inverted tet: dNdx ~ 1/det would poison the accumulator
    // with inf/NaN; skip rather than scatter garbage.
    if (!(V > RealType(0))) return;

    // integral_e u = (V/4) * (u_0+u_1+u_2+u_3), exact for linear u.
    RealType usx = 0, usy = 0, usz = 0;
    for (int j = 0; j < NPE; ++j) {
        usx += vx[n[j]];
        usy += vy[n[j]];
        usz += vz[n[j]];
    }

    RealType quarterV = V * RealType(0.25);
    for (int i = 0; i < NPE; ++i) {
        RealType gdotu = dNdx[i][0] * usx + dNdx[i][1] * usy + dNdx[i][2] * usz;
        atomicAdd(&divAccNode[n[i]], -quarterV * gdotu);
    }
}

// Adjoint gradient accumulator: scatters (V/4)*(grad phi)_e to each of the 4
// nodes. The caller normalizes by d_massNode and handles the reverse-halo +
// periodic folding, mirroring the SCS gradient pipeline.
template<typename KeyType, typename RealType>
__global__ void computeFemGradientTetKernel(
    const KeyType* c0, const KeyType* c1, const KeyType* c2, const KeyType* c3,
    const RealType* phi,
    const RealType* nodeX, const RealType* nodeY, const RealType* nodeZ,
    RealType* gxAcc, RealType* gyAcc, RealType* gzAcc,
    size_t startElem, size_t numLocal)
{
    size_t k = blockIdx.x * blockDim.x + threadIdx.x;
    if (k >= numLocal) return;
    size_t e = startElem + k;

    constexpr int NPE = ElemTraits<TetTag>::NodesPerElem;   // 4
    const KeyType* cc[4] = {c0, c1, c2, c3};
    KeyType n[NPE];
    for (int i = 0; i < NPE; ++i) n[i] = cc[i][e];

    RealType coords[4][3];
    for (int i = 0; i < NPE; ++i) {
        coords[i][0] = nodeX[n[i]];
        coords[i][1] = nodeY[n[i]];
        coords[i][2] = nodeZ[n[i]];
    }
    RealType det, dNdx[4][3];
    Tet4CVFEM::jacobian_and_dNdx<RealType>(coords, det, dNdx);
    RealType V = det / RealType(6);
    if (!(V > RealType(0))) return;

    // constant element gradient of phi (linear tet)
    RealType gx = 0, gy = 0, gz = 0;
    for (int i = 0; i < NPE; ++i) {
        RealType ph = phi[n[i]];
        gx += dNdx[i][0] * ph;
        gy += dNdx[i][1] * ph;
        gz += dNdx[i][2] * ph;
    }

    RealType quarterV = V * RealType(0.25);
    for (int j = 0; j < NPE; ++j) {
        atomicAdd(&gxAcc[n[j]], quarterV * gx);
        atomicAdd(&gyAcc[n[j]], quarterV * gy);
        atomicAdd(&gzAcc[n[j]], quarterV * gz);
    }
}
