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

// Node-driven assembly of the EXACT conjugate operator A_fem = D_fem M^-1 D_fem^T,
// the assembled form of the weak div/grad projection pair above. The K-path
// (--pressure-k) corrector applies M^-1 D_fem^T phi (computeFemGradientTetKernel
// + normalize-by-massNode); for the projection D_fem u^{n+1} = 0 to hold exactly
// the SOLVED operator must be D_fem M^-1 D_fem^T, NOT the Galerkin stiffness K
// (K != D_fem M^-1 D_fem^T -- they differ by a per-node valence factor on
// non-uniform meshes), so solving K leaves a residual divergence the corrector
// cannot remove. This kernel assembles A_fem.
//
// D_fem entry (element e, node a): a_a^e = -(V_e/4) * dNdx_a^e (a 3-vector,
//   independent of the velocity node it multiplies -- (D_fem u)_a = a_a^e . sum_j u_j).
// A_fem = D_fem M^-1 D_fem^T, M = lumped mass diag (massNode = sum_e V_e/4):
//   A_fem[a,b] = sum_i (1/M_i) * S_a^(i) . S_b^(i),
//   S_p^(i) = sum_{e containing both i and p} a_p^e.
// So at intermediate node i (the M^-1 node), gather the 1-ring of i, accumulate
// per-partner S_p, then scatter the OUTER PRODUCT (1/M_i) S_a.S_b over all pairs
// (a,b) in {i} U 1-ring(i). This is a pure sum of squares -> SPD, Laplacian-like,
// no SCS area-vector cancellation (the verified DDT artifact); the diagonal
// stays healthy (~ K/valence scale, NOT 1e-12).
//
// Sparsity REUSE: A_fem couples (a,b) iff both lie in some node i's 1-ring, the
// same two-hop-through-shared-node pattern buildDDTSparsityTet emits (every other
// tet corner is an SCS-edge neighbour of i), so the DDT tet sparsity is a correct
// superset -- atomicAddSparseEntry finds every column this kernel writes.
//
// Halo fold: one thread per node i over ALL nodes (owned + ghost incident
// elements via the node->element CSR). Writes gate on ownership[*]==1 + valid
// owned DOF, exactly like assembleDDTPerNodeKernelTet, so ghost rows are dropped
// and only owned rows are assembled (no separate reverse-halo needed -- the CSR
// is owned-row only and every owning rank visits its own intermediate nodes).
template<typename KeyType, typename RealType>
__global__ void assembleFemGramPerNodeKernelTet(
    const KeyType* c0, const KeyType* c1, const KeyType* c2, const KeyType* c3,
    const RealType* nodeX, const RealType* nodeY, const RealType* nodeZ,
    const KeyType* nodeToElemOffsets,
    const KeyType* nodeToElemList,
    const int* nodeToDof,
    const uint8_t* ownership,
    const RealType* lumpedMassNode,
    const int* rowPtr,
    const int* colInd,
    int numOwnedDofs,
    RealType* values,
    size_t numNodes)
{
    size_t i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i >= numNodes) return;
    if (ownership[i] != 1) return;
    int dofI = nodeToDof[i];
    if (dofI < 0 || dofI >= numOwnedDofs) return;

    RealType mi = lumpedMassNode[i];
    if (!(mi > RealType(0))) return;
    RealType invMi = RealType(1) / mi;

    constexpr int NPE = ElemTraits<TetTag>::NodesPerElem;   // 4

    // SCATTER-DIRECT, NO PARTNER BUFFER. A_fem[a,b] = (1/M_i) S_a.S_b with
    // S_p = sum_{e at i, p in e} a_p^e and a_p^e = -(V_e/4) dNdx_p^e. Expanding
    // the outer product of the two element-sums:
    //   (1/M_i) S_a.S_b = (1/M_i) sum_e sum_f (a_a^e . a_b^f),
    // a double sum over the elements e,f incident to i (a a corner of e, b a
    // corner of f). atomicAdd accumulates these element-pair contributions into
    // the SAME CSR entry (dofA,dofB), so the result is bit-for-bit the assembled
    // (1/M_i) S_a.S_b -- identical to the old buffered version, but WITHOUT a
    // fixed-size partner array. The CSR sparsity (built by buildDDTSparsityTet)
    // is the authoritative column set; atomicAddSparseEntry drops anything not
    // in it. So this kernel has NO per-node neighbour cap and CANNOT silently
    // overflow regardless of valence -- the only cap left is the sparsity one,
    // which is now counted+reported. Cost is O(m^2) element pairs per node (m =
    // incident elements), heavier than the buffered O(valence^2), but this runs
    // ONCE at setup, not per timestep, so it is amortized to nothing.
    KeyType eStart = nodeToElemOffsets[i];
    KeyType eEnd   = nodeToElemOffsets[i + 1];
    for (KeyType ep = eStart; ep < eEnd; ++ep)
    {
        KeyType e = nodeToElemList[ep];
        KeyType en[NPE] = {c0[e], c1[e], c2[e], c3[e]};
        bool eHasI = false;
        #pragma unroll
        for (int k = 0; k < NPE; ++k) if (en[k] == (KeyType)i) eHasI = true;
        if (!eHasI) continue;

        RealType coordsE[NPE][3];
        #pragma unroll
        for (int k = 0; k < NPE; ++k) {
            coordsE[k][0] = nodeX[en[k]];
            coordsE[k][1] = nodeY[en[k]];
            coordsE[k][2] = nodeZ[en[k]];
        }
        RealType detE, dNdxE[NPE][3];
        Tet4CVFEM::jacobian_and_dNdx<RealType>(coordsE, detE, dNdxE);
        RealType VE = detE / RealType(6);
        if (!(VE > RealType(0))) continue;
        RealType nqVE = -VE * RealType(0.25);   // -(V_e/4)

        // Inner loop over the second element f (also incident to i). a comes
        // from e (using dNdxE), b comes from f (using dNdxF). When f==e the two
        // geometries coincide -> reuse dNdxE, no recompute.
        for (KeyType fp = eStart; fp < eEnd; ++fp)
        {
            KeyType f = nodeToElemList[fp];
            KeyType fn[NPE] = {c0[f], c1[f], c2[f], c3[f]};
            bool fHasI = false;
            #pragma unroll
            for (int k = 0; k < NPE; ++k) if (fn[k] == (KeyType)i) fHasI = true;
            if (!fHasI) continue;

            RealType nqVF;
            const RealType (*dNdxF)[3];
            RealType dNdxFbuf[NPE][3];
            if (f == e)
            {
                nqVF  = nqVE;
                dNdxF = dNdxE;
            }
            else
            {
                RealType coordsF[NPE][3];
                #pragma unroll
                for (int k = 0; k < NPE; ++k) {
                    coordsF[k][0] = nodeX[fn[k]];
                    coordsF[k][1] = nodeY[fn[k]];
                    coordsF[k][2] = nodeZ[fn[k]];
                }
                RealType detF;
                Tet4CVFEM::jacobian_and_dNdx<RealType>(coordsF, detF, dNdxFbuf);
                RealType VF = detF / RealType(6);
                if (!(VF > RealType(0))) continue;
                nqVF  = -VF * RealType(0.25);   // -(V_f/4)
                dNdxF = dNdxFbuf;
            }

            // For each corner a of e (row) and corner b of f (col), scatter
            // (1/M_i) (a_a^e . a_b^f) into A_fem[dofA,dofB]. a_p = nqV * dNdx_p.
            #pragma unroll
            for (int ka = 0; ka < NPE; ++ka)
            {
                int dofA = nodeToDof[en[ka]];
                if (dofA < 0 || dofA >= numOwnedDofs) continue;  // owned rows only
                int rs = rowPtr[dofA], re = rowPtr[dofA + 1];
                RealType ax = nqVE * dNdxE[ka][0];
                RealType ay = nqVE * dNdxE[ka][1];
                RealType az = nqVE * dNdxE[ka][2];
                #pragma unroll
                for (int kb = 0; kb < NPE; ++kb)
                {
                    int dofB = nodeToDof[fn[kb]];
                    if (dofB < 0) continue;
                    RealType bx = nqVF * dNdxF[kb][0];
                    RealType by = nqVF * dNdxF[kb][1];
                    RealType bz = nqVF * dNdxF[kb][2];
                    RealType dot = ax * bx + ay * by + az * bz;
                    fem::atomicAddSparseEntry(values, colInd, rs, re, dofB, invMi * dot);
                }
            }
        }
    }
}
