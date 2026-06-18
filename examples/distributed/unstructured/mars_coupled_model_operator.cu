// Coupled collocated solver — Phase 0 model operator (see
// internal-notes/plan_coupled_collocated_solver.md).
//
// This standalone driver builds, validates, and (eventually) directly solves a
// small coupled [u,v,w,p] block operator for the paper benchmarks (lid-driven
// skew cavity, backward-facing step), as the throwaway oracle that de-risks the
// Phase-1 preconditioner bake-off.
//
// INCREMENT 1 (this file, now): the foundational geometry gate only.
//   Verify G = -D^T (discrete integration by parts) on the weak FEM div/grad
//   operators the coupled off-diagonal blocks are built from. If the geometry /
//   connectivity is sound, <phi, D v> = -<v, G phi> holds to round-off for ANY
//   fields phi, v -- so a failure here is a mesh/connectivity bug, caught before
//   any coupled physics is written.
// NEXT INCREMENTS (not yet here): momentum block a^uu (advection+diffusion),
//   the Rhie-Chow continuity coefficient a^pp/a^pu, block-CSR assembly, the
//   direct-solve oracle, and the published-profile comparison.

#include "backend/distributed/unstructured/fem/mars_ns_solver.hpp"  // pulls domain, AMR, CVFEM, ElemTraits, SparseMatrix

using namespace mars;
using namespace mars::fem;
using namespace mars::amr;

// Must follow the usings: the projection kernels live in the global namespace
// and resolve Tet4CVFEM / ElemTraits unqualified (header is not standalone).
#include "backend/distributed/unstructured/fem/mars_fem_projection.hpp"
#include "backend/distributed/unstructured/solvers/mars_acm_coarsen.hpp"  // ACM Stage 3 coarse-op + transfers
#include "backend/distributed/unstructured/solvers/mars_acm_vcycle.hpp"   // ACM Stage 4a V-cycle
#include "backend/distributed/unstructured/solvers/mars_gpu_acm_preconditioner.hpp"  // ACM Stage 4b precond
#include "backend/distributed/unstructured/solvers/mars_gmres_solver.hpp"            // native (Flex)GMRES

#include <thrust/device_vector.h>
#include <thrust/inner_product.h>
#include <thrust/extrema.h>
#include <thrust/transform.h>
#include <thrust/reduce.h>
#include <thrust/functional.h>
#include <thrust/fill.h>
#include <thrust/copy.h>
#include <cusolverSp.h>
#include <cusparse.h>
#include <mpi.h>
#include <iostream>
#include <iomanip>
#include <string>
#include <vector>
#include <set>
#include <map>
#include <utility>
#include <cmath>
#include <cstdlib>

// ---- Increment 2a: coupled Stokes block assembly (advection OFF) ----
// Assembles the interleaved (dof = 4*node + comp) CSR with the coefficients
// host-replica-validated in scripts/coupled_assembly_replica.py:
//   a^up = +(V/4) dNdx[b][d]   (momentum row, pressure col)
//   a^pu = -(V/4) dNdx[a][d]   (continuity row, velocity col) = -(a^up)^T
//   a^uu = nu*K (component-diagonal)   a^pp = tau*K   with K_ab = V*(dNdx[a].dNdx[b])
// One thread per element scatters its 4-node x 4-comp contributions into the CSR.
template<typename RealType>
__device__ inline void csrAdd(RealType* vals, const int* col, int rs, int re, int c, RealType v)
{
    for (int k = rs; k < re; ++k)
        if (col[k] == c) { atomicAdd(&vals[k], v); return; }
}

template<typename KeyType, typename RealType>
__global__ void assembleCoupledStokesKernel(
    const KeyType* c0, const KeyType* c1, const KeyType* c2, const KeyType* c3,
    const RealType* nodeX, const RealType* nodeY, const RealType* nodeZ,
    const int* rowOff, const int* colInd, RealType* vals,
    RealType nu, RealType tau,
    const RealType* avx, const RealType* avy, const RealType* avz,  // frozen advecting vel; null = OFF
    size_t startElem, size_t numLocal)
{
    size_t k = blockIdx.x * blockDim.x + threadIdx.x;
    if (k >= numLocal) return;
    size_t e = startElem + k;

    const KeyType* cc[4] = {c0, c1, c2, c3};
    int n[4];
    for (int i = 0; i < 4; ++i) n[i] = static_cast<int>(cc[i][e]);

    RealType coords[4][3];
    for (int i = 0; i < 4; ++i) {
        coords[i][0] = nodeX[n[i]];
        coords[i][1] = nodeY[n[i]];
        coords[i][2] = nodeZ[n[i]];
    }
    RealType det, dNdx[4][3];
    Tet4CVFEM::jacobian_and_dNdx<RealType>(coords, det, dNdx);
    RealType V = det / RealType(6);
    if (!(V > RealType(0))) return;
    RealType qV = V * RealType(0.25);

    // frozen element-mean advecting velocity (consistent-Galerkin convection)
    RealType ax = 0, ay = 0, az = 0;
    if (avx) {
        for (int i = 0; i < 4; ++i) { ax += avx[n[i]]; ay += avy[n[i]]; az += avz[n[i]]; }
        ax *= RealType(0.25); ay *= RealType(0.25); az *= RealType(0.25);
    }

    for (int a = 0; a < 4; ++a) {
        int ra3 = 4 * n[a] + 3, rs3 = rowOff[ra3], re3 = rowOff[ra3 + 1];
        for (int b = 0; b < 4; ++b) {
            RealType Kab = V * (dNdx[a][0] * dNdx[b][0]
                              + dNdx[a][1] * dNdx[b][1]
                              + dNdx[a][2] * dNdx[b][2]);
            RealType Cab = avx ? qV * (ax * dNdx[b][0] + ay * dNdx[b][1] + az * dNdx[b][2])
                               : RealType(0);                          // C_ij = (V/4)(a.gradN_j)
            csrAdd<RealType>(vals, colInd, rs3, re3, 4 * n[b] + 3, tau * Kab);   // a^pp
            for (int d = 0; d < 3; ++d)
                csrAdd<RealType>(vals, colInd, rs3, re3, 4 * n[b] + d, -qV * dNdx[a][d]); // a^pu
            for (int d = 0; d < 3; ++d) {
                int rad = 4 * n[a] + d, rs = rowOff[rad], re = rowOff[rad + 1];
                csrAdd<RealType>(vals, colInd, rs, re, 4 * n[b] + d, nu * Kab + Cab);     // a^uu diff+adv
                csrAdd<RealType>(vals, colInd, rs, re, 4 * n[b] + 3, qV * dNdx[b][d]);    // a^up
            }
        }
    }
}

// ---- Increment 3a: Dirichlet BC elimination + direct reference solve ----
// Identity-row BC: for a Dirichlet dof, zero its CSR row, set its diagonal to 1,
// rhs = prescribed value. Correct for a DIRECT solve (the interior rows still see
// the pinned value through their column coupling); no column elimination needed.
template<typename RealType>
__global__ void applyBCKernel(const uint8_t* bcFlag, const RealType* bcVal,
                              const int* rowOff, const int* colInd,
                              RealType* vals, RealType* rhs, int ND)
{
    int r = blockIdx.x * blockDim.x + threadIdx.x;
    if (r >= ND || !bcFlag[r]) return;
    for (int k = rowOff[r]; k < rowOff[r + 1]; ++k)
        vals[k] = (colInd[k] == r) ? RealType(1) : RealType(0);
    rhs[r] = bcVal[r];
}

template<typename RealType>
__global__ void extractCompKernel(const RealType* x, RealType* out, int stride, int comp, int N)
{
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i >= N) return;
    out[i] = x[stride * i + comp];
}

int main(int argc, char** argv)
{
    MPI_Init(&argc, &argv);
    int rank = 0, numRanks = 1;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &numRanks);

    // TET AmrManager init needs a bound device (see memory: tet drivers must
    // set the device themselves, unlike the hex tgv driver).
    int deviceCount = 0;
    cudaGetDeviceCount(&deviceCount);
    if (deviceCount > 0) cudaSetDevice(rank % deviceCount);

    using KeyType  = uint64_t;
    using RealType = double;

    std::string meshFile;
    int blockSize  = 256;
    int bucketSize = 64;
    bool doAssemble = false;     // --assemble: also run the coupled Stokes block assembly + gates
    bool doAdvect   = false;     // --advect: add frozen-velocity convection into a^uu (implies --assemble)
    bool doSolve    = false;     // --solve: BC elimination + cusolverSp QR direct reference solve (implies --assemble)
    bool doHypre    = false;     // --hypre: Phase-1 1a -- Hypre point-block BoomerAMG GMRES iters vs cusolver ref (implies --assemble)
    bool doAcm3     = false;      // --acm-stage3: GPU coarse-operator + transfers vs host P^T A P (implies --assemble)
    bool doAcm4     = false;      // --acm-stage4: GPU multilevel V-cycle vs host V-cycle replica (implies --assemble)
    bool doAcm4b    = false;      // --acm-stage4b: ACM-preconditioned FlexGMRES iters vs 1a baseline + cusolver ref
    double beta     = 45.0;      // --beta: skew-cavity angle (deg), used only for the --solve BC geometry
    int    Re       = 100;       // --Re: Reynolds number for the lid-driven solve (3b)
    int    picard   = 30;        // --picard: max Picard outer iterations
    for (int i = 1; i < argc; ++i)
    {
        std::string a = argv[i];
        if      (a.rfind("--mesh=", 0) == 0)        meshFile   = a.substr(7);
        else if (a.rfind("--block-size=", 0) == 0)  blockSize  = std::stoi(a.substr(13));
        else if (a.rfind("--bucket-size=", 0) == 0) bucketSize = std::stoi(a.substr(14));
        else if (a == "--assemble")                 doAssemble = true;
        else if (a == "--advect")                 { doAssemble = true; doAdvect = true; }
        else if (a == "--solve")                  { doAssemble = true; doSolve  = true; }
        else if (a == "--hypre")                  { doAssemble = true; doHypre  = true; }
        else if (a == "--acm-stage3")             { doAssemble = true; doAcm3   = true; }
        else if (a == "--acm-stage4")             { doAssemble = true; doAcm4   = true; }
        else if (a == "--acm-stage4b")            { doAssemble = true; doAcm4b  = true; }
        else if (a.rfind("--beta=", 0) == 0)        beta       = std::stod(a.substr(7));
        else if (a.rfind("--Re=", 0) == 0)          Re         = std::stoi(a.substr(5));
        else if (a.rfind("--picard=", 0) == 0)      picard     = std::stoi(a.substr(9));
        else if (a[0] != '-' && meshFile.empty())   meshFile   = a;
    }
    if (meshFile.empty())
    {
        if (rank == 0)
            std::cout << "Usage: " << argv[0] << " --mesh=FILE.exo [--block-size=N]\n";
        MPI_Finalize();
        return 1;
    }

    // Increment 1 is a single-rank geometry check: the div/grad kernels scatter
    // to ghosts too, so a correct multi-rank check needs reverseExchangeNodeHaloAdd
    // before the dot products. That fold is a later increment; refuse >1 rank now
    // rather than report a wrong (ghost-double-counted) residual.
    if (numRanks > 1)
    {
        if (rank == 0)
            std::cout << "[phase0] increment-1 adjointness gate is single-rank only "
                         "(multi-rank halo-fold is a later increment). Run with 1 rank.\n";
        MPI_Finalize();
        return 1;
    }

    AmrManager<TetTag, KeyType, RealType>::Config cfg;
    cfg.maxLevels  = 0;          // frozen mesh
    cfg.blockSize  = blockSize;
    cfg.bucketSize = bucketSize;
    AmrManager<TetTag, KeyType, RealType> amr(cfg);
    amr.initialize(meshFile, rank, numRanks);
    auto& domain = amr.domain();

    const size_t nNodes   = domain.getNodeCount();
    const size_t startEl  = domain.startIndex();
    const size_t numLocal = domain.localElementCount();
    if (rank == 0)
        std::cout << "[phase0] mesh: elements=" << domain.getElementCount()
                  << " nodes=" << nNodes << " localElems=" << numLocal << "\n";

    const auto& conn = domain.getElementToNodeConnectivity();
    const KeyType* c0 = std::get<0>(conn).data();
    const KeyType* c1 = std::get<1>(conn).data();
    const KeyType* c2 = std::get<2>(conn).data();
    const KeyType* c3 = std::get<3>(conn).data();
    const RealType* nx = domain.getNodeX().data();
    const RealType* ny = domain.getNodeY().data();
    const RealType* nz = domain.getNodeZ().data();

    // Deterministic smooth test fields from coordinates. Any phi, v validate the
    // adjointness identity; smooth polynomials make a failure easy to read.
    std::vector<RealType> hx(nNodes), hy(nNodes), hz(nNodes);
    thrust::copy(thrust::device_pointer_cast(nx), thrust::device_pointer_cast(nx + nNodes), hx.begin());
    thrust::copy(thrust::device_pointer_cast(ny), thrust::device_pointer_cast(ny + nNodes), hy.begin());
    thrust::copy(thrust::device_pointer_cast(nz), thrust::device_pointer_cast(nz + nNodes), hz.begin());

    std::vector<RealType> hphi(nNodes), hvx(nNodes), hvy(nNodes), hvz(nNodes);
    for (size_t i = 0; i < nNodes; ++i)
    {
        const RealType x = hx[i], y = hy[i], z = hz[i];
        hphi[i] = 1.0 + 2.0 * x - 3.0 * y + 0.7 * z + 0.4 * x * y - 0.2 * y * z;
        hvx[i]  = 0.5 + x - 0.3 * z + 0.1 * x * x;
        hvy[i]  = -0.2 + 0.4 * y + 0.6 * z;
        hvz[i]  = 0.9 - 0.5 * x + 0.3 * y * y;
    }

    thrust::device_vector<RealType> phi(hphi.begin(), hphi.end());
    thrust::device_vector<RealType> vx(hvx.begin(), hvx.end());
    thrust::device_vector<RealType> vy(hvy.begin(), hvy.end());
    thrust::device_vector<RealType> vz(hvz.begin(), hvz.end());
    thrust::device_vector<RealType> divAcc(nNodes, RealType(0));
    thrust::device_vector<RealType> gx(nNodes, RealType(0));
    thrust::device_vector<RealType> gy(nNodes, RealType(0));
    thrust::device_vector<RealType> gz(nNodes, RealType(0));

    const int eBlocks = static_cast<int>((numLocal + blockSize - 1) / blockSize);

    // D v  (weak divergence, scattered to nodes)
    computeFemDivergenceTetKernel<KeyType, RealType><<<eBlocks, blockSize>>>(
        c0, c1, c2, c3,
        thrust::raw_pointer_cast(vx.data()),
        thrust::raw_pointer_cast(vy.data()),
        thrust::raw_pointer_cast(vz.data()),
        nx, ny, nz,
        thrust::raw_pointer_cast(divAcc.data()),
        startEl, numLocal);

    // G phi  (un-normalized weak gradient accumulator -- this is exactly D^T's
    // partner; we do NOT divide by the lumped mass here, since the adjoint
    // identity holds between the raw operators).
    computeFemGradientTetKernel<KeyType, RealType><<<eBlocks, blockSize>>>(
        c0, c1, c2, c3,
        thrust::raw_pointer_cast(phi.data()),
        nx, ny, nz,
        thrust::raw_pointer_cast(gx.data()),
        thrust::raw_pointer_cast(gy.data()),
        thrust::raw_pointer_cast(gz.data()),
        startEl, numLocal);

    cudaError_t err = cudaDeviceSynchronize();
    if (err != cudaSuccess)
    {
        std::cerr << "[phase0] kernel launch failed: " << cudaGetErrorString(err) << "\n";
        MPI_Finalize();
        return 2;
    }

    // <phi, D v> and <v, G phi>; the discrete integration-by-parts identity is
    // <phi, D v> + <v, G phi> = 0.
    const RealType dDiv = thrust::inner_product(phi.begin(), phi.end(), divAcc.begin(), RealType(0));
    const RealType dGx  = thrust::inner_product(vx.begin(),  vx.end(),  gx.begin(),     RealType(0));
    const RealType dGy  = thrust::inner_product(vy.begin(),  vy.end(),  gy.begin(),     RealType(0));
    const RealType dGz  = thrust::inner_product(vz.begin(),  vz.end(),  gz.begin(),     RealType(0));
    const RealType dGrad = dGx + dGy + dGz;

    const RealType residual = dDiv + dGrad;
    const RealType scale    = std::max({std::abs(dDiv), std::abs(dGrad), RealType(1e-300)});
    const RealType relErr   = std::abs(residual) / scale;

    const RealType tol = 1e-12;
    const bool pass = relErr < tol;

    std::cout << std::scientific << std::setprecision(6);
    std::cout << "[phase0][G=-D^T gate] <phi,Dv> = " << dDiv
              << "  <v,Gphi> = " << dGrad
              << "  rel|sum| = " << relErr
              << "  -> " << (pass ? "PASS" : "FAIL") << "\n";

    // ---- Increment 2a: assemble the coupled Stokes block + matrix-consistency gates ----
    bool assemblePass = true;
    if (doAssemble)
    {
        const RealType nu = 1e-3, tau = 1.0 / 24.0;   // values irrelevant to the gates
        const int ND = 4 * static_cast<int>(nNodes);

        // connectivity to host (single rank: numLocal == elementCount, startEl == 0)
        std::vector<KeyType> h0(numLocal), h1(numLocal), h2(numLocal), h3(numLocal);
        thrust::copy(thrust::device_pointer_cast(c0), thrust::device_pointer_cast(c0 + numLocal), h0.begin());
        thrust::copy(thrust::device_pointer_cast(c1), thrust::device_pointer_cast(c1 + numLocal), h1.begin());
        thrust::copy(thrust::device_pointer_cast(c2), thrust::device_pointer_cast(c2 + numLocal), h2.begin());
        thrust::copy(thrust::device_pointer_cast(c3), thrust::device_pointer_cast(c3 + numLocal), h3.begin());

        // node 1-ring adjacency (nodes sharing a tet), then interleaved 4N block sparsity
        std::vector<std::set<int>> ring(nNodes);
        for (size_t e = 0; e < numLocal; ++e) {
            int nn[4] = {(int)h0[e], (int)h1[e], (int)h2[e], (int)h3[e]};
            for (int a = 0; a < 4; ++a)
                for (int b = 0; b < 4; ++b) ring[nn[a]].insert(nn[b]);
        }
        std::vector<int> h_rowOff(ND + 1, 0);
        for (size_t i = 0; i < nNodes; ++i) {
            int deg = 4 * static_cast<int>(ring[i].size());
            for (int a = 0; a < 4; ++a) h_rowOff[4 * i + a + 1] = deg;
        }
        for (int r = 0; r < ND; ++r) h_rowOff[r + 1] += h_rowOff[r];
        const int nnz = h_rowOff[ND];
        std::vector<int> h_colInd(nnz);
        for (size_t i = 0; i < nNodes; ++i) {
            std::vector<int> rs(ring[i].begin(), ring[i].end());   // sorted (std::set)
            for (int a = 0; a < 4; ++a) {
                int p = h_rowOff[4 * i + a];
                for (int j : rs)
                    for (int b = 0; b < 4; ++b) h_colInd[p++] = 4 * j + b;
            }
        }

        thrust::device_vector<int>      d_rowOff(h_rowOff.begin(), h_rowOff.end());
        thrust::device_vector<int>      d_colInd(h_colInd.begin(), h_colInd.end());
        thrust::device_vector<RealType> d_vals(nnz, RealType(0));

        // frozen advecting velocity (divergence-free rotation about the mesh centre)
        thrust::device_vector<RealType> d_avx, d_avy, d_avz;
        const RealType *avxp = nullptr, *avyp = nullptr, *avzp = nullptr;
        if (doAdvect) {
            RealType xc = 0, yc = 0;
            for (size_t i = 0; i < nNodes; ++i) { xc += hx[i]; yc += hy[i]; }
            xc /= nNodes; yc /= nNodes;
            std::vector<RealType> hax(nNodes), hay(nNodes), haz(nNodes, RealType(0));
            for (size_t i = 0; i < nNodes; ++i) { hax[i] = -(hy[i] - yc); hay[i] = (hx[i] - xc); }
            d_avx.assign(hax.begin(), hax.end());
            d_avy.assign(hay.begin(), hay.end());
            d_avz.assign(haz.begin(), haz.end());
            avxp = thrust::raw_pointer_cast(d_avx.data());
            avyp = thrust::raw_pointer_cast(d_avy.data());
            avzp = thrust::raw_pointer_cast(d_avz.data());
        }

        assembleCoupledStokesKernel<KeyType, RealType><<<eBlocks, blockSize>>>(
            c0, c1, c2, c3, nx, ny, nz,
            thrust::raw_pointer_cast(d_rowOff.data()),
            thrust::raw_pointer_cast(d_colInd.data()),
            thrust::raw_pointer_cast(d_vals.data()),
            nu, tau, avxp, avyp, avzp, startEl, numLocal);
        cudaError_t aerr = cudaDeviceSynchronize();
        if (aerr != cudaSuccess) {
            std::cerr << "[phase0][assemble] kernel failed: " << cudaGetErrorString(aerr) << "\n";
            MPI_Finalize();
            return 4;
        }

        // copy CSR back; gates are host-side (Phase-0 model operator, small mesh)
        std::vector<RealType> hv(nnz);
        thrust::copy(d_vals.begin(), d_vals.end(), hv.begin());
        std::map<std::pair<int, int>, RealType> A;
        for (int r = 0; r < ND; ++r)
            for (int k = h_rowOff[r]; k < h_rowOff[r + 1]; ++k) A[{r, h_colInd[k]}] = hv[k];
        auto get = [&](int r, int c) -> RealType {
            auto it = A.find({r, c}); return it == A.end() ? RealType(0) : it->second;
        };

        RealType e1 = 0, e2 = 0, e3 = 0, e4 = 0, e5 = 0;
        for (size_t i = 0; i < nNodes; ++i)
            for (int j : ring[i])
                for (int d = 0; d < 3; ++d) {
                    e1 = std::max(e1, std::abs(get(4 * i + 3, 4 * j + d) + get(4 * j + d, 4 * i + 3))); // a^pu==-(a^up)^T
                    e4 = std::max(e4, std::abs(get(4 * i + d, 4 * j + d) - get(4 * j + d, 4 * i + d))); // a^uu symmetry residual
                }
        for (size_t i = 0; i < nNodes; ++i) {
            RealType spp = 0, sdc = 0, suu[3] = {0, 0, 0};
            for (int j : ring[i]) {
                spp += get(4 * i + 3, 4 * j + 3);
                for (int d = 0; d < 3; ++d) { sdc += get(4 * i + 3, 4 * j + d); suu[d] += get(4 * i + d, 4 * j + d); }
            }
            e2 = std::max(e2, std::abs(spp));                     // a^pp . 1 == 0
            e3 += sdc;                                            // global sum(D.const)
            for (int d = 0; d < 3; ++d) e5 = std::max(e5, std::abs(suu[d])); // a^uu . 1 (conservation)
        }
        e3 = std::abs(e3);

        const RealType atol = 1e-12;
        auto pf = [&](bool ok) { return ok ? "PASS" : "FAIL"; };
        std::cout << std::scientific << std::setprecision(3);
        std::cout << "[phase0][assemble] DOFs=" << ND << " nnz=" << nnz
                  << (doAdvect ? "  (advection ON)" : "  (Stokes)") << "\n";
        std::cout << "[phase0][assemble] a_pu==-a_up^T : " << e1 << "  " << pf(e1 < atol) << "\n";
        std::cout << "[phase0][assemble] a_pp.1==0     : " << e2 << "  " << pf(e2 < atol) << "\n";
        std::cout << "[phase0][assemble] sum(D.const)  : " << e3 << "  " << pf(e3 < atol) << "\n";
        if (doAdvect) {
            std::cout << "[phase0][assemble] a_uu nonsymm  : " << e4 << "  " << pf(e4 > atol) << "\n";
            std::cout << "[phase0][assemble] a_uu conserves: " << e5 << "  " << pf(e5 < atol) << "\n";
            assemblePass = (e1 < atol) && (e2 < atol) && (e3 < atol) && (e4 > atol) && (e5 < atol);
        } else {
            std::cout << "[phase0][assemble] a_uu symmetric: " << e4 << "  " << pf(e4 < atol) << "\n";
            assemblePass = (e1 < atol) && (e2 < atol) && (e3 < atol) && (e4 < atol);
        }
        std::cout << "[phase0][assemble] -> " << (assemblePass ? "ALL PASS" : "FAIL") << "\n";

        // ---- ACM Stage 3: GPU coarse-operator (sum-of-fine == P^T A P) + grid transfers ----
        if (doAcm3)
        {
            // host-once aggregation (node-granular, greedy over the 1-ring). The map is an INPUT to
            // Stage 3 -- any valid map validates the GPU coarse-operator recipe; the GPU directional
            // aggregation is Stage 5.
            const int KMAX = 4;
            std::vector<int> agg(nNodes, -1);
            int nCoarse = 0;
            for (size_t i = 0; i < nNodes; ++i) {
                if (agg[i] != -1) continue;
                agg[i] = nCoarse; int cnt = 1;
                for (int j : ring[i]) {
                    if (cnt >= KMAX) break;
                    if (j != (int)i && agg[j] == -1) { agg[j] = nCoarse; ++cnt; }
                }
                ++nCoarse;
            }
            const int NDc = 4 * nCoarse;

            thrust::device_vector<int> d_agg(agg.begin(), agg.end());
            thrust::device_vector<int> d_crowOff, d_ccolInd;
            thrust::device_vector<RealType> d_cvals;
            buildCoarseOperator<RealType>(d_rowOff, d_colInd, d_vals, d_agg,
                                          (int)nNodes, nCoarse, d_crowOff, d_ccolInd, d_cvals);

            // host reference: the same map applied to the fine CSR by an independent accumulation
            auto cdof = [&](int r) { return 4 * agg[r >> 2] + (r & 3); };
            std::map<std::pair<int, int>, RealType> Aref;
            for (int r = 0; r < ND; ++r)
                for (int k = h_rowOff[r]; k < h_rowOff[r + 1]; ++k)
                    Aref[{cdof(r), cdof(h_colInd[k])}] += hv[k];

            const int ccnnz = (int)d_cvals.size();
            std::vector<int> hcrow(NDc + 1), hccol(ccnnz);
            std::vector<RealType> hcval(ccnnz);
            thrust::copy(d_crowOff.begin(), d_crowOff.end(), hcrow.begin());
            thrust::copy(d_ccolInd.begin(), d_ccolInd.end(), hccol.begin());
            thrust::copy(d_cvals.begin(),   d_cvals.end(),   hcval.begin());
            std::map<std::pair<int, int>, RealType> Agpu;
            for (int r = 0; r < NDc; ++r)
                for (int k = hcrow[r]; k < hcrow[r + 1]; ++k) Agpu[{r, hccol[k]}] = hcval[k];
            RealType dRecipe = 0;
            for (auto& kv : Aref) { auto it = Agpu.find(kv.first); dRecipe = std::max(dRecipe, std::abs(kv.second - (it == Agpu.end() ? RealType(0) : it->second))); }
            for (auto& kv : Agpu) { auto it = Aref.find(kv.first); dRecipe = std::max(dRecipe, std::abs(kv.second - (it == Aref.end() ? RealType(0) : it->second))); }

            // Galerkin consistency: P^T (A (P x)) == A_coarse x for a random coarse x (host)
            std::vector<RealType> xc(NDc), xf(ND, 0), yf(ND, 0), yc(NDc, 0), zc(NDc, 0), rc(NDc, 0);
            for (int i = 0; i < NDc; ++i) xc[i] = std::sin(RealType(0.7) * i + RealType(1));
            for (size_t i = 0; i < nNodes; ++i) for (int c = 0; c < 4; ++c) xf[4 * i + c] = xc[4 * agg[i] + c];   // P xc
            for (size_t i = 0; i < nNodes; ++i) for (int c = 0; c < 4; ++c) rc[4 * agg[i] + c] += xf[4 * i + c];  // P^T xf (restrict-kernel ref)
            for (int r = 0; r < ND; ++r) { RealType s = 0; for (int k = h_rowOff[r]; k < h_rowOff[r + 1]; ++k) s += hv[k] * xf[h_colInd[k]]; yf[r] = s; }  // A xf
            for (size_t i = 0; i < nNodes; ++i) for (int c = 0; c < 4; ++c) yc[4 * agg[i] + c] += yf[4 * i + c];  // P^T A xf
            for (auto& kv : Aref) zc[kv.first.first] += kv.second * xc[kv.first.second];                          // A_coarse xc
            RealType dGal = 0; for (int i = 0; i < NDc; ++i) dGal = std::max(dGal, std::abs(yc[i] - zc[i]));

            // GPU transfer kernels vs host (injection prolongation + sum restriction)
            thrust::device_vector<RealType> d_xc(xc.begin(), xc.end()), d_xf(ND, RealType(0)), d_yc(NDc, RealType(0));
            int nb = ((int)nNodes + blockSize - 1) / blockSize;
            acmInjectKernel<RealType><<<nb, blockSize>>>(thrust::raw_pointer_cast(d_xc.data()),
                thrust::raw_pointer_cast(d_xf.data()), thrust::raw_pointer_cast(d_agg.data()), (int)nNodes);
            acmRestrictAddKernel<RealType><<<nb, blockSize>>>(thrust::raw_pointer_cast(d_xf.data()),
                thrust::raw_pointer_cast(d_yc.data()), thrust::raw_pointer_cast(d_agg.data()), (int)nNodes);
            cudaDeviceSynchronize();
            std::vector<RealType> gxf(ND), gyc(NDc);
            thrust::copy(d_xf.begin(), d_xf.end(), gxf.begin());
            thrust::copy(d_yc.begin(), d_yc.end(), gyc.begin());
            RealType dInj = 0, dRes = 0;
            for (int i = 0; i < ND;  ++i) dInj = std::max(dInj, std::abs(gxf[i] - xf[i]));
            for (int i = 0; i < NDc; ++i) dRes = std::max(dRes, std::abs(gyc[i] - rc[i]));

            const RealType ctol = 1e-10;
            auto pf2 = [&](bool ok) { return ok ? "PASS" : "FAIL"; };
            std::cout << "[phase0][acm3] fine DOFs=" << ND << " -> coarse DOFs=" << NDc
                      << " (aggregates=" << nCoarse << ", coarse nnz=" << ccnnz << ")\n";
            std::cout << "[phase0][acm3] GPU coarse == host P^T A P : " << dRecipe << "  " << pf2(dRecipe < ctol) << "\n";
            std::cout << "[phase0][acm3] Galerkin P^T A P x == Ac x : " << dGal    << "  " << pf2(dGal    < ctol) << "\n";
            std::cout << "[phase0][acm3] GPU inject == host          : " << dInj    << "  " << pf2(dInj    < ctol) << "\n";
            std::cout << "[phase0][acm3] GPU restrict == host        : " << dRes    << "  " << pf2(dRes    < ctol) << "\n";
            std::cout << "[phase0][acm3] -> "
                      << ((dRecipe < ctol && dGal < ctol && dInj < ctol && dRes < ctol) ? "ALL PASS" : "FAIL") << "\n";
        }

        // ---- ACM Stage 4a: GPU multilevel V-cycle vs host V-cycle replica (one cycle, exact match) ----
        if (doAcm4)
        {
            std::vector<AcmLevel<RealType>> levels;
            acmBuildHierarchy<RealType>(d_rowOff, d_colInd, d_vals, (int)nNodes,
                                        /*kmax=*/4, /*maxCoarseND=*/64, levels);

            std::vector<RealType> hb(ND);
            for (int i = 0; i < ND; ++i) hb[i] = std::sin(RealType(0.3) * i + RealType(0.5));   // deterministic test RHS

            // GPU: one V-cycle on (b, x=0)
            // coarsest = 10 Jacobi sweeps (the paper's value). NOTE: Jacobi is non-contractive on the
            // indefinite coarse operator (it mildly amplifies) -- a Stage-4a stopgap; Algorithm A's
            // direct coarsest solve (Stage 6) is the proper treatment, and GMRES (Stage 4b) makes the
            // cycle effective regardless.
            thrust::copy(hb.begin(), hb.end(), levels[0].bvec.begin());
            thrust::fill(levels[0].xvec.begin(), levels[0].xvec.end(), RealType(0));
            acmVcycleGpu<RealType>(levels, 0, /*pre=*/2, /*post=*/1, /*coarse=*/10, RealType(0.7));
            cudaDeviceSynchronize();
            std::vector<RealType> gx(ND);
            thrust::copy(levels[0].xvec.begin(), levels[0].xvec.end(), gx.begin());

            // host replica: same hierarchy + same arithmetic, one V-cycle on (b, x=0)
            std::vector<AcmLevelHost<RealType>> hLevels;
            acmHierarchyToHost<RealType>(levels, hLevels);
            std::vector<RealType> hx(ND, RealType(0));
            acmHostVcycle<RealType>(hLevels, 0, hx, hb, 2, 1, 10, RealType(0.7));

            // GPU SpMV/Jacobi contract to FMA, the host does not -> compare RELATIVE to ||x|| (a logic
            // bug gives O(1), not ~1e-9; this gate validates kernel correctness, not FP-bit-equality).
            RealType dCycle = 0, mhx = 0;
            for (int i = 0; i < ND; ++i) { dCycle = std::max(dCycle, std::abs(gx[i] - hx[i])); mhx = std::max(mhx, std::abs(hx[i])); }
            const RealType dRel = dCycle / std::max(RealType(1e-30), mhx);

            std::cout << "[phase0][acm4] V-cycle levels=" << levels.size() << "  sizes(ND):";
            for (auto& L : levels) std::cout << " " << L.ND;
            std::cout << "\n";
            const RealType ctol = 1e-7;
            std::cout << "[phase0][acm4] GPU V-cycle == host replica : abs=" << dCycle << " rel=" << dRel
                      << "  " << (dRel < ctol ? "PASS" : "FAIL") << "\n";
            std::cout << "[phase0][acm4] -> " << (dRel < ctol ? "ALL PASS" : "FAIL") << "\n";
        }

        // ---- Phase 1 1a (Hypre point-block) + ACM Stage 4b (FlexGMRES): shared BC-eliminated
        //      Stokes operator + cusolverSp QR reference, then each solver's own benchmark ----
        if (doHypre || doAcm4b)
        {
            const RealType nuS  = RealType(1) / static_cast<RealType>(Re);
            const RealType hpar = RealType(1) / std::sqrt(static_cast<RealType>(nNodes) * RealType(0.5));
            const RealType tauS = hpar * hpar / (RealType(4) * nuS);
            // re-assemble Stokes (advection off; matches the host pyamg study)
            thrust::fill(d_vals.begin(), d_vals.end(), RealType(0));
            assembleCoupledStokesKernel<KeyType, RealType><<<eBlocks, blockSize>>>(
                c0, c1, c2, c3, nx, ny, nz,
                thrust::raw_pointer_cast(d_rowOff.data()), thrust::raw_pointer_cast(d_colInd.data()),
                thrust::raw_pointer_cast(d_vals.data()), nuS, tauS, nullptr, nullptr, nullptr,
                startEl, numLocal);
            cudaDeviceSynchronize();
            // BC tags + eliminate (lid u=1 / walls / w=0 / pressure pin)
            const double bb = beta * M_PI / 180.0;
            const RealType sinb = std::sin(bb), cotb = std::cos(bb) / std::sin(bb), eps = 1e-6;
            std::vector<uint8_t> bcF(ND, 0); std::vector<RealType> bcV(ND, RealType(0));
            for (size_t i = 0; i < nNodes; ++i) {
                RealType y0 = hy[i] / sinb, x0 = hx[i] - hy[i] * cotb;
                bcF[4 * i + 2] = 1;
                if (y0 > 1.0 - eps) { bcF[4 * i + 0] = 1; bcV[4 * i + 0] = RealType(1); bcF[4 * i + 1] = 1; }
                else if (x0 < eps || x0 > 1.0 - eps || y0 < eps) { bcF[4 * i + 0] = 1; bcF[4 * i + 1] = 1; }
            }
            bcF[3] = 1;
            thrust::device_vector<uint8_t> d_bcF(bcF.begin(), bcF.end());
            thrust::device_vector<RealType> d_bcV(bcV.begin(), bcV.end());
            thrust::device_vector<RealType> d_rhs(ND, RealType(0));
            const int dblk = (ND + blockSize - 1) / blockSize;
            applyBCKernel<RealType><<<dblk, blockSize>>>(
                thrust::raw_pointer_cast(d_bcF.data()), thrust::raw_pointer_cast(d_bcV.data()),
                thrust::raw_pointer_cast(d_rowOff.data()), thrust::raw_pointer_cast(d_colInd.data()),
                thrust::raw_pointer_cast(d_vals.data()), thrust::raw_pointer_cast(d_rhs.data()), ND);
            cudaDeviceSynchronize();
            // cusolver reference solution
            thrust::device_vector<RealType> d_xref(ND, RealType(0));
            { cusolverSpHandle_t cs = nullptr; cusolverSpCreate(&cs);
              cusparseMatDescr_t de = nullptr; cusparseCreateMatDescr(&de);
              cusparseSetMatType(de, CUSPARSE_MATRIX_TYPE_GENERAL);
              cusparseSetMatIndexBase(de, CUSPARSE_INDEX_BASE_ZERO); int sg = -2;
              cusolverSpDcsrlsvqr(cs, ND, nnz, de,
                  thrust::raw_pointer_cast(d_vals.data()), thrust::raw_pointer_cast(d_rowOff.data()),
                  thrust::raw_pointer_cast(d_colInd.data()), thrust::raw_pointer_cast(d_rhs.data()),
                  1e-12, 1, thrust::raw_pointer_cast(d_xref.data()), &sg);
              cudaDeviceSynchronize(); cusparseDestroyMatDescr(de); cusolverSpDestroy(cs); }
            // wrap the BC-eliminated CSR into a SparseMatrix (shared by the Hypre + ACM solves)
            SparseMatrix<int, RealType, cstone::GpuTag> Am; Am.allocate(ND, ND, nnz);
            thrust::copy(d_rowOff.begin(), d_rowOff.end(), thrust::device_pointer_cast(Am.rowOffsetsPtr()));
            thrust::copy(d_colInd.begin(), d_colInd.end(), thrust::device_pointer_cast(Am.colIndicesPtr()));
            thrust::copy(d_vals.begin(),   d_vals.end(),   thrust::device_pointer_cast(Am.valuesPtr()));

            if (doHypre)
            {
            // Hypre point-block BoomerAMG + GMRES
            cstone::DeviceVector<RealType> bH, xH; bH.resize(ND); xH.resize(ND);
            thrust::copy(d_rhs.begin(), d_rhs.end(), thrust::device_pointer_cast(bH.data()));
            HypreGMRESSolver<RealType, int, cstone::GpuTag> hg(
                MPI_COMM_WORLD, 1000, 1e-8,
                HypreGMRESSolver<RealType, int, cstone::GpuTag>::BOOMERAMG, 50);
            hg.setVerbose(false);
            hg.setPointBlock(4);
            hg.solve(Am, bH, xH, 0, ND, 0, ND, std::vector<int>{});
            const int iters = hg.getLastIterations();
            // compare to cusolver reference (guards the false-x=0 convergence mode)
            thrust::device_vector<RealType> d_xh(ND), d_df(ND);
            thrust::copy(thrust::device_pointer_cast(xH.data()),
                         thrust::device_pointer_cast(xH.data() + ND), d_xh.begin());
            thrust::transform(d_xh.begin(), d_xh.end(), d_xref.begin(), d_df.begin(),
                [] __device__ (RealType a, RealType c) { return fabs(a - c); });
            RealType dmax = thrust::reduce(d_df.begin(), d_df.end(), RealType(0), thrust::maximum<RealType>());
            thrust::transform(d_xref.begin(), d_xref.end(), d_df.begin(),
                [] __device__ (RealType c) { return fabs(c); });
            RealType rmax = thrust::reduce(d_df.begin(), d_df.end(), RealType(0), thrust::maximum<RealType>());
            RealType rel = dmax / (rmax > 0 ? rmax : RealType(1));
            // rel guards the false-x=0 collapse (rel~1 = AMG returned ~0). A real
            // GMRES solve to a 1e-8 RESIDUAL tol on a kappa-growing saddle operator
            // legitimately differs from the direct solve by ~kappa*1e-8, so allow a
            // generous band; it is NOT a precision check.
            bool hyOk = (iters > 0) && (rel < 5e-2);
            std::cout << std::scientific << std::setprecision(3);
            std::cout << "[phase0][hypre] Re=" << Re << " point-block BoomerAMG GMRES iters=" << iters
                      << " nnz=" << nnz << " |x_hypre-x_qr|rel=" << rel
                      << "  -> " << (hyOk ? "PASS (converged + matches reference)" : "CHECK") << "\n";
            assemblePass = assemblePass && hyOk;
            }  // if (doHypre)

            // ---- ACM Stage 4b: ACM-preconditioned native FlexGMRES vs the 1a baseline + cusolver ref ----
            if (doAcm4b)
            {
                using Vec = cstone::DeviceVector<RealType>;
                Vec bvec; bvec.resize(ND);
                thrust::copy(d_rhs.begin(), d_rhs.end(), thrust::device_pointer_cast(bvec.data()));

                // max|x_xref| for the relative-error denominator
                thrust::device_vector<RealType> d_absref(ND);
                thrust::transform(d_xref.begin(), d_xref.end(), d_absref.begin(),
                    [] __device__ (RealType c) { return fabs(c); });
                const RealType rmax = thrust::reduce(d_absref.begin(), d_absref.end(), RealType(0), thrust::maximum<RealType>());
                auto relErr = [&](Vec& xv) -> RealType {
                    thrust::device_vector<RealType> xt(ND);
                    thrust::copy(thrust::device_pointer_cast(xv.data()),
                                 thrust::device_pointer_cast(xv.data() + ND), xt.begin());
                    thrust::transform(xt.begin(), xt.end(), d_xref.begin(), xt.begin(),
                        [] __device__ (RealType a, RealType c) { return fabs(a - c); });
                    RealType dmax = thrust::reduce(xt.begin(), xt.end(), RealType(0), thrust::maximum<RealType>());
                    return dmax / (rmax > 0 ? rmax : RealType(1));
                };

                std::cout << std::scientific << std::setprecision(3);

                // GATE 1 (Z-basis correctness): with an IDENTITY preconditioner, flexible GMRES must
                // reduce EXACTLY to plain GMRES (Z_j=V_j, x=x0+Z y = x0+V y) -- a convergence-independent
                // algebraic identity. A Z-basis reconstruction bug breaks this even if neither converges.
                // (Jacobi is too weak to converge on this indefinite saddle operator, so it cannot be the
                // oracle; identity can.)
                Vec xp; xp.resize(ND);
                GMRESSolver<RealType, int, cstone::GpuTag> gp(2000, 1e-8, 30);
                gp.setVerbose(false); gp.solve(Am, bvec, xp, false);            // plain (unpreconditioned)
                const int itP = gp.getLastIterations();
                IdentityPreconditioner<RealType, int, cstone::GpuTag> id;
                Vec xi; xi.resize(ND);
                GMRESSolver<RealType, int, cstone::GpuTag> gi(2000, 1e-8, 30);
                gi.setVerbose(false); gi.setPreconditioner(&id); gi.setFlexible(true);
                gi.solve(Am, bvec, xi, true);                                   // flexible, M = I
                const int itI = gi.getLastIterations();
                thrust::device_vector<RealType> xpv(ND), xiv(ND);
                thrust::copy(thrust::device_pointer_cast(xp.data()), thrust::device_pointer_cast(xp.data() + ND), xpv.begin());
                thrust::copy(thrust::device_pointer_cast(xi.data()), thrust::device_pointer_cast(xi.data() + ND), xiv.begin());
                thrust::transform(xpv.begin(), xpv.end(), xiv.begin(), xiv.begin(),
                    [] __device__ (RealType a, RealType b) { return fabs(a - b); });
                const RealType ddiff = thrust::reduce(xiv.begin(), xiv.end(), RealType(0), thrust::maximum<RealType>());
                thrust::transform(xpv.begin(), xpv.end(), xpv.begin(), [] __device__ (RealType a) { return fabs(a); });
                const RealType dpmax = thrust::reduce(xpv.begin(), xpv.end(), RealType(0), thrust::maximum<RealType>());
                const RealType relZ = ddiff / (dpmax > 0 ? dpmax : RealType(1));
                const bool gate1 = (itP == itI) && (relZ < 1e-10);
                std::cout << "[phase0][acm4b][gate1] FGMRES(I)==plain GMRES (Z-basis): iters " << itP << "/" << itI
                          << " rel=" << relZ << "  -> " << (gate1 ? "PASS" : "FAIL") << "\n";

                // GATE 2 (ACM quality): ACM-preconditioned FlexGMRES converges to the reference.
                GpuAcmPreconditioner<RealType, int, cstone::GpuTag> acm;
                Vec xa; xa.resize(ND);
                GMRESSolver<RealType, int, cstone::GpuTag> ga(2000, 1e-8, 30);
                ga.setVerbose(false); ga.setPreconditioner(&acm); ga.setFlexible(true);
                ga.solve(Am, bvec, xa, true);
                const int itA = ga.getLastIterations();
                const RealType relA = relErr(xa);
                const bool gate2 = (itA > 0) && (itA < 1980) && (relA < 5e-2);   // converged (< maxIter) + correct
                std::cout << "[phase0][acm4b] Re=" << Re << " ACM-FlexGMRES iters=" << itA
                          << " (Hypre 1a baseline n16/32/64 = 13/19/39) levels=" << acm.numLevels()
                          << " |x_acm-x_qr|rel=" << relA
                          << "  -> " << (gate1 && gate2 ? "PASS" : "CHECK") << "\n";
                assemblePass = assemblePass && gate1 && gate2;
            }
        }

        // ---- 3b: Picard loop + cusolverSp QR direct reference solve + benchmark match ----
        if (doSolve)
        {
            // Erturk-Dursun (arXiv:physics/0505121) Table 5, alpha=45: u along A-B, Re=100.
            const RealType ED_s[17] = {0,32,64,96,128,160,192,224,256,288,320,352,384,416,448,480,512};
            const RealType ED_u[17] = {0.0,-2.389e-3,-8.901e-3,-1.949e-2,-3.426e-2,-5.373e-2,
                                       -7.846e-2,-1.080e-1,-1.390e-1,-1.634e-1,-1.669e-1,-1.315e-1,
                                       -3.945e-2,1.217e-1,3.544e-1,6.477e-1,1.0};
            auto edInterp = [&](RealType s) -> RealType {
                if (s <= 0) return ED_u[0];
                for (int k = 1; k < 17; ++k)
                    if (s <= ED_s[k]) {
                        RealType t = (s - ED_s[k-1]) / (ED_s[k] - ED_s[k-1]);
                        return ED_u[k-1] + t * (ED_u[k] - ED_u[k-1]);
                    }
                return ED_u[16];
            };

            // physics: nu = 1/Re, PSPG tau = h^2/(4 nu), h ~ in-plane node spacing
            const RealType nuS  = RealType(1) / static_cast<RealType>(Re);
            const RealType hpar = RealType(1) / std::sqrt(static_cast<RealType>(nNodes) * RealType(0.5));
            const RealType tauS = hpar * hpar / (RealType(4) * nuS);

            // BC tags (one-time host setup): lid u=1/v=0, no-slip walls, w=0 (thin slab), pressure pin
            const double bb = beta * M_PI / 180.0;
            const RealType sinb = std::sin(bb), cotb = std::cos(bb) / std::sin(bb), eps = 1e-6;
            std::vector<uint8_t> bcF(ND, 0);
            std::vector<RealType> bcV(ND, RealType(0));
            for (size_t i = 0; i < nNodes; ++i) {
                RealType y0 = hy[i] / sinb, x0 = hx[i] - hy[i] * cotb;
                bcF[4 * i + 2] = 1;
                if (y0 > 1.0 - eps) { bcF[4 * i + 0] = 1; bcV[4 * i + 0] = RealType(1); bcF[4 * i + 1] = 1; }
                else if (x0 < eps || x0 > 1.0 - eps || y0 < eps) { bcF[4 * i + 0] = 1; bcF[4 * i + 1] = 1; }
            }
            bcF[3] = 1;
            thrust::device_vector<uint8_t>  d_bcF(bcF.begin(), bcF.end());
            thrust::device_vector<RealType> d_bcV(bcV.begin(), bcV.end());
            thrust::device_vector<RealType> d_avx(nNodes, RealType(0)), d_avy(nNodes, RealType(0)),
                                            d_avz(nNodes, RealType(0)), d_uprev(nNodes, RealType(0)),
                                            d_tmp(nNodes);
            thrust::device_vector<RealType> d_rhs(ND, RealType(0)), d_x(ND, RealType(0));
            const int dblk = (ND + blockSize - 1) / blockSize;
            const int nblk = (static_cast<int>(nNodes) + blockSize - 1) / blockSize;

            cusolverSpHandle_t cs = nullptr; cusolverSpCreate(&cs);
            cusparseMatDescr_t descr = nullptr; cusparseCreateMatDescr(&descr);
            cusparseSetMatType(descr, CUSPARSE_MATRIX_TYPE_GENERAL);
            cusparseSetMatIndexBase(descr, CUSPARSE_INDEX_BASE_ZERO);

            int it = 0, singularity = -2; RealType du = 0;
            cusolverStatus_t st = CUSOLVER_STATUS_SUCCESS;
            for (it = 0; it < picard; ++it) {
                // re-assemble with the frozen (previous-iterate) velocity, proper nu/tau
                thrust::fill(d_vals.begin(), d_vals.end(), RealType(0));
                assembleCoupledStokesKernel<KeyType, RealType><<<eBlocks, blockSize>>>(
                    c0, c1, c2, c3, nx, ny, nz,
                    thrust::raw_pointer_cast(d_rowOff.data()), thrust::raw_pointer_cast(d_colInd.data()),
                    thrust::raw_pointer_cast(d_vals.data()), nuS, tauS,
                    thrust::raw_pointer_cast(d_avx.data()), thrust::raw_pointer_cast(d_avy.data()),
                    thrust::raw_pointer_cast(d_avz.data()), startEl, numLocal);
                cudaDeviceSynchronize();
                thrust::fill(d_rhs.begin(), d_rhs.end(), RealType(0));
                applyBCKernel<RealType><<<dblk, blockSize>>>(
                    thrust::raw_pointer_cast(d_bcF.data()), thrust::raw_pointer_cast(d_bcV.data()),
                    thrust::raw_pointer_cast(d_rowOff.data()), thrust::raw_pointer_cast(d_colInd.data()),
                    thrust::raw_pointer_cast(d_vals.data()), thrust::raw_pointer_cast(d_rhs.data()), ND);
                cudaDeviceSynchronize();
                st = cusolverSpDcsrlsvqr(cs, ND, nnz, descr,
                    thrust::raw_pointer_cast(d_vals.data()), thrust::raw_pointer_cast(d_rowOff.data()),
                    thrust::raw_pointer_cast(d_colInd.data()), thrust::raw_pointer_cast(d_rhs.data()),
                    1e-12, 1, thrust::raw_pointer_cast(d_x.data()), &singularity);
                cudaDeviceSynchronize();
                if (st != CUSOLVER_STATUS_SUCCESS || singularity >= 0) break;
                extractCompKernel<RealType><<<nblk, blockSize>>>(thrust::raw_pointer_cast(d_x.data()),
                    thrust::raw_pointer_cast(d_avx.data()), 4, 0, (int)nNodes);
                extractCompKernel<RealType><<<nblk, blockSize>>>(thrust::raw_pointer_cast(d_x.data()),
                    thrust::raw_pointer_cast(d_avy.data()), 4, 1, (int)nNodes);
                extractCompKernel<RealType><<<nblk, blockSize>>>(thrust::raw_pointer_cast(d_x.data()),
                    thrust::raw_pointer_cast(d_avz.data()), 4, 2, (int)nNodes);
                cudaDeviceSynchronize();
                thrust::transform(d_avx.begin(), d_avx.end(), d_uprev.begin(), d_tmp.begin(),
                    [] __device__ (RealType a, RealType c) { return fabs(a - c); });
                du = thrust::reduce(d_tmp.begin(), d_tmp.end(), RealType(0), thrust::maximum<RealType>());
                thrust::copy(d_avx.begin(), d_avx.end(), d_uprev.begin());
                if (du < 1e-8) { ++it; break; }
            }

            // centerline u (line A-B: x0~0.5, one z-layer) vs Erturk-Dursun -- host post-process diagnostic
            std::vector<RealType> hu(nNodes);
            thrust::copy(d_avx.begin(), d_avx.end(), hu.begin());
            RealType zmin = *std::min_element(hz.begin(), hz.end());
            std::vector<std::pair<RealType, RealType>> prof;
            for (size_t i = 0; i < nNodes; ++i) {
                RealType x0 = hx[i] - hy[i] * cotb, y0 = hy[i] / sinb;
                if (std::abs(x0 - 0.5) < 1e-6 && std::abs(hz[i] - zmin) < 1e-6)
                    prof.push_back({y0, hu[i]});
            }
            std::sort(prof.begin(), prof.end());
            RealType emax = 0, el2 = 0;
            for (auto& pr : prof) {
                RealType e = std::abs(pr.second - edInterp(pr.first * RealType(512)));
                emax = std::max(emax, e); el2 += e * e;
            }
            int np = static_cast<int>(prof.size());
            el2 = np > 0 ? std::sqrt(el2 / np) : 0;
            bool solveOk = (st == CUSOLVER_STATUS_SUCCESS) && (singularity < 0) && (np > 2) && (emax < 6e-2);
            std::cout << std::scientific << std::setprecision(3);
            std::cout << "[phase0][solve] Re=" << Re << " picard=" << it << " d(u)=" << du
                      << " nu=" << nuS << " tau=" << tauS << " status=" << static_cast<int>(st)
                      << " sing=" << singularity << "\n";
            std::cout << "[phase0][solve] centerline u vs Erturk-Dursun(b45,Re100): n=" << np
                      << " max|err|=" << emax << " L2=" << el2 << "  -> "
                      << (solveOk ? "PASS (matches benchmark)" : "CHECK") << "\n";
            assemblePass = assemblePass && solveOk;
            cusparseDestroyMatDescr(descr);
            cusolverSpDestroy(cs);
        }
    }

    MPI_Finalize();
    return (pass && assemblePass) ? 0 : 3;
}
