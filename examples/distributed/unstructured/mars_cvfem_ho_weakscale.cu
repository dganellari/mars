// Multi-rank (multi-GPU) WEAK-SCALING driver for the high-order matrix-free CVFEM
// APPLY (Knaus Alg 2, mars_cvfem_ho_matfree.hpp) at p=1. This is the first wiring
// of a DISTRIBUTED matrix-free matvec in MARS. It mirrors the prior PASC deck's
// weak-scaling methodology (which scaled the p=1 ASSEMBLY) but for the matrix-free
// apply: hold ~2M elements/rank and scale rank count, reporting total MDOF/s and
// parallel efficiency vs the 1-rank baseline.
//
// THE DISTRIBUTED MATVEC y = A*u (the gap this driver fills over compare.cu):
//   1. exchangeNodeHalo(u)         FORWARD: ghost-DOF slots <- owners' u.
//   2. memset(y, 0)
//   3. ho_cvfem_apply_launch<1>    per-element apply+atomicAdd over OWNED elements
//      ONLY ([startIndex, startIndex+localElementCount)). The scatter writes both
//      owned-DOF rows AND ghost-DOF rows (a ghost DOF here is owned elsewhere).
//   4. reverseExchangeNodeHaloAdd(y) REVERSE-ADD: this rank's ghost-DOF partials
//      are SUMMED into the owners' owned-DOF slots on the owning ranks.
//   After step 4, y[0..numDofs) on every rank is the COMPLETE owned-row result.
//
// WHY owned-elements-only (not owned+halo): the apply's atomicAdd scatter has no
// ownership mask, so running it over halo elements too would double-count the
// seam-shared element contributions. Owned-only + reverse-add is the unique
// parity-exact scheme (sum is commutative; only FP reduction order differs from
// single-rank, ~1e-13). This matches the NS-solver conservative-scatter pattern
// (mars_ns_channel_solver.hpp:3185 -- owned scatter -> reverseAdd -> forward).
//
// PARITY GATES (partition-invariant, hold on any rank count):
//   A*1 = 0       constant null space of a pure-Neumann diffusion operator.
//   A*linear = 0  on INTERIOR owned rows. A linear field has constant gradient, so
//                 over a CLOSED sub-control-volume the net flux is zero. Excluded:
//                 geometric domain-boundary nodes (open control volume; the bare
//                 apply uses natural zero-flux, so a linear field's flux is
//                 unbalanced -> nonzero there) and partition-seam nodes (stencil
//                 split). The interior mask is global-box + halo-send-list based,
//                 NOT ownership: buildNodeOwnership emits only owned/ghost, every
//                 owned node (boundary, seam included) is ownership==1. Linear field
//                 is keyed by node coordinate so it is identical across partitions.
// Optionally --dump-parity=<file> writes (globalSfcKey, y) for every owned DOF so a
// 1-rank and an N-rank run on the SAME GLOBAL MESH (identical --ncells, or the same
// .exo; same global bounding box) can be diffed bit-for-bit externally (the true
// single-vs-multi parity). NOTE this is SEPARATE from the weak-scaling ladder, which
// deliberately uses DIFFERENT --ncells per rung -- do NOT diff dumps across rungs.
//
// Mesh: --mesh=<dir/.exo> (one file, cstone self-partitions on N ranks) OR
//       --ncells=N (procedural per-rank cube; pick N so N^3/numRanks ~ 2M for a
//       held-per-rank weak-scaling ladder, e.g. 126/252/400 on 1/8/32 ranks).
//
// This driver is its OWN executable because mars_cvfem_ho_matfree.hpp defines
// __constant__ banks that may be included by EXACTLY ONE translation unit.

#include "mars.hpp"
#include "backend/distributed/unstructured/domain.hpp"
#include "backend/distributed/unstructured/fem/mars_cvfem_utils.hpp"
#include "backend/distributed/unstructured/fem/mars_cvfem_ho_basis.hpp"
#include "backend/distributed/unstructured/fem/mars_cvfem_ho_matfree.hpp"
#include "backend/distributed/unstructured/utils/mars_generate_cube.hpp"

#include <cuda_runtime.h>
#include <mpi.h>
#include <thrust/extrema.h>      // min_element/max_element for the global box (Gate 2)
#include <thrust/for_each.h>     // overlap: scatter recv-ghost flags + pack on a stream
#include <thrust/copy.h>         // overlap: copy_if to partition interior/boundary lists
#include <thrust/iterator/counting_iterator.h>
#include <thrust/execution_policy.h>
#include <thrust/system/cuda/execution_policy.h>   // thrust::cuda::par.on(stream)
#include <algorithm>
#include <cstdio>
#include <cmath>
#include <string>
#include <vector>

using namespace mars;
using namespace mars::fem;

#define CK(call) do { cudaError_t e=(call); if(e!=cudaSuccess){ \
  printf("CUDA error %s at %s:%d\n", cudaGetErrorString(e), __FILE__, __LINE__); \
  MPI_Abort(MPI_COMM_WORLD, 1); } } while(0)

// MARS hex corner c -> HO local DOF index l = i*4 + j*2 + k (i slowest). Same
// reconciliation of the two corner conventions as mars_cvfem_ho_compare.cu; the
// metric kernel reads d_corners in MARS corner order, only elemDof is permuted.
__constant__ int c_hexCornerToHoLocal_ws[8] = {0, 4, 6, 2, 1, 5, 7, 3};

// One thread per OWNED element: gather 8 corner coords (MARS order) into
// d_corners[(e-startElem)*24 + c*3 + d] and 8 node DOFs into elemDof (HO order).
// elemBase is the global element index of the first OWNED element so corners/dof
// are packed densely over the owned slice [0, localElementCount).
template<typename KeyType, typename RealType>
__global__ void ho_p1_extract_owned_kernel(
    const KeyType* c0, const KeyType* c1, const KeyType* c2, const KeyType* c3,
    const KeyType* c4, const KeyType* c5, const KeyType* c6, const KeyType* c7,
    const RealType* __restrict__ d_x, const RealType* __restrict__ d_y, const RealType* __restrict__ d_z,
    const int* __restrict__ d_nodeToDof, size_t elemBase, size_t numLocal,
    RealType* __restrict__ d_corners, int* __restrict__ d_elemDof)
{
    size_t le = (size_t)blockIdx.x * blockDim.x + threadIdx.x;
    if (le >= numLocal) return;
    const size_t e = elemBase + le;
    const KeyType node[8] = {c0[e], c1[e], c2[e], c3[e], c4[e], c5[e], c6[e], c7[e]};
    #pragma unroll
    for (int c = 0; c < 8; ++c) {
        const KeyType nd = node[c];
        d_corners[le * 24 + c * 3 + 0] = d_x[nd];
        d_corners[le * 24 + c * 3 + 1] = d_y[nd];
        d_corners[le * 24 + c * 3 + 2] = d_z[nd];
        d_elemDof[le * 8 + c_hexCornerToHoLocal_ws[c]] = d_nodeToDof[nd];
    }
}

// Fill a node-indexed scalar field f(x,y,z) for the A*linear gate. We use a
// coordinate-defined field so the SAME value lands on a node regardless of which
// rank owns it -- a prerequisite for the gate to be partition-invariant.
template<typename RealType>
__global__ void fill_linear_field_kernel(
    const RealType* __restrict__ d_x, const RealType* __restrict__ d_y, const RealType* __restrict__ d_z,
    const int* __restrict__ d_nodeToDof, size_t nodeCount,
    RealType cx, RealType cy, RealType cz, RealType* __restrict__ d_u)
{
    size_t n = (size_t)blockIdx.x * blockDim.x + threadIdx.x;
    if (n >= nodeCount) return;
    int dof = d_nodeToDof[n];
    if (dof < 0) return;
    d_u[dof] = cx * d_x[n] + cy * d_y[n] + cz * d_z[n];
}

// A*linear=0 holds ONLY on rows whose full CVFEM stencil is assembled AND whose
// sub-control-volume surface is CLOSED. ownership tells us only owned(1)/ghost(0)
// -- it is NOT an interior flag (buildNodeOwnership never emits a "2"; every owned
// node, including domain-boundary and partition-seam nodes, is ownership==1). A
// domain-boundary node has an OPEN control volume (the bare matfree apply imposes
// natural zero-flux there), so a linear field's flux is unbalanced and A*linear!=0;
// a seam node's stencil is split across ranks. We therefore build the interior mask
// in three steps below: (1) mark every owned DOF, (2) un-mark geometric
// domain-boundary nodes (coord at the GLOBAL bounding-box min/max in any axis),
// (3) un-mark seam nodes (any node in the halo send list). This kernel does steps
// 1+2; step 3 is a separate scatter over the send-id list.
template<typename RealType>
__global__ void mark_interior_owned_kernel(
    const uint8_t* __restrict__ d_nodeOwnership, const int* __restrict__ d_nodeToDof,
    const RealType* __restrict__ d_x, const RealType* __restrict__ d_y, const RealType* __restrict__ d_z,
    size_t nodeCount, int numDofs,
    RealType bxlo, RealType bylo, RealType bzlo, RealType bxhi, RealType byhi, RealType bzhi,
    RealType tol, uint8_t* __restrict__ d_isInterior)
{
    size_t n = (size_t)blockIdx.x * blockDim.x + threadIdx.x;
    if (n >= nodeCount) return;
    int dof = d_nodeToDof[n];
    if (dof < 0 || dof >= numDofs) return;
    if (d_nodeOwnership[n] != 1) { d_isInterior[dof] = 0; return; }   // ghost
    const RealType x = d_x[n], y = d_y[n], z = d_z[n];
    const bool onBoundary =
        fabs((double)(x - bxlo)) < tol || fabs((double)(x - bxhi)) < tol ||
        fabs((double)(y - bylo)) < tol || fabs((double)(y - byhi)) < tol ||
        fabs((double)(z - bzlo)) < tol || fabs((double)(z - bzhi)) < tol;
    d_isInterior[dof] = onBoundary ? 0 : 1;
}

// Step 3: a node in the halo SEND list is a partition-seam node (its stencil is
// split across ranks) -> not interior. Scatter zeros over the send-id list.
__global__ void unmark_seam_kernel(
    const int* __restrict__ d_sendNodeIds, int sendCount, const int* __restrict__ d_nodeToDof,
    int numDofs, uint8_t* __restrict__ d_isInterior)
{
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i >= sendCount) return;
    int dof = d_nodeToDof[d_sendNodeIds[i]];
    if (dof >= 0 && dof < numDofs) d_isInterior[dof] = 0;
}

// --overlap classification. One thread per OWNED element le in [0,numLocal).
// Element e = elemBase + le. An element is BOUNDARY if ANY of its 8 corner nodes
// is a RECV-ghost node (its input u needs a halo value from another rank); else
// INTERIOR. We read raw connectivity C(c)[e] (MARS corner order) -- classification
// keys on node IDENTITY only, so the HO-local corner permutation is irrelevant.
// d_isRecvGhost is a per-node bool (1 at recvNodeIds_). flag=1 -> boundary.
template<typename KeyType>
__global__ void classify_owned_elem_kernel(
    const KeyType* c0, const KeyType* c1, const KeyType* c2, const KeyType* c3,
    const KeyType* c4, const KeyType* c5, const KeyType* c6, const KeyType* c7,
    const uint8_t* __restrict__ d_isRecvGhost, size_t elemBase, size_t numLocal,
    uint8_t* __restrict__ d_elemIsBoundary)
{
    size_t le = (size_t)blockIdx.x * blockDim.x + threadIdx.x;
    if (le >= numLocal) return;
    const size_t e = elemBase + le;
    const KeyType node[8] = {c0[e], c1[e], c2[e], c3[e], c4[e], c5[e], c6[e], c7[e]};
    uint8_t boundary = 0;
    #pragma unroll
    for (int c = 0; c < 8; ++c)
        boundary |= d_isRecvGhost[node[c]];   // OR: any recv-ghost corner -> boundary
    d_elemIsBoundary[le] = boundary;
}

int main(int argc, char** argv)
{
    MPI_Init(&argc, &argv);
    int rank = 0, numRanks = 1;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &numRanks);
    // Bind one GPU per rank. Wrap in CK so a device-binding failure (wrong GPU count
    // from srun) reports HERE, not at a confusing downstream CUDA call.
    int devCount = 0; CK(cudaGetDeviceCount(&devCount));
    if (devCount > 0) CK(cudaSetDevice(rank % devCount));

    using KeyType = uint64_t; using RealType = double; using ElemTag = HexTag;
    using Domain  = ElementDomain<ElemTag, RealType, KeyType, cstone::GpuTag>;

    std::string mesh; size_t ncells = 0; int iters = 50; std::string dumpFile;
    bool overlap = false;   // --overlap: hide the forward halo behind the interior apply
    // --irregular: opt-in genuinely-unstructured procedural mesh (warp+deform) so the
    // scaling test is not flattered by the perfect cube's free SFC locality. The deform
    // is DOMAIN-scale (amplitude + wavelength fixed in domain units), so the halo bump
    // is INDEPENDENT of --ncells -> it still holds at trillion scale (a perturbation
    // tied to local h would be cosmetic there). Default OFF -> byte-identical perfect
    // cube, existing runs unchanged.
    CubeIrregularity irr{};
    for (int i = 1; i < argc; ++i) { std::string a = argv[i];
        if (a.rfind("--mesh=", 0) == 0) mesh = a.substr(7);
        else if (a.rfind("--ncells=", 0) == 0) ncells = std::stoull(a.substr(9));
        else if (a.rfind("--iters=", 0) == 0) iters = std::stoi(a.substr(8));
        else if (a == "--overlap") overlap = true;
        else if (a == "--irregular") { irr.warp = true; irr.deform = true; }
        else if (a.rfind("--dump-parity=", 0) == 0) dumpFile = a.substr(14); }
    if (mesh.empty() && ncells == 0) {
        if (rank == 0) printf("need --mesh=<dir/.exo> OR --ncells=<N> [--iters=N] [--dump-parity=<file>]\n");
        MPI_Finalize(); return 1;
    }

    // --- build ElementDomain: file (self-partitions on N ranks) or procedural cube ---
    Domain* domainPtr = nullptr;
    if (!mesh.empty()) {
        domainPtr = new Domain(mesh, rank, numRanks, true, 64, 8u);
    } else {
        auto [genNodes, genElems, gx, gy, gz, lconn] =
            generateCubeElementPartition<RealType, KeyType>(ncells, rank, numRanks, irr);
        (void)genNodes; (void)genElems;
        typename Domain::HostCoordsTuple h_coords{std::move(gx), std::move(gy), std::move(gz)};
        typename Domain::HostConnectivityTuple h_conn{
            std::move(lconn[0]), std::move(lconn[1]), std::move(lconn[2]), std::move(lconn[3]),
            std::move(lconn[4]), std::move(lconn[5]), std::move(lconn[6]), std::move(lconn[7])};
        domainPtr = new Domain(h_coords, h_conn, rank, numRanks, 64, false, 8u);
    }
    Domain& domain = *domainPtr;

    const auto& d_nodeOwnership = domain.getNodeOwnershipMap();
    size_t nodeCount    = domain.getNodeCount();
    size_t startElem    = domain.startIndex();
    size_t numLocal     = domain.localElementCount();   // OWNED elements only

    // --- SHARED p=1 node DOF map (same builder as the assembly path). Owned nodes
    //     -> [0,numDofs), ghost nodes -> [numDofs, nodeCount). numTotalDofs sizes
    //     the u/y vectors so ghost-DOF slots exist for the halo round-trip. ---
    cstone::DeviceVector<int> d_nodeToDof(nodeCount);
    int numDofs      = buildDofMappingGpu<KeyType>(d_nodeOwnership.data(), d_nodeToDof.data(), nodeCount);
    int numTotalDofs = (int)nodeCount;

    const auto& d_conn = domain.getElementToNodeConnectivity();
    auto C = [&](int i)->const KeyType* {
        return (i==0?std::get<0>(d_conn):i==1?std::get<1>(d_conn):i==2?std::get<2>(d_conn):i==3?std::get<3>(d_conn):
                i==4?std::get<4>(d_conn):i==5?std::get<5>(d_conn):i==6?std::get<6>(d_conn):std::get<7>(d_conn)).data(); };

    domain.cacheNodeCoordinates();
    const auto& d_x = domain.getNodeX(); const auto& d_y = domain.getNodeY(); const auto& d_z = domain.getNodeZ();

    // --- reference operators -> constant memory (once) ---
    HoCvfemOperators op = buildHoCvfemOperators(1);
    // The hex-corner -> HO-local map is only valid if GLL nodes are ascending
    // (zeta[0]=-1 < zeta[1]=+1); assert so a future basis reorder trips loudly.
    if (!(op.zeta[0] < op.zeta[1])) {
        if (rank == 0) printf("FATAL: GLL zeta not ascending (%.3f,%.3f); corner map invalid.\n",
                              op.zeta[0], op.zeta[1]);
        MPI_Finalize(); return 1;
    }
    CK(ho_cvfem_upload_operators(1, op.Btil.data(), op.Dtil.data(), op.D.data(),
                                 op.W.data(), op.xi.data(), op.zeta.data()));

    // --- extract corners (MARS order) + elemDof (HO order) over OWNED elements ---
    cstone::DeviceVector<RealType> d_corners(numLocal * 24);
    cstone::DeviceVector<int>      d_elemDof(numLocal * 8);
    if (numLocal > 0) {
        int blk = 256, grid = (int)((numLocal + blk - 1) / blk);
        ho_p1_extract_owned_kernel<KeyType, RealType><<<grid, blk>>>(
            C(0),C(1),C(2),C(3),C(4),C(5),C(6),C(7),
            d_x.data(), d_y.data(), d_z.data(), d_nodeToDof.data(), startElem, numLocal,
            d_corners.data(), d_elemDof.data());
        CK(cudaGetLastError()); cudaDeviceSynchronize();
    }

    // --- per-element metric, precomputed ONCE over the owned-element slice ---
    constexpr int kHoG = 3 * 1 * (1 + 1) * (1 + 1) * 3;   // 36 doubles/element at P=1
    cstone::DeviceVector<RealType> d_G(numLocal * kHoG);
    if (numLocal > 0) {
        CK((ho_cvfem_metric_perpoint_launch<RealType, 1>(d_corners.data(), d_G.data(), numLocal)));
        cudaDeviceSynchronize();
    }

    // ====================== the distributed matvec y = A*u ======================
    // d_u holds authoritative OWNED-DOF values; ghost-DOF slots are filled by the
    // forward halo each apply. d_out is the result (owned rows complete after the
    // reverse-add). NOTE: result is d_out, NOT d_y -- d_y is the node-Y COORDINATE
    // vector (domain.getNodeY(), line above) and the two MUST stay distinct, or the
    // A*linear gate would build its field from the result buffer instead of coords.
    cstone::DeviceVector<RealType> d_u(numTotalDofs, 0.0), d_out(numTotalDofs, 0.0);

    // The loop holds d_u constant across calls, so the leading forward-exchange each
    // call refreshes ghosts from the unchanged owned u and a trailing forward is not
    // needed here. A real CG that feeds y back as the next u MUST add a forward
    // exchange on the new u (or a trailing one here) before the next apply, or the
    // ghost slots stay stale -- see the NS conservative-scatter pattern in the header.
    auto matvec = [&]() {
        // 1. FORWARD: ghost-DOF slots <- owners' u (no-op when numRanks==1).
        domain.exchangeNodeHalo(d_u, d_nodeToDof.data());
        // 2. zero y (additive atomic scatter).
        cudaMemset(d_out.data(), 0, numTotalDofs * sizeof(RealType));
        // 3. per-element apply+scatter over OWNED elements only.
        if (numLocal > 0)
            ho_cvfem_apply_launch<RealType, 1>(d_u.data(), d_out.data(), d_elemDof.data(),
                                               d_G.data(), numLocal);
        // 4. REVERSE-ADD: ghost-DOF partials summed into owners (no-op when 1 rank).
        //    After this, only d_out[0,numDofs) is authoritative; ghost slots
        //    [numDofs,numTotalDofs) hold this rank's stale partials (re-published
        //    only by a follow-up forward exchange, which this driver does not need).
        domain.reverseExchangeNodeHaloAdd(d_out, d_nodeToDof.data());
    };

    // ====================== --overlap setup (once, mesh-static) ======================
    // Hide the FORWARD halo exchange behind the interior apply. We replicate
    // exchangeNodeHalo's pack/Isend/Irecv/Waitall/unpack split in two (post / wait)
    // DRIVER-SIDE using the exposed NodeHaloTopology, leaving the shared domain code
    // (used by the NS solvers) untouched. Interior elements (no recv-ghost input)
    // are applied while the exchange is in flight; boundary elements after the wait.
    //
    // Single-rank guard: numRanks==1 has an empty/peerless topology (and fetching it
    // can segfault, see Gate-2 note), so --overlap falls back to the plain full-range
    // apply -- a no-op equal to the blocking path.
    cudaStream_t applyStream = nullptr, packStream = nullptr;
    cstone::DeviceVector<int> d_interiorElems, d_boundaryElems;
    size_t nInterior = 0, nBoundary = 0;
    // Raw topology pointers/sizes cached once; replicate the exchange members.
    const int* h_sendOffsets = nullptr; const int* h_recvOffsets = nullptr;
    const int* d_sendNodeIds = nullptr; const int* d_recvNodeIds = nullptr;
    RealType*  d_sendBuf = nullptr;     RealType*  d_recvBuf = nullptr;
    const std::vector<int>* peersPtr = nullptr;
    int* epochPtr = nullptr;             // mutable epoch_ -- bump exactly once/forward
    std::vector<MPI_Request> fwdReqs;    // live across the interior apply window

    if (overlap && numRanks > 1) {
        CK(cudaStreamCreate(&applyStream));
        CK(cudaStreamCreate(&packStream));
        const auto& topo = domain.getNodeHaloTopology();
        peersPtr      = &topo.peers_;
        h_sendOffsets = topo.sendOffsets_.data();
        h_recvOffsets = topo.recvOffsets_.data();
        d_sendNodeIds = thrust::raw_pointer_cast(topo.sendNodeIds_.data());
        d_recvNodeIds = thrust::raw_pointer_cast(topo.recvNodeIds_.data());
        d_sendBuf     = thrust::raw_pointer_cast(topo.sendBuf_.data());
        d_recvBuf     = thrust::raw_pointer_cast(topo.recvBuf_.data());
        epochPtr      = const_cast<int*>(&topo.epoch_);

        // recv-ghost set: mark every LOCAL node id this rank receives (the INPUT
        // dependency set -- NOT the send list used by the A*linear interior mask).
        size_t recvTotal = topo.recvOffsets_.empty() ? 0 : (size_t)topo.recvOffsets_.back();
        cstone::DeviceVector<uint8_t> d_isRecvGhost(nodeCount, (uint8_t)0);
        if (recvTotal > 0) {
            uint8_t* rg = thrust::raw_pointer_cast(d_isRecvGhost.data());
            const int* rn = d_recvNodeIds;
            thrust::for_each(thrust::device,
                thrust::counting_iterator<size_t>(0),
                thrust::counting_iterator<size_t>(recvTotal),
                [rg, rn] __device__ (size_t i) { rg[rn[i]] = 1; });
            cudaDeviceSynchronize();
        }
        // classify owned elements: boundary iff any corner is a recv-ghost node.
        cstone::DeviceVector<uint8_t> d_elemIsBoundary(numLocal > 0 ? numLocal : 1, (uint8_t)0);
        if (numLocal > 0) {
            int blk = 256, grid = (int)((numLocal + blk - 1) / blk);
            classify_owned_elem_kernel<KeyType><<<grid, blk>>>(
                C(0),C(1),C(2),C(3),C(4),C(5),C(6),C(7),
                thrust::raw_pointer_cast(d_isRecvGhost.data()), startElem, numLocal,
                thrust::raw_pointer_cast(d_elemIsBoundary.data()));
            CK(cudaGetLastError()); cudaDeviceSynchronize();
        }
        // partition [0,numLocal) into the two index lists. copy_if over the predicate
        // and its negation is EXHAUSTIVE: every le lands in exactly one list.
        d_interiorElems.resize(numLocal > 0 ? numLocal : 1);
        d_boundaryElems.resize(numLocal > 0 ? numLocal : 1);
        if (numLocal > 0) {
            const uint8_t* flag = thrust::raw_pointer_cast(d_elemIsBoundary.data());
            auto begin = thrust::counting_iterator<int>(0);
            auto end   = thrust::counting_iterator<int>((int)numLocal);
            auto intEnd = thrust::copy_if(thrust::device, begin, end,
                thrust::device_pointer_cast(d_interiorElems.data()),
                [flag] __device__ (int le) { return flag[le] == 0; });
            auto bndEnd = thrust::copy_if(thrust::device, begin, end,
                thrust::device_pointer_cast(d_boundaryElems.data()),
                [flag] __device__ (int le) { return flag[le] != 0; });
            nInterior = (size_t)(intEnd - thrust::device_pointer_cast(d_interiorElems.data()));
            nBoundary = (size_t)(bndEnd - thrust::device_pointer_cast(d_boundaryElems.data()));
            cudaDeviceSynchronize();
        }
        if (nInterior + nBoundary != numLocal) {
            printf("FATAL rank %d: overlap split not exhaustive (%zu+%zu != %zu)\n",
                   rank, nInterior, nBoundary, numLocal);
            MPI_Abort(MPI_COMM_WORLD, 1);
        }
    }

    // postForward: pack the send buffer on packStream, sync ONLY that stream (cheaper
    // than a full device sync, lets the interior apply on applyStream proceed), then
    // post Irecv/Isend for every peer with nonzero count. Mirrors exchangeNodeHalo
    // (domain.hpp:882) exactly, incl. the tag base 0x4d52 + epoch_.
    auto postForward = [&]() {
        // INVARIANT: d_u OWNED slots must have no in-flight writer when this pack
        // runs. Guaranteed here by the per-iteration trailing applyStream sync +
        // reverseExchangeNodeHaloAdd device sync (matvecOverlap step e), so syncing
        // ONLY packStream below is sufficient. A future caller that mutates d_u on a
        // different stream must sync that stream before calling postForward.
        size_t sendTotal = (peersPtr->empty()) ? 0 : (size_t)h_sendOffsets[peersPtr->size()];
        if (sendTotal > 0) {
            RealType* arr = thrust::raw_pointer_cast(d_u.data());
            const int* snds = d_sendNodeIds; RealType* sbuf = d_sendBuf;
            const int* n2d = d_nodeToDof.data();
            thrust::for_each(thrust::cuda::par.on(packStream),
                thrust::counting_iterator<size_t>(0),
                thrust::counting_iterator<size_t>(sendTotal),
                [arr, snds, sbuf, n2d] __device__ (size_t i) {
                    int n = snds[i];
                    int idx = (n2d != nullptr) ? n2d[n] : n;
                    sbuf[i] = (idx >= 0) ? arr[idx] : arr[0];
                });
            CK(cudaStreamSynchronize(packStream));   // send buffer fully packed before Isend
        }
        constexpr int nodeHaloTagBase = 0x4d52;   // forward tag, identical to domain.hpp:926
        const int tag = nodeHaloTagBase + *epochPtr;
        fwdReqs.clear();
        fwdReqs.reserve(2 * peersPtr->size());
        for (size_t p = 0; p < peersPtr->size(); ++p) {
            int peer = (*peersPtr)[p];
            int rcnt = h_recvOffsets[p+1] - h_recvOffsets[p];
            if (rcnt > 0) {
                MPI_Request r;
                MPI_Irecv(d_recvBuf + h_recvOffsets[p], rcnt, MPI_DOUBLE, peer,
                          tag, MPI_COMM_WORLD, &r);
                fwdReqs.push_back(r);
            }
        }
        for (size_t p = 0; p < peersPtr->size(); ++p) {
            int peer = (*peersPtr)[p];
            int scnt = h_sendOffsets[p+1] - h_sendOffsets[p];
            if (scnt > 0) {
                MPI_Request r;
                MPI_Isend(d_sendBuf + h_sendOffsets[p], scnt, MPI_DOUBLE, peer,
                          tag, MPI_COMM_WORLD, &r);
                fwdReqs.push_back(r);
            }
        }
    };

    // waitUnpackForward: drain the posted reqs, bump epoch_ exactly once (tag-symmetry
    // with the separate reverse base 0x4d53), then scatter recvBuf into the ghost u
    // slots on applyStream so the unpack is ORDERED before the boundary apply.
    auto waitUnpackForward = [&]() {
        if (!fwdReqs.empty())
            MPI_Waitall((int)fwdReqs.size(), fwdReqs.data(), MPI_STATUSES_IGNORE);
        ++(*epochPtr);
        size_t recvTotal = (peersPtr->empty()) ? 0 : (size_t)h_recvOffsets[peersPtr->size()];
        if (recvTotal > 0) {
            RealType* arr = thrust::raw_pointer_cast(d_u.data());
            const int* rnds = d_recvNodeIds; RealType* rbuf = d_recvBuf;
            const int* n2d = d_nodeToDof.data();
            thrust::for_each(thrust::cuda::par.on(applyStream),
                thrust::counting_iterator<size_t>(0),
                thrust::counting_iterator<size_t>(recvTotal),
                [arr, rnds, rbuf, n2d] __device__ (size_t i) {
                    int n = rnds[i];
                    int idx = (n2d != nullptr) ? n2d[n] : n;
                    if (idx >= 0) arr[idx] = rbuf[i];
                });
        }
    };

    // Overlapped matvec. SINGLE-IN-FLIGHT invariant: sendBuf_/recvBuf_ are shared
    // staging buffers reused (roles swapped) by reverseExchangeNodeHaloAdd, so the
    // forward MUST fully drain (waitUnpackForward) before the reverse-add starts.
    // Result is bit-identical to the blocking matvec up to atomicAdd FP reordering
    // (interior+boundary scatter in a different global order than one full launch).
    auto matvecOverlap = [&]() {
        if (!(overlap && numRanks > 1)) { matvec(); return; }
        // (a) post async forward (pack on packStream, Isend/Irecv, return).
        postForward();
        // (b) zero y + apply INTERIOR elements on applyStream (overlaps in-flight MPI).
        CK(cudaMemsetAsync(d_out.data(), 0, numTotalDofs * sizeof(RealType), applyStream));
        if (nInterior > 0)
            ho_cvfem_apply_launch<RealType, 1>(d_u.data(), d_out.data(), d_elemDof.data(),
                d_G.data(), numLocal, applyStream,
                thrust::raw_pointer_cast(d_interiorElems.data()), nInterior);
        // (c) wait + unpack ghost u (unpack enqueued on applyStream, after interior).
        waitUnpackForward();
        // (d) apply BOUNDARY elements on applyStream (ghost u now present + ordered).
        if (nBoundary > 0)
            ho_cvfem_apply_launch<RealType, 1>(d_u.data(), d_out.data(), d_elemDof.data(),
                d_G.data(), numLocal, applyStream,
                thrust::raw_pointer_cast(d_boundaryElems.data()), nBoundary);
        CK(cudaStreamSynchronize(applyStream));   // y complete before reverse-add reads it
        // (e) reverse-add (blocking as today; overlaps naturally inside a CG solve).
        domain.reverseExchangeNodeHaloAdd(d_out, d_nodeToDof.data());
    };

    // ====================== PARITY GATE 1: A*1 = 0 ======================
    // Constant null space of pure-Neumann diffusion. Holds on owned rows on ANY
    // rank count; this is the first proof the forward/reverse halo pair is
    // conservative for the pure matvec (it was only exercised in the NS solver).
    thrust::fill(thrust::device_pointer_cast(d_u.data()),
                 thrust::device_pointer_cast(d_u.data() + numTotalDofs), RealType(1));
    matvec(); CK(cudaGetLastError()); cudaDeviceSynchronize();
    std::vector<RealType> yconst(numDofs);
    if (numDofs > 0)
        CK(cudaMemcpy(yconst.data(), d_out.data(), numDofs * sizeof(RealType), cudaMemcpyDeviceToHost));
    double locMaxConst = 0;
    for (int i = 0; i < numDofs; ++i) locMaxConst = std::max(locMaxConst, std::abs((double)yconst[i]));
    double gMaxConst = 0;
    MPI_Allreduce(&locMaxConst, &gMaxConst, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);

    // ====================== PARITY GATE 2: A*linear = 0 (interior) ======================
    // A linear field has constant gradient; over a CLOSED sub-control-volume the net
    // flux is zero (divergence theorem), so the CVFEM diffusion apply annihilates it
    // on genuine interior rows ONLY. We exclude (a) ghosts, (b) geometric
    // domain-boundary nodes -- their control volume is open and the bare apply uses
    // natural zero-flux, leaving an unbalanced linear flux -> A*linear!=0, and (c)
    // partition-seam nodes -- stencil split across ranks. ownership alone cannot do
    // this (it is owned/ghost only), so we use the global box + the halo send list.
    // First the GLOBAL bounding box (must match the mesh's; for procedural cube it is
    // [0,1]^3, computed here so the .exo path works too).
    RealType locLo[3], locHi[3], gLo[3], gHi[3];
    if (nodeCount > 0) {
        auto px = thrust::device_pointer_cast(d_x.data());
        auto py = thrust::device_pointer_cast(d_y.data());
        auto pz = thrust::device_pointer_cast(d_z.data());
        locLo[0] = *thrust::min_element(px, px + nodeCount);
        locHi[0] = *thrust::max_element(px, px + nodeCount);
        locLo[1] = *thrust::min_element(py, py + nodeCount);
        locHi[1] = *thrust::max_element(py, py + nodeCount);
        locLo[2] = *thrust::min_element(pz, pz + nodeCount);
        locHi[2] = *thrust::max_element(pz, pz + nodeCount);
    } else {
        for (int d = 0; d < 3; ++d) { locLo[d] =  1e300; locHi[d] = -1e300; }
    }
    MPI_Allreduce(locLo, gLo, 3, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
    MPI_Allreduce(locHi, gHi, 3, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
    // tol relative to the global extent so the face test is scale-independent.
    double ext = std::max({gHi[0] - gLo[0], gHi[1] - gLo[1], gHi[2] - gLo[2], 1e-300});
    RealType bndTol = (RealType)(1e-6 * ext);

    cstone::DeviceVector<uint8_t> d_isInterior(numDofs > 0 ? numDofs : 1, 0);
    {
        int blk = 256, grid = (int)((nodeCount + blk - 1) / blk);
        mark_interior_owned_kernel<RealType><<<grid, blk>>>(
            d_nodeOwnership.data(), d_nodeToDof.data(),
            d_x.data(), d_y.data(), d_z.data(), nodeCount, numDofs,
            gLo[0], gLo[1], gLo[2], gHi[0], gHi[1], gHi[2], bndTol, d_isInterior.data());
        CK(cudaGetLastError()); cudaDeviceSynchronize();
    }
    // step 3: un-mark partition-seam nodes (the halo send list). Single-rank has
    // no seams AND no halo topology -- skip entirely (forcing getNodeHaloTopology()
    // at numRanks==1 builds an empty/peerless topology and segfaults).
    if (numRanks > 1) {
        const auto& topo = domain.getNodeHaloTopology();
        int sendCount = topo.sendOffsets_.empty() ? 0 : (int)topo.sendOffsets_.back();
        if (sendCount > 0) {
            int blk = 256, grid = (sendCount + blk - 1) / blk;
            unmark_seam_kernel<<<grid, blk>>>(
                thrust::raw_pointer_cast(topo.sendNodeIds_.data()), sendCount,
                d_nodeToDof.data(), numDofs, d_isInterior.data());
            CK(cudaGetLastError()); cudaDeviceSynchronize();
        }
    }
    {
        int blk = 256, grid = (int)((nodeCount + blk - 1) / blk);
        // arbitrary non-degenerate linear field; coords identical across partitions.
        fill_linear_field_kernel<RealType><<<grid, blk>>>(
            d_x.data(), d_y.data(), d_z.data(), d_nodeToDof.data(), nodeCount,
            RealType(0.37), RealType(-1.11), RealType(0.53), d_u.data());
        CK(cudaGetLastError()); cudaDeviceSynchronize();
    }
    matvec(); CK(cudaGetLastError()); cudaDeviceSynchronize();
    std::vector<RealType> ylin(numDofs);
    std::vector<uint8_t>  isInt(numDofs);
    if (numDofs > 0) {
        CK(cudaMemcpy(ylin.data(), d_out.data(), numDofs * sizeof(RealType), cudaMemcpyDeviceToHost));
        CK(cudaMemcpy(isInt.data(), d_isInterior.data(), numDofs * sizeof(uint8_t), cudaMemcpyDeviceToHost));
    }
    double locMaxLin = 0; long long locInterior = 0;
    for (int i = 0; i < numDofs; ++i)
        if (isInt[i]) { locMaxLin = std::max(locMaxLin, std::abs((double)ylin[i])); ++locInterior; }
    double gMaxLin = 0; long long gInterior = 0;
    MPI_Allreduce(&locMaxLin, &gMaxLin, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
    MPI_Allreduce(&locInterior, &gInterior, 1, MPI_LONG_LONG, MPI_SUM, MPI_COMM_WORLD);

    // ====================== PARITY GATE 3: overlap == blocking ======================
    // Run BOTH paths on the SAME random u and assert max|overlap-blocking|/max|blocking|
    // is at round-off. NOT bitwise-equal: interior+boundary scatter in a different
    // global order than the single full-range launch, so the only difference is
    // atomicAdd FP reordering (already non-deterministic in the baseline) -> ~1e-13
    // relative, well under 1e-12. Skipped when --overlap is off or single-rank (the
    // overlap path is then a no-op alias of blocking and parity is trivially 0).
    // NOTE: single-rank aliases the blocking path (matvecOverlap line ~478), so the
    // subset-kernel split is ONLY exercised at numRanks>1 -- this gate must run
    // multi-rank (e.g. -n 8) to actually cover the interior/boundary partition.
    double gRelDiff = 0.0;
    if (overlap && numRanks > 1) {
        // reproducible per-node field so u is identical to what every rank computes.
        {
            int blk = 256, grid = (int)((nodeCount + blk - 1) / blk);
            fill_linear_field_kernel<RealType><<<grid, blk>>>(
                d_x.data(), d_y.data(), d_z.data(), d_nodeToDof.data(), nodeCount,
                RealType(0.91), RealType(0.13), RealType(-0.47), d_u.data());
            CK(cudaGetLastError()); cudaDeviceSynchronize();
        }
        // blocking reference.
        matvec(); CK(cudaGetLastError()); cudaDeviceSynchronize();
        std::vector<RealType> yBlock(numDofs > 0 ? numDofs : 1);
        if (numDofs > 0)
            CK(cudaMemcpy(yBlock.data(), d_out.data(), numDofs * sizeof(RealType), cudaMemcpyDeviceToHost));
        // re-fill u (the blocking matvec leaves ghost slots overwritten, but u OWNED
        // slots are unchanged; refill anyway so the overlap path sees the same input).
        {
            int blk = 256, grid = (int)((nodeCount + blk - 1) / blk);
            fill_linear_field_kernel<RealType><<<grid, blk>>>(
                d_x.data(), d_y.data(), d_z.data(), d_nodeToDof.data(), nodeCount,
                RealType(0.91), RealType(0.13), RealType(-0.47), d_u.data());
            CK(cudaGetLastError()); cudaDeviceSynchronize();
        }
        matvecOverlap(); CK(cudaGetLastError()); cudaDeviceSynchronize();
        std::vector<RealType> yOver(numDofs > 0 ? numDofs : 1);
        if (numDofs > 0)
            CK(cudaMemcpy(yOver.data(), d_out.data(), numDofs * sizeof(RealType), cudaMemcpyDeviceToHost));
        double locMaxDiff = 0, locMaxRef = 0;
        for (int i = 0; i < numDofs; ++i) {
            locMaxDiff = std::max(locMaxDiff, std::abs((double)(yOver[i] - yBlock[i])));
            locMaxRef  = std::max(locMaxRef,  std::abs((double)yBlock[i]));
        }
        double gMaxDiff = 0, gMaxRef = 0;
        MPI_Allreduce(&locMaxDiff, &gMaxDiff, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
        MPI_Allreduce(&locMaxRef,  &gMaxRef,  1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
        gRelDiff = (gMaxRef > 0) ? gMaxDiff / gMaxRef : gMaxDiff;
    }

    // ====================== OPTIONAL single-vs-multi parity dump ======================
    // Write (globalSfcKey, y) for every owned DOF. The SFC key is a stable global
    // node id, so a 1-rank and an N-rank run on the SAME global mesh produce the
    // same (key,value) set -- diff them externally for bit-exact parity. We do a
    // reproducible coordinate-keyed field again so the input is identical.
    if (!dumpFile.empty()) {
        {
            int blk = 256, grid = (int)((nodeCount + blk - 1) / blk);
            fill_linear_field_kernel<RealType><<<grid, blk>>>(
                d_x.data(), d_y.data(), d_z.data(), d_nodeToDof.data(), nodeCount,
                RealType(1.0), RealType(2.0), RealType(3.0), d_u.data());
            CK(cudaGetLastError()); cudaDeviceSynchronize();
        }
        matvec(); CK(cudaGetLastError()); cudaDeviceSynchronize();

        const auto& d_sfc = domain.getLocalToGlobalSfcMap();   // size nodeCount
        std::vector<KeyType> h_sfc(nodeCount);
        std::vector<int>     h_dof(nodeCount);
        std::vector<RealType> h_y(numDofs > 0 ? numDofs : 1);
        CK(cudaMemcpy(h_sfc.data(), d_sfc.data(), nodeCount * sizeof(KeyType), cudaMemcpyDeviceToHost));
        CK(cudaMemcpy(h_dof.data(), d_nodeToDof.data(), nodeCount * sizeof(int), cudaMemcpyDeviceToHost));
        if (numDofs > 0)
            CK(cudaMemcpy(h_y.data(), d_out.data(), numDofs * sizeof(RealType), cudaMemcpyDeviceToHost));
        // owned node -> (key, y[dof]); rank-suffixed file so ranks don't clobber.
        std::string fn = dumpFile + ".r" + std::to_string(rank);
        FILE* f = std::fopen(fn.c_str(), "w");
        if (f) {
            for (size_t n = 0; n < nodeCount; ++n) {
                int dof = h_dof[n];
                if (dof >= 0 && dof < numDofs)
                    std::fprintf(f, "%llu %.17g\n", (unsigned long long)h_sfc[n], (double)h_y[dof]);
            }
            std::fclose(f);
        }
    }

    // ====================== timed weak-scaling loop ======================
    // Time the WHOLE matvec (forward halo + apply + reverse-add) as one unit, plus
    // a compute-only number (apply alone) so comm overhead is visible. The halo
    // calls each cudaDeviceSynchronize + MPI_Waitall internally (blocking), so the
    // per-iter wall time already folds comm in. We barrier before/after so the
    // reported max-over-ranks time reflects the true synchronized cost.
    cstone::DeviceVector<RealType> d_uTime(numTotalDofs, 1.0);
    thrust::copy(thrust::device_pointer_cast(d_uTime.data()),
                 thrust::device_pointer_cast(d_uTime.data() + numTotalDofs),
                 thrust::device_pointer_cast(d_u.data()));
    matvec(); cudaDeviceSynchronize();   // warmup (build halo topology, JIT, etc.)

    MPI_Barrier(MPI_COMM_WORLD);
    cudaEvent_t t0, t1; cudaEventCreate(&t0); cudaEventCreate(&t1);

    // full matvec (compute + comm). We take BOTH a cudaEvent number (stream timeline)
    // and an MPI_Wtime wall-clock number. The halo calls block on the host
    // (cudaDeviceSync + MPI_Waitall), so the cudaEvent capture of comm is fragile;
    // the MPI_Wtime number is the honest weak-scaling figure and is what we report.
    MPI_Barrier(MPI_COMM_WORLD);
    double wt0 = MPI_Wtime();
    cudaEventRecord(t0);
    for (int it = 0; it < iters; ++it) matvec();
    cudaEventRecord(t1); cudaEventSynchronize(t1);
    double msFullWall = (MPI_Wtime() - wt0) * 1e3 / iters;
    float msFull = 0; cudaEventElapsedTime(&msFull, t0, t1); msFull /= iters;

    // compute-only: forward halo once (ghosts valid), then time apply alone.
    domain.exchangeNodeHalo(d_u, d_nodeToDof.data()); cudaDeviceSynchronize();
    cudaEventRecord(t0);
    for (int it = 0; it < iters; ++it) {
        cudaMemset(d_out.data(), 0, numTotalDofs * sizeof(RealType));
        if (numLocal > 0)
            ho_cvfem_apply_launch<RealType, 1>(d_u.data(), d_out.data(), d_elemDof.data(),
                                               d_G.data(), numLocal);
    }
    cudaEventRecord(t1); cudaEventSynchronize(t1);
    float msApply = 0; cudaEventElapsedTime(&msApply, t0, t1); msApply /= iters;
    CK(cudaGetLastError());
    MPI_Barrier(MPI_COMM_WORLD);

    // --- OVERLAP matvec timing (wall clock). Same methodology as the full blocking
    //     matvec: barrier, MPI_Wtime over the loop. The overlap hides the FORWARD
    //     exchange behind the interior apply; the reverse-add stays exposed (it
    //     cannot be hidden in a standalone matvec, but overlaps inside a CG). ---
    double msOverlapWall = 0.0;
    if (overlap && numRanks > 1) {
        // restore u (the parity check left a linear field; timing wants the all-1s
        // input, matching the blocking timer for an apples-to-apples comparison).
        thrust::copy(thrust::device_pointer_cast(d_uTime.data()),
                     thrust::device_pointer_cast(d_uTime.data() + numTotalDofs),
                     thrust::device_pointer_cast(d_u.data()));
        matvecOverlap(); cudaDeviceSynchronize();   // warmup
        MPI_Barrier(MPI_COMM_WORLD);
        double wo0 = MPI_Wtime();
        for (int it = 0; it < iters; ++it) matvecOverlap();
        cudaDeviceSynchronize();
        msOverlapWall = (MPI_Wtime() - wo0) * 1e3 / iters;
        MPI_Barrier(MPI_COMM_WORLD);
    }

    // --- global totals + slowest-rank times (the weak-scaling metric) ---
    unsigned long long locElems = numLocal, gElems = 0;
    unsigned long long locDofs  = (unsigned long long)numDofs, gDofs = 0;
    MPI_Reduce(&locElems, &gElems, 1, MPI_UNSIGNED_LONG_LONG, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(&locDofs,  &gDofs,  1, MPI_UNSIGNED_LONG_LONG, MPI_SUM, 0, MPI_COMM_WORLD);
    // Per-rank owned-element min/max so an imbalanced or wrong invocation (same
    // --ncells on 1 and N ranks -> per-rank work drops 1/N) is caught instead of
    // silently inflating efficiency. A held ladder must have min~max~target/rank.
    unsigned long long minElems = 0, maxElems = 0;
    MPI_Reduce(&locElems, &minElems, 1, MPI_UNSIGNED_LONG_LONG, MPI_MIN, 0, MPI_COMM_WORLD);
    MPI_Reduce(&locElems, &maxElems, 1, MPI_UNSIGNED_LONG_LONG, MPI_MAX, 0, MPI_COMM_WORLD);
    float maxFull = 0, maxApply = 0; double maxFullWall = 0, maxOverlapWall = 0;
    MPI_Reduce(&msFull,        &maxFull,        1, MPI_FLOAT,  MPI_MAX, 0, MPI_COMM_WORLD);
    MPI_Reduce(&msApply,       &maxApply,       1, MPI_FLOAT,  MPI_MAX, 0, MPI_COMM_WORLD);
    MPI_Reduce(&msFullWall,    &maxFullWall,    1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
    MPI_Reduce(&msOverlapWall, &maxOverlapWall, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);

    // --- parallel efficiency vs 1-rank baseline. Each rank holds ~equal work, so
    //     the weak-scaling baseline is the per-rank MDOF/s on 1 rank. We report
    //     per-rank-equivalent throughput (gDofs/numRanks / time) so efficiency is
    //     (perRankThroughput_N / perRankThroughput_1). The 1-rank run prints its
    //     own per-rank number; the user reads efficiency across the ladder runs.
    if (rank == 0) {
        // Wall clock (MPI_Wtime) is the reported weak-scaling figure; cudaEvents are
        // a compute-side cross-check only.
        double fullDofPerS = (double)gDofs / (maxFullWall * 1e-3);  // total throughput
        double appDofPerS  = (double)gDofs / (maxApply    * 1e-3);
        double perRankFull = fullDofPerS / numRanks;               // weak-scaling unit
        printf("HO CVFEM matrix-free DISTRIBUTED weak-scaling apply (Knaus Alg 2, p=1)\n");
        printf("ranks=%d | global elems=%llu  global DOFs=%llu | ~%.2fM elems/rank\n",
               numRanks, gElems, gDofs, (double)gElems / numRanks / 1e6);
        printf("per-rank owned elems: min=%llu max=%llu  [%s -- must be ~constant for a held ladder]\n",
               minElems, maxElems,
               (maxElems > 0 && (maxElems - minElems) <= maxElems / 20) ? "balanced" : "IMBALANCED");
        printf("\n-- PARITY (partition-invariant gates) --\n");
        printf("  A*1=0           : max|y| = %.3e   [%s]\n", gMaxConst,
               gMaxConst < 1e-9 ? "PASS" : "FAIL");
        printf("  A*linear=0 (int): max|y| = %.3e over %lld interior owned rows   [%s]\n",
               gMaxLin, gInterior, gMaxLin < 1e-8 ? "PASS" : "FAIL");
        if (overlap && numRanks > 1)
            printf("  overlap==block  : rel diff = %.3e   [%s]\n", gRelDiff,
                   gRelDiff < 1e-12 ? "PASS" : "FAIL");
        printf("\n-- THROUGHPUT (slowest rank) --\n");
        printf("  full matvec  : %.4f ms wall | %.1f MDOF/s total | %.1f MDOF/s/rank\n",
               maxFullWall, fullDofPerS / 1e6, perRankFull / 1e6);
        printf("  (cudaEvent full matvec: %.4f ms -- compute-side cross-check)\n", maxFull);
        printf("  apply only   : %.4f ms | %.1f MDOF/s total\n", maxApply, appDofPerS / 1e6);
        printf("  comm fraction: %.1f%% (wall_full - apply)/wall_full\n",
               maxFullWall > 0 ? 100.0 * (maxFullWall - maxApply) / maxFullWall : 0.0);
        if (overlap && numRanks > 1) {
            double overDofPerS = (double)gDofs / (maxOverlapWall * 1e-3);
            printf("  -- OVERLAP (forward halo hidden behind interior apply) --\n");
            printf("  overlap matvec: %.4f ms wall | %.1f MDOF/s total | %.1f MDOF/s/rank\n",
                   maxOverlapWall, overDofPerS / 1e6, overDofPerS / numRanks / 1e6);
            printf("  speedup vs block: %.2fx | residual comm fraction: %.1f%% (reverse-add exposed)\n",
                   maxOverlapWall > 0 ? maxFullWall / maxOverlapWall : 0.0,
                   maxOverlapWall > 0 ? 100.0 * (maxOverlapWall - maxApply) / maxOverlapWall : 0.0);
        }
        printf("\n  WEAK-SCALING: record this per-rank MDOF/s across the ladder.\n");
        printf("  efficiency_N = perRankMDOFs(N) / perRankMDOFs(1). Run -n 1,8,32\n");
        printf("  holding ~2M elems/rank (procedural --ncells so N^3/ranks ~ 2M).\n");
    }

    cudaEventDestroy(t0); cudaEventDestroy(t1);
    if (applyStream) cudaStreamDestroy(applyStream);
    if (packStream)  cudaStreamDestroy(packStream);
    delete domainPtr;
    MPI_Finalize();
    return 0;
}
