#ifndef MARS_NONCONFORMAL_HPP
#define MARS_NONCONFORMAL_HPP

// GPU-native non-conformal (DG-IP) CVFEM interface -- the MARS-native mirror of Sierra/Nalu's
// NonConformalManager/DgInfo, per Domino JCP 359 (2018) 331-351. Couples two mesh surfaces that are
// geometrically coincident but whose nodes/DOFs do not match, via the DG interior-penalty flux
//   T_hat(a,b) = 1/2 (q^a.n^a - q^b.n^b) + lambda (phi^a - phi^b),   lambda = 1/2 (Gamma^a/L^a + Gamma^b/L^b)
// assembled DUAL-PASS with the FULL implicit cross-interface stencil (dropping the opposing columns
// diverges -- Nalu's hard-won fact). Design doc: internal-notes/NONCONFORMAL_INTERFACE_DESIGN.md.
//
// MARS-native wins over the STK/Trilinos path (design doc S8):
//  - coarse search = cstone traversal/collisions.hpp (box-vs-octree overlap) ON-DEVICE, reusing the
//    domain octree we already build -- replaces stk::search::coarse_search (host-side in Nalu).
//  - ghosting = the node-halo machinery + an interface ghost list -- replaces STK custom ghosting.
//  - fine search = closed-form barycentric at p=1 (planar tet faces) -- no Newton.
//
// STAGING (this file starts P0; P1-P4 add search/ghost/momentum-continuity/sliding):
//  P0 [here]  scalar DG-IP flux + identity pairing (case c: matching mesh, discontinuous DOFs). No search,
//             no ghosting. GATE: MMS design order O(h^2) across the interface.
//  P1         cstone-collision coarse search + barycentric fine search -> non-matching static blocks.
//  P2         MPI cross-rank opposing elements via the interface ghost list.
//  P3         momentum + continuity (C_hat with PSPG inside the interface mdot -- the pump pressure system).
//  P4         sliding rotor-stator (re-pair on-device each step).

#include <string>
#include <cstdint>
#include <thrust/device_vector.h>
#include "backend/distributed/unstructured/domain.hpp"

namespace mars {
namespace fem {

// Per-integration-point pairing record (mirrors Nalu DgInfo, device SoA). One entry per MASTER-side
// sub-control-surface integration point on the interface. For P0 (matching mesh) the pairing is identity
// (opposingNode == the coincident master node), so the search-derived fields are trivially filled.
template<typename RealType>
struct DgInfoSoA {
    thrust::device_vector<int>      currentElem;    // master element (local id)
    thrust::device_vector<int>      currentFace;    // master face ordinal within the element
    thrust::device_vector<int>      opposingElem;   // opposing (slave) element; == current for identity pairing
    thrust::device_vector<RealType> curIsoCoords;   // master-side isoparametric/barycentric coords (nDim per ip)
    thrust::device_vector<RealType> oppIsoCoords;   // opposing-side coords (== cur for matching mesh)
    thrust::device_vector<RealType> areaVec;        // interface SCS area vector at the ip (nDim per ip)
    thrust::device_vector<uint8_t>  opposingGhosted;// 1 = opposing element is off-rank (P2+); 0 for P0/P1 local
    size_t numIps = 0;

    void resize(size_t n, int nDim) {
        numIps = n;
        currentElem.resize(n); currentFace.resize(n); opposingElem.resize(n);
        curIsoCoords.resize(n * nDim); oppIsoCoords.resize(n * nDim);
        areaVec.resize(n * nDim); opposingGhosted.assign(n, 0);
    }
};

// The DG interior-penalty scalar flux at one interface ip (Domino Eq. 6-7). Gamma = diffusivity;
// L^a/L^b are the shape-function normal length scales L = (sum_k dN_k/dx_j n_j)^-1 (Eq. 24). NO tunable
// penalty constant. Returned as the numerical flux per unit area; the caller scales by |area| and
// assembles dual-pass. HOST_DEVICE so the P0 kernel and the MMS host oracle share it.
template<typename RealType>
__host__ __device__ inline RealType dgScalarFlux(
    RealType qDotN_a, RealType qDotN_b,   // q^s . n^s = -Gamma^s dphi^s/dx_j n_j^s, each side's own normal
    RealType phi_a,   RealType phi_b,
    RealType Gamma_a, RealType Gamma_b, RealType Linv_a, RealType Linv_b)
{
    RealType lambda = RealType(0.5) * (Gamma_a * Linv_a + Gamma_b * Linv_b);   // Eq. 7
    return RealType(0.5) * (qDotN_a - qDotN_b) + lambda * (phi_a - phi_b);     // Eq. 6
}

// P0 interface object: register two side sets, build the (identity, for matching mesh) pairing + DgInfo,
// and assemble the scalar DG-IP flux into the CSR. P1 replaces buildPairingIdentity() with the cstone
// collision search; the flux kernel + dual-pass assembly are reused unchanged.
template<typename ElementTag, typename RealType, typename KeyType>
class NonConformalInterface {
public:
    using Domain = ElementDomain<ElementTag, RealType, KeyType, cstone::GpuTag>;

    // masterBlock / slaveBlock are the element-block ids of the two interface sides (fluid vs solid,
    // or the two sides of a non-conformal cut). They disambiguate the coincident interface nodes,
    // which share an SFC key but differ in block. Default -1 = single-block (the two sides are already
    // distinct, e.g. matching-mesh discontinuous DOFs); the caller/driver supplies the ids because a
    // side set carries no block of its own (side sets store node keys only).
    NonConformalInterface(Domain& domain, std::string masterSideSet, std::string slaveSideSet,
                          int masterBlock = -1, int slaveBlock = -1)
        : domain_(domain), master_(std::move(masterSideSet)), slave_(std::move(slaveSideSet)),
          masterBlock_(masterBlock), slaveBlock_(slaveBlock) {}

    // P0: matching mesh -> pair each master interface node to the SLAVE node at the same SFC key (coincident
    // coordinates). Reuses the domain's side-set key storage (sideSetNodeKeys) + the SFC map, exactly like
    // ElementDomain::setSideSetFromKeys resolves keys->local. No geometric search. Fills DgInfoSoA identity.
    void buildPairingIdentity();   // TODO(P0): implement against sideSetNodeKeys(master_/slave_) + SFC map

    // P1: cstone traversal/collisions.hpp coarse search (slave query boxes vs the master octree) + closed-
    // form barycentric fine search -> DgInfoSoA for non-matching meshes. P2 flags off-rank opposing as ghost.
    void buildPairingSearch();     // TODO(P1)

    // Assemble the scalar DG-IP flux (dgScalarFlux) into the coupled CSR, dual-pass, with the full implicit
    // cross-interface stencil (master rows carry opposing-element columns and vice versa). Gamma = the scalar
    // diffusivity. rowOff/colInd/vals are the interleaved coupled CSR the pump already builds.
    void assembleScalarFlux(const int* rowOff, const int* colInd, RealType* vals,
                            RealType Gamma, int blockSize);   // TODO(P0): the dual-pass kernel

    const DgInfoSoA<RealType>& dgInfo() const { return dg_; }

private:
    Domain&      domain_;
    std::string  master_, slave_;
    int          masterBlock_ = -1, slaveBlock_ = -1;  // element-block ids of the two interface sides
    DgInfoSoA<RealType> dg_;
};

}  // namespace fem
}  // namespace mars

#endif  // MARS_NONCONFORMAL_HPP
