#pragma once

// Face (sub-control-surface) fluxes for the pump solver's projection, TET4 and HEX8.
//
// The projection keeps a volumetric flux F per (element, sub-control surface), F = u . A on the
// surface between nodes L and R (A from L to R, d_areaVec's layout e * NSCS + ip), and makes it,
// not the nodal velocity, the field continuity holds for:
//
//   F*      = ubar** . A  [ - c (beta_L beta_R)^m ( (grad p^n)^K_f . A - 1/2 (g_L + g_R) . A ) ]
//   K_beta phi = -(rho/dtEff) sum_s sigma_s F*_s                                         pressure
//   F^{n+1} = F* - (dtEff/rho) beta_f (grad phi)^K_f . A                                 corrector
//
// (grad q)^K_f . A is the face gradient the pressure operator K is assembled from, so
// sum_s sigma_s F^{n+1}_s = 0 in every control volume whose pressure row is free, to the
// solver's tolerance:
//   TET4: K_ij = V gamma_e grad N_i . grad N_j, the median-dual face form of P1, with the element
//         gradient and gamma_e the element mean (mars_cvfem_tet_kernel_full.hpp);
//   HEX8: gamma_f grad q(x_ip) . A at the shifted integration point, gamma_f = (gamma_L + gamma_R)/2
//         (mars_cvfem_hex_kernel_tensor_faces.hpp).
// The nodal velocities keep MARS's predictor and corrector. The bracket is an optional Rhie-Chow term
// (c = 0 by default): g the gradient the predictor used, so it swaps exactly what u** carries for the
// compact face gradient. Measured 2026-11-14 on LDP24H/LDP24HX at 1 bar: with K's average beta it
// carried 13-19 % of the delivered volume through the penalised plates in 10 ms (the pressure jumps
// across a plate); weighting it by (beta_L beta_R)^m cuts that, but the projection does not need it.
//
// G (faceFluxNodalGradient, faceFluxElementGradient) is the volume-weighted average of the element
// gradients over each node's sub-control volumes: exact for linear fields, no boundary closure; on
// HEX8 evaluated at the sub-control volumes' centres, so it sees the corner (hourglass) modes. It is
// the reconstruction the convective deferred correction needs. (Used as the predictor's and
// corrector's gradient it made the nodal velocities noisier than MARS's own, and the solid's Newton
// failed on LDP24HX at 7.9 ms.)
//
// The leaf functions (upwind switch, limited increments) follow the SFEM CVFEM spike
// (sfem/spikes/cvfem/src/kernels/microkernels): same formulas, element-agnostic here.

#include "mars_element_traits.hpp"
#include "mars_cvfem_kernel.hpp"
#include "mars_cvfem_hex_kernel_tensor_faces.hpp"
#include "mars_hex_geometry.hpp"

#include <cmath>

namespace mars
{
namespace fem
{

// ---------------------------------------------------------------------------------------------
// Leaf functions (the SFEM spike's cvfem_upwind_abs and limiter arms)

// |m|_eps: |m| for |m| >= eps, m^2 (2 eps - |m|) / eps^2 below: value and slope match at eps,
// zero with zero slope at 0, never above |m|. eps = 0 is the hard switch.
template<typename RealType>
__host__ __device__ inline RealType faceFluxUpwindAbs(RealType m, RealType eps)
{
    const RealType am = m > RealType(0) ? m : -m;
    if (eps > RealType(0) && am < eps) return m * m * (RealType(2) * eps - am) / (eps * eps);
    return am;
}

// The face value base + inc confined to [lo, hi] (the edge's two nodes): returns the increment.
template<typename RealType>
__host__ __device__ inline RealType faceFluxClipInc(RealType base, RealType inc, RealType lo, RealType hi)
{
    RealType f = base + inc;
    f          = f < lo ? lo : (f > hi ? hi : f);
    return f - base;
}

// Venkatakrishnan's limiter on an increment, eps2 dimensional (eps2 = 0: the plain bound).
template<typename RealType>
__host__ __device__ inline RealType faceFluxVenkatakrishnanInc(RealType base, RealType inc, RealType lo, RealType hi,
                                                               RealType eps2)
{
    const RealType dp  = (inc >= RealType(0)) ? (hi - base) : (lo - base);
    const RealType num = dp * dp + RealType(2) * inc * dp + eps2;
    const RealType den = dp * dp + RealType(2) * inc * inc + inc * dp + eps2;
    return (den != RealType(0)) ? (num / den) * inc : inc;
}

// Darwish-Moukalled: a virtual upwind node and van Leer; graddotd = grad phi_C . (x_face - x_C).
template<typename RealType>
__host__ __device__ inline RealType faceFluxDarwishMoukalledInc(RealType phiC, RealType phiD, RealType graddotd)
{
    const RealType a   = RealType(2) * graddotd;
    const RealType b   = phiD - phiC;
    const RealType aa  = a < RealType(0) ? -a : a;
    const RealType ab  = b < RealType(0) ? -b : b;
    const RealType den = RealType(2) * (aa + ab);
    return (den != RealType(0)) ? (a * ab + aa * b) / den : RealType(0);
}

// ---------------------------------------------------------------------------------------------
// Element geometry

template<typename RealType>
__host__ __device__ inline RealType faceFluxInverse3(const RealType J[3][3], RealType inv[3][3])
{
    const RealType det = J[0][0] * (J[1][1] * J[2][2] - J[1][2] * J[2][1]) - J[0][1] * (J[1][0] * J[2][2] - J[1][2] * J[2][0]) +
                         J[0][2] * (J[1][0] * J[2][1] - J[1][1] * J[2][0]);
    const RealType id = RealType(1) / det;
    inv[0][0]         = (J[1][1] * J[2][2] - J[1][2] * J[2][1]) * id;
    inv[0][1]         = (J[0][2] * J[2][1] - J[0][1] * J[2][2]) * id;
    inv[0][2]         = (J[0][1] * J[1][2] - J[0][2] * J[1][1]) * id;
    inv[1][0]         = (J[1][2] * J[2][0] - J[1][0] * J[2][2]) * id;
    inv[1][1]         = (J[0][0] * J[2][2] - J[0][2] * J[2][0]) * id;
    inv[1][2]         = (J[0][2] * J[1][0] - J[0][0] * J[1][2]) * id;
    inv[2][0]         = (J[1][0] * J[2][1] - J[1][1] * J[2][0]) * id;
    inv[2][1]         = (J[0][1] * J[2][0] - J[0][0] * J[2][1]) * id;
    inv[2][2]         = (J[0][0] * J[1][1] - J[0][1] * J[1][0]) * id;
    return det;
}

// Physical gradients of the trilinear basis at (s1, s2, s3) in the tensor kernels' reference cube
// [-1/2, 1/2]^3 (N_n = prod (1/2 +- s), Exodus corner order).
template<typename RealType>
__host__ __device__ inline void faceFluxHexGrad(const RealType x[8][3], RealType s1, RealType s2, RealType s3,
                                                RealType g[8][3])
{
    const RealType f1[2] = {RealType(0.5) - s1, RealType(0.5) + s1};
    const RealType f2[2] = {RealType(0.5) - s2, RealType(0.5) + s2};
    const RealType f3[2] = {RealType(0.5) - s3, RealType(0.5) + s3};
    RealType dN[8][3];
    for (int n = 0; n < 8; ++n)
    {
        const int a = (n & 3) == 1 || (n & 3) == 2, b = (n & 3) >= 2, c = n >= 4;
        const RealType sa = a ? RealType(1) : RealType(-1), sb = b ? RealType(1) : RealType(-1),
                       sc = c ? RealType(1) : RealType(-1);
        dN[n][0] = sa * (f2[b] * f3[c]);
        dN[n][1] = sb * (f1[a] * f3[c]);
        dN[n][2] = sc * (f1[a] * f2[b]);
    }
    RealType J[3][3] = {{0, 0, 0}, {0, 0, 0}, {0, 0, 0}};
    for (int n = 0; n < 8; ++n)
        for (int a = 0; a < 3; ++a)
            for (int b = 0; b < 3; ++b)
                J[a][b] += dN[n][b] * x[n][a];
    RealType inv[3][3];
    faceFluxInverse3(J, inv);
    // grad_x N = J^-T grad_s N
    for (int n = 0; n < 8; ++n)
        for (int a = 0; a < 3; ++a)
            g[n][a] = inv[0][a] * dN[n][0] + inv[1][a] * dN[n][1] + inv[2][a] * dN[n][2];
}

// An element's nodes, coordinates and nodal values of up to two fields.
template<typename KeyType, typename RealType, typename ElementTag>
struct FaceFluxElement
{
    static constexpr int NPE  = ElemTraits<ElementTag>::NodesPerElem;
    static constexpr int NSCS = ElemTraits<ElementTag>::ScsPerElem;
    KeyType n[NPE];
    RealType x[NPE][3];

    __device__ void load(const KeyType* const* cc, size_t e, const RealType* X, const RealType* Y, const RealType* Z)
    {
        for (int i = 0; i < NPE; ++i)
        {
            n[i]    = cc[i][e];
            x[i][0] = X[n[i]];
            x[i][1] = Y[n[i]];
            x[i][2] = Z[n[i]];
        }
    }
};

// An element's contribution to G at each of its nodes, and its weight there. TET4: the constant
// gradient, weight V/4. HEX8: the gradient at the centre of the node's sub-control volume (the octant
// towards the corner, s = +-1/4 in the tensor kernels' cube), weight the sub-control volume. Exact for
// linear fields either way; on the hex, unlike the gradient at the element centre, it sees the corner
// (hourglass) modes of a trilinear field, so G couples the pressure's checkerboard back into the
// velocity as it does on tets -- with the centre gradient that mode grew without the Rhie-Chow bracket.
template<typename KeyType, typename RealType, typename ElementTag>
__device__ inline void faceFluxElementGradient(const FaceFluxElement<KeyType, RealType, ElementTag>& el,
                                               const RealType* phi,
                                               RealType grad[ElemTraits<ElementTag>::NodesPerElem][3],
                                               RealType weight[ElemTraits<ElementTag>::NodesPerElem])
{
    constexpr int NPE = ElemTraits<ElementTag>::NodesPerElem;
    RealType g[NPE][3];
    if constexpr (std::is_same_v<ElementTag, TetTag>)
    {
        RealType det;
        Tet4CVFEM::jacobian_and_dNdx<RealType>(el.x, det, g);
        const RealType quarterV = fabs(det) / RealType(24);
        RealType ge[3] = {0, 0, 0};
        for (int i = 0; i < 4; ++i)
        {
            const RealType p = phi[el.n[i]];
            ge[0] += g[i][0] * p;
            ge[1] += g[i][1] * p;
            ge[2] += g[i][2] * p;
        }
        for (int i = 0; i < 4; ++i)
        {
            weight[i]  = quarterV;
            grad[i][0] = ge[0];
            grad[i][1] = ge[1];
            grad[i][2] = ge[2];
        }
    }
    else
    {
        for (int c = 0; c < 8; ++c)
        {
            const int a = (c & 3) == 1 || (c & 3) == 2, b = (c & 3) >= 2, z = c >= 4;
            const RealType q = RealType(0.25);
            faceFluxHexGrad<RealType>(el.x, a ? q : -q, b ? q : -q, z ? q : -q, g);
            grad[c][0] = grad[c][1] = grad[c][2] = RealType(0);
            for (int n = 0; n < 8; ++n)
            {
                const RealType p = phi[el.n[n]];
                grad[c][0] += g[n][0] * p;
                grad[c][1] += g[n][1] * p;
                grad[c][2] += g[n][2] * p;
            }
            weight[c] = hexSubVolume<RealType>(el.x, c);
        }
    }
}

// K's face gradient (grad phi)^K_f . A on sub-control surface ip, and the face's gamma (beta)
// as K averages it. tetGrad: the element gradient (TET4), computed once by the caller.
template<typename KeyType, typename RealType, typename ElementTag>
__device__ inline RealType faceFluxKGradDotArea(const FaceFluxElement<KeyType, RealType, ElementTag>& el,
                                                const RealType* phi, int ip, const RealType A[3],
                                                const RealType tetGrad[3])
{
    if constexpr (std::is_same_v<ElementTag, TetTag>)
    {
        (void)el;
        (void)phi;
        (void)ip;
        return tetGrad[0] * A[0] + tetGrad[1] * A[1] + tetGrad[2] * A[2];
    }
    else
    {
        (void)tetGrad;
        double s1, s2, s3;
        tensorFacesPoint(ip, s1, s2, s3);
        RealType g[8][3];
        faceFluxHexGrad<RealType>(el.x, RealType(s1), RealType(s2), RealType(s3), g);
        RealType r = 0;
        for (int n = 0; n < 8; ++n) r += phi[el.n[n]] * (g[n][0] * A[0] + g[n][1] * A[1] + g[n][2] * A[2]);
        return r;
    }
}

template<typename KeyType, typename RealType, typename ElementTag>
__device__ inline RealType faceFluxGamma(const FaceFluxElement<KeyType, RealType, ElementTag>& el, const RealType* gamma,
                                         int L, int R)
{
    if (!gamma) return RealType(1);
    if constexpr (std::is_same_v<ElementTag, TetTag>)
        return RealType(0.25) * (gamma[el.n[0]] + gamma[el.n[1]] + gamma[el.n[2]] + gamma[el.n[3]]);
    else
        return RealType(0.5) * (gamma[el.n[L]] + gamma[el.n[R]]);
}

template<typename KeyType, typename RealType, typename ElementTag>
__device__ inline void faceFluxTetGrad(const FaceFluxElement<KeyType, RealType, ElementTag>& el, const RealType* phi,
                                       RealType grad[3])
{
    grad[0] = grad[1] = grad[2] = RealType(0);
    if constexpr (std::is_same_v<ElementTag, TetTag>)
    {
        RealType det, g[4][3];
        Tet4CVFEM::jacobian_and_dNdx<RealType>(el.x, det, g);
        for (int i = 0; i < 4; ++i)
        {
            const RealType p = phi[el.n[i]];
            grad[0] += g[i][0] * p;
            grad[1] += g[i][1] * p;
            grad[2] += g[i][2] * p;
        }
    }
}

// ---------------------------------------------------------------------------------------------
// Kernels. Element-local formulas: run over all local elements (halo included) where the result is
// per element (F), over the owned elements where it is scattered to nodes (then reverse-added).

// G, accumulated: sum_e w_en grad_e phi (and sum_e w_en when wAcc is given), owned elements.
template<typename KeyType, typename RealType, typename ElementTag>
__global__ void faceFluxNodalGradientKernel(const KeyType* c0, const KeyType* c1, const KeyType* c2, const KeyType* c3,
                                            const KeyType* c4, const KeyType* c5, const KeyType* c6, const KeyType* c7,
                                            const RealType* X, const RealType* Y, const RealType* Z, const RealType* phi,
                                            RealType* gx, RealType* gy, RealType* gz, RealType* wAcc, size_t startElem,
                                            size_t numLocal)
{
    const size_t k = blockIdx.x * size_t(blockDim.x) + threadIdx.x;
    if (k >= numLocal) return;
    const KeyType* cc[8] = {c0, c1, c2, c3, c4, c5, c6, c7};
    FaceFluxElement<KeyType, RealType, ElementTag> el;
    el.load(cc, startElem + k, X, Y, Z);
    constexpr int NPE = ElemTraits<ElementTag>::NodesPerElem;
    RealType grad[NPE][3], w[NPE];
    faceFluxElementGradient<KeyType, RealType, ElementTag>(el, phi, grad, w);
    for (int i = 0; i < NPE; ++i)
    {
        atomicAdd(&gx[el.n[i]], w[i] * grad[i][0]);
        atomicAdd(&gy[el.n[i]], w[i] * grad[i][1]);
        atomicAdd(&gz[el.n[i]], w[i] * grad[i][2]);
        if (wAcc) atomicAdd(&wAcc[el.n[i]], w[i]);
    }
}

// g /= w on the owned nodes (the others are set by the forward exchange that follows).
template<typename RealType>
__global__ void faceFluxNormalizeKernel(const uint8_t* own, const RealType* w, RealType* gx, RealType* gy, RealType* gz,
                                        size_t n)
{
    const size_t i = blockIdx.x * size_t(blockDim.x) + threadIdx.x;
    if (i >= n || own[i] != 1) return;
    const RealType iw = w[i] > RealType(0) ? RealType(1) / w[i] : RealType(0);
    gx[i] *= iw;
    gy[i] *= iw;
    gz[i] *= iw;
}

// F* = ubar . A - coefRC (beta_L beta_R)^m [ (grad p)^K_f . A - 1/2 (g_L + g_R) . A ], all local
// elements (weight 1 without the penalised projection); with coefRC = 0, F* = ubar . A.
template<typename KeyType, typename RealType, typename ElementTag>
__global__ void faceFluxPredictKernel(const KeyType* c0, const KeyType* c1, const KeyType* c2, const KeyType* c3,
                                      const KeyType* c4, const KeyType* c5, const KeyType* c6, const KeyType* c7,
                                      const RealType* X, const RealType* Y, const RealType* Z, const RealType* u,
                                      const RealType* v, const RealType* w, const RealType* p, const RealType* Gx,
                                      const RealType* Gy, const RealType* Gz, const RealType* beta, const RealType* Ax,
                                      const RealType* Ay, const RealType* Az, RealType coefRC, int rcBetaPower,
                                      RealType* F, size_t numElements)
{
    const size_t e = blockIdx.x * size_t(blockDim.x) + threadIdx.x;
    if (e >= numElements) return;
    constexpr int NSCS   = ElemTraits<ElementTag>::ScsPerElem;
    const KeyType* cc[8] = {c0, c1, c2, c3, c4, c5, c6, c7};
    FaceFluxElement<KeyType, RealType, ElementTag> el;
    el.load(cc, e, X, Y, Z);
    RealType tg[3];
    faceFluxTetGrad<KeyType, RealType, ElementTag>(el, p, tg);
    for (int ip = 0; ip < NSCS; ++ip)
    {
        int L, R;
        scsLR<ElementTag>(ip, L, R);
        const KeyType nL = el.n[L], nR = el.n[R];
        const size_t f   = e * NSCS + ip;
        const RealType A[3] = {Ax[f], Ay[f], Az[f]};
        RealType flux = RealType(0.5) * ((u[nL] + u[nR]) * A[0] + (v[nL] + v[nR]) * A[1] + (w[nL] + w[nR]) * A[2]);
        if (coefRC != RealType(0))
        {
            const RealType compact = faceFluxKGradDotArea<KeyType, RealType, ElementTag>(el, p, ip, A, tg);
            const RealType smooth  = RealType(0.5) * ((Gx[nL] + Gx[nR]) * A[0] + (Gy[nL] + Gy[nR]) * A[1] +
                                                     (Gz[nL] + Gz[nR]) * A[2]);
            RealType weight = RealType(1);
            if (beta)
            {
                const RealType bb = beta[nL] * beta[nR];
                for (int k = 0; k < rcBetaPower; ++k) weight *= bb;
            }
            flux -= coefRC * weight * (compact - smooth);
        }
        F[f] = flux;
    }
}

// F -= coef beta_f (grad phi)^K_f . A, all local elements.
template<typename KeyType, typename RealType, typename ElementTag>
__global__ void faceFluxCorrectKernel(const KeyType* c0, const KeyType* c1, const KeyType* c2, const KeyType* c3,
                                      const KeyType* c4, const KeyType* c5, const KeyType* c6, const KeyType* c7,
                                      const RealType* X, const RealType* Y, const RealType* Z, const RealType* phi,
                                      const RealType* beta, const RealType* Ax, const RealType* Ay, const RealType* Az,
                                      RealType coef, RealType* F, size_t numElements)
{
    const size_t e = blockIdx.x * size_t(blockDim.x) + threadIdx.x;
    if (e >= numElements) return;
    constexpr int NSCS   = ElemTraits<ElementTag>::ScsPerElem;
    const KeyType* cc[8] = {c0, c1, c2, c3, c4, c5, c6, c7};
    FaceFluxElement<KeyType, RealType, ElementTag> el;
    el.load(cc, e, X, Y, Z);
    RealType tg[3];
    faceFluxTetGrad<KeyType, RealType, ElementTag>(el, phi, tg);
    for (int ip = 0; ip < NSCS; ++ip)
    {
        int L, R;
        scsLR<ElementTag>(ip, L, R);
        const size_t f      = e * NSCS + ip;
        const RealType A[3] = {Ax[f], Ay[f], Az[f]};
        F[f] -= coef * faceFluxGamma<KeyType, RealType, ElementTag>(el, beta, L, R) *
                faceFluxKGradDotArea<KeyType, RealType, ElementTag>(el, phi, ip, A, tg);
    }
}

// acc[L] += F, acc[R] -= F over the owned elements: each control volume's net outflow through its
// interior faces (the projection's divergence; reverse-add afterwards).
template<typename KeyType, typename ElementTag, typename RealType>
__global__ void faceFluxScatterKernel(const KeyType* c0, const KeyType* c1, const KeyType* c2, const KeyType* c3,
                                      const KeyType* c4, const KeyType* c5, const KeyType* c6, const KeyType* c7,
                                      const RealType* F, RealType* acc, size_t startElem, size_t numLocal)
{
    const size_t k = blockIdx.x * size_t(blockDim.x) + threadIdx.x;
    if (k >= numLocal) return;
    constexpr int NPE    = ElemTraits<ElementTag>::NodesPerElem;
    constexpr int NSCS   = ElemTraits<ElementTag>::ScsPerElem;
    const KeyType* cc[8] = {c0, c1, c2, c3, c4, c5, c6, c7};
    const size_t e       = startElem + k;
    KeyType n[NPE];
    for (int i = 0; i < NPE; ++i) n[i] = cc[i][e];
    for (int ip = 0; ip < NSCS; ++ip)
    {
        int L, R;
        scsLR<ElementTag>(ip, L, R);
        const RealType f = F[e * NSCS + ip];
        atomicAdd(&acc[n[L]], f);
        atomicAdd(&acc[n[R]], -f);
    }
}

// K_gamma phi from the face form, over the owned elements (the setup self-check against the
// assembled K): acc[L] -= gamma_f (grad phi)^K_f . A, acc[R] += the same.
template<typename KeyType, typename RealType, typename ElementTag>
__global__ void faceFluxLaplacianKernel(const KeyType* c0, const KeyType* c1, const KeyType* c2, const KeyType* c3,
                                        const KeyType* c4, const KeyType* c5, const KeyType* c6, const KeyType* c7,
                                        const RealType* X, const RealType* Y, const RealType* Z, const RealType* phi,
                                        const RealType* gamma, const RealType* Ax, const RealType* Ay,
                                        const RealType* Az, RealType* acc, size_t startElem, size_t numLocal)
{
    const size_t k = blockIdx.x * size_t(blockDim.x) + threadIdx.x;
    if (k >= numLocal) return;
    constexpr int NSCS   = ElemTraits<ElementTag>::ScsPerElem;
    const KeyType* cc[8] = {c0, c1, c2, c3, c4, c5, c6, c7};
    const size_t e       = startElem + k;
    FaceFluxElement<KeyType, RealType, ElementTag> el;
    el.load(cc, e, X, Y, Z);
    RealType tg[3];
    faceFluxTetGrad<KeyType, RealType, ElementTag>(el, phi, tg);
    for (int ip = 0; ip < NSCS; ++ip)
    {
        int L, R;
        scsLR<ElementTag>(ip, L, R);
        const size_t f      = e * NSCS + ip;
        const RealType A[3] = {Ax[f], Ay[f], Az[f]};
        const RealType q    = faceFluxGamma<KeyType, RealType, ElementTag>(el, gamma, L, R) *
                           faceFluxKGradDotArea<KeyType, RealType, ElementTag>(el, phi, ip, A, tg);
        atomicAdd(&acc[el.n[L]], -q);
        atomicAdd(&acc[el.n[R]], q);
    }
}

// ---------------------------------------------------------------------------------------------
// Pressure outlets with continuity (no Dirichlet rows). At an outlet node i with free velocity the
// control volume's boundary face carries F_b = u_i . A_b (A_b the node's outward opening area), and
// the nodal gradient MARS corrects with, (1/V_i) sum_f 1/2 (phi_j - phi_i) sigma_f A_f, is closed by
// the boundary term (phi_b - phi_i) A_b / V_i with phi_b the prescribed value (Green-Gauss). With
// u^{n+1} = u** - (dtEff/rho) beta_i (closed gradient of phi), continuity in the outlet control volume,
// sum_f sigma_f F_f + F_b = 0, is the row
//   (K phi)_i - (beta_i/V_i) [ sum_f 1/2 (phi_j - phi_i) sigma_f (A_f . A_b) + (phi_b - phi_i) |A_b|^2 ]
//       = -(rho/dtEff) (sum_f sigma_f F*_f + u**_i . A_b)
// (phi_b = 0: the predictor's gradient is closed with the prescribed pressure p_bar itself). The
// added part has row sum beta_i |A_b|^2 / V_i > 0, so the outlets fix the pressure level.

// The face terms of those rows: for each face of an outlet node's control volume, c = (beta_i/V_i)
// 1/2 sigma (A_f . A_b,i) added to (i, i) and subtracted from (i, j). All local elements, owned rows.
template<typename KeyType, typename RealType, typename ElementTag>
__global__ void faceFluxOutletRowsKernel(const KeyType* c0, const KeyType* c1, const KeyType* c2, const KeyType* c3,
                                         const KeyType* c4, const KeyType* c5, const KeyType* c6, const KeyType* c7,
                                         const uint8_t* outlet, const RealType* beta, const RealType* V,
                                         const RealType* Abx, const RealType* Aby, const RealType* Abz,
                                         const RealType* Ax, const RealType* Ay, const RealType* Az, const int* nodeToDof,
                                         const int* rowPtr, const int* colInd, RealType* values, int numOwnedDofs,
                                         size_t numElements)
{
    const size_t e = blockIdx.x * size_t(blockDim.x) + threadIdx.x;
    if (e >= numElements) return;
    constexpr int NSCS   = ElemTraits<ElementTag>::ScsPerElem;
    const KeyType* cc[8] = {c0, c1, c2, c3, c4, c5, c6, c7};
    for (int ip = 0; ip < NSCS; ++ip)
    {
        int L, R;
        scsLR<ElementTag>(ip, L, R);
        const KeyType nL = cc[L][e], nR = cc[R][e];
        const size_t f   = e * NSCS + ip;
        for (int side = 0; side < 2; ++side)
        {
            const KeyType i = side == 0 ? nL : nR, j = side == 0 ? nR : nL;
            if (!outlet[i]) continue;
            const int r = nodeToDof[i];
            if (r < 0 || r >= numOwnedDofs) continue;
            const RealType sigma = side == 0 ? RealType(1) : RealType(-1);
            const RealType b     = beta ? beta[i] : RealType(1);
            const RealType c     = b / V[i] * RealType(0.5) * sigma * (Ax[f] * Abx[i] + Ay[f] * Aby[i] + Az[f] * Abz[i]);
            atomicAddSparseEntry(values, colInd, rowPtr[r], rowPtr[r + 1], r, c);
            atomicAddSparseEntry(values, colInd, rowPtr[r], rowPtr[r + 1], nodeToDof[j], -c);
        }
    }
}

// ... and the boundary term's diagonal, (beta_i/V_i) |A_b,i|^2, one thread per node.
template<typename RealType>
__global__ void faceFluxOutletDiagonalKernel(const uint8_t* outlet, const RealType* beta, const RealType* V,
                                             const RealType* Abx, const RealType* Aby, const RealType* Abz,
                                             const int* nodeToDof, const int* diagPtr, RealType* values,
                                             int numOwnedDofs, size_t n)
{
    const size_t i = blockIdx.x * size_t(blockDim.x) + threadIdx.x;
    if (i >= n || !outlet[i]) return;
    const int r = nodeToDof[i];
    if (r < 0 || r >= numOwnedDofs) return;
    const RealType b = beta ? beta[i] : RealType(1);
    values[diagPtr[r]] += b / V[i] * (Abx[i] * Abx[i] + Aby[i] * Aby[i] + Abz[i] * Abz[i]);
}

// g_i += (target_i - phi_i) A_b,i / V_i at the outlet nodes: the boundary closure of MARS's nodal
// gradient (target = p_bar for the predictor's grad p^n, 0 for the corrector's grad phi).
template<typename RealType>
__global__ void faceFluxOutletGradientKernel(const uint8_t* outlet, const RealType* phi, const RealType* targetDof,
                                             const int* nodeToDof, const RealType* V, const RealType* Abx,
                                             const RealType* Aby, const RealType* Abz, RealType* gx, RealType* gy,
                                             RealType* gz, int numOwnedDofs, size_t n)
{
    const size_t i = blockIdx.x * size_t(blockDim.x) + threadIdx.x;
    if (i >= n || !outlet[i]) return;
    const int r = nodeToDof[i];
    if (r < 0 || r >= numOwnedDofs) return;
    const RealType t = targetDof ? targetDof[r] : RealType(0);
    const RealType s = (t - phi[i]) / V[i];
    gx[i] += s * Abx[i];
    gy[i] += s * Aby[i];
    gz[i] += s * Abz[i];
}

// acc_i += u_i . A_b,i at the outlet nodes (the boundary face's flux).
template<typename RealType>
__global__ void faceFluxOutletFluxKernel(const uint8_t* outlet, const RealType* u, const RealType* v, const RealType* w,
                                         const RealType* Abx, const RealType* Aby, const RealType* Abz, RealType* acc,
                                         size_t n)
{
    const size_t i = blockIdx.x * size_t(blockDim.x) + threadIdx.x;
    if (i >= n || !outlet[i]) return;
    acc[i] += u[i] * Abx[i] + v[i] * Aby[i] + w[i] * Abz[i];
}

// ---------------------------------------------------------------------------------------------
// Convection on the stored flux (the SFEM spike's interior flux, cvfem_hex8_scs_convection and
// cvfem_hex8_scs_defcor, on MARS's geometry). Per face f from L to R with flux F (from the last
// projection, so sum F = 0 per control volume): F+ = (F + |F|_eps)/2, F- = (F - |F|_eps)/2, and the
// momentum flux F+ u_L + F- u_R, implicit (first order, the matrix) plus the deferred correction
// F+ D_L + F- D_R, explicit (the right-hand side), D_C the limited increment from node C to the face
// (the edge midpoint, as MARS's own CVFEM advection reconstructs) along G u_C. At an outlet node's
// boundary face max(F_b, 0) u_i: outflow carries the node's momentum out, inflow brings none (the
// spike's backflow guard; the energy flux through the face is then outgoing only).

// The first-order matrix: all local elements, owned non-Dirichlet rows (as the Galerkin convection
// kernels), Dirichlet columns kept.
template<typename KeyType, typename RealType, typename ElementTag>
__global__ void faceFluxConvectionKernel(const KeyType* c0, const KeyType* c1, const KeyType* c2, const KeyType* c3,
                                         const KeyType* c4, const KeyType* c5, const KeyType* c6, const KeyType* c7,
                                         const RealType* F, RealType eps, const int* nodeToDof, const uint8_t* isBdryDof,
                                         const int* rowPtr, const int* colInd, RealType* values, int numOwnedDofs,
                                         size_t numElements)
{
    const size_t e = blockIdx.x * size_t(blockDim.x) + threadIdx.x;
    if (e >= numElements) return;
    constexpr int NSCS   = ElemTraits<ElementTag>::ScsPerElem;
    const KeyType* cc[8] = {c0, c1, c2, c3, c4, c5, c6, c7};
    for (int ip = 0; ip < NSCS; ++ip)
    {
        int L, R;
        scsLR<ElementTag>(ip, L, R);
        const RealType f = F[e * NSCS + ip];
        if (f == RealType(0)) continue;
        const RealType a  = faceFluxUpwindAbs(f, eps);
        const RealType fp = RealType(0.5) * (f + a), fm = RealType(0.5) * (f - a);
        const int dL = nodeToDof[cc[L][e]], dR = nodeToDof[cc[R][e]];
        if (dL >= 0 && dL < numOwnedDofs && !isBdryDof[dL])
        {
            atomicAddSparseEntry(values, colInd, rowPtr[dL], rowPtr[dL + 1], dL, fp);
            atomicAddSparseEntry(values, colInd, rowPtr[dL], rowPtr[dL + 1], dR, fm);
        }
        if (dR >= 0 && dR < numOwnedDofs && !isBdryDof[dR])
        {
            atomicAddSparseEntry(values, colInd, rowPtr[dR], rowPtr[dR + 1], dL, -fp);
            atomicAddSparseEntry(values, colInd, rowPtr[dR], rowPtr[dR + 1], dR, -fm);
        }
    }
}

// The outlets' boundary faces: max(u_i . A_b, 0) on the diagonal (non-Dirichlet rows).
template<typename RealType>
__global__ void faceFluxOutletConvectionKernel(const uint8_t* outlet, const RealType* u, const RealType* v,
                                               const RealType* w, const RealType* Abx, const RealType* Aby,
                                               const RealType* Abz, const int* nodeToDof, const uint8_t* isBdryDof,
                                               const int* diagPtr, RealType* values, int numOwnedDofs, size_t n)
{
    const size_t i = blockIdx.x * size_t(blockDim.x) + threadIdx.x;
    if (i >= n || !outlet[i]) return;
    const int r = nodeToDof[i];
    if (r < 0 || r >= numOwnedDofs || isBdryDof[r]) return;
    const RealType fb = u[i] * Abx[i] + v[i] * Aby[i] + w[i] * Abz[i];
    if (fb > RealType(0)) values[diagPtr[r]] += fb;
}

// The deferred correction of one velocity component q, as each control volume's net outflow:
// acc[L] += F+ D_L + F- D_R, acc[R] -= the same, owned elements (reverse-add afterwards).
// limiter: 0 none, 1 bounded-face clip, 2 Venkatakrishnan (eps2 = venkatC |x_R - x_L|^3),
// 3 Darwish-Moukalled.
template<typename KeyType, typename RealType, typename ElementTag>
__global__ void faceFluxDeferredCorrectionKernel(const KeyType* c0, const KeyType* c1, const KeyType* c2,
                                                 const KeyType* c3, const KeyType* c4, const KeyType* c5,
                                                 const KeyType* c6, const KeyType* c7, const RealType* X,
                                                 const RealType* Y, const RealType* Z, const RealType* F, RealType eps,
                                                 const RealType* q, const RealType* gx, const RealType* gy,
                                                 const RealType* gz, int limiter, RealType venkatC, RealType* acc,
                                                 size_t startElem, size_t numLocal)
{
    const size_t k = blockIdx.x * size_t(blockDim.x) + threadIdx.x;
    if (k >= numLocal) return;
    constexpr int NSCS   = ElemTraits<ElementTag>::ScsPerElem;
    const KeyType* cc[8] = {c0, c1, c2, c3, c4, c5, c6, c7};
    const size_t e       = startElem + k;
    for (int ip = 0; ip < NSCS; ++ip)
    {
        int L, R;
        scsLR<ElementTag>(ip, L, R);
        const RealType f = F[e * NSCS + ip];
        if (f == RealType(0)) continue;
        const KeyType nL = cc[L][e], nR = cc[R][e];
        const RealType dx = X[nR] - X[nL], dy = Y[nR] - Y[nL], dz = Z[nR] - Z[nL];
        // node -> edge midpoint: +d/2 from L, -d/2 from R
        RealType incL = RealType(0.5) * (gx[nL] * dx + gy[nL] * dy + gz[nL] * dz);
        RealType incR = -RealType(0.5) * (gx[nR] * dx + gy[nR] * dy + gz[nR] * dz);
        const RealType qL = q[nL], qR = q[nR];
        const RealType lo = qL < qR ? qL : qR, hi = qL < qR ? qR : qL;
        if (limiter == 1)
        {
            incL = faceFluxClipInc(qL, incL, lo, hi);
            incR = faceFluxClipInc(qR, incR, lo, hi);
        }
        else if (limiter == 2)
        {
            const RealType h2 = dx * dx + dy * dy + dz * dz, eps2 = venkatC * h2 * sqrt(h2);
            incL              = faceFluxVenkatakrishnanInc(qL, incL, lo, hi, eps2);
            incR              = faceFluxVenkatakrishnanInc(qR, incR, lo, hi, eps2);
        }
        else if (limiter == 3)
        {
            incL = faceFluxDarwishMoukalledInc(qL, qR, incL);
            incR = faceFluxDarwishMoukalledInc(qR, qL, incR);
        }
        const RealType a    = faceFluxUpwindAbs(f, eps);
        const RealType corr = RealType(0.5) * (f + a) * incL + RealType(0.5) * (f - a) * incR;
        atomicAdd(&acc[nL], corr);
        atomicAdd(&acc[nR], -corr);
    }
}

// y = A x on the owned rows of a CSR matrix over DOFs (the self-check's reference product).
template<typename RealType>
__global__ void faceFluxCsrApplyKernel(const int* rowPtr, const int* colInd, const RealType* vals, const RealType* x,
                                       RealType* y, int nRows)
{
    const int r = blockIdx.x * blockDim.x + threadIdx.x;
    if (r >= nRows) return;
    RealType s = 0;
    for (int k = rowPtr[r]; k < rowPtr[r + 1]; ++k) s += vals[k] * x[colInd[k]];
    y[r] = s;
}

} // namespace fem
} // namespace mars
