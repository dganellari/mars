#pragma once

// Incompressible Navier-Stokes on a triply periodic hex mesh (CVFEM, Q1).
//
// Every field -- the three velocity components and the pressure -- lives in
// the reduced periodic space of mars_periodic_space.hpp: one DOF per periodic
// point. The state vectors hold one value per node slot and are always kept
// prolonged (every slot of a periodic point holds the DOF value). Every
// operator is restrict(element scatter(prolonged input)) = P^T A P:
//
//   M    lumped mass                     (diagonal, per DOF)
//   D    divergence    (D u)_i = sum_f +-A_f . (u_L + u_R) / 2
//   G    gradient      G = D^T, so M^-1 G p ~ -grad p
//   K    CVFEM viscous stiffness         (symmetric on the periodic box)
//   N(u) explicit advection, skew-symmetric or upwind
//
// One time step, BDF2 + EXT2 incremental pressure correction (BDF1 on the
// first step), with dtEff = dt (BDF1) or 2 dt / 3 (BDF2):
//
//   predictor   u*  = BDF extrapolation + dtEff (M^-1 N_ext + M^-1 G p^n / rho)
//   viscous     (M / dtEff + nu K) u** = (M / dtEff) u*
//   projection  A phi = -(rho / dtEff) D u**      with A = D M^-1 G
//   corrector   u^{n+1} = u** + (dtEff / rho) M^-1 G phi,   p^{n+1} = p^n + phi
//
// Because D, G and M all come from the same P, D u^{n+1} = D u** + (dtEff/rho) A phi
// = 0 holds exactly on every rank count, and the corrected velocity is
// single-valued at every periodic point without any extra copy.

#include "backend/distributed/unstructured/domain.hpp"
#include "backend/distributed/unstructured/fem/mars_cvfem_matfree.hpp"
#include "backend/distributed/unstructured/fem/mars_cvfem_utils.hpp"
#include "backend/distributed/unstructured/fem/mars_periodic_space.hpp"

#include <thrust/execution_policy.h>
#include <thrust/functional.h>
#include <thrust/iterator/counting_iterator.h>
#include <thrust/transform_reduce.h>

#include <algorithm>
#include <cmath>
#include <mpi.h>

namespace mars
{
namespace fem
{

// Owned elements of a hex domain: kernels loop over k in [0, count), read
// connectivity at element first + k and per-element geometry at k.
template<typename KeyType>
struct OwnedHexElements
{
    const KeyType* node[8];
    size_t first;
    size_t count;
};

// ---------------------------------------------------------------------------
// Element scatters over owned elements. They write every corner slot (owned or
// ghost); PeriodicSpace::restrict then completes and reduces the sums.
// ---------------------------------------------------------------------------

// Each corner takes 1/8 of the element volume. The bounding-box volume is
// exact for the axis-aligned hexes of a periodic box.
template<typename KeyType, typename RealType>
__global__ void pnsLumpedMassKernel(OwnedHexElements<KeyType> hex, const RealType* x, const RealType* y,
                                    const RealType* z, RealType* mass)
{
    size_t k = size_t(blockIdx.x) * blockDim.x + threadIdx.x;
    if (k >= hex.count) return;
    size_t e = hex.first + k;
    KeyType n0 = hex.node[0][e];
    RealType lo[3] = {x[n0], y[n0], z[n0]};
    RealType hi[3] = {lo[0], lo[1], lo[2]};
    for (int c = 1; c < 8; ++c)
    {
        KeyType n    = hex.node[c][e];
        RealType p[3] = {x[n], y[n], z[n]};
        for (int d = 0; d < 3; ++d)
        {
            lo[d] = fmin(lo[d], p[d]);
            hi[d] = fmax(hi[d], p[d]);
        }
    }
    RealType share = (hi[0] - lo[0]) * (hi[1] - lo[1]) * (hi[2] - lo[2]) * RealType(0.125);
    for (int c = 0; c < 8; ++c)
        atomicAdd(&mass[hex.node[c][e]], share);
}

// D u: the flux through each sub-control face leaves node L and enters node R.
template<typename KeyType, typename RealType>
__global__ void pnsDivergenceKernel(OwnedHexElements<KeyType> hex, const RealType* ax, const RealType* ay,
                                    const RealType* az, const RealType* u, const RealType* v, const RealType* w,
                                    RealType* div)
{
    size_t k = size_t(blockIdx.x) * blockDim.x + threadIdx.x;
    if (k >= hex.count) return;
    size_t e = hex.first + k;
    KeyType n[8];
    for (int c = 0; c < 8; ++c) n[c] = hex.node[c][e];
#pragma unroll
    for (int ip = 0; ip < 12; ++ip)
    {
        KeyType L = n[hexLRSCV[2 * ip]], R = n[hexLRSCV[2 * ip + 1]];
        size_t f      = k * 12 + ip;
        RealType flow = RealType(0.5) * ((u[L] + u[R]) * ax[f] + (v[L] + v[R]) * ay[f] + (w[L] + w[R]) * az[f]);
        atomicAdd(&div[L], flow);
        atomicAdd(&div[R], -flow);
    }
}

// G p = D^T p: the transpose of the scatter above, so it adds the same
// (p_L - p_R) / 2 * A_f to both ends of the face.
template<typename KeyType, typename RealType>
__global__ void pnsDivergenceTransposeKernel(OwnedHexElements<KeyType> hex, const RealType* ax, const RealType* ay,
                                             const RealType* az, const RealType* p, RealType* gx, RealType* gy,
                                             RealType* gz)
{
    size_t k = size_t(blockIdx.x) * blockDim.x + threadIdx.x;
    if (k >= hex.count) return;
    size_t e = hex.first + k;
    KeyType n[8];
    for (int c = 0; c < 8; ++c) n[c] = hex.node[c][e];
#pragma unroll
    for (int ip = 0; ip < 12; ++ip)
    {
        KeyType L = n[hexLRSCV[2 * ip]], R = n[hexLRSCV[2 * ip + 1]];
        size_t f    = k * 12 + ip;
        RealType dp = RealType(0.5) * (p[L] - p[R]);
        atomicAdd(&gx[L], dp * ax[f]);
        atomicAdd(&gy[L], dp * ay[f]);
        atomicAdd(&gz[L], dp * az[f]);
        atomicAdd(&gx[R], dp * ax[f]);
        atomicAdd(&gy[R], dp * ay[f]);
        atomicAdd(&gz[R], dp * az[f]);
    }
}

// Advection of the three velocity components by the face mass flux mdot.
// Skew: the Verstappen skew-symmetric flux. Its kinetic-energy production is
// sum_i q_i^2 (D u)_i, zero once the projection makes D u = 0, so advection
// neither creates nor destroys energy. Upwind: first order, dissipative.
template<typename KeyType, typename RealType>
__global__ void pnsAdvectionKernel(OwnedHexElements<KeyType> hex, const RealType* ax, const RealType* ay,
                                   const RealType* az, const RealType* u, const RealType* v, const RealType* w,
                                   int skew, RealType* au, RealType* av, RealType* aw)
{
    size_t k = size_t(blockIdx.x) * blockDim.x + threadIdx.x;
    if (k >= hex.count) return;
    size_t e = hex.first + k;
    KeyType n[8];
    for (int c = 0; c < 8; ++c) n[c] = hex.node[c][e];
    const RealType* q[3] = {u, v, w};
    RealType* out[3]     = {au, av, aw};
#pragma unroll
    for (int ip = 0; ip < 12; ++ip)
    {
        KeyType L = n[hexLRSCV[2 * ip]], R = n[hexLRSCV[2 * ip + 1]];
        size_t f      = k * 12 + ip;
        RealType mdot = RealType(0.5) * ((u[L] + u[R]) * ax[f] + (v[L] + v[R]) * ay[f] + (w[L] + w[R]) * az[f]);
        for (int c = 0; c < 3; ++c)
        {
            RealType qL = q[c][L], qR = q[c][R];
            if (skew)
            {
                atomicAdd(&out[c][L], -RealType(0.5) * mdot * (RealType(2) * qL + qR));
                atomicAdd(&out[c][R], RealType(0.5) * mdot * (qL + RealType(2) * qR));
            }
            else
            {
                RealType flux = mdot * (mdot > RealType(0) ? qL : qR);
                atomicAdd(&out[c][L], -flux);
                atomicAdd(&out[c][R], flux);
            }
        }
    }
}

// K q with the element matrices stored once at setup (8x8 per owned element).
template<typename KeyType, typename RealType>
__global__ void pnsStiffnessKernel(OwnedHexElements<KeyType> hex, const RealType* Ke, const RealType* q,
                                   RealType* out)
{
    size_t k = size_t(blockIdx.x) * blockDim.x + threadIdx.x;
    if (k >= hex.count) return;
    size_t e = hex.first + k;
    KeyType n[8];
    RealType qe[8];
    for (int c = 0; c < 8; ++c)
    {
        n[c]  = hex.node[c][e];
        qe[c] = q[n[c]];
    }
    const RealType* K = Ke + k * 64;
    for (int i = 0; i < 8; ++i)
    {
        RealType y = 0;
#pragma unroll
        for (int j = 0; j < 8; ++j)
            y += K[i * 8 + j] * qe[j];
        atomicAdd(&out[n[i]], y);
    }
}

template<typename KeyType, typename RealType>
__global__ void pnsStiffnessDiagonalKernel(OwnedHexElements<KeyType> hex, const RealType* Ke, RealType* diag)
{
    size_t k = size_t(blockIdx.x) * blockDim.x + threadIdx.x;
    if (k >= hex.count) return;
    size_t e = hex.first + k;
    for (int i = 0; i < 8; ++i)
        atomicAdd(&diag[hex.node[i][e]], Ke[k * 64 + i * 9]);
}

// Jacobi diagonal of A = D M^-1 G keeping each face's own coupling only:
// 0.25 |A_f|^2 (1/M_L + 1/M_R) on both ends. Positive, and constant on a
// uniform mesh, so it only helps on graded (AMR) meshes. mass is prolonged.
template<typename KeyType, typename RealType>
__global__ void pnsProjectionDiagonalKernel(OwnedHexElements<KeyType> hex, const RealType* ax, const RealType* ay,
                                            const RealType* az, const RealType* mass, RealType* diag)
{
    size_t k = size_t(blockIdx.x) * blockDim.x + threadIdx.x;
    if (k >= hex.count) return;
    size_t e = hex.first + k;
    KeyType n[8];
    for (int c = 0; c < 8; ++c) n[c] = hex.node[c][e];
#pragma unroll
    for (int ip = 0; ip < 12; ++ip)
    {
        KeyType L = n[hexLRSCV[2 * ip]], R = n[hexLRSCV[2 * ip + 1]];
        size_t f     = k * 12 + ip;
        RealType a2  = ax[f] * ax[f] + ay[f] * ay[f] + az[f] * az[f];
        RealType val = RealType(0.25) * a2 * (RealType(1) / mass[L] + RealType(1) / mass[R]);
        atomicAdd(&diag[L], val);
        atomicAdd(&diag[R], val);
    }
}

// ---------------------------------------------------------------------------
// Per-DOF kernels. Slots that are not DOFs are set to 0; prolong fills them.
// ---------------------------------------------------------------------------

// y = a + alpha * b
template<typename RealType>
__global__ void pnsAxpyKernel(const uint8_t* isDof, size_t n, const RealType* a, RealType alpha, const RealType* b,
                              RealType* y)
{
    size_t i = size_t(blockIdx.x) * blockDim.x + threadIdx.x;
    if (i >= n) return;
    y[i] = isDof[i] ? a[i] + alpha * b[i] : RealType(0);
}

template<typename RealType>
__global__ void pnsJacobiKernel(const uint8_t* isDof, size_t n, const RealType* r, const RealType* diag, RealType* z)
{
    size_t i = size_t(blockIdx.x) * blockDim.x + threadIdx.x;
    if (i >= n) return;
    z[i] = isDof[i] ? r[i] / diag[i] : RealType(0);
}

// out = c * M * q
template<typename RealType>
__global__ void pnsMassTimesKernel(const uint8_t* isDof, size_t n, RealType c, const RealType* mass,
                                   const RealType* q, RealType* out)
{
    size_t i = size_t(blockIdx.x) * blockDim.x + threadIdx.x;
    if (i >= n) return;
    out[i] = isDof[i] ? c * mass[i] * q[i] : RealType(0);
}

// On entry Kq holds the restricted K q; on exit (c M + nu K) q.
template<typename RealType>
__global__ void pnsViscousCombineKernel(const uint8_t* isDof, size_t n, RealType c, RealType nu,
                                        const RealType* mass, const RealType* q, RealType* Kq)
{
    size_t i = size_t(blockIdx.x) * blockDim.x + threadIdx.x;
    if (i >= n) return;
    Kq[i] = isDof[i] ? nu * Kq[i] + c * mass[i] * q[i] : RealType(0);
}

template<typename RealType>
__global__ void pnsViscousDiagonalKernel(const uint8_t* isDof, size_t n, RealType c, RealType nu,
                                         const RealType* mass, const RealType* diagK, RealType* diag)
{
    size_t i = size_t(blockIdx.x) * blockDim.x + threadIdx.x;
    if (i >= n) return;
    diag[i] = isDof[i] ? c * mass[i] + nu * diagK[i] : RealType(1);
}

template<typename RealType>
__global__ void pnsInverseMassKernel(const uint8_t* isDof, size_t n, const RealType* mass, RealType* gx,
                                     RealType* gy, RealType* gz)
{
    size_t i = size_t(blockIdx.x) * blockDim.x + threadIdx.x;
    if (i >= n) return;
    RealType s = isDof[i] ? RealType(1) / mass[i] : RealType(0);
    gx[i] *= s;
    gy[i] *= s;
    gz[i] *= s;
}

// u* from u^n (and u^{n-1} for BDF2), advection a = N(u^n) and a_prev =
// N(u^{n-1}) extrapolated to 2 a - a_prev, and g = M^-1 G p^n = -grad p^n.
// Writes u^n into the history slot for the next step.
template<typename RealType>
__global__ void pnsPredictorKernel(const uint8_t* isDof, size_t n, int bdf2, RealType dt, RealType invRho,
                                   const RealType* mass, const RealType* u, const RealType* v, const RealType* w,
                                   RealType* um1, RealType* vm1, RealType* wm1, const RealType* au,
                                   const RealType* av, const RealType* aw, const RealType* aum1,
                                   const RealType* avm1, const RealType* awm1, const RealType* gx,
                                   const RealType* gy, const RealType* gz, RealType* us, RealType* vs, RealType* ws)
{
    size_t i = size_t(blockIdx.x) * blockDim.x + threadIdx.x;
    if (i >= n) return;
    if (!isDof[i])
    {
        us[i] = vs[i] = ws[i] = RealType(0);
        return;
    }
    RealType invM            = RealType(1) / mass[i];
    const RealType q[3]      = {u[i], v[i], w[i]};
    RealType* qm1[3]         = {um1, vm1, wm1};
    const RealType a[3]      = {au[i], av[i], aw[i]};
    const RealType* am1[3]   = {aum1, avm1, awm1};
    const RealType g[3]      = {gx[i], gy[i], gz[i]};
    RealType* out[3]         = {us, vs, ws};
    for (int c = 0; c < 3; ++c)
    {
        if (bdf2)
        {
            RealType ext = RealType(2) * a[c] - am1[c][i];
            out[c][i]    = (RealType(4) * q[c] - qm1[c][i]) / RealType(3) +
                        RealType(2) * dt / RealType(3) * (ext * invM + g[c] * invRho);
        }
        else
        {
            out[c][i] = q[c] + dt * (a[c] * invM + g[c] * invRho);
        }
        qm1[c][i] = q[c];
    }
}

// u = u** + s * g
template<typename RealType>
__global__ void pnsCorrectorKernel(const uint8_t* isDof, size_t n, RealType s, const RealType* gx,
                                   const RealType* gy, const RealType* gz, const RealType* us, const RealType* vs,
                                   const RealType* ws, RealType* u, RealType* v, RealType* w)
{
    size_t i = size_t(blockIdx.x) * blockDim.x + threadIdx.x;
    if (i >= n) return;
    bool dof = isDof[i];
    u[i]     = dof ? us[i] + s * gx[i] : RealType(0);
    v[i]     = dof ? vs[i] + s * gy[i] : RealType(0);
    w[i]     = dof ? ws[i] + s * gz[i] : RealType(0);
}

// Per-DOF scale of a velocity a divergence could have: |u| times a face area M^(2/3).
template<typename RealType>
__global__ void pnsFluxScaleKernel(const uint8_t* isDof, size_t n, const RealType* mass, const RealType* u,
                                   const RealType* v, const RealType* w, RealType* out)
{
    size_t i = size_t(blockIdx.x) * blockDim.x + threadIdx.x;
    if (i >= n) return;
    out[i] = isDof[i] ? sqrt(u[i] * u[i] + v[i] * v[i] + w[i] * w[i]) * cbrt(mass[i] * mass[i]) : RealType(0);
}

template<typename RealType>
__global__ void pnsKineticEnergyKernel(const uint8_t* isDof, size_t n, const RealType* mass, const RealType* u,
                                       const RealType* v, const RealType* w, RealType* out)
{
    size_t i = size_t(blockIdx.x) * blockDim.x + threadIdx.x;
    if (i >= n) return;
    out[i] = isDof[i] ? RealType(0.5) * mass[i] * (u[i] * u[i] + v[i] * v[i] + w[i] * w[i]) : RealType(0);
}

// Adds grad(q_c) to omega = curl(u). g = M^-1 G q_c = -grad q_c.
template<typename RealType>
__global__ void pnsCurlKernel(const uint8_t* isDof, size_t n, int c, const RealType* gx, const RealType* gy,
                              const RealType* gz, RealType* ox, RealType* oy, RealType* oz)
{
    size_t i = size_t(blockIdx.x) * blockDim.x + threadIdx.x;
    if (i >= n || !isDof[i]) return;
    if (c == 0) { oy[i] -= gz[i]; oz[i] += gy[i]; }       // + du/dz, - du/dy
    if (c == 1) { oz[i] -= gx[i]; ox[i] += gz[i]; }       // + dv/dx, - dv/dz
    if (c == 2) { ox[i] -= gy[i]; oy[i] += gx[i]; }       // + dw/dy, - dw/dx
}

template<typename RealType>
__global__ void pnsMagnitudeKernel(const uint8_t* isDof, size_t n, const RealType* x, const RealType* y,
                                   const RealType* z, RealType* out)
{
    size_t i = size_t(blockIdx.x) * blockDim.x + threadIdx.x;
    if (i >= n) return;
    out[i] = isDof[i] ? sqrt(x[i] * x[i] + y[i] * y[i] + z[i] * z[i]) : RealType(0);
}

template<typename KeyType, typename RealType>
class PeriodicNavierStokes
{
public:
    using Domain = ElementDomain<HexTag, RealType, KeyType, cstone::execution::Gpu>;
    using Space  = PeriodicSpace<KeyType, RealType, Domain>;
    using Vector = cstone::DeviceVector<RealType>;

    struct Params
    {
        RealType nu        = RealType(1) / RealType(1600);
        RealType rho       = 1;
        RealType dt        = RealType(1e-3);
        bool skewAdvection = true;
        int maxIter        = 1000;
        RealType tolerance = RealType(1e-10);
        int blockSize      = 256;
    };

    // State, one value per node slot, kept prolonged between calls.
    Vector u, v, w, p;

    PeriodicNavierStokes(Domain& domain, const PeriodicMap<KeyType, RealType>& map, const Params& params)
        : domain_(domain)
        , space_(domain, map, params.blockSize)
        , prm_(params)
        , n_(domain.getNodeCount())
    {
        const auto& conn = domain.getElementToNodeConnectivity();
        hex_ = {{std::get<0>(conn).data(), std::get<1>(conn).data(), std::get<2>(conn).data(),
                 std::get<3>(conn).data(), std::get<4>(conn).data(), std::get<5>(conn).data(),
                 std::get<6>(conn).data(), std::get<7>(conn).data()},
                domain.startIndex(),
                domain.localElementCount()};

        for (Vector* f : {&u, &v, &w, &p, &um1_, &vm1_, &wm1_, &us_, &vs_, &ws_, &au_, &av_, &aw_, &aum1_, &avm1_,
                          &awm1_, &gx_, &gy_, &gz_, &phi_, &rhs_, &mass_, &diagK_, &diagP_, &diagV_, &r_, &z_,
                          &d_, &Ad_, &scratch_})
        {
            f->resize(n_);
            zero(*f);
        }
        buildGeometry();
        buildMassAndDiagonals();
    }

    // Call after writing u, v, w, p: make them single-valued at periodic
    // points and fix the pressure gauge.
    void start()
    {
        space_.prolong(u);
        space_.prolong(v);
        space_.prolong(w);
        space_.removeMean(p);
        space_.prolong(p);
        numSteps_ = 0;
    }

    // One time step. Returns false (on every rank) if a linear solve failed.
    bool step()
    {
        const bool bdf2        = numSteps_ > 0;
        const RealType dt      = prm_.dt;
        const RealType dtEff   = bdf2 ? RealType(2) * dt / RealType(3) : dt;
        const RealType invRho  = RealType(1) / prm_.rho;
        const uint8_t* dof     = space_.isDof();

        // 1. Predictor: explicit advection and the old pressure gradient.
        gradientOf(p);
        advectionOf();
        pnsPredictorKernel<RealType><<<grid(), bs()>>>(dof, n_, int(bdf2), dt, invRho, mass_.data(), u.data(),
                                                       v.data(), w.data(), um1_.data(), vm1_.data(), wm1_.data(),
                                                       au_.data(), av_.data(), aw_.data(), aum1_.data(),
                                                       avm1_.data(), awm1_.data(), gx_.data(), gy_.data(),
                                                       gz_.data(), us_.data(), vs_.data(), ws_.data());
        cudaCheckError();
        au_.swap(aum1_);
        av_.swap(avm1_);
        aw_.swap(awm1_);

        // 2. Implicit viscous step, one CG per component, u** overwrites u*.
        const RealType c = RealType(1) / dtEff;
        pnsViscousDiagonalKernel<RealType><<<grid(), bs()>>>(dof, n_, c, prm_.nu, mass_.data(), diagK_.data(),
                                                             diagV_.data());
        cudaCheckError();
        auto viscous = [this, c](Vector& q, Vector& out) { applyViscous(c, q, out); };
        Vector* star[3] = {&us_, &vs_, &ws_};
        for (int k = 0; k < 3; ++k)
        {
            pnsMassTimesKernel<RealType><<<grid(), bs()>>>(dof, n_, c, mass_.data(), star[k]->data(), rhs_.data());
            cudaCheckError();
            velocityIters_[k] = pcg(viscous, rhs_, *star[k], diagV_, RealType(0), false);
            if (velocityIters_[k] < 0) return false;
        }

        // 3. Projection: A phi = -(rho / dtEff) D u**.
        divergenceOf(us_, vs_, ws_, rhs_);
        scale(rhs_, -prm_.rho / dtEff);
        space_.removeMean(rhs_);
        pnsFluxScaleKernel<RealType><<<grid(), bs()>>>(dof, n_, mass_.data(), us_.data(), vs_.data(), ws_.data(),
                                                       scratch_.data());
        cudaCheckError();
        // The increment D u** shrinks as the pressure converges in time; a CG
        // stopping rule relative to it alone would ask for a residual below
        // roundoff. Also accept tol times the divergence this velocity could carry.
        const RealType bRef = prm_.rho / dtEff * std::sqrt(space_.dot(scratch_, scratch_));
        zero(phi_);
        auto projection = [this](Vector& x, Vector& out) { applyProjection(x, out); };
        pressureIters_ = pcg(projection, rhs_, phi_, diagP_, bRef, true);
        if (pressureIters_ < 0) return false;
        space_.removeMean(phi_);

        // 4. Corrector and pressure update.
        gradientOf(phi_);
        pnsCorrectorKernel<RealType><<<grid(), bs()>>>(dof, n_, dtEff * invRho, gx_.data(), gy_.data(), gz_.data(),
                                                       us_.data(), vs_.data(), ws_.data(), u.data(), v.data(),
                                                       w.data());
        pnsAxpyKernel<RealType><<<grid(), bs()>>>(dof, n_, p.data(), RealType(1), phi_.data(), p.data());
        cudaCheckError();
        space_.prolong(u);
        space_.prolong(v);
        space_.prolong(w);
        space_.removeMean(p);
        space_.prolong(p);

        ++numSteps_;
        return true;
    }

    RealType kineticEnergy()
    {
        pnsKineticEnergyKernel<RealType><<<grid(), bs()>>>(space_.isDof(), n_, mass_.data(), u.data(), v.data(),
                                                           w.data(), scratch_.data());
        cudaCheckError();
        return space_.sum(scratch_);
    }

    // max over DOFs of |D u| / M: the pointwise divergence the projection removed.
    RealType maxDivergence()
    {
        divergenceOf(u, v, w, scratch_);
        const uint8_t* dof   = space_.isDof();
        const RealType* div  = scratch_.data();
        const RealType* mass = mass_.data();
        RealType local       = thrust::transform_reduce(
            thrust::device, thrust::counting_iterator<size_t>(0), thrust::counting_iterator<size_t>(n_),
            [dof, div, mass] __device__(size_t i) { return dof[i] ? fabs(div[i] / mass[i]) : RealType(0); },
            RealType(0), thrust::maximum<RealType>());
        RealType global = 0;
        MPI_Allreduce(&local, &global, 1, mpiDatatype<RealType>(), MPI_MAX, MPI_COMM_WORLD);
        return global;
    }

    // |curl u| per node slot, for output.
    void vorticityMagnitude(Vector& out)
    {
        Vector* omega[3] = {&omegaX_, &omegaY_, &omegaZ_};
        for (Vector* o : omega)
        {
            o->resize(n_);
            zero(*o);
        }
        const Vector* comp[3] = {&u, &v, &w};
        for (int c = 0; c < 3; ++c)
        {
            gradientOf(*comp[c]);
            pnsCurlKernel<RealType><<<grid(), bs()>>>(space_.isDof(), n_, c, gx_.data(), gy_.data(), gz_.data(),
                                                      omegaX_.data(), omegaY_.data(), omegaZ_.data());
            cudaCheckError();
        }
        out.resize(n_);
        pnsMagnitudeKernel<RealType><<<grid(), bs()>>>(space_.isDof(), n_, omegaX_.data(), omegaY_.data(),
                                                       omegaZ_.data(), out.data());
        cudaCheckError();
        space_.prolong(out);
    }

    const Space& space() const { return space_; }
    int velocityIterations(int component) const { return velocityIters_[component]; }
    int pressureIterations() const { return pressureIters_; }

private:
    // gx, gy, gz = M^-1 G q = -grad q at DOF slots. q must be prolonged.
    void gradientOf(const Vector& q)
    {
        zero(gx_);
        zero(gy_);
        zero(gz_);
        pnsDivergenceTransposeKernel<KeyType, RealType><<<elemGrid(), bs()>>>(
            hex_, areaX_.data(), areaY_.data(), areaZ_.data(), q.data(), gx_.data(), gy_.data(), gz_.data());
        cudaCheckError();
        space_.restrict(gx_);
        space_.restrict(gy_);
        space_.restrict(gz_);
        pnsInverseMassKernel<RealType><<<grid(), bs()>>>(space_.isDof(), n_, mass_.data(), gx_.data(), gy_.data(),
                                                         gz_.data());
        cudaCheckError();
    }

    // out = D (a, b, c) per DOF. The inputs must be prolonged.
    void divergenceOf(const Vector& a, const Vector& b, const Vector& c, Vector& out)
    {
        zero(out);
        pnsDivergenceKernel<KeyType, RealType><<<elemGrid(), bs()>>>(
            hex_, areaX_.data(), areaY_.data(), areaZ_.data(), a.data(), b.data(), c.data(), out.data());
        cudaCheckError();
        space_.restrict(out);
    }

    void advectionOf()
    {
        zero(au_);
        zero(av_);
        zero(aw_);
        pnsAdvectionKernel<KeyType, RealType><<<elemGrid(), bs()>>>(
            hex_, areaX_.data(), areaY_.data(), areaZ_.data(), u.data(), v.data(), w.data(),
            int(prm_.skewAdvection), au_.data(), av_.data(), aw_.data());
        cudaCheckError();
        space_.restrict(au_);
        space_.restrict(av_);
        space_.restrict(aw_);
    }

    // out = (c M + nu K) q
    void applyViscous(RealType c, Vector& q, Vector& out)
    {
        space_.prolong(q);
        zero(out);
        pnsStiffnessKernel<KeyType, RealType><<<elemGrid(), bs()>>>(hex_, Ke_.data(), q.data(), out.data());
        cudaCheckError();
        space_.restrict(out);
        pnsViscousCombineKernel<RealType><<<grid(), bs()>>>(space_.isDof(), n_, c, prm_.nu, mass_.data(), q.data(),
                                                            out.data());
        cudaCheckError();
    }

    // out = D M^-1 G x
    void applyProjection(Vector& x, Vector& out)
    {
        space_.prolong(x);
        gradientOf(x);
        space_.prolong(gx_);
        space_.prolong(gy_);
        space_.prolong(gz_);
        divergenceOf(gx_, gy_, gz_, out);
    }

    // Jacobi-preconditioned CG in the reduced space. Stops when the residual
    // is tolerance * max(|b|, bRef). Returns the iteration count, or -1 on
    // breakdown or when maxIter is reached; x is returned prolonged.
    template<class Operator>
    int pcg(Operator&& A, const Vector& b, Vector& x, const Vector& diag, RealType bRef, bool zeroGuess)
    {
        const uint8_t* dof = space_.isDof();
        if (zeroGuess)
            copy(b, r_);
        else
        {
            A(x, Ad_);
            pnsAxpyKernel<RealType><<<grid(), bs()>>>(dof, n_, b.data(), RealType(-1), Ad_.data(), r_.data());
        }
        pnsJacobiKernel<RealType><<<grid(), bs()>>>(dof, n_, r_.data(), diag.data(), z_.data());
        cudaCheckError();
        copy(z_, d_);

        auto [rz, bb]       = space_.dot2(r_, z_, b, b);
        RealType rr         = space_.dot(r_, r_);
        const RealType stop = prm_.tolerance * std::max(std::sqrt(bb), bRef);
        int it              = 0;
        while (std::sqrt(rr) > stop)
        {
            if (it == prm_.maxIter) return -1;
            A(d_, Ad_);
            RealType dAd = space_.dot(d_, Ad_);
            if (!(dAd > RealType(0))) return -1;
            RealType alpha = rz / dAd;
            pnsAxpyKernel<RealType><<<grid(), bs()>>>(dof, n_, x.data(), alpha, d_.data(), x.data());
            pnsAxpyKernel<RealType><<<grid(), bs()>>>(dof, n_, r_.data(), -alpha, Ad_.data(), r_.data());
            pnsJacobiKernel<RealType><<<grid(), bs()>>>(dof, n_, r_.data(), diag.data(), z_.data());
            cudaCheckError();
            auto [rrNew, rzNew] = space_.dot2(r_, r_, r_, z_);
            pnsAxpyKernel<RealType><<<grid(), bs()>>>(dof, n_, z_.data(), rzNew / rz, d_.data(), d_.data());
            cudaCheckError();
            rr = rrNew;
            rz = rzNew;
            ++it;
        }
        space_.prolong(x);
        return it;
    }

    void buildGeometry()
    {
        const size_t m = hex_.count;
        areaX_.resize(12 * m);
        areaY_.resize(12 * m);
        areaZ_.resize(12 * m);
        Ke_.resize(64 * m);
        if (m == 0) return;
        domain_.cacheNodeCoordinates();
        const auto& x = domain_.getNodeX();
        const auto& y = domain_.getNodeY();
        const auto& z = domain_.getNodeZ();
        const KeyType* c[8];
        for (int i = 0; i < 8; ++i)
            c[i] = hex_.node[i] + hex_.first;
        precomputeAreaVectorsGpu<KeyType, RealType>(c[0], c[1], c[2], c[3], c[4], c[5], c[6], c[7], m, x.data(),
                                                    y.data(), z.data(), areaX_.data(), areaY_.data(), areaZ_.data());
        // CVFEM viscous element matrices with unit diffusivity, the same
        // operator the assembled solvers use.
        Vector ones(n_, RealType(1));
        constexpr int buildBlock = 256;
        cvfem_hex_matfree_build_lhs<KeyType, RealType, buildBlock>
            <<<int((m + buildBlock - 1) / buildBlock), buildBlock>>>(
                c[0], c[1], c[2], c[3], c[4], c[5], c[6], c[7], m, x.data(), y.data(), z.data(), ones.data(),
                areaX_.data(), areaY_.data(), areaZ_.data(), Ke_.data());
        cudaCheckError();
        cudaDeviceSynchronize();
    }

    void buildMassAndDiagonals()
    {
        pnsLumpedMassKernel<KeyType, RealType><<<elemGrid(), bs()>>>(
            hex_, domain_.getNodeX().data(), domain_.getNodeY().data(), domain_.getNodeZ().data(), mass_.data());
        pnsStiffnessDiagonalKernel<KeyType, RealType><<<elemGrid(), bs()>>>(hex_, Ke_.data(), diagK_.data());
        cudaCheckError();
        space_.restrict(mass_);
        space_.prolong(mass_);
        space_.restrict(diagK_);
        pnsProjectionDiagonalKernel<KeyType, RealType><<<elemGrid(), bs()>>>(
            hex_, areaX_.data(), areaY_.data(), areaZ_.data(), mass_.data(), diagP_.data());
        cudaCheckError();
        space_.restrict(diagP_);
    }

    void zero(Vector& a) { cudaMemsetAsync(a.data(), 0, a.size() * sizeof(RealType)); }
    void copy(const Vector& from, Vector& to)
    {
        cudaMemcpyAsync(to.data(), from.data(), n_ * sizeof(RealType), cudaMemcpyDeviceToDevice);
    }
    void scale(Vector& a, RealType s)
    {
        pnsAxpyKernel<RealType><<<grid(), bs()>>>(space_.isDof(), n_, a.data(), s - RealType(1), a.data(), a.data());
        cudaCheckError();
    }

    // At least one block: a rank may own no elements, and an empty grid is a launch error.
    int bs() const { return prm_.blockSize; }
    int grid() const { return std::max(1, int((n_ + bs() - 1) / bs())); }
    int elemGrid() const { return std::max(1, int((hex_.count + bs() - 1) / bs())); }

    Domain& domain_;
    Space space_;
    Params prm_;
    size_t n_;
    OwnedHexElements<KeyType> hex_{};
    int numSteps_ = 0;
    int velocityIters_[3] = {0, 0, 0};
    int pressureIters_    = 0;

    Vector um1_, vm1_, wm1_;          // velocity at the previous step (BDF2)
    Vector us_, vs_, ws_;             // u* then u**
    Vector au_, av_, aw_;             // advection N(u^n), restricted
    Vector aum1_, avm1_, awm1_;       // advection N(u^{n-1})
    Vector gx_, gy_, gz_;             // M^-1 G q
    Vector phi_, rhs_;
    Vector mass_, diagK_, diagP_, diagV_;
    Vector r_, z_, d_, Ad_, scratch_; // CG work vectors
    Vector omegaX_, omegaY_, omegaZ_;
    Vector areaX_, areaY_, areaZ_;    // sub-control face area vectors, 12 per owned element
    Vector Ke_;                       // 8x8 viscous matrix per owned element
};

} // namespace fem
} // namespace mars
