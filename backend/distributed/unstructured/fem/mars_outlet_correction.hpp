#pragma once

#include "mars_outlet_iteration.hpp"

// Included by mars_ns_pump_solver.hpp after its flux and pressure helpers.
// The trace, reconstructed gradient, and momentum coefficients stay frozen here.

struct OutletResidualMoments
{
    double square, volume, sum, absolute;
};

struct OutletResidualPlus
{
    __host__ __device__ OutletResidualMoments operator()(OutletResidualMoments a,
                                                       OutletResidualMoments b) const
    {
        return {a.square + b.square, a.volume + b.volume,
                a.sum + b.sum, a.absolute + b.absolute};
    }
};

struct OutletResidualNorm
{
    double rms, maximum, sum, absolute, volume;
};

template<typename KeyType, typename RealType, typename ElementTag>
void require_outlet_correction(NSStepper<KeyType, RealType, ElementTag>& s,
                             bool local_ok, const char* message)
{
    int local_bad = local_ok ? 0 : 1, bad = 0;
    MPI_Allreduce(&local_bad, &bad, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
    if (bad)
    {
        if (s.rank == 0) std::cerr << "Average-pressure outlet: " << message << '\n';
        MPI_Abort(MPI_COMM_WORLD, 1);
        std::abort();
    }
}

template<typename KeyType, typename RealType, typename ElementTag>
OutletResidualNorm outlet_residual_norm(NSStepper<KeyType, RealType, ElementTag>& s,
                                      const cstone::DeviceVector<RealType>& residual)
{
    const auto* own = s.ownershipMap().data();
    const auto* dof = s.d_node_to_dof.data();
    const auto* mass = s.d_massNode.data();
    const auto* r = residual.data();
    const int n_owned = s.numOwnedDofs;
    auto first = thrust::counting_iterator<size_t>(0);
    OutletResidualMoments local = thrust::transform_reduce(thrust::device, first, first + s.nodeCount,
        [own, dof, mass, r, n_owned] __device__(size_t i) -> OutletResidualMoments {
            if (own[i] != 1 || dof[i] < 0 || dof[i] >= n_owned) return {0, 0, 0, 0};
            double q = double(r[i]), v = double(mass[i]);
            if (!(v > 0) || !isfinite(v) || !isfinite(q)) return {INFINITY, 0, 0, 0};
            return {q * q / v, v, q, fabs(q)};
        }, OutletResidualMoments{0, 0, 0, 0}, OutletResidualPlus{});
    double local_max = thrust::transform_reduce(thrust::device, first, first + s.nodeCount,
        [own, dof, mass, r, n_owned] __device__(size_t i) -> double {
            if (own[i] != 1 || dof[i] < 0 || dof[i] >= n_owned) return 0;
            return mass[i] > RealType(0) ? fabs(double(r[i]) / double(mass[i])) : INFINITY;
        }, 0.0, thrust::maximum<double>());
    double input[4] = {local.square, local.volume, local.sum, local.absolute}, global[4] = {};
    double maximum = 0;
    MPI_Allreduce(input, global, 4, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    MPI_Allreduce(&local_max, &maximum, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
    return {global[1] > 0 ? std::sqrt(global[0] / global[1]) : INFINITY,
            maximum, global[2], global[3], global[1]};
}

template<typename KeyType, typename RealType, typename ElementTag>
void run_outlet_pressure_correction(NSStepper<KeyType, RealType, ElementTag>& s, RealType dt, RealType rho)
{
    using Stepper = NSStepper<KeyType, RealType, ElementTag>;
    require_outlet_correction(s,
        std::is_same_v<ElementTag, TetTag> && s.bcKind == Stepper::BCKind::Pump
        && s.solverKind == SolverKind::Hypre && s.pressureSolve == PressureSolveKind::DDT
        && s.useVMSStab && s.useRcImplicit && !s.useRcOnly && !s.useRcBlend
        && !s.useFemProjection && !s.usePSPG && !s.useOpeningFluxSource
        && !s.fluxNeumann && !s.useOpenFaceNormalProj && !s.rotationalPressureCorrection
        && s.pumpDp == RealType(0) && s.relaxMass == RealType(1) && s.outletDoNothing
        && s.pressurePinDof < 0 && s.outletBeta <= RealType(1)
        && std::isfinite(s.outletBeta) && std::isfinite(s.outletPRef)
        && !s.implicitAdvection && s.relaxU > RealType(0) && s.relaxU <= RealType(1)
        && std::isfinite(s.outlet_relative_tolerance) && s.outlet_relative_tolerance > RealType(0)
        && s.outlet_relative_tolerance < RealType(1)
        && std::isfinite(s.outlet_divergence_tolerance) && s.outlet_divergence_tolerance > RealType(0)
        && std::isfinite(s.outlet_flux_tolerance) && s.outlet_flux_tolerance > RealType(0)
        && std::isfinite(s.outlet_max_damping) && s.outlet_max_damping > RealType(0)
        && s.outlet_max_damping <= RealType(1)
        && std::getenv("MARS_HYPRE_USE_DDT") == nullptr
        && std::getenv("MARS_FEMGRAM_SOLVE") == nullptr
        && s.outlet_max_corrections > 0 && std::isfinite(dt) && dt > RealType(0)
        && std::isfinite(rho) && rho > RealType(0),
        "unsupported correction configuration");
    require_outlet_correction(s, s.lastUIters >= 0 && s.lastVIters >= 0 && s.lastWIters >= 0,
                             "momentum solve failed before pressure correction");
    buildVmsFluxCtx(s, dt, rho, s.vmsCtx);
    validate_outlet_continuity(s);
    require_outlet_correction(s, s.vmsCtx.valid && s.vmsCtx.keepSmooth,
                             "the full frozen VMS gradient difference is required");
    refresh_outlet_pressure_operator(s, rho);
    const RealType h = s.vmsCtx.dtEff / rho;
    const auto* own = s.ownershipMap().data();
    const auto* dof = s.d_node_to_dof.data();
    const auto* fixed = s.d_isBdryDof.data();
    const int n_owned = s.numOwnedDofs;
    const int blocks = int((s.nodeCount + s.blockSize - 1) / s.blockSize);
    auto copy = [](const auto& source, auto& target) {
        target.resize(source.size());
        if (!source.empty())
            thrust::copy(thrust::device, thrust::device_pointer_cast(source.data()),
                         thrust::device_pointer_cast(source.data() + source.size()),
                         thrust::device_pointer_cast(target.data()));
    };
    auto publish = [&] {
        s.domain.exchangeNodeHalo(s.d_u);
        s.domain.exchangeNodeHalo(s.d_v);
        s.domain.exchangeNodeHalo(s.d_w);
        s.domain.exchangeNodeHalo(s.d_p);
    };
    copy(s.d_uStarStar, s.d_u); copy(s.d_vStarStar, s.d_v); copy(s.d_wStarStar, s.d_w);
    publish();
    s.d_outlet_rhs.resize(s.numOwnedDofs);
    s.d_outlet_solution.resize(s.numTotalDofs);
    assemble_outlet_continuity(s, s.d_u.data(), s.d_v.data(), s.d_w.data(), s.d_outlet_residual);
    auto norm = outlet_residual_norm(s, s.d_outlet_residual);
    const double rms_limit = std::max(double(s.outlet_divergence_tolerance),
                                     double(s.outlet_relative_tolerance) * norm.rms);
    const double max_limit = std::max(double(s.outlet_divergence_tolerance),
                                     double(s.outlet_relative_tolerance) * norm.maximum);
    s.last_outlet_corrections = 0;
    s.last_outlet_damping = 0;
    s.lastPressureIters = 0;
    s.lastGradPhiRms = 0;
    if (s.nodeCount > 0)
        thrust::fill(thrust::device, thrust::device_pointer_cast(s.d_phi.data()),
                     thrust::device_pointer_cast(s.d_phi.data() + s.nodeCount), RealType(0));

    for (int k = 0; ; ++k)
    {
        double q_in = 0, q_out = 0;
        boundaryMassBalance(s, &s.vmsCtx, q_in, q_out);
        const double flux_scale = std::abs(q_in) + std::abs(q_out);
        const double roundoff = std::max(1e-9, 512.0 * double(std::numeric_limits<RealType>::epsilon()));
        require_outlet_correction(s, std::isfinite(norm.rms) && std::isfinite(norm.maximum)
            && std::isfinite(q_in) && std::isfinite(q_out)
            && std::abs(norm.sum - q_in - q_out) <= double(s.outlet_flux_tolerance)
                + roundoff * (norm.absolute + flux_scale),
            "continuity sum disagrees with boundary flux, or the residual is nonfinite");
        s.last_outlet_residual_rms = RealType(norm.rms);
        s.last_outlet_residual_max = RealType(norm.maximum);
        s.last_outlet_balance = RealType(q_in + q_out);
        if (norm.rms <= rms_limit && norm.maximum <= max_limit
            && std::abs(q_in + q_out) <= double(s.outlet_flux_tolerance)
                + double(s.outlet_relative_tolerance) * flux_scale) break;
        require_outlet_correction(s, k < s.outlet_max_corrections,
                                 "full continuity residual did not converge within --outlet-max-corrections");

        thrust::fill(thrust::device, thrust::device_pointer_cast(s.d_outlet_rhs.data()),
                     thrust::device_pointer_cast(s.d_outlet_rhs.data() + s.d_outlet_rhs.size()), RealType(0));
        thrust::fill(thrust::device, thrust::device_pointer_cast(s.d_outlet_solution.data()),
                     thrust::device_pointer_cast(s.d_outlet_solution.data() + s.d_outlet_solution.size()), RealType(0));
        if (blocks > 0)
            buildPressureRhsKernel<RealType><<<blocks, s.blockSize>>>(
                s.d_outlet_residual.data(), dof, own, RealType(1) / h,
                s.d_outlet_rhs.data(), s.nodeCount, n_owned);
        int iterations = solveOneComponent(s, s.d_outlet_rhs, s.d_outlet_solution,
                                          s.d_phi, s.Apre, KrylovHint::GMRES);
        require_outlet_correction(s, iterations >= 0, "Hypre failed to solve the correction system");
        s.lastPressureIters += iterations;
        s.domain.exchangeNodeHalo(s.d_phi);
        compute_pressure_increment_gradient(s);
        copy(s.d_u, s.d_outlet_base_u); copy(s.d_v, s.d_outlet_base_v);
        copy(s.d_w, s.d_outlet_base_w); copy(s.d_p, s.d_outlet_base_p);
        auto trial = [&](RealType omega) {
            auto* u = s.d_u.data(); auto* v = s.d_v.data(); auto* w = s.d_w.data(); auto* p = s.d_p.data();
            const auto* bu = s.d_outlet_base_u.data(); const auto* bv = s.d_outlet_base_v.data();
            const auto* bw = s.d_outlet_base_w.data(); const auto* bp = s.d_outlet_base_p.data();
            const auto* gx = s.d_gradPhix.data(); const auto* gy = s.d_gradPhiy.data();
            const auto* gz = s.d_gradPhiz.data(); const auto* phi = s.d_phi.data();
            thrust::for_each(thrust::device, thrust::counting_iterator<size_t>(0),
                thrust::counting_iterator<size_t>(s.nodeCount),
                [=] __device__(size_t i) {
                    if (own[i] != 1 || dof[i] < 0 || dof[i] >= n_owned) return;
                    const RealType a = fixed[dof[i]] ? RealType(0) : omega * h;
                    u[i] = bu[i] - a * gx[i]; v[i] = bv[i] - a * gy[i]; w[i] = bw[i] - a * gz[i];
                    p[i] = bp[i] + omega * phi[i];
                });
            publish();
            assemble_outlet_continuity(s, s.d_u.data(), s.d_v.data(), s.d_w.data(), s.d_outlet_trial_residual);
        };
        // R(omega) is affine for this frozen step. Minimize its volume-weighted square norm.
        trial(RealType(1));
        const auto* r = s.d_outlet_residual.data();
        const auto* t = s.d_outlet_trial_residual.data();
        const auto* mass = s.d_massNode.data();
        auto first = thrust::counting_iterator<size_t>(0);
        auto products = thrust::transform_reduce(thrust::device, first, first + s.nodeCount,
            [=] __device__(size_t i) -> OutletResidualMoments {
                if (own[i] != 1 || dof[i] < 0 || dof[i] >= n_owned) return {0, 0, 0, 0};
                double delta = double(t[i]) - double(r[i]), volume = double(mass[i]);
                return {double(r[i]) * delta / volume, delta * delta / volume, 0, 0};
            }, OutletResidualMoments{0, 0, 0, 0}, OutletResidualPlus{});
        double local[2] = {products.square, products.volume}, global[2] = {};
        MPI_Allreduce(local, global, 2, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
        RealType omega = RealType(outlet_correction_damping(global[0], global[1], double(s.outlet_max_damping)));
        require_outlet_correction(s, omega > RealType(0),
            "compact correction has no descent direction; a true-J Krylov solve is required");
        bool accepted = false;
        OutletResidualNorm next{};
        for (int backtrack = 0; backtrack < 16; ++backtrack)
        {
            trial(omega);
            next = outlet_residual_norm(s, s.d_outlet_trial_residual);
            if (outlet_correction_contracts(norm.rms, next.rms, double(omega)))
            { accepted = true; break; }
            omega *= RealType(0.5);
        }
        require_outlet_correction(s, accepted, "damped correction failed the measured contraction gate");
        s.last_outlet_damping = omega;
        s.last_outlet_corrections = k + 1;
        s.lastGradPhiRms *= omega;
        if (std::getenv("MARS_SOLVE_TRACE") && s.rank == 0)
            std::cout << "  [outlet-correction] k=" << k + 1 << " omega=" << omega
                      << " contraction=" << next.rms / norm.rms << " div-rms=" << next.rms << '\n';
        s.d_outlet_residual.swap(s.d_outlet_trial_residual);
        if (s.nodeCount > 0)
            thrust::transform(thrust::device, thrust::device_pointer_cast(s.d_phi.data()),
                              thrust::device_pointer_cast(s.d_phi.data() + s.nodeCount),
                              thrust::device_pointer_cast(s.d_phi.data()),
                              [omega] __device__(RealType value) { return omega * value; });
        norm = next;
    }
    // Legacy reports remain distinct from the full boundary-inclusive acceptance norm above.
    divMaxAndRmsOwned(s, s.d_u, s.d_v, s.d_w, s.lastDivMax, s.lastDivRms);
    divMaxVmsOwned(s, dt, rho, s.lastDivRC, s.lastDivRCRms);
    const cudaError_t error = cudaDeviceSynchronize();
    require_outlet_correction(s, error == cudaSuccess, "CUDA failure in pressure correction");
}
