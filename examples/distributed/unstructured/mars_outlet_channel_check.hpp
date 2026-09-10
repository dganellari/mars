#pragma once

// GPT/Codex, 2026-09-10. Opt-in checks for the public channel integration fixture.
struct OutletChannelMoments
{
    double volume, u, v, w, p, u2, p2, trace, area;
};

struct OutletChannelPlus
{
    __host__ __device__ OutletChannelMoments operator()(const OutletChannelMoments& a,
                                                       const OutletChannelMoments& b) const
    {
        return {a.volume+b.volume, a.u+b.u, a.v+b.v, a.w+b.w, a.p+b.p,
                a.u2+b.u2, a.p2+b.p2, a.trace+b.trace, a.area+b.area};
    }
};

template<class KeyType, class RealType>
struct OutletChannelCheck
{
    cstone::DeviceVector<RealType> previous_u, previous_v, previous_w, scratch, residual;

    static void copy(const cstone::DeviceVector<RealType>& source, cstone::DeviceVector<RealType>& target)
    {
        target.resize(source.size());
        if (!source.empty())
            thrust::copy(thrust::device, thrust::device_pointer_cast(source.data()),
                         thrust::device_pointer_cast(source.data()+source.size()),
                         thrust::device_pointer_cast(target.data()));
    }

    void capture(NSStepper<KeyType, RealType, TetTag>& s)
    {
        copy(s.d_u, previous_u); copy(s.d_v, previous_v); copy(s.d_w, previous_w);
    }

    void verify(NSStepper<KeyType, RealType, TetTag>& s, int step, RealType dt)
    {
        const RealType expected_dt = step == 1 ? dt : RealType(2)*dt/RealType(3);
        require_outlet_correction(s, s.vmsCtx.bdf2 == (step > 1)
            && std::abs(s.vmsCtx.dtEff-expected_dt) <= 16*std::numeric_limits<RealType>::epsilon()*dt,
            "channel gate: wrong BDF startup coefficient");
        require_outlet_correction(s, s.d_u_nm1.size() == s.nodeCount
            && s.d_v_nm1.size() == s.nodeCount && s.d_w_nm1.size() == s.nodeCount,
            "channel gate: missing BDF history");
        const auto* own = s.ownershipMap().data();
        auto difference = [&](const cstone::DeviceVector<RealType>& a,
                              const cstone::DeviceVector<RealType>& b, bool owned_only) {
            const auto* ap = a.data(); const auto* bp = b.data();
            double local = thrust::transform_reduce(thrust::device,
                thrust::counting_iterator<size_t>(0), thrust::counting_iterator<size_t>(s.nodeCount),
                [=] __device__(size_t i) -> double {
                    if (owned_only && own[i] != 1) return 0.;
                    return isfinite(ap[i]) && isfinite(bp[i]) ? fabs(double(ap[i]-bp[i])) : INFINITY;
                }, 0., thrust::maximum<double>());
            double global = 0;
            MPI_Allreduce(&local, &global, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
            return global;
        };
        double history_error = difference(previous_u, s.d_u_nm1, true);
        history_error = std::max(history_error, difference(previous_v, s.d_v_nm1, true));
        history_error = std::max(history_error, difference(previous_w, s.d_w_nm1, true));
        require_outlet_correction(s, history_error == 0., "channel gate: BDF history advanced incorrectly");

        // Re-publish copies to check that the solved fields already have current ghost values.
        double halo_error = 0;
        for (const auto* field : {&s.d_u, &s.d_v, &s.d_w, &s.d_p})
        {
            copy(*field, scratch);
            s.domain.exchangeNodeHalo(scratch);
            halo_error = std::max(halo_error, difference(scratch, *field, false));
        }
        require_outlet_correction(s, halo_error == 0., "channel gate: stale solved-field halo");
        assemble_outlet_continuity(s, s.d_u.data(), s.d_v.data(), s.d_w.data(), residual);
        const auto norm = outlet_residual_norm(s, residual);
        double q_in = 0, q_out = 0;
        boundaryMassBalance(s, &s.vmsCtx, q_in, q_out);
        require_outlet_correction(s, norm.rms <= 1e-7 && norm.maximum <= 1e-7
            && std::abs(q_in+q_out) <= 1e-8 && std::abs(norm.sum-q_in-q_out) <= 1e-10
            && std::abs(q_in+double(s.Uinf)) <= 1e-10 && q_out > 0,
            "channel gate: continuity, source, or boundary balance failed");

        const auto* mass = s.d_massNode.data();
        const auto* u = s.d_u.data(); const auto* v = s.d_v.data(); const auto* w = s.d_w.data();
        const auto* p = s.d_p.data(); const auto* trace = s.d_pTraceOutlet.data();
        const auto* area = s.d_outletAreaScalar.data();
        const auto local = thrust::transform_reduce(thrust::device,
            thrust::counting_iterator<size_t>(0), thrust::counting_iterator<size_t>(s.nodeCount),
            [=] __device__(size_t i) -> OutletChannelMoments {
                if (own[i] != 1) return {};
                const double m = mass[i], a = area[i];
                return {m, m*u[i], m*v[i], m*w[i], m*p[i],
                    m*(u[i]*u[i]+v[i]*v[i]+w[i]*w[i]), m*p[i]*p[i], a*trace[i], a};
            }, OutletChannelMoments{}, OutletChannelPlus{});
        double values[9] = {local.volume, local.u, local.v, local.w, local.p,
                            local.u2, local.p2, local.trace, local.area}, sums[9] = {};
        MPI_Allreduce(values, sums, 9, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
        bool finite = true;
        for (double value : sums) finite &= std::isfinite(value);
        // Tet geometry is decoded from SFC keys, while opening areas come from the Exodus triangles.
        // Compare volume to the encoded/decoded public box, not to its pre-quantization volume.
        const auto& box = s.domain.getBoundingBox();
        using SfcKey = cstone::SfcKind<KeyType>;
        const KeyType lo_key = cstone::sfc3D<SfcKey>(RealType(0), RealType(0), RealType(0), box);
        const KeyType hi_key = cstone::sfc3D<SfcKey>(RealType(4), RealType(1), RealType(1), box);
        const auto [lx, ly, lz] = s.domain.sfcToPhysicalCoordinate(lo_key);
        const auto [hx, hy, hz] = s.domain.sfcToPhysicalCoordinate(hi_key);
        const double expected_volume = double(hx-lx)*double(hy-ly)*double(hz-lz);
        require_outlet_correction(s, finite && std::isfinite(expected_volume) && expected_volume > 0
            && std::abs(sums[0]-expected_volume) <= 1e-10
            && std::abs(sums[8]-1.) <= 1e-10
            && std::abs(sums[7]-double(s.outletPRef)*sums[8]) <= 1e-10,
            "channel gate: volume, outlet area, or frozen trace mean failed");
        if (s.rank == 0)
        {
            const auto flags = std::cout.flags(); const auto precision = std::cout.precision();
            std::cout << std::scientific << std::setprecision(16)
                << "[outlet-channel] step=" << step << " ranks=" << s.numRanks
                << " bdf=" << (s.vmsCtx.bdf2 ? 2 : 1) << " dt_eff=" << s.vmsCtx.dtEff
                << " rms=" << norm.rms << " max=" << norm.maximum
                << " q_in=" << q_in << " q_out=" << q_out << " residual_sum=" << norm.sum
                << " u_mean=" << sums[1]/sums[0] << " v_mean=" << sums[2]/sums[0]
                << " w_mean=" << sums[3]/sums[0] << " p_mean=" << sums[4]/sums[0]
                << " u_rms=" << std::sqrt(sums[5]/sums[0]) << " p_rms=" << std::sqrt(sums[6]/sums[0])
                << " trace_mean=" << sums[7]/sums[8]
                << " history_error=" << history_error << " halo_error=" << halo_error << '\n';
            std::cout.flags(flags); std::cout.precision(precision);
        }
    }
};
