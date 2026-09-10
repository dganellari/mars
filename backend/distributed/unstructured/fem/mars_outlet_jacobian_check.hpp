#pragma once

// Opt-in public fixture gate. Snapshots stay on device; only scalar errors are reported.
template<typename KeyType, typename RealType, typename ElementTag>
void check_outlet_jacobian(OutletKrylovOps<KeyType, RealType, ElementTag>& ops)
{
    auto& s = ops.s;
    using Vector = cstone::DeviceVector<RealType>;
    std::array<Vector*, 12> fields{&s.d_u, &s.d_v, &s.d_w, &s.d_p, &s.d_pTraceOutlet,
        &s.vmsCtx.Gx, &s.vmsCtx.Gy, &s.vmsCtx.Gz, &s.vmsCtx.tauNode,
        &s.d_u_nm1, &s.d_v_nm1, &s.d_w_nm1};
    std::array<Vector, 12> saved;
    auto copy = [&](const Vector& source, Vector& target) {
        target.resize(source.size());
        if (!source.empty()) ops.check(cudaMemcpy(target.data(), source.data(),
            source.size()*sizeof(RealType), cudaMemcpyDeviceToDevice));
    };
    for (size_t j = 0; j < fields.size(); ++j) copy(*fields[j], saved[j]);
    auto state_error = [&]() {
        double local = 0;
        for (size_t j = 0; j < fields.size(); ++j)
        {
            const auto* a = fields[j]->data();
            const auto* b = saved[j].data();
            const double error = thrust::transform_reduce(thrust::device,
                thrust::counting_iterator<size_t>(0), thrust::counting_iterator<size_t>(saved[j].size()),
                [=] __device__(size_t i) -> double {
                    return isfinite(a[i]) && isfinite(b[i]) ? fabs(double(a[i])-double(b[i])) : INFINITY;
                }, 0.0, thrust::maximum<double>());
            local = std::max(local, error);
        }
        double global = 0;
        MPI_Allreduce(&local, &global, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
        return global;
    };
    ops.zero(ops.solution());
    ops.apply(ops.solution(), ops.work());
    require_outlet_correction(s, ops.norm(ops.work()) == 0, "Jacobian gate: J*0 is nonzero");
    const auto* own = s.ownershipMap().data();
    const auto* dof = s.d_node_to_dof.data();
    const auto* x = s.domain.getNodeX().data();
    const auto* y = s.domain.getNodeY().data();
    const auto* z = s.domain.getNodeZ().data();
    const int n_owned = s.numOwnedDofs;
    auto* direction = ops.solution();
    thrust::for_each(thrust::device, thrust::counting_iterator<size_t>(0),
        thrust::counting_iterator<size_t>(s.nodeCount), [=] __device__(size_t i) {
            const int row = dof[i];
            if (own[i] == 1 && row >= 0 && row < n_owned)
                direction[row] = 1 + RealType(.2)*x[i] - RealType(.1)*y[i]
                    + RealType(.3)*z[i] + sin(x[i] + y[i]*z[i]);
        });
    ops.apply(direction, ops.basis(0));
    const double action_norm = ops.norm(ops.basis(0));
    require_outlet_correction(s, std::isfinite(action_norm) && action_norm > 0,
                             "Jacobian gate: nonfinite or zero test action");
    ops.copy(direction, ops.direction(0));
    ops.scale(2, ops.direction(0));
    ops.apply(ops.direction(0), ops.work());
    ops.axpy(-2, ops.basis(0), ops.work());
    require_outlet_correction(s, ops.norm(ops.work())/action_norm < 1e-12,
                             "Jacobian gate: homogeneous action is not linear");
    ops.apply(direction, ops.basis(0));
    require_outlet_correction(s, state_error() == 0, "Jacobian gate: apply changed frozen or physical state");

    auto* u = s.d_u.data(); auto* v = s.d_v.data(); auto* w = s.d_w.data(); auto* p = s.d_p.data();
    const auto* bu = saved[0].data(); const auto* bv = saved[1].data();
    const auto* bw = saved[2].data(); const auto* bp = saved[3].data();
    const auto* du = s.d_outlet_delta_u.data(); const auto* dv = s.d_outlet_delta_v.data();
    const auto* dw = s.d_outlet_delta_w.data(); const auto* phi = s.d_phi.data();
    const auto* mass = s.d_massNode.data();
    const auto* base_residual = s.d_outlet_residual.data();
    const auto* expected = ops.basis(0);
    const RealType h = ops.h;
    auto* error = ops.work();
    for (double value : {1e-3, 1e-4, 1e-5})
    {
        const RealType epsilon = RealType(value);
        thrust::for_each(thrust::device, thrust::counting_iterator<size_t>(0),
            thrust::counting_iterator<size_t>(s.nodeCount), [=] __device__(size_t i) {
                if (own[i] != 1) return;
                u[i] = bu[i] + epsilon*du[i]; v[i] = bv[i] + epsilon*dv[i];
                w[i] = bw[i] + epsilon*dw[i]; p[i] = bp[i] + epsilon*phi[i];
            });
        s.domain.exchangeNodeHalo(s.d_u); s.domain.exchangeNodeHalo(s.d_v);
        s.domain.exchangeNodeHalo(s.d_w); s.domain.exchangeNodeHalo(s.d_p);
        assemble_outlet_continuity(s, u, v, w, s.d_outlet_trial_residual);
        const auto* trial = s.d_outlet_trial_residual.data();
        ops.zero(error);
        thrust::for_each(thrust::device, thrust::counting_iterator<size_t>(0),
            thrust::counting_iterator<size_t>(s.nodeCount), [=] __device__(size_t i) {
                const int row = dof[i];
                if (own[i] == 1 && row >= 0 && row < n_owned)
                    error[row] = (trial[i]-base_residual[i])/(epsilon*h*sqrt(mass[i])) - expected[row];
            });
        const double relative_error = ops.norm(error)/action_norm;
        require_outlet_correction(s, relative_error < 1e-8,
                                 "Jacobian gate: action disagrees with the full correction finite difference");
        if (s.rank == 0)
            std::cout << "[outlet-jacobian] bdf=" << (s.vmsCtx.bdf2 ? 2 : 1)
                      << " epsilon=" << value << " relative_error=" << relative_error << '\n';
    }
    for (int j = 0; j < 4; ++j) copy(saved[j], *fields[j]);
    require_outlet_correction(s, state_error() == 0, "Jacobian gate: perturbation changed frozen state");
}
