#pragma once

#include "mars_outlet_fgmres.hpp"

// Included after the stepper's gradient, continuity, and collective failure helpers.
template<typename RealType>
__global__ void outlet_krylov_products_kernel(const RealType* basis, const RealType* vector,
                                              size_t stride, int n, double* products)
{
    __shared__ double partial[256];
    const auto* column = basis + size_t(blockIdx.x)*stride;
    double value = 0;
    for (int i = threadIdx.x; i < n; i += blockDim.x) value += double(column[i])*double(vector[i]);
    partial[threadIdx.x] = value;
    __syncthreads();
    for (int width = 128; width > 0; width /= 2)
    {
        if (threadIdx.x < width) partial[threadIdx.x] += partial[threadIdx.x + width];
        __syncthreads();
    }
    if (threadIdx.x == 0) products[blockIdx.x] = partial[0];
}

template<typename KeyType, typename RealType, typename ElementTag>
struct OutletKrylovOps
{
    using Stepper = NSStepper<KeyType, RealType, ElementTag>;
    Stepper& s;
    RealType h;
    int restart, n;
    size_t stride;
    std::vector<double> coefficients;
#ifdef MARS_ENABLE_HYPRE
    using PreparedSolver = mars::fem::HypreGMRESSolver<RealType, int, cstone::GpuTag>;
    std::unique_ptr<PreparedSolver> prepared_solver;
#endif

    OutletKrylovOps(Stepper& stepper, RealType step_scale, int depth)
        : s(stepper), h(step_scale), restart(depth), n(s.numOwnedDofs),
          stride(std::max(size_t(1), size_t(n))), coefficients(depth + 1)
    {
        s.d_outlet_krylov.resize(stride*size_t(2*restart + 5));
        s.d_outlet_krylov_products.resize(restart + 1);
        s.d_outlet_delta_u.resize(s.nodeCount);
        s.d_outlet_delta_v.resize(s.nodeCount);
        s.d_outlet_delta_w.resize(s.nodeCount);
#ifdef MARS_ENABLE_HYPRE
        if (s.outlet_preconditioner != 0)
        {
            const char* flex = std::getenv("MARS_HYPRE_FLEXGMRES");
            const char* precond = std::getenv("MARS_HYPRE_PRECOND");
            const bool jacobi = precond && std::string(precond) == "jacobi";
            require_outlet_correction(s, (!flex || std::string(flex) == "0")
                && (s.outlet_preconditioner != 2 || !jacobi),
                "prepared outlet requires plain inner GMRES; amg-cycle requires BoomerAMG");
            prepared_solver = std::make_unique<PreparedSolver>(MPI_COMM_WORLD, s.maxIter, s.tolerance,
                jacobi ? PreparedSolver::JACOBI : PreparedSolver::BOOMERAMG);
            prepared_solver->setVerbose(false);
            prepared_solver->enable_reuse(s.outlet_preconditioner == 2);
            if (s.outlet_profile.enabled()) prepared_solver->set_profile(&s.outlet_profile);
        }
#endif
    }

    RealType* solution() { return s.d_outlet_krylov.data(); }
    RealType* rhs() { return solution() + stride; }
    RealType* residual() { return solution() + 2*stride; }
    RealType* work() { return solution() + 3*stride; }
    RealType* basis(int i) { return solution() + size_t(4 + i)*stride; }
    RealType* direction(int i) { return solution() + size_t(5 + restart + i)*stride; }

    void check(cudaError_t error)
    {
        if (error == cudaSuccess) return;
        std::cerr << "ERROR: outlet Krylov: " << cudaGetErrorString(error) << '\n';
        MPI_Abort(MPI_COMM_WORLD, 1);
        std::abort();
    }

    void zero(RealType* target)
    {
        if (n > 0) check(cudaMemsetAsync(target, 0, size_t(n)*sizeof(RealType)));
    }

    void copy(const RealType* source, RealType* target)
    {
        if (n > 0) check(cudaMemcpyAsync(target, source, size_t(n)*sizeof(RealType), cudaMemcpyDeviceToDevice));
    }

    void scale(double alpha, RealType* vector)
    {
        const RealType a = RealType(alpha);
        thrust::for_each(thrust::device, thrust::counting_iterator<int>(0),
            thrust::counting_iterator<int>(n), [=] __device__(int i) { vector[i] *= a; });
        check(cudaGetLastError());
    }

    void axpy(double alpha, const RealType* x, RealType* y)
    {
        const RealType a = RealType(alpha);
        thrust::for_each(thrust::device, thrust::counting_iterator<int>(0),
            thrust::counting_iterator<int>(n), [=] __device__(int i) { y[i] += a*x[i]; });
        check(cudaGetLastError());
    }

    double norm(const RealType* vector)
    {
        double local = thrust::transform_reduce(thrust::device, thrust::counting_iterator<int>(0),
            thrust::counting_iterator<int>(n),
            [=] __device__(int i) -> double { return double(vector[i])*double(vector[i]); },
            0.0, thrust::plus<double>());
        check(cudaGetLastError());
        double global = 0;
        MPI_Allreduce(&local, &global, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
        return std::sqrt(global);
    }

    void orthogonalize(RealType* vector, int count, double* total)
    {
        std::fill(total, total + count, 0.0);
        const auto* v = basis(0);
        auto* products = s.d_outlet_krylov_products.data();
        const size_t column_stride = stride;
        for (int pass = 0; pass < 2; ++pass)
        {
            // Even an empty rank writes zero partials and participates in both reductions.
            outlet_krylov_products_kernel<RealType><<<count, 256>>>(v, vector, stride, n, products);
            check(cudaGetLastError());
            check(cudaMemcpy(coefficients.data(), products, size_t(count)*sizeof(double), cudaMemcpyDeviceToHost));
            MPI_Allreduce(MPI_IN_PLACE, coefficients.data(), count, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
            for (int j = 0; j < count; ++j) total[j] += coefficients[j];
            check(cudaMemcpy(products, coefficients.data(), size_t(count)*sizeof(double), cudaMemcpyHostToDevice));
            thrust::for_each(thrust::device, thrust::counting_iterator<int>(0),
                thrust::counting_iterator<int>(n), [=] __device__(int i) {
                    double value = vector[i];
                    for (int j = 0; j < count; ++j) value -= products[j]*double(v[size_t(j)*column_stride + i]);
                    vector[i] = RealType(value);
                });
            check(cudaGetLastError());
        }
    }

    void set_pressure(const RealType* owned_pressure)
    {
        auto* phi = s.d_phi.data();
        const auto* own = s.ownershipMap().data();
        const auto* dof = s.d_node_to_dof.data();
        const int n_owned = n;
        thrust::for_each(thrust::device, thrust::counting_iterator<size_t>(0),
            thrust::counting_iterator<size_t>(s.nodeCount), [=] __device__(size_t i) {
                const int row = dof[i];
                phi[i] = own[i] == 1 && row >= 0 && row < n_owned ? owned_pressure[row] : RealType(0);
            });
        check(cudaGetLastError());
        s.domain.exchangeNodeHalo(s.d_phi);
    }

    void build_rhs()
    {
        const auto* flux = s.d_outlet_residual.data();
        const auto* mass = s.d_massNode.data();
        const auto* own = s.ownershipMap().data();
        const auto* dof = s.d_node_to_dof.data();
        const int n_owned = n;
        const RealType step_scale = h;
        auto* b = rhs();
        zero(b);
        thrust::for_each(thrust::device, thrust::counting_iterator<size_t>(0),
            thrust::counting_iterator<size_t>(s.nodeCount), [=] __device__(size_t i) {
                const int row = dof[i];
                if (own[i] == 1 && row >= 0 && row < n_owned)
                    b[row] = -flux[i]/(step_scale*sqrt(mass[i]));
            });
        check(cudaGetLastError());
    }

    // Solve Apre*z = S^-1*v, S_ii=1/sqrt(V_i). Hypre may vary its iteration count.
    bool precondition(const RealType* v, RealType* z)
    {
        SolverProfile::Scope profile(s.outlet_profile, SolverProfile::Preconditioner);
        auto* b = s.d_outlet_rhs.data();
        const auto* mass = s.d_massNode.data();
        const auto* own = s.ownershipMap().data();
        const auto* dof = s.d_node_to_dof.data();
        const int n_owned = n;
        zero(b);
        thrust::for_each(thrust::device, thrust::counting_iterator<size_t>(0),
            thrust::counting_iterator<size_t>(s.nodeCount), [=] __device__(size_t i) {
                const int row = dof[i];
                if (own[i] == 1 && row >= 0 && row < n_owned) b[row] = sqrt(mass[i])*v[row];
            });
        check(cudaGetLastError());
        if (!s.d_outlet_solution.empty())
            check(cudaMemsetAsync(s.d_outlet_solution.data(), 0, s.d_outlet_solution.size()*sizeof(RealType)));
        int iterations = -2;
#ifdef MARS_ENABLE_HYPRE
        if (prepared_solver)
        {
            bool usable = prepared_solver->solve(s.Apre, s.d_outlet_rhs, s.d_outlet_solution,
                static_cast<int>(s.globalRowStart), static_cast<int>(s.globalRowEnd),
                0, static_cast<int>(s.numInteriorGlobal), s.d_localToGlobalDof);
            // Preserve the existing GMRES fallback; a direct cycle has no inner residual target.
            if (s.outlet_preconditioner == 1)
            {
                RealType accept_res = RealType(1e-6);
                if (const char* value = std::getenv("MARS_HYPRE_ACCEPT_RES"))
                { const double parsed = std::atof(value); if (parsed > 0) accept_res = RealType(parsed); }
                usable = usable || (prepared_solver->getLastFinalResidual() < accept_res
                                     && !prepared_solver->lastReturnedNullSolution());
            }
            iterations = usable ? prepared_solver->getLastIterations() : -2;
        }
        else
#endif
            iterations = solveOneComponent(s, s.d_outlet_rhs, s.d_outlet_solution,
                                           s.d_phi, s.Apre, KrylovHint::GMRES);
        require_outlet_correction(s, iterations >= 0, "Hypre failed inside the true-J preconditioner");
        s.lastPressureIters += iterations;
        copy(s.d_outlet_solution.data(), z);
        return true;
    }

    // S*J*phi = S*R(-h Q G_v phi, phi; delta Gbar=0, delta trace=0)/h.
    // Homogeneous flux evaluation avoids subtracting two nonzero residuals.
    void apply(const RealType* pressure, RealType* action)
    {
        set_pressure(pressure);
        compute_pressure_increment_gradient(s, false);
        const auto* own = s.ownershipMap().data();
        const auto* dof = s.d_node_to_dof.data();
        const auto* fixed = s.d_isBdryDof.data();
        const int n_owned = n;
        const RealType step_scale = h;
        auto* u = s.d_outlet_delta_u.data();
        auto* v = s.d_outlet_delta_v.data();
        auto* w = s.d_outlet_delta_w.data();
        const auto* gx = s.d_gradPhix.data();
        const auto* gy = s.d_gradPhiy.data();
        const auto* gz = s.d_gradPhiz.data();
        thrust::for_each(thrust::device, thrust::counting_iterator<size_t>(0),
            thrust::counting_iterator<size_t>(s.nodeCount), [=] __device__(size_t i) {
                const int row = dof[i];
                const bool free = own[i] == 1 && row >= 0 && row < n_owned && !fixed[row];
                u[i] = free ? -step_scale*gx[i] : RealType(0);
                v[i] = free ? -step_scale*gy[i] : RealType(0);
                w[i] = free ? -step_scale*gz[i] : RealType(0);
            });
        check(cudaGetLastError());
        s.domain.exchangeNodeHalo(s.d_outlet_delta_u);
        s.domain.exchangeNodeHalo(s.d_outlet_delta_v);
        s.domain.exchangeNodeHalo(s.d_outlet_delta_w);
        assemble_outlet_continuity_fields(s, u, v, w, s.d_phi.data(),
            static_cast<const RealType*>(nullptr), RealType(0), false, s.d_outlet_action);
        const auto* flux = s.d_outlet_action.data();
        const auto* mass = s.d_massNode.data();
        zero(action);
        thrust::for_each(thrust::device, thrust::counting_iterator<size_t>(0),
            thrust::counting_iterator<size_t>(s.nodeCount), [=] __device__(size_t i) {
                const int row = dof[i];
                if (own[i] == 1 && row >= 0 && row < n_owned)
                    action[row] = flux[i]/(step_scale*sqrt(mass[i]));
            });
        check(cudaGetLastError());
    }

    void report(int iterations, double relative_residual)
    {
        if (std::getenv("MARS_SOLVE_TRACE") && s.rank == 0)
            std::cout << "[outlet-krylov] iterations=" << iterations
                      << " true_relative_residual=" << relative_residual << '\n';
    }
};
