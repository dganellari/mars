#ifndef MARS_ST_CG_HPP
#define MARS_ST_CG_HPP

#include <err.h>
#include <algorithm>
#include <cassert>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <map>
#include <numeric>
#include <string>

#ifdef WITH_KOKKOS_KERNELS
#include <KokkosBlas1_axpby.hpp>
#include <KokkosBlas1_dot.hpp>
#include <KokkosBlas1_mult.hpp>
#include <KokkosBlas1_nrm2.hpp>
#include <KokkosBlas1_nrm2_squared.hpp>
#include <KokkosBlas1_scal.hpp>
#include <KokkosSparse_spmv.hpp>
#include <Kokkos_Core.hpp>
#include "KokkosSparse_CrsMatrix.hpp"
#include "KokkosSparse_spmv_spec.hpp"

#include "mars_matrix_free_operator.hpp"

// Note on Blas1 Mult. We want to do z[i] = b*z[i] + a*x[i]*y[i] --> use as
// KokkosBlas::mult(b, z, a, x, y);
namespace mars {
    // using SparseMatrix = KokkosSparse::CrsMatrix<Real, Integer, Kokkos::Serial>;
    using VecType = mars::ViewVectorType<Real>;

    template <class OperatorA, class OperatorP>
    bool bcg_stab(OperatorA &A, OperatorP &P, const VecType &b, Integer max_iter, VecType &x_0, Integer &num_iter) {
        const auto &comm = A.comm();
        const size_t N = x_0.extent(0);

        std::cout << "Max iter: " << max_iter << std::endl;

        // FIXME can any of these buffers be avoided?
        VecType r_0("r_0", N);
        VecType r_1("r_1", N);
        VecType p_0("p_0", N);
        VecType p_1("p_1", N);
        VecType v_0("v_0", N);
        VecType v_1("v_1", N);
        VecType x_1("x_1", N);
        VecType r_hat("r_hat", N);
        VecType y("y", N);
        VecType z("z", N);
        VecType s("s", N);
        VecType h("h", N);
        VecType t("t", N);
        VecType tP("tP", N);

        Real rho_0 = 1.0;
        Real rho_1;
        Real alpha_0 = 1.0;
        Real beta;
        Real omega_0 = 1.0;
        Real omega_1;

        Real TOL = 1.0e-8;

        // Compute initial residual
        A.apply(x_0, r_0);                     // r_0 =  A*x_0
        KokkosBlas::axpby(1.0, b, -1.0, r_0);  // r_0 = b - A*x_0

        // Compute initial r_hat
        Kokkos::deep_copy(r_hat, r_0);
        double eps2_d = std::sqrt(comm.sum(KokkosBlas::nrm2_squared(b)));
        assert(!std::isnan(eps2_d));

        long count = 0;
        while (count < max_iter) {
            count++;
            rho_1 = comm.sum(KokkosBlas::dot(r_hat, r_0));  // rho = (r_hat)' * r_0
            assert(!std::isnan(rho_1));

            if (rho_1 == 0) {
                std::cerr << "METHOD FAILS: rho = 0" << std::endl;
                break;
            }

            if (count == 1.0) {
                Kokkos::deep_copy(p_1, r_0);
            } else {
                beta = (rho_1 / rho_0) * (alpha_0 / omega_0);
                assert(!std::isnan(beta));

                KokkosBlas::scal(p_1, omega_0, v_0);
                KokkosBlas::axpby(1.0, p_0, -1.0, p_1);
                KokkosBlas::axpby(1.0, r_0, beta, p_1);
            }

            P.precondition_right(p_1, y);  // y = P*p_0
            A.apply(y, v_1);

            Real alpha_1 = 0.0;
            Real dot_rv = comm.sum(KokkosBlas::dot(r_hat, v_1));

            assert(!std::isnan(dot_rv));

            alpha_1 = rho_1 / dot_rv;

            assert(!std::isnan(alpha_1));

            KokkosBlas::scal(h, alpha_1, y);
            KokkosBlas::axpby(1.0, x_0, 1.0, h);

            // compute s
            KokkosBlas::scal(s, alpha_1, v_1);
            KokkosBlas::axpby(1.0, r_0, -1.0, s);

            Real norm_s = std::sqrt(comm.sum(KokkosBlas::nrm2_squared(s)));

            if (norm_s < TOL) {
                Kokkos::deep_copy(x_0, h);

                break;
            }

            // P.apply(s, z);  // z = P*s;
            P.precondition_left(s, z);  // z = P*s;
            A.apply(z, t);              // t = A*z;

            P.precondition_right(t, tP);

            // FIXME make one reduce instead of two
            omega_1 = comm.sum(KokkosBlas::dot(tP, s)) / comm.sum(KokkosBlas::dot(tP, t));

            assert(!std::isnan(omega_1));

            // compute new x
            KokkosBlas::scal(x_1, omega_1, z);
            KokkosBlas::axpby(1.0, h, 1.0, x_1);

            KokkosBlas::scal(r_1, omega_1, t);
            KokkosBlas::axpby(1.0, s, -1.0, r_1);

            Real stop_criteria = std::sqrt(comm.sum(KokkosBlas::nrm2_squared(r_1))) / eps2_d;

            assert(!std::isnan(stop_criteria));

            if (count == 1 || count % 100 == 0) {
                std::cout << count << " " << stop_criteria << std::endl;
            }

            if (stop_criteria < TOL) {
                Kokkos::deep_copy(x_0, x_1);
                break;
            }

            Kokkos::deep_copy(r_hat, r_1);
            Kokkos::deep_copy(p_0, p_1);
            Kokkos::deep_copy(r_0, r_1);
            Kokkos::deep_copy(v_0, v_1);
            Kokkos::deep_copy(x_0, x_1);

            // std::swap(r_hat, r_1);
            // std::swap(p_0, p_1);
            // std::swap(r_0, r_1);
            // std::swap(v_0, v_1);
            // std::swap(x_0, x_1);

            alpha_0 = alpha_1;
            omega_0 = omega_1;
            rho_0 = rho_1;
        }

        // std::swap(x_0, x_1);

        num_iter = count;
        A.apply(x_0, r_0);               // r_0 =  A*x_0
        KokkosBlas::axpy(-1.0, b, r_0);  // r_0 = b - A*x_0
        std::cout << "norm(Residual): " << std::sqrt(comm.sum(KokkosBlas::nrm2_squared(r_0))) << std::endl;
        std::cout << "Iter: " << count << std::endl;

        return num_iter < max_iter;
    }

}  // namespace mars

#endif  // MARS_ST_CG_HPP
#endif