#ifndef MARS_ST_CG_HPP
#define MARS_ST_CG_HPP

#include <algorithm>
#include <cassert>
#include <cstdlib>
#include <err.h>
#include <fstream>
#include <iostream>
#include <map>
#include <numeric>
#include <string>

#include "KokkosSparse_CrsMatrix.hpp"
#include "KokkosSparse_spmv_spec.hpp"
#include <KokkosBlas1_axpby.hpp>
#include <KokkosBlas1_dot.hpp>
#include <KokkosBlas1_mult.hpp>
#include <KokkosBlas1_nrm2.hpp>
#include <KokkosBlas1_scal.hpp>
#include <KokkosSparse_spmv.hpp>
#include <Kokkos_Core.hpp>

#include "mars_matrix_free_operator.hpp"

// Note on Blas1 Mult. We want to do z[i] = b*z[i] + a*x[i]*y[i] --> use as
// KokkosBlas::mult(b, z, a, x, y);
namespace mars {
// using SparseMatrix = KokkosSparse::CrsMatrix<Real, Integer, Kokkos::Serial>;
using KokkosVector = Kokkos::View<Real *>;

void BiCGSTAB(Operator &A, Operator &P, const KokkosVector &b,
              KokkosVector &x_0, const size_t &N, Integer max_iter,
              double &omega_precon, KokkosVector &result, Real &num_iter

) {

  // CSVWriter csv;
  // std::string path_error =
  // "/Users/liudmilakaragyaur/code/space_time_error_estimator/results/csv/error_bcgstab/rel_err_init_zero.csv";
  // std::string path_error = "err_70.csv";
  // std::string iterations = "rel_error_zero.csv";
  std::cout << "Max iter: " << max_iter << std::endl;

  KokkosVector r_0("r_0", N);
  KokkosVector r_1("r_1", N);
  KokkosVector p_0("p_0", N);
  KokkosVector p_1("p_1", N);
  KokkosVector v_0("v_0", N);
  KokkosVector v_1("v_1", N);
  KokkosVector x_1("x_1", N);
  KokkosVector r_hat("r_hat", N);
  KokkosVector y("y", N);
  KokkosVector z("z", N);
  KokkosVector s("s", N);
  KokkosVector h("h", N);
  KokkosVector t("t", N);
  KokkosVector tP("tP", N);

  Real rho_0 = 1.0;
  Real rho_1;
  Real alpha_0 = 1.0;
  Real beta;
  Real omega_0 = 1.0;
  Real omega_1;

  Real TOL = 1.0e-10;

  // Compute initial residual

  A.apply(x_0, r_0);                    // r_0 =  A*x_0
  KokkosBlas::axpby(1.0, b, -1.0, r_0); // r_0 = b - A*x_0

  // Compute initial r_hat
  Kokkos::deep_copy(r_hat, r_0);
  double eps2_d = KokkosBlas::nrm2(b);

  std::vector<Real> rel_error;

  double count = 0.0;
  while (count < max_iter) {
    count++;
    rho_1 = KokkosBlas::dot(r_hat, r_0); // rho = (r_hat)' * r_0
    if (rho_1 == 0) {
      std::cerr << "METHOD FAILS: rho = 0" << std::endl;
      break;
    }

    if (count == 1.0) {
      Kokkos::deep_copy(p_1, r_0);
    } else {
      beta = 0.0;
      beta = (rho_1 / rho_0) * (alpha_0 / omega_0);
      KokkosBlas::scal(p_1, omega_0, v_0);
      KokkosBlas::axpby(1.0, p_0, -1.0, p_1);
      KokkosBlas::axpby(1.0, r_0, beta, p_1);
    }

    P.apply(p_1, y); // y = P*p_0
    A.apply(y, v_1);

    Real alpha_1 = 0.0;
    alpha_1 = rho_1 / KokkosBlas::dot(r_hat, v_1);

    KokkosBlas::scal(h, alpha_1, y);
    KokkosBlas::axpby(1.0, x_0, 1.0, h);

    // compute s
    KokkosBlas::scal(s, alpha_1, v_1);
    KokkosBlas::axpby(1.0, r_0, -1.0, s);

    rel_error.push_back(KokkosBlas::nrm2(s));

    // if (max_iter == 6719) {

    //   FILE *file =
    //       fopen("/Users/liudmilakaragyaur/code/space_time_error_estimator/"
    //             "test_3D/bicgstab/rel_error_close.csv",
    //             "a");
    //   fprintf(file, "%e \n", KokkosBlas::nrm2(s));
    //   fclose(file);
    // }

    if (KokkosBlas::nrm2(s) < TOL) {

      Kokkos::deep_copy(x_0, h);

      break;
    }

    P.apply(s, z); // z = P*s;
    A.apply(z, t); // t = A*z;

    P.apply(t, tP);
    omega_1 = KokkosBlas::dot(tP, s) / KokkosBlas::dot(tP, t);

    // compute new x
    KokkosBlas::scal(x_1, omega_1, z);
    KokkosBlas::axpby(1.0, h, 1.0, x_1);

    KokkosBlas::scal(r_1, omega_1, t);
    KokkosBlas::axpby(1.0, s, -1.0, r_1);

    Real stop_criteria = KokkosBlas::nrm2(r_1) / eps2_d;
    // rel_error.push_back(stop_criteria);

    if (stop_criteria < TOL) {

      Kokkos::deep_copy(x_0, x_1);

      break;
    }

    Kokkos::deep_copy(r_hat, r_1);
    Kokkos::deep_copy(p_0, p_1);
    Kokkos::deep_copy(r_0, r_1);
    Kokkos::deep_copy(v_0, v_1);
    Kokkos::deep_copy(x_0, x_1);

    alpha_0 = alpha_1;
    omega_0 = omega_1;
    rho_0 = rho_1;
  }
  Kokkos::deep_copy(result, x_0);

  num_iter = count;
  A.apply(result, r_0);                 // r_0 =  A*x_0
  KokkosBlas::axpby(1.0, b, -1.0, r_0); // r_0 = b - A*x_0
  std::cout << "RESIDUE: " << KokkosBlas::nrm2(r_0) << std::endl;
  // if (max_iter == 6719){
  //     csv.write_value(iterations,rel_error);
  // }
  std::cout << "Iter: " << count << std::endl;
}

// void preconCG(
// const SparseMatrix &A,
// const SparseMatrix &M,
// const KokkosVector &b,
// std::vector<bool> &node_is_boundary,
// const Real &dt,
// KokkosVector &init_guess,
// KokkosVector &result,
// Real &num_iter) {

//     Integer N =  A.numRows();
//     // Kokkos::View<Real*> x_0("x_0", N );
//     Kokkos::View<Real*> a("diag_A", N );
//     Kokkos::View<Real*> Precon("precond", N );

//     // ImplicitEulerOperator S(A,M,dt,N);
//     MatrixWrapper S(A);

//     for(Integer i = 0; i < N; ++i) {
//         auto row = A.row(i);
//         auto n_vals = row.length;
//         for(Integer k = 0; k < n_vals; ++k) {
//             auto j = row.colidx(k);
//             if (i==j){
//                 a(i) = 1./row.value(k);
//                 // a(i) = row.value(k);
//                 // std::cout << a(i) <<  std::endl;
//             }
//         }
//     }

//     // for(int i=0;i<N;i++){
//     //     if (node_is_boundary[i]  == 1){
//     //         init_guess(i)=b(i);
//     //     } else {
//     //         init_guess(i)=1.0;
//     //     }
//     // }

//     // display_vector(init_guess, N);

//     Integer  max_iter = N;
//     VectorWrapper P(a);

//     double o1;
//     P.apply2(o1);
//     // matrix_free_CG(S, P, b, x_0, N, result, num_iter);

//     // BiCGSTAB(S, P, b, init_guess, N, max_iter, o1, result, num_iter);
//     BiCGSTAB_new(S, P, b, init_guess, N, max_iter, o1, result, num_iter);

// }

} // namespace mars

#endif // MARS_ST_CG_HPP