#ifndef ST_OPERATOR_HPP
#define ST_OPERATOR_HPP
#include "mars_base.hpp"

#ifdef WITH_MPI
#include <mpi.h>
#endif

#include <err.h>

#ifdef WITH_KOKKOS_KERNELS
#include <KokkosBlas1_axpby.hpp>
#include <KokkosBlas1_dot.hpp>
#include <KokkosBlas1_mult.hpp>
#include <KokkosBlas1_nrm1.hpp>
#include <KokkosBlas1_nrm2.hpp>
#include <KokkosBlas1_scal.hpp>
#include <KokkosSparse_spadd.hpp>
#include <KokkosSparse_spmv.hpp>
#include <Kokkos_Core.hpp>
#include "KokkosSparse_CrsMatrix.hpp"
#include "KokkosSparse_spmv_spec.hpp"

#include <algorithm>
#include <cassert>
#include <cstdlib>
#include <iostream>
#include <map>
#include <numeric>

#include "mars_context.hpp"
#include "mars_fedm_values.hpp"
#include "mars_matrix.hpp"
#include "mars_utils.hpp"
#include "mars_vector.hpp"

namespace mars {
    // using KokkosVector = Kokkos::View<Real *>;
    // using SparseMatrix = KokkosSparse::CrsMatrix<Real, Integer, Kokkos::Serial>;

    template <class DM, class FEM>
    class Operator {
    public:
        virtual ~Operator(){};

        virtual void apply(const ViewVectorType<Real> &x, ViewVectorType<Real> &op_x) = 0;
        virtual void init() = 0;

        Operator(const context &ctx, DM &dm, FEM &fe) : ctx_(ctx), values_(dm, fe) {}

        const context &ctx() const { return ctx_; }

        const distributed_context &comm() const {
            auto commctx = ctx_->distributed;
            return *commctx;
        }

        MARS_INLINE_FUNCTION FEDMValues<DM, FEM> &values() { return values_; }
        MARS_INLINE_FUNCTION const FEDMValues<DM, FEM> &values() const { return values_; }

    private:
        const context &ctx_;
        FEDMValues<DM, FEM> values_;
    };

    // class ImplicitEulerOperator final : public Operator {
    // public:
    //   ImplicitEulerOperator(const SparseMatrix &stiffness_matrix,
    //                         const SparseMatrix &mass_matrix, const Real &dt,
    //                         const size_t &N)
    //       : stiffness_matrix(stiffness_matrix), mass_matrix(mass_matrix), dt(dt),
    //         N(N) {}

    //   void apply(const KokkosVector &a, KokkosVector &b) override {
    //     Kokkos::View<Real *> b0("b0", N);
    //     KokkosSparse::spmv("N", 1.0, stiffness_matrix, a, 0.0, b0); // A*x = b0

    //     KokkosSparse::spmv("N", 1.0, mass_matrix, a, 0.0, b); // b = M*x
    //     KokkosBlas::axpy(dt, b0, b);                          // b = b + dt*b1
    //   }

    // private:
    //   const SparseMatrix &stiffness_matrix;
    //   const SparseMatrix &mass_matrix;
    //   const Real &dt;
    //   const size_t &N;
    // };

    // class MatrixWrapper final : public Operator {
    // public:
    //   // storare la matrice, sparse matrix vector multiplication
    //   MatrixWrapper(const SparseMatrix &mat) : mat(mat) {}

    //   void apply(const KokkosVector &a, KokkosVector &b) override {
    //     KokkosSparse::spmv("N", 1.0, mat, a, 0.0, b); // mat*a = b
    //   }

    // private:
    //   const SparseMatrix &mat;
    // };

    // class VectorWrapper final : public Operator {
    // public:
    //   // storare la matrice, sparse matrix vector multiplication
    //   VectorWrapper(const KokkosVector &v) : v(v) {}

    //   void apply(const KokkosVector &a, KokkosVector &b) override {
    //     KokkosBlas::mult(0.0, b, 1.0, a, v);
    //   }

    //   void apply2(double &nrm) { nrm = KokkosBlas::nrm1(v); }

    // private:
    //   const KokkosVector &v;
    // };
}  // namespace mars

#endif  // ST_OPERATOR_HPP
#endif