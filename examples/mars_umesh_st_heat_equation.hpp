#ifndef MARS_UMESH_ST_HEAT_EQUATION_HPP
#define MARS_UMESH_ST_HEAT_EQUATION_HPP

#include <memory>

#include "mars_base.hpp"
#include "mars_fe_values.hpp"
#include "mars_globals.hpp"
#include "mars_identity_operator.hpp"
#include "mars_simplex_laplacian.hpp"
#include "mars_umesh_operator.hpp"

namespace mars {

    template <class Mesh>
    class SpaceTimeMixed {
    public:
        using Elem = typename Mesh::Elem;
        using Point = typename Mesh::Point;
        static constexpr int Dim = Mesh::Dim;
        static constexpr int NFuns = Mesh::Dim + 1;

        MARS_INLINE_FUNCTION static void one_thread_eval_diag_add(const Real *J_inv, const Real &det_J, Real *val) {
            // Real g_ref[Dim], g[Dim];

            // for (int d = 0; d < Dim; ++d) {
            //     g_ref[d] = -1;
            // }

            // m_t_v_mult(J_inv, g_ref, g);
            // val[0] = dot(g, g) * det_J;

            // for (int d = 0; d < Dim; ++d) {
            //     g_ref[d] = 0;
            // }

            // for (int d = 0; d < Dim; ++d) {
            //     g_ref[d] = 1;
            //     m_t_v_mult(J_inv, g_ref, g);

            //     val[d + 1] = dot(g, g) * det_J;

            //     g_ref[d] = 0;
            // }
        }

        MARS_INLINE_FUNCTION static void one_thread_eval_add(const Real *J_inv,
                                                             const Real &det_J,
                                                             const Real *u,
                                                             Real *val) {
            Real g_ref[Dim], g[Dim], g_fe[Dim];

            ///////////////////////////////////////////////////////////////////
            ////////////////// Gradient with local basis function /////////////

            for (int d = 0; d < Dim; ++d) {
                g_ref[d] = -1 * u[0];
            }

            for (int i = 1; i < NFuns; ++i) {
                g_ref[i - 1] += u[i];
            }

            ///////////////////////////////////////////////////////////////////
            ////////////////// Transform gradient to physical coordinates //////

            Algebra<Dim>::m_t_v_mult(J_inv, g_ref, g);

            Real ut = g[Dim - 1];

            for (int i = 0; i < NFuns; ++i) {
                val[i] += ut * det_J * 1. / NFuns;
            }
        }
    };

    template <class Mesh>
    class UMeshSTHeatEquation final : public UMeshOperator<Mesh> {
    public:
        using Elem = typename Mesh::Elem;
        using SideElem = typename Mesh::SideElem;
        using Super = mars::UMeshOperator<Mesh>;

        static constexpr int Dim = Mesh::Dim;
        static constexpr int NFuns = Mesh::Dim + 1;

        UMeshSTHeatEquation(Mesh &mesh) : Super(mesh) {}

        void init() override {
            this->values().init();
            this->set_precontitioner(std::make_shared<JacobiPreconditioner>());
            this->preconditioner()->init(this->values());
        }

        void apply(const ViewVectorType<Real> &x, ViewVectorType<Real> &op_x) override {
            // For Kokkos-Cuda
            auto values = this->values();
            auto det_J = values.det_J();
            auto J_inv = values.J_inv();
            auto mesh = values.mesh();

            ViewMatrixType<Integer> elems = mesh.get_view_elements();

            const Integer n_nodes = mesh.n_nodes();

            Kokkos::parallel_for(
                n_nodes, MARS_LAMBDA(const Integer i) { op_x(i) = 0.0; });

            Kokkos::parallel_for(
                "UMeshSTHeatEquation::apply", mesh.n_elements(), MARS_LAMBDA(const Integer i) {
                    Real u[NFuns];
                    Real Au[NFuns];
                    Integer idx[NFuns];
                    Real J_inv_e[Dim * Dim];

                    for (int k = 0; k < (Dim * Dim); ++k) {
                        J_inv_e[k] = J_inv(i, k);
                    }

                    for (Integer k = 0; k < NFuns; ++k) {
                        idx[k] = elems(i, k);
                        u[k] = x(idx[k]);
                    }

                    SimplexLaplacian<Mesh>::one_thread_eval(J_inv_e, det_J(i), u, Au);
                    SpaceTimeMixed<Mesh>::one_thread_eval_add(J_inv_e, det_J(i), u, Au);

                    for (Integer k = 0; k < NFuns; ++k) {
                        Kokkos::atomic_add(&op_x(idx[k]), Au[k]);
                    }
                });

            this->identity()->apply(x, op_x);
        }

        class JacobiPreconditioner final : public UMeshJacobiPreconditioner<Mesh> {
        public:
            void init(FEValues<Mesh> &values) override {
                auto mesh = values.mesh();
                ViewMatrixType<Integer> elems = values.mesh().get_view_elements();
                auto det_J = values.det_J();
                auto J_inv = values.J_inv();

                ViewVectorType<Real> inv_diag("inv_diag", mesh.n_nodes());

                Kokkos::parallel_for(
                    "JacobiPreconditioner::init", mesh.n_elements(), MARS_LAMBDA(const Integer i) {
                        Integer idx[NFuns];
                        Real val[NFuns];
                        Real J_inv_e[Dim * Dim];
                        const Real det_J_e = det_J(i);

                        assert(det_J_e > 0.0);

                        for (Integer k = 0; k < (Dim * Dim); ++k) {
                            J_inv_e[k] = J_inv(i, k);
                        }

                        for (Integer k = 0; k < NFuns; ++k) {
                            idx[k] = elems(i, k);
                        }

                        SimplexLaplacian<Mesh>::one_thread_eval_diag(J_inv_e, det_J_e, val);
                        SpaceTimeMixed<Mesh>::one_thread_eval_diag_add(J_inv_e, det_J_e, val);

                        for (Integer k = 0; k < NFuns; ++k) {
                            assert(val[k] != 0.0);
                            Real inv_val = 1. / val[k];

                            assert(inv_val == inv_val);

                            Kokkos::atomic_add(&inv_diag(idx[k]), inv_val);
                        }
                    });

                this->inv_diag_ = inv_diag;
            }
        };
    };

}  // namespace mars

#endif  // MARS_UMESH_ST_HEAT_EQUATION_HPP
