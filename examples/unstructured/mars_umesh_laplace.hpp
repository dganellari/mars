#ifndef MARS_UMESH_LAPLACE_HPP
#define MARS_UMESH_LAPLACE_HPP

#include <memory>

#include "mars_base.hpp"
#include "mars_fe_values.hpp"
#include "mars_globals.hpp"
#include "mars_identity_operator.hpp"
#include "mars_simplex_laplacian.hpp"
#include "mars_tensor_laplacian.hpp"
#include "mars_umesh_operator.hpp"

namespace mars {

    template <class Mesh>
    class UMeshLaplace final : public UMeshOperator<Mesh> {
    public:
        using Elem = typename Mesh::Elem;
        using SideElem = typename Mesh::SideElem;
        using Super = mars::UMeshOperator<Mesh>;
        using Space = ViewVectorType<Real>::execution_space;

        static constexpr int Dim = Mesh::Dim;
        static constexpr int NFuns = Elem::NNodes;

        UMeshLaplace(Mesh &mesh) : Super(mesh) {}

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
            auto active = mesh.get_view_active();

            const Integer n_nodes = mesh.n_nodes();

            Kokkos::parallel_for(
                n_nodes, MARS_LAMBDA(const Integer i) { op_x(i) = 0.0; });

            // using PolicyType = TeamPolicy<Space>;

            // parallel_for(
            //     policy_type(N, TEAM_SIZE).set_scratch_size(PerTeam(0, 4096)),
            //     KOKKOS_LAMBDA(const typename policy_type::member_type &team_handle) {
            //       int ts = team_handle.team_size();     // TEAM_SIZE
            //       int tid = team_handle.team_rank();    // between 0 and TEAM_SIZE
            //       int ls = team_handle.league_size();   // returns N
            //       int lid = team_handle.league_rank();  // between 0 and N

            //       int value = tid * 5;
            //       team_handle.team_broadcast(value, 3);
            //       // value==15 on every thread
            //       value += tid;
            //       team_handle.team_broadcast([&](int &var) { var *= 2 }, value, 2);
            //       // value==34 on every thread
            //       int global;
            //       int scan = team_handle.team_scan(tid + 1, &global);
            //       // scan == tid*(tid+1)/2 on every thread
            //       // global == ts*(ts-1)/2 on every thread
            //       Kokkos::View<int *,
            //       policy_type::execution_space::scratch_memory_type>
            //           a(team_handle.team_scratch(0), 1024);
            //     });

            Laplacian<Elem> lapl;

            Kokkos::parallel_for(
                "UMeshLaplace::apply", mesh.n_elements(), MARS_LAMBDA(const Integer i) {
                    Real u[NFuns];
                    Real Au[NFuns];
                    Integer idx[NFuns];
                    Real J_inv_e[Dim * Dim];

                    if (!active(i)) return;  // ACTIVE

                    for (int k = 0; k < (Dim * Dim); ++k) {
                        J_inv_e[k] = J_inv(i, k);
                    }

                    for (Integer k = 0; k < NFuns; ++k) {
                        idx[k] = elems(i, k);
                        u[k] = x(idx[k]);
                    }

                    lapl.one_thread_eval(J_inv_e, det_J(i), u, Au);

                    for (Integer k = 0; k < NFuns; ++k) {
                        Kokkos::atomic_add(&op_x(idx[k]), Au[k]);
                    }
                });

            if (this->identity()) {
                this->identity()->apply(x, op_x);
            }
        }

        class JacobiPreconditioner final : public UMeshJacobiPreconditioner<Mesh> {
        public:
            void init(FEValues<Mesh> &values) override {
                auto mesh = values.mesh();
                ViewMatrixType<Integer> elems = values.mesh().get_view_elements();
                auto det_J = values.det_J();
                auto J_inv = values.J_inv();

                ViewVectorType<Real> inv_diag("inv_diag", mesh.n_nodes());
                auto active = mesh.get_view_active();

                Laplacian<Elem> lapl;

                Kokkos::parallel_for(
                    "JacobiPreconditioner::init", mesh.n_elements(), MARS_LAMBDA(const Integer i) {
                        Integer idx[NFuns];
                        Real val[NFuns];
                        Real J_inv_e[Dim * Dim];

                        if (!active(i)) return;  // ACTIVE

                        const Real det_J_e = det_J(i);

                        assert(det_J_e > 0.0);

                        for (Integer k = 0; k < (Dim * Dim); ++k) {
                            J_inv_e[k] = J_inv(i, k);
                        }

                        for (Integer k = 0; k < NFuns; ++k) {
                            idx[k] = elems(i, k);
                        }

                        lapl.one_thread_eval_diag(J_inv_e, det_J_e, val);

                        for (Integer k = 0; k < NFuns; ++k) {
                            assert(val[k] != 0.0);
                            Real inv_val = val[k];

                            assert(inv_val == inv_val);

                            Kokkos::atomic_add(&inv_diag(idx[k]), inv_val);
                        }
                    });

                Kokkos::parallel_for(
                    mesh.n_nodes(), MARS_LAMBDA(const Integer d) { inv_diag(d) = 1. / inv_diag(d); });

                this->inv_diag_ = inv_diag;
            }
        };
    };

}  // namespace mars

#endif  // MARS_UMESH_LAPLACE_HPP
