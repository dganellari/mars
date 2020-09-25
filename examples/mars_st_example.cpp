#include <err.h>
#include <algorithm>
#include <cassert>
#include <cstdlib>
#include <iostream>
#include <numeric>

#include <KokkosBlas1_nrm1.hpp>
#include <KokkosBlas1_nrminf.hpp>

#include "mars_benchmark.hpp"
#include "mars_bisection.hpp"
#include "mars_lagrange_element.hpp"
#include "mars_lepp_benchmark.hpp"
#include "mars_mesh.hpp"
#include "mars_mesh_partition.hpp"
#include "mars_partitioned_bisection.hpp"
#include "mars_prelepp_benchmark.hpp"
#include "mars_quality.hpp"
#include "mars_simplex.hpp"
#include "mars_utils.hpp"
#include "mars_vtk_writer.hpp"

#include "mars_longest_edge.hpp"
#include "mars_memory.hpp"
#include "mars_mesh_reader.hpp"
#include "mars_mesh_writer.hpp"
#include "mars_oldest_edge.hpp"
#include "mars_ranked_edge.hpp"
#include "mars_test.hpp"

#include "mars_mesh_generation.hpp"

#ifdef WITH_KOKKOS
#include "mars_lepp_benchmark_kokkos.hpp"
#include "mars_test_kokkos.hpp"
#endif  // WITH_KOKKOS

#include "mars_matrix_free_operator.hpp"
#include "mars_poisson.hpp"

#include "mars_laplace_ex.hpp"
#include "mars_precon_conjugate_grad.hpp"
#include "vtu_writer.hpp"

namespace mars {

    template <class Mesh>
    class FEValues {
    public:
        using Elem = typename Mesh::Elem;
        using Point = typename Mesh::Point;
        static const int Dim = Mesh::Dim;
        static const int ManifoldDim = Mesh::ManifoldDim;
        static const int NFuns = Elem::ElemType;
        static const int NQPoints = 1;
        // template <class Quadrature>
        FEValues(const Mesh &mesh) : mesh_(mesh) {}

        void init() {
            ViewMatrixType<Integer> elems = mesh_.get_view_elements();
            ViewMatrixType<Real> points = mesh_.get_view_points();

            auto ne = elems.extent(0);
            auto nen = elems.extent(1);

            det_J_ = ViewVectorType<Real>("det_J", mesh_.n_elements());
            J_inv_ = ViewMatrixType<Real>("J_inv", mesh_.n_elements(), Dim * Dim);

            auto det_J = det_J_;
            auto J_inv = J_inv_;

            const Integer n_nodes = mesh_.n_nodes();

            Kokkos::parallel_for(
                "FEValues::init", mesh_.n_elements(), MARS_LAMBDA(const Integer i) {
                    Integer idx[NFuns];
                    Real J[Dim * Dim];
                    Real p0[Dim], pk[Dim];

                    for (int k = 0; k < NFuns; ++k) {
                        idx[k] = elems(i, k);
                        assert(idx[k] < n_nodes);
                    }

                    for (int d = 0; d < Dim; ++d) {
                        p0[d] = points(idx[0], d);
                    }

                    for (int k = 1; k < NFuns; ++k) {
                        const int km1 = k - 1;

                        for (int d = 0; d < Dim; ++d) {
                            pk[d] = points(idx[k], d);
                        }

                        for (int d = 0; d < Dim; ++d) {
                            J[d * Dim + km1] = pk[d] - p0[d];
                        }
                    }

                    Real e_det_J = det(J);
                    invert(J, &J_inv(i, 0), e_det_J);
                    det_J(i) = Kokkos::ArithTraits<Real>::abs(e_det_J);

                    assert(e_det_J == e_det_J);
                    assert(e_det_J != 0.0);
                    assert(e_det_J > 0.0);

                    if (e_det_J < 0) {
                        static const int dim = Dim;
                        printf("found element with wrong orientation\n");
                        Integer temp = elems(i, dim);
                        elems(i, dim) = elems(i, dim - 1);
                        elems(i, dim - 1) = temp;
                    }
                });

            //     Kokkos::parallel_for(
            //         "print_det_J", mesh_.n_elements(), MARS_LAMBDA(const Integer i) { printf("%g\n", det_J(i)); });

            Real measure = KokkosBlas::nrm1(det_J_);

            std::cout << "measure: " << measure << std::endl;
        }

        MARS_INLINE_FUNCTION Real det1(const Real *m) { return m[0]; }

        MARS_INLINE_FUNCTION Real det2(const Real *m) { return m[0] * m[3] - m[2] * m[1]; }

        MARS_INLINE_FUNCTION Real det3(const Real *m) {
            return m[0] * m[4] * m[8] + m[1] * m[5] * m[6] + m[2] * m[3] * m[7] - m[0] * m[5] * m[7] -
                   m[1] * m[3] * m[8] - m[2] * m[4] * m[6];
        }

        MARS_INLINE_FUNCTION Real det(const Real *m) {
            switch (Dim) {
                case 1:
                    return det1(m);
                case 2:
                    return det2(m);
                case 3:
                    return det3(m);
                default:
                    return 0.0;
            }
        }

        MARS_INLINE_FUNCTION static bool invert(const Real *mat, Real *mat_inv, const Real &det) {
            switch (Dim) {
                case 2: {
                    return invert2(mat, mat_inv, det);
                }

                case 3: {
                    return invert3(mat, mat_inv, det);
                }

                default:
                    return false;
            }
        }

        MARS_INLINE_FUNCTION static bool invert2(const Real *mat, Real *mat_inv, const Real &det) {
            mat_inv[0] = mat[3] / det;
            mat_inv[1] = -mat[1] / det;
            mat_inv[2] = -mat[2] / det;
            mat_inv[3] = mat[0] / det;
            return true;
        }

        MARS_INLINE_FUNCTION static bool invert3(const Real *mat, Real *mat_inv, const Real &det) {
            assert(det != 0.);

            if (det == 0.) {
                return false;
            }

            mat_inv[0] = (mat[4] * mat[8] - mat[5] * mat[7]) / det;
            mat_inv[1] = (mat[2] * mat[7] - mat[1] * mat[8]) / det;
            mat_inv[2] = (mat[1] * mat[5] - mat[2] * mat[4]) / det;
            mat_inv[3] = (mat[5] * mat[6] - mat[3] * mat[8]) / det;
            mat_inv[4] = (mat[0] * mat[8] - mat[2] * mat[6]) / det;
            mat_inv[5] = (mat[2] * mat[3] - mat[0] * mat[5]) / det;
            mat_inv[6] = (mat[3] * mat[7] - mat[4] * mat[6]) / det;
            mat_inv[7] = (mat[1] * mat[6] - mat[0] * mat[7]) / det;
            mat_inv[8] = (mat[0] * mat[4] - mat[1] * mat[3]) / det;
            return true;
        }

        MARS_INLINE_FUNCTION Mesh &mesh() { return mesh_; }
        MARS_INLINE_FUNCTION const Mesh &mesh() const { return mesh_; }

        // MARS_INLINE_FUNCTION Real det_J(const Integer &i) const { return det_J_(i); }
        // MARS_INLINE_FUNCTION Real *J_inv_ptr(const Integer &i) const { return &J_inv_(i, 0); }

        MARS_INLINE_FUNCTION ViewVectorType<Real> det_J() const { return det_J_; }
        MARS_INLINE_FUNCTION ViewMatrixType<Real> J_inv() const { return J_inv_; }

    private:
        Mesh mesh_;

    public:
        ViewMatrixType<Real> J_inv_;
        ViewVectorType<Real> det_J_;
    };

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

            m_t_v_mult(J_inv, g_ref, g);

            Real ut = g[Dim - 1];

            for (int i = 0; i < NFuns; ++i) {
                val[i] += ut * det_J * 1. / NFuns;
            }
        }
    };

    template <class Mesh>
    class SimplexLaplacian {
    public:
        using Elem = typename Mesh::Elem;
        using Point = typename Mesh::Point;
        static constexpr int Dim = Mesh::Dim;
        static constexpr int NFuns = Mesh::Dim + 1;

        MARS_INLINE_FUNCTION static void m_t_v_mult(const Real *A, const Real *x, Real *y) {
            for (int d1 = 0; d1 < Dim; ++d1) {
                y[d1] = 0;

                for (int d2 = 0; d2 < Dim; ++d2) {
                    y[d1] += A[d1 + d2 * Dim] * x[d2];
                }
            }
        }

        MARS_INLINE_FUNCTION static Real dot(const Real *l, const Real *r) {
            Real ret = 0.0;
            for (Integer i = 0; i < Dim; ++i) {
                ret += l[i] * r[i];
            }

            return ret;
        }

        MARS_INLINE_FUNCTION static void one_thread_eval_diag(const Real *J_inv, const Real &det_J, Real *val) {
            Real g_ref[Dim], g[Dim];

            for (int d = 0; d < Dim; ++d) {
                g_ref[d] = -1;
            }

            m_t_v_mult(J_inv, g_ref, g);
            val[0] = dot(g, g) * det_J;

            for (int d = 0; d < Dim; ++d) {
                g_ref[d] = 0;
            }

            for (int d = 0; d < Dim; ++d) {
                g_ref[d] = 1;
                m_t_v_mult(J_inv, g_ref, g);

                val[d + 1] = dot(g, g) * det_J;

                g_ref[d] = 0;
            }
        }

        MARS_INLINE_FUNCTION static void one_thread_eval(const Real *J_inv,
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

            m_t_v_mult(J_inv, g_ref, g);

            ///////////////////////////////////////////////////////////////////
            ////////////////// evaluate bilinear form ////////////////////////

            ///////////// Evaluate for Phi_0(x) = 1 - x_0 - ... x_n

            for (int d = 0; d < Dim; ++d) {
                g_ref[d] = -1;
            }

            m_t_v_mult(J_inv, g_ref, g_fe);

            val[0] = dot(g, g_fe) * det_J;

            for (int d = 0; d < Dim; ++d) {
                g_ref[d] = 0;
            }

            ///////////// Evaluate for Phi_i(x) = x_{i-1}
            for (int i = 1; i < NFuns; ++i) {
                int d = i - 1;

                g_ref[d] = 1.0;

                m_t_v_mult(J_inv, g_ref, g_fe);

                val[i] = dot(g, g_fe) * det_J;

                g_ref[d] = 0.0;
            }
        }
    };

    template <int Dim>
    MARS_INLINE_FUNCTION bool is_boundary_of_unit_cube(const Real *p) {
        bool ret = false;
        for (int d = 0; d < Dim; ++d) {
            if (p[d] <= 1e-14 || p[d] >= 1 - 1e-14) {
                ret = true;
                break;
            }
        }

        return ret;
    }

    template <class Mesh>
    class ZeroDirchletOnUnitCube {
    public:
        static const int Dim = Mesh::Dim;

        MARS_INLINE_FUNCTION void operator()(const Real *p, Real &val) const {
            if (is_boundary(p)) {
                val = 0.0;
            }
        }

        MARS_INLINE_FUNCTION static bool is_boundary(const Real *p) { return is_boundary_of_unit_cube<Dim>(p); }
    };

    class Example1Dirichlet {
    public:
        MARS_INLINE_FUNCTION void operator()(const Real *p, Real &val) const {
            if (is_boundary(p)) {
                val = ex1_exact(p);
            }
        }

        MARS_INLINE_FUNCTION static bool is_boundary(const Real *p) { return is_boundary_of_unit_cube<2>(p); }
    };

    class Example1RHS {
    public:
        MARS_INLINE_FUNCTION Real operator()(const Real *p) const { return ex1_laplacian(p); }
    };

    class Example2Dirichlet {
    public:
        MARS_INLINE_FUNCTION void operator()(const Real *p, Real &val) const {
            if (is_boundary(p)) {
                val = ex2_exact(p);
            }
        }

        MARS_INLINE_FUNCTION static bool is_boundary(const Real *p) { return is_boundary_of_unit_cube<2>(p); }
    };

    class Example2Analitcal {
    public:
        MARS_INLINE_FUNCTION Real operator()(const Real *p) const { return ex2_exact(p); }
    };

    class Example2RHS {
    public:
        MARS_INLINE_FUNCTION Real operator()(const Real *p) const { return ex2_laplacian(p); }
    };

    ////////////////////////////////

    class Example3Dirichlet {
    public:
        MARS_INLINE_FUNCTION void operator()(const Real *p, Real &val) const {
            if (is_boundary(p)) {
                val = ex3_exact(p);
            }
        }

        MARS_INLINE_FUNCTION static bool is_boundary(const Real *p) { return is_boundary_of_unit_cube<2>(p); }
    };

    class Example3RHS {
    public:
        MARS_INLINE_FUNCTION Real operator()(const Real *p) const { return ex3_laplacian(p); }
    };

    template <class Mesh>
    class One {
    public:
        static const int Dim = Mesh::Dim;
        MARS_INLINE_FUNCTION Real operator()(const Real *) const { return 1.0; }
    };

    template <class Mesh>
    class Norm2Squared {
    public:
        static const int Dim = Mesh::Dim;
        MARS_INLINE_FUNCTION Real operator()(const Real *p) const {
            Real ret = 0.0;

            for (int i = 0; i < Dim; ++i) {
                const Real x = p[i];
                ret += x * x;
            }

            return ret;
        }
    };

    template <class Mesh>
    class Interpolate {
    public:
        Interpolate(Mesh &mesh) : mesh_(mesh) {}

        template <typename F>
        void apply(ViewVectorType<Real> &v, F fun) {
            static const int Dim = Mesh::Dim;

            ViewMatrixType<Real> points = mesh_.get_view_points();

            Kokkos::parallel_for(
                "Interpolate::apply", mesh_.n_nodes(), MARS_LAMBDA(const Integer i) {
                    Real p[Dim];

                    for (int d = 0; d < Dim; ++d) {
                        p[d] = points(i, d);
                    }

                    v(i) = fun(p);
                });
        }

    private:
        Mesh &mesh_;
    };

    class IdentityOperator {
    public:
        template <class Mesh, class BC>
        void init(Mesh &mesh, BC bc) {
            static const int Dim = Mesh::Dim;

            auto points = mesh.points();
            ViewVectorType<bool> is_boundary("is_boundary", mesh.n_nodes());

            Kokkos::parallel_for(
                "IdentityOperator::init", mesh.n_nodes(), MARS_LAMBDA(const Integer i) {
                    Real p[Dim];

                    for (int d = 0; d < Dim; ++d) {
                        p[d] = points(i, d);
                    }
                    is_boundary(i) = bc.is_boundary(p);
                });

            is_boundary_ = is_boundary;
        }

        void apply(const ViewVectorType<Real> &input, ViewVectorType<Real> &x) {
            auto is_boundary = is_boundary_;

            Kokkos::parallel_for(
                "IdentityOperator::apply", is_boundary_.extent(0), MARS_LAMBDA(const Integer i) {
                    x(i) = x(i) * (!is_boundary(i)) + input(i) * is_boundary(i);
                });
        }

    private:
        ViewVectorType<bool> is_boundary_;
    };

    template <class Mesh>
    class BoundaryConditions {
    public:
        static const int Dim = Mesh::Dim;

        BoundaryConditions(Mesh &mesh) : mesh_(mesh) {}

        template <typename F>
        void apply(ViewVectorType<Real> &v, F fun) {
            ViewMatrixType<Real> points = mesh_.get_view_points();

            Kokkos::parallel_for(
                "BoundaryConditions::apply", mesh_.n_nodes(), MARS_LAMBDA(const Integer i) {
                    Real p[Dim];

                    for (int d = 0; d < Dim; ++d) {
                        p[d] = points(i, d);
                    }

                    fun(p, v(i));
                });
        }

    private:
        Mesh &mesh_;
    };

    template <class Mesh>
    class UMeshLaplace final {
    public:
        using Elem = typename Mesh::Elem;
        using SideElem = typename Mesh::SideElem;

        static constexpr int Dim = Mesh::Dim;
        static constexpr int NFuns = Mesh::Dim + 1;

        UMeshLaplace(Mesh &mesh) : values_(mesh) {}

        class FakeComm {
        public:
            static Real sum(const Real &v) { return v; }
        };

        FakeComm comm() { return FakeComm(); }

        void init() {
            values_.init();
            preconditioner_.init(values_);
        }

        void apply(const ViewVectorType<Real> &x, ViewVectorType<Real> &op_x) {
            // For Kokkos-Cuda
            auto det_J = values_.det_J();
            auto J_inv = values_.J_inv();
            auto mesh = values_.mesh();

            ViewMatrixType<Integer> elems = mesh.get_view_elements();

            const Integer n_nodes = mesh.n_nodes();

            Kokkos::parallel_for(
                n_nodes, MARS_LAMBDA(const Integer i) { op_x(i) = 0.0; });

            Kokkos::parallel_for(
                "UMeshLaplace::apply", mesh.n_elements(), MARS_LAMBDA(const Integer i) {
                    Real u[NFuns];
                    Real Au[NFuns];
                    Integer idx[NFuns];

                    for (Integer k = 0; k < NFuns; ++k) {
                        idx[k] = elems(i, k);
                        u[k] = x(idx[k]);
                    }

                    SimplexLaplacian<Mesh>::one_thread_eval(&J_inv(i, 0), det_J(i), u, Au);

                    for (Integer k = 0; k < NFuns; ++k) {
                        Kokkos::atomic_add(&op_x(idx[k]), Au[k]);
                    }
                });

            id_->apply(x, op_x);
        }

        template <class F>
        void assemble_rhs(ViewVectorType<Real> &rhs, F f) {
            auto det_J = values_.det_J();

            ViewMatrixType<Integer> elems = values_.mesh().get_view_elements();
            ViewMatrixType<Real> points = values_.mesh().get_view_points();

            auto mesh = values_.mesh();

            Kokkos::parallel_for(
                "UMeshLaplace::assemble_rhs", mesh.n_elements(), MARS_LAMBDA(const Integer i) {
                    Integer idx[NFuns];
                    Real p[Dim];
                    Real det_J_e = det_J(i);

                    for (Integer k = 0; k < NFuns; ++k) {
                        idx[k] = elems(i, k);
                    }

                    for (Integer k = 0; k < NFuns; ++k) {
                        for (Integer d = 0; d < Dim; ++d) {
                            p[d] = points(idx[k], d);
                        }

                        const Real val = f(p);
                        const Real scaled_val = val * det_J_e / NFuns;
                        Kokkos::atomic_add(&rhs(idx[k]), scaled_val);
                    }
                });
        }

        class JacobiPreconditioner {
        public:
            void init(FEValues<Mesh> &values) {
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

                        for (Integer k = 0; k < Dim * Dim; ++k) {
                            J_inv_e[k] = J_inv(i, k);
                        }

                        for (Integer k = 0; k < NFuns; ++k) {
                            idx[k] = elems(i, k);
                        }

                        SimplexLaplacian<Mesh>::one_thread_eval_diag(J_inv_e, det_J_e, val);

                        for (Integer k = 0; k < NFuns; ++k) {
                            assert(val[k] != 0.0);
                            Real inv_val = 1. / val[k];

                            assert(inv_val == inv_val);

                            Kokkos::atomic_add(&inv_diag(idx[k]), inv_val);
                        }
                    });

                inv_diag_ = inv_diag;
            }

            void apply(const ViewVectorType<Real> &x, ViewVectorType<Real> &op_x) {
                auto n = inv_diag_.extent(0);
                auto inv_diag = inv_diag_;

                Kokkos::parallel_for(
                    "JacobiPreconditioner::apply", n, MARS_LAMBDA(const Integer i) { op_x(i) = inv_diag(i) * x(i); });

                // Kokkos::parallel_for(
                //     "JacobiPreconditioner::copy", n, MARS_LAMBDA(const Integer i) { op_x(i) = x(i); });

                id_->apply(x, op_x);
            }

            void set_identity(std::shared_ptr<IdentityOperator> id) { id_ = id; }

        private:
            ViewVectorType<Real> inv_diag_;
            std::shared_ptr<IdentityOperator> id_;
        };

        void set_identity(const std::shared_ptr<IdentityOperator> &id) {
            id_ = id;
            preconditioner_.set_identity(id);
        }

        inline JacobiPreconditioner &preconditioner() { return preconditioner_; }

        FEValues<Mesh> values_;
        JacobiPreconditioner preconditioner_;
        std::shared_ptr<IdentityOperator> id_;
    };
}  // namespace mars

int main(int argc, char *argv[]) {
    using namespace mars;
    Kokkos::initialize(argc, argv);

#ifdef MARS_USE_CUDA
    cudaDeviceSetLimit(cudaLimitStackSize,
                       32768);  // set stack to 32KB only for cuda since it is
                                // not yet supported in kokkos.
#endif

    {
        // using PMesh = ParallelMesh2;
        // using SMesh = Mesh2;

        using PMesh = ParallelMesh3;
        using SMesh = Mesh3;

        using Elem = typename PMesh::Elem;
        using SideElem = typename PMesh::SideElem;

        Integer nx = 6, ny = 6, nz = (PMesh::Dim > 2) ? 6 : 0;
        if (argc > 1) {
            nx = atol(argv[1]);
            ny = atol(argv[1]);
            nz = atol(argv[1]);

            if (PMesh::Dim <= 2) nz = 0;
        }

        PMesh mesh;
        generate_cube(mesh, nx, ny, nz);

        // ParallelBisection<PMesh> bisection(&mesh);

        // Integer n_marked = 3;
        // ViewVectorType<Integer> marked("marked", n_marked);

        // Kokkos::parallel_for(
        //     n_marked, MARS_LAMBDA(const Integer i) { marked(i) = i; });

        // bisection.refine(marked);

        ZeroDirchletOnUnitCube<PMesh> bc_fun;
        One<PMesh> rhs_fun;
        One<PMesh> an_fun;  // FIXME

        // ZeroDirchletOnUnitCube<PMesh> bc_fun;
        // Norm2Squared<PMesh> rhs_fun;

        // Example1Dirichlet bc_fun;
        // Example1RHS rhs_fun;

        // Example2Dirichlet bc_fun;
        // Example2RHS rhs_fun;
        // Example2Analitcal an_fun;

        const Integer n_nodes = mesh.n_nodes();

        ViewVectorType<Real> x("X", n_nodes);
        ViewVectorType<Real> rhs("rhs", n_nodes);
        ViewVectorType<Real> Ax("Ax", n_nodes);

        UMeshLaplace<PMesh> op(mesh);
        BoundaryConditions<PMesh> bc(mesh);
        op.init();

        auto id = std::make_shared<IdentityOperator>();
        id->init(mesh, bc_fun);

        op.set_identity(id);
        op.assemble_rhs(rhs, rhs_fun);

        Real nrm_rhs = KokkosBlas::nrm1(rhs);

        std::cout << "nrm_rhs : " << nrm_rhs << std::endl;

        bc.apply(rhs, bc_fun);
        bc.apply(x, bc_fun);

        auto prec = op.preconditioner();

        Integer num_iter = 0;
        bcg_stab(op, prec, rhs, rhs.extent(0), x, num_iter);

        // Compute Error
        ViewVectorType<Real> x_exact("X_exact", n_nodes);
        ViewVectorType<Real> diff("Diff", n_nodes);

        Interpolate<PMesh> interp(mesh);
        interp.apply(x_exact, an_fun);

        Kokkos::deep_copy(diff, x_exact);

        KokkosBlas::axpy(-1.0, x, diff);
        Real err = KokkosBlas::nrminf(diff);
        std::cout << "err : " << err << std::endl;

        ///////////////////////////////////////////////////////////////////////////

        ViewVectorType<Real>::HostMirror x_host("x_host", n_nodes);
        ViewVectorType<Real>::HostMirror rhs_host("rhs_host", n_nodes);
        Kokkos::deep_copy(x_host, x);
        Kokkos::deep_copy(rhs_host, rhs);

        SMesh serial_mesh;
        convert_parallel_mesh_to_serial(serial_mesh, mesh);

        std::cout << "n_active_elements: " << serial_mesh.n_active_elements() << std::endl;
        std::cout << "n_nodes:           " << serial_mesh.n_nodes() << std::endl;

        VTUMeshWriter<SMesh> w;
        w.write("mesh.vtu", serial_mesh, x_host);
        w.write("mesh_rhs.vtu", serial_mesh, rhs_host);

        Kokkos::deep_copy(x_host, x_exact);
        w.write("analitic.vtu", serial_mesh, x_host);
    }

    Kokkos::finalize();
}
