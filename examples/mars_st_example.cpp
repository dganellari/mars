#include <err.h>
#include <algorithm>
#include <cassert>
#include <cstdlib>
#include <iostream>
#include <numeric>

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

namespace mars {

    template <class Mesh>
    class FEValues {
    public:
        using Elem = typename Mesh::Elem;
        using Point = typename Mesh::Point;
        static constexpr int Dim = Mesh::Dim;
        static constexpr int NFuns = Mesh::Dim + 1;
        static constexpr int NQPoints = 1;
        // template <class Quadrature>
        FEValues(const Mesh &mesh)
            : mesh_(mesh), det_J("det_J", mesh.n_elements()), J_inv("J_inv", mesh.n_elements(), Dim * Dim) {}

        void init() {
            Mesh mesh = mesh_;

            Kokkos::parallel_for(
                mesh.n_elements(), MARS_LAMBDA(const Integer i) {
                    Elem e = mesh.elem(i);
                    Point p0 = mesh.point(e.nodes[0]);
                    Point pk;

                    Real J[Dim * Dim];

                    for (int k = 1; k < e.n_nodes(); ++k) {
                        const int km1 = k - 1;
                        pk = mesh.point(e.nodes[k]);

                        for (int d = 0; d < Dim; ++d) {
                            J[d * Dim + km1] = pk[d] - p0[d];
                        }
                    }

                    Real e_det_J = det(J);
                    invert(J, &J_inv(i, 0), e_det_J);
                    det_J(i) = e_det_J;
                });
        }

        MARS_INLINE_FUNCTION Real det1(const Real *m) { return m[0]; }

        MARS_INLINE_FUNCTION Real det2(const Real *m) { return m[0] * m[3] - m[2] * m[1]; }

        MARS_INLINE_FUNCTION Real det3(const Real *m) {
            return m[0] * m[4] * m[8] + m[1] * m[5] * m[6] + m[2] * m[3] * m[4] - m[0] * m[5] * m[4] -
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

    private:
        Mesh mesh_;

    public:
        ViewMatrixType<Real> J_inv;
        ViewVectorType<Real> det_J;
    };

    template <class Mesh>
    class SimplexLaplacian {
    public:
        using Elem = typename Mesh::Elem;
        using Point = typename Mesh::Point;
        static constexpr int Dim = Mesh::Dim;
        static constexpr int NFuns = Mesh::Dim + 1;

        MARS_INLINE_FUNCTION static void one_thread_eval(const Real *J_inv,
                                                         const Real &det_J,
                                                         const Real *u,
                                                         Real *val) {
            Real g_ref[Dim], g[Dim];

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

            for (int d1 = 0; d1 < Dim; ++d1) {
                g[d1] = 0;

                for (int d2 = 0; d2 < Dim; ++d2) {
                    g[d1] += J_inv[d1 + d2 * Dim] * g_ref[d2];
                }
            }

            ///////////////////////////////////////////////////////////////////
            ////////////////// evaluate bilinear form ////////////////////////

            ///////////// Evaluate for Phi_0(x) = 1 - x_0 - ... x_n
            val[0] = 0.0;
            for (int d1 = 0; d1 < Dim; ++d1) {
                Real temp = 0.0;

                for (int d2 = 0; d2 < Dim; ++d2) {
                    temp += J_inv[d1 + d2 * Dim];
                }

                val[0] += g[d1] * -temp;
            }

            ///////////// Evaluate for Phi_i(x) = x_{i-1}
            for (int i = 1; i < NFuns; ++i) {
                int d = i - 1;
                val[i] = 0.0;

                for (int d1 = 0; d1 < Dim; ++d1) {
                    val[i] += J_inv[d * Dim + d1] * g[d1];
                }

                val[i] *= det_J;
            }
        }
    };

    template <class Mesh>
    class UMeshLaplace final {
    public:
        using Elem = typename Mesh::Elem;
        using SideElem = typename Mesh::SideElem;

        static constexpr int Dim = Mesh::Dim;
        static constexpr int NFuns = Mesh::Dim + 1;

        UMeshLaplace(Mesh &mesh) : values_(mesh) { values_.init(); }

        void apply(const ViewVectorType<Real> &x, ViewVectorType<Real> &op_x) {
            // For Kokkos-Cuda
            FEValues<Mesh> values = values_;

            Mesh mesh = values.mesh();
            const Integer n_nodes = mesh.n_nodes();

            Kokkos::parallel_for(
                n_nodes, MARS_LAMBDA(const Integer i) { op_x(i) = 0.0; });

            Kokkos::parallel_for(
                mesh.n_elements(), MARS_LAMBDA(const Integer i) {
                    Real u[NFuns];
                    Real Au[NFuns];

                    Elem e = mesh.elem(i);
                    for (Integer i = 0; i < e.n_nodes(); ++i) {
                        u[i] = x(e.nodes[i]);
                    }

                    SimplexLaplacian<Mesh>::one_thread_eval(&values.J_inv(i, 0), values.det_J(i), u, Au);

                    for (Integer i = 0; i < e.n_nodes(); ++i) {
                        Kokkos::atomic_add(&op_x(e.nodes[i]), Au[i]);
                    }
                });
        }

        FEValues<Mesh> values_;
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
        using PMesh = ParallelMesh2;
        using Elem = typename PMesh::Elem;
        using SideElem = typename PMesh::SideElem;

        using SMesh = Mesh2;

        Integer nx = 6, ny = 6, nz = 0;
        PMesh mesh;
        generate_cube(mesh, nx, ny, nz);

        ParallelBisection<PMesh> bisection(&mesh);

        Integer n_marked = 3;
        ViewVectorType<Integer> marked("marked", n_marked);

        Kokkos::parallel_for(
            n_marked, MARS_LAMBDA(const Integer i) { marked(i) = i; });

        bisection.refine(marked);

        FEValues<PMesh> values(mesh);
        values.init();

        // Kokkos::parallel_for(
        //     mesh.n_elements(), MARS_LAMBDA(const Integer i) { printf("%g\n", values.det_J[i]); });

        const Integer n_nodes = mesh.n_nodes();

        ViewVectorType<Real> x("X", n_nodes);
        ViewVectorType<Real> Ax("Ax", n_nodes);

        Kokkos::parallel_for(
            n_nodes, MARS_LAMBDA(const Integer i) { x(i) = 1.0; });

        UMeshLaplace<PMesh> op(mesh);
        op.apply(x, Ax);

        ViewVectorType<Real>::HostMirror Ax_host("Ax_host", n_nodes);
        Kokkos::deep_copy(Ax_host, Ax);

        for (Integer i = 0; i < n_nodes; ++i) {
            std::cout << Ax_host(i) << std::endl;
        }

        /////////////////////////////////////////////////////////////////////////////

        SMesh serial_mesh;
        convert_parallel_mesh_to_serial(serial_mesh, mesh);

        std::cout << "n_active_elements: " << serial_mesh.n_active_elements() << std::endl;
        std::cout << "n_nodes:           " << serial_mesh.n_nodes() << std::endl;

        VTKMeshWriter<SMesh> w;
        w.write("mesh.vtu", serial_mesh);
    }

    Kokkos::finalize();
}
