#include <algorithm>
#include <cassert>
#include <cstdlib>
#include <err.h>
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
#endif // WITH_KOKKOS

#include "mars_matrix_free_operator.hpp"
#include "mars_poisson.hpp"

namespace mars {

template <class Mesh> class FEValues {
public:
  static constexpr int NFuns = Mesh::Dim + 1;

  template <class Quadrature>
  FEValues(const Mesh &mesh, const Quadrature &q)
      : mesh_(mesh), q(q),
        grad(mesh.n_elements() * q.n_points() * NFuns, Mesh::Dim),
        J_inv(mesh.n_elements()) {}

  Mesh mesh;
  Quadrature q;
  ViewMatrixType<Real> grad;
  ViewMatrixType<Real> J_inv;
};

template <class Mesh> class UMeshLaplace final {
public:
  using Elem = typename Mesh::Elem;
  using SideElem = typename Mesh::SideElem;

  UMeshLaplace(Mesh &mesh) : mesh_(mesh) {}

  void apply(const ViewVectorType<Real> &x, ViewVectorType<Real> &op_x) {
    const Integer n_nodes = mesh_.n_nodes();

    // Avoid capturing this
    Mesh mesh = mesh_;

    Kokkos::parallel_for(
        n_nodes, MARS_LAMBDA(const Integer i) { op_x(i) = 0.0; });

    Kokkos::parallel_for(
        mesh.n_elements(), MARS_LAMBDA(const Integer i) {
          Elem e = mesh.elem(i);
          for (Integer i = 0; i < e.n_nodes(); ++i) {
            Kokkos::atomic_add(&op_x(e.nodes[i]), 1.0);
          }
        });
  }

  Mesh mesh_;
};
} // namespace mars

int main(int argc, char *argv[]) {
  using namespace mars;
  Kokkos::initialize(argc, argv);

#ifdef MARS_USE_CUDA
  cudaDeviceSetLimit(cudaLimitStackSize,
                     32768); // set stack to 32KB only for cuda since it is
                             // not yet supported in kokkos.
#endif

  {
    using PMesh = ParallelMesh2;
    using Elem = typename PMesh::Elem;
    using SideElem = typename PMesh::SideElem;

    using SMesh = Mesh2;

    Integer nx = 20, ny = 20, nz = 0;
    PMesh mesh;
    generate_cube(mesh, nx, ny, nz);

    ParallelBisection<PMesh> bisection(&mesh);

    ViewVectorType<Integer> marked("marked", 1);

    bisection.refine(marked);

    const Integer n_nodes = mesh.n_nodes();

    ViewVectorType<Real> x("X", n_nodes);
    ViewVectorType<Real> Ax("Ax", n_nodes);

    Kokkos::parallel_for(
        n_nodes, MARS_LAMBDA(const Integer i) { x(i) = 0.0; });

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

    std::cout << "n_active_elements: " << serial_mesh.n_active_elements()
              << std::endl;
    std::cout << "n_nodes:           " << serial_mesh.n_nodes() << std::endl;

    VTKMeshWriter<SMesh> w;
    w.write("mesh.vtu", serial_mesh);
  }

  Kokkos::finalize();
}
