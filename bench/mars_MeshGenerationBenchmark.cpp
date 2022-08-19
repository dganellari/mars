#include "mars.hpp"

#include <benchmark/benchmark.h>
#include <array>

using namespace mars;

// FIXME (depends on serial)
//  template <Integer Dim>
//  static void BM_MeshGeneration(benchmark::State& state) {
//      std::array<Integer, 4> dims{0, 0, 0, 0};

//     constexpr Integer n_x_dim = 100;
//     for (Integer d = 0; d < Dim; ++d) {
//         dims[d] = n_x_dim;
//     }

//     for (auto _ : state) {
//         Mesh<Dim> mesh;
//         generate_cube(mesh, dims[0], dims[1], dims[2]);

//         benchmark::DoNotOptimize(mesh);
//     }
// }

// BENCHMARK_TEMPLATE(BM_MeshGeneration, 1);
// BENCHMARK_TEMPLATE(BM_MeshGeneration, 2);
// BENCHMARK_TEMPLATE(BM_MeshGeneration, 3);

// FIXME
// BENCHMARK_TEMPLATE(BM_MeshGeneration, 4);