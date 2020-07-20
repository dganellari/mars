#include "mars.hpp"

#include <benchmark/benchmark.h>

using namespace mars;

static void BM_MeshGeneration4(benchmark::State& state) {
    for (auto _ : state) {
        Mesh4 mesh4;
    }
}

BENCHMARK(BM_MeshGeneration4);