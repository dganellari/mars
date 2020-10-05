#include <benchmark/benchmark.h>
#include "mars_instance.hpp"

// BENCHMARK_MAIN();

int main(int argc, char** argv) {
    using namespace mars;

    MARS::init(argc, argv);

    ::benchmark::Initialize(&argc, argv);
    ::benchmark::RunSpecifiedBenchmarks();

    return MARS::finalize();
}
