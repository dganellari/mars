#pragma once

#include <iostream>
#include <iomanip>
#include <vector>
#include <algorithm>
#include <cuda_runtime.h>

namespace mars {
namespace fem {

// CUDA event-based GPU timer
class GpuTimer {
public:
    GpuTimer() {
        cudaEventCreate(&start_);
        cudaEventCreate(&stop_);
    }

    ~GpuTimer() {
        cudaEventDestroy(start_);
        cudaEventDestroy(stop_);
    }

    void start(cudaStream_t stream = 0) {
        cudaEventRecord(start_, stream);
    }

    void stop(cudaStream_t stream = 0) {
        cudaEventRecord(stop_, stream);
    }

    float elapsedMs() {
        cudaEventSynchronize(stop_);
        float ms = 0.0f;
        cudaEventElapsedTime(&ms, start_, stop_);
        return ms;
    }

private:
    cudaEvent_t start_, stop_;
};

// Performance counter structure for detailed timing and metrics
struct PerfCounters {
    // Timing (in milliseconds)
    float meshLoadTime = 0.0f;
    float sparsityBuildTime = 0.0f;
    float fieldInitTime = 0.0f;
    float assemblyWarmupTime = 0.0f;
    float assemblyTime = 0.0f;       // Total time for all iterations
    float assemblyMinTime = 1e9f;
    float assemblyMaxTime = 0.0f;
    int   assemblyIterations = 0;

    // Data sizes (in bytes)
    size_t inputFieldBytes = 0;      // Field data read per element
    size_t outputMatrixBytes = 0;    // Matrix/RHS data written
    size_t connectivityBytes = 0;    // Connectivity data read
    size_t coordinateBytes = 0;      // Coordinate data read

    // Counts
    size_t numElements = 0;
    size_t numNodes = 0;
    size_t numDofs = 0;
    size_t matrixNnz = 0;

    // Domain decomposition load balance metrics (element counts)
    double ddOwnedImbalancePct = 0.0;  // (max-min)/mean * 100 for owned elements
    double ddGhostImbalancePct = 0.0;  // (max-min)/mean * 100 for ghost elements
    double ddTotalImbalancePct = 0.0;  // (max-min)/mean * 100 for total elements
    double ddOwnedMin = 0.0, ddOwnedMax = 0.0, ddOwnedMean = 0.0;
    double ddGhostMin = 0.0, ddGhostMax = 0.0, ddGhostMean = 0.0;
    double ddTotalMin = 0.0, ddTotalMax = 0.0, ddTotalMean = 0.0;

    // Per-iteration times for detailed output
    std::vector<float> iterTimes;

    // Computed metrics
    double avgAssemblyTimeMs() const {
        return assemblyIterations > 0 ? assemblyTime / assemblyIterations : 0.0;
    }

    size_t totalBytesPerAssembly() const {
        return inputFieldBytes + outputMatrixBytes + connectivityBytes + coordinateBytes;
    }

    double bandwidthGBs() const {
        double timeS = avgAssemblyTimeMs() / 1000.0;
        return timeS > 0 ? (totalBytesPerAssembly() / 1e9) / timeS : 0.0;
    }

    double elementsPerSecond() const {
        double timeS = avgAssemblyTimeMs() / 1000.0;
        return timeS > 0 ? numElements / timeS : 0.0;
    }

    // FLOP estimate for CVFEM hex assembly (approximate)
    // Per element: 12 SCS * (shape functions + Jacobian + diffusion) ~ 2000 FLOPs
    double estimatedGFLOPs() const {
        const double flopsPerElement = 2000.0;
        double timeS = avgAssemblyTimeMs() / 1000.0;
        return timeS > 0 ? (numElements * flopsPerElement / 1e9) / timeS : 0.0;
    }

    // Record an iteration time
    void recordIteration(float timeMs) {
        iterTimes.push_back(timeMs);
        assemblyTime += timeMs;
        assemblyMinTime = std::min(assemblyMinTime, timeMs);
        assemblyMaxTime = std::max(assemblyMaxTime, timeMs);
        assemblyIterations++;
    }

    // Print detailed performance summary
    void print(int rank) const {
        if (rank != 0) return;

        std::cout << "\n";
        std::cout << "================================================================================\n";
        std::cout << "                         PERFORMANCE SUMMARY\n";
        std::cout << "================================================================================\n";
        std::cout << std::fixed << std::setprecision(3);

        std::cout << "\n--- Problem Size (rank 0, owned only) ---\n";
        std::cout << "  Elements:              " << std::setw(12) << numElements << "\n";
        std::cout << "  Nodes:                 " << std::setw(12) << numNodes << "\n";
        std::cout << "  DOFs:                  " << std::setw(12) << numDofs << "\n";
        std::cout << "  Matrix NNZ:            " << std::setw(12) << matrixNnz << "\n";
        std::cout << "  Avg NNZ/row:           " << std::setw(12) << std::setprecision(1)
                  << (numDofs > 0 ? static_cast<double>(matrixNnz) / numDofs : 0.0) << "\n";

        std::cout << "\n--- Data Transfer (per assembly) ---\n";
        std::cout << std::setprecision(2);
        std::cout << "  Input fields:          " << std::setw(12) << inputFieldBytes / 1024.0 / 1024.0 << " MB\n";
        std::cout << "  Coordinates:           " << std::setw(12) << coordinateBytes / 1024.0 / 1024.0 << " MB\n";
        std::cout << "  Connectivity:          " << std::setw(12) << connectivityBytes / 1024.0 / 1024.0 << " MB\n";
        std::cout << "  Output (matrix+RHS):   " << std::setw(12) << outputMatrixBytes / 1024.0 / 1024.0 << " MB\n";
        std::cout << "  Total:                 " << std::setw(12) << totalBytesPerAssembly() / 1024.0 / 1024.0 << " MB\n";

        std::cout << "\n--- Timing Breakdown ---\n";
        std::cout << std::setprecision(3);
        std::cout << "  Mesh loading:          " << std::setw(12) << meshLoadTime << " ms\n";
        std::cout << "  Sparsity build:        " << std::setw(12) << sparsityBuildTime << " ms\n";
        std::cout << "  Field initialization:  " << std::setw(12) << fieldInitTime << " ms\n";
        std::cout << "  Assembly warmup:       " << std::setw(12) << assemblyWarmupTime << " ms\n";

        std::cout << "\n--- Assembly Performance (" << assemblyIterations << " iterations) ---\n";
        std::cout << "  Average time:          " << std::setw(12) << avgAssemblyTimeMs() << " ms\n";
        std::cout << "  Min time:              " << std::setw(12) << assemblyMinTime << " ms\n";
        std::cout << "  Max time:              " << std::setw(12) << assemblyMaxTime << " ms\n";
        std::cout << "  Bandwidth:             " << std::setw(12) << bandwidthGBs() << " GB/s\n";
        std::cout << "  Throughput:            " << std::setw(12) << elementsPerSecond() / 1e6 << " M elem/s\n";
        std::cout << "  Est. compute:          " << std::setw(12) << estimatedGFLOPs() << " GFLOP/s\n";
        std::cout << "\n--- Domain Decomposition Imbalance ---\n";
        std::cout << std::setprecision(1);
        std::cout << "                    " << std::setw(12) << "min" << std::setw(12) << "mean" << std::setw(12) << "max" << std::setw(10) << "imbalance" << "   (ghost: spread/total_mean)\n";
        std::cout << "  Owned nodes:      "
                  << std::setw(12) << static_cast<size_t>(ddOwnedMin)
                  << std::setw(12) << static_cast<size_t>(ddOwnedMean)
                  << std::setw(12) << static_cast<size_t>(ddOwnedMax)
                  << std::setw(9)  << std::setprecision(2) << ddOwnedImbalancePct << " %\n";
        double ghostOverheadPct = ddTotalMean > 0.0 ? ddGhostMean / ddTotalMean * 100.0 : 0.0;
        std::cout << "  Ghost nodes:      "
                  << std::setw(12) << static_cast<size_t>(ddGhostMin)
                  << std::setw(12) << static_cast<size_t>(ddGhostMean)
                  << std::setw(12) << static_cast<size_t>(ddGhostMax)
                  << std::setw(9)  << std::setprecision(2) << ddGhostImbalancePct << " %"
                  << "   (" << std::setprecision(1) << ghostOverheadPct << "% of rank workload)\n";
        std::cout << "  Total nodes:      "
                  << std::setw(12) << static_cast<size_t>(ddTotalMin)
                  << std::setw(12) << static_cast<size_t>(ddTotalMean)
                  << std::setw(12) << static_cast<size_t>(ddTotalMax)
                  << std::setw(9)  << std::setprecision(2) << ddTotalImbalancePct << " %\n";

        std::cout << "================================================================================\n";
    }
};

} // namespace fem
} // namespace mars
