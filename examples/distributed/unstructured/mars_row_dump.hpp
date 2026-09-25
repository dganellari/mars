#pragma once

// Debug aid for checking multi-rank assembly row by row. With MARS_ROW_DUMP=<prefix> set, each
// rank writes its owned rows as "sfc_key diag abs_row_sum rhs" to <prefix>.rank<r>.txt. The SFC key
// identifies the node independently of the rank count, so tests/release/compare_rows.py can check
// that every node is owned exactly once and that its row matches a 1-rank run. Global norms cannot
// show that: a missing contribution in one row can hide behind an extra one in another.

#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <string>
#include <vector>
#include <cuda_runtime.h>

namespace mars
{

template<class KeyType, class RealType>
void dumpOwnedRowsIfRequested(int rank, size_t nodeCount, int numOwnedDofs, int nnz, const int* d_nodeToDof,
                              const uint8_t* d_ownership, const KeyType* d_nodeSfcKey, const int* d_rowPtr,
                              const RealType* d_values, const int* d_diagPtr, const RealType* d_rhs)
{
    const char* prefix = std::getenv("MARS_ROW_DUMP");
    if (!prefix || !*prefix) return;

    std::vector<int> nodeToDof(nodeCount), rowPtr(numOwnedDofs + 1), diagPtr(numOwnedDofs);
    std::vector<uint8_t> own(nodeCount);
    std::vector<KeyType> key(nodeCount);
    std::vector<RealType> values(nnz), rhs(numOwnedDofs);
    cudaDeviceSynchronize();
    cudaMemcpy(nodeToDof.data(), d_nodeToDof, nodeCount * sizeof(int), cudaMemcpyDeviceToHost);
    cudaMemcpy(own.data(), d_ownership, nodeCount * sizeof(uint8_t), cudaMemcpyDeviceToHost);
    cudaMemcpy(key.data(), d_nodeSfcKey, nodeCount * sizeof(KeyType), cudaMemcpyDeviceToHost);
    cudaMemcpy(rowPtr.data(), d_rowPtr, (numOwnedDofs + 1) * sizeof(int), cudaMemcpyDeviceToHost);
    cudaMemcpy(diagPtr.data(), d_diagPtr, numOwnedDofs * sizeof(int), cudaMemcpyDeviceToHost);
    cudaMemcpy(values.data(), d_values, nnz * sizeof(RealType), cudaMemcpyDeviceToHost);
    cudaMemcpy(rhs.data(), d_rhs, numOwnedDofs * sizeof(RealType), cudaMemcpyDeviceToHost);

    std::string path = std::string(prefix) + ".rank" + std::to_string(rank) + ".txt";
    FILE* f = std::fopen(path.c_str(), "w");
    if (!f)
    {
        std::fprintf(stderr, "MARS_ROW_DUMP: cannot open %s\n", path.c_str());
        return;
    }
    for (size_t i = 0; i < nodeCount; ++i)
    {
        int dof = nodeToDof[i];
        if (own[i] == 0 || dof < 0 || dof >= numOwnedDofs) continue;
        double absSum = 0.0;
        for (int j = rowPtr[dof]; j < rowPtr[dof + 1]; ++j)
            absSum += values[j] < 0 ? -double(values[j]) : double(values[j]);
        double diag = diagPtr[dof] >= 0 ? double(values[diagPtr[dof]]) : 0.0;
        std::fprintf(f, "%llu %.17g %.17g %.17g\n", static_cast<unsigned long long>(key[i]), diag, absSum,
                     double(rhs[dof]));
    }
    std::fclose(f);
}

} // namespace mars
