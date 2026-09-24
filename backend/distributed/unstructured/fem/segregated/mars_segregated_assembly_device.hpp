#pragma once
#include <cuda_runtime.h>
#include "mars_segregated_assembly.hpp"
#include "mars_segregated_geometry_device.hpp"
#include <thrust/binary_search.h>
#include <thrust/device_vector.h>
#include <thrust/iterator/counting_iterator.h>
#include <thrust/iterator/transform_iterator.h>
#include <thrust/sort.h>
#include <thrust/unique.h>
#include <algorithm>
#include <cstdint>
#include <limits>
#include <stdexcept>

namespace mars::segregated {
inline void assembly_cuda_check(cudaError_t error)
{
    if (error != cudaSuccess) throw std::runtime_error(cudaGetErrorString(error));
}
template<class T> T* device_data(thrust::device_vector<T>& vector)
{
    return thrust::raw_pointer_cast(vector.data());
}
struct AssemblyRowKey {
    __host__ __device__ std::uint64_t operator()(int row) const { return std::uint64_t(row)<<32; }
};

template<class KeyType>
__global__ void assembly_graph_keys(TetMeshView<KeyType,double> mesh, std::uint64_t* keys)
{
    const int i = blockIdx.x*blockDim.x+threadIdx.x;
    if (i < 16*mesh.element_count) {
        const int e = i/16, local = i%16;
        keys[i] = (std::uint64_t(mesh.nodes[local/4][e])<<32)|std::uint64_t(mesh.nodes[local%4][e]);
    }
    if (i < mesh.node_count) keys[16*mesh.element_count+i] = (std::uint64_t(i)<<32)|std::uint64_t(i);
}
__global__ void assembly_graph_columns(const std::uint64_t* keys, int count, int* columns)
{
    const int i = blockIdx.x*blockDim.x+threadIdx.x;
    if (i < count) columns[i] = int(keys[i]&0xffffffffu);
}

// Build once per mesh after checking native geometry/connectivity status. Only the
// unique count crosses to the host inside Thrust; connectivity and keys stay on device.
struct DeviceBlockGraph {
    int nodes = 0;
    thrust::device_vector<int> offsets, columns;
    template<class KeyType> void build(TetMeshView<KeyType,double> mesh) {
        if (mesh.node_count <= 0 || mesh.element_count <= 0
            || mesh.element_count > (std::numeric_limits<int>::max()-mesh.node_count)/16)
            throw std::runtime_error("invalid or oversized segregated graph");
        nodes = mesh.node_count;
        const int count = 16*mesh.element_count+nodes;
        thrust::device_vector<std::uint64_t> keys(count);
        const int work = std::max(16*mesh.element_count,nodes);
        assembly_graph_keys<<<(work+127)/128,128>>>(mesh,device_data(keys));
        assembly_cuda_check(cudaGetLastError());
        thrust::sort(keys.begin(),keys.end());
        const int unique = int(thrust::unique(keys.begin(),keys.end())-keys.begin());
        keys.resize(unique); columns.resize(unique); offsets.resize(nodes+1);
        auto rows = thrust::make_transform_iterator(thrust::make_counting_iterator(0),AssemblyRowKey{});
        thrust::lower_bound(keys.begin(),keys.end(),rows,rows+nodes+1,offsets.begin());
        assembly_graph_columns<<<(unique+127)/128,128>>>(device_data(keys),unique,device_data(columns));
        assembly_cuda_check(cudaGetLastError());
        // keys' lifetime ends here; its release and Thrust's host return complete construction.
    }
    template<int Components> BlockCsrView<Components> view(double* values, double* rhs) {
        return {nodes,device_data(offsets),device_data(columns),values,rhs};
    }
};

template<int Components, class KeyType>
__global__ void assemble_interior_blocks(TetMeshView<KeyType,double> mesh,
    const TetGeometry<double>* geometry, const TetInteriorInput* inputs,
    const double* velocity_gradient, const double* pressure_gradient, const double* influence,
    int begin, int end, BlockCsrView<Components> matrix, int* error)
{
    const int e = begin+blockIdx.x*blockDim.x+threadIdx.x;
    if (e >= end) return;
    int nodes[4]; for (int n = 0; n < 4; ++n) nodes[n] = int(mesh.nodes[n][e]);
    auto input = inputs[e];
    if (input.stage != (Components == 3 ? 1 : 0)) { atomicExch(error,1); return; }
    native_interior(input,geometry[e],nodes,velocity_gradient,pressure_gradient,influence);
    TetInteriorOutput output;
    tet_interior(input,output);
    if (!scatter_block(matrix,nodes,4,output.lhs,output.rhs)) atomicExch(error,1);
}

template<int Components, class KeyType>
__global__ void assemble_boundary_blocks(TetMeshView<KeyType,double> mesh,
    const TetGeometry<double>* geometry, const BoundaryAssemblyInput* inputs, int count,
    const double* pressure_gradient, const double* influence, BlockCsrView<Components> matrix, int* error)
{
    const int i = blockIdx.x*blockDim.x+threadIdx.x;
    if (i >= count) return;
    const auto& input = inputs[i];
    if (input.element < 0 || input.element >= mesh.element_count || input.face < 0 || input.face >= 4
        || input.values.stage < 0 || input.values.stage > 5 || (input.values.stage >= 3) != (Components == 3)) {
        atomicExch(error,1); return;
    }
    int nodes[4]; for (int n = 0; n < 4; ++n) nodes[n] = int(mesh.nodes[n][input.element]);
    auto values = input.values;
    if (!native_boundary(values,input,geometry[input.element],nodes,pressure_gradient,influence)) {
        atomicExch(error,1); return;
    }
    BoundaryOutput output;
    boundary_block(values,output);
    if (!scatter_block(matrix,input.nodes,values.stage == 5 ? 3 : 4,output.lhs,output.rhs)) atomicExch(error,1);
}

__global__ void assemble_momentum_nodes(const SteadyMomentumNode* inputs, int begin, int end,
    const double* volumes, const double* pressure_gradient, BlockCsrView<3> matrix, int* error)
{
    const int node = begin+blockIdx.x*blockDim.x+threadIdx.x;
    if (node >= end) return;
    auto input = inputs[node]; input.volume = volumes[node];
    for (int j = 0; j < 3; ++j) input.pressure_gradient[j] = pressure_gradient[3*node+j];
    double lhs[9],rhs[3]; steady_momentum_node(input,lhs,rhs);
    if (!scatter_block(matrix,&node,1,lhs,rhs)) atomicExch(error,1);
}
__global__ void finish_momentum_rows(BlockCsrView<3> matrix, const double* volumes,
    const double* alpha, const double* boundary_factor, bool consistent, int begin, int end,
    double* influence, double* influence_tilde, int* error)
{
    const int node = begin+blockIdx.x*blockDim.x+threadIdx.x;
    if (node >= end) return;
    if (!finish_momentum_row(matrix,node,volumes[node],alpha[node],boundary_factor[node],consistent,
                             influence+3*node,influence_tilde+3*node)) atomicExch(error,1);
}
} // namespace mars::segregated
