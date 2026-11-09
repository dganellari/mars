#pragma once

#include "tet_checks.hpp"
#include "backend/distributed/unstructured/domain.hpp"

namespace amr_test
{
using Domain = mars::ElementDomain<mars::TetTag, RealType, KeyType, cstone::execution::Gpu>;
using Values = cstone::DeviceVector<RealType>;

template<class T>
std::vector<T> read_device(const T* data, size_t size)
{
    std::vector<T> result(size);
    if (size) cudaCheckCall(cudaMemcpy(result.data(), data, size * sizeof(T), cudaMemcpyDeviceToHost));
    return result;
}

template<class T>
cstone::DeviceVector<T> upload(const std::vector<T>& data)
{
    cstone::DeviceVector<T> result(data.size());
    if (!data.empty())
        cudaCheckCall(cudaMemcpy(result.data(), data.data(), data.size() * sizeof(T), cudaMemcpyHostToDevice));
    return result;
}

std::unique_ptr<Domain> make_cube(int rank, int ranks);
std::vector<Tet> snapshot(const std::array<const KeyType*, 4>& conn,
                          size_t start,
                          size_t count,
                          const RealType* x,
                          const RealType* y,
                          const RealType* z,
                          size_t nodes,
                          const KeyType* keys,
                          int rank);
std::vector<Tet> snapshot(Domain& domain, int rank);
std::array<Values, 2> initial_fields(Domain& domain);
void check_fields(const RealType* x,
                  const RealType* y,
                  const RealType* z,
                  size_t nodes,
                  const std::array<Values, 2>& fields,
                  RealType tolerance,
                  const std::vector<uint8_t>* current = nullptr);
std::vector<uint8_t> current_nodes(Domain& domain);
void check_sparse_device_exchange(int rank, int ranks);
void check_domain_layout(Domain& domain, int rank, int ranks);
void poison_ghost_fields(Domain& domain, std::array<Values, 2>& fields);
std::vector<Node> check_identity(Domain& domain, int rank);
void check_forced_migration(Domain& domain, int rank, int ranks);
RealType affine_quantization_bound(const Domain& domain);
void check_coarsening_ignored(Domain& domain, int rank);
void check_single_request(Domain& domain, int rank, int ranks, bool staged);
void check_default_red(Domain& domain, int rank, int ranks);
void label_boundary(Domain& domain, int rank);
void check_boundary(Domain& domain, const std::vector<Tet>& global);
void check_manager(std::unique_ptr<Domain> domain, int rank, int ranks, bool staged);
void check_dependent_device_lepp(int rank, int ranks, bool staged);
} // namespace amr_test
