#include "tet_gpu_checks.hpp"

namespace amr_test
{
void run_gpu_checks(int rank, int ranks, bool staged)
{
    check_sparse_device_exchange(rank, ranks);
    check_dependent_device_lepp(rank, ranks, staged);
    auto domain = make_cube(rank, ranks);
    domain->cacheNodeCoordinates();
    check_domain_layout(*domain, rank, ranks);
    (void)domain->getNodeOwnershipMap();
    (void)domain->getElementToNodeConnectivity();
    check_forced_migration(*domain, rank, ranks);
    check_coarsening_ignored(*domain, rank);
    check_default_red(*domain, rank, ranks);
    check_single_request(*domain, rank, ranks, staged);
    check_manager(std::move(domain), rank, ranks, staged);
    cudaCheckCall(cudaDeviceSynchronize());
}
} // namespace amr_test
