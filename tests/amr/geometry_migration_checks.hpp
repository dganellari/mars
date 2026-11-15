#pragma once

#include "field_migration_checks.hpp"
#include "public_tet_mesh.hpp"
#include <cstone/sfc/sfc.hpp>
#include <cmath>

namespace amr_test
{
struct GeometryOwner
{
    int ranks, shift;
    HOST_DEVICE_FUN int operator()(uint64_t key) const { return int((key % ranks + shift) % ranks); }
};
inline cstone::Box<double> geometry_box(const std::array<std::vector<double>, 3>& xyz)
{
    double lo[3], hi[3];
    for (int d = 0; d < 3; ++d)
    {
        auto [a,b] = std::minmax_element(xyz[d].begin(), xyz[d].end());
        double padding = .05 * (*b - *a);
        lo[d] = *a - padding;
        hi[d] = *b + padding;
    }
    return {lo[0], hi[0], lo[1], hi[1], lo[2], hi[2]};
}
inline double geometry_area(const std::array<std::vector<double>, 3>& xyz, std::array<uint64_t,3> face)
{
    double a[3], b[3], cross[3];
    for (int d = 0; d < 3; ++d)
    {
        a[d] = xyz[d][face[1]] - xyz[d][face[0]];
        b[d] = xyz[d][face[2]] - xyz[d][face[0]];
    }
    for (int d = 0; d < 3; ++d) cross[d] = a[(d+1)%3]*b[(d+2)%3] - a[(d+2)%3]*b[(d+1)%3];
    return .5 * std::sqrt(cross[0]*cross[0] + cross[1]*cross[1] + cross[2]*cross[2]);
}

template<class Execution>
void check_geometry_migration(MPI_Comm comm, bool staged)
{
    using Key = uint64_t;
    int rank, ranks;
    MPI_Comm_rank(comm, &rank);
    MPI_Comm_size(comm, &ranks);
    auto mesh = read_public_tet(MARS_PUBLIC_TET_MESH);
    const auto original = mesh.mesh.xyz;
    const size_t n = original[0].size();
    auto box = geometry_box(original);
    std::vector<Key> keys(n);
    for (size_t i = 0; i < n; ++i)
        keys[i] = cstone::sfc3D<cstone::SfcKind<Key>>(original[0][i], original[1][i], original[2][i], box);
    auto d_keys = field_upload<Execution>(keys);
    auto physical = original;
    for (int round = 0; round < 20; ++round)
    {
        GeometryOwner old_owner{ranks, round % ranks}, owner{ranks, (round+1) % ranks};
        std::vector<uint8_t> valid(n), owned(n);
        for (size_t i = 0; i < n; ++i)
        {
            valid[i] = old_owner(keys[i]) == rank;
            owned[i] = owner(keys[i]) == rank;
        }
        auto d_valid = field_upload<Execution>(valid), d_owned = field_upload<Execution>(owned);
        mars::amr::AmrFieldMigration<Execution, Key, double> migration(
            Execution::data(d_keys), n, Execution::data(d_valid), Execution::data(d_keys),
            Execution::data(d_owned), n, owner, comm, staged, n, true);
        for (int axis = 0; axis < 3; ++axis)
        {
            auto values = physical[axis];
            for (size_t i = 0; i < n; ++i) if (!valid[i]) values[i] = -1.e30;
            auto d_values = field_upload<Execution>(values);
            typename Execution::template Vector<double> output;
            migration.apply(Execution::data(d_values), n, output, true);
            Execution::download(Execution::data(output), n, values.data());
            mars::amr::require_collectively(values == original[axis], "Coordinate migration changed geometry", comm);
            physical[axis] = std::move(values);
        }
    }
    // The old refit/decode route must fail even with unchanged connectivity and no refinement.
    auto decoded = original;
    std::map<int,double> initial;
    double drift = 0;
    for (int round = 0; round < 20; ++round)
    {
        auto current = geometry_box(decoded);
        auto bits = current.getBoxDimBits(cstone::maxTreeLevel<Key>{});
        double lo[3] = {current.xmin(), current.ymin(), current.zmin()};
        double length[3] = {current.lx(), current.ly(), current.lz()};
        for (size_t i = 0; i < n; ++i)
        {
            auto key = cstone::sfc3D<cstone::SfcKind<Key>>(decoded[0][i], decoded[1][i], decoded[2][i], current);
            auto [ix,iy,iz] = cstone::decodeSfc(key, bits);
            unsigned index[3] = {unsigned(ix), unsigned(iy), unsigned(iz)};
            for (int d = 0; d < 3; ++d)
                decoded[d][i] = lo[d] + index[d] * (1. / ((1u<<bits[d])-1)) * length[d];
        }
        for (const auto& [tag, faces] : mesh.boundary)
        {
            double area = 0;
            for (auto face : faces) area += geometry_area(decoded, face);
            if (!round) initial[tag] = area;
            drift = std::max(drift, std::abs(area/initial[tag]-1));
        }
    }
    mars::amr::require_collectively(drift > 1.e-5, "Old coordinate path did not reproduce surface drift", comm);
    if (rank == 0) std::printf("geometry migration rounds=20 staged=%d exact=1 old_drift=%.17g PASS\n", staged, drift);
}
} // namespace amr_test
