#pragma once

#include <gtest/gtest.h>
#include "mars.hpp"
#include "mars_structured_gtest.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <iostream>
#include <set>
#include <vector>

namespace mars::test
{

using Cell = std::array<Integer, 3>;

inline int rank()
{
    int value = 0;
    MPI_Comm_rank(MPI_COMM_WORLD, &value);
    return value;
}

inline std::vector<int> sum(const std::vector<int>& local)
{
    std::vector<int> global(local.size());
    EXPECT_EQ(
        MPI_Allreduce(local.data(), global.data(), static_cast<int>(local.size()), MPI_INT, MPI_SUM, MPI_COMM_WORLD),
        MPI_SUCCESS);
    return global;
}

inline bool all_valid(bool local)
{
    int value  = local ? 1 : 0;
    int global = 0;
    EXPECT_EQ(MPI_Allreduce(&value, &global, 1, MPI_INT, MPI_MIN, MPI_COMM_WORLD), MPI_SUCCESS);
    return global != 0;
}

template<Integer Type, class SfcKey>
std::vector<Cell> cells(const ViewVectorType<typename SfcKey::ValueType>& keys)
{
    ViewMatrixType<Integer> coordinates("test_cell_coordinates", keys.extent(0), 3);
    Kokkos::parallel_for(
        "test_decode_cells", keys.extent(0), KOKKOS_LAMBDA(const Integer i) {
            const auto cell   = get_octant_from_sfc<Type, SfcKey>(keys(i));
            coordinates(i, 0) = cell.x;
            coordinates(i, 1) = cell.y;
            coordinates(i, 2) = Type == ElementType::Hex8 ? cell.z : 0;
        });
    const auto host = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), coordinates);
    std::vector<Cell> result(keys.extent(0));
    for (std::size_t i = 0; i < result.size(); ++i)
        result[i] = {host(i, 0), host(i, 1), host(i, 2)};
    return result;
}

inline bool valid_cell(const Cell& cell, const Cell& shape)
{
    for (int d = 0; d < 3; ++d)
    {
        if (cell[d] < 0 || cell[d] >= shape[d]) return false;
    }
    return true;
}

inline Integer cell_index(const Cell& cell, const Cell& shape)
{
    return cell[0] + shape[0] * (cell[1] + shape[1] * cell[2]);
}

template<Integer Type, class SfcKey>
void check_mesh(Integer nx, Integer ny, Integer nz, bool periodic = false)
{
    auto context = make_context(proc_allocation{}, MPI_COMM_WORLD);
    DistributedMesh<SfcKey, Type> mesh(context);
    if (periodic) mesh.set_periodic();
    generate_distributed_cube(mesh, nx, ny, nz);
    const Cell shape{nx, ny, Type == ElementType::Hex8 ? nz : 1};
    const auto local  = cells<Type, SfcKey>(mesh.get_view_sfc());
    const auto ghosts = cells<Type, SfcKey>(mesh.get_view_ghost());
    EXPECT_EQ(local.size(), mesh.get_chunk_size());
    bool valid = true;
    for (const auto& cell : local)
        valid = valid && valid_cell(cell, shape);
    for (const auto& cell : ghosts)
        valid = valid && valid_cell(cell, shape);
    ASSERT_TRUE(all_valid(valid));

    const Integer count = shape[0] * shape[1] * shape[2];
    const int ranks = num_ranks(context), rank_id = rank(context);
    const Integer expected_begin = count / ranks * rank_id + std::min<Integer>(rank_id, count % ranks);
    const Integer expected_end = count / ranks * (rank_id + 1) + std::min<Integer>(rank_id + 1, count % ranks);
    EXPECT_EQ(local.size(), expected_end - expected_begin);
    EXPECT_EQ(mesh.partition_generated_keys(), expected_end - expected_begin);
    EXPECT_EQ(mesh.partition_received_keys(), expected_end - expected_begin);
    using Key = typename SfcKey::ValueType;
    std::vector<Key> reference_keys;
    for (Integer z = 0; z < shape[2]; ++z)
        for (Integer y = 0; y < ny; ++y)
            for (Integer x = 0; x < nx; ++x) {
                if constexpr (Type == ElementType::Quad4) reference_keys.push_back(encode_sfc_2D<SfcKey>(x, y));
                else reference_keys.push_back(encode_sfc_3D<SfcKey>(x, y, z));
            }
    std::sort(reference_keys.begin(), reference_keys.end());
    const auto actual_keys = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace{}, mesh.get_view_sfc());
    for (Integer i = expected_begin; i < expected_end; ++i)
        EXPECT_EQ(actual_keys(i - expected_begin), reference_keys[i]);
    std::vector<int> ownership(count, 0);
    for (const auto& cell : local)
        ++ownership[cell_index(cell, shape)];
    const auto global = sum(ownership);
    for (Integer i = 0; i < count; ++i)
        EXPECT_EQ(global[i], 1) << "Cartesian cell " << i;

    // The halo contains all remote cells sharing a vertex with a local cell.
    std::set<Cell> expected;
    const int z_radius = Type == ElementType::Hex8 ? 1 : 0;
    for (const auto& cell : local)
    {
        for (int dz = -z_radius; dz <= z_radius; ++dz)
        {
            for (int dy = -1; dy <= 1; ++dy)
            {
                for (int dx = -1; dx <= 1; ++dx)
                {
                    Cell neighbor{cell[0] + dx, cell[1] + dy, cell[2] + dz};
                    if (periodic)
                    {
                        for (int d = 0; d < 3; ++d)
                            neighbor[d] = (neighbor[d] + shape[d]) % shape[d];
                    }
                    if (valid_cell(neighbor, shape) && !ownership[cell_index(neighbor, shape)])
                        expected.insert(neighbor);
                }
            }
        }
    }
    const std::set<Cell> actual(ghosts.begin(), ghosts.end());
    EXPECT_EQ(actual.size(), ghosts.size()) << "Duplicate ghost cell";
    EXPECT_EQ(actual, expected);

    UserData<decltype(mesh), double, int> cell_data(mesh);
    const auto mesh_view = mesh.device_view();
    const auto cell_values = cell_data.template get_data<0>();
    const auto cell_integer = cell_data.template get_data<1>();
    for (int step = 0; step < 3; ++step) {
        Kokkos::parallel_for("cell_field_refresh", mesh.get_chunk_size(), KOKKOS_LAMBDA(Integer i) {
            const auto cell = mesh_view.get_octant(i);
            const Integer linear = cell.x + nx * (cell.y + ny * cell.z);
            cell_values(i) = linear + step + 0.5;
            cell_integer(i) = linear + step;
        });
        cell_data.exchange_ghost_data(context);
        const auto ghost_values = cell_data.template get_ghost_data<0>();
        const auto ghost_integer = cell_data.template get_ghost_data<1>();
        int errors = 0;
        Kokkos::parallel_reduce("cell_field_halo_check", mesh.get_ghost_size(),
            KOKKOS_LAMBDA(Integer i, int& n) {
                const auto cell = mesh_view.get_ghost_octant(i);
                const Integer linear = cell.x + nx * (cell.y + ny * cell.z);
                n += ghost_values(i) != linear + step + 0.5;
                n += ghost_integer(i) != linear + step;
            }, errors);
        EXPECT_EQ(errors, 0);
    }
    if (periodic) return;
    EXPECT_TRUE(mesh.template generate_points<Type>());
    EXPECT_TRUE(mesh.template generate_elements<Type>());
    const auto points   = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), mesh.get_view_points());
    const auto elements = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), mesh.get_view_elements());
    const auto active   = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), mesh.get_view_active());
    ASSERT_EQ(elements.extent(0), local.size());
    ASSERT_EQ(elements.extent(1), Type);
    const int dim = Type == ElementType::Hex8 ? 3 : 2;
    // Reference corners run counterclockwise on z=0, then in the same order on z=1.
    const std::array<Cell, 8> reference{{{0, 0, 0}, {1, 0, 0}, {1, 1, 0}, {0, 1, 0},
                                       {0, 0, 1}, {1, 0, 1}, {1, 1, 1}, {0, 1, 1}}};
    const double expected_jacobian = 1.0 / nx / ny / (dim == 3 ? nz : 1);
    for (std::size_t i = 0; i < local.size(); ++i)
    {
        EXPECT_TRUE(active(i));
        std::set<Cell> vertices;
        std::array<std::array<double, 3>, 8> physical{};
        for (Integer j = 0; j < Type; ++j)
        {
            const auto vertex = elements(i, j);
            ASSERT_GE(vertex, 0);
            ASSERT_LT(vertex, points.extent(0));
            Cell coordinate{0, 0, 0};
            for (int d = 0; d < dim; ++d)
            {
                physical[j][d]    = points(vertex, d);
                const auto scaled = physical[j][d] * shape[d];
                coordinate[d]     = std::llround(scaled);
                EXPECT_NEAR(scaled, coordinate[d], 1e-12);
                EXPECT_EQ(coordinate[d], local[i][d] + reference[j][d]) << "Element " << i << ", corner " << j;
            }
            vertices.insert(coordinate);
        }
        EXPECT_EQ(vertices.size(), Type);
        // The affine map from [0,1]^dim has columns along edges 0->1, 0->3 and 0->4.
        std::array<double, 3> dx{}, dy{}, dz{};
        for (int d = 0; d < dim; ++d)
        {
            dx[d] = physical[1][d] - physical[0][d];
            dy[d] = physical[3][d] - physical[0][d];
            if (dim == 3) dz[d] = physical[4][d] - physical[0][d];
        }
        const double jacobian = dim == 2 ? dx[0] * dy[1] - dx[1] * dy[0] :
            dx[0] * (dy[1] * dz[2] - dy[2] * dz[1]) -
            dx[1] * (dy[0] * dz[2] - dy[2] * dz[0]) +
            dx[2] * (dy[0] * dz[1] - dy[1] * dz[0]);
        EXPECT_GT(jacobian, 0.0);
        EXPECT_NEAR(jacobian, expected_jacobian, 1e-12 * expected_jacobian);
    }
}

template<class Handler>
auto dof_snapshot(const Handler& owner)
{
    const auto handler = owner.device_view();
    ViewMatrixType<Integer> indices("test_dof_indices", handler.get_dof_size(), 6);
    ViewMatrixType<double> coordinates("test_dof_coordinates", handler.get_dof_size(), 3);
    Kokkos::parallel_for(
        "test_dof_snapshot", handler.get_dof_size(), KOKKOS_LAMBDA(const Integer i) {
            double point[3] = {0.0, 0.0, 0.0};
            handler.get_local_dof_coordinates(i, point);
            for (int d = 0; d < 3; ++d)
                coordinates(i, d) = point[d];
            const auto key       = handler.local_to_sfc(i);
            const auto component = handler.get_component(i);
            indices(i, 0)        = handler.local_to_global(i);
            indices(i, 1)        = handler.local_to_global_proc(i);
            indices(i, 2)        = handler.is_owned(i);
            indices(i, 3)        = handler.sfc_to_local(key, component);
            const auto owned     = handler.local_to_owned(i);
            indices(i, 4)        = handler.is_owned(i) ? handler.get_owned_dof(owned) : -1;
            indices(i, 5)        = component;
        });
    return std::make_pair(Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), indices),
                          Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), coordinates));
}

template<Integer Type, Integer Degree, class SfcKey, Integer Block = 1>
void check_dofs(Integer nx, Integer ny, Integer nz, bool halo = false, Integer runtime_block = 0)
{
    auto context   = make_context(proc_allocation{}, MPI_COMM_WORLD);
    using MeshType = DistributedMesh<SfcKey, Type>;
    MeshType mesh(context);
    generate_distributed_cube(mesh, nx, ny, nz);
    using Handler = DofHandler<MeshType, Degree, Block>;
    Handler handler(mesh);
    if constexpr (Block == 0) handler.set_block(runtime_block);
    handler.enumerate_dofs();
    const Integer block = handler.get_block();
    const Cell lattice{Degree * nx + 1, Degree * ny + 1, Type == ElementType::Hex8 ? Degree * nz + 1 : 1};
    const Integer expected_count = lattice[0] * lattice[1] * lattice[2] * block;
    EXPECT_EQ(handler.get_global_dof_size(), expected_count);
    const auto snapshot       = dof_snapshot(handler);
    const auto& indices       = snapshot.first;
    const auto& coordinates   = snapshot.second;
    const Integer local_count = handler.get_dof_size();
    int ranks                 = 0;
    MPI_Comm_size(MPI_COMM_WORLD, &ranks);
    std::vector<int> owner_count(expected_count, 0), gid_count(expected_count, 0), presence(expected_count, 0);
    std::vector<int> owner_rank(expected_count, 0), owner_gid(expected_count, 0);
    std::vector<Integer> oracle_id(local_count, -1);
    bool valid = true;
    for (Integer i = 0; i < local_count; ++i)
    {
        Cell node{0, 0, 0};
        for (int d = 0; d < 3; ++d)
        {
            const auto scaled = coordinates(i, d) * (lattice[d] - 1);
            node[d]           = std::llround(scaled);
            EXPECT_NEAR(scaled, node[d], 1e-11);
        }
        valid = valid && valid_cell(node, lattice) && indices(i, 0) >= 0 && indices(i, 0) < expected_count &&
                indices(i, 1) >= 0 && indices(i, 1) < ranks && indices(i, 5) >= 0 && indices(i, 5) < block;
        if (!valid) continue;
        const Integer id = block * cell_index(node, lattice) + indices(i, 5);
        oracle_id[i]     = id;
        ++presence[id];
        EXPECT_EQ(indices(i, 3), i);
        if (indices(i, 2))
        {
            EXPECT_EQ(indices(i, 4), i);
            EXPECT_EQ(indices(i, 1), rank());
            ++owner_count[id];
            ++gid_count[indices(i, 0)];
            owner_rank[id] = rank() + 1;
            owner_gid[id]  = indices(i, 0) + 1;
        }
    }
    ASSERT_TRUE(all_valid(valid));
    for (const auto count : presence)
        EXPECT_LE(count, 1) << "Duplicate local coordinate/component";
    const auto global_owners     = sum(owner_count);
    const auto global_ids        = sum(gid_count);
    const auto global_owner_rank = sum(owner_rank);
    const auto global_owner_gid  = sum(owner_gid);
    for (Integer i = 0; i < expected_count; ++i)
    {
        EXPECT_EQ(global_owners[i], 1) << "Lattice DOF " << i;
        EXPECT_EQ(global_ids[i], 1) << "Global ID " << i;
    }
    for (Integer i = 0; i < local_count; ++i)
    {
        EXPECT_EQ(indices(i, 0) + 1, global_owner_gid[oracle_id[i]]);
        EXPECT_EQ(indices(i, 1) + 1, global_owner_rank[oracle_id[i]]);
    }

    if (halo)
    {
        DM<Handler, double, double> data(handler);
        auto values = data.template get_dof_data<0>();
        auto second = data.template get_dof_data<1>();
        EXPECT_EQ(values.extent(0), local_count);
        EXPECT_EQ(second.extent(0), local_count);
        ASSERT_TRUE(all_valid(values.extent(0) == local_count && second.extent(0) == local_count));
        const auto& const_data = data;
        const auto device_handler = handler.device_view();
        const_data.template data_iterate<0>(KOKKOS_LAMBDA(double& value) { value = -1234.0; });
        const_data.template data_iterate<1>(KOKKOS_LAMBDA(double& value) { value = -1234.0; });
        const_data.template owned_data_iterate<0>(KOKKOS_LAMBDA(const Integer i, double& value) {
                double point[3] = {0.0, 0.0, 0.0};
                device_handler.get_local_dof_coordinates(i, point);
                value = 2.0 + point[0] + 3.0 * point[1] + 7.0 * point[2] + 11.0 * device_handler.get_component(i);
            });
        const_data.template owned_data_iterate<1>(KOKKOS_LAMBDA(const Integer i, double& value) {
                double point[3] = {0.0, 0.0, 0.0};
                device_handler.get_local_dof_coordinates(i, point);
                value = -5.0 + 2.0 * point[0] - point[1] + 4.0 * point[2] - 13.0 * device_handler.get_component(i);
            });
        auto expected_first = [&](Integer i) {
            return 2.0 + coordinates(i, 0) + 3.0 * coordinates(i, 1) + 7.0 * coordinates(i, 2) +
                   11.0 * indices(i, 5);
        };
        auto expected_second = [&](Integer i) {
            return -5.0 + 2.0 * coordinates(i, 0) - coordinates(i, 1) + 4.0 * coordinates(i, 2) -
                   13.0 * indices(i, 5);
        };
        data.template start_gather<0>();
        EXPECT_TRUE(data.exchange_plan().active());
        EXPECT_THROW(data.template start_gather<0>(), std::logic_error);
        data.wait();
        EXPECT_FALSE(data.exchange_plan().active());
        EXPECT_EQ(data.exchange_plan().buffer_allocations(), 1);
        auto host_values = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), values);
        auto host_second = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), second);
        for (Integer i = 0; i < local_count; ++i)
        {
            EXPECT_NEAR(host_values(i), expected_first(i), 1e-12);
            EXPECT_NEAR(host_second(i), indices(i, 2) ? expected_second(i) : -1234.0, 1e-12);
        }
        data.gather_ghost_data();
        host_second = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), second);
        for (Integer i = 0; i < local_count; ++i) EXPECT_NEAR(host_second(i), expected_second(i), 1e-12);

        EXPECT_EQ(data.exchange_plan().buffer_allocations(), 2);
        ViewVectorType<double> independent("independent_halo_work", local_count);
        for (int step = 1; step <= 4; ++step) {
            Kokkos::parallel_for("refresh_halo_source", local_count, KOKKOS_LAMBDA(Integer i) {
                double point[3] = {};
                device_handler.get_local_dof_coordinates(i, point);
                values(i) = device_handler.is_owned(i)
                    ? 2.0 + point[0] + 3.0 * point[1] + 7.0 * point[2] +
                      11.0 * device_handler.get_component(i) + step : -1234.0;
            });
            data.template start_gather<0>();
            Kokkos::parallel_for("independent_halo_work", local_count, KOKKOS_LAMBDA(Integer i) {
                independent(i) = i + step;
            });
            data.wait();
            auto refreshed = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace{}, values);
            for (Integer i = 0; i < local_count; ++i) EXPECT_NEAR(refreshed(i), expected_first(i) + step, 1e-12);
            EXPECT_EQ(data.exchange_plan().buffer_allocations(), 2);
        }
        Kokkos::parallel_for("restore_halo_values", local_count, KOKKOS_LAMBDA(Integer i) { values(i) -= 4; });

        ViewVectorType<double> owned("test_owned_values", handler.get_owned_dof_size());
        ViewVectorType<Integer> owned_indices("test_owned_indices", handler.get_owned_dof_size());
        Kokkos::parallel_for("test_owned_indices", handler.get_owned_dof_size(), KOKKOS_LAMBDA(const Integer i) {
            owned_indices(i) = device_handler.get_owned_dof(i);
        });
        data.template get_locally_owned_data<0>(owned);
        const auto host_owned = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), owned);
        const auto host_owned_indices = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), owned_indices);
        for (Integer i = 0; i < handler.get_owned_dof_size(); ++i)
            EXPECT_NEAR(host_owned(i), expected_first(host_owned_indices(i)), 1e-12);
        Kokkos::parallel_for("test_modify_owned", owned.extent(0), KOKKOS_LAMBDA(const Integer i) { owned(i) += 17.0; });
        data.template set_locally_owned_data<0>(owned);
        host_values = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), values);
        for (Integer i = 0; i < local_count; ++i)
        {
            EXPECT_NEAR(host_values(i), expected_first(i) + (indices(i, 2) ? 17.0 : 0.0), 1e-12);
        }

        // Each local copy contributes once; the owner must receive the exact rank sum.
        const int contribution = rank() + 1;
        std::vector<int> local_contributions(expected_count, 0);
        for (Integer i = 0; i < local_count; ++i)
            local_contributions[oracle_id[i]] = contribution * (indices(i, 5) + 1);
        const auto expected_sum = sum(local_contributions);
        Kokkos::parallel_for("test_scatter_values", local_count, KOKKOS_LAMBDA(const Integer i) {
            values(i) = contribution * (device_handler.get_component(i) + 1);
            second(i) = 2.0 * values(i);
        });
        const auto selected = data.template scatter_ghost_data<1>();
        EXPECT_EQ(std::get<0>(selected).extent(0), 0);
        EXPECT_EQ(std::get<1>(selected).extent(0), handler.get_boundary_dof_size());
        host_values = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), values);
        host_second = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), second);
        for (Integer i = 0; i < local_count; ++i)
        {
            EXPECT_DOUBLE_EQ(host_values(i), contribution * (indices(i, 5) + 1));
            EXPECT_DOUBLE_EQ(host_second(i), 2 * contribution * (indices(i, 5) + 1));
        }
        data.template scatter_add_ghost_data<1>();
        data.template gather_ghost_data<1>();
        host_values = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), values);
        host_second = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), second);
        for (Integer i = 0; i < local_count; ++i)
        {
            EXPECT_DOUBLE_EQ(host_values(i), contribution * (indices(i, 5) + 1));
            EXPECT_DOUBLE_EQ(host_second(i), 2 * expected_sum[oracle_id[i]]);
        }
        Kokkos::parallel_for("test_reset_scatter_values", local_count, KOKKOS_LAMBDA(const Integer i) {
            values(i) = contribution * (device_handler.get_component(i) + 1);
            second(i) = 2.0 * values(i);
        });
        data.scatter_add_ghost_data();
        data.gather_ghost_data();
        host_values = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), values);
        host_second = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), second);
        for (Integer i = 0; i < local_count; ++i)
        {
            EXPECT_DOUBLE_EQ(host_values(i), expected_sum[oracle_id[i]]);
            EXPECT_DOUBLE_EQ(host_second(i), 2 * expected_sum[oracle_id[i]]);
        }
    }
}

} // namespace mars::test
