#pragma once

#include <mars.hpp>
#include <type_traits>

namespace mars::test
{

inline ViewMatrixType<Real> make_preflight_points()
{
    ViewMatrixType<Real> values("preflight_points", 2, 2);
    Kokkos::parallel_for("preflight_point_fill", 2, KOKKOS_LAMBDA(const int i) {
        values(0, i) = 3 + i;
        values(1, i) = 1 + i;
    });
    return values;
}

template <Integer Type>
inline int check_structured_element_wrapper(const context& owner)
{
    constexpr Integer dim = Type == ElementType::Quad4 ? 2 : 3;
    using Key = MortonKey<Unsigned>;
    DMesh<dim, dim, Type, Key> mesh(owner);
    generate_distributed_cube<Key, dim, dim, Type, false>(mesh, 4, 3, dim == 2 ? 0 : 2);
    const auto connectivity = mesh.get_view_elements();
    auto first = mesh.elem(0);
    const auto& const_mesh = mesh;
    auto last = const_mesh.elem(mesh.n_elements() - 1);
    const Integer last_id = mesh.n_elements() - 1;
    int errors = 0;
    Kokkos::parallel_reduce("preflight_element_wrapper", Type,
        KOKKOS_LAMBDA(const int i, int& failures) {
            failures += first.get_elem_id() != 0 || last.get_elem_id() != last_id;
            failures += first.node(i) != connectivity(0, i);
            failures += last.node(i) != connectivity(last_id, i);
        }, errors);
    return errors;
}

inline int check_structured_view_lifetimes(MPI_Comm comm)
{
    using Mesh = DistributedQuad4Mesh<MortonKey<Unsigned>>;
    using Data = UserData<Mesh, double>;
    static_assert(std::is_constructible_v<Mesh, context&>);
    static_assert(!std::is_constructible_v<Mesh, context&&>);
    static_assert(!std::is_constructible_v<Mesh, const context&&>);
    static_assert(std::is_same_v<decltype(std::declval<const Data&>().get_view_boundary()), ViewVectorType<Unsigned>>);
    static_assert(std::is_same_v<decltype(std::declval<const Data&>().get_view_ghost()), ViewVectorType<Unsigned>>);
    static_assert(std::is_same_v<decltype(std::declval<const Data&>().get_view_scan_ghost()), ViewVectorType<Integer>>);

    // Both source View handles disappear before their wrappers reach the kernel.
    Point<Real, 2> point(make_preflight_points(), 0);
    auto row = [] {
        auto values = make_preflight_points();
        return SubView<Real, 2>(&values, 1);
    }();
    auto row_alias = row;
    int errors = 0;
    Kokkos::parallel_reduce("preflight_point_lifetime", 1,
        KOKKOS_LAMBDA(const int, int& failures) {
            failures += point[0] != 3 || point[1] != 4 || point.squared_norm() != 25;
            failures += row[0] != 1 || row[1] != 2 || dot(row, row) != 5;
            auto writable = row;
            writable += row_alias;
            writable *= 3;
            writable /= 2;
            failures += writable[0] != 3 || writable[1] != 6;
        }, errors);

    ViewVectorType<Unsigned> boundary, ghosts;
    ViewVectorType<Integer> scan;
    {
        auto owner = make_context(proc_allocation{}, comm);
        Mesh mesh(owner);
        mesh.set_chunk_size(0);
        mesh.set_view_boundary(ViewVectorType<Unsigned>("preflight_boundary", 3));
        mesh.reserve_ghost(2);
        mesh.reserve_scan_ghost(3);
        const auto mesh_boundary = mesh.get_view_boundary();
        const auto mesh_ghosts = mesh.get_view_ghost();
        const auto mesh_scan = mesh.get_view_scan_ghost();
        Kokkos::parallel_for("preflight_userdata_fill", 3, KOKKOS_LAMBDA(const int i) {
            mesh_boundary(i) = 7 + i;
            mesh_scan(i) = i;
            if (i < 2) mesh_ghosts(i) = 11 + i;
        });
        Data data(mesh);
        boundary = data.get_view_boundary();
        ghosts = data.get_view_ghost();
        scan = data.get_view_scan_ghost();

        ViewVectorType<Unsigned> duplicate_keys("preflight_duplicate_keys", 2);
        Kokkos::deep_copy(duplicate_keys, Unsigned(1));
        Mesh duplicate_mesh(owner);
        duplicate_mesh.set_chunk_size(2);
        duplicate_mesh.set_view_sfc(duplicate_keys);
        errors += duplicate_mesh.generate_sfc_to_local_map();
        SFC<ElementType::Quad4, MortonKey<Unsigned>> duplicate_dofs(2, 2, 0);
        duplicate_dofs.set_elements(duplicate_keys);
        errors += duplicate_dofs.generate_sfc_to_local_map();

        errors += check_structured_element_wrapper<ElementType::Quad4>(owner);
        errors += check_structured_element_wrapper<ElementType::Hex8>(owner);
    }
    int getter_errors = 0;
    Kokkos::parallel_reduce("preflight_userdata_lifetime", 3,
        KOKKOS_LAMBDA(const int i, int& failures) {
            failures += boundary(i) != Unsigned(7 + i) || scan(i) != i;
            if (i < 2) failures += ghosts(i) != Unsigned(11 + i);
        }, getter_errors);
    return errors + getter_errors;
}

} // namespace mars::test
