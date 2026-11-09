#include "mars.hpp"
#include "core/mars_simplex_bisection.hpp"

#include <array>
#include <cmath>
#include <iostream>
#include <stdexcept>

void check_single_request_rule()
{
    const unsigned parents[2][4] = {{0, 1, 2, 3}, {0, 2, 1, 4}};
    unsigned children[4][4];
    for (int e = 0; e < 2; ++e)
    {
        int b = e == 0 ? 1 : 2;
        mars::bisection::child_connectivity<4>(parents[e], b, 5u, children[2 * e]);
        mars::bisection::child_connectivity<4>(parents[e], 0, 5u, children[2 * e + 1]);
    }
    int longer_edge_copies = 0;
    for (const auto& child : children)
    {
        bool has_b = false, has_e = false;
        for (unsigned node : child)
        {
            has_b |= node == 1;
            has_e |= node == 4;
        }
        longer_edge_copies += has_b && has_e;
    }
    if (longer_edge_copies != 1) throw std::runtime_error("Requested-edge control changed the longer edge");
    std::cout << "PASS requested-edge control: 4 active tets, 6 nodes; edge (1,4) remains unsplit\n";
}

int main(int argc, char** argv)
{
    Kokkos::ScopeGuard guard(argc, argv);
    check_single_request_rule();
    using Mesh = mars::Mesh<3, 3, mars::KokkosImplementation, mars::Simplex<3, 3, mars::KokkosImplementation>>;
    Mesh mesh;
    mesh.reserve_elements(2);
    mesh.reserve_points(5);
    const double xyz[5][3]  = {{0, 0, 0}, {2, 0, 0}, {1, .5, 0}, {1, .2, .3}, {-.5, 0, -1}};
    const int corners[2][4] = {{0, 1, 2, 3}, {0, 2, 1, 4}};
    auto h_points           = Kokkos::create_mirror_view(mesh.get_view_points());
    auto h_elements         = Kokkos::create_mirror_view(mesh.get_view_elements());
    for (int n = 0; n < 5; ++n)
        for (int d = 0; d < 3; ++d)
            h_points(n, d) = xyz[n][d];
    for (int e = 0; e < 2; ++e)
        for (int c = 0; c < 4; ++c)
            h_elements(e, c) = corners[e][c];
    Kokkos::deep_copy(mesh.get_view_points(), h_points);
    Kokkos::deep_copy(mesh.get_view_elements(), h_elements);
    Kokkos::deep_copy(mesh.get_view_active(), true);

    mars::ViewVectorType<mars::Integer> marked("marked", 1);
    Kokkos::deep_copy(marked, mars::Integer(0));
    mars::ParallelBisection<Mesh> refiner(&mesh);
    refiner.refine(marked);
    Kokkos::fence();
    auto points      = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace{}, mesh.get_view_points());
    auto elements    = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace{}, mesh.get_view_elements());
    auto active      = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace{}, mesh.get_view_active());
    int active_count = 0;
    double volume    = 0;
    for (mars::Integer e = 0; e < mesh.n_elements(); ++e)
    {
        if (!active(e)) continue;
        ++active_count;
        double a[3][3];
        for (int c = 0; c < 3; ++c)
            for (int d = 0; d < 3; ++d)
                a[c][d] = points(elements(e, c + 1), d) - points(elements(e, 0), d);
        double determinant = a[0][0] * (a[1][1] * a[2][2] - a[1][2] * a[2][1]) -
                             a[0][1] * (a[1][0] * a[2][2] - a[1][2] * a[2][0]) +
                             a[0][2] * (a[1][0] * a[2][1] - a[1][1] * a[2][0]);
        if (!(determinant > 0)) throw std::runtime_error("Nonpositive orientation");
        volume += determinant / 6;
    }
    bool terminal_midpoint = false, requested_midpoint = false;
    for (mars::Integer n = 5; n < mesh.n_nodes(); ++n)
    {
        std::cout << "midpoint " << n << ": " << points(n, 0) << ' ' << points(n, 1) << ' ' << points(n, 2) << '\n';
        terminal_midpoint |= points(n, 0) == .75 && points(n, 1) == 0 && points(n, 2) == -.5;
        requested_midpoint |= points(n, 0) == 1 && points(n, 1) == 0 && points(n, 2) == 0;
    }
    if (active_count != 5 || mesh.n_nodes() != 7 || !terminal_midpoint || !requested_midpoint ||
        std::abs(volume - 13. / 60) > 1.e-13)
        throw std::runtime_error("Legacy terminal-edge witness failed");
    std::cout << "PASS legacy LEPP: 5 active tets, 7 nodes, positive volume " << volume
              << "; terminal edge (1,4) before requested edge (0,1)\n";
}
