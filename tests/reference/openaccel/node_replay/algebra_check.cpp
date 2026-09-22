#include "mars_segregated_node.hpp"
#include <algorithm>
#include <cmath>
#include <iostream>
#include <limits>
#include <stdexcept>

int checks = 0;
void check(bool value) { ++checks; if (!value) throw std::runtime_error("node algebra check failed"); }
void close(double got, double expected) { check(std::isfinite(got) && std::abs(got-expected) <= 1e-13*std::max(1.0, std::abs(expected))); }
int main()
{
    using namespace mars::segregated;
    try {
        // rho V/dt=12 kg/s. Prescribed integrated forces and divergence are independent.
        SteadyMomentumNode x{2, 3, .5, -4, {1, -2, 3}, {2, 4, -1}, {1, 2, 3}, {-1, 1, 2}};
        double lhs[9], rhs[3];
        steady_momentum_node(x, lhs, rhs);
        for (int i = 0; i < 9; ++i) close(lhs[i], i%4 == 0 ? 16 : 0);
        close(rhs[0], -10); close(rhs[1], 5); close(rhs[2], 6);
        x.mass_divergence = 4;
        steady_momentum_node(x, lhs, rhs);
        close(lhs[0], 12); close(rhs[0], -2); close(rhs[1], -11); close(rhs[2], 30);
        x.mass_divergence = 0;
        steady_momentum_node(x, lhs, rhs);
        close(lhs[0], 12); close(rhs[0], -6);
        x.pseudo_dt = .25;
        steady_momentum_node(x, lhs, rhs);
        close(lhs[0], 24); close(rhs[0], -6); // pseudo-time has no steady RHS
        double block[9] = {4, 900, -700, 300, 8, 400, -800, 500, 12};
        relax_momentum_diagonal(block, .25);
        const double relaxed[9] = {16, 900, -700, 300, 32, 400, -800, 500, 48};
        for (int i = 0; i < 9; ++i) close(block[i], relaxed[i]);
        // Diagonal is the middle block. Cross-component neighbor entries must be ignored.
        double rows[27] = {-2, 1e8, -1e8, 1e8, -4, 1e8, 1e8, 1e8, -6,
                          16, 900, -700, 300, 32, 400, -800, 500, 48,
                          -2, -1e8, 1e8, -1e8, -4, -1e8, -1e8, -1e8, -6};
        double d[3], dt[3];
        check(momentum_influence(8, rows, 3, 1, true, d, dt));
        close(d[0], .5); close(d[1], .25); close(d[2], 1.0/6);
        close(dt[0], 2.0/3); close(dt[1], 1.0/3); close(dt[2], 2.0/9);
        check(momentum_influence(8, rows, 3, 1, false, d, dt));
        for (double v : dt) close(v, 0);
        // Only one alpha: d=alpha*V/raw_diagonal, not alpha squared.
        close(d[0], .25*8/4);
        rows[9] = -16;
        check(momentum_influence(8, rows, 3, 1, true, d, dt));
        close(d[0], -.5); close(dt[0], -.4);
        rows[9] = -std::numeric_limits<double>::epsilon();
        check(!momentum_influence(8, rows, 3, 1, false, d, dt));
        rows[9] = std::numeric_limits<double>::infinity();
        check(!momentum_influence(8, rows, 3, 1, false, d, dt));
        rows[9] = std::numeric_limits<double>::quiet_NaN();
        check(!momentum_influence(8, rows, 3, 1, false, d, dt));
        rows[9] = 16;
        check(!momentum_influence(8, rows, 3, 3, false, d, dt));
        check(!momentum_influence(8, nullptr, 3, 1, false, d, dt));
        check(!momentum_influence(std::numeric_limits<double>::infinity(), rows, 3, 1, false, d, dt));
        rows[0] = std::numeric_limits<double>::infinity();
        check(!momentum_influence(8, rows, 3, 1, true, d, dt));
        double boundary[3] = {8, -4, 0};
        relax_boundary_rhs(boundary, .75);
        close(boundary[0], 6); close(boundary[1], -3); close(boundary[2], 0);
        std::cout << "PASS: " << checks << " independent node algebra checks\n";
    } catch (const std::exception& e) { std::cerr << e.what() << '\n'; return 1; }
}
